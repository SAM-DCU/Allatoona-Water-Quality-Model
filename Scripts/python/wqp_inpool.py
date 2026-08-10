#!/usr/bin/env python3
"""Pull and summarize EPA Water Quality Portal in-pool data for both reservoirs.

Downloads every lake result within a bounding box around each reservoir for the
characteristics that matter to this study, saves the raw portal CSV as
``inpool_<lake>.csv`` in the output directory, and prints a per-characteristic
summary: result count, number of stations and dates, date range, depth range, and
units. The saved files are the input to forebay_profiles.py and
tailrace_and_stations.py, so give all three the same directory.

The portal CSV is written exactly as delivered; nothing is filtered, converted, or
renamed on the way to disk, so the raw record stays intact and every later step is a
documented transformation of it. Only the printed summary converts anything, and the
one conversion it makes is depth to metres.

Provenance
----------
Source: EPA Water Quality Portal Result service, ``resultPhysChem`` profile, queried by
bounding box with ``siteType=Lake, Reservoir, Impoundment``. The portal aggregates
state and federal data; at these two reservoirs the profile records come from GA EPD
under the organization identifier ``21GAEPD_WQX``. The full query URL and the retrieval
time are printed on every run, because the portal has no way to re-request a past
snapshot and the returned record grows and is revised over time.

Usage:
    /opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3 \
        Scripts/python/wqp_inpool.py <output_dir>
"""
import argparse
import csv
import datetime as dt
import io
import os
import statistics
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from typing import Optional

# Bounding boxes, decimal degrees WGS84, as west,south,east,north. Each was drawn to
# enclose the reservoir pool and no more; they are query extents, not project boundaries.
LAKES = {
    "Allatoona": "-84.78,34.09,-84.58,34.23",
    "Lanier":    "-84.10,34.13,-83.80,34.45",
}
CHARS = [
    "Dissolved oxygen (DO)",
    "Temperature, water",
    "Phosphorus",
    "Nitrogen",
    "Chlorophyll a",
    "Depth, Secchi disk depth",
]
BASE = "https://www.waterqualitydata.us/data/Result/search"
TIMEOUT_S = 300

# Depth arrives in the portal's own unit codes. Metres dominate, but a few hundred Lanier
# activity depths are coded feet (585 in the current pull, none of them at either forebay
# station), so the summary converts before it prints a range in metres. Any other code is
# left unconverted and counted, so a new one cannot enter the range silently.
# 0.3048 m/ft is exact by definition.
FT_TO_M = 0.3048
DEPTH_UNITS_M = {"m", "meters", "metres"}
DEPTH_UNITS_FT = {"ft", "feet"}

DEPTH_COLS = [("ResultDepthHeightMeasure/MeasureValue",
               "ResultDepthHeightMeasure/MeasureUnitCode"),
              ("ActivityDepthHeightMeasure/MeasureValue",
               "ActivityDepthHeightMeasure/MeasureUnitCode")]


def build_url(bbox: str) -> str:
    """Return the portal query URL for one bounding box.

    Parameters
    ----------
    bbox : str
        Bounding box as "west,south,east,north" in decimal degrees.

    Returns
    -------
    str
        Fully encoded Water Quality Portal Result service URL.
    """
    params = [
        ("bBox", bbox),
        ("siteType", "Lake, Reservoir, Impoundment"),
        ("mimeType", "csv"),
        ("zip", "no"),
        ("dataProfile", "resultPhysChem"),
    ] + [("characteristicName", c) for c in CHARS]
    return BASE + "?" + urllib.parse.urlencode(params)


def fetch(lake: str, bbox: str, out_dir: str) -> str:
    """Download one lake's in-pool results and save the portal CSV verbatim.

    Parameters
    ----------
    lake : str
        Lake name, used only in the output file name and log lines.
    bbox : str
        Bounding box as "west,south,east,north" in decimal degrees.
    out_dir : str
        Directory to write ``inpool_<lake>.csv`` into.

    Returns
    -------
    str
        Path to the saved CSV.

    Raises
    ------
    urllib.error.URLError
        If the request fails or times out.
    ValueError
        If the response is not the expected result CSV, which is how the portal
        signals an error page or a changed profile rather than by an HTTP status.
    """
    url = build_url(bbox)
    path = os.path.join(out_dir, f"inpool_{lake}.csv")
    print(f"[{lake}] GET {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "data-inventory"})
    with urllib.request.urlopen(req, timeout=TIMEOUT_S) as r:
        data = r.read()
    head = data[:2048].decode("utf-8", "replace")
    if "CharacteristicName" not in head:
        raise ValueError(f"{lake}: response is not a resultPhysChem CSV; "
                         f"first bytes were: {head[:200]!r}")
    with open(path, "wb") as f:
        f.write(data)
    print(f"[{lake}] saved {len(data):,} bytes -> {path} "
          f"(retrieved {dt.datetime.now().astimezone().isoformat(timespec='seconds')})",
          flush=True)
    return path


def num(x: Optional[str]) -> Optional[float]:
    """Return x as a float, or None if it is blank or not numeric.

    Portal values can be blank, a detection-limit qualifier, or free text, so a
    non-numeric entry is a normal outcome rather than an error.
    """
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def depth_m(row: dict) -> Optional[float]:
    """Sample depth in metres from a portal result row.

    Parameters
    ----------
    row : dict
        One row of the portal result CSV.

    Returns
    -------
    float or None
        Depth below surface in metres, taken from the result depth if present and
        otherwise from the activity depth. A blank unit code is taken as metres. None if
        neither is populated, or if the unit code is one this function does not convert,
        in which case the caller counts the row and reports it.
    """
    for vcol, ucol in DEPTH_COLS:
        v = num(row.get(vcol, ""))
        if v is None:
            continue
        u = (row.get(ucol, "") or "").strip().lower()
        if u in DEPTH_UNITS_FT:
            return v * FT_TO_M
        if u in DEPTH_UNITS_M or u == "":
            return v
        return None
    return None


def summarize(lake: str, path: str) -> None:
    """Print a per-characteristic inventory of one saved portal CSV.

    Reports result count, distinct stations and dates, date range, depth range in
    metres, and the reported unit codes, then the vertical structure of the dissolved
    oxygen profiles. Characteristics that were requested but returned nothing are named
    explicitly, and so is any depth whose unit code could not be converted, so that a
    silent absence is never mistaken for a real absence of data.
    """
    rows = list(csv.DictReader(io.open(path, encoding="utf-8", errors="replace")))
    print(f"\n===== {lake}: {len(rows):,} result rows =====")
    if not rows:
        return
    by = defaultdict(lambda: {"n": 0, "stns": set(), "dates": set(), "depths": [],
                              "profiles": set(), "units": set(), "depth_units": set(),
                              "depth_unconverted": 0})
    for r in rows:
        c = r.get("CharacteristicName", "?")
        d = by[c]
        d["n"] += 1
        stn = r.get("MonitoringLocationIdentifier", "")
        date = r.get("ActivityStartDate", "")
        d["stns"].add(stn)
        d["dates"].add(date)
        d["profiles"].add((stn, date))
        u = r.get("ResultMeasure/MeasureUnitCode", "")
        if u:
            d["units"].add(u)
        z = depth_m(r)
        if z is not None:
            d["depths"].append(z)
        # Only the depth column that depth_m would have used is credited, so the reported
        # unit codes describe the depths in the range above rather than every code in the
        # file. The break keeps the two in step.
        for vcol, ucol in DEPTH_COLS:
            if num(r.get(vcol, "")) is not None:
                d["depth_units"].add((r.get(ucol, "") or "").strip())
                if z is None:
                    d["depth_unconverted"] += 1
                break
    print(f"{'characteristic':<28}{'n':>7}{'stns':>6}{'dates':>7}  {'date range':<24}"
          f"{'depth(m) min/max':<18}units")
    for c in sorted(by, key=lambda k: -by[k]["n"]):
        d = by[c]
        ds = sorted(x for x in d["dates"] if x)
        dr = f"{ds[0]}..{ds[-1]}" if ds else "-"
        if d["depths"]:
            dep = f"{min(d['depths']):.1f}/{max(d['depths']):.1f} (n={len(d['depths'])})"
        else:
            dep = "no depth"
        units = ",".join(sorted(d["units"]))[:14]
        print(f"{c:<28}{d['n']:>7}{len(d['stns']):>6}{len(d['dates']):>7}  {dr:<24}"
              f"{dep:<18}{units}")
        mixed = {u for u in d["depth_units"] if u}
        if len(mixed) > 1:
            print(f"    depth unit codes present: {sorted(mixed)} (converted to metres)")
        if d["depth_unconverted"]:
            print(f"    {d['depth_unconverted']} depths dropped from the range above: "
                  f"unit code not recognized")
    absent = [c for c in CHARS if c not in by]
    if absent:
        print(f"  requested characteristics with no results: {', '.join(absent)}")

    # Profile depth structure for DO: how many depths a typical cast carries.
    do = by.get("Dissolved oxygen (DO)")
    if do and do["depths"]:
        prof = defaultdict(int)
        for r in rows:
            if r.get("CharacteristicName") == "Dissolved oxygen (DO)":
                if depth_m(r) is not None:
                    prof[(r.get("MonitoringLocationIdentifier"),
                          r.get("ActivityStartDate"))] += 1
        sizes = sorted(prof.values())
        if sizes:
            print(f"  DO vertical profiles: {len(sizes)} station-dates, "
                  f"depths/profile median={statistics.median(sizes):.0f} max={max(sizes)}")


def main() -> int:
    """Pull and summarize both lakes. Returns a process exit status."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out_dir", nargs="?", default=".",
                    help="directory for the raw portal CSVs (default: current dir)")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    failed = []
    for lake, bbox in LAKES.items():
        try:
            path = fetch(lake, bbox, args.out_dir)
        except (urllib.error.URLError, OSError, ValueError) as e:
            # Reported rather than raised so the second lake is still attempted, but the
            # exit status is non-zero so a partial pull cannot pass for a complete one.
            print(f"[{lake}] ERROR: {e}", flush=True)
            failed.append(lake)
            continue
        summarize(lake, path)
    if failed:
        print(f"\nFAILED: {', '.join(failed)}; downstream scripts will read a stale or "
              f"missing file for those lakes", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
