#!/usr/bin/env python3
"""Pull the Buford Dam tailrace daily records and rank the near-dam in-pool stations.

Two steps. The first downloads USGS daily values for the Chattahoochee River at Buford
Dam (02334430) and writes each one to the directory the analysis scripts read it from,
in the same column layout the dataRetrieval scripts in Scripts/ produce, so the two
routes are interchangeable. The second ranks the Water Quality Portal in-pool stations
at each lake by distance from the dam and by the size of their DO profile record,
which is how the two forebay stations used throughout the study were chosen.

Step two reads the in-pool result files pulled by wqp_inpool.py.

Units as delivered
------------------
The four daily-value series are written with the value column named for its USGS
parameter code and no unit in the file, matching the dataRetrieval layout. Units, taken
from the parameter descriptions the service returns in the file header and echoed on
every run, are: 00060 discharge in ft3/s, 00010 water temperature in deg C, 00300
dissolved oxygen in mg/L, and 00095 specific conductance in uS/cm at 25 deg C. All four
are daily means (statistic code 00003).

Time convention
---------------
A USGS daily value is a calendar-day statistic in the station's own local standard time.
The daily-value service returns no time-zone column, so the dates are used as delivered
and every downstream daily join is on that same clock.

Provisional data
----------------
The ``_cd`` qualifier column is written but never filtered downstream. A code of A means
approved for publication and P means provisional and subject to revision; composite
codes such as "A:e" add an estimation qualifier. The current qualifier counts are
printed on every run, because a rerun does not reproduce an earlier pull: the archived
files carry values that were provisional when pulled and that USGS has since revised.
The end date below is pinned for that reason, so that at least the date range of a
rerun matches the archived files.

Usage, with the directory holding the wqp_inpool.py output:
    /opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3 \
        Scripts/python/tailrace_and_stations.py <inpool_dir>
"""
import argparse
import csv
import io
import math
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ---- 1. Buford Dam tailrace daily values -> CSVs mirroring repo convention ----
# The four Buford parameters are split across two directories: discharge and DO are
# kept in Data/Lanier/ and temperature and specific conductance in Data/. The split is
# historical, from an earlier file move, and each file is written back to the location
# the analysis scripts read it from.
SITE = "02334430"
PARAMS = {"00060": ("discharge", ("Data", "Lanier")),
          "00010": ("temp", ("Data",)),
          "00300": ("do", ("Data", "Lanier")),
          "00095": ("spcond", ("Data",))}
STAT_MEAN = "00003"              # USGS statistic code for the daily mean
# Retrieval window. The end date is pinned to the date the archived files were pulled so
# that a rerun covers the same span; it is part of the file name as well.
START, END = "2005-01-01", "2026-06-29"

# ---- 2. In-pool profile stations near each dam ----
# Dam positions in decimal degrees (latitude, longitude) with the query bounding box.
# These coordinates are approximate dam-axis positions carried in this script from the
# start of the study and their source is not recorded; they are used only to order
# stations by distance, and the ordering is insensitive to an error of about a
# kilometre because the nearest station is under 1 km from the dam at both projects and
# the next nearest is more than 7 km away. Station coordinates themselves come from the
# portal Station service, and the surveyed gage and forebay coordinates used in the maps
# are in Data/station_coordinates.csv with their sources.
DAMS = {"Allatoona": (34.1631, -84.7277, "-84.78,34.09,-84.58,34.23"),
        "Lanier":    (34.1606, -84.0756, "-84.10,34.13,-83.80,34.45")}

EARTH_RADIUS_KM = 6371.0         # mean radius, adequate for a distance ranking


def get(url: str, timeout: int = 120) -> str:
    """Fetch a URL and return the decoded body.

    Parameters
    ----------
    url : str
        Request URL.
    timeout : int, optional
        Socket timeout in seconds.

    Returns
    -------
    str
        Response body decoded as UTF-8, replacing undecodable bytes.

    Raises
    ------
    urllib.error.URLError
        If the request fails or times out. Left to propagate: a partial pull that
        silently overwrites an archived file would be worse than a stopped run.
    """
    req = urllib.request.Request(url, headers={"User-Agent": "data-inventory"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def rdb_provenance(txt: str) -> list:
    """Extract the provenance lines from an RDB comment block.

    Parameters
    ----------
    txt : str
        Full RDB response, comment lines included.

    Returns
    -------
    list of str
        The retrieval timestamp, the site line, the parameter description (which is
        where the units are stated), and the data-value qualification legend. These
        lines are dropped from the written CSV, so they are echoed to stdout instead of
        being lost.
    """
    keep = []
    in_codes = False
    # The parameter-description line is found by its leading time-series identifier, which
    # is why the test below is on the digits "36" and "33": those are the ts_id prefixes
    # this site's four series happen to carry (36275, 36276, 333743, 333746). A new series
    # at another site would need another prefix, and the only symptom would be a missing
    # units line in the log.
    for line in txt.splitlines():
        if not line.startswith("#"):
            break
        body = line[1:].strip()
        if body.startswith("retrieved:") or body.startswith("USGS " + SITE):
            keep.append(body)
        elif body.startswith("Data-value qualification codes"):
            in_codes = True
        elif in_codes and body:
            keep.append(body)
        elif in_codes and not body:
            in_codes = False
        elif body.startswith(("36", "33")) and STAT_MEAN in body:
            keep.append(" ".join(body.split()))
    return keep


def fetch_daily_values(dest_root: str) -> None:
    """Download the four Buford tailrace daily-value series and write them to CSV.

    Parameters
    ----------
    dest_root : str
        Repository root to write under. Each series goes to the directory the analysis
        scripts read it from, either ``Data/`` or ``Data/Lanier/``.

    Notes
    -----
    Rows are written exactly as the service returns them, including days whose value is
    blank, so the written file is a faithful copy of the service response reshaped to
    the dataRetrieval column layout. Blank-value days and the qualifier mix are counted
    and printed rather than being removed here.
    """
    print("== Buford Dam tailrace (USGS 02334430) daily pulls ==")
    for pcode, (label, destparts) in PARAMS.items():
        url = ("https://waterservices.usgs.gov/nwis/dv/?"
               f"sites={SITE}&parameterCd={pcode}&statCd={STAT_MEAN}"
               f"&startDT={START}&endDT={END}&format=rdb")
        txt = get(url)
        lines = [ln for ln in txt.splitlines() if ln and not ln.startswith("#")]
        if len(lines) < 3:
            print(f"  {label} ({pcode}): no data")
            continue
        hdr = lines[0].split("\t")
        # Value and qualifier columns are named <ts_id>_<pcode>_<stat> by the service.
        vcol = next((i for i, h in enumerate(hdr)
                     if h.endswith(f"_{pcode}_{STAT_MEAN}")), None)
        ccol = next((i for i, h in enumerate(hdr)
                     if h.endswith(f"_{pcode}_{STAT_MEAN}_cd")), None)
        if vcol is None:
            raise ValueError(f"{label} ({pcode}): no value column in RDB header {hdr}")
        dcol = hdr.index("datetime")
        rows = []
        # Line 1 is the RDB field-width line, so data starts at line 2.
        for ln in lines[2:]:
            f = ln.split("\t")
            if len(f) <= vcol:
                continue
            rows.append((f[dcol], f[vcol], f[ccol] if ccol is not None else ""))
        dest = os.path.join(dest_root, *destparts)
        os.makedirs(dest, exist_ok=True)
        out = os.path.join(dest, f"{label}_daily_{SITE}_{START[:4]}_{END[:4]}.csv")
        with open(out, "w", newline="") as fo:
            w = csv.writer(fo, quoting=csv.QUOTE_NONNUMERIC)
            w.writerow(["Date", f"X_{pcode}_{STAT_MEAN}", f"X_{pcode}_{STAT_MEAN}_cd"])
            for d, v, c in rows:
                w.writerow([d, v, c])
        dates = [r[0] for r in rows if r[1] != ""]
        nblank = sum(1 for r in rows if r[1] == "")
        qual = defaultdict(int)
        for r in rows:
            qual[r[2]] += 1
        print(f"  {label:10s} ({pcode}): {sum(1 for r in rows if r[1] != ''):5d} daily values "
              f"{dates[0] if dates else '-'}..{dates[-1] if dates else '-'} -> "
              f"{os.path.basename(out)}")
        print(f"    {len(rows)} rows written, {nblank} with no value; "
              f"qualifier codes {dict(sorted(qual.items()))}")
        for line in rdb_provenance(txt):
            print(f"    {line}")


def haversine(a: tuple, b: tuple) -> float:
    """Great-circle distance in kilometres between two (lat, lon) pairs in degrees."""
    (la1, lo1), (la2, lo2) = a, b
    p1, p2 = math.radians(la1), math.radians(la2)
    dp = math.radians(la2 - la1)
    dl = math.radians(lo2 - lo1)
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * EARTH_RADIUS_KM * math.asin(math.sqrt(h))


def rank_stations(inpool_dir: str) -> int:
    """Rank each lake's in-pool DO profile stations by distance from the dam.

    Parameters
    ----------
    inpool_dir : str
        Directory holding ``inpool_<lake>.csv`` from wqp_inpool.py.

    Returns
    -------
    int
        0 on success, 1 if any lake's input file is missing.

    Notes
    -----
    Prints one line per station: distance in km, identifier, number of DO results, the
    date range of those results, the deepest reported sample, and the station name. The
    maximum depth is taken from whichever depth column is populated and is not unit
    converted, so a station reporting depth in feet would show a misleadingly large
    figure; every station in this pull reports metres, which is checked and reported.
    """
    status = 0
    for lake, (dlat, dlon, bbox) in DAMS.items():
        inpool = os.path.join(inpool_dir, f"inpool_{lake}.csv")
        if not os.path.exists(inpool):
            print(f"\n== {lake}: ERROR missing {inpool}; run wqp_inpool.py into that "
                  f"directory first")
            status = 1
            continue
        # Station coordinates and names come from a second portal service: the Result
        # profile carries neither, so the distance ranking cannot be built from the
        # in-pool file alone.
        surl = ("https://www.waterqualitydata.us/data/Station/search?"
                f"bBox={urllib.parse.quote(bbox)}"
                f"&siteType={urllib.parse.quote('Lake, Reservoir, Impoundment')}"
                "&mimeType=csv&zip=no")
        smeta = {}
        for r in csv.DictReader(io.StringIO(get(surl, 180))):
            smeta[r["MonitoringLocationIdentifier"]] = (
                r.get("MonitoringLocationName", ""),
                r.get("LatitudeMeasure", ""), r.get("LongitudeMeasure", ""))
        # DO profile counts per station from the in-pool file pulled by wqp_inpool.py
        stn = defaultdict(lambda: {"n": 0, "dates": set(), "dmax": 0.0})
        depth_units = set()
        for r in csv.DictReader(io.open(inpool, encoding="utf-8", errors="replace")):
            if r.get("CharacteristicName") != "Dissolved oxygen (DO)":
                continue
            s = r["MonitoringLocationIdentifier"]
            stn[s]["n"] += 1
            stn[s]["dates"].add(r.get("ActivityStartDate", ""))
            for vc, uc in (("ResultDepthHeightMeasure/MeasureValue",
                            "ResultDepthHeightMeasure/MeasureUnitCode"),
                           ("ActivityDepthHeightMeasure/MeasureValue",
                            "ActivityDepthHeightMeasure/MeasureUnitCode")):
                try:
                    v = float(r.get(vc, ""))
                except (TypeError, ValueError):
                    continue
                stn[s]["dmax"] = max(stn[s]["dmax"], v)
                depth_units.add((r.get(uc, "") or "").strip())
                break
        print(f"\n== {lake}: in-pool DO-profile stations (nearest dam first) ==")
        print(f"  depth unit codes in this pull: {sorted(depth_units)}")
        rows = []
        nocoord = []
        for s, d in stn.items():
            name, lat, lon = smeta.get(s, ("?", "", ""))
            try:
                dist = haversine((dlat, dlon), (float(lat), float(lon)))
            except (TypeError, ValueError):
                dist = float("nan")
                nocoord.append(s)
            ds = sorted(x for x in d["dates"] if x)
            rows.append((dist, s, name, d["n"], ds[0] if ds else "-",
                         ds[-1] if ds else "-", d["dmax"]))
        # NaN distances sort last rather than being dropped, and are named below.
        for dist, s, name, n, d0, d1, dmax in sorted(
                rows, key=lambda x: (x[0] if x[0] == x[0] else 9e9)):
            print(f"  {dist:5.1f} km  {s:22s} n={n:5d} {d0}..{d1} maxZ={dmax:.0f}m  "
                  f"{name[:40]}")
        if nocoord:
            print(f"  no coordinates in the Station service, distance unknown: "
                  f"{', '.join(nocoord)}")
    return status


def main() -> int:
    """Pull the daily values, then rank the in-pool stations."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("inpool_dir", nargs="?", default=".",
                    help="directory holding inpool_<lake>.csv (default: current dir)")
    ap.add_argument("--dest-root", default=REPO,
                    help="repository root to write the daily-value CSVs under "
                         "(default: the repository this script lives in)")
    args = ap.parse_args()
    try:
        fetch_daily_values(args.dest_root)
    except (urllib.error.URLError, OSError, ValueError) as e:
        print(f"ERROR pulling USGS daily values: {e}")
        return 1
    return rank_stations(args.inpool_dir)


if __name__ == "__main__":
    sys.exit(main())
