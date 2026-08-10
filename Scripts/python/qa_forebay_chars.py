#!/usr/bin/env python3
"""List every water-quality characteristic measured at the two forebay stations.

Queries the Water Quality Portal for the complete result set at each station and
prints one line per characteristic: result count, the share of results that carry a
sample depth, the date range, and the reported units. This is the check that resolved
the naming and unit questions behind the report's forebay constituent table, including
that GA EPD files chlorophyll under "Chlorophyll a, uncorrected for pheophytin" and
reports forebay temperature in degrees C throughout.

Nothing is written to disk; the output is an audit listing for reading. Because the
unit column is truncated to fit the table, every characteristic reported under more
than one unit is listed again in full underneath, which is the case this script exists
to catch.

Provenance
----------
Source: EPA Water Quality Portal Result service, ``resultPhysChem`` profile, queried by
station identifier with no characteristic or date filter, so the listing is the whole
record the portal holds for that station at the time of the run. Percent-with-depth is
computed over the rows returned and is floor-rounded to a whole percent.

Run with the clearwater conda env:
    /opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3 \
        Scripts/python/qa_forebay_chars.py
"""
import argparse
import csv
import io
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict

STATIONS = {
    "Allatoona": "21GAEPD_WQX-LK_14_4494",
    "Lanier":    "21GAEPD_WQX-LK_12_4028",
}
BASE = "https://www.waterqualitydata.us/data/Result/search"
TIMEOUT_S = 180

DEPTH_COLS = ("ResultDepthHeightMeasure/MeasureValue",
              "ActivityDepthHeightMeasure/MeasureValue")


def get(url: str, timeout: int = TIMEOUT_S) -> str:
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
        Response body decoded as UTF-8, replacing undecodable bytes so that a stray byte
        in a free-text portal field cannot stop the audit.

    Raises
    ------
    urllib.error.URLError
        If the request fails or times out.
    """
    req = urllib.request.Request(url, headers={"User-Agent": "qa"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def fetch_station(sid: str) -> list:
    """Return every portal result row for one station identifier.

    Parameters
    ----------
    sid : str
        MonitoringLocationIdentifier, for example ``21GAEPD_WQX-LK_14_4494``.

    Returns
    -------
    list of dict
        Rows of the ``resultPhysChem`` CSV.

    Raises
    ------
    ValueError
        If the response is not the expected result CSV, which is how the portal signals
        an error page or a changed profile rather than by an HTTP status.
    """
    url = (f"{BASE}?siteid={urllib.parse.quote(sid)}"
           "&dataProfile=resultPhysChem&mimeType=csv&zip=no")
    txt = get(url)
    if "CharacteristicName" not in txt[:2048]:
        raise ValueError(f"{sid}: response is not a resultPhysChem CSV; "
                         f"first bytes were: {txt[:200]!r}")
    return list(csv.DictReader(io.StringIO(txt)))


def report(lake: str, sid: str, rows: list) -> None:
    """Print the characteristic listing for one station.

    Parameters
    ----------
    lake : str
        Lake name, used only in the heading.
    sid : str
        MonitoringLocationIdentifier the rows were fetched for.
    rows : list of dict
        Rows of the ``resultPhysChem`` CSV from :func:`fetch_station`.

    Notes
    -----
    A result counts as carrying a depth if either depth column holds a non-blank string;
    no conversion or range check is applied, so the depth share answers whether the
    characteristic is profiled at all, not whether the depths are usable. The share is
    floor-divided to a whole percent, so 99% covers everything from 99.0 to 99.9.

    The units column is truncated to 16 characters to hold the table width, which is why
    every characteristic carrying more than one unit is repeated in full underneath.
    """
    agg = defaultdict(lambda: {"n": 0, "dates": set(), "units": set(), "depth": 0})
    for r in rows:
        c = r.get("CharacteristicName", "?")
        a = agg[c]
        a["n"] += 1
        a["dates"].add(r.get("ActivityStartDate", ""))
        u = r.get("ResultMeasure/MeasureUnitCode", "")
        if u:
            a["units"].add(u)
        for dc in DEPTH_COLS:
            if r.get(dc, "").strip():
                a["depth"] += 1
                break
    print(f"\n===== {lake} forebay {sid}: {len(rows):,} results, "
          f"{len(agg)} characteristics =====")
    print(f"{'characteristic':<46}{'n':>6}{'depth%':>7}  {'date range':<24}units")
    for c in sorted(agg, key=lambda k: -agg[k]["n"]):
        a = agg[c]
        ds = sorted(x for x in a["dates"] if x)
        dr = f"{ds[0]}..{ds[-1]}" if ds else "-"
        depthpct = f"{100 * a['depth'] // max(a['n'], 1)}%"
        units = ",".join(sorted(a["units"]))[:16]
        print(f"{c[:45]:<46}{a['n']:>6}{depthpct:>7}  {dr:<24}{units}")
    mixed = {c: sorted(a["units"]) for c, a in agg.items() if len(a["units"]) > 1}
    if mixed:
        print("  characteristics reported under more than one unit "
              "(units column above is truncated):")
        for c, us in sorted(mixed.items()):
            print(f"    {c}: {us}")


def main() -> int:
    """Audit both forebay stations. Returns a process exit status."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--station", action="append", metavar="NAME=ID",
                    help="audit an extra station instead of the two defaults; "
                         "repeatable")
    args = ap.parse_args()
    stations = dict(STATIONS)
    if args.station:
        stations = dict(s.split("=", 1) for s in args.station)

    status = 0
    for lake, sid in stations.items():
        try:
            rows = fetch_station(sid)
        except (urllib.error.URLError, OSError, ValueError) as e:
            # Reported rather than raised so the second station is still attempted, but
            # the exit status is non-zero so a partial audit cannot pass for a full one.
            print(f"{lake} ({sid}): ERROR {e}")
            status = 1
            continue
        report(lake, sid, rows)
    return status


if __name__ == "__main__":
    sys.exit(main())
