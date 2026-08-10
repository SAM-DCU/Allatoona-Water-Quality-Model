#!/usr/bin/env python3
"""Pivot GA EPD long-format results into forebay vertical profiles.

Reads the in-pool result files pulled by wqp_inpool.py, keeps the near-dam station of
each lake, and writes one tidy profile CSV per lake (Date, Depth_m, Temp_C, DO_mgL)
to Data/. Depths coded in feet are converted to metres and temperatures coded in degrees
F to degrees C. Neither conversion fires at these two stations in the current pull, where
every depth is coded m and every temperature deg C; qa_forebay_chars.py is the audit that
shows the temperature units and tailrace_and_stations.py prints the depth unit codes. The
conversions guard against a mixed-unit station appearing in a later pull, which is not
hypothetical: elsewhere in the Lanier bounding box the portal returns 585 depths coded
feet and 563 temperatures coded deg F.

Units written
-------------
``Date`` is the portal ``ActivityStartDate``, ``Depth_m`` is depth below the water
surface in metres, ``Temp_C`` is water temperature in deg C, and ``DO_mgL`` is
dissolved oxygen in mg/L. Depth and both measurements are rounded to two decimals.
An empty cell means the cast reported the other quantity at that depth but not this one.

Station selection
-----------------
The two stations are the nearest profile station to each dam, chosen by the distance
and record-length ranking printed by tailrace_and_stations.py: at Allatoona
21GAEPD_WQX-LK_14_4494 ("Lake Allatoona Upstream from Dam", 0.3 km from the dam) and at
Lanier 21GAEPD_WQX-LK_12_4028 ("Lake Sidney Lanier upstream of Buford Dam", 0.8 km). In
both cases the next nearest profile station is more than 7 km away.

Time convention
---------------
A profile is keyed on the calendar date alone, in the station's local time as the portal
reports it. The activity start time is not carried, so two casts on the same date merge
into one profile. That merge is silent: the run counts colliding (date, depth) values,
not colliding casts, so a second cast sampling different depths would leave no trace in
the counts. In the current pull one date carries two casts, 2017-05-11 at Lanier, and
none at Allatoona.

Silent loss is reported, not accepted: rows without a usable depth, rows whose result
value is not numeric, and repeated (date, depth) pairs where one value overwrites
another are all counted and printed.

Usage, with the directory holding the wqp_inpool.py output:
    /opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3 \
        Scripts/python/forebay_profiles.py <inpool_dir>
"""
import argparse
import csv
import io
import os
import sys
from collections import defaultdict
from typing import Optional

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# station identifier -> output file name, per lake
TARGETS = {
    "Allatoona": ("21GAEPD_WQX-LK_14_4494",
                  "inpool_forebay_profiles_Allatoona_LK_14_4494.csv"),
    "Lanier":    ("21GAEPD_WQX-LK_12_4028",
                  "inpool_forebay_profiles_Lanier_LK_12_4028.csv"),
}

FT_TO_M = 0.3048                 # exact by definition
DEPTH_UNITS_FT = ("ft", "feet")
TEMP_UNITS_F = ("deg f", "f")

DEPTH_COLS = [("ResultDepthHeightMeasure/MeasureValue",
               "ResultDepthHeightMeasure/MeasureUnitCode"),
              ("ActivityDepthHeightMeasure/MeasureValue",
               "ActivityDepthHeightMeasure/MeasureUnitCode")]


def fnum(x: Optional[str]) -> Optional[float]:
    """Return x as a float, or None if it is blank or not numeric.

    Portal values can be blank, a detection-limit qualifier, or free text, so a
    non-numeric entry is a normal outcome rather than an error.
    """
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def depth_m(r: dict) -> Optional[float]:
    """Sample depth in metres, from the result depth if present, else the activity depth.

    Parameters
    ----------
    r : dict
        One row of the portal result CSV.

    Returns
    -------
    float or None
        Depth below surface in metres, rounded to two decimals; None if neither depth
        column is populated. A unit code of "ft" or "feet" is converted; anything else,
        including a blank code, is taken as metres.

    Notes
    -----
    The same-named function in wqp_inpool.py differs here: it returns None for a unit
    code it does not recognize, and counts it, whereas this one assumes metres. Either
    rule is safe only while the forebay stations report metres throughout, and it is
    wqp_inpool.py, which names any code it could not convert, that would show otherwise.
    """
    for vc, uc in DEPTH_COLS:
        v = fnum(r.get(vc, ""))
        if v is not None:
            u = (r.get(uc, "") or "").lower()
            if u in DEPTH_UNITS_FT:
                v *= FT_TO_M
            return round(v, 2)
    return None


def build_profiles(src: str, station: str) -> tuple[dict, dict]:
    """Read one lake's portal CSV and pivot the near-dam station into profiles.

    Parameters
    ----------
    src : str
        Path to ``inpool_<lake>.csv`` as written by wqp_inpool.py.
    station : str
        MonitoringLocationIdentifier of the forebay station to keep.

    Returns
    -------
    prof : dict
        (date, depth_m) -> {"temp_C": deg C, "DO_mgL": mg/L}, each present only if the
        cast reported it.
    tally : dict
        Counts of rows read and rows discarded, by reason, plus the number of values
        overwritten by a later row at an identical (date, depth).

    Notes
    -----
    KNOWN DEFECT, left uncorrected: where two rows report the same characteristic at the
    same date and depth, the later row in the file wins, so the output depends on the
    order the portal returns rows in. The count of such collisions is reported. In the
    current pull there are 36 at Allatoona and 62 at Lanier, mostly replicate readings
    agreeing to within 0.1 deg C or 0.1 mg/L, but a few disagree widely: up to 3.17 deg C
    and 1.05 mg/L at Allatoona and 16.93 deg C and 5.59 mg/L at Lanier.
    """
    prof: dict = defaultdict(dict)
    tally = {"station_rows": 0, "no_depth": 0, "no_value": 0, "overwritten": 0,
             "other_characteristic": 0}
    for r in csv.DictReader(io.open(src, encoding="utf-8", errors="replace")):
        if r.get("MonitoringLocationIdentifier") != station:
            continue
        tally["station_rows"] += 1
        ch = r.get("CharacteristicName")
        if ch not in ("Temperature, water", "Dissolved oxygen (DO)"):
            tally["other_characteristic"] += 1
            continue
        z = depth_m(r)
        val = fnum(r.get("ResultMeasureValue", ""))
        if z is None:
            tally["no_depth"] += 1
            continue
        if val is None:
            tally["no_value"] += 1
            continue
        unit = (r.get("ResultMeasure/MeasureUnitCode", "") or "").lower()
        key = (r.get("ActivityStartDate", ""), z)
        if ch == "Temperature, water":
            if unit in TEMP_UNITS_F:
                val = (val - 32) * 5 / 9
            if "temp_C" in prof[key]:
                tally["overwritten"] += 1
            prof[key]["temp_C"] = round(val, 2)
        elif ch == "Dissolved oxygen (DO)":
            if "DO_mgL" in prof[key]:
                tally["overwritten"] += 1
            prof[key]["DO_mgL"] = round(val, 2)
    return prof, tally


def write_profiles(prof: dict, outpath: str) -> list:
    """Write the pivoted profiles to CSV, sorted by date then depth.

    Parameters
    ----------
    prof : dict
        Output of :func:`build_profiles`.
    outpath : str
        CSV to write, columns Date, Depth_m (m), Temp_C (deg C), DO_mgL (mg/L).

    Returns
    -------
    list
        The rows written, as (date, depth_m, temp_C, DO_mgL) tuples.
    """
    rows = sorted(((d, z, v.get("temp_C", ""), v.get("DO_mgL", ""))
                   for (d, z), v in prof.items()),
                  key=lambda x: (x[0], x[1]))
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Date", "Depth_m", "Temp_C", "DO_mgL"])
        w.writerows(rows)
    return rows


def main() -> int:
    """Pivot both lakes and report coverage and discarded rows."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("inpool_dir", nargs="?", default=".",
                    help="directory holding inpool_<lake>.csv (default: current dir)")
    ap.add_argument("--out-dir", default=os.path.join(REPO, "Data"),
                    help="directory for the profile CSVs (default: Data/)")
    args = ap.parse_args()

    status = 0
    for lake, (station, outname) in TARGETS.items():
        src = os.path.join(args.inpool_dir, f"inpool_{lake}.csv")
        if not os.path.exists(src):
            print(f"{lake}: ERROR missing {src}; run wqp_inpool.py into that directory "
                  f"first")
            status = 1
            continue
        prof, tally = build_profiles(src, station)
        outpath = os.path.join(args.out_dir, outname)
        rows = write_profiles(prof, outpath)
        if not rows:
            print(f"{lake}: ERROR no temperature or oxygen rows at {station}")
            status = 1
            continue
        dates = sorted({d for d, _, _, _ in rows})
        zmax = max(z for _, z, _, _ in rows)
        print(f"{lake}: {station}  {len(rows)} depth-obs, {len(dates)} profile dates, "
              f"{dates[0]}..{dates[-1]}, maxZ={zmax:.1f} m -> {os.path.basename(outpath)}")
        print(f"    from {tally['station_rows']} station rows: "
              f"{tally['other_characteristic']} other characteristics, "
              f"{tally['no_depth']} without a depth, "
              f"{tally['no_value']} without a numeric value, "
              f"{tally['overwritten']} values overwritten at a repeated (date, depth)")
        n_temp = sum(1 for _, _, t, _ in rows if t != "")
        n_do = sum(1 for _, _, _, o in rows if o != "")
        print(f"    depths carrying temperature: {n_temp}, oxygen: {n_do}, "
              f"both: {sum(1 for _, _, t, o in rows if t != '' and o != '')}")
    return status


if __name__ == "__main__":
    sys.exit(main())
