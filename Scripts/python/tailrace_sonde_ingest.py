#!/usr/bin/env python3
"""Ingest and QA the USACE-SAM Allatoona tailrace water-quality sonde spreadsheets.

TJ (USACE Mobile District) supplied quarterly Hydrolab DataSonde downloads for the
Allatoona Dam tailrace (below-dam continuous monitoring station co-located with USGS
02394000). This script parses the two spreadsheet layouts present in
``Data/Allatoona_from_TJ_2026-07-07/``, harmonizes units, applies gross-range and
instrument-failure quality control, and writes a tidy record to
``Data/allatoona_tailrace_sonde_2011_2019.csv``.

Two layouts are handled:
  A. Continuous sub-hourly logger export with a two-row column header, parameter names
     then units, beginning with a cell literally equal to ``Date``. A DataSonde metadata
     block precedes the header in the three 2011 files and the 2014 file, follows it in
     2012 Q3 (which then repeats the header), and is absent from the rest.
     Parameter columns: Temp (C), pH, ORP, SpCond (mS/cm), LDO (% sat), LDO (mg/L).
  B. Sparse grab-reading export (2019 file) with a two-row header on rows 0-1 and no
     metadata block. Temperature is reported in degrees F; DO is column ``DO`` (mg/L).

Units as delivered and as written
---------------------------------
The exports carry temperature in deg C, deg F, and deg K side by side, pH in standard
units, ORP in mV, specific conductance in both mS/cm and uS/cm, and oxygen as percent
saturation and mg/L. The tidy CSV writes T_C in deg C, pH in standard units, ORP_mV in
mV, SpCond_mScm in mS/cm as delivered, DO_pctsat in percent, and DO_mgL and DO_mgL_qa
in mg/L. Only the mS/cm conductance column is carried, so no conductance conversion is
applied.

Time convention
---------------
The logger writes a local clock time and the export names no time zone, so timestamps
are carried through as recorded and every daily aggregate downstream is on that same
logger clock. Whether the logger was held on standard time year round is not documented
in the source. The logging interval is 60 min in the first two 2011 quarters and 30 min
from 2011 Q3 onward.

File coverage note: the file named ``2014- 4th Qtr WQ Data- Allatoona.xlsx`` contains
2013-10-01 through 2013-12-31, so the continuous record spans 2011 through 2013 with no
2014 data. The 2019 grab file runs past its title as well, to 2020-02-18.

Provenance and location were established by cross-validating the 2012 third-quarter
sonde temperature against the published USGS 02394000 tailrace daily temperature
(mean absolute difference 0.18 C, Pearson r = 0.995), confirming the same below-dam
station.

Run with the clearwater conda env:
    /opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3 \
        Scripts/python/tailrace_sonde_ingest.py
"""
import argparse
import glob
import os
import warnings
import datetime as dt
from typing import Optional

import pandas as pd
import numpy as np

# pandas falls back to dateutil, one element at a time, wherever a date column holds
# mixed types; that is expected here because the metadata block is text. Suppress only
# that message so real warnings still reach the terminal.
warnings.filterwarnings("ignore", message="Could not infer format",
                        category=UserWarning)

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC_DIR = os.path.join(REPO, "Data", "Allatoona_from_TJ_2026-07-07")
OUT_CSV = os.path.join(REPO, "Data", "allatoona_tailrace_sonde_2011_2019.csv")

# Gross-range physical limits for the Allatoona tailrace. Unlike Buford, the Allatoona
# release warms through the season, so the temperature ceiling is the physical limit for
# a Georgia tailwater rather than a site-specific one; it removes nothing in this record.
T_MIN, T_MAX = 0.0, 35.0         # deg C
DO_MIN, DO_MAX = 0.0, 15.0       # mg/L; > ~14.6 exceeds freshwater saturation at 0 C
DO_SUPERSAT_PCT = 140.0          # % saturation above which a reading is a sensor error
PH_PHYSICAL_MIN = 4.0            # circumneutral soft water; pH < 4 indicates sensor failure

# 2011 deployment is instrument-compromised: the pH sensor reported exactly 0.00 on 687
# readings and below 4 on 3,613 of the year's 6,257, and the DO sensor sat near zero for
# most of the year, with monthly medians of 0.15, 0.14, 0.10, and 0.09 mg/L in January,
# February, May, and June against 10.45, 9.39, 7.10, and 5.10 mg/L in the same months of
# 2012-2013. Temperature that year remains usable, but 2011 DO is excluded from
# quantitative DO statistics.
DO_INSTRUMENT_FAIL_YEARS = {2011}

# Median spacing at or below this many hours marks a continuous logger deployment; above
# it, the file is intermittent grab readings. The continuous files log at 30 to 60 min
# and the 2019 grab file at a median of six days, so the classification is not sensitive
# to where in that gap the cut falls.
CONTINUOUS_MAX_GAP_H = 6.0

# The DataSonde export writes a marker character in the column immediately to the right
# of a parameter value. These markers are read only to report how many oxygen readings
# carry one; they do not enter the QA flags, because the export supplies no legend and
# the meaning of "*" against "#" has not been established from Hydrolab documentation.
INSTRUMENT_MARKERS = ("*", "#")


def _combine_datetime(date_cell, time_cell) -> pd.Timestamp:
    """Combine a date cell and a time-of-day cell into a timestamp.

    A missing or unparsable time of day yields midnight of the date, which is what the
    2019 grab file needs on the days it records no reading time.
    """
    d = pd.to_datetime(date_cell, errors="coerce")
    if pd.isna(d):
        return pd.NaT
    if isinstance(time_cell, dt.time):
        return pd.Timestamp(d.year, d.month, d.day, time_cell.hour,
                            time_cell.minute, time_cell.second)
    t = pd.to_datetime(time_cell, errors="coerce")
    if pd.isna(t):
        return pd.Timestamp(d.year, d.month, d.day)
    return pd.Timestamp(d.year, d.month, d.day, t.hour, t.minute, t.second)


def _unit(v) -> str:
    """Return a spreadsheet cell as a stripped string, with missing cells as ''."""
    return str(v).strip() if pd.notna(v) else ""


def _build_column_map(param_row: list, unit_row: list) -> dict:
    """Map physical quantities to column indices using the two header rows.

    Parameters
    ----------
    param_row : list
        Header row cells naming each parameter (Temp, pH, ORP, LDO, ...).
    unit_row : list
        The row beneath it, naming each unit (deg C, mg/l, ...).

    Returns
    -------
    dict
        Keys among ``T_C``, ``T_F``, ``pH``, ``ORP_mV``, ``SpCond_mScm``,
        ``DO_pctsat``, ``DO_mgL``, each mapped to a zero-based column index.

    Notes
    -----
    Column positions are stable across the continuous files, but resolving by label also
    covers the 2019 grab file, whose oxygen columns are named ``DO``/``% O2`` rather
    than ``LDO``/``LDO%``. Both temperature and conductance appear more than once with
    different units, so the unit cell, not the parameter name, decides which column is
    taken: the deg K temperature and the uS/cm conductance fall through unmapped. Where
    a parameter name did repeat with the same unit the later column would win, which is
    why the resolved map is printed for every file.

    The degree symbol in the unit cells arrives as mojibake from the original encoding,
    so only the trailing letter of the unit is tested.
    """
    cmap = {}
    for c in range(len(param_row)):
        p = _unit(param_row[c])
        u = _unit(unit_row[c])
        if p == "Temp" and u.endswith("C"):
            cmap["T_C"] = c
        elif p == "Temp" and u.endswith("F"):
            cmap["T_F"] = c
        elif p == "pH":
            cmap["pH"] = c
        elif p == "ORP":
            cmap["ORP_mV"] = c
        elif p == "SpCond" and "mS" in u:
            cmap["SpCond_mScm"] = c
        elif p in ("LDO%", "% O2"):
            cmap["DO_pctsat"] = c
        elif p in ("LDO", "DO"):
            cmap["DO_mgL"] = c
    return cmap


def count_marked_oxygen(data: pd.DataFrame, cmap: dict) -> int:
    """Count oxygen readings the export marks with an instrument character.

    Parameters
    ----------
    data : pandas.DataFrame
        Data rows of one file, positionally indexed as in the sheet.
    cmap : dict
        Column map from :func:`_build_column_map`.

    Returns
    -------
    int
        Number of rows with a numeric DO value whose adjacent marker cell holds one of
        INSTRUMENT_MARKERS. Zero if the file has no DO column or no marker column.
    """
    c = cmap.get("DO_mgL")
    if c is None or c + 1 >= data.shape[1]:
        return 0
    do = pd.to_numeric(data.iloc[:, c], errors="coerce")
    mark = data.iloc[:, c + 1].astype(str).str.strip().isin(INSTRUMENT_MARKERS)
    return int((do.notna() & mark).sum())


def parse_file(path: str) -> Optional[pd.DataFrame]:
    """Parse one sonde spreadsheet into harmonized columns and units.

    Parameters
    ----------
    path : str
        Path to a quarterly ``.xls`` or ``.xlsx`` export. Sheet1 is read.

    Returns
    -------
    pandas.DataFrame or None
        Columns ``datetime``, ``T_C`` (deg C), ``pH``, ``ORP_mV`` (mV),
        ``SpCond_mScm`` (mS/cm), ``DO_pctsat`` (percent), ``DO_mgL`` (mg/L),
        ``sampling``, and ``source_file``. Rows carrying no measurement at all are
        dropped. None if no header row is found, which is reported by the caller.

    Notes
    -----
    Data rows are those whose first cell parses as a calendar date, which skips the
    DataSonde metadata block wherever it is positioned relative to the header and also
    steps over the second header block some files repeat below the metadata. Where a
    file carries more than one header row the first is used; in this data set the
    repeated headers are identical to the first, and that is checked on every run.
    """
    fname = os.path.basename(path)
    df = pd.read_excel(path, sheet_name="Sheet1", header=None)
    date_rows = [r for r in range(min(40, len(df)))
                 if str(df.iloc[r, 0]).strip() == "Date"]
    if not date_rows:
        return None
    hdr = date_rows[0]
    cmap = _build_column_map(df.iloc[hdr].tolist(), df.iloc[hdr + 1].tolist())
    # Compared as text because a direct cell comparison would report every empty
    # spacer column as a difference: NaN is not equal to NaN.
    header_text = [_unit(v) for v in df.iloc[hdr].tolist()]
    for extra in date_rows[1:]:
        if [_unit(v) for v in df.iloc[extra].tolist()] != header_text:
            print(f"  {fname}: WARNING repeated header at row {extra} differs from the "
                  f"header at row {hdr}; columns are mapped from row {hdr}")

    body = df.iloc[hdr + 1:].reset_index(drop=True)
    col0_dt = pd.to_datetime(body.iloc[:, 0], errors="coerce")
    data = body[col0_dt.notna()].reset_index(drop=True)
    if not len(data):
        return None

    ts = [_combine_datetime(data.iloc[i, 0], data.iloc[i, 1]) for i in range(len(data))]

    def col(name: str) -> pd.Series:
        """Numeric values of a mapped column; all-NaN if the file lacks it."""
        if name in cmap:
            return pd.to_numeric(data.iloc[:, cmap[name]], errors="coerce")
        return pd.Series(np.nan, index=data.index)

    tC = col("T_C")
    if "T_C" not in cmap and "T_F" in cmap:   # no Celsius column at all
        tC = (col("T_F") - 32.0) * 5.0 / 9.0
    elif "T_F" in cmap:
        # Both columns exist but the Celsius one can be empty, which is how the 2019
        # grab file reports: a labelled but blank deg C column and readings in deg F.
        tF = col("T_F")
        tC = tC.where(tC.notna(), (tF - 32.0) * 5.0 / 9.0)

    out = pd.DataFrame({
        "datetime": ts,
        "T_C": tC,
        "pH": col("pH"),
        "ORP_mV": col("ORP_mV"),
        "SpCond_mScm": col("SpCond_mScm"),
        "DO_pctsat": col("DO_pctsat"),
        "DO_mgL": col("DO_mgL"),
    })
    meas = ["T_C", "pH", "ORP_mV", "SpCond_mScm", "DO_pctsat", "DO_mgL"]
    kept = out[meas].notna().any(axis=1)
    out = out[kept].copy()
    # The sampling label is derived rather than read from the file, because nothing in the
    # export says whether a deployment logged continuously or was visited for grab
    # readings, and downstream statistics weight the two differently.
    tsv = pd.to_datetime(out["datetime"]).sort_values()
    med_gap_h = tsv.diff().dt.total_seconds().median() / 3600.0 if len(tsv) > 2 else 24.0
    out["sampling"] = "continuous" if med_gap_h <= CONTINUOUS_MAX_GAP_H else "grab"
    out["source_file"] = fname
    print(f"  {fname}: rows {len(body)} -> {len(data)} dated -> {len(out)} with a "
          f"measurement; median gap {med_gap_h:.2f} h ({out['sampling'].iloc[0]}); "
          f"DO readings carrying an instrument marker: "
          f"{count_marked_oxygen(data, cmap)}; columns {sorted(cmap)}")
    return out


def apply_qa(rec: pd.DataFrame) -> pd.DataFrame:
    """Apply gross-range and instrument-failure quality control.

    Parameters
    ----------
    rec : pandas.DataFrame
        Concatenated parsed files, one row per observation.

    Returns
    -------
    pandas.DataFrame
        Time-sorted, with ``year``, ``month``, ``do_flag``, ``DO_mgL_qa`` (mg/L), and
        ``ph_flag`` added. Temperature outside [T_MIN, T_MAX] is set to NaN in place.

    Notes
    -----
    ``do_flag`` takes one of three values, assigned in this order so that later
    assignments win: "range" for a reading at or below DO_MIN, above DO_MAX, or above
    DO_SUPERSAT_PCT percent saturation; "instrument_fail_2011" for any reading in a year
    listed in DO_INSTRUMENT_FAIL_YEARS; "good" otherwise. ``DO_mgL_qa`` carries the
    value only where the flag is "good". There is no "missing" flag, so a row that
    reports another quantity but no oxygen is left "good" with an empty ``DO_mgL_qa``:
    count DO_mgL_qa, not the flag, to get the number of usable oxygen readings. One row
    in the current record takes that path.

    The luminescent (optical) DO sensor is independent of the pH probe, so pH-sensor
    failures do not invalidate DO; the pH flag is reported separately and removes
    nothing. A reading of exactly 0.000 mg/L is treated as an instrument artifact rather
    than a measurement, which is why the lower bound is exclusive.
    """
    rec = rec.dropna(subset=["datetime"]).copy()
    rec["datetime"] = pd.to_datetime(rec["datetime"])
    rec = rec.sort_values("datetime").reset_index(drop=True)
    rec["year"] = rec["datetime"].dt.year
    rec["month"] = rec["datetime"].dt.month

    # Temperature QA. The count is reported because this edit is silent in the CSV:
    # a screened temperature is written as an empty cell with no accompanying flag.
    t_out = ((rec["T_C"] < T_MIN) | (rec["T_C"] > T_MAX)).sum()
    rec.loc[(rec["T_C"] < T_MIN) | (rec["T_C"] > T_MAX), "T_C"] = np.nan
    print(f"Temperature outside {T_MIN:.0f} to {T_MAX:.0f} deg C set to missing: {t_out}")

    do_flag = pd.Series("good", index=rec.index)
    do_flag[(rec["DO_mgL"] <= DO_MIN) | (rec["DO_mgL"] > DO_MAX)] = "range"
    do_flag[rec["DO_pctsat"] > DO_SUPERSAT_PCT] = "range"
    do_flag[rec["year"].isin(DO_INSTRUMENT_FAIL_YEARS)] = "instrument_fail_2011"
    rec["do_flag"] = do_flag
    rec["DO_mgL_qa"] = rec["DO_mgL"].where(do_flag == "good", np.nan)

    # Separate pH reliability flag: the pH probe drifted to non-physical values
    # (readings of 0 to about 4) during several deployments. pH is not used
    # quantitatively; this flag documents the affected records.
    rec["ph_flag"] = np.where(rec["pH"] < PH_PHYSICAL_MIN, "suspect", "ok")
    rec.loc[rec["pH"].isna(), "ph_flag"] = "missing"
    return rec


def summarize(rec: pd.DataFrame) -> None:
    """Print coverage, QA disposition, and the seasonal and summer DO statistics."""
    print("=" * 78)
    print("PARSED RECORD SUMMARY")
    print("=" * 78)
    print(f"Total records: {len(rec):,}")
    print(f"Date span: {rec['datetime'].min()} -> {rec['datetime'].max()}")
    # Repeated timestamps are carried through as separate observations rather than merged.
    # All 37 rows involved are in 2011: 28 of them are the 14 hours of 2011-04-01 that the
    # first and second quarter files both export, and the other 9 are repeats inside the
    # third quarter file on 2011-08-27, where three timestamps appear twice and one
    # appears three times.
    ndup = int(rec["datetime"].duplicated().sum())
    if ndup:
        dup_years = rec[rec["datetime"].duplicated(keep=False)]["year"].value_counts()
        print(f"Repeated timestamps: {ndup} rows share a timestamp with an earlier row "
              f"(rows involved by year: {dup_years.to_dict()})")
    print("\nRecords and DO-QA disposition by year:")
    g = rec.groupby("year").agg(
        n=("datetime", "size"),
        sampling=("sampling", lambda s: s.iloc[0]),
        T_med=("T_C", "median"),
        DO_good=("DO_mgL_qa", lambda s: s.notna().sum()),
        DO_med_qa=("DO_mgL_qa", "median"),
    )
    print(g.to_string())
    print(f"\nDO-flag counts:\n{rec['do_flag'].value_counts().to_string()}")

    good = rec[rec["do_flag"] == "good"].copy()
    print("\n" + "=" * 78)
    print("QUALITY-CONTROLLED CONTINUOUS DO (good flag only)")
    print("=" * 78)
    print(f"Good DO records: {len(good):,}  "
          f"({good['datetime'].min()} -> {good['datetime'].max()})")
    seas = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
            6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}
    good["season"] = good["month"].map(seas)
    print("\nSeasonal DO (mg/L), good continuous data:")
    s = good.groupby("season")["DO_mgL_qa"].agg(["size", "median", "min", "max"])
    print(s.reindex(["DJF", "MAM", "JJA", "SON"]).to_string())

    summer = good[good["month"].isin([6, 7, 8, 9])]
    print("\nSummer (Jun-Sep) tailrace DO, good continuous data:")
    print(f"  n = {len(summer):,} hourly obs")
    print(f"  median = {summer['DO_mgL_qa'].median():.2f} mg/L")
    print(f"  mean   = {summer['DO_mgL_qa'].mean():.2f} mg/L")
    print(f"  min    = {summer['DO_mgL_qa'].min():.2f} mg/L")
    for thr in (2, 4, 5):
        pct = 100.0 * (summer["DO_mgL_qa"] < thr).mean()
        print(f"  % summer hours < {thr} mg/L: {pct:.1f}%")

    # Diel amplitude, as the summer daily maximum minus the daily minimum. The 12-reading
    # floor keeps out the three grab readings that fall in summer, but the continuous
    # files log at 30 min, so it still admits a day covering as little as six hours, whose
    # amplitude understates the true diel range.
    summer_daily = summer.set_index("datetime")["DO_mgL_qa"].resample("D")
    diel = (summer_daily.max() - summer_daily.min()).dropna()
    diel = diel[summer.set_index("datetime")["DO_mgL_qa"].resample("D").count() >= 12]
    print(f"\n  Median summer diel DO amplitude (days with >=12 hourly obs): "
          f"{diel.median():.2f} mg/L (n={len(diel)} days)")

    # Temperature is summarized over every year, 2011 included: the instrument-failure
    # flag removes 2011 oxygen only, and the 2011 temperature record is sound.
    print("\n" + "=" * 78)
    print("TAILRACE TEMPERATURE (all valid records, all years)")
    print("=" * 78)
    tvalid = rec.dropna(subset=["T_C"])
    tsum = tvalid[tvalid["month"].isin([6, 7, 8, 9])]
    print(f"  Summer (Jun-Sep) median T = {tsum['T_C'].median():.1f} C "
          f"(range {tvalid['T_C'].min():.1f}-{tvalid['T_C'].max():.1f} C over full record)")


def main() -> None:
    """Parse every workbook, apply QA, write the tidy CSV, and print the summary."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--src", default=SRC_DIR,
                    help="directory of sonde workbooks "
                         "(default: Data/Allatoona_from_TJ_2026-07-07)")
    ap.add_argument("--out", default=OUT_CSV,
                    help="tidy CSV to write (default: "
                         "Data/allatoona_tailrace_sonde_2011_2019.csv)")
    args = ap.parse_args()

    paths = sorted(glob.glob(os.path.join(args.src, "*.xls*")))
    if not paths:
        raise SystemExit(f"no .xls/.xlsx workbooks found in {args.src}")
    print(f"Reading {len(paths)} workbooks from {args.src}")
    frames = []
    skipped = []
    for p in paths:
        rec = parse_file(p)
        if rec is not None and len(rec):
            frames.append(rec)
        else:
            skipped.append(os.path.basename(p))
    if skipped:
        print(f"Skipped, no 'Date' header row or no dated rows: {', '.join(skipped)}")
    if not frames:
        raise SystemExit(f"no data rows parsed from {args.src}")
    allrec = pd.concat(frames, ignore_index=True)
    allrec = apply_qa(allrec)
    cols = ["datetime", "year", "month", "sampling", "T_C", "pH", "ph_flag", "ORP_mV",
            "SpCond_mScm", "DO_pctsat", "DO_mgL", "do_flag", "DO_mgL_qa", "source_file"]
    allrec[cols].to_csv(args.out, index=False)
    print(f"Wrote {args.out} ({len(allrec):,} records)\n")
    summarize(allrec)


if __name__ == "__main__":
    main()
