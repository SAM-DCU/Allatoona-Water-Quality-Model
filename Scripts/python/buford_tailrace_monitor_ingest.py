#!/usr/bin/env python3
"""Ingest and QA the Buford Dam tailrace water-quality monitor record, 2002-2008.

USACE Mobile District supplied ``Data/Lanier_2002-2010.xls``, a workbook with one sheet
per year holding sub-hourly dissolved oxygen, temperature, pH, and conductivity from the
Corps continuous water-quality monitoring station on the Chattahoochee River below
Buford Dam. The Buford Water Control Manual (ACF Appendix B, p. 5-02) records that the
Corps operated this station from 1981 to 2008, measuring dissolved oxygen, temperature,
pH, and conductivity, with the data held by the Mobile District Planning Division,
Inland Environment (PD-EI) Office. The workbook is therefore a subset of a longer record
that remains to be requested.

This script harmonizes the three spreadsheet layouts present in the workbook, applies
gross-range, supersaturation, and flatline quality control, and writes a tidy record to
``Data/buford_tailrace_monitor_2002_2008.csv``.

Three layouts are handled:
  A. 2002: hourly export with columns ``Recorded Date``, ``D.O.``, ``Temp`` (deg F),
     ``pH``, ``COND.``; no separate time column.
  B. 2004-2005: ``Date`` and ``Time`` columns with both ``Temp deg C`` and
     ``Temp deg F``, and dissolved oxygen labelled ``DO mg/L``.
  C. 2006-2008: ``Date`` and ``Time`` columns with a separate unit row, oxygen labelled
     ``DO`` or ``LDO``, and percent saturation in the following column.

Despite the file name, the workbook carries sheets for 2002 and 2004 through 2008 only;
there is no 2003, 2009, or 2010 sheet.

Units as delivered and as written
---------------------------------
The workbook reports temperature in deg C, deg F, or both; dissolved oxygen in mg/L;
percent saturation in percent; conductivity in uS/cm except on the 2004 sheet, which
reports mS/cm. The tidy CSV writes T_C in deg C, DO_mgL and DO_mgL_qa in mg/L,
DO_pctsat and DO_pctsat_est in percent, SpCond_uScm in uS/cm, and pH in standard units.

Time convention
---------------
The workbook records a local clock time from the Corps logger and names no time zone,
so timestamps are carried through exactly as recorded and every daily or seasonal
aggregate downstream is on that same logger clock. Whether the logger was held on
standard time year round is not documented in the source and has not been established.
The sampling interval changes within the record: 60 min in 2002, 30 min in 2004, and
15 min in 2005 through 2008.

Station identity was established by cross-validating the daily mean monitor temperature
against the published USGS 02334430 tailrace daily temperature: mean absolute difference
0.44 C over 681 paired days, and by year 0.50 C in 2005 (Pearson r = 0.991), 0.53 C in
2006 (r = 0.970), 0.39 C in 2007 (r = 0.952), and 0.34 C in 2008 (r = 0.950), confirming
the same below-dam station. That comparison is owned by tailwater_do_statistics.py and
its current values are in analysis/tailwater_do_statistics.txt.

Two defects in the source workbook require handling:
  1. In the 2008 sheet the Celsius and Fahrenheit columns exchange positions without a
     change of header, three times, at the deployment breaks on 2008-03-12, 2008-06-10,
     and 2008-08-05. Temperature is therefore resolved per row from the internal
     consistency of the two columns rather than from the header.
  2. The oxygen sensor reports exactly 0.000 mg/L continuously from 2008-06-27 to
     2008-07-24 and returns isolated values above the freshwater saturation ceiling in
     several months. Both are removed by the quality control below.

Every row the parser discards is counted and reported to stdout, including the
deployment notes the logger writes inline in the data block ("Power loss from ...",
"DEAD BATTERIES RESULTED IN DATA LOSS FROM 1/23 - 2/7 2006"), so no data or gap
annotation leaves the ingest unaccounted for.

Run with the clearwater conda env:
    /opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3 \
        Scripts/python/buford_tailrace_monitor_ingest.py
"""
import argparse
import os
import warnings
import datetime as dt
from typing import Optional

import numpy as np
import pandas as pd

# pandas falls back to dateutil, one element at a time, wherever a date column holds
# mixed types; that is expected here because the annotation rows are text. Suppress only
# that message so real warnings still reach the terminal.
warnings.filterwarnings("ignore", message="Could not infer format",
                        category=UserWarning)

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC_XLS = os.path.join(REPO, "Data", "Lanier_2002-2010.xls")
OUT_CSV = os.path.join(REPO, "Data", "buford_tailrace_monitor_2002_2008.csv")

# Gross-range limits for the Buford tailrace. The upper temperature limit is
# site-specific rather than merely physical: the Buford release is drawn from the deep
# hypolimnion and never exceeds 16.4 C anywhere in the 21-year USGS 02334430 record
# (2005-2026), so a 17.0 C ceiling flags sensor excursions without touching valid data.
# The screen removes 40 records. The 11 of them that fall on days the USGS record covers
# exceed the concurrent USGS daily mean by 4.2 to 13.5 C, which identifies them as
# instrument error rather than real warm releases; the other 29 are in May 2002, before
# the USGS record begins. A further 2,563 records, all but 9 of them in the
# failed-instrument window of June and July 2008, satisfy neither Celsius/Fahrenheit
# assignment and carry no temperature at all.
T_MIN, T_MAX = 0.0, 17.0         # deg C
DO_MIN, DO_MAX = 0.0, 15.0       # mg/L; > ~14.6 exceeds freshwater saturation at 0 C
DO_SUPERSAT_PCT = 140.0          # % saturation above which a reading is a sensor error
PH_PHYSICAL_MIN = 4.0            # circumneutral soft water; pH < 4 indicates sensor failure

# Flatline screen, in samples rather than in time. Because the logging interval changes
# within the record, 24 consecutive identical values span 24 h in 2002, 12 h in 2004,
# and 6 h from 2005 on. In the archived record the screen fires only inside the known
# failed-instrument window of June and July 2008 (2,554 records), where the sensor
# reports exactly 0.000 mg/L, so the varying time span costs nothing here; a record with
# a longer stable period at 15 min spacing would need the threshold set per interval.
FLATLINE_N = 24

# Celsius/Fahrenheit agreement tolerance, deg F. The two temperature columns are rounded,
# at most to hundredths in Celsius and thousandths in Fahrenheit, so a consistent pair
# differs by at most 0.06 F anywhere in this workbook, while the wrong orientation misses
# by at least 106 F over this record's temperature range. 0.6 F sits in that gap, so the
# classification does not depend on the exact value chosen.
CF_TOLERANCE_F = 0.6

# Minimum length, in characters, of a first-column string worth echoing as a logger
# annotation. The rows without a parseable date carry the deployment log (power loss,
# late probe start, dead batteries), which is reported rather than dropped in silence,
# mixed with blank separator rows; the threshold keeps the blanks out of the report.
ANNOTATION_MIN_CHARS = 3

# Sensor-decay episodes, as half-open [start, end) intervals in the record's own clock.
#
# 2004-06-17 to 06-24: the deployment opens at 6.82 mg/L and falls monotonically to
# 0.01 mg/L over seven days, then returns to 11.51 mg/L (97.7 percent saturation) in a
# single 30 min step at 14:00 on 06-24. Temperature holds at 8.1 to 9.4 deg C, specific
# conductance at 33 uS/cm, and pH at 6.6 to 6.7 across both the decay and the recovery,
# so no property of the water changes when the oxygen reading does. A hypolimnetic
# release cannot fall to 0.01 mg/L and recover to saturation within half an hour, and
# the concurrent forebay casts show no such source water. The pattern is a fouled
# membrane or a depleted electrolyte, restored at a service visit.
#
# These readings passed every gross-range test because they are individually physical,
# which is why the screen is a named window rather than a threshold. The window supplies
# the record minimum, so its removal changes the reported minimum DO.
SENSOR_DECAY_WINDOWS = [("2004-06-17 00:00", "2004-06-24 14:00")]


def do_saturation(t_c: pd.Series | np.ndarray) -> np.ndarray:
    """Freshwater dissolved-oxygen saturation at 1 atm.

    Parameters
    ----------
    t_c : pandas.Series or numpy.ndarray
        Water temperature, deg C.

    Returns
    -------
    numpy.ndarray
        Saturation dissolved-oxygen concentration, mg/L.

    Notes
    -----
    Benson and Krause as given in APHA Standard Methods, evaluated at 1 atm and zero
    salinity. Used only to screen non-physical supersaturation, because the 2002-2005
    sheets do not report percent saturation.

    No barometric or elevation correction is applied. The gage datum at the co-located
    USGS station 02334430 is 912.01 ft NAVD88, about 278 m (USGS NWIS site service),
    where mean barometric pressure is roughly 3 percent below 1 atm. Saturation scales
    with pressure, so evaluating at 1 atm overstates the ceiling and understates percent
    saturation by about the same 3 percent; the fixed 140 percent screen below is
    therefore conservative, flagging slightly fewer records than an elevation-corrected
    screen would.
    """
    tk = np.asarray(t_c, dtype=float) + 273.15
    with np.errstate(invalid="ignore", divide="ignore"):
        ln_c = (-139.34411 + 1.575701e5 / tk - 6.642308e7 / tk ** 2
                + 1.243800e10 / tk ** 3 - 8.621949e11 / tk ** 4)
    return np.exp(ln_c)


def _combine_datetime(date_cell, time_cell) -> pd.Timestamp:
    """Combine a date cell and an optional time-of-day cell into a timestamp.

    A missing or unparsable time of day yields midnight of the date, which is what the
    2002 sheet needs: it carries the hour inside the date cell and has no time column.
    """
    d = pd.to_datetime(date_cell, errors="coerce")
    if pd.isna(d):
        return pd.NaT
    if time_cell is None:
        return d
    if isinstance(time_cell, dt.time):
        return pd.Timestamp(d.year, d.month, d.day, time_cell.hour,
                            time_cell.minute, time_cell.second)
    t = pd.to_datetime(time_cell, errors="coerce")
    if pd.isna(t):
        return pd.Timestamp(d.year, d.month, d.day)
    return pd.Timestamp(d.year, d.month, d.day, t.hour, t.minute, t.second)


def _label(v) -> str:
    """Return a spreadsheet cell as a stripped string, with missing cells as ''."""
    return str(v).strip() if pd.notna(v) else ""


def _resolve_temperature(col_a: pd.Series, col_b: pd.Series) -> pd.Series:
    """Resolve temperature in deg C from a Celsius/Fahrenheit pair of unknown order.

    Parameters
    ----------
    col_a, col_b : pandas.Series
        The two temperature columns as read from the sheet, in either order. One holds
        deg C and the other deg F, but the 2008 sheet exchanges them mid-year without
        changing the header.

    Returns
    -------
    pandas.Series
        Temperature in deg C, NaN where neither orientation is consistent.

    Notes
    -----
    Each row is tested against F = C * 9/5 + 32 in both orientations, within
    CF_TOLERANCE_F, and the consistent orientation is used. Rows consistent with neither
    are set to missing: 2,563 such rows occur in 2008, almost all inside the June and
    July failed-instrument window. Where only one column carries a value the pair test
    cannot run, so the value is treated as deg F if it exceeds 35 and deg C otherwise;
    no row in this workbook takes that path.
    """
    a = pd.to_numeric(col_a, errors="coerce")
    b = pd.to_numeric(col_b, errors="coerce")
    a_is_celsius = (b - (a * 9.0 / 5.0 + 32.0)).abs() < CF_TOLERANCE_F
    b_is_celsius = (a - (b * 9.0 / 5.0 + 32.0)).abs() < CF_TOLERANCE_F
    t_c = pd.Series(np.nan, index=a.index)
    t_c[a_is_celsius] = a[a_is_celsius]
    t_c[b_is_celsius & ~a_is_celsius] = b[b_is_celsius & ~a_is_celsius]
    lone = t_c.isna() & a.notna() & b.isna()
    t_c[lone] = np.where(a[lone] > 35.0, (a[lone] - 32.0) * 5.0 / 9.0, a[lone])
    return t_c


def parse_sheet(xls: pd.ExcelFile, sheet: str) -> Optional[pd.DataFrame]:
    """Parse one year sheet into harmonized columns and units.

    Parameters
    ----------
    xls : pandas.ExcelFile
        The open workbook.
    sheet : str
        Sheet name, one of the year sheets in the workbook.

    Returns
    -------
    pandas.DataFrame or None
        Columns ``datetime``, ``T_C`` (deg C), ``pH``, ``SpCond_uScm`` (uS/cm),
        ``DO_pctsat`` (percent), ``DO_mgL`` (mg/L), and ``source_sheet``. Rows carrying
        no measurement at all are dropped. None if the sheet holds no data rows.

    Raises
    ------
    ValueError
        If no header row containing "date" is found in the first 10 rows, which means
        the sheet layout is one this parser does not recognize.

    Notes
    -----
    Data rows are selected by requiring a parseable date in the first column rather than
    by slicing to the end of the block. That is what lets the parser step over the logger
    annotations inside the 2005, 2006, and 2007 data and the repeated header block inside
    the 2006 data. Discarded rows are counted and reported by :func:`report_skipped`.
    """
    raw = xls.parse(sheet, header=None)
    hdr = next((r for r in range(min(10, len(raw)))
                if "date" in _label(raw.iloc[r, 0]).lower()), None)
    if hdr is None:
        raise ValueError(f"sheet {sheet!r}: no header row containing 'date' in the "
                         f"first {min(10, len(raw))} rows")
    labels = [_label(v) for v in raw.iloc[hdr].tolist()]
    # Layout C carries a unit row under the header, identified by its date format cell.
    has_unit_row = _label(raw.iloc[hdr + 1, 0]).upper().startswith("M/D/")
    body = raw.iloc[hdr + (2 if has_unit_row else 1):].reset_index(drop=True)
    body = body[pd.to_datetime(body.iloc[:, 0], errors="coerce").notna()].reset_index(drop=True)
    if not len(body):
        return None

    lower = [c.lower() for c in labels]
    has_time = len(lower) > 1 and lower[1] == "time"
    ts = [_combine_datetime(body.iloc[i, 0], body.iloc[i, 1] if has_time else None)
          for i in range(len(body))]

    def by_label(*names: str) -> pd.Series:
        """First column whose header matches one of names, as numbers; NaN if absent."""
        for n in names:
            if n in lower:
                return pd.to_numeric(body.iloc[:, lower.index(n)], errors="coerce")
        return pd.Series(np.nan, index=body.index)

    # Temperature: the sheets carry either a single Fahrenheit column (2002) or a
    # Celsius/Fahrenheit pair whose order is unreliable in 2008.
    temp_cols = [i for i, c in enumerate(lower) if c.startswith("temp")]
    if len(temp_cols) >= 2:
        t_c = _resolve_temperature(body.iloc[:, temp_cols[0]], body.iloc[:, temp_cols[1]])
    elif len(temp_cols) == 1:
        v = pd.to_numeric(body.iloc[:, temp_cols[0]], errors="coerce")
        # A median above 35 can only be Fahrenheit for a hypolimnetic release.
        t_c = (v - 32.0) * 5.0 / 9.0 if v.median() > 35.0 else v
    else:
        t_c = pd.Series(np.nan, index=body.index)

    do_mgl = by_label("d.o.", "do mg/l", "ldo", "do")
    pct = by_label("dosat %", "ldo%", "do%")
    # Conductivity is reported as uS/cm except in 2004, where it is mS/cm. The 2002 header
    # ("COND.") states no unit and is taken as uS/cm; its median of 48.6 sits inside the
    # 33 to 53 uS/cm range of the labelled uS/cm sheets, which supports that reading.
    cond_idx = next((i for i, c in enumerate(lower) if c.startswith(("cond", "spcond"))), None)
    if cond_idx is None:
        cond = pd.Series(np.nan, index=body.index)
    else:
        cond = pd.to_numeric(body.iloc[:, cond_idx], errors="coerce")
        if "ms" in lower[cond_idx]:
            cond = cond * 1000.0        # mS/cm -> uS/cm

    out = pd.DataFrame({
        "datetime": ts,
        "T_C": t_c,
        "pH": by_label("ph"),
        "SpCond_uScm": cond,
        "DO_pctsat": pct,
        "DO_mgL": do_mgl,
        "source_sheet": sheet,
    })
    meas = ["T_C", "pH", "SpCond_uScm", "DO_pctsat", "DO_mgL"]
    return out[out[meas].notna().any(axis=1)].copy()


def report_skipped(xls: pd.ExcelFile, sheet: str) -> None:
    """Print the rows the sheet parser discards, so no loss is silent.

    Two classes are reported: rows whose first cell does not parse as a date, which are
    the logger annotations and the repeated header block, and rows that parse as dates
    but carry no measurement at all.
    """
    raw = xls.parse(sheet, header=None)
    hdr = next((r for r in range(min(10, len(raw)))
                if "date" in _label(raw.iloc[r, 0]).lower()), None)
    if hdr is None:
        return
    has_unit_row = _label(raw.iloc[hdr + 1, 0]).upper().startswith("M/D/")
    body = raw.iloc[hdr + (2 if has_unit_row else 1):].reset_index(drop=True)
    nodate = body[pd.to_datetime(body.iloc[:, 0], errors="coerce").isna()]
    dated = body[pd.to_datetime(body.iloc[:, 0], errors="coerce").notna()]
    blank = 0
    if len(dated):
        vals = dated.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        blank = int((~vals.notna().any(axis=1)).sum())
    # A repeated header block sits inside the 2006 data; its cells are skipped the same
    # way the annotations are, but they carry no information worth echoing.
    repeated_header = {_label(raw.iloc[hdr, 0]).lower(), "m/d/yy", "m/d/yyyy"}
    notes = [_label(v) for v in nodate.iloc[:, 0].tolist()
             if len(_label(v)) >= ANNOTATION_MIN_CHARS
             and _label(v).lower() not in repeated_header]
    print(f"  sheet {sheet}: {len(dated) - blank:,} rows kept, "
          f"{len(nodate)} without a parseable date, {blank} dated but empty")
    for n in notes:
        print(f"      note: {n}")


def flag_flatline(rec: pd.DataFrame) -> pd.Series:
    """Flag runs of identical consecutive dissolved-oxygen values.

    Parameters
    ----------
    rec : pandas.DataFrame
        Time-sorted record carrying ``DO_mgL`` in mg/L.

    Returns
    -------
    pandas.Series of bool
        True where the value belongs to a run of at least FLATLINE_N identical
        consecutive readings and is not missing.

    Notes
    -----
    A stuck sensor is the target: the 2008 summer deployment reports exactly 0.000 mg/L
    continuously from 2008-06-27 to 2008-07-24. Runs are counted over the concatenated,
    time-sorted record, so a run is not interrupted by a gap in time or by a change of
    source sheet. Missing values break a run because NaN never equals NaN.
    """
    v = rec["DO_mgL"]
    run_id = (v != v.shift()).cumsum()
    run_len = run_id.map(run_id.value_counts())
    return (run_len >= FLATLINE_N) & v.notna()


def in_sensor_decay(when: pd.Series) -> pd.Series:
    """Flag timestamps inside a SENSOR_DECAY_WINDOWS interval.

    Parameters
    ----------
    when : pandas.Series
        Timestamps, in the record's own clock.

    Returns
    -------
    pandas.Series of bool
        True where the timestamp falls in a half-open [start, end) decay window.
    """
    hit = pd.Series(False, index=when.index)
    for start, end in SENSOR_DECAY_WINDOWS:
        hit |= (when >= pd.Timestamp(start)) & (when < pd.Timestamp(end))
    return hit


def collapse_repeated_timestamps(rec: pd.DataFrame) -> pd.DataFrame:
    """Collapse rows sharing a timestamp to one row, averaging the measured channels.

    The 2004 sheet writes most of its deployment on two consecutive rows per reading, so
    without this step that deployment enters every observation-weighted statistic at
    roughly double weight while the other years enter at single weight.

    Parameters
    ----------
    rec : pandas.DataFrame
        Parsed record, already sorted by ``datetime``.

    Returns
    -------
    pandas.DataFrame
        One row per timestamp, with ``n_source_rows`` giving how many rows were
        combined. Measured channels are averaged; ``source_sheet`` is taken from the
        first row of the group.

    Notes
    -----
    Most repeated timestamps carry identical readings, for which the mean is the value
    itself. Where they disagree, both members are genuine samples logged against the
    same 30 min stamp during a period when the oxygen sensor was oscillating, so the
    mean is a better estimate of the interval than either member, and taking it makes
    the result independent of row order. The count and the largest disagreement are
    printed by the run summary, and ``n_source_rows`` keeps the collapse auditable in
    the archived CSV.
    """
    if not rec["datetime"].duplicated().any():
        rec["n_source_rows"] = 1
        return rec
    channels = ["T_C", "pH", "SpCond_uScm", "DO_pctsat", "DO_mgL"]
    grouped = rec.groupby("datetime", sort=True)
    out = grouped[channels].mean()
    out["source_sheet"] = grouped["source_sheet"].first()
    out["n_source_rows"] = grouped.size()
    return out.reset_index()


def apply_qa(rec: pd.DataFrame) -> pd.DataFrame:
    """Apply gross-range, supersaturation, and flatline quality control.

    Parameters
    ----------
    rec : pandas.DataFrame
        Concatenated parsed sheets, one row per observation.

    Returns
    -------
    pandas.DataFrame
        Time-sorted, with ``year``, ``month``, ``DO_pctsat_est`` (percent),
        ``do_flag``, ``DO_mgL_qa`` (mg/L), and ``ph_flag`` added. Temperature outside
        [T_MIN, T_MAX] is set to NaN in place.

    Notes
    -----
    ``do_flag`` takes one of five values, assigned in this order so that later
    assignments win: "range" for a reading at or below DO_MIN, above DO_MAX, or above
    DO_SUPERSAT_PCT percent saturation; "flatline" for a stuck sensor; "sensor_decay"
    for a reading inside a SENSOR_DECAY_WINDOWS interval; "missing" for no reading;
    "good" otherwise. ``DO_mgL_qa`` carries the value only where the flag is "good", so
    every downstream script sees the same screened series and the screening stays
    auditable in the archived CSV.

    Rows sharing a timestamp are collapsed before any screen is applied, so a flag is
    assigned once per instant rather than once per source row.

    A reading of exactly 0.000 mg/L is treated as an instrument artifact rather than a
    measurement, which is why the lower bound is exclusive.
    """
    rec = rec.dropna(subset=["datetime"]).copy()
    rec["datetime"] = pd.to_datetime(rec["datetime"])
    # A stable sort keeps rows that share a timestamp in the order the sheets supply
    # them, so the written record does not depend on the sort implementation.
    rec = rec.sort_values("datetime", kind="stable").reset_index(drop=True)
    rec = collapse_repeated_timestamps(rec)
    rec["year"] = rec["datetime"].dt.year
    rec["month"] = rec["datetime"].dt.month

    # Temperature QA. Unlike oxygen, a screened temperature carries no flag into the CSV:
    # it is written as an empty cell, indistinguishable from a temperature the sheet never
    # reported, and the count is not echoed by the run summary either. The 40 records this
    # removes are described at T_MAX above; tailrace_sonde_ingest.py prints the equivalent
    # count for the Allatoona record.
    rec.loc[(rec["T_C"] < T_MIN) | (rec["T_C"] > T_MAX), "T_C"] = np.nan

    # Percent saturation: use the reported value where present, otherwise compute it
    # from the paired temperature so that the supersaturation screen applies to every
    # year of the record.
    sat = do_saturation(rec["T_C"])
    computed_pct = rec["DO_mgL"] / sat * 100.0
    rec["DO_pctsat_est"] = rec["DO_pctsat"].where(rec["DO_pctsat"].notna(), computed_pct)

    do_flag = pd.Series("good", index=rec.index)
    do_flag[(rec["DO_mgL"] <= DO_MIN) | (rec["DO_mgL"] > DO_MAX)] = "range"
    do_flag[rec["DO_pctsat_est"] > DO_SUPERSAT_PCT] = "range"
    do_flag[flag_flatline(rec)] = "flatline"
    do_flag[in_sensor_decay(rec["datetime"])] = "sensor_decay"
    do_flag[rec["DO_mgL"].isna()] = "missing"
    rec["do_flag"] = do_flag
    rec["DO_mgL_qa"] = rec["DO_mgL"].where(do_flag == "good", np.nan)

    # Separate pH reliability flag. pH is not used quantitatively in this study; the
    # flag documents the affected records rather than removing them.
    rec["ph_flag"] = np.where(rec["pH"] < PH_PHYSICAL_MIN, "suspect", "ok")
    rec.loc[rec["pH"].isna(), "ph_flag"] = "missing"
    return rec


def summarize(rec: pd.DataFrame) -> None:
    """Print the coverage, QA disposition, and summer DO statistics of the record."""
    good = rec[rec["do_flag"] == "good"]
    print("=" * 78)
    print("BUFORD TAILRACE MONITOR, 2002-2008")
    print("=" * 78)
    print(f"Total records: {len(rec):,}   "
          f"{rec['datetime'].min()} -> {rec['datetime'].max()}")
    step = rec["datetime"].diff().dt.total_seconds().median() / 60.0
    print(f"Median sampling interval: {step:.0f} min")
    print(f"\nDO-flag counts:\n{rec['do_flag'].value_counts().to_string()}")

    merged = rec[rec["n_source_rows"] > 1]
    if len(merged):
        by_sheet = merged["source_sheet"].value_counts()
        print(f"\nRepeated timestamps collapsed: {len(merged):,} timestamps built from "
              f"{int(merged['n_source_rows'].sum()):,} source rows "
              f"(by sheet: {by_sheet.to_dict()})")

    decayed = rec[rec["do_flag"] == "sensor_decay"]
    if len(decayed):
        print(f"Sensor-decay windows excluded: {len(decayed):,} records, "
              f"{decayed['datetime'].min()} to {decayed['datetime'].max()}, "
              f"DO {decayed['DO_mgL'].min():.2f} to {decayed['DO_mgL'].max():.2f} mg/L")

    print("\nRecords and DO-QA disposition by year:")
    g = rec.groupby("year").agg(
        n=("datetime", "size"),
        T_med=("T_C", "median"),
        DO_good=("DO_mgL_qa", lambda s: s.notna().sum()),
        DO_med_qa=("DO_mgL_qa", "median"),
    )
    print(g.to_string())

    summer = good[good["month"].isin([6, 7, 8, 9])]
    print("\n" + "=" * 78)
    print("QUALITY-CONTROLLED SUMMER (Jun-Sep) TAILRACE DO")
    print("=" * 78)
    print(f"  n = {len(summer):,} sub-hourly observations over "
          f"{summer['datetime'].dt.date.nunique()} days")
    print(f"  median = {summer['DO_mgL_qa'].median():.2f} mg/L")
    print(f"  mean   = {summer['DO_mgL_qa'].mean():.2f} mg/L")
    print(f"  min    = {summer['DO_mgL_qa'].min():.2f} mg/L")
    for thr in (2.0, 5.0, 6.0):
        print(f"  % summer observations < {thr:.1f} mg/L: "
              f"{100.0 * (summer['DO_mgL_qa'] < thr).mean():.1f}%")

    daily = good.set_index("datetime")["DO_mgL_qa"].resample("D").mean().dropna()
    ds = daily[daily.index.month.isin([6, 7, 8, 9])]
    print(f"\n  Daily means, summer: n = {len(ds)} days, median "
          f"{ds.median():.2f} mg/L; {100.0 * (ds < 6.0).mean():.1f}% below the "
          f"6.0 mg/L daily-average trout criterion")

    # Diel range needs a well-sampled day: 20 observations is under a quarter of a 15 min
    # day but five sixths of a 2002 hourly day, so the threshold admits every year of the
    # record rather than dropping the coarse early years.
    cnt = good.set_index("datetime")["DO_mgL_qa"].resample("D").count()
    rng = (good.set_index("datetime")["DO_mgL_qa"].resample("D").max()
           - good.set_index("datetime")["DO_mgL_qa"].resample("D").min())[cnt >= 20]
    rs = rng[rng.index.month.isin([6, 7, 8, 9])]
    print(f"  Summer diel DO range: median {rs.median():.2f} mg/L "
          f"(90th percentile {rs.quantile(0.9):.2f}, n = {len(rs)} days)")

    seas = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
            6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}
    s = good.assign(season=good["month"].map(seas)).groupby("season")["DO_mgL_qa"].median()
    print("\n  Seasonal DO medians (mg/L): "
          + ", ".join(f"{k} {s[k]:.2f}" for k in ("DJF", "MAM", "JJA", "SON") if k in s))

    # Sustained near-anoxic readings survive the gross-range screen by design, because
    # an anoxic hypolimnetic release is physically possible here. They are reported so
    # that a reader can see how much of the low tail rests on a few episodes.
    low = good[good["DO_mgL_qa"] < 1.0]
    if len(low):
        months = (low.groupby([low["year"], low["month"]]).size()
                  .rename("n").reset_index())
        spans = ", ".join(f"{int(r.year)}-{int(r.month):02d} n={int(r.n)}"
                          for r in months.itertuples())
        print(f"\n  Retained readings below 1.0 mg/L: {len(low):,} over "
              f"{low['datetime'].dt.date.nunique()} days ({spans})")


def main() -> None:
    """Parse every sheet, apply QA, write the tidy CSV, and print the summary."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--src", default=SRC_XLS,
                    help="source workbook (default: Data/Lanier_2002-2010.xls)")
    ap.add_argument("--out", default=OUT_CSV,
                    help="tidy CSV to write (default: "
                         "Data/buford_tailrace_monitor_2002_2008.csv)")
    args = ap.parse_args()

    if not os.path.exists(args.src):
        raise SystemExit(f"source workbook not found: {args.src}")

    xls = pd.ExcelFile(args.src)
    print(f"Reading {args.src}: sheets {', '.join(xls.sheet_names)}")
    for sheet in xls.sheet_names:
        report_skipped(xls, sheet)
    frames = [f for f in (parse_sheet(xls, s) for s in xls.sheet_names)
              if f is not None and len(f)]
    if not frames:
        raise SystemExit(f"no data rows parsed from {args.src}")
    rec = apply_qa(pd.concat(frames, ignore_index=True))
    cols = ["datetime", "year", "month", "T_C", "pH", "ph_flag", "SpCond_uScm",
            "DO_pctsat", "DO_pctsat_est", "DO_mgL", "do_flag", "DO_mgL_qa",
            "n_source_rows", "source_sheet"]
    rec[cols].to_csv(args.out, index=False)
    print(f"Wrote {args.out} ({len(rec):,} records)\n")
    summarize(rec)


if __name__ == "__main__":
    main()
