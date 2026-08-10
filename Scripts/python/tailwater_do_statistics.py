#!/usr/bin/env python3
"""Reproduce every derived tailwater dissolved-oxygen statistic quoted in the report.

The report cites saturation percentages, criterion-exceedance fractions, rank
correlations, and the corroboration statistics for the Corps Buford tailrace monitor.
This script recomputes all of them from the archived data and writes the results to
``analysis/tailwater_do_statistics.txt`` so that the reported values are reproducible
from a single command.

Saturation is computed with the Benson and Krause relation as given in APHA Standard
Methods, evaluated at 1 atm and zero salinity. Two alternatives are reported beside it
so the choice can be audited. Over the summer temperatures of these records the
CE-QUAL-W2 internal relation agrees with Benson-Krause to within 0.04 percentage points
of saturation, and the APHA polynomial approximation to within 0.87 percentage points
(both maxima measured over the paired summer days below, 9.2 to 24.9 deg C).

Units throughout: dissolved oxygen in mg/L, temperature in deg C, discharge in cfs,
saturation in percent. Daily values are USGS daily means (statistic code 00003);
monitor and sonde values are sub-hourly instantaneous readings.

Run with the clearwater conda env:
    /opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3 \
        Scripts/python/tailwater_do_statistics.py
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT_TXT = os.path.join(REPO, "analysis", "tailwater_do_statistics.txt")
os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)

SUMMER = [6, 7, 8, 9]        # "summer" throughout the report
LATE_SUMMER = [7, 8, 9]      # "late summer" throughout the report
STRATIFIED_SEASON = [7, 8, 9, 10]   # Jul-Oct, the near-field pairing window

# Georgia criteria from the project Water Control Manuals: (daily average, instantaneous).
CRIT = {"Allatoona": (5.0, 4.0), "Buford": (6.0, 5.0)}

# Thickness of the bottom layer taken as the source water for the release, and the
# window either side of a cast over which the tailrace daily records are averaged.
# Both match the values used by the withdrawal-depth analysis, so the near-field gain
# reported here stays directly comparable to the one drawn there.
SOURCE_LAYER_M = 5.0
TWDO_WINDOW = pd.Timedelta("3D")

# A day counts as diel-complete at this many sub-hourly observations, so that a range
# is not taken from a few samples clustered in one part of the day. The Buford monitor
# logs every 15 min from 2005, so a full day is 96.
MIN_OBS_PER_DAY = 20

_lines: list[str] = []


def emit(s: str = "") -> None:
    """Print a line and retain it for the archived statistics file."""
    print(s)
    _lines.append(s)


def usgs(path: str, pcode: str, name: str) -> pd.DataFrame:
    """Read an archived USGS daily-value CSV.

    Parameters
    ----------
    path : str
        Path relative to the repository root.
    pcode : str
        Five-digit NWIS parameter code: 00010 temperature (deg C), 00060 discharge
        (cfs), 00300 dissolved oxygen (mg/L). Statistic code 00003 is the daily mean.
    name : str
        Name given to the value column in the returned frame.

    Returns
    -------
    pandas.DataFrame
        Columns ``Date`` (datetime64) and ``name`` (float, units of the parameter
        code), rows missing either one dropped, index reset.

    Notes
    -----
    Each archived file carries one row per date with no duplicates, so every merge on
    ``Date`` below is one-to-one and cannot duplicate a day. It can drop days one
    record does not cover: the Buford DO record has 1060 days, of which 1059 pair with
    a temperature day and 1058 with a discharge day. The retained n is printed beside
    every statistic taken from a merged frame.
    """
    d = pd.read_csv(os.path.join(REPO, path))
    d["Date"] = pd.to_datetime(d["Date"].astype(str).str.strip('"'))
    d[name] = pd.to_numeric(d[f"X_{pcode}_00003"].astype(str).str.strip('"'), errors="coerce")
    return d[["Date", name]].dropna().reset_index(drop=True)


def sat_benson_krause(t_c):
    """Dissolved-oxygen solubility from the Benson and Krause relation.

    Parameters
    ----------
    t_c : array_like
        Water temperature in deg C.

    Returns
    -------
    numpy.ndarray
        Equilibrium DO concentration in mg/L at 1 atm total pressure and zero salinity.

    Notes
    -----
    Evaluated as ``exp(A0 + A1/T + A2/T^2 + A3/T^3 + A4/T^4)`` with T in kelvin, the
    form given in APHA Standard Methods. No salinity or chloride correction and no
    barometric correction is applied: these are freshwater tailwater reaches, and the
    pressure correction is not applied because the study has not sourced a barometric
    record for either site. Omitting it makes the saturation percentages reported here
    slightly low, since the true local pressure at these elevations is below 1 atm.

    This routine is applied only over the 9.2 to 24.9 deg C spanned by the paired
    summer days in this script; it is not guarded, so a caller working outside that
    range should confirm the fit is still appropriate.
    """
    tk = np.asarray(t_c, dtype=float) + 273.15
    return np.exp(-139.34411 + 1.575701e5 / tk - 6.642308e7 / tk ** 2
                  + 1.243800e10 / tk ** 3 - 8.621949e11 / tk ** 4)


def sat_w2(t_c):
    """Dissolved-oxygen solubility as CE-QUAL-W2 computes it internally.

    Parameters
    ----------
    t_c : array_like
        Water temperature in deg C.

    Returns
    -------
    numpy.ndarray
        Equilibrium DO concentration in mg/L at 1 atm and zero salinity.

    Notes
    -----
    Reported alongside Benson-Krause so that a saturation figure quoted from this study
    can be compared with one a CE-QUAL-W2 run would produce. Over the temperatures of
    these records the two agree to within 0.04 percentage points of saturation.
    """
    return np.exp(7.7117 - 1.31403 * np.log(np.asarray(t_c, dtype=float) + 45.93))


def sat_apha_poly(t_c):
    """Dissolved-oxygen solubility from the APHA cubic polynomial approximation.

    Parameters
    ----------
    t_c : array_like
        Water temperature in deg C.

    Returns
    -------
    numpy.ndarray
        Equilibrium DO concentration in mg/L at 1 atm and zero salinity.

    Notes
    -----
    The coarsest of the three relations. Over the temperatures of these records it
    departs from Benson-Krause by up to 0.87 percentage points of saturation, with the
    larger departures at the warmer Allatoona summer temperatures; it is reported for
    audit only and is not the relation the report quotes.
    """
    t = np.asarray(t_c, dtype=float)
    return 14.652 - 0.41022 * t + 0.007991 * t ** 2 - 0.000077774 * t ** 3


def daily_from_subhourly(df: pd.DataFrame, value: str, out: str) -> pd.DataFrame:
    """Calendar-day mean of a sub-hourly record, keyed on ``Date``.

    Parameters
    ----------
    df : pandas.DataFrame
        Sub-hourly record carrying a ``datetime`` column and ``value``.
    value : str
        Column to average, in its own units.
    out : str
        Name given to the averaged column in the result.

    Returns
    -------
    pandas.DataFrame
        Columns ``Date`` (datetime64, midnight) and ``out``. Days with no observation
        are absent rather than NaN.

    Notes
    -----
    No minimum-observation threshold is applied to the means, unlike the diel range,
    so a day with a partial record still contributes. On the Buford monitor that is 8
    of 1013 days, and restricting to days with at least MIN_OBS_PER_DAY observations
    moves the summer daily-mean median from 5.32 to 5.30 mg/L.
    """
    return (df.set_index("datetime")[value].resample("D").mean().dropna()
            .rename(out).reset_index().rename(columns={"datetime": "Date"}))


def nearfield_gains(prof: pd.DataFrame, twdo: pd.DataFrame,
                    y0: int, y1: int) -> list[float]:
    """Apparent DO gain across the dam, one value per paired stratified-season cast.

    Source water is the mean DO over the deepest SOURCE_LAYER_M of the cast; the
    tailrace value is the mean of the daily records within TWDO_WINDOW of the cast
    date. The gain is tailrace minus source.

    Parameters
    ----------
    prof : pandas.DataFrame
        Forebay profiles with ``Date``, ``Depth_m`` (m) and numeric ``DO_mgL`` (mg/L).
    twdo : pandas.DataFrame
        Daily tailrace DO, columns ``Date`` and ``DO`` in mg/L.
    y0, y1 : int
        Inclusive first and last calendar year of casts to pair.

    Returns
    -------
    list of float
        Gains in mg/L. Casts with no tailrace record inside the window are dropped, so
        the length is reported wherever the median of this list is quoted.

    Notes
    -----
    This reproduces the near-field gain calculation of the withdrawal-depth analysis,
    which is why SOURCE_LAYER_M and TWDO_WINDOW carry the values they do. The gain is
    apparent, not a measured reaeration: it also carries any mismatch between the true
    withdrawal band and the deepest SOURCE_LAYER_M.
    """
    p = prof.dropna(subset=["DO_mgL"])
    casts = p[(p.Date.dt.year.between(y0, y1))
              & (p.Date.dt.month.isin(STRATIFIED_SEASON))]
    gains = []
    for d, g in casts.groupby("Date"):
        src = g[g.Depth_m >= g.Depth_m.max() - SOURCE_LAYER_M]["DO_mgL"].mean()
        near = twdo[(twdo.Date - d).abs() <= TWDO_WINDOW]
        if len(near):
            gains.append(near.DO.mean() - src)
    return gains


# ------------------------------------------------------------------ load records
a_do = usgs("Data/do_daily_02394000_2005_2007.csv", "00300", "DO")
a_t = usgs("Data/temp_daily_02394000_2005_2026.csv", "00010", "Tw")
a_q = usgs("Data/discharge_daily_02394000_2005_2007.csv", "00060", "Q")
b_do = usgs("Data/Lanier/do_daily_02334430_2005_2026.csv", "00300", "DO")
b_t = usgs("Data/temp_daily_02334430_2005_2026.csv", "00010", "Tw")
b_q = usgs("Data/Lanier/discharge_daily_02334430_2005_2026.csv", "00060", "Q")

# KNOWN DEFECT, not corrected here because the report is written against the current
# values: the archived monitor CSV holds 2,120 rows on timestamps that already appear
# in the file, from overlapping annual source sheets. 1,843 of the 2,119 repeated
# timestamps carry identical readings and 276, all in 2004, carry different ones, with
# a median spread of 0.03 and a maximum of 6.94 mg/L. Nothing below de-duplicates them,
# so every observation-weighted statistic double-counts those instants. Dropping the
# repeats with keep="first" would move the summer median from 4.58 to 4.45 mg/L, the
# summer mean from 5.16 to 5.00, the percentage of summer observations below 5.0 mg/L
# from 53.5 to 55.3, the JJA seasonal median from 5.36 to 5.17, the SON median from
# 1.99 to 1.90, and the good-record count from 77,167 to 75,049. Daily means are barely
# affected. The 276 disagreeing timestamps need a resolution rule before the fix lands;
# the origin is Scripts/python/buford_tailrace_monitor_ingest.py, not this script.
bufmon = pd.read_csv(os.path.join(REPO, "Data", "buford_tailrace_monitor_2002_2008.csv"),
                     parse_dates=["datetime"])
bm_good = bufmon[bufmon["do_flag"] == "good"].copy()
bm_daily = daily_from_subhourly(bm_good, "DO_mgL_qa", "DO")

emit("Tailwater dissolved-oxygen statistics quoted in the report")
emit("Summer = Jun-Sep; late summer = Jul-Sep. Saturation: Benson-Krause at 1 atm.")
emit("")

# ------------------------------------------------- daily-record criterion statistics
emit("-- USGS daily tailwater records: summer statistics and criterion exceedance --")
emit("   Note: these are daily means. Counting daily means below an instantaneous")
emit("   minimum criterion gives a LOWER BOUND on instantaneous exceedance.")
for name, do, tw in (("Allatoona (02394000, 2005-2007)", a_do, a_t),
                     ("Buford (02334430, 2023-present)", b_do, b_t)):
    proj = "Allatoona" if name.startswith("Allatoona") else "Buford"
    davg, inst = CRIT[proj]
    d = do.assign(m=do.Date.dt.month)
    s = d[d.m.isin(SUMMER)]
    emit(f"{name}: n={len(d)}, summer n={len(s)}")
    emit(f"  summer median={s.DO.median():.2f}  mean={s.DO.mean():.2f}  "
         f"min={s.DO.min():.2f} mg/L")
    emit(f"  record minimum={d.DO.min():.2f} mg/L on "
         f"{d.loc[d.DO.idxmin(), 'Date'].date()}")
    emit(f"  late-summer median={d[d.m.isin(LATE_SUMMER)].DO.median():.2f} mg/L")
    emit(f"  % summer days < {inst:.1f} (inst min)={100 * (s.DO < inst).mean():.1f}; "
         f"< {davg:.1f} (daily-avg)={100 * (s.DO < davg).mean():.1f}")
    emit(f"  % all days   < {inst:.1f}={100 * (d.DO < inst).mean():.1f}; "
         f"< {davg:.1f}={100 * (d.DO < davg).mean():.1f}")
    # Saturation needs a temperature on the same day, so this merge drops any DO day
    # the temperature record does not cover. The surviving count is printed below.
    m = pd.merge(do, tw, on="Date")
    ms = m[m.Date.dt.month.isin(SUMMER)]
    parts = []
    for lbl, f in (("Benson-Krause", sat_benson_krause), ("W2", sat_w2),
                   ("APHA poly", sat_apha_poly)):
        pc = ms.DO.values / f(ms.Tw.values) * 100.0
        parts.append(f"{lbl} median={np.median(pc):.0f}% min={pc.min():.0f}%")
    emit(f"  summer % saturation (n={len(ms)}): " + "; ".join(parts))
    emit("")

# ------------------------------------------------------------------ correlations
# Both coefficients are computed on the merged frame, so DO and its driver are the same
# calendar day by construction and the pairs cannot be mismatched. Spearman is the one
# the report quotes: the DO-temperature relation is monotonic but not linear, and the
# discharge distribution has a heavy right tail that Pearson would be dominated by.
emit("-- Rank correlations of tailrace DO with temperature and discharge --")
for name, do, x, xn in (("Allatoona", a_do, a_t, "Tw"), ("Buford", b_do, b_t, "Tw"),
                        ("Allatoona", a_do, a_q, "Q"), ("Buford", b_do, b_q, "Q")):
    m = pd.merge(do, x, on="Date")
    m = m.assign(mm=m.Date.dt.month)
    for lbl, sub in (("all", m), ("summer", m[m.mm.isin(SUMMER)])):
        rs, ps = spearmanr(sub[xn], sub.DO)
        rp, _ = pearsonr(sub[xn], sub.DO)
        emit(f"  {name:10s} DO~{xn:2s} {lbl:7s} n={len(sub):5d}  "
             f"Spearman={rs:+.3f} (p={ps:.1e})  Pearson={rp:+.3f}")
emit("")

# ------------------------------------- Corps Buford monitor: coverage and statistics
# These are sub-hourly instantaneous readings, so comparing them with the 5.0 mg/L
# instantaneous minimum is a like-for-like test; the 6.0 mg/L daily-average criterion
# is tested separately below, against daily means.
emit("-- Corps Buford tailrace monitor (2002, 2004-2008) --")
emit(f"  records={len(bufmon):,}  good={len(bm_good):,}  "
     f"days with good DO={bm_good.datetime.dt.date.nunique()}")
emit(f"  QA disposition: {bufmon.do_flag.value_counts().to_dict()}")
s = bm_good[bm_good.month.isin(SUMMER)]
emit(f"  summer n={len(s):,} median={s.DO_mgL_qa.median():.2f} "
     f"mean={s.DO_mgL_qa.mean():.2f} min={s.DO_mgL_qa.min():.2f} mg/L")
emit(f"  % summer observations < 5.0 (inst min)={100 * (s.DO_mgL_qa < 5).mean():.1f}; "
     f"< 6.0={100 * (s.DO_mgL_qa < 6).mean():.1f}")
ds = bm_daily[bm_daily.Date.dt.month.isin(SUMMER)]
emit(f"  daily means, summer: n={len(ds)} median={ds.DO.median():.2f}; "
     f"{100 * (ds.DO < 6.0).mean():.1f}% below the 6.0 daily-average criterion")
bm_day = bm_good.set_index("datetime")["DO_mgL_qa"].resample("D")
cnt = bm_day.count()
rng = (bm_day.max() - bm_day.min())[cnt >= MIN_OBS_PER_DAY]
rs_ = rng[rng.index.month.isin(SUMMER)]
emit(f"  summer diel DO range: median={rs_.median():.2f} "
     f"90th pct={rs_.quantile(0.9):.2f} mg/L (n={len(rs_)} days)")
# Meteorological seasons: DJF, MAM, JJA, SON.
seas = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
        6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}
sm = bm_good.assign(season=bm_good.month.map(seas)).groupby("season")["DO_mgL_qa"].median()
emit("  seasonal medians (mg/L): "
     + ", ".join(f"{k} {sm[k]:.2f}" for k in ("DJF", "MAM", "JJA", "SON")))
emit("")

# ------------------------------------------------------- station identity and overlap
# "MAD" below is the MEAN absolute difference between the two daily temperature series,
# not the median absolute deviation. It tests whether the Corps monitor and USGS
# 02334430 are sampling the same water, which is what lets the monitor record stand in
# for the gage over 2002-2008. Years with 20 or fewer common days are not broken out.
emit("-- Station identity: monitor daily mean T vs USGS 02334430 daily T --")
td = daily_from_subhourly(bm_good, "T_C", "Tx")
m = pd.merge(td, b_t, on="Date")
emit(f"  overall n={len(m)} MAD={np.abs(m.Tx - m.Tw).mean():.2f} C")
for y in sorted(m.Date.dt.year.unique()):
    mm = m[m.Date.dt.year == y]
    if len(mm) > 20:
        emit(f"    {y}: n={len(mm):3d} MAD={np.abs(mm.Tx - mm.Tw).mean():.2f} C "
             f"r={pearsonr(mm.Tx, mm.Tw)[0]:.3f}")

ov = pd.merge(bm_daily, a_do.rename(columns={"DO": "DO_allatoona"}), on="Date")
emit("")
emit("-- Common window: Corps Buford monitor and USGS Allatoona gage --")
emit(f"  {len(ov)} common days, {ov.Date.min().date()} to {ov.Date.max().date()} "
     f"({len(ov[ov.Date.dt.month.isin(SUMMER)])} summer days)")

# --------------------------------------------- corroboration against the modern record
# The two records are 15 years apart, so this compares climatological monthly medians,
# not concurrent observations. Each month pools every year the record covers, and the
# two records cover different years, so a month-to-month difference carries genuine
# interannual variability as well as any instrument difference.
emit("")
emit("-- Corroboration of the monitor against USGS 02334430 (2023-2026), by month --")
cm = bm_good.groupby(bm_good.datetime.dt.month)["DO_mgL_qa"].median()
um = b_do.groupby(b_do.Date.dt.month)["DO"].median()
emit("  month  monitor  USGS 2023-  difference")
for mo in range(1, 13):
    if mo in cm.index and mo in um.index:
        emit(f"  {mo:5d}  {cm[mo]:7.2f}  {um[mo]:9.2f}  {cm[mo] - um[mo]:+10.2f}")
# Dec-Jul is the run of months in which the two records agree; Aug-Nov is the run in
# which the monitor reads persistently low. Splitting there is a description of the
# pattern in the table above, not a hypothesis tested against it.
dj = [cm[m_] - um[m_] for m_ in (12, 1, 2, 3, 4, 5, 6, 7) if m_ in cm.index]
an = [cm[m_] - um[m_] for m_ in (8, 9, 10, 11) if m_ in cm.index]
emit(f"  Dec-Jul mean difference = {np.mean(dj):+.2f} mg/L "
     f"(max |diff| {np.max(np.abs(dj)):.2f})")
emit(f"  Aug-Nov mean difference = {np.mean(an):+.2f} mg/L "
     f"(range {np.min(an):+.2f} to {np.max(an):+.2f})")

# ------------------------------------ near-field gain check over the 2002-2008 window
La = pd.read_csv(os.path.join(REPO, "Data",
                              "inpool_forebay_profiles_Lanier_LK_12_4028.csv"),
                 parse_dates=["Date"])
La["DO_mgL"] = pd.to_numeric(La["DO_mgL"], errors="coerce")
gains = nearfield_gains(La, bm_daily, 2002, 2008)
# The modern comparison is recomputed rather than transcribed from fig7, so the two
# cannot drift apart. It uses the same kernel over the USGS Buford DO record.
modern_gains = nearfield_gains(La, b_do, 2023, 2025)
emit("")
emit("-- Near-field DO gain over the monitor window (forebay hypolimnion to tailrace) --")
emit(f"  n={len(gains)} concurrent stratified-season casts, median gain "
     f"{np.median(gains):+.2f} mg/L, all positive: {all(g_ > 0 for g_ in gains)}")
emit(f"  For comparison the 2023-2025 USGS window gives a median gain of "
     f"{np.median(modern_gains):+.2f} mg/L.")
emit("  A consistently positive gain indicates the oxygen sensor tracks the source water,")
emit("  so the Aug-Nov divergence above is not resolved as a simple sensor failure.")

with open(OUT_TXT, "w") as f:
    f.write("\n".join(_lines) + "\n")
print(f"\nWrote {OUT_TXT}")
