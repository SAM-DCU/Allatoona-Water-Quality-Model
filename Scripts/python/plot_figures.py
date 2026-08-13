#!/usr/bin/env python3
"""Forebay and tailwater figures for the Lanier/Allatoona DO study.

Figures are numbered to their order in the report, so this script writes them out of
sequence: the two forebay figures and the data-availability timeline belong to the
study-area chapter and the four tailwater figures to the model-capabilities chapter.

  fig1_forebay_profiles       representative late-summer forebay profiles
  fig2_forebay_do_seasonal    forebay DO by month and depth, full period of record
  fig8_data_timeline          periods of record by data category
  fig9_tailwater_do_windows   the four continuous tailwater DO records
  fig10_tailrace_do_climatology   Allatoona sonde monthly DO and diel range
  fig11_tailrace_subhourly    Allatoona sonde at native 30-minute resolution
  fig12_tailrace_do_drivers   DO against temperature and discharge, and saturation

Units throughout: dissolved oxygen in mg/L, temperature in deg C, depth in m below the
water surface, discharge in cfs, saturation in percent. Sonde and monitor timestamps
are as delivered by the source records and are not shifted between zones; a "day" is a
calendar day in that same convention.

The numerical kernels (select_representative_cast, monthly_depth_mean, daily_mean,
diel_range, break_gaps) take and return plain frames and series, so any value quoted
below can be recomputed without drawing a figure.

Also writes analysis/tailrace_sonde_stats.txt, which it overwrites in full: no other
script may append to that file. Figures go to analysis/figures/ as PDF (for LaTeX) and
PNG (preview); copy the PDFs into report/latex/src/images/ when they change.

Run with the clearwater conda env:
    /opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3 \
        Scripts/python/plot_figures.py
"""
import os
from typing import Sequence

import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Patch

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIG = os.path.join(REPO, "analysis", "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25,
                     "figure.dpi": 130, "savefig.bbox": "tight"})

# Month groupings used across this study, 1-based, and the meteorological seasons.
SUMMER = [6, 7, 8, 9]        # Jun-Sep, "summer" in the report
LATE_SUMMER = [7, 8, 9]      # Jul-Sep, "late summer" in the report
# fig1's representative cast is drawn from August alone. September was included until a
# seasonal breakdown showed the two lakes diverge within late summer: the median peak
# temperature gradient at Allatoona falls from 1.13 deg C/m in August to 0.57 in
# September as the thermocline erodes, while Lanier holds 3.34 and 3.15. Pooling the two
# months gives Allatoona a median of 0.80, which describes neither of them.
PEAK_STRATIFICATION = [8]

# Eligibility floor for the representative cast: a cast must carry at least this many
# depth observations complete in both temperature and oxygen, and must reach at least
# this fraction of the deepest cast in the selection months.
MIN_CAST_POINTS = 15
MIN_CAST_DEPTH_FRAC = 0.75
SEASONS = {"DJF": [12, 1, 2], "MAM": [3, 4, 5], "JJA": [6, 7, 8], "SON": [9, 10, 11]}

# A day counts as diel-complete at this many sub-hourly observations. The Allatoona
# sonde logs every 30 min, so a full day is 48; 20 keeps days with a partial deployment
# while rejecting days that would report a spurious within-day range from a few samples
# clustered in one part of the day.
MIN_OBS_PER_DAY = 20

FOREBAY = {
    "Lake Lanier (Buford Dam forebay)":
        os.path.join(REPO, "Data", "inpool_forebay_profiles_Lanier_LK_12_4028.csv"),
    "Lake Allatoona (Allatoona Dam forebay)":
        os.path.join(REPO, "Data", "inpool_forebay_profiles_Allatoona_LK_14_4494.csv"),
}


def load_profiles(p: str) -> pd.DataFrame:
    """Read a forebay profile CSV and add the calendar month.

    Parameters
    ----------
    p : str
        Absolute path to a tidy forebay profile CSV written by forebay_profiles.py.

    Returns
    -------
    pandas.DataFrame
        Columns ``Date`` (datetime64), ``Depth_m`` (m below surface), ``Temp_C``
        (deg C), ``DO_mgL`` (mg/L), ``month`` (1-based). Non-numeric temperature and
        oxygen entries become NaN rather than raising, so a cast that reports only one
        of the two survives and is filtered per figure by the caller.
    """
    df = pd.read_csv(p, parse_dates=["Date"])
    for c in ("Temp_C","DO_mgL"): df[c] = pd.to_numeric(df[c], errors="coerce")
    df["month"] = df["Date"].dt.month
    return df

prof = {k: load_profiles(v) for k, v in FOREBAY.items()}


def save(fig: Figure, name: str) -> None:
    """Write a figure as PDF (for LaTeX) and PNG (preview) to analysis/figures/."""
    for ext in ("pdf","png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig); print(f"  wrote {name}.pdf/.png")


def q(s, p: float) -> float:
    """Percentile p (0-100) of s, ignoring NaN. Units are those of s."""
    return float(np.nanpercentile(s, p))


def select_representative_cast(df: pd.DataFrame,
                               months: Sequence[int] = PEAK_STRATIFICATION) -> pd.Timestamp:
    """Date of the cast whose stratification is most typical of the given months.

    Eligibility is settled first: a cast must carry at least MIN_CAST_POINTS depth
    observations complete in both temperature and oxygen, and must reach
    MIN_CAST_DEPTH_FRAC of the deepest cast in those months. Among the eligible casts,
    the one returned is that whose peak temperature gradient lies closest to the median
    peak temperature gradient of the eligible set.

    An earlier version ranked on completeness alone and returned the deepest, most
    finely sampled cast. Completeness is not typicality. At Allatoona that rule returned
    a profile with no mixed layer, whose peak gradient of 0.80 deg C/m sits below the
    August interquartile range of 0.92 to 1.40, so the figure understated the thermocline
    that the surrounding text describes and made the sharp oxycline look unsupported by
    the temperature structure.

    Parameters
    ----------
    df : pandas.DataFrame
        Forebay profiles as returned by :func:`load_profiles`.
    months : sequence of int, optional
        Calendar months to choose from, 1-based. Defaults to PEAK_STRATIFICATION.

    Returns
    -------
    pandas.Timestamp
        Date of the selected cast.

    Notes
    -----
    The peak gradient is the largest single-interval temperature decrease per metre, so
    it is sensitive to the vertical sampling interval: a sparsely sampled cast has wide
    intervals and a smoothed gradient. The completeness floor keeps such casts out of
    the comparison rather than letting them depress the median.

    If no cast clears the floor the function falls back to the completeness ranking
    rather than raising, so a thin record still yields a figure.

    forebay_profile_gallery.py mirrors this selection and prints the cast it would
    choose, so the gallery and the report figure stay traceable to each other.
    """
    late = df[df.month.isin(months)].dropna(subset=["Temp_C","DO_mgL"])
    z_deepest = late.Depth_m.max()
    dates: list[pd.Timestamp] = []
    peaks: list[float] = []
    for date, g in late.groupby("Date"):
        g = g.sort_values("Depth_m")
        if len(g) < MIN_CAST_POINTS or g.Depth_m.max() < MIN_CAST_DEPTH_FRAC * z_deepest:
            continue
        z, t = g.Depth_m.values, g.Temp_C.values
        dates.append(date)
        peaks.append(float((-np.diff(t) / np.diff(z)).max()))
    if not dates:
        return (late.groupby("Date")
                    .agg(n=("Depth_m","size"), zmax=("Depth_m","max"))
                    .sort_values(["n","zmax"]).index[-1])
    p = np.asarray(peaks)
    return dates[int(np.argmin(np.abs(p - np.median(p))))]


def monthly_depth_mean(df: pd.DataFrame, value: str = "DO_mgL") -> pd.DataFrame:
    """Mean of ``value`` on a 1 m depth by calendar month grid, over the whole record.

    Depths are binned by taking the floor, so bin k holds every observation from k to
    k+1 m. Bins with no observation in a given month stay NaN and are drawn as gaps
    rather than interpolated across.

    Parameters
    ----------
    df : pandas.DataFrame
        Forebay profiles as returned by :func:`load_profiles`.
    value : str, optional
        Column to average. Defaults to ``DO_mgL`` (mg/L).

    Returns
    -------
    pandas.DataFrame
        Rows are 1 m depth bins from 0 to the deepest observed bin, columns are the
        months present in the record, values are means in the units of ``value``.
    """
    d = df.dropna(subset=[value]).copy()
    d["zbin"] = np.floor(d.Depth_m).astype(int)
    piv = d.pivot_table(index="zbin", columns="month", values=value, aggfunc="mean")
    return piv.reindex(index=range(0, int(d.zbin.max())+1))


# Georgia water quality dissolved-oxygen criteria for the tailwater reaches, taken
# from the project Water Control Manuals. Etowah River below Allatoona Dam: daily
# average not less than 5.0 mg/L and not less than 4.0 mg/L at all times. Chattahoochee
# River below Buford Dam (secondary trout stream): daily average 6.0 mg/L and not less
# than 5.0 mg/L at all times.
DO_CRIT = {"Allatoona": {"inst": 4.0, "davg": 5.0},
           "Buford":    {"inst": 5.0, "davg": 6.0}}
TEMP_MAX_C = {"Allatoona": (90.0 - 32.0) * 5.0 / 9.0}  # Etowah: not to exceed 90 deg F

def do_criteria(ax: Axes, project: str, legend: bool = True,
                loc: str = "upper right", fs: float = 6.8) -> None:
    """Draw the daily-average and instantaneous-minimum DO criteria for a project.

    Both lines are reference lines in mg/L. On a panel of sub-daily observations the
    daily-average line is the standard the daily mean of those observations must meet,
    not a threshold each observation must clear, and its legend entry says so.
    """
    c = DO_CRIT[project]
    ax.axhline(c["davg"], color="red", ls="-", lw=0.9,
               label=f"{c['davg']:.0f} mg/L, daily-avg min")
    ax.axhline(c["inst"], color="red", ls="--", lw=0.9,
               label=f"{c['inst']:.0f} mg/L, min at all times")
    if legend:
        ax.legend(loc=loc, fontsize=fs)

# ---- Fig 1: representative late-summer stratification profiles ----
# One August or September cast per lake stands for the stratified-season condition.
# The complete set of casts, in this style, is rendered by forebay_profile_gallery.py,
# so an alternative cast can be checked before it is substituted here.
print("Figure 1: representative forebay profiles")
fig, axes = plt.subplots(1, 2, figsize=(7.6, 5.8), sharey=True, constrained_layout=True)
# Both panels share one depth axis, scaled to the deeper of the two casts.
zmax = 0.0
for ax, (lake, df) in zip(axes, prof.items()):
    d = select_representative_cast(df)
    cast = df[df.Date == d].sort_values("Depth_m")
    zmax = max(zmax, cast.Depth_m.max())
    ax.plot(cast.DO_mgL, cast.Depth_m, "-o", ms=3, color="#1f77b4")
    ax.set_xlabel("Dissolved oxygen (mg/L)", color="#1f77b4")
    ax.tick_params(axis="x", labelcolor="#1f77b4")
    ax.set_xlim(0, 14)
    axt = ax.twiny()
    axt.plot(cast.Temp_C, cast.Depth_m, "-s", ms=3, color="#d62728")
    axt.set_xlabel("Temperature (°C)", color="#d62728"); axt.set_xlim(0, 32)
    axt.tick_params(axis="x", labelcolor="#d62728"); axt.grid(False)
    # 4 mg/L is drawn as an orientation line for the in-pool profile. It is the Etowah
    # tailwater instantaneous minimum, not an in-pool standard, so it carries no label.
    ax.axvline(4.0, color="gray", ls="--", lw=0.8)
    # lake + date label placed clear of both top axes, inside the panel
    ax.text(0.5, 0.045, f"{lake}\n{pd.Timestamp(d).date()}", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=8.3,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9))
# surface at top, depth increasing downward (set once on the shared y-axis)
axes[0].set_ylim(np.ceil(zmax / 5) * 5, 0)
axes[0].set_ylabel("Depth (m)")
save(fig, "fig1_forebay_profiles")

# ---- Fig 2: forebay DO seasonal heatmap (month x depth) ----
print("Figure 2: forebay DO seasonal pattern")
fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.6), sharey=True)
for ax, (lake, df) in zip(axes, prof.items()):
    piv = monthly_depth_mean(df)
    im = ax.pcolormesh(piv.columns, piv.index, piv.values, cmap="RdYlBu",
                       vmin=0, vmax=12, shading="nearest")
    ax.set_xlabel("Month"); ax.set_title(lake, fontsize=8.5)
    ax.set_xticks(range(1,13)); ax.set_xticklabels(list("JFMAMJJASOND"))
# Surface at top, capped at 50 m. Four points on the 2005-04-13 Lanier cast reach 50-65 m,
# which is below the dam foundation and not physical; every other cast ends by about 47 m.
axes[0].set_ylim(50, 0)
axes[0].set_ylabel("Depth (m)")
cb = fig.colorbar(im, ax=axes, fraction=0.04, pad=0.02); cb.set_label("Mean DO (mg/L)")
fig.suptitle("Forebay dissolved-oxygen seasonal cycle (monthly mean by depth, full period of record)",
             fontsize=9.5)
save(fig, "fig2_forebay_do_seasonal")

# ---- Fig 2b: forebay temperature seasonal heatmap, companion to fig2 ----
# Same grid and depth axis as the oxygen figure, so the two can be read against each
# other: the thermocline in this figure sits at the depth where the oxygen figure turns
# anoxic. Temperature is reported on every cast at both stations, where oxygen is
# missing from one Lanier cast, so this grid rests on slightly more data.
#
# Scale: 5 to 30 deg C spans the monthly means at both lakes (7.4 to 28.8 at Lanier,
# 8.7 to 29.7 at Allatoona) and is held common to both panels so the colors compare.
# A warm sequential map matches the report convention of red for temperature, and
# reversing RdYlBu keeps this figure in the same color family as its oxygen companion.
print("Figure 2b: forebay temperature seasonal pattern")
fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.6), sharey=True)
for ax, (lake, df) in zip(axes, prof.items()):
    piv = monthly_depth_mean(df, value="Temp_C")
    im = ax.pcolormesh(piv.columns, piv.index, piv.values, cmap="RdYlBu_r",
                       vmin=5, vmax=30, shading="nearest")
    ax.set_xlabel("Month"); ax.set_title(lake, fontsize=8.5)
    ax.set_xticks(range(1,13)); ax.set_xticklabels(list("JFMAMJJASOND"))
axes[0].set_ylim(50, 0)
axes[0].set_ylabel("Depth (m)")
cb = fig.colorbar(im, ax=axes, fraction=0.04, pad=0.02)
cb.set_label("Mean temperature (°C)")
fig.suptitle("Forebay temperature seasonal cycle (monthly mean by depth, full period of record)",
             fontsize=9.5)
save(fig, "fig2b_forebay_temp_seasonal")

# ---- Fig 9: tailwater DO records (calibration windows) ----
print("Figure 9: tailwater DO records")
def tw(path: str, col: str) -> pd.DataFrame:
    """Read an archived USGS daily-value CSV as Date plus a value column named ``v``.

    Parameters
    ----------
    path : str
        Path relative to the repository root.
    col : str
        Full NWIS column name, ``X_<pcode>_00003`` for the daily mean: 00060 discharge
        (cfs), 00300 dissolved oxygen (mg/L).

    Returns
    -------
    pandas.DataFrame
        Columns ``Date`` (datetime64) and ``v`` (float, units of the parameter code),
        rows with a missing value dropped. Every archived file carries one row per date
        with no duplicates, so a later merge on ``Date`` is one-to-one.
    """
    df = pd.read_csv(os.path.join(REPO, path))
    df["Date"] = pd.to_datetime(df["Date"]); df["v"] = pd.to_numeric(df[col], errors="coerce")
    return df.dropna(subset=["v"])
allat_do = tw("Data/do_daily_02394000_2005_2007.csv", "X_00300_00003")
buford_do = tw("Data/Lanier/do_daily_02334430_2005_2026.csv", "X_00300_00003")


def daily_mean(df: pd.DataFrame, value: str = "DO_mgL_qa",
               out: str = "v") -> pd.DataFrame:
    """Calendar-day mean of a sub-hourly record, keyed on ``Date``.

    Parameters
    ----------
    df : pandas.DataFrame
        Sub-hourly record carrying a ``datetime`` column and ``value``.
    value : str, optional
        Column to average. Defaults to the quality-controlled oxygen column
        ``DO_mgL_qa`` (mg/L).
    out : str, optional
        Name given to the averaged column in the result.

    Returns
    -------
    pandas.DataFrame
        Columns ``Date`` (datetime64, midnight) and ``out``. Days with no observation
        are absent rather than NaN, so a plotted line would join across a gap unless
        :func:`break_gaps` is applied first.

    Notes
    -----
    No minimum-observation threshold is applied, so a day with a partial record still
    contributes a mean weighted by whatever hours it holds. On the Allatoona sonde that
    is 2 of 692 days; on the Buford monitor 8 of 1013. The diel-range statistics do
    apply MIN_OBS_PER_DAY, because a range from a few clustered samples is meaningless
    where a mean from them is merely imprecise.
    """
    s = (df.set_index("datetime")[value].resample("D").mean().dropna()
         .rename(out).reset_index().rename(columns={"datetime": "Date"}))
    return s


def diel_range(df: pd.DataFrame, value: str = "DO_mgL_qa",
               min_obs: int = MIN_OBS_PER_DAY) -> pd.Series:
    """Within-day range (daily maximum minus daily minimum) of a sub-hourly record.

    Parameters
    ----------
    df : pandas.DataFrame
        Sub-hourly record carrying a ``datetime`` column and ``value``.
    value : str, optional
        Column to range. Defaults to ``DO_mgL_qa`` (mg/L).
    min_obs : int, optional
        Days carrying fewer than this many observations are dropped, so the range is
        not taken from a handful of samples clustered in one part of the day.

    Returns
    -------
    pandas.Series
        Daily range in the units of ``value``, indexed by date, for qualifying days
        only. The number of days retained is reported wherever this is summarized.
    """
    d = df.set_index("datetime")[value]
    day_max = d.resample("D").max(); day_min = d.resample("D").min()
    day_n = d.resample("D").count()
    return (day_max - day_min)[day_n >= min_obs].dropna()


# USACE-SAM Allatoona tailrace sonde (30-min). Continuous rows flagged good become the
# daily-mean DO series; the good grab readings are held separately and drawn as points,
# because they are single instants months apart rather than part of the daily series.
sam = pd.read_csv(os.path.join(REPO, "Data", "allatoona_tailrace_sonde_2011_2019.csv"),
                  parse_dates=["datetime"])
sam_good = sam[sam["do_flag"] == "good"].copy()
sam_cont = sam_good[sam_good["sampling"] == "continuous"].copy()
sam_daily = daily_mean(sam_cont)
sam_grab = sam_good[sam_good["sampling"] == "grab"].dropna(subset=["DO_mgL_qa"])

# Corps Buford Dam tailrace monitor (2002-2008; hourly in 2002, 30-min in 2004, 15-min
# from 2005) - daily mean of quality-controlled DO. Ingest and QA:
# Scripts/python/buford_tailrace_monitor_ingest.py.
bufmon = pd.read_csv(os.path.join(REPO, "Data", "buford_tailrace_monitor_2002_2008.csv"),
                     parse_dates=["datetime"])
bufmon_daily = daily_mean(bufmon[bufmon["do_flag"] == "good"])


def break_gaps(df: pd.DataFrame, max_gap_days: int = 7) -> pd.DataFrame:
    """Insert missing rows so that plotted lines break across deployment gaps rather
    than interpolating across them.

    Parameters
    ----------
    df : pandas.DataFrame
        Columns ``Date`` and ``v``.
    max_gap_days : int, optional
        A step longer than this many days is treated as a deployment gap. Seven days
        separates the routine multi-day outages in these records from the seasonal
        redeployment breaks the figure should show as gaps.

    Returns
    -------
    pandas.DataFrame
        The input with one NaN-valued row added one day before each gap, sorted by
        date. The marker dates fall inside the gap and so cannot collide with a real
        observation; on the Buford monitor record this adds 19 rows to 1013.
    """
    d = df.sort_values("Date").reset_index(drop=True)
    gap = d["Date"].diff().dt.days > max_gap_days
    fill = pd.DataFrame({"Date": d.loc[gap, "Date"] - pd.Timedelta(days=1), "v": np.nan})
    return pd.concat([d, fill]).sort_values("Date").reset_index(drop=True)


bufmon_daily = break_gaps(bufmon_daily)

# The Corps Buford monitor and the USGS Allatoona gage overlap from 2005-02-10 to
# 2006-08-09, the only period in which both projects carry concurrent tailrace DO.
COMMON = (pd.Timestamp("2005-02-10"), pd.Timestamp("2006-08-09"))

fig, axes = plt.subplots(4, 1, figsize=(8.0, 8.6))
axes[0].plot(bufmon_daily.Date, bufmon_daily.v, lw=0.7, color="#2ca02c")
axes[0].set_title("Buford Dam tailrace (Corps monitor) — daily-mean DO from "
                  "sub-hourly record, 2002 and 2004--2008", fontsize=8.5)
axes[1].plot(allat_do.Date, allat_do.v, lw=0.7, color="#1f77b4")
axes[1].set_title("Allatoona Dam tailrace (USGS 02394000) — daily DO, 2005--2007",
                  fontsize=8.5)
axes[2].plot(sam_daily.Date, sam_daily.v, lw=0.7, color="#ff7f0e")
axes[2].scatter(pd.to_datetime(sam_grab.datetime), sam_grab.DO_mgL_qa, s=9,
                color="#8c564b", zorder=3, label="Grab readings (2019--2020)")
axes[2].set_title("Allatoona Dam tailrace (USACE-SAM sonde) — daily-mean DO from "
                  "30-min record, 2012--2013", fontsize=8.5)
axes[2].legend(loc="lower right", fontsize=6.5)
axes[3].plot(buford_do.Date, buford_do.v, lw=0.7, color="#2ca02c")
axes[3].set_title("Buford Dam tailrace (USGS 02334430) — daily DO, 2023 to present",
                  fontsize=8.5)
# Project-specific GA DO criteria: Buford (panels 1, 4) 6.0/5.0; Allatoona 5.0/4.0.
for ax, proj in zip(axes, ["Buford", "Allatoona", "Allatoona", "Buford"]):
    do_criteria(ax, proj, legend=False)
    ax.set_ylabel("Tailwater DO (mg/L)"); ax.set_ylim(0, 16)
# Mark the common window on the two panels that share it.
for ax in axes[:2]:
    ax.axvspan(*COMMON, color="gold", alpha=0.20, zorder=0)
    ax.text(COMMON[0], 15.2, " common window", fontsize=6.2, color="#7a6000", va="top")
axes[0].legend(loc="lower center", ncol=2, fontsize=6.2)
axes[1].legend(loc="lower center", ncol=2, fontsize=6.2)
fig.suptitle("Continuous tailwater DO records with the Georgia DO criteria "
             "(Allatoona 5.0/4.0 mg/L; Buford 6.0/5.0 mg/L)", fontsize=9.2, y=0.995)
fig.tight_layout(rect=[0,0,1,0.975]); save(fig, "fig9_tailwater_do_windows")

# ---- Fig 8: data-availability timeline (gap analysis) ----
print("Figure 8: data availability timeline")


def decimal_year(t: pd.Timestamp) -> float:
    """Timestamp as a decimal year, so a date can be drawn on a year axis."""
    start = pd.Timestamp(year=t.year, month=1, day=1)
    return t.year + (t - start).days / (366.0 if t.is_leap_year else 365.0)


def coverage_spans(dates, max_gap_days: float = 7.0) -> list:
    """Measured coverage of a record, as (start, end) decimal-year segments.

    A break longer than max_gap_days opens a new segment, so the bars show the gaps
    in the record rather than a single span from first observation to last. This is
    the figure the gap analysis rests on, so the segments are measured from the
    archived data rather than transcribed.

    Parameters
    ----------
    dates : sequence of datetime-like
        Observation timestamps, in any order.
    max_gap_days : float, optional
        Longest break drawn as continuous coverage.

    Returns
    -------
    list of (float, float)
        One (start, end) pair per segment, ascending. A segment covering a single
        instant is widened to one day so that it remains visible on a decade axis.
    """
    d = pd.Series(pd.to_datetime(pd.Series(dates).dropna().unique())).sort_values()
    if d.empty:
        return []
    breaks = d.diff() > pd.Timedelta(days=max_gap_days)
    out = []
    for _, seg in d.groupby(breaks.cumsum()):
        lo, hi = decimal_year(seg.iloc[0]), decimal_year(seg.iloc[-1])
        out.append((lo, max(hi, lo + 1.0 / 365.0)))
    return out


# Discharge and temperature spans are transcribed from the NWIS site record, which
# reaches back further than the slices archived here. Every dissolved-oxygen and
# profile span is measured from the archived data by coverage_spans, because those are
# the records whose gaps the analysis turns on. Profile rows carry the depth-resolved
# casts, which is what the study uses; the sparse pre-2000 casts at these stations
# report no depth and are not plotted.
prof_dates = {lake: df.Date for lake, df in prof.items()}
rows = [  # (project, label, [(start, end), ...], category)
 ("Allatoona","Discharge (02394000)",[(1938,2026)],"discharge"),
 ("Allatoona","Tailwater temp (02394000)",[(2005,2026)],"temperature"),
 ("Allatoona","Tailwater DO (USGS 02394000)",coverage_spans(allat_do.Date, 31),"do"),
 ("Allatoona","Tailwater DO (USACE-SAM sonde)",coverage_spans(sam_daily.Date, 31),"do"),
 ("Allatoona","Forebay profiles (GA EPD, depth-resolved)",
  coverage_spans(prof_dates["Lake Allatoona (Allatoona Dam forebay)"], 366),"profiles"),
 ("Lanier","Discharge (02334430)",[(1942,2026)],"discharge"),
 ("Lanier","Tailwater temp (02334430)",[(1975,2026)],"temperature"),
 ("Lanier","Tailwater DO (Corps monitor)",coverage_spans(bufmon.datetime),"do"),
 ("Lanier","Tailwater DO (02334430)",coverage_spans(buford_do.Date, 31),"do"),
 ("Lanier","Forebay profiles (GA EPD, depth-resolved)",
  coverage_spans(prof_dates["Lake Lanier (Buford Dam forebay)"], 366),"profiles"),
]
for _p, _label, _spans, _cat in rows:
    if len(_spans) > 1:
        print(f"  {_label}: {len(_spans)} segments")
CAT_COLOR = {"do": "#2ca02c", "discharge": "#4c72b0", "temperature": "#d62728",
             "profiles": "#9467bd"}
fig, ax = plt.subplots(figsize=(8.5, 4.4))
y = 0; yt=[]; yl=[]
for proj in ["Lanier","Allatoona"]:
    for p,label,spans,cat in [r for r in rows if r[0]==proj]:
        color = CAT_COLOR[cat]
        for s, e in spans:
            ax.barh(y, e-s, left=s, height=0.6, color=color, alpha=0.9)
        yt.append(y); yl.append(f"{label}"); y+=1
    y+=0.6
ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=7.5)
ax.set_xlim(1935, 2028); ax.set_xlabel("Year")
# The 2005-02 to 2006-08 overlap is the only period with concurrent tailrace DO at
# both projects; the remaining DO windows are single-project.
ax.axvspan(2005.11, 2006.60, color="gold", alpha=0.30)
ax.axvspan(2006.60, 2007.05, color="gold", alpha=0.14)
ax.axvspan(2012.0, 2014.0, color="gold", alpha=0.14)
ax.axvspan(2023.55, 2026, color="lightgreen", alpha=0.25)
ax.set_title("Periods of record by data category: the 2005--2006 overlap (gold) is the only\n"
             "period with concurrent tailrace DO at both projects", fontsize=9)
ax.legend(handles=[Patch(color="#2ca02c",label="Tailwater DO"),
                   Patch(color="#9467bd",label="Forebay T/DO profiles"),
                   Patch(color="#4c72b0",label="Discharge"),
                   Patch(color="#d62728",label="Temperature")],
          loc="upper left", fontsize=7.5)
ax.annotate("common\nwindow", xy=(2005.9, len(yt) - 0.2), ha="center", va="bottom",
            fontsize=6.4, color="#7a6000", annotation_clip=False)
save(fig, "fig8_data_timeline")

# ---- Summary statistics for the narrative ----
# Console only: nothing in this block reaches a statistics file, yet the report quotes
# it. The cast counts and the anoxic fraction appear in
# report/latex/src/sections/02_study_area.tex (169 and 144 casts; 17 percent of
# Jul-Sep depth observations below 2 mg/L at Lanier), so a change of format or
# precision here changes text written against it.
print("\n===== SUMMARY STATISTICS =====")
for lake, df in prof.items():
    late = df[df.month.isin(LATE_SUMMER)].dropna(subset=["DO_mgL"])
    ncast = late.Date.nunique()
    # anoxic fraction in late-summer casts (DO < 2 mg/L)
    anox = (late.DO_mgL < 2).mean()*100
    print(f"{lake}: {df.Date.nunique()} casts total; {ncast} Jul-Sep casts; "
          f"{anox:.0f}% of Jul-Sep depth-obs < 2 mg/L; maxZ {df.Depth_m.max():.0f} m")
# Each project is counted against its own instantaneous minimum: 4.0 mg/L on the Etowah
# below Allatoona, 5.0 mg/L on the Chattahoochee below Buford, which is a secondary
# trout stream. These are daily means counted against an instantaneous criterion, so
# each count is a lower bound on the true exceedance.
for name, df, proj in [("Allatoona tailrace 2005-07", allat_do, "Allatoona"),
                       ("Buford tailrace 2023-26", buford_do, "Buford")]:
    s = df[(df.Date.dt.month.isin(SUMMER))]
    inst = DO_CRIT[proj]["inst"]
    print(f"{name}: n={len(df)}; summer(JJAS) median DO={s.v.median():.1f} mg/L, "
          f"min={df.v.min():.1f} mg/L; days<{inst:.0f}mg/L overall={(df.v < inst).sum()} "
          f"({(df.v < inst).mean()*100:.0f}%)")

# ---- Fig 10: Allatoona tailrace sonde DO climatology and diel range (2012-2013) ----
print("\nFigure 10: Allatoona tailrace sonde DO climatology (USACE-SAM, 2012-2013)")
sc = sam_cont.copy()
sc["DO"] = sc["DO_mgL_qa"]
# One box per calendar month, pooling every 30-min observation of both years, so the
# distribution is observation-weighted and a month with fuller coverage carries more
# weight than one interrupted by a deployment gap.
by_month = [sc.loc[sc.month == m, "DO"].dropna().values for m in range(1, 13)]
# diel amplitude by month, over days carrying at least MIN_OBS_PER_DAY observations
amp = diel_range(sc, "DO")
amp_month = amp.groupby(amp.index.month).median().reindex(range(1, 13))

fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2))
ax = axes[0]
# Whiskers at the 5th and 95th percentiles, outliers suppressed: at 30-minute
# resolution the outlier cloud would obscure the boxes.
bp = ax.boxplot(by_month, positions=range(1, 13), widths=0.6, showfliers=False,
                whis=(5, 95), patch_artist=True, medianprops=dict(color="black"))
for patch in bp["boxes"]:
    patch.set(facecolor="#ff7f0e", alpha=0.55)
do_criteria(ax, "Allatoona", loc="upper right", fs=7.0)
ax.set_xticks(range(1, 13)); ax.set_xticklabels(list("JFMAMJJASOND"))
ax.set_ylabel("Tailrace DO (mg/L)"); ax.set_xlabel("Month"); ax.set_ylim(0, 14)
ax.set_title("Monthly DO distribution", fontsize=9)
ax = axes[1]
ax.bar(range(1, 13), amp_month.values, color="#1f77b4", alpha=0.75)
ax.set_xticks(range(1, 13)); ax.set_xticklabels(list("JFMAMJJASOND"))
ax.set_ylabel("Median diel DO range (mg/L)"); ax.set_xlabel("Month")
ax.set_title("Within-day DO variability", fontsize=9)
fig.suptitle("Allatoona tailrace DO from the USACE-SAM 30-min sonde record (2012--2013): "
             "a summer sag with the largest diel swings", fontsize=9.2)
fig.tight_layout(rect=[0, 0, 1, 0.95]); save(fig, "fig10_tailrace_do_climatology")

# ---- Write a machine-checked statistics file for the report narrative ----
STATS = os.path.join(REPO, "analysis", "tailrace_sonde_stats.txt")
summer = sc[sc.month.isin(SUMMER)]["DO"].dropna()
jas = sc[sc.month.isin(LATE_SUMMER)]["DO"].dropna()
# Temperature is taken from every valid row of the sonde record, continuous and grab,
# and is not restricted to do_flag == good: the 2011 oxygen sensor failed but its
# thermistor did not, which is why this line spans 2011-2020 while the DO lines above
# span 2012-2013 only.
tvalid = sam.dropna(subset=["T_C"])
tsum = tvalid[tvalid.datetime.dt.month.isin(SUMMER)]["T_C"]
with open(STATS, "w") as f:
    f.write("Allatoona tailrace USACE-SAM sonde (30-min) — statistics for the report\n")
    f.write("Source: Data/allatoona_tailrace_sonde_2011_2019.csv (do_flag == good)\n")
    # Provenance line. These two figures are fixed text, not recomputed on this run:
    # they come from the sonde ingest, which cross-validated the 2012 third-quarter
    # sonde temperature against USGS 02394000 daily temperature (see the header of
    # Scripts/python/tailrace_sonde_ingest.py). "MAD" here is the mean absolute
    # difference, not the median absolute deviation. tailwater_do_statistics.py runs
    # the same kind of check for the Buford monitor against USGS 02334430; no script
    # recomputes this one.
    f.write("Provenance confirmed vs USGS 02394000 (2012 JAS temp: MAD 0.18 C, r 0.995)\n\n")
    f.write(f"Continuous good-DO record: {sam_cont.datetime.min():%Y-%m-%d} to "
            f"{sam_cont.datetime.max():%Y-%m-%d}; {sam_cont.datetime.dt.date.nunique()} days; "
            f"n={len(sam_cont):,} sub-hourly obs; median interval 30 min\n")
    f.write(f"2011: temperature usable; DO instrument-compromised (excluded). "
            f"2019-09 to 2020-02: {len(sam_grab)} grab readings.\n\n")
    f.write(f"Annual DO median: 2012 = {sc[sc.year==2012].DO.median():.2f}, "
            f"2013 = {sc[sc.year==2013].DO.median():.2f} mg/L\n")
    f.write("Seasonal DO median (mg/L): "
            f"DJF {sc[sc.month.isin(SEASONS['DJF'])].DO.median():.2f}, "
            f"MAM {sc[sc.month.isin(SEASONS['MAM'])].DO.median():.2f}, "
            f"JJA {sc[sc.month.isin(SEASONS['JJA'])].DO.median():.2f}, "
            f"SON {sc[sc.month.isin(SEASONS['SON'])].DO.median():.2f}\n\n")
    f.write(f"Summer Jun-Sep DO: n={len(summer):,}, median={summer.median():.2f}, "
            f"mean={summer.mean():.2f}, min={summer.min():.2f} mg/L\n")
    f.write(f"  %hours <2 mg/L={100*(summer<2).mean():.1f}; "
            f"<4 mg/L={100*(summer<4).mean():.1f}; <5 mg/L={100*(summer<5).mean():.1f}\n")
    f.write(f"Late summer Jul-Sep DO: median={jas.median():.2f}, "
            f"%hours <4={100*(jas<4).mean():.1f}\n")
    amp_summer = amp[amp.index.month.isin(SUMMER)]
    f.write(f"Summer (Jun-Sep) diel DO range: median={amp_summer.median():.2f}, "
            f"mean={amp_summer.mean():.2f}, 90th pct={q(amp_summer,90):.2f} mg/L "
            f"(n={len(amp_summer)} days >=20 obs)\n\n")
    f.write(f"Tailrace temperature (2011-2020 valid): summer Jun-Sep median={tsum.median():.1f} C "
            f"(10-90 pct {q(tsum,10):.1f}-{q(tsum,90):.1f}), full range "
            f"{tvalid.T_C.min():.1f}-{tvalid.T_C.max():.1f} C\n")
print(f"  wrote {STATS}")

# ---- Fig 11: raw 30-min DO and temperature signal, 2012-2013, with a summer diel zoom ----
print("\nFigure 11: raw 30-min tailrace DO and temperature (USACE-SAM, 2012-2013)")
DO_BLUE = "#1f77b4"; T_RED = "#d62728"
sig = sam_cont.sort_values("datetime").set_index("datetime")
fig, axes = plt.subplots(3, 1, figsize=(8.4, 7.2))
# Panel A: full-record DO at native 30-min resolution, with the GA Allatoona DO criteria
axes[0].plot(sig.index, sig["DO_mgL_qa"], lw=0.3, color=DO_BLUE)
do_criteria(axes[0], "Allatoona", loc="upper right", fs=6.8)
axes[0].set_ylabel("Tailrace DO (mg/L)"); axes[0].set_ylim(0, 14)
axes[0].set_title("Dissolved oxygen at native 30-minute resolution, 2012--2013", fontsize=8.5)
# Panel B: full-record temperature at native 30-min resolution, with the 90 deg F maximum
axes[1].plot(sig.index, sig["T_C"], lw=0.3, color=T_RED)
axes[1].axhline(TEMP_MAX_C["Allatoona"], color="gray", ls="--", lw=0.9,
                label="90 °F (32.2 °C) maximum")
axes[1].set_ylabel("Tailrace temperature (°C)"); axes[1].set_ylim(0, 34)
axes[1].legend(loc="upper right", fontsize=6.8)
axes[1].set_title("Temperature at native 30-minute resolution, 2012--2013", fontsize=8.5)
# Panel C: a representative late-summer window showing the diel cycle (DO + temperature)
win = sig.loc["2012-08-01":"2012-08-22"]
axc = axes[2]
lnA = axc.plot(win.index, win["DO_mgL_qa"], lw=0.9, color=DO_BLUE, label="DO")
lnC = axc.axhline(DO_CRIT["Allatoona"]["inst"], color="red", ls="--", lw=0.9,
                  label="DO min, 4.0 mg/L at all times")
axc.set_ylabel("DO (mg/L)", color=DO_BLUE); axc.tick_params(axis="y", labelcolor=DO_BLUE)
axc.set_ylim(0, 8)
axt = axc.twinx()
lnB = axt.plot(win.index, win["T_C"], lw=0.9, color=T_RED, label="Temperature")
axt.set_ylabel("Temperature (°C)", color=T_RED); axt.tick_params(axis="y", labelcolor=T_RED)
axt.grid(False)
axc.legend(handles=[lnA[0], lnB[0], lnC], loc="upper right", fontsize=6.6, ncol=3)
axc.set_title("Detail: 1--21 August 2012, showing the diel (within-day) cycle", fontsize=8.5)
axc.set_xlabel("Date")
fig.suptitle("Allatoona tailrace at 30-minute resolution (USACE-SAM sonde): the seasonal cycle "
             "and the diel swings that daily means hide", fontsize=9.0)
fig.tight_layout(rect=[0, 0, 1, 0.96]); save(fig, "fig11_tailrace_subhourly")

# ---- Fig 12: DO drivers and measured saturation for the 2012-2013 window ----
print("Figure 12: tailrace DO drivers and measured saturation (2012-2013)")
allat_q_new = tw("Data/discharge_daily_02394000_2011_2014.csv", "X_00060_00003")
d = sam_cont.copy()
d["season"] = d.month.map({m: ("Jun-Sep" if m in SUMMER else "Oct-May") for m in range(1, 13)})
# Daily means for the discharge pairing. Discharge is published as a daily mean only, so
# DO and temperature are reduced to daily means to match it rather than the reverse.
dd = d.set_index("datetime")
daily = pd.DataFrame({"DO": dd["DO_mgL_qa"].resample("D").mean(),
                      "T": dd["T_C"].resample("D").mean()}).dropna().reset_index()
daily["Date"] = daily["datetime"].dt.normalize()
# Inner merge on Date. The discharge file spans 2011-2014 with one row per date and no
# duplicates, so this cannot duplicate a day, and it retains all 692 sonde days.
doq = pd.merge(daily, allat_q_new.rename(columns={"v": "Q"})[["Date", "Q"]], on="Date").dropna()
# Spearman on the merged frame, so DO and Q are the same day by construction. The two
# sub-daily correlations use the same rows of one frame, so DO and temperature are the
# same instant, and pandas drops any pair where either side is missing.
rho_tq = doq["DO"].corr(doq["Q"], method="spearman")
rho_tt = d["DO_mgL_qa"].corr(d["T_C"], method="spearman")
rho_tt_sum = d.loc[d.season=="Jun-Sep", "DO_mgL_qa"].corr(d.loc[d.season=="Jun-Sep", "T_C"], method="spearman")

fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.5))
CJ = {"Jun-Sep": "#d62728", "Oct-May": "#4c72b0"}
# (a) DO vs temperature, 30-min, colored by season
for s, g in d.groupby("season"):
    axes[0].scatter(g["T_C"], g["DO_mgL_qa"], s=1.5, c=CJ[s], alpha=0.15, label=s)
do_criteria(axes[0], "Allatoona", legend=False)
axes[0].set_xlabel("Temperature (°C)"); axes[0].set_ylabel("DO (mg/L)")
axes[0].set_title(f"DO vs temperature\nSpearman $\\rho$ = {rho_tt:.2f}", fontsize=8.5)
axes[0].set_ylim(0, 14)
lg = axes[0].legend(loc="upper right", fontsize=6.5, markerscale=4)
for lh in lg.legend_handles: lh.set_alpha(1.0)
# (b) daily DO vs daily mean discharge (log x), colored by season
doq["season"] = doq["datetime"].dt.month.map(
    {m: "Jun-Sep" for m in SUMMER}).fillna("Oct-May")
for s, g in doq.groupby("season"):
    axes[1].scatter(g["Q"], g["DO"], s=8, c=CJ[s], alpha=0.6, label=s)
axes[1].set_xscale("log")
do_criteria(axes[1], "Allatoona", legend=False)
axes[1].set_xlabel("Daily mean discharge (cfs)"); axes[1].set_ylabel("Daily mean DO (mg/L)")
axes[1].set_title(f"DO vs discharge\nSpearman $\\rho$ = {rho_tq:.2f}", fontsize=8.5)
axes[1].set_ylim(0, 14)
# (c) DO percent saturation as logged by the sonde, not recomputed from DO and
# temperature, so it carries the instrument's own barometric compensation
by_month_sat = [d.loc[d.month == m, "DO_pctsat"].dropna().values for m in range(1, 13)]
bp = axes[2].boxplot(by_month_sat, positions=range(1, 13), widths=0.6, showfliers=False,
                     whis=(5, 95), patch_artist=True, medianprops=dict(color="black"))
for patch in bp["boxes"]:
    patch.set(facecolor="#2ca02c", alpha=0.5)
axes[2].axhline(100.0, color="gray", ls=":", lw=0.8)   # equilibrium with the atmosphere
axes[2].set_xticks(range(1, 13)); axes[2].set_xticklabels(list("JFMAMJJASOND"))
axes[2].set_ylabel("Measured DO saturation (%)"); axes[2].set_xlabel("Month")
axes[2].set_title("Measured DO saturation", fontsize=8.5); axes[2].set_ylim(0, 120)
fig.suptitle("Allatoona tailrace DO drivers and directly measured saturation, USACE-SAM sonde "
             "2012--2013", fontsize=9.0)
fig.tight_layout(rect=[0, 0, 1, 0.94]); save(fig, "fig12_tailrace_do_drivers")

# ---- Append driver/saturation statistics to the stats file ----
sat_sum = d.loc[d.season == "Jun-Sep", "DO_pctsat"].dropna()
with open(STATS, "a") as f:
    f.write("\n-- Drivers and measured saturation (2012-2013 continuous, good DO) --\n")
    f.write(f"Measured DO percent saturation, summer Jun-Sep: median={sat_sum.median():.0f}%, "
            f"min={sat_sum.min():.0f}%, 5th pct={q(sat_sum,5):.0f}%\n")
    f.write(f"Spearman rho(DO, temperature), all data = {rho_tt:.2f}; summer = {rho_tt_sum:.2f}\n")
    f.write(f"Spearman rho(daily DO, daily discharge) = {rho_tq:.2f} "
            f"(n={len(doq)} paired days; discharge USGS 02394000)\n")
print("  appended driver/saturation stats")
