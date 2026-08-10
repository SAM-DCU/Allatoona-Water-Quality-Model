#!/usr/bin/env python3
"""Release regime and tailrace-temperature analysis for the Lanier/Allatoona DO study.

Characterizes the release regime at both dams and the thermal signature of
hypolimnetic withdrawal, then pairs tailrace DO with concurrent discharge and
temperature over each project's DO window. Rank correlations are printed for the
report narrative; the two figures are written to analysis/figures/ as PDF and PNG:

  fig4_release_temperature_regime  release flow-duration and monthly tailrace
                                   temperature, both projects
  fig5_do_vs_release_temperature   tailrace DO against discharge and against
                                   temperature, by project and season

Inputs are the archived USGS daily-value CSVs in Data/ and Data/Lanier/. Every value
here is a daily mean as published by NWIS (statistic code 00003): discharge in cfs,
temperature in deg C, dissolved oxygen in mg/L.

The numerical kernels (discharge_summary, temperature_summary, pair_records,
rank_correlations, exceedance) return plain values and frames, so any printed number
can be recomputed without drawing a figure; the print_* helpers only format them.

Unlike the other analysis scripts, this one archives nothing to analysis/*.txt: its
results reach the report through the console. See the note on
:func:`print_temperature_summary`.

Run with the clearwater conda env:
    /opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3 \
        Scripts/python/discharge_temp_analysis.py
"""
import os
from typing import Any

import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy.stats import spearmanr

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIG = os.path.join(REPO, "analysis", "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25,
                     "figure.dpi": 130, "savefig.bbox": "tight"})

# Cubic feet per second to cubic metres per second. The international foot is exactly
# 0.3048 m (NIST Handbook 44, Appendix C), so 1 cfs is exactly 0.028316846592 m3/s;
# the value below is that rounded to six significant figures, which is finer than the
# one decimal place these conversions are printed to.
CFS_TO_CMS = 0.0283168
SUMMER = [6, 7, 8, 9]  # JJAS, consistent with plot_figures.py
WINTER = [12, 1, 2]    # DJF, the meteorological winter


def load_dv(path: str, pcode: str) -> pd.DataFrame:
    """Load a USGS daily-values CSV (Date, X_<pcode>_00003, _cd).

    Parameters
    ----------
    path : str
        Absolute path to an archived NWIS daily-values export.
    pcode : str
        Five-digit NWIS parameter code: 00010 temperature (deg C), 00060 discharge
        (cfs), 00300 dissolved oxygen (mg/L). Statistic code 00003 is the daily mean.

    Returns
    -------
    pandas.DataFrame
        Columns ``Date`` (datetime64), ``v`` (float, units of the parameter code),
        ``month`` (1-based), ``summer`` (bool, month in SUMMER), sorted ascending by
        date with the index reset. Rows whose date or value will not parse are dropped;
        on the archived files that is none of them.

    Notes
    -----
    The file is read as ``str`` and quotes are stripped, because some archived exports
    are quoted and some are not. Each file carries one row per date with no duplicates,
    so :func:`pair_records` cannot duplicate a day.

    Dropping unparseable rows means consecutive rows are not necessarily consecutive
    days; :func:`discharge_summary` documents where that matters.
    """
    df = pd.read_csv(path, dtype=str)
    df["Date"] = pd.to_datetime(df["Date"].str.strip('"'), errors="coerce")
    vcol = f"X_{pcode}_00003"
    df["v"] = pd.to_numeric(df[vcol].astype(str).str.strip('"'), errors="coerce")
    df = df.dropna(subset=["Date", "v"]).sort_values("Date").reset_index(drop=True)
    df["month"] = df["Date"].dt.month
    df["summer"] = df["month"].isin(SUMMER)
    return df[["Date", "v", "month", "summer"]]


def save(fig: Figure, name: str) -> None:
    """Write a figure as PDF (for LaTeX) and PNG (preview) to analysis/figures/."""
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig)
    print(f"  wrote {name}.pdf/.png")


# ------------------------------------------------------------------ load
allat_q = load_dv(os.path.join(REPO, "Data", "discharge_daily_02394000_2005_2007.csv"), "00060")
allat_do = load_dv(os.path.join(REPO, "Data", "do_daily_02394000_2005_2007.csv"), "00300")
buf_q = load_dv(os.path.join(REPO, "Data", "Lanier", "discharge_daily_02334430_2005_2026.csv"), "00060")
buf_do = load_dv(os.path.join(REPO, "Data", "Lanier", "do_daily_02334430_2005_2026.csv"), "00300")
buf_t = load_dv(os.path.join(REPO, "Data", "temp_daily_02334430_2005_2026.csv"), "00010")
allat_t = load_dv(os.path.join(REPO, "Data", "temp_daily_02394000_2005_2026.csv"), "00010")


def discharge_summary(df: pd.DataFrame) -> dict[str, Any]:
    """Flow statistics describing the release regime at one project.

    Parameters
    ----------
    df : pandas.DataFrame
        Daily-mean discharge as returned by :func:`load_dv`, with ``v`` in cfs.

    Returns
    -------
    dict
        ``n`` (count of days), ``start`` and ``end`` (dates), ``min``, ``p10``,
        ``median``, ``mean``, ``p90``, ``max`` (all cfs), ``peaking_ratio`` (p90/p10,
        dimensionless), and ``d_day`` (median absolute change between consecutive rows,
        cfs).

    Notes
    -----
    ``d_day`` differences consecutive *rows*, which are consecutive days only where the
    record has no missing dates. The archived Buford discharge record has one step
    longer than a day in its 7,848 rows and the Allatoona record none, and restricting
    the difference to true one-day steps leaves the median at 170 cfs for Buford, so the
    distinction does not affect the reported value. It would matter for a record with
    real gaps.

    Percentiles use the numpy default linear interpolation between order statistics.
    """
    v = df["v"].values
    p10, p50, p90 = np.percentile(v, [10, 50, 90])
    dday = np.abs(np.diff(v))
    return {"n": len(v), "start": df.Date.min().date(), "end": df.Date.max().date(),
            "min": v.min(), "p10": p10, "median": p50, "mean": v.mean(),
            "p90": p90, "max": v.max(), "peaking_ratio": p90 / p10,
            "d_day": np.median(dday)}


def temperature_summary(df: pd.DataFrame) -> dict[str, Any]:
    """Tailrace temperature statistics, including the monthly-median seasonal cycle.

    Parameters
    ----------
    df : pandas.DataFrame
        Daily-mean temperature as returned by :func:`load_dv`, with ``v`` in deg C.

    Returns
    -------
    dict
        ``n``, ``start``, ``end``; ``min``, ``median``, ``mean``, ``max`` (deg C);
        ``mo_min`` and ``mo_max`` (the coldest and warmest monthly medians, deg C) with
        their months ``mo_min_month`` and ``mo_max_month``; ``amplitude`` (deg C, the
        difference between them); and the ``summer_mean``, ``summer_median`` and
        ``winter_mean`` (deg C).

    Notes
    -----
    The seasonal amplitude is the spread of monthly medians over the whole record, not
    a within-year range, so it is a climatological amplitude and is damped relative to
    any single year. This is the value the report quotes as the thermal signature of
    hypolimnetic withdrawal.

    Every month is pooled across all years of the record; months with uneven year
    coverage therefore carry different sample sizes into their medians.
    """
    v = df["v"].values
    mc = df.groupby("month")["v"].median()
    summer = df[df.summer]["v"]
    winter = df[df.month.isin(WINTER)]["v"]
    return {"n": len(v), "start": df.Date.min().date(), "end": df.Date.max().date(),
            "min": v.min(), "median": np.median(v), "mean": v.mean(), "max": v.max(),
            "mo_min": mc.min(), "mo_min_month": mc.idxmin(),
            "mo_max": mc.max(), "mo_max_month": mc.idxmax(),
            "amplitude": mc.max() - mc.min(),
            "summer_mean": summer.mean(), "summer_median": summer.median(),
            "winter_mean": winter.mean()}


def print_discharge_summary(name: str, df: pd.DataFrame) -> None:
    """Print :func:`discharge_summary` with units on every line."""
    s = discharge_summary(df)
    print(f"{name}: n={s['n']}, {s['start']}..{s['end']}")
    print(f"  Q cfs: min={s['min']:.0f} p10={s['p10']:.0f} median={s['median']:.0f} "
          f"mean={s['mean']:.0f} p90={s['p90']:.0f} max={s['max']:.0f}")
    print(f"  Q m3/s: median={s['median']*CFS_TO_CMS:.1f} p10={s['p10']*CFS_TO_CMS:.1f} "
          f"p90={s['p90']*CFS_TO_CMS:.1f}")
    print(f"  peaking: p90/p10={s['peaking_ratio']:.1f}x; median |dQ day-to-day|="
          f"{s['d_day']:.0f} cfs ({s['d_day']/s['median']*100:.0f}% of median)")


def print_temperature_summary(name: str, df: pd.DataFrame) -> None:
    """Print :func:`temperature_summary` with units on every line.

    The seasonal amplitude on the third line is quoted directly in the report
    (report/latex/src/sections/02_study_area.tex) and in analysis/data_inventory.md, so
    this console block is a report source even though it is not archived to one of the
    analysis/*.txt statistics files. Changing its wording or precision changes text
    written against it.
    """
    s = temperature_summary(df)
    print(f"{name}: n={s['n']}, {s['start']}..{s['end']}")
    print(f"  T C: min={s['min']:.1f} median={s['median']:.1f} mean={s['mean']:.1f} "
          f"max={s['max']:.1f}")
    print(f"  monthly-median range: {s['mo_min']:.1f} (mo {s['mo_min_month']}) .. "
          f"{s['mo_max']:.1f} (mo {s['mo_max_month']}); "
          f"seasonal amplitude={s['amplitude']:.1f} C")
    print(f"  summer(JJAS) mean={s['summer_mean']:.1f} (median={s['summer_median']:.1f}); "
          f"winter(DJF) mean={s['winter_mean']:.1f}")


# Smallest paired sample a rank correlation is reported for. Below about ten pairs the
# Spearman coefficient is dominated by individual points and its p value is unreliable.
MIN_PAIRS = 10


def pair_records(do: pd.DataFrame, other: pd.DataFrame) -> pd.DataFrame:
    """Inner-join a DO record with a concurrent driver record on the calendar date.

    Parameters
    ----------
    do : pandas.DataFrame
        Daily-mean DO as returned by :func:`load_dv`; ``v`` is renamed to ``do`` (mg/L).
    other : pandas.DataFrame
        Daily-mean driver as returned by :func:`load_dv`; ``v`` is renamed to ``x``,
        in cfs for discharge or deg C for temperature.

    Returns
    -------
    pandas.DataFrame
        Columns ``Date``, ``do``, ``x``, plus ``month`` and ``summer`` carried from the
        DO record. One row per date present in both records.

    Notes
    -----
    Both inputs carry one row per date with no duplicates, so the join is one-to-one:
    it can drop days that only one record covers, but it cannot duplicate a day. The
    retained count is printed as ``n_overlap`` next to the date range it spans, so the
    silent drop is always visible.

    The pairing is by date only, so ``do`` and ``x`` are the same calendar day by
    construction and a rank correlation over the result cannot mismatch its pairs.
    """
    return pd.merge(do.rename(columns={"v": "do"}),
                    other.rename(columns={"v": "x"})[["Date", "x"]],
                    on="Date", how="inner")


def rank_correlations(m: pd.DataFrame) -> list[dict[str, Any]]:
    """Spearman rank correlation of DO against the paired driver, by season.

    Parameters
    ----------
    m : pandas.DataFrame
        Paired records as returned by :func:`pair_records`.

    Returns
    -------
    list of dict
        One entry per subset (all days, summer JJAS, non-summer), each carrying
        ``label``, ``n``, and either ``too_few`` when fewer than MIN_PAIRS days remain
        or ``rho`` (dimensionless), ``p`` (two-sided) and ``do_median`` (mg/L).

    Notes
    -----
    Spearman is used rather than Pearson because the DO-discharge relation is monotonic
    but far from linear, and because rank correlation is insensitive to the heavy right
    tail of a peaking release record.
    """
    out = []
    for lab, sub in (("all", m), ("summer JJAS", m[m.summer]), ("non-summer", m[~m.summer])):
        if len(sub) < MIN_PAIRS:
            out.append({"label": lab, "n": len(sub), "too_few": True})
            continue
        rho, p = spearmanr(sub["x"], sub["do"])
        out.append({"label": lab, "n": len(sub), "too_few": False,
                    "rho": rho, "p": p, "do_median": sub["do"].median()})
    return out


def pair(name: str, do: pd.DataFrame, other: pd.DataFrame, olabel: str) -> pd.DataFrame:
    """Pair a DO record with a driver, print the rank correlations, return the pairs."""
    m = pair_records(do, other)
    print(f"{name}: n_overlap={len(m)}, {m.Date.min().date()}..{m.Date.max().date()}")
    for r in rank_correlations(m):
        if r["too_few"]:
            print(f"  {r['label']:12s}: n={r['n']} (too few)"); continue
        print(f"  {r['label']:12s}: n={r['n']:4d}  Spearman {olabel}~DO rho={r['rho']:+.2f} "
              f"(p={r['p']:.1e})  DO median={r['do_median']:.1f}")
    return m


print("\n===== DISCHARGE REGIME =====")
print_discharge_summary("Buford 02334430 (Lanier tailrace)", buf_q)
print_discharge_summary("Allatoona 02394000 tailrace", allat_q)

print("\n===== TAILRACE TEMPERATURE =====")
print_temperature_summary("Buford 02334430 temp", buf_t)
print_temperature_summary("Allatoona 02394000 temp", allat_t)

print("\n===== DO vs DISCHARGE (pairing over DO windows) =====")
buf_do_q = pair("Buford DO~Q (2023-)", buf_do, buf_q, "Q")
allat_do_q = pair("Allatoona DO~Q (2005-07)", allat_do, allat_q, "Q")

print("\n===== DO vs TEMPERATURE =====")
buf_do_t = pair("Buford DO~T (2023-)", buf_do, buf_t, "T")
allat_do_t = pair("Allatoona DO~T (2005-07)", allat_do, allat_t, "T")


# --------------------------------------------- fig4: release regime and temperature
def exceedance(v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Flow-duration curve: discharge sorted descending against plotting position.

    Parameters
    ----------
    v : numpy.ndarray
        Daily mean discharge in cfs. Order is irrelevant; the routine sorts.

    Returns
    -------
    p : numpy.ndarray
        Exceedance probability in percent of days, from the Weibull plotting position
        ``i / (n + 1)`` for rank i counted from the largest value.
    s : numpy.ndarray
        The same discharges in cfs, sorted descending, aligned with ``p``.

    Notes
    -----
    The Weibull position is unbiased for the exceedance probability of the ranked
    sample and never returns exactly 0 or 100 percent, which keeps both ends of the
    curve on a logarithmic discharge axis.
    """
    s = np.sort(v)[::-1]
    p = np.arange(1, len(s) + 1) / (len(s) + 1) * 100
    return p, s


fig, ax = plt.subplots(1, 2, figsize=(10, 4.0))
for df, lab, col in ((buf_q, "Buford (02334430), 2005-2026", "#4c72b0"),
                     (allat_q, "Allatoona (02394000), 2005-2007", "#dd8452")):
    p, s = exceedance(df["v"].values)
    ax[0].plot(p, s, color=col, lw=1.6, label=lab)
ax[0].set_yscale("log"); ax[0].set_xlabel("Exceedance (% of days)")
ax[0].set_ylabel("Daily mean discharge (cfs)")
ax[0].set_title("Release flow-duration (daily mean)"); ax[0].legend(fontsize=7.5, loc="upper right")

mo = np.arange(1, 13)
# Band is the 10th to 90th percentile of daily means within each month, pooled over the
# whole record, so it shows the year-to-year plus within-month spread, not a single year.
for tdf, lab, col in ((buf_t, "Buford (02334430)", "#4c72b0"),
                      (allat_t, "Allatoona (02394000)", "#dd8452")):
    tm = tdf.groupby("month")["v"]
    med = tm.median().reindex(mo)
    lo = tm.quantile(0.1).reindex(mo); hi = tm.quantile(0.9).reindex(mo)
    ax[1].fill_between(mo, lo, hi, color=col, alpha=0.18)
    ax[1].plot(mo, med, color=col, marker="o", ms=4, lw=1.6, label=lab)
ax[1].set_xticks(mo); ax[1].set_xticklabels(list("JFMAMJJASOND"))
ax[1].set_xlabel("Month"); ax[1].set_ylabel("Tailrace temperature (°C)")
ax[1].set_title("Tailrace temperature, monthly median and 10-90th pct (2005-2026)", fontsize=8.5)
ax[1].legend(fontsize=7.5, loc="upper left")
fig.suptitle("Release regime and tailrace thermal signature of hypolimnetic withdrawal", fontsize=10)
save(fig, "fig4_release_temperature_regime")


# ------------------------------------------------- fig5: DO against its two drivers
# Georgia DO criteria (from the project Water Control Manuals): daily-average minimum
# and instantaneous ("at all times") minimum. Etowah below Allatoona: 5.0 / 4.0 mg/L.
# Chattahoochee below Buford (secondary trout stream): 6.0 / 5.0 mg/L.
DO_CRIT = {"Allatoona": (5.0, 4.0), "Buford": (6.0, 5.0)}

def scatter_do(ax: Axes, m: pd.DataFrame, xcol: str, xlabel: str, title: str,
               project: str, xlog: bool = False) -> None:
    """Scatter paired daily-mean DO against a driver, split by season.

    Both criterion lines are drawn. The plotted DO is a daily mean, so the solid
    daily-average line is the directly comparable standard; the dashed instantaneous
    minimum is shown for context, and a daily mean above it does not establish that
    every instant in that day cleared it.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    m : pandas.DataFrame
        Paired records as returned by :func:`pair_records`.
    xcol : str
        Column of ``m`` to place on the x axis, normally ``x``.
    xlabel, title : str
        Axis label (carrying the driver's units) and panel title.
    project : str
        ``Allatoona`` or ``Buford``, selecting the criteria drawn.
    xlog : bool, optional
        Use a logarithmic x axis, wanted for discharge and not for temperature.
    """
    s = m[m.summer]; w = m[~m.summer]
    ax.scatter(w[xcol], w["do"], s=10, c="#4c72b0", alpha=0.5, label="Oct-May", edgecolors="none")
    ax.scatter(s[xcol], s["do"], s=12, c="#c44e52", alpha=0.6, label="Jun-Sep", edgecolors="none")
    davg, inst = DO_CRIT[project]
    ax.axhline(davg, ls="-", color="red", lw=1.0, label=f"{davg:.0f} mg/L, daily-avg min")
    ax.axhline(inst, ls="--", color="red", lw=1.0, label=f"{inst:.0f} mg/L, min at all times")
    if xlog: ax.set_xscale("log")
    ax.set_xlabel(xlabel); ax.set_ylabel("Tailrace DO (mg/L)"); ax.set_title(title, fontsize=9)
    ax.legend(fontsize=6.5, loc="lower right")


fig, ax = plt.subplots(2, 2, figsize=(9.6, 7.4))
scatter_do(ax[0, 0], allat_do_q, "x", "Daily mean discharge (cfs)",
           "Allatoona DO vs discharge, 2005-2007", "Allatoona", xlog=True)
scatter_do(ax[0, 1], buf_do_q, "x", "Daily mean discharge (cfs)",
           "Buford DO vs discharge, 2023-2026", "Buford", xlog=True)
scatter_do(ax[1, 0], allat_do_t, "x", "Tailrace temperature (°C)",
           "Allatoona DO vs temperature, 2005-2007", "Allatoona")
scatter_do(ax[1, 1], buf_do_t, "x", "Tailrace temperature (°C)",
           "Buford DO vs temperature, 2023-2026", "Buford")
fig.suptitle("Tailrace dissolved oxygen versus release and temperature", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.96])
save(fig, "fig5_do_vs_release_temperature")
