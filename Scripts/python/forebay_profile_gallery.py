#!/usr/bin/env python3
"""Every forebay temperature/DO profile, plotted one cast per PNG.

The report shows a single representative late-summer cast per lake (fig1, produced by
plot_figures.py). This script renders the complete set behind that figure: one PNG for
each GA EPD cast at the two near-dam stations, in the same style as fig1 (DO in blue on
the lower axis, temperature in red on the upper axis, depth increasing downward, and the
4 mg/L reference line), so that any cast can be substituted for the one in the report.

Outputs, under analysis/figures/forebay_profiles/:
  Lanier/lanier_forebay_YYYY-MM-DD.png       -- one file per cast
  Allatoona/allatoona_forebay_YYYY-MM-DD.png
  contact_sheets/<lake>_YYYY.png             -- all of that year's casts on one sheet
  profile_index.csv                          -- one row per cast (date, depth, surface and
                                                bottom values, file name) for cast selection

Axis limits are held fixed within a lake so casts are directly comparable across dates.

Units: depth in m below the water surface, temperature in deg C, dissolved oxygen in
mg/L. The index file names each measured column with its unit (z_max_m, T_surface_C,
DO_bottom_mgL, and so on).

Run with the clearwater conda env:
    /opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3 \
        Scripts/python/forebay_profile_gallery.py
"""
import os
from typing import Any

import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, "analysis", "figures", "forebay_profiles")
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25,
                     "figure.dpi": 150, "savefig.bbox": "tight"})

DO_BLUE, T_RED = "#1f77b4", "#d62728"      # matched to fig1 and fig11
DO_LIM, T_LIM = (0, 14), (0, 32)           # mg/L and deg C, fixed across all casts, as in fig1
# 4 mg/L is an orientation line for reading the in-pool profile. It is the Etowah
# tailwater instantaneous minimum from the Allatoona Water Control Manual, not an
# in-pool standard, which is why it carries no label here or in fig1.
DO_REF = 4.0                               # mg/L

# Depth axis floor per lake, in m. Casts deeper than the floor extend the axis to the
# next multiple of 5 m rather than being clipped. On the present record that happens for
# one cast only, Lanier 2005-04-13, which reaches 65 m; the next deepest Lanier cast ends
# at 47.2 m and no Allatoona cast passes its 45 m floor, so the axis is otherwise fixed.
FOREBAY = {
    "Lanier": dict(lake="Lake Lanier (Buford Dam forebay)", station="LK_12_4028",
                   csv="Data/inpool_forebay_profiles_Lanier_LK_12_4028.csv", zlim=50.0),
    "Allatoona": dict(lake="Lake Allatoona (Allatoona Dam forebay)", station="LK_14_4494",
                      csv="Data/inpool_forebay_profiles_Allatoona_LK_14_4494.csv", zlim=45.0),
}


def load_profiles(path: str) -> pd.DataFrame:
    """Read a forebay profile CSV, sorted by cast date then increasing depth.

    Parameters
    ----------
    path : str
        Path relative to the repository root.

    Returns
    -------
    pandas.DataFrame
        Columns ``Date`` (datetime64), ``Depth_m`` (m below surface), ``Temp_C``
        (deg C), ``DO_mgL`` (mg/L), sorted by ``Date`` then ``Depth_m``.

    Notes
    -----
    The depth sort is what makes "first row of a cast" mean the shallowest sample in
    :func:`cast_stats`. Non-numeric temperature and oxygen entries become NaN rather
    than raising, so a cast that reports only one of the two still plots the other.
    """
    df = pd.read_csv(os.path.join(REPO, path), parse_dates=["Date"])
    for c in ("Temp_C", "DO_mgL"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.sort_values(["Date", "Depth_m"])


def draw_profile(ax: Axes, cast: pd.DataFrame, zlim: float,
                 labels: bool = True, ms: float = 3) -> Axes:
    """Draw one cast: DO on ax (blue, lower axis), temperature on a twin (red, upper axis).

    Returns the twin axis so callers can adjust its tick labels on a contact sheet.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes for the DO trace; its x axis is fixed to DO_LIM in mg/L.
    cast : pandas.DataFrame
        One cast, already sorted by increasing depth.
    zlim : float
        Depth axis floor in m. A cast reaching deeper extends the axis to the next
        multiple of 5 m so that no observation is clipped out of the figure.
    labels : bool, optional
        Draw the axis labels. Suppressed on contact sheets, where the suptitle carries
        the variables and their ranges instead.
    ms : float, optional
        Marker size in points.

    Notes
    -----
    Temperature and DO are dropped independently, so a depth at which only one was
    recorded still contributes to the other trace. Missing depths inside a cast leave
    the line joined across them, which is the usual reading of a profile plot.
    """
    do = cast.dropna(subset=["DO_mgL"])
    ax.plot(do.DO_mgL, do.Depth_m, "-o", ms=ms, lw=1.0, color=DO_BLUE)
    ax.set_xlim(*DO_LIM)
    ax.tick_params(axis="x", labelcolor=DO_BLUE)
    ax.axvline(DO_REF, color="gray", ls="--", lw=0.8)

    axt = ax.twiny()
    t = cast.dropna(subset=["Temp_C"])
    axt.plot(t.Temp_C, t.Depth_m, "-s", ms=ms, lw=1.0, color=T_RED)
    axt.set_xlim(*T_LIM)
    axt.tick_params(axis="x", labelcolor=T_RED)
    axt.grid(False)

    ax.set_ylim(max(zlim, np.ceil(cast.Depth_m.max() / 5) * 5), 0)
    if labels:
        ax.set_xlabel("Dissolved oxygen (mg/L)", color=DO_BLUE)
        ax.set_ylabel("Depth (m)")
        axt.set_xlabel("Temperature (°C)", color=T_RED)
    return axt


def cast_stats(cast: pd.DataFrame) -> dict[str, Any]:
    """Surface and bottom values used for the index file (surface = shallowest sample).

    Parameters
    ----------
    cast : pandas.DataFrame
        One cast, sorted by increasing depth as :func:`load_profiles` leaves it.

    Returns
    -------
    dict
        ``n_points`` (rows in the cast, including any that carry only one variable),
        ``z_max_m`` (m), ``T_surface_C`` and ``T_bottom_C`` (deg C), ``DO_surface_mgL``
        and ``DO_bottom_mgL`` (mg/L), ``DO_min_mgL`` (mg/L).

    Notes
    -----
    Surface and bottom are the shallowest and deepest samples that reported the
    variable in question, not a fixed depth, and temperature and oxygen are taken
    independently. A cast whose deepest DO sample sits above its deepest temperature
    sample therefore reports the two bottom values from different depths.

    ``n_points`` counts every row of the cast. plot_figures.py selects fig1's cast on
    the count of rows complete in both variables, so the selection this script prints
    at the end can in principle differ from the figure's. On the present record they
    agree: 2005-08-10 at Lanier, 2025-08-18 at Allatoona.
    """
    do = cast.dropna(subset=["DO_mgL"])
    t = cast.dropna(subset=["Temp_C"])

    def first_or_nan(d: pd.DataFrame, col: str, i: int) -> float:
        return float(d[col].iloc[i]) if len(d) else np.nan

    return dict(n_points=len(cast), z_max_m=float(cast.Depth_m.max()),
                T_surface_C=first_or_nan(t, "Temp_C", 0),
                T_bottom_C=first_or_nan(t, "Temp_C", -1),
                DO_surface_mgL=first_or_nan(do, "DO_mgL", 0),
                DO_bottom_mgL=first_or_nan(do, "DO_mgL", -1),
                DO_min_mgL=(float(do.DO_mgL.min()) if len(do) else np.nan))


index_rows = []
for key, meta in FOREBAY.items():
    df = load_profiles(meta["csv"])
    lake_dir = os.path.join(OUT, key)
    os.makedirs(lake_dir, exist_ok=True)
    dates = sorted(df.Date.unique())
    print(f"{meta['lake']} ({meta['station']}): {len(dates)} casts, "
          f"{pd.Timestamp(dates[0]):%Y-%m-%d} to {pd.Timestamp(dates[-1]):%Y-%m-%d}")

    # ---- one PNG per cast ----
    # File names are the cast date, so the output paths are deterministic and a rerun
    # overwrites in place rather than accumulating.
    for d in dates:
        cast = df[df.Date == d]
        day = pd.Timestamp(d)
        fig, ax = plt.subplots(figsize=(4.2, 5.8), constrained_layout=True)
        draw_profile(ax, cast, meta["zlim"])
        ax.text(0.5, 0.045, f"{meta['lake']}\n{meta['station']}  {day:%Y-%m-%d}",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=8.3,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9))
        name = f"{key.lower()}_forebay_{day:%Y-%m-%d}.png"
        fig.savefig(os.path.join(lake_dir, name))
        plt.close(fig)
        index_rows.append(dict(lake=key, station=meta["station"], date=f"{day:%Y-%m-%d}",
                               month=day.month, year=day.year, **cast_stats(cast),
                               file=os.path.join(key, name)))

    # ---- contact sheet per year, for scanning the record quickly ----
    sheet_dir = os.path.join(OUT, "contact_sheets")
    os.makedirs(sheet_dir, exist_ok=True)
    for year, group in pd.Series(dates).groupby(pd.DatetimeIndex(dates).year):
        cols = min(4, len(group))
        rows = int(np.ceil(len(group) / cols))
        fig, axes = plt.subplots(rows, cols, figsize=(2.5 * cols, 3.4 * rows),
                                 squeeze=False, constrained_layout=True)
        flat = axes.ravel()
        for ax, d in zip(flat, group):
            cast = df[df.Date == d]
            axt = draw_profile(ax, cast, meta["zlim"], labels=False, ms=1.8)
            ax.set_title(f"{pd.Timestamp(d):%Y-%m-%d}", fontsize=8)
            axt.set_xticklabels([])
            ax.tick_params(labelsize=7)
        # Trailing cells of the last row are hidden rather than left as empty axes.
        for ax in flat[len(group):]:
            ax.set_visible(False)
        for r in range(rows):
            axes[r][0].set_ylabel("Depth (m)", fontsize=8)
        fig.suptitle(f"{meta['lake']} — forebay profiles, {year}. "
                     f"DO (blue, 0--14 mg/L), temperature (red, 0--32 °C)", fontsize=9.5)
        fig.savefig(os.path.join(sheet_dir, f"{key.lower()}_{year}.png"))
        plt.close(fig)
    print(f"  wrote {len(dates)} cast PNGs to {os.path.relpath(lake_dir, REPO)}/ "
          f"and {pd.DatetimeIndex(dates).year.nunique()} contact sheets")

# Sorted by lake then date so the file is stable across runs and diffs cleanly.
idx = pd.DataFrame(index_rows).sort_values(["lake", "date"])
idx_path = os.path.join(OUT, "profile_index.csv")
idx.to_csv(idx_path, index=False, float_format="%.2f")
print(f"\nwrote {os.path.relpath(idx_path, REPO)} ({len(idx)} casts)")

# The report's fig1 selects the most complete late-summer (August-September) cast per lake;
# list the current selection so the gallery and the report figure stay traceable to each other.
# See cast_stats on why this count and plot_figures.py's can differ in principle.
late = idx[idx.month.isin([8, 9])]
for key in FOREBAY:
    sel = late[late.lake == key].sort_values(["n_points", "z_max_m"])
    if len(sel):
        r = sel.iloc[-1]
        print(f"fig1 late-summer selection, {key}: {r.date} "
              f"({int(r.n_points)} points, {r.z_max_m:.0f} m)")
