# Allatoona Water Quality Model

A U.S. Army Corps of Engineers (USACE) Mobile District (SAM) repository holding the data, retrieval and analysis scripts, and reference documents for a model of in-pool and downstream water quality at Allatoona Dam, focused initially on dissolved oxygen (DO), discharge, and reservoir elevation in the Etowah River tailwater.

That remains the aim of the work. The repository has since grown to cover Lake Sidney Lanier and Buford Dam as well, whose Chattahoochee River tailwater is the companion case, so the datasets and scripts here are organized for the two projects side by side.

The repository currently holds data, scripts, and references. It does not contain models, model results, or a written report.

## 1. Directory layout

| Directory | Contents |
|---|---|
| `Data/` | 32 data files: continuous tailwater records for both projects, in-pool forebay profiles, U.S. Geological Survey (USGS) daily values, earlier Mobile District files, and map inputs |
| `Scripts/` | 5 R retrieval scripts at the top level and 11 Python scripts in `Scripts/python/` |
| `references/` | Local copies of the Allatoona and Buford Water Control Manuals (WCMs) |

Analysis scripts write their figures, statistics files, and maps into `analysis/` and `report/` directories that they create on demand. Those outputs are not tracked in this repository.

## 2. Data

Continuous tailwater records, one per project, produced by the ingest scripts from the raw sources listed beneath them:

- `Data/buford_tailrace_monitor_2002_2008.csv`, from the USACE below-dam project
  monitor at Buford. Raw source: `Data/Lanier_2002-2010.xls`.
- `Data/allatoona_tailrace_sonde_2011_2019.csv`, from the USACE Mobile District
  Hydrolab sonde in the Allatoona tailrace. Raw source: the quarterly workbooks in
  `Data/Allatoona_from_TJ_2026-07-07/`.

In-pool forebay temperature and DO profiles from the Georgia Environmental Protection Division (GA EPD), retrieved through the Water Quality Portal:

- `Data/inpool_forebay_profiles_Allatoona_LK_14_4494.csv`, Lake Allatoona upstream from
  the dam.
- `Data/inpool_forebay_profiles_Lanier_LK_12_4028.csv`, Lake Sidney Lanier upstream of
  the Buford Dam forebay.

USGS daily values, Etowah River at Allatoona Dam (gage 02394000): discharge over two windows (`discharge_daily_02394000_2005_2007.csv` and `discharge_daily_02394000_2011_2014.csv`), dissolved oxygen (`do_daily_02394000_2005_2007.csv`), and water temperature (`temp_daily_02394000_2005_2026.csv`).

USGS daily values, Chattahoochee River at Buford Dam (gage 02334430): discharge and dissolved oxygen under `Data/Lanier/`, water temperature and specific conductance at the top of `Data/`. The split is deliberate. Each file is written where the analysis scripts already read it from, and the retrieval scripts document that.

Earlier Mobile District files, retained unchanged:

- `Data/Allatoona Elevation_ from Troy.xlsx`
- `Data/allatoona_elevation_daily_2005_2007.csv`
- `Data/allatoona_wqp_water_DO_temp_wqx3.csv`
- `Data/GAEPD_AllatoonaLake_nutrients_wide.csv`
- `Data/GAEPD_EtowahRiver_DO_Temp_wide.csv`

Map inputs: `Data/station_coordinates.csv` (forebay and tailrace station locations for both projects) and `Data/river_label_anchors.csv` (inflow and outflow anchor points used to label rivers on the station maps).

## 3. Scripts

The R scripts retrieve USGS daily values through the `dataRetrieval` package and are run from the repository root. The Python scripts handle ingest, quality control, statistics, figures, and station maps, and resolve the repository root from their own location, so they can be run from anywhere.

`Scripts/README.md` documents the run order and each script's inputs and outputs. Start there rather than from this page.

The Python scripts need a scientific Python stack: pandas, numpy, scipy, matplotlib, folium, selenium, and Pillow. `Scripts/README.md` names the specific conda environment the project uses. The station map script also drives headless Chrome for its PNG export.

## 4. Reference documents

Both WCMs are the source of the Georgia water quality criteria applied in the analysis.  Local copies are in `references/`.

- Allatoona Water Control Manual: https://www.sam.usace.army.mil/Portals/46/docs/planning_environmental/act/docs/ACR/Appendix%20A3%20-%20Allatoona%20WCM_April%202022.pdf?ver=EdlbhSSSwkuZlmysKB6VBg%3d%3d
- Lanier Water Control Manual: https://www.sam.usace.army.mil/Portals/46/docs/planning_environmental/acf/docs/ACF%20Buford_Final_Mar%202017.pdf?ver=2017-04-17-120546-503
