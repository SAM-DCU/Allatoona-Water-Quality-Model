# Scripts/

Data retrieval, quality control, analysis, and figure generation for the Lanier and
Allatoona tailwater dissolved-oxygen study. Every derived number and figure in the
report comes from one of these scripts.

Python scripts live in `Scripts/python/` and run under the `clearwater` conda
environment (`/opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3`),
which carries pandas, numpy, scipy, matplotlib, folium, selenium, and Pillow. Each
resolves the repository root from its own location, so they can be run from anywhere.
The R scripts document the team's `dataRetrieval` workflow and are run from the
repository root; the archived CSVs were pulled with their Python equivalents.

## Order of operations

Retrieval and ingest first, then the analysis and figure scripts, which read the tidy
CSVs the ingest steps produce.

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3
SCRATCH=$(mktemp -d)                                        # holds the raw portal pulls

$PY Scripts/python/wqp_inpool.py $SCRATCH                   # in-pool results, both lakes
$PY Scripts/python/forebay_profiles.py $SCRATCH             # -> Data/inpool_forebay_profiles_*.csv
$PY Scripts/python/tailrace_and_stations.py $SCRATCH        # Buford tailrace CSVs; station ranking
$PY Scripts/python/qa_forebay_chars.py                      # characteristic and unit audit
$PY Scripts/python/tailrace_sonde_ingest.py                 # -> Data/allatoona_tailrace_sonde_2011_2019.csv
$PY Scripts/python/buford_tailrace_monitor_ingest.py        # -> Data/buford_tailrace_monitor_2002_2008.csv

$PY Scripts/python/plot_figures.py                          # fig1, fig2, fig8-fig12
$PY Scripts/python/discharge_temp_analysis.py               # fig4, fig5
$PY Scripts/python/forebay_profile_gallery.py               # every forebay cast as a PNG
$PY Scripts/python/tailwater_do_statistics.py               # all derived tailwater DO statistics
$PY Scripts/python/station_map.py                           # three station maps, HTML and PNG
```

## What each script does

| Script | Reads | Writes |
|---|---|---|
| `wqp_inpool.py` | Water Quality Portal | `inpool_<lake>.csv` in the given directory |
| `forebay_profiles.py` | `inpool_<lake>.csv` | `Data/inpool_forebay_profiles_*.csv` |
| `tailrace_and_stations.py` | USGS NWIS, WQP, `inpool_<lake>.csv` | `Data/Lanier/*_daily_02334430_*.csv`; station ranking to stdout |
| `qa_forebay_chars.py` | Water Quality Portal | characteristic and unit listing to stdout |
| `tailrace_sonde_ingest.py` | `Data/Allatoona_from_TJ_2026-07-07/*.xls[x]` | `Data/allatoona_tailrace_sonde_2011_2019.csv` |
| `buford_tailrace_monitor_ingest.py` | `Data/Lanier_2002-2010.xls` | `Data/buford_tailrace_monitor_2002_2008.csv` |
| `plot_figures.py` | forebay profiles, tailrace records | fig1, fig2, fig8-fig12; `analysis/tailrace_sonde_stats.txt` |
| `discharge_temp_analysis.py` | USGS daily values | fig4, fig5 |
| `forebay_profile_gallery.py` | forebay profiles | `analysis/figures/forebay_profiles/` (one PNG per cast, contact sheets, index) |
| `tailwater_do_statistics.py` | tailrace records, forebay profiles | `analysis/tailwater_do_statistics.txt` |
| `station_map.py` | `Data/station_coordinates.csv`, `Data/river_label_anchors.csv` | `analysis/maps/`; report copies in `report/latex/src/images/` |

Figures are written to `analysis/figures/` as PDF for LaTeX and PNG for preview, and
are numbered to their order in the report rather than to the order they are generated.
The numbering therefore has gaps here: fig3, fig6, and fig7 come from the
withdrawal-depth analysis, which is not part of this repository.
Copy a changed PDF into `report/latex/src/images/` before rebuilding the report.
`station_map.py` is the exception: it writes its report copies directly, because those
are PNG captures rather than vector figures.

## Conventions

- **Statistics files are owned by one script each.** `plot_figures.py` rewrites
  `analysis/tailrace_sonde_stats.txt` in full on every run and
  `tailwater_do_statistics.py` rewrites `analysis/tailwater_do_statistics.txt`, so
  nothing may be appended to either from elsewhere.
- **Numbers in the report trace to a statistics file.** If a value changes, regenerate
  the file that owns it and update the report from that file rather than editing the
  number in the LaTeX.
- **Quality control is applied at ingest, not at analysis time.** The tidy CSVs carry
  both the raw column and a `_qa` column plus a flag, so every downstream script sees
  the same screened record and the screening remains auditable.
- **Water quality criteria are the Georgia standards named in the project Water
  Control Manuals**: Etowah River below Allatoona Dam 5.0 mg/L daily average and
  4.0 mg/L at all times, 90 deg F maximum; Chattahoochee River below Buford Dam, a
  secondary trout stream, 6.0 and 5.0 mg/L.
- **Project colours are consistent across figures and maps**: Buford and Lake Lanier
  blue, Allatoona orange. In tailrace plots that carry both, DO is blue and
  temperature red.

`station_map.py` needs `folium`, `selenium`, and Pillow; its PNG export drives
headless Chrome, whose driver Selenium Manager resolves automatically.
