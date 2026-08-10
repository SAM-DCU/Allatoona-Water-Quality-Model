# Daily mean dissolved oxygen at the Chattahoochee River at Buford Dam (USGS 02334430),
# the Lake Lanier tailrace. Writes Data/Lanier/do_daily_02334430_2005_2026.csv.
# Run from the repository root.
#
# NWIS parameter 00300 (dissolved oxygen), statistic 00003 (daily mean), delivered in
# mg/L. This is one of the two tailwater DO records the study calibrates against at
# Buford. It is read by Scripts/python/plot_figures.py (fig9),
# Scripts/python/discharge_temp_analysis.py (fig5), and
# Scripts/python/tailwater_do_statistics.py.
#
# The output directory is not interchangeable with Data/. The four Buford tailrace
# parameters are split across two directories, discharge (00060) and DO (00300) under
# Data/Lanier/ and temperature (00010) and specific conductance (00095) under Data/,
# and each file is written where the analysis scripts already read it from.
#
# The archived CSVs were pulled with Scripts/python/tailrace_and_stations.py, which
# issues the same query and writes the same file. Every downstream reader indexes the
# archived CSV by the column names Date, X_00300_00003, and X_00300_00003_cd.

library(dataRetrieval)
library(dplyr)

site <- "02334430"          # Chattahoochee River at Buford Dam (Lake Lanier tailrace)
# The window matches the discharge and temperature pulls at this gage so the three
# records line up. The end date is the day the study's NWIS and Water Quality Portal
# sources were queried; holding it fixed keeps a rerun reproducing the archived file
# rather than extending it.
start_date <- "2005-01-01"
end_date <- "2026-06-29"

# NOTE: continuous DO at the Buford Dam tailrace begins 2023-07-27 (none for 2005-2007).
# The archived file therefore holds 1,060 daily values, 2023-07-27 to 2026-06-29, out of
# a request spanning 2005 to 2026. Buford tailrace DO for 2002 and 2004-2008, which is
# what covers the study's 2005 to 2006 calibration window, comes instead from the Corps
# below-dam monitor in Data/buford_tailrace_monitor_2002_2008.csv.
do_dv <- readNWISdv(
  siteNumbers = site,
  parameterCd = "00300",   # dissolved oxygen, mg/L
  startDate = start_date,
  endDate = end_date,
  statCd = "00003"         # daily mean; the _00003 column suffix comes from this code
)

# Substitutes dataRetrieval's readable column names for the raw
# X_<parameter>_<statistic> names NWIS returns.
do_dv <- do_dv %>%
  renameNWISColumns()

write.csv(do_dv, "Data/Lanier/do_daily_02334430_2005_2026.csv", row.names = FALSE)
