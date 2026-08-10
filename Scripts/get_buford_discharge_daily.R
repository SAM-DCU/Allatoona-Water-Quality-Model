# Daily mean discharge at the Chattahoochee River at Buford Dam (USGS 02334430),
# the Lake Lanier tailrace. Writes Data/Lanier/discharge_daily_02334430_2005_2026.csv.
# Run from the repository root.
#
# NWIS parameter 00060 (discharge), statistic 00003 (daily mean), delivered in cubic
# feet per second. This is the release record the tailrace water quality records are
# read against: Scripts/python/discharge_temp_analysis.py (fig4, fig5) and
# Scripts/python/tailwater_do_statistics.py pair it day for day with tailrace DO and
# temperature.
#
# The output directory is not interchangeable with Data/. The four Buford tailrace
# parameters are split across two directories, discharge (00060) and DO (00300) under
# Data/Lanier/ and temperature (00010) and specific conductance (00095) under Data/,
# and each file is written where the analysis scripts already read it from.
#
# The archived CSVs were pulled with Scripts/python/tailrace_and_stations.py, which
# issues the same query and writes the same file. Every downstream reader indexes the
# archived CSV by the column names Date, X_00060_00003, and X_00060_00003_cd.

library(dataRetrieval)
library(dplyr)

site <- "02334430"          # Chattahoochee River at Buford Dam (Lake Lanier tailrace)
# Discharge at this gage runs from 1942, so the 2005 start is a deliberate truncation
# to the study window shared with the Allatoona pulls. The end date is the day the
# study's NWIS and Water Quality Portal sources were queried; holding it fixed keeps a
# rerun reproducing the archived file rather than extending it.
start_date <- "2005-01-01"
end_date <- "2026-06-29"

q_dv <- readNWISdv(
  siteNumbers = site,
  parameterCd = "00060",   # discharge, cubic feet per second
  startDate = start_date,
  endDate = end_date,
  statCd = "00003"         # daily mean; the _00003 column suffix comes from this code
)

# Substitutes dataRetrieval's readable column names for the raw
# X_<parameter>_<statistic> names NWIS returns.
q_dv <- q_dv %>%
  renameNWISColumns()

write.csv(q_dv, "Data/Lanier/discharge_daily_02334430_2005_2026.csv", row.names = FALSE)
