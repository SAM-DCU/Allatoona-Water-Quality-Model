# Daily mean water temperature at the Chattahoochee River at Buford Dam (USGS 02334430),
# the Lake Lanier tailrace. Writes Data/temp_daily_02334430_2005_2026.csv.
# Run from the repository root.
#
# NWIS parameter 00010 (water temperature), statistic 00003 (daily mean), delivered in
# degrees Celsius. The record is the thermal signature of hypolimnetic withdrawal and
# is read by Scripts/python/discharge_temp_analysis.py (fig4, fig5) and
# Scripts/python/tailwater_do_statistics.py. It is also the record the withdrawal-depth
# inversion is driven by, in the wider study. Unlike the DO record it is nearly
# complete over the request, 7,818 daily values, 2005-01-01 to 2026-06-29.
#
# The output directory is not interchangeable with Data/Lanier/. The four Buford
# tailrace parameters are split across two directories, discharge (00060) and DO
# (00300) under Data/Lanier/ and temperature (00010) and specific conductance (00095)
# under Data/, and each file is written where the analysis scripts already read it
# from.
#
# The archived CSVs were pulled with Scripts/python/tailrace_and_stations.py, which
# issues the same query and writes the same file. Every downstream reader indexes the
# archived CSV by the column names Date, X_00010_00003, and X_00010_00003_cd.

library(dataRetrieval)
library(dplyr)

site <- "02334430"          # Chattahoochee River at Buford Dam (Lake Lanier tailrace)
# Temperature at this gage runs from 1975, so the 2005 start is a deliberate truncation
# to the study window shared with the Allatoona pulls. The end date is the day the
# study's NWIS and Water Quality Portal sources were queried; holding it fixed keeps a
# rerun reproducing the archived file rather than extending it.
start_date <- "2005-01-01"
end_date <- "2026-06-29"

temp_dv <- readNWISdv(
  siteNumbers = site,
  parameterCd = "00010",   # water temperature, degrees Celsius
  startDate = start_date,
  endDate = end_date,
  statCd = "00003"         # daily mean; the _00003 column suffix comes from this code
)

# Substitutes dataRetrieval's readable column names for the raw
# X_<parameter>_<statistic> names NWIS returns.
temp_dv <- temp_dv %>%
  renameNWISColumns()

write.csv(temp_dv, "Data/temp_daily_02334430_2005_2026.csv", row.names = FALSE)
