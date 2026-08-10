# Daily mean discharge at the Etowah River at Allatoona Dam (USGS 02394000), the Lake
# Allatoona tailrace. Writes Data/discharge_daily_02394000_2005_2007.csv.
# Run from the repository root.
#
# NWIS parameter 00060 (discharge), statistic 00003 (daily mean), delivered in cubic
# feet per second. Pulled over the same window as the Allatoona tailrace DO record so
# the two pair day for day: Scripts/python/discharge_temp_analysis.py (fig4, fig5) and
# Scripts/python/tailwater_do_statistics.py read them together.
#
# Allatoona files live in the repository-wide Data/ directory; only the Lanier project
# has a project directory of its own, and only two of its four Buford parameters are
# kept there.
#
# The archived CSVs were pulled with the Python route in Scripts/python/, which issues
# the same query and writes the same file. Every downstream reader indexes the archived
# CSV by the column names Date, X_00060_00003, and X_00060_00003_cd.

library(dataRetrieval)
library(dplyr)

site <- "02394000"          # Etowah River at Allatoona Dam (Lake Allatoona tailrace)
# The window brackets the tailrace DO record at this gage, which runs 2005-01 to
# 2007-01, so the request is padded past the end of that record rather than cut to it.
# Discharge itself goes back to 1938. The concurrent discharge for the second Allatoona
# DO window, the 2012-2013 USACE-SAM sonde record, is a separate pull archived as
# Data/discharge_daily_02394000_2011_2014.csv.
start_date <- "2005-01-01"
end_date <- "2007-12-31"

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

write.csv(q_dv, "Data/discharge_daily_02394000_2005_2007.csv", row.names = FALSE)
