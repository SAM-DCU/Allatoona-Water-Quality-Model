# Daily mean dissolved oxygen at the Etowah River at Allatoona Dam (USGS 02394000), the
# Lake Allatoona tailrace. Writes Data/do_daily_02394000_2005_2007.csv.
# Run from the repository root.
#
# NWIS parameter 00300 (dissolved oxygen), statistic 00003 (daily mean), delivered in
# mg/L. This is the first of the two Allatoona tailrace DO records, and the one that
# overlaps the Corps Buford monitor on the 309-day common calibration window
# (2005-02-10 to 2006-08-09). Read by Scripts/python/plot_figures.py (fig9),
# Scripts/python/discharge_temp_analysis.py (fig5), and
# Scripts/python/tailwater_do_statistics.py. The second Allatoona window, 2012-2013,
# is the USACE-SAM sonde record in Data/allatoona_tailrace_sonde_2011_2019.csv, not a
# USGS pull.
#
# Allatoona files live in the repository-wide Data/ directory; only the Lanier project
# has a project directory of its own, and only two of its four Buford parameters are
# kept there.
#
# The archived CSVs were pulled with the Python route in Scripts/python/, which issues
# the same query and writes the same file. Every downstream reader indexes the archived
# CSV by the column names Date, X_00300_00003, and X_00300_00003_cd.

library(dataRetrieval)
library(dplyr)

site <- "02394000"          # Etowah River at Allatoona Dam (Lake Allatoona tailrace)
# The window brackets the daily DO record at this gage rather than being cut to it: the
# request is padded to the end of 2007, and the archived file holds 613 daily values
# from 2005-01-27 to 2007-01-05. analysis/data_inventory.md records that span as the
# full extent of daily DO at this gage.
start_date <- "2005-01-01"
end_date <- "2007-12-31"

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

write.csv(do_dv, "Data/do_daily_02394000_2005_2007.csv", row.names = FALSE)
