library(dataRetrieval)
library(dplyr)

site <- "02394000"
start_date <- "2005-01-01"
end_date <- "2007-12-31"

do_dv <- readNWISdv(
  siteNumbers = site,
  parameterCd = "00300",   # dissolved oxygen
  startDate = start_date,
  endDate = end_date,
  statCd = "00003"         # daily mean
)

do_dv <- do_dv %>%
  renameNWISColumns()

write.csv(do_dv, "data/do_daily_02394000_2005_2007.csv", row.names = FALSE)