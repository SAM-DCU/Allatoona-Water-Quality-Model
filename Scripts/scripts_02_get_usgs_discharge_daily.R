library(dataRetrieval)
library(dplyr)

site <- "02394000"
start_date <- "2005-01-01"
end_date <- "2007-12-31"

q_dv <- readNWISdv(
  siteNumbers = site,
  parameterCd = "00060",   # discharge
  startDate = start_date,
  endDate = end_date,
  statCd = "00003"         # daily mean
)

q_dv <- q_dv %>%
  renameNWISColumns()

write.csv(q_dv, "data/discharge_daily_02394000_2005_2007.csv", row.names = FALSE)
