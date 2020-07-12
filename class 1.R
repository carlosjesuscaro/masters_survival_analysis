#Title:COVIDAnalysisbasedoncrowdsourceddata
#Objective:TODO
#Createdby:carloscaro
#Createdon:2020-07-12

library(tidyverse)
d_raw<-read.csv('covid19.csv', col_types=cols(id='c',
                               case_in='c',
                               age='d',
                               if_onset_approximated='l',
                               international_traveler='l',
                               domestic_traveler='l',
                               traveler='l',
                               `visitingWuhan`='l',
                               `fromWuhan`='l',
                               .default='c'))
# l -> logical, c -> character, d -> decimal

as_date <- function(x) as.Date(x, format = "%m/%d/%y")
d <-
  d_raw %>%
    mutate(reporting_date = as_date(reporting_date),
           hosp_visit_date = as_date(hosp_visit_date),
           exposure_start = as_date(exposure_start),
           exposure_end = as_date(exposure_end),
           symptom_onset = as_date(symptom_onset),
           death_status = death != "0",
           death_date = as.Date(ifelse(!death %in% c("0", "1"), as.Date(death, format = "%m/%d/%y", origin = "1970-01-01"), NA), origin = "1970-01-01"),
           gender = factor(gender, levels = c("female", "male")))