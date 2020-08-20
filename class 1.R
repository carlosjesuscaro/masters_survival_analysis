 #Title:COVID Analysis based on crowd sourced data
#Objective:TODO
#Createdby:carloscaro
#Createdon:2020-07-12

library(tidyverse)

d_raw <- read.csv('covid19.csv', col_types=cols(id='c',
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

# Data transformation and date manipulation
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

# Binary outcome: alive/dead

## Sex impact
### Frequency tables and independence tests
# Absolute frequency:
with(d, table(death_status, gender))
# Relatibe frequency (%):
with(d, round(100 * prop.table(table(death_status, gender), 2), 1))
# Test of independence: Chi square (this is an asymptotic test)
with(d, chisq.test(table(death_status, gender)))
# Test of independence: Fisher test (this is an exact test)
with(d, fisher.test(table(death_status, gender)))

### Logistic regression: Death vs Gender
# glm: generic linear model
summary(glm(death_status ~ gender, data = d, family = "binomial"))
# The outcome on Pr(>z) is the Wald test
# An alternative:
confint(glm(death_status ~ gender, data = d, family = "binomial"))

# Logistic regression: Death vs age
summary(glm(death_status ~ age, data = d, family = "binomial"))
fit <- glm(death_status ~ I(age/10), data = d, family = "binomial")
summary(fit)
exp(coef(fit)[2])