# Title     : "Nonparametric methods for Survival Analysis"
# Objective : TODO
# Created by: carloscaro
# Created on: 2020-07-16

  # One sample
library(survival)
library(tidyverse)

## Entering right-censored data in R

dat <- data.frame(ratID = paste0("rat", 1:5),
                  time = c(55, 50, 70, 120, 110),
                  status = c(0, 1, 1, 0, 1))

## Kaplan-Meyer estimator
fit.KM <- survfit(Surv(time, status) ~ 1, data = dat)
summary(fit.KM)

plot(fit.KM, mark.time = TRUE,
     main = "Kaplan-Meier estimator",
     ylab = "Survival probability",
     xlab = "time (seconds)")

# Question: what is the median survival time?
fit.KM

## Nelson-AAlen estimator
fit.NA <- survfit(Surv(time, status) ~ 1, data = dat, type = "fh")
summary(fit.NA)

fit.NA

###########################################################################################
###########################################################################################
## Case study: the Xelox trial

library(survival)
library(tidyverse)
library(asaur)
dat <- gastricXelox

# How many events, how many censored data points?
table(dat$delta)

summary(dat)

# How the Progress Free Survival times data looks like (ignoring censoring info)?
hist(dat$timeWeeks * 7 / 365.25)
# No conclusion should be driven from this histogram, any conclusion would be incorrect

### Check censoring coding
with(dat, Surv(timeWeeks, delta))

### Kaplan-Meyer estimator
fit.KM <- survfit(Surv(timeWeeks, delta) ~ 1, data = dat)
fit.KM
summary(fit.KM)

plot(fit.KM)

# Time in weeks might be cumbersome to read: we can re-express it in years
dat <- mutate(dat, timeYears = timeWeeks * 7 / 365.25)
fit.KM <- survfit(Surv(timeYears, delta) ~ 1, data = dat, conf.type = "log-log")
summary(fit.KM)
plot(fit.KM)

### Median survival
# Question: what is the median survival time?
fit.KM

# Note that the definition of censoring depends on what's the quantity of interest.
# If we're interested in measuring the follow-up time, delta is to be 'inverted':
dat <- mutate(dat, delta_followUp = 1 - delta)
fit.followUp <- survfit(Surv(timeYears, delta_followUp) ~ 1, data = dat, conf.type = "log-log")
fit.followUp

###########################################################################################
###########################################################################################

# Nonparametric comparison of two samples
## Entering right-censored data in R
dat <- data.frame(ratID = paste0("rat", 1:5),
                  time = c(55, 50, 70, 120, 110),
                  status = c(0, 1, 1, 0, 1),
                  group = c(0, 1, 0, 1, 1))

## The logrank test
fit.logrank <- survdiff(Surv(time, status) ~ group, data = dat)
fit.logrank

###########################################################################################
###########################################################################################

## Case study: the pancreatic dataset
library(asaur)

dat <- pancreatic
head(dat)

# * M: metastatic
# * LA: locally advanced

# This dataset requires some preprocessing before proper survival analysis.

# 1. parse 'onstudy', 'progression' and 'death' dates correctly
# 2. compute progression free survival times and overall survival times (this dataset has no censored data)

### step 1: parse dates
# Check the manual page of 'as.Date'
fmt <- "%m/%d/%Y"
dat <- mutate(dat,
  onstudy = as.Date(as.character(onstudy), format = fmt),
  progression = as.Date(as.character(progression), format = fmt),
  death = as.Date(as.character(death), format = fmt)
)
head(dat)


### step 2: compute survival times
dat <- mutate(dat,
  OS = difftime(death, onstudy, units = "days"),
  PFS = ifelse(!is.na(progression), difftime(progression, onstudy, units = "days"), OS)
)
# OS: overall survival
# PFS: progression free survival

# Note: OS and PFS are expressed in days. We want them in months:

dat <- mutate(dat,
  OS = as.numeric(OS) / 30.5,
  PFS = as.numeric(PFS) / 30.5
)

### compare PFS in the 2 disease groups
# As we have no censoring, we can produce use simple boxplots:
library(ggplot2)

ggplot(dat, aes(stage, PFS)) +
  geom_boxplot() +
  theme_bw()

# more generally, Kaplan-Meier estimates:

fit.KM <- survfit(Surv(PFS) ~ stage, data = dat, conf.type = "log-log")
plot(fit.KM, col = 1:2)

fit.KM

### The logrank test
survdiff(Surv(PFS) ~ stage, data = dat)

# What's the estimated probability of not experiencing a cancer progression for (at least) 1 year?
summary(fit.KM, time = 12)

# It is similar in the 2 groups, namely between 13% and 15%.
# Said otherwise, chances are high that the cancer is going to make a comeback
# within one year.

# Can you repeat the analysis above, this time for OS?


###########################################################################################
###########################################################################################


## Stratified logrank test: pharmacoSmoking dataset

### The data
dat <- pharmacoSmoking
head(dat)

summary(dat)

## Kaplan-Meyer estimator
fit.KM <- survfit(Surv(ttr,relapse)~grp, data=dat)
fit.KM
plot(fit.KM, col=1:2)
# white line is the combination therapy
# red line is the patch only therapy


# Question: do the 2 treatment group differ significantly in terms of survival to relapse?
# calculating the p-value with the Logrank test
# p-value is very small so it is very significant
survdiff(Surv(ttr, relapse) ~ grp, data = dat)

# Critique: the 2 groups have different age distribution, which might confound our results.
# Lets investigate:
with(dat, prop.table(table(grp, ageGroup2), 1))

with(dat, mosaicplot(table(grp, ageGroup2)))

### stratified logrank test
survdiff(Surv(ttr, relapse) ~ grp + strata(ageGroup2), data = dat)
# based on the p-value (0.8%), we are still confident that on the relevance of the
# experiment results. The strata command essentially performs the test for the 2 different age
# groups and then combines the results into only one p-value

### extra
fit.4 <- survfit(Surv(ttr, relapse) ~ grp + employment, data = dat)
fit.4

plot(fit.4, col = 1:6)
legend("topright", lty = 1, col = 1:6, legend = names(fit.4$strata))


# The 3 'combination' curves seem all higher than the 3 'patchOnly' curves. Lets make a stratified test:
survdiff(Surv(ttr, relapse) ~ grp + strata(employment), data = dat)



