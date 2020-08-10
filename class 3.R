# Title     : "Cox regression"
# Objective : TODO
# Created by: carloscaro
# Created on: 2020-07-21

# A manually worked out, simple example: two groups
## Load libraries

library(tidyverse)
library(maxLik)
library(survival)

## Data definition
dat <- data.frame(ratID = paste0("rat", 1:5),
                  time = c(55, 50, 70, 120, 110),
                  failure = c(0, 1, 1, 0, 1),
                  group = c(0, 1, 0, 1, 1))

# Total number of failures D:
sum(dat$failure)

# For convenience, rename 'group' to 'x':
dat <- rename(dat, x = group)
dat

# We also define an auxiliary data.frame containing events only:
dat.events <- subset(dat, failure == 1)

## Partial log-likelihood function

# Lets define the partial (log-)likelihood function
pLogLik <- function(beta) {
  numerator <- with(dat.events, x * beta)
  denominator <- rep(NA_real_, length(numerator))
  for(j in seq_along(denominator)) {
    t_j <- dat.events[j, "time"]
    risk_set <- subset(dat, time >= t_j)
    theta_j <- with(risk_set, exp(x * beta))
    denominator[j] <- log(sum(theta_j))
  }
  return(sum(numerator - denominator))
}

debugonce(pLogLik)
pLogLik(3)

# We can plot it:
f <- Vectorize(pLogLik)
curve(f, from = -4, to = 4)

## Maximum partial-Likelihood estimation
fit.ML <- maxLik(pLogLik, start = c(beta = 0))
summary(fit.ML)
# The p-value is 0.698 which means we accept the null hypothesis:
# there isnt significant difference among the groups

# With the `coxph` function:
fit.cph <- coxph(Surv(time, failure) ~ x, data = dat)
confint(fit.cph)
summary(fit.cph)

# We can reproduce the Likelihood-ratio test:
LRT <- 2 * (fit.ML$maximum - pLogLik(0))
data.frame(LRT = LRT,
           pvalue = pchisq(LRT, df = 1, lower.tail = FALSE))

# The Wald test is already in the `maxLik` summary output.

######################################################################################
######################################################################################

# A manually worked out, simple example: one continuous covariate

dat <- data.frame(time = c(6, 7, 10, 15, 19, 25),
                  event = c(1, 0, 1, 1, 0, 1),
                  age = c(67, 62, 34, 41, 46, 28))

fit <- coxph(Surv(time, event) ~ age, data = dat)
summary(fit)

pred <- survfit(fit, newdata = data.frame(age = c(20, 40, 60)))
pred
plot(pred, col = 1:3)

# We might express age in decades:
dat <- mutate(dat, age_dec = age / 10)
summary(coxph(Surv(time, event) ~ age_dec, data = dat))

######################################################################################
######################################################################################

# Case study: the pharmacoSmoking dataset

## Load the data
library(asaur)
dat <- pharmacoSmoking
head(dat)

## Fit the Cox model
fit <- coxph(Surv(ttr, relapse) ~ grp + age + gender + priorAttempts, data = dat)
summary(fit)

# We can change the contrasts as we see fit:
dat <- mutate(dat, grp = relevel(grp, ref = "patchOnly"))
fit <- update(fit)
summary(fit)

dat <- mutate(dat, employment = relevel(employment, ref = "pt"))
fit <- coxph(Surv(ttr, relapse) ~ grp + employment, data = dat)
summary(fit)

######################################################################################
######################################################################################

# Case study: pharmakoSmoking, reloaded
# Nicotine addiction

# A. Model for TTR given: TRT, age, employment, gender, race
# - Table of estimates
# - with commentary

# B. Risk stratification: based on the model, identify subjects at high risk
# Top 10 subjects at risk

# C. Predictions
# - Median TTR
# - S ('6 months' covariates)

# D. Compare TTR efficacy in full time vs non full time emplotyment

# Loading the data
d_raw <- pharmacoSmoking
head(d_raw)
table(d_raw$employment, useNA = "always")

# Data preparation
dat <-
  mutate(d_raw,
    employment = ifelse(employment == "ft", "ft", "other"),
    grp = relevel(grp, ref = "patchOnly"),
    race = ifelse(race == "white", "white",
           ifelse(race == "black", "black",
                  "other")))
table(dat$employment)
table(dat$race)

## A. Model fit
fit <- coxph(Surv(ttr, relapse) ~ grp + age + gender + employment + race, data = dat)
summary(fit)
# this model is comparing against patch only

# Based on the outcome:
# the most significant variables are: grpcombination, age and employment
# For example, a person not working full time is 2x more likely to relapse than a
# person working full time

# Before exporting the outcome of summary(fit), we can use the following command
# for better formating:
broom::tidy(fit) %>%
  write_csv("coefficients_table.csv")

## B. Data segmentation
# This is the creation of a new ('fiction') data set amnd it will be used for the
# prediction based on the model
d_new <- select(dat, -ttr, -relapse)

# Running the model from Point A with the 'new' data and creating a new columm called
# risk score
d_segmented <-
  d_new %>%
  mutate(risk_score = predict(fit, newdata = d_new, type = "lp"))
head(d_segmented)

# Sorting the data based on the top 10 with the highest risk
d_segmented %>%
  arrange(desc(risk_score)) %>%
  head(10)

## C. Predicting median ttr and Surv(6 months | covariates)
# Getting only one sample
d0 <- d_new[1,]
p_S <- survfit(fit, newdata = d0)
summary(p_S)

plot(p_S)
p_S
# The TTR (time to relapse) is 21 days

summary(p_S, time=10)
# 65% chances to survive 10 days
summary(p_S, time=180)
# 15.4% chances to survive 180 days

# Making it as a function
make_individual_prediction <- function(grp, age, gender, employment, race) {
  ## LEFT AS AN EXERCISE
  tibble(median_ttr = ..., surv_6_months = ...)
}

# Homework - D
library(tidyverse)
data_ft <- filter(dat, employment == 'ft')
data_ot <- filter(dat, employment == 'other')
model_ft = coxph(Surv(ttr, relapse) ~ grp + age + gender + race, data = data_ft)
summary(model_ft)
model_ot = coxph(Surv(ttr, relapse) ~ grp + age + gender + race, data = data_ot)
summary(model_ot)

#p_ft <- survfit(model_emp, newdata = data_ft)
#p_ft
#p_ot <- survfit(model_emp, newdata = data_ot)
#p_ot


######################################################################################
######################################################################################

# Case study: the lung cancer dataset

## Load the data
library(survival)

dat <- lung
dat$sex <- factor(dat$sex)

## Nelson-AAlen estimators
pred.NA <- survfit(Surv(time, status) ~ sex, data = dat, type = "fh")
plot(pred.NA, col = 1:2)

## Cox regression: predictions
fit.cph <- coxph(Surv(time, status) ~ sex, data = dat)

pred.cph <- survfit(fit.cph, newdata = data.frame(sex = factor(1:2)),
                    type = "aalen")

plot(pred.cph, col = 1:2)

# How does the proportional hazards assumption hold?
plot(pred.NA, fun = "cloglog", col = 1:2)

plot(pred.cph, fun = "cloglog", col = 1:2)