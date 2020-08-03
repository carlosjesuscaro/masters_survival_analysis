# Title     : "Cox Model Building and Diagnostics"
# Objective : TODO
# Created by: carloscaro
# Created on: 2020-08-02

  # Model building
  ## Load the data

library(tidyverse)
library(survival)
library(asaur)

dat <- pharmacoSmoking

## The 4 candidate models
M0 <- coxph(Surv(ttr, relapse) ~ 1, data = dat)
MA <- coxph(Surv(ttr, relapse) ~ ageGroup4, data = dat)
MB <- coxph(Surv(ttr, relapse) ~ employment, data = dat)
MC <- coxph(Surv(ttr, relapse) ~ ageGroup4 + employment, data = dat)

summary(MC)

d <- mutate(dat, employment = ifelse(employment == "ft", "ft", "other"))
M_additive <- coxph(Surv(ttr, relapse) ~ grp + employment, data = d)
M_int <- coxph(Surv(ttr, relapse) ~ grp * employment, data = d)
summary(M_int)
anova(M_additive, M_int)

d$race <- relevel(d$race, ref = "other")
fit <- coxph(Surv(ttr, relapse) ~ grp + employment + gender + race + age, data = d)
summary(fit)
fit_noRaceNoGnd <- coxph(Surv(ttr, relapse) ~ grp + employment + age, data = d)
anova(fit_noRaceNoGnd, fit)

## Comparing nested models: LRT
anova(MA, MC)

## Comparing non-nested models: AIC
fits <- list(MA = MA, MB = MB, MC = MC)
sapply(fits, AIC)

sapply(list(ageOnly = MA, emplOnly = MB, full = MC), AIC)

## Automatic model selection based on AIC
Mfull <- coxph(Surv(ttr, relapse) ~ grp + gender + race +
                 employment + yearsSmoking + levelSmoking +
                 ageGroup4 + priorAttempts + longestNoSmoke,
               data = dat)

MAIC <- step(Mfull)

summary(MAIC)

## Predictive power: concordance index
summary(MA)

summary(MAIC)

## Predictive power: AUC

library(survivalROC)
data(mayo)

head(mayo)

plot(survfit(Surv(time / 365.25, censor) ~ 1, data = mayo))

ROC.4 <- survivalROC(Stime = mayo$time,
                     status = mayo$censor,
                     marker = mayo$mayoscore4,
                     predict.time = 365.25 * 5,
                     method="KM")
ROC.5 <- survivalROC(Stime = mayo$time,
                     status = mayo$censor,
                     marker = mayo$mayoscore5,
                     predict.time = 365.25 * 5,
                     method = "KM")

ROC <- list(mayo4 = ROC.4, mayo5 = ROC.5)
map_dbl(ROC, "AUC")

dfl <- map(ROC, ~ with(., tibble(cutoff = cut.values, FP, TP)))
for(nm in names(dfl)) {
  dfl[[ nm ]]$marker <- nm
}
dat <- do.call(rbind, dfl)

dat

ggplot(dat, aes(FP, TP, color = marker)) +
  geom_line() +
  theme_bw(base_size = 9)


cutoff <- min(filter(dat, marker == "mayo5", FP <= 0.1)$cutoff)


mayo$prediction <-
  ifelse(mayo$mayoscore5 <= cutoff,
         "low_risk", "high_risk")

plot(survfit(Surv(time/365, censor) ~ prediction, data = mayo),
     col = c("red", "blue"))


# Model diagnostics
## Martingale residuals

library(survival)
library(asaur) ## dataset

data(pharmacoSmoking)
dat <- pharmacoSmoking

fit <- coxph(Surv(ttr, relapse) ~ grp + age + employment, data = dat)
dat$residual <- residuals(fit, type = "martingale")


par(mfrow = c(1, 3), mar = c(4.2, 2, 2, 2))
with(dat, {

  plot(age, residual)
  lines(lowess(age, residual), lwd = 2)

  plot(residual ~ grp)

  plot(residual ~ employment)

})

dfbetas <- residuals(fit, type = 'dfbetas')
dat$dfbetas <- sqrt(rowSums(dfbetas^2))

plot(dat$dfbetas, type = 'h')
abline(h = 0)


## Proportionality of hazards

# Pancreatic cancer dataset

library(survival)
library(asaur) ## dataset
library(plyr)
library(ggplot2)

fmt <- "%m/%d/%Y"
dat <- as.tibble(pancreatic) %>%
  mutate(
  onstudy = as.Date(as.character(onstudy), format = fmt),
  progression = as.Date(as.character(progression), format = fmt),
  death = as.Date(as.character(death), format = fmt),
  OS = death - onstudy,
  PFS = ifelse(is.na(progression), OS, pmin(progression - onstudy, OS))) %>%
  mutate(
  PFS = Surv(as.numeric(PFS / 30.5)),
  OS = Surv(as.numeric(OS / 30.5))
  )
dat

fit <- coxph(PFS ~ stage, data = dat)
summary(fit)

fit.KM <- survfit(PFS ~ stage, data = dat)
plot(fit.KM, fun= "cloglog", col = 1:2)

fit.KM <- survfit(Surv(ttr, relapse) ~ grp, data = pharmacoSmoking)

plot(fit.KM, fun = "cloglog", col = 1:2)


### Schoenfeld residuals
fit <- coxph(PFS ~ stage, data = dat)
residual.sch <- cox.zph(fit)
residual.sch


plot(residual.sch)

# Dealing with assumptions violations

## Stratification

library(asaur)
d <- pharmacoSmoking
d$employment <- ifelse(d$employment == "ft", "ft", "other")

table(d$employment)

# Stratified Cox model:
fit <- coxph(Surv(ttr, relapse) ~ grp + strata(employment), data = d)

summary(fit)

# Note how there is no estimate associated with 'employment'.

## Truncation
library(asaur)
library(survival)
d <- pancreatic2

plot(survfit(Surv(pfs, status) ~ stage, data = d), col = 1:2)

# THIS IS *NOT* HOW IT IS DONE:
d_WRONG <- subset(d, pfs <= 180)

plot(survfit(Surv(pfs, status) ~ stage, data = d_WRONG), col = 1:2)

# Here is how you do it:
d_RIGHT <- within(d, {
  status_truncated <- ifelse(pfs > 180, 0, status)
  pfs_truncated <- ifelse(pfs > 180, 180, pfs)
})

plot(survfit(Surv(pfs, status) ~ stage, data = d_RIGHT),
     col = 1:2)

plot(survfit(Surv(pfs, status) ~ stage, data = d_RIGHT),
     fun = "cloglog",
     col = 1:2)

summary(coxph(Surv(pfs, status) ~ stage, data = d_RIGHT))
