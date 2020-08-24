---
  title: "Penalized Cox Regression"
---

  # Simulated data

library(survival)
library(glmnet)


## Data generation

set.seed(1234)

N <- 1000
p <- 30  # total num. features
nzc <- p/3 # num. 'true' predictors

X <- matrix(rnorm(N * p), nrow = N, ncol = p)
beta <- rnorm(nzc)
y <- local({
  linear_predictor <- X[, seq_len(nzc)] %*% beta
  hazard <- exp(linear_predictor)
  y_time_event <- rexp(N, rate = hazard)
  y_time_censoring <- rexp(N, rate = hazard * 0.5)
  y_time <- pmin(y_time_event, y_time_censoring)
  y_event <- y_time_event <= y_time_censoring
  Surv(y_time, y_event)
})


# We should have about 1/3 of the data points censored:

table(y[,2])


## Fit the model
fit <- glmnet(X, y, family = "cox")
plot(fit)

## Selecting a threshold through cross-validation
set.seed(1234)
fit.cv10 <- cv.glmnet(X, y, family = "cox")
plot(fit.cv10)

str(fit.cv10)


# Estimated coefficients:
coef(fit.cv10, s = "lambda.1se")


# Compare estimated values with 'true' values:
beta.true <- c(round(beta, 2), rep(0, ncol(X) - length(beta)))
beta.est <- round(coef(fit.cv10, s = "lambda.1se"), 2)
cbind(beta.true, beta.est)

## Make predictions
# Using built-in functions:
predict(fit.cv10, newx = X[1:5, ], s = "lambda.1se")

# 'Manually':
b <- coef(fit.cv10, s = "lambda.1se")
b.i <- which(b != 0)
bnz <- b[b.i]
y0 <- X[1:5, b.i, drop = FALSE] %*% bnz
print(y0)

# Case study nr. 1

## Load data and packages
library(survival)
library(glmnet)

load('data/LymphomaData.rda')

str(patient.data)

gex <- t(patient.data$x)
y <- Surv(patient.data$time, patient.data$status)

## Exploratory analysis

plot(survfit(y ~ 1))

hist(as.vector(gex))

hist(colMeans(gex))

hist(apply(gex, 2, sd))


## Split the data randomly into a training and a testing set

i.training <- sample.int(nrow(gex), size = 160, replace = FALSE)
i.testing <- setdiff(seq_len(nrow(gex)), i.training)

gex.training <- gex[i.training,, drop = FALSE]
y.training <- y[i.training,, drop = FALSE]

gex.testing <- gex[i.testing,, drop = FALSE]
y.testing <- y[i.testing,, drop = FALSE]

## Train the model
fit.cv10 <- cv.glmnet(gex.training, y.training, family = "cox")

plot(fit.cv10)

b <- coef(fit.cv10, s = "lambda.min")
sum(b != 0)

round(b[b != 0], digits = 3)


## Test the model
score.testing <- predict(fit.cv10, newx = gex.testing, s = "lambda.min")
hist(score.testing)

# Question: how well is the score predicting survival?

# A continuous predictor vs a right-censored time-to-failure outcome: Cox regression!
summary(coxph(y.testing ~ score.testing))

# A more interpretable scale: by IQR variation:
score_scaled.testing <- score.testing / IQR(score.testing)
summary(coxph(y.testing ~ score_scaled.testing))

# We can split the scores into 2 categories, and compare patients with 'low' vs 'high' score:
gex_risk <- ifelse(score.testing <= median(score.testing), "low", "high")
table(gex_risk)

fit.KM <- survfit(y.testing ~ gex_risk, conf.type = "log-log")
plot(fit.KM)

fit.KM

# Question: how is the 6 months survival for patients classified as low risk, compared to patients classified as high risk?
summary(fit.KM, time = 6)

library(survivalROC)
ROC <- survivalROC(Stime = y.testing[, 1],
                   status = y.testing[, 2],
                   marker = score.testing,
                   cut.values = quantile(score.testing, prob = 0:10/10),
                   predict.time = 10,
                   method = "KM")
ROC$AUC

with(ROC, plot(TP ~ FP, type = "l", xlim = c(0, 1), ylim = c(0, 1)))
