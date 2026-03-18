##########################
# The Imputation Estimator
##########################
rm(list = ls())
library(fect)
data(fect)

set.seed(1234)

##########################
# From: ./02-fect.Rmd
##########################

## Simulated data
library(panelView)
panelview(Y ~ D, data = simdata1, index = c("id","time"),
  axis.lab = "time", xlab = "Time", ylab = "Unit",
  gridOff = TRUE, by.timing = TRUE,
  background = "white", main = "Simulated Data: Treatment Status")

panelview(Y ~ D, data = simdata1, index = c("id","time"),
  axis.lab = "time", xlab = "Time", ylab = "Unit",
  theme.bw = TRUE, type = "outcome",
  main = "Simulated Data: Outcome")

##########################
# The imputation estimator (FEct)
##########################

out.fect <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time"),
                 method = "fe", force = "two-way")

plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y",
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## Uncertainty estimates
out.fect <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time"),
  method = "fe", force = "two-way", se = TRUE,
  cores = 8, parallel = TRUE, nboots = 1000)

plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y",
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8, stats = "F.p")

## Save estimates
print(out.fect)

out.fect$est.att
out.fect$est.avg
out.fect$beta

## Leave-one-out approach
out.fect.loo <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time"),
  method = "fe", force = "two-way", se = TRUE, loo = TRUE,
  cores = 8, parallel = TRUE, nboots = 200)

plot(out.fect.loo, main = "Estimated ATT (FEct) -- LOO",
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# Weighting treatment effects
##########################

## Balanced treated sample
out.bal <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time"),
  balance.period = c(-3, 4), force = "two-way", method = "ife",
  CV = FALSE, r = 2, se = TRUE, nboots = 200, parallel = TRUE)

plot(out.bal, main = "Estimated ATT (Balanced Sample)")

plot(out.bal, main = "Estimated ATT (Balanced Sample)",
  color = "red", count.color = "blue")

## Weighted treatment effect
simdata1$Weight <- abs(rnorm(n = dim(simdata1)[1]))
out.w <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time"),
  force = "two-way", method = "ife", W = 'Weight',
  CV = FALSE, r = 2, se = TRUE, nboots = 200, parallel = TRUE)

plot(out.w, main = "Estimated Weighted ATT")

##########################
# Effect heterogeneity
##########################

## Box plots
plot(out.fect, type = "box", xlim = c(-15, 10))

## CATT by calendar time
plot(out.fect, type = "calendar", xlim = c(1, 35))

## CATT by a covariate
plot(out.fect, type = "hte", covariate = "X1")

simdata1$X3 <- sample(1:3, size = nrow(simdata1), replace = TRUE)
out.fect.X3 <- fect(Y ~ D + X1 + X2 + X3, data = simdata1, index = c("id","time"),
                   method = "fe", se = TRUE, seed = 123,
                   cores = 8, nboots = 1000, parallel = TRUE)

plot(out.fect.X3, type="hte", covariate = "X3",
     xlab = "", ylab = "Effet of D on Y",
     covariate.labels = c("USA", "China", "UK"),
     ylim = c(-2, 6))

##########################
# Quick diagnostic check
##########################
plot(out.fect, type = "equiv", ylim = c(-4,4),
     cex.legend = 0.6, main = "Testing Pre-Trend (FEct)", cex.text = 0.8)

##########################
# Other estimators: CFE demo
##########################
simdata1[,"FE3"] <- sample(c(1,2,3,4,5), size = dim(simdata1)[1], replace = TRUE)
out.cfe <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time","FE3"),
  method = "cfe", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200,
  Q.type = "linear")
plot(out.cfe)

##########################
# Other options
##########################

## Remove count bar
plot(out.fect, show.count = FALSE)

## Counterfactual plot
plot(out.fect, type = "counterfactual")

## Status plot
plot(out.fect, type = 'status', axis.lab = "time", cex.axis  = 0.6)

##########################
# Average cohort treatment effect
##########################
panelview(Y ~ D, data = simdata1, index = c("id","time"), by.timing = TRUE,
  axis.lab = "time", xlab = "Time", ylab = "Unit",
  background = "white", main = "Simulated Data: Treatment Status")

simdata1.cohort <- get.cohort(data = simdata1, D = 'D', index = c("id","time"))
print(table(simdata1.cohort[,'Cohort']))

simdata1.cohort2 <- get.cohort(data = simdata1, D = 'D', index = c("id","time"),
                               entry.time = list(c(21,27),c(30,33)))
print(table(simdata1.cohort2[,'Cohort']))

out.fe.g <- fect(Y ~ D + X1 + X2, data = simdata1.cohort, index = c("id","time"),
          force = "two-way", method = "fe",
          se = TRUE, nboots = 200, parallel = TRUE, group = 'Cohort')

plot(out.fe.g, show.group = "Cohort:22",
          xlim = c(-15, 10), ylim = c(-10, 10))
