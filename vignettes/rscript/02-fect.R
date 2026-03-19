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
# The FEct estimator
##########################

## Estimation
out.fect <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time"),
                 method = "fe", force = "two-way")

plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y",
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## Inferences (bootstrap)
out.fect <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time"),
  method = "fe", force = "two-way", se = TRUE,
  cores = 8, parallel = TRUE, nboots = 1000)

plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y",
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8, stats = "F.p")

## Exiting the treatment
plot(out.fect, type = "exit", ylim = c(-2.5,4.5), main = "Exit Plot (FEct)")

## Save estimates
print(out.fect)

out.fect$est.att
out.fect$est.avg
out.fect$beta
out.fect$eff.boot

##########################
# Diagnostics
##########################

## Placebo test
out.fect.placebo <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time"),
  force = "two-way", method = "fe",
  se = TRUE, cores = 8, nboots = 200, parallel = TRUE,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.fect.placebo, cex.text = 0.8)

## Carryover test
out.fect.carry <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time"),
  force = "two-way", method = "fe",
  se = TRUE, cores = 8, nboots = 200, parallel = TRUE,
  carryoverTest = TRUE, carryover.period = c(1, 3))
plot(out.fect.carry, type = "exit", cex.text = 0.8, main = "Carryover Effects (FEct)")

## Leave-one-out approach
out.fect.loo <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time"),
  method = "fe", force = "two-way", se = TRUE, loo = TRUE,
  cores = 8, parallel = TRUE, nboots = 200)

plot(out.fect.loo, main = "Estimated ATT (FEct) -- LOO",
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# Other estimands
##########################

## Balanced treated sample
out.bal <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time"),
  balance.period = c(-3, 4), force = "two-way", method = "ife",
  CV = FALSE, r = 2, se = TRUE, nboots = 200, parallel = TRUE)

plot(out.bal, main = "Estimated ATT (Balanced Sample)")

plot(out.bal, main = "Estimated ATT (Balanced Sample)",
  color = "red", count.color = "blue")

## Average cohort treatment effect
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

## User-supplied weights
simdata1$Weight <- abs(rnorm(n = dim(simdata1)[1]))
out.w <- fect(Y ~ D + X1 + X2, data = simdata1, index = c("id","time"),
  force = "two-way", method = "ife", W = 'Weight',
  CV = FALSE, r = 2, se = TRUE, nboots = 200, parallel = TRUE)

plot(out.w, main = "Estimated Weighted ATT")
