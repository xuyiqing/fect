##########################
# Complex Fixed Effects (CFE)
##########################
rm(list = ls())
library(fect)
data(fect)

set.seed(1234)

##########################
# Basic CFE: additional additive fixed effects
##########################
simdata[,"FE3"] <- sample(1:5, size = nrow(simdata), replace = TRUE)

out.cfe.basic <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time", "FE3"),
  method = "cfe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200)
plot(out.cfe.basic, main = "CFE with Additional Fixed Effect",
     ylab = "Effect of D on Y",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# Unit-specific time trends via Q.type
##########################

# Linear time trend
out.cfe.linear <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "linear",
  se = TRUE, parallel = TRUE, nboots = 200)
plot(out.cfe.linear, main = "CFE with Unit-Specific Linear Trend",
     ylab = "Effect of D on Y",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

# Quadratic time trend
out.cfe.quad <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = c("linear", "quadratic"),
  se = TRUE, parallel = TRUE, nboots = 200)
plot(out.cfe.quad, main = "CFE with Unit-Specific Quadratic Trend",
     ylab = "Effect of D on Y",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

# B-spline time trend
out.cfe.bs <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "bspline", Q.bspline.degree = 3,
  se = TRUE, parallel = TRUE, nboots = 200)
plot(out.cfe.bs, main = "CFE with B-spline Time Trend",
     ylab = "Effect of D on Y",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# Time-invariant covariates with time-varying coefficients (Z/gamma)
##########################
set.seed(42)
unit_values <- rnorm(length(unique(simdata$id)))
names(unit_values) <- unique(simdata$id)
simdata$Z_baseline <- unit_values[as.character(simdata$id)]
simdata$time_group <- ifelse(simdata$time <= 17, 1, 2)

out.cfe.z <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "Z_baseline", gamma = "time_group",
  se = TRUE, parallel = TRUE, nboots = 200)
plot(out.cfe.z, main = "CFE with Time-Invariant Z",
     ylab = "Effect of D on Y",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# CFE with interactive fixed effects
##########################
out.cfe.ife <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "linear",
  r = 2,
  se = TRUE, parallel = TRUE, nboots = 200)
plot(out.cfe.ife,
     main = "CFE with Linear Trend + 2 Factors",
     ylab = "Effect of D on Y",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# Combining multiple components
##########################
out.cfe.full <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time", "FE3"),
  method = "cfe", force = "two-way",
  Q.type = c("linear", "quadratic"),
  r = 2,
  se = TRUE, parallel = TRUE, nboots = 200)
plot(out.cfe.full,
     main = "CFE: Full Model",
     ylab = "Effect of D on Y",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# Diagnostic tests
##########################

# Pre-trend test
plot(out.cfe.linear, type = "equiv", ylim = c(-4, 4),
     cex.legend = 0.6, main = "Pre-Trend Test (CFE)", cex.text = 0.8)

# Placebo test
out.cfe.p <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "linear",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.cfe.p, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "Placebo Test (CFE)")

# Carryover test
out.cfe.c <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "linear",
  se = TRUE, parallel = TRUE, nboots = 200,
  carryoverTest = TRUE, carryover.period = c(1, 3))
plot(out.cfe.c, type = "exit", ylim = c(-2.5, 4.5),
     cex.text = 0.8, main = "Carryover Test (CFE)")

