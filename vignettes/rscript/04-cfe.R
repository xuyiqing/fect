##########################
# Complex Fixed Effects (CFE)
##########################
rm(list = ls())
library(fect)
data(fect)

set.seed(1234)

##########################
# 4.2 Additional fixed effects
##########################

## sim_region ships with the package (see data-raw/sim_region.R for DGP)
head(sim_region)

## FE only
out.fe.only <- fect(Y ~ D, data = sim_region,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.fe.only, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## With region FE
out.cfe.region <- fect(Y ~ D, data = sim_region,
  index = c("id", "time", "region"),
  method = "cfe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.cfe.region, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with Region FE — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# 4.3 Time-invariant covariates with time-varying coefficients
##########################

## FE baseline (simdata)
out.fe.base <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.fe.base, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only (simdata) — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## With Z = L1
simdata$gamma_t <- simdata$time

out.cfe.z <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.cfe.z, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with Z = L1 — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# 4.4 Unit-specific time trends
##########################

## sim_trend ships with the package (see data-raw/sim_trend.R for DGP)
head(sim_trend)

## FE only
out.fe.trend <- fect(Y ~ D, data = sim_trend,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.fe.trend, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only (sim_trend) — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## Linear trend
out.cfe.lin <- fect(Y ~ D, data = sim_trend,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "linear",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.cfe.lin, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with Linear Trend — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## B-spline trend
out.cfe.bs <- fect(Y ~ D, data = sim_trend,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "bspline",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.cfe.bs, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with B-spline Trend — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# 4.5 CFE with factors
##########################
simdata$gamma_t <- simdata$time

## Model comparison
out.fe <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "fe", force = "two-way", se = FALSE)

out.cfe.z.only <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  se = FALSE)

out.cfe.z.f1 <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  r = 1, se = FALSE)

out.ife.r2 <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "ife", force = "two-way",
  r = 2, se = FALSE)

mspe.out <- fect_mspe(
  list(FE = out.fe,
       CFE_Z = out.cfe.z.only,
       CFE_Z_F1 = out.cfe.z.f1,
       IFE_r2 = out.ife.r2),
  seed = 1234)
print(mspe.out$summary[, c("Model", "MSPE")])

## Best model placebo test
out.cfe.best <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  r = 1,
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.cfe.best, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE (Z + 1 Factor) — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)
