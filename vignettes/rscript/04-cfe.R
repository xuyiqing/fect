##############################
# 04-cfe.R
# Generated from 04-cfe.Rmd
##############################
rm(list = ls())
set.seed(1234)

library(fect)
data(fect)

## --- cfe-42-load ---
head(sim_region)

## --- cfe-42-fe-only ---
out.fe.only <- fect(Y ~ D, data = sim_region,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))

## --- cfe-42-fe-only-plot ---
plot(out.fe.only, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only \u2014 Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## --- cfe-42-with-region ---
out.cfe.region <- fect(Y ~ D, data = sim_region,
  index = c("id", "time", "region_time"),
  method = "cfe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))

## --- cfe-42-with-region-plot ---
plot(out.cfe.region, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with Region\u00d7Time FE \u2014 Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## --- cfe-43-fe-baseline ---
out.fe.base <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))

## --- cfe-43-fe-baseline-plot ---
plot(out.fe.base, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only (simdata) \u2014 Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## --- cfe-43-gamma-setup ---
simdata$gamma_t <- simdata$time

## --- cfe-43-with-z ---
out.cfe.z <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))

## --- cfe-43-with-z-plot ---
plot(out.cfe.z, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with Z = L1 \u2014 Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## --- cfe-44-linear-load ---
head(sim_linear)

## --- cfe-44-lin-fe-only ---
out.fe.lin <- fect(Y ~ D, data = sim_linear,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))

## --- cfe-44-lin-fe-only-plot ---
plot(out.fe.lin, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only (Linear Trend DGP) \u2014 Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## --- cfe-44-lin-cfe ---
out.cfe.lin <- fect(Y ~ D, data = sim_linear,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "linear",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))

## --- cfe-44-lin-cfe-plot ---
plot(out.cfe.lin, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with Linear Trend \u2014 Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## --- cfe-44-sin-load ---
head(sim_trend)

## --- cfe-44-sin-fe-only ---
out.fe.trend <- fect(Y ~ D, data = sim_trend,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))

## --- cfe-44-sin-fe-only-plot ---
plot(out.fe.trend, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only (Sin Trend DGP) \u2014 Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## --- cfe-44-sin-bspline ---
out.cfe.bs <- fect(Y ~ D, data = sim_trend,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "bspline",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))

## --- cfe-44-sin-bspline-plot ---
plot(out.cfe.bs, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with B-spline Trend \u2014 Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## --- cfe-45-gamma-setup ---
simdata$gamma_t <- simdata$time

## --- cfe-45-fit-models ---
# Model 1: FE only
out.fe <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "fe", force = "two-way", se = FALSE)

# Model 2: CFE with Z = L1 only
out.cfe.z.only <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  se = FALSE)

# Model 3: CFE with Z = L1 + 1 factor
out.cfe.z.f1 <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  r = 1, se = FALSE)

# Model 4: IFE with 2 factors
out.ife.r2 <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "ife", force = "two-way",
  r = 2, se = FALSE)

## --- cfe-45-mspe ---
mspe.out <- fect_mspe(
  list(FE = out.fe,
       CFE_Z = out.cfe.z.only,
       CFE_Z_F1 = out.cfe.z.f1,
       IFE_r2 = out.ife.r2),
  seed = 1234)
print(mspe.out$summary[, c("Model", "MSPE", "RMSE", "MAD")])

## --- cfe-45-best-placebo ---
out.cfe.best <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  r = 1,
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))

## --- cfe-45-best-placebo-plot ---
plot(out.cfe.best, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE (Z + 1 Factor) \u2014 Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## --- cfe-46-zparam-example ---
## # Example with Z.param (not run --- requires appropriate data)
## # out <- fect(Y ~ D, data = mydata,
## #   index = c("unit", "time"),
## #   method = "cfe", force = "two-way",
## #   Z = c("baseline_gdp", "baseline_pop"),
## #   gamma = c("decade", "political_era"),
## #   Z.param = list(decade = "baseline_gdp",
## #                  political_era = "baseline_pop"))
