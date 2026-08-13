## ----.common, include = FALSE-------------------------------------------------
source("_common.R")


## ----setup-cfe, echo = FALSE--------------------------------------------------
set.seed(1234)


## ----load-packages-cfe, message = FALSE, warning = FALSE----------------------
data(simdata)
data(sim_region)
data(sim_linear)
data(sim_trend)


## ----cfe-42-load, eval = TRUE-------------------------------------------------
head(sim_region)


## ----cfe-42-fe-only, eval = TRUE, cache = TRUE, message = FALSE, warning = FALSE, results = 'hide'----
out.fe.only <- fect(Y ~ D, data = sim_region,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, cores = 16, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))


## ----cfe-42-fe-only-plot, fig.width = 6, fig.height = 4.5---------------------
plot(out.fe.only, 
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only — Placebo Test")


## ----cfe-42-with-region, eval = TRUE, cache = TRUE, message = FALSE, warning = FALSE, results = 'hide'----
out.cfe.region <- fect(Y ~ D, data = sim_region,
  index = c("id", "time", "region_time"),
  method = "cfe", force = "two-way",
  se = TRUE, parallel = TRUE, cores = 16, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))


## ----cfe-42-with-region-plot, fig.width = 6, fig.height = 4.5-----------------
plot(out.cfe.region, 
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with Region×Time FE — Placebo Test")


## ----cfe-43-fe-baseline, eval = TRUE, cache = TRUE, message = FALSE, warning = FALSE, results = 'hide'----
out.fe.base <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, cores = 16, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))


## ----cfe-43-fe-baseline-plot, fig.width = 6, fig.height = 4.5-----------------
plot(out.fe.base, 
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only (simdata) — Placebo Test")


## ----cfe-43-gamma-setup, eval = TRUE------------------------------------------
simdata$gamma_t <- simdata$time


## ----cfe-43-with-z, eval = TRUE, cache = TRUE, message = FALSE, warning = FALSE, results = 'hide'----
out.cfe.z <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  se = TRUE, parallel = TRUE, cores = 16, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))


## ----cfe-43-with-z-plot, fig.width = 6, fig.height = 4.5----------------------
plot(out.cfe.z, 
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with Z = L1 — Placebo Test")


## ----cfe-44-linear-load, eval = TRUE------------------------------------------
head(sim_linear)


## ----cfe-44-lin-fe-only, eval = TRUE, cache = TRUE, message = FALSE, warning = FALSE, results = 'hide'----
out.fe.lin <- fect(Y ~ D, data = sim_linear,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, cores = 16, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))


## ----cfe-44-lin-fe-only-plot, fig.width = 6, fig.height = 4.5-----------------
plot(out.fe.lin, 
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only (Linear Trend DGP) — Placebo Test")


## ----cfe-44-lin-cfe, eval = TRUE, cache = TRUE, message = FALSE, warning = FALSE, results = 'hide'----
out.cfe.lin <- fect(Y ~ D, data = sim_linear,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "linear",
  se = TRUE, parallel = TRUE, cores = 16, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))


## ----cfe-44-lin-cfe-plot, fig.width = 6, fig.height = 4.5---------------------
plot(out.cfe.lin, 
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with Linear Trend — Placebo Test")


## ----cfe-44-sin-load, eval = TRUE---------------------------------------------
head(sim_trend)


## ----cfe-44-sin-fe-only, eval = TRUE, cache = TRUE, message = FALSE, warning = FALSE, results = 'hide'----
out.fe.trend <- fect(Y ~ D, data = sim_trend,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, cores = 16, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))


## ----cfe-44-sin-fe-only-plot, fig.width = 6, fig.height = 4.5-----------------
plot(out.fe.trend, 
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only (Sin Trend DGP) — Placebo Test")


## ----cfe-44-sin-bspline, eval = TRUE, cache = TRUE, message = FALSE, warning = FALSE, results = 'hide'----
out.cfe.bs <- fect(Y ~ D, data = sim_trend,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "bspline",
  se = TRUE, parallel = TRUE, cores = 16, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))


## ----cfe-44-sin-bspline-plot, fig.width = 6, fig.height = 4.5-----------------
plot(out.cfe.bs, 
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with B-spline Trend — Placebo Test")


## ----cfe-45-gamma-setup, eval = TRUE------------------------------------------
simdata$gamma_t <- simdata$time


## ----cfe-45-fit-models, eval = TRUE, cache = TRUE, message = FALSE, warning = FALSE, results = 'hide'----
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


## ----cfe-45-mspe, eval = TRUE, cache = TRUE-----------------------------------
mspe.out <- fect_mspe(list(FE = out.fe,
       CFE_Z = out.cfe.z.only,
       CFE_Z_F1 = out.cfe.z.f1,
       IFE_r2 = out.ife.r2),
  seed = 1234)
print(mspe.out$summary[, c("Model", "MSPE", "RMSE", "MAD")])


## ----cfe-45-best-placebo, eval = TRUE, cache = TRUE, message = FALSE, warning = FALSE, results = 'hide'----
out.cfe.best <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  r = 1,
  se = TRUE, parallel = TRUE, cores = 16, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))


## ----cfe-45-best-placebo-plot, fig.width = 6, fig.height = 4.5----------------
plot(out.cfe.best, 
     stats = c("placebo.p", "equiv.p"),
     main = "CFE (Z + 1 Factor) — Placebo Test")


## ----cfe-46-zparam-example, eval = FALSE--------------------------------------
# # Example with Z.param (not run — requires appropriate data)
# # out <- fect(Y ~ D, data = mydata,
# #   index = c("unit", "time"),
# #   method = "cfe", force = "two-way",
# #   Z = c("baseline_gdp", "baseline_pop"),
# #   gamma = c("decade", "political_era"),
# #   Z.param = list(decade = "baseline_gdp",
# #                  political_era = "baseline_pop"))


## ----cfe-cv-dispatcher, eval = FALSE------------------------------------------
# fit.cfe <- fect(Y ~ D, data = sim_region,
#   index   = c("id", "time", "region_time"),
#   method  = "cfe", force = "two-way",
#   CV      = TRUE, r = c(0, 3),
#   cv.method = "rolling",
#   cv.buffer = 1, cv.nobs = 3, k = 20, cv.prop = 0.1,
#   cv.rule = "1se",
#   se = TRUE, parallel = TRUE, cores = 16, nboots = 200
# )
# fit.cfe$r.cv  # selected r

