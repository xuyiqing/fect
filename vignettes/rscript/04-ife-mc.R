## ----.common, include = FALSE-------------------------------------------------
source("_common.R")


## ----setup-ife-mc, echo = FALSE, message = FALSE, warning = FALSE-------------
set.seed(1234)
data(simdata)


## ----panelview-treatment-ifemc, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5----
library(panelView)
panelview(Y ~ D, data = simdata, index = c("id", "time"),
  axis.lab = "time", xlab = "Time", ylab = "Unit",
  gridOff = TRUE, by.timing = TRUE,
  background = "white", main = "simdata: Treatment Status")


## ----panelview-outcome-ifemc, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5----
panelview(Y ~ D, data = simdata, index = c("id", "time"),
  axis.lab = "time", xlab = "Time", ylab = "Outcome",
  theme.bw = TRUE, type = "outcome", by.group = FALSE,
  main = "simdata: Outcome")


## ----simdata_ife, eval=TRUE, cache = TRUE-------------------------------------
out.ife <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
          force = "two-way", method = "ife", CV = TRUE, r = c(0, 5),
          se = TRUE, nboots = 200, parallel = TRUE, cores = 16)
print(out.ife)


## ----plot-att-ife,  fig.width = 6, fig.height = 4.5---------------------------
plot(out.ife, main = "Estimated ATT (IFEct)")


## ----plot-factors-ifemc, fig.width = 6, fig.height = 4------------------------
plot(out.ife, type = "factors", main = "Estimated Factors")


## ----plot-loadings-ifemc, fig.width = 6, fig.height = 5-----------------------
plot(out.ife, type = "loadings", main = "Factor Loadings")


## ----plot-loading-overlap-ifemc, fig.width = 6, fig.height = 5----------------
plot(out.ife, type = "loading.overlap")


## ----simdata_mc, eval=TRUE, cache = TRUE--------------------------------------
out.mc <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
          force = "two-way", method = "mc", CV = TRUE,
          se = TRUE, nboots = 200, parallel = TRUE, cores = 16)

print(out.mc)


## ----plot-att-mc, fig.width = 6, fig.height = 4.5-----------------------------
plot(out.mc, main = "Estimated ATT (MC)")


## ----cv_ife_demo, eval=TRUE, cache=TRUE, message=FALSE, results='hide'--------
out.cv <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
               method = "ife", CV = TRUE, r = c(0, 5),
               se = FALSE, parallel = TRUE, cores = 16)


## ----print-cv-selected-r------------------------------------------------------
cat("Selected r:", out.cv$r.cv, "\n")


## ----cv_method_compare, eval=TRUE, cache=TRUE, message=FALSE, results='hide'----
out.roll <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                 method = "ife", CV = TRUE, r = c(0, 5),
                 cv.method = "rolling", se = FALSE, parallel = TRUE, cores = 16)

out.block <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                  method = "ife", CV = TRUE, r = c(0, 5),
                  cv.method = "block", se = FALSE, parallel = TRUE, cores = 16)


## ----print-cv-method-compare--------------------------------------------------
cat("cv.method = 'rolling': r.cv =", out.roll$r.cv, "\n")
cat("cv.method = 'block':   r.cv =", out.block$r.cv, "\n")


## ----cv-strategies-fig, echo=FALSE, out.width='100%', fig.cap='Block CV (left) versus rolling-window CV (right) on a synthetic 40-unit panel with staggered treatment timing. Both panels use the same random missing pattern. Block CV (panel a) drops random anchors and masks contiguous holdouts (red) flanked by donut buffer (orange). Rolling CV (panel b) samples a fraction of eligible units per fold; for each sampled unit it masks a past-side buffer (orange), a scored holdout (red), and drops everything from the holdout onward (purple). Treated post-treatment cells (dark blue) are never masked.'----
knitr::include_graphics("fig/cv-strategies.png")


## ----rcv_dispatcher_demo, eval=FALSE------------------------------------------
# fit <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
#             method = "ife", force = "two-way",
#             CV = TRUE, r = c(0, 5),
#             cv.method = "rolling",
#             cv.buffer = 1, cv.nobs = 3, k = 20, cv.prop = 0.1,
#             cv.rule = "1se", se = TRUE)
# fit$r.cv      # selected r


## ----criterion_compare, eval=TRUE, cache=TRUE, message=FALSE, results='hide'----
out.mspe <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                 method = "ife", CV = TRUE, r = c(0, 5),
                 criterion = "mspe", se = FALSE, parallel = TRUE, cores = 16)

out.pc <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
               method = "ife", CV = TRUE, r = c(0, 5),
               criterion = "gmspe", se = FALSE, parallel = TRUE, cores = 16)


## ----print-criterion-compare--------------------------------------------------
cat("criterion = 'mspe': r.cv =", out.mspe$r.cv, "\n")
cat("criterion = 'gmspe': r.cv =", out.pc$r.cv, "\n")


## ----placebo_ife, eval = TRUE, cache = TRUE, message = FALSE, results='hide'----
out.ife.p <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "ife",  r = 2, CV = 0,
  parallel = TRUE, cores = 16, se = TRUE,
  nboots = 200, placeboTest = TRUE, placebo.period = c(-2, 0))

out.mc.p <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "mc",  lambda = out.mc$lambda.cv,
  CV = 0, parallel = TRUE, cores = 16, se = TRUE,
  nboots = 200, placeboTest = TRUE, placebo.period = c(-2, 0))


## ----placebo_ife_plot, eval = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.ife.p, ylab = "Effect of D on Y", main = "Estimated ATT (IFE)",
     stats = c("placebo.p","equiv.p"))


## ----placebo_mc_plot, eval = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.mc.p, stats = c("placebo.p","equiv.p"),
     main = "Estimated ATT (MC)")


## ----simdata_ife_loo, eval=TRUE, cache = TRUE, message = FALSE, results = 'hide'----
out.ife.loo <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
  method = "ife", force = "two-way", se = TRUE, parallel = TRUE, cores = 16, nboots = 200, loo = TRUE)
out.mc.loo <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
  method = "mc", force = "two-way", se = TRUE, parallel = TRUE, cores = 16, nboots = 200, loo = TRUE)


## ----pretrend_ife, eval = TRUE, fig.width = 6, fig.height = 4.5, warning = FALSE----
plot(out.ife.loo, type = "equiv", ylim = c(-4,4), loo = TRUE,
     main = "Testing Pre-Trend (IFEct)")


## ----pretrend_mc, eval = TRUE, fig.width = 6, fig.height = 4.5, warning = FALSE----
plot(out.mc.loo, type = "equiv", ylim = c(-4,4), loo = TRUE,
     main = "Testing Pre-Trend (MC)")


## ----carryover_ife, eval = TRUE, cache = TRUE, message = FALSE, results='hide'----
out.ife.c <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "ife", r = 2, CV = 0,
  parallel = TRUE, cores = 16, se = TRUE,
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))

out.mc.c <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "mc",  lambda = out.mc$lambda.cv,
  CV = 0, parallel = TRUE, cores = 16, se = TRUE,
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))


## ----carryover_ife_plot, eval = TRUE, warning = FALSE, fig.width = 6, fig.height = 5----
plot(out.ife.c, type = "exit", main = "Carryover Effects (IFE)")


## ----carryover_mc_plot, eval = TRUE, warning = FALSE, fig.width = 6, fig.height = 5----
plot(out.mc.c, type = "exit", ylim = c(-2.5,4.5),
          main = "Carryover Effects (MC)")


## ----carryover_rm_fit, eval = TRUE, cache = TRUE, message = FALSE, results='hide'----
out.ife.rm.test <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "ife", r = 2, CV = 0,
  parallel = TRUE, cores = 16, se = TRUE,  carryover.rm = 3,
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))# remove three periods


## ----carryover_rm, eval = TRUE, fig.width = 6, fig.height = 4.5---------------
plot(out.ife.rm.test)


## ----carryover_rm_only_test, eval = TRUE, fig.width = 6, fig.height = 4.5-----
plot(out.ife.rm.test, highlight = "carryover",
     main = "Highlight only the carryover-test periods (blue diamonds)")


## ----carryover_rm_only_removed, eval = TRUE, fig.width = 6, fig.height = 4.5----
plot(out.ife.rm.test, highlight = "carryover.rm",
     main = "Highlight only the removed periods (orange triangles)")


## ----carryover_rm_fill, eval = TRUE, fig.width = 6, fig.height = 4.5----------
plot(out.ife.rm.test, highlight.fill = TRUE,
     main = "Both test types highlighted, with rectangles")

