## ----.common, include = FALSE-------------------------------------------------
source("_common.R")


## ----setup-ife-mc, echo = FALSE, message = FALSE, warning = FALSE-------------
set.seed(1234)
data(simdata)


## ----simdata_ife, eval=TRUE, cache = TRUE-------------------------------------
out.ife <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
          force = "two-way", method = "ife", CV = TRUE, r = c(0, 5),
          se = TRUE, nboots = 200, parallel = TRUE, cores = 16)
print(out.ife)


## ----plot-att-ife,  fig.width = 6, fig.height = 4.5---------------------------
plot(out.ife, main = "Estimated ATT (IFEct)")


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
out.all <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                method = "ife", CV = TRUE, r = c(0, 5),
                cv.method = "all_units", se = FALSE, parallel = TRUE, cores = 16)

out.tr <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
               method = "ife", CV = TRUE, r = c(0, 5),
               cv.method = "treated_units", se = FALSE, parallel = TRUE, cores = 16)


## ----print-cv-method-compare--------------------------------------------------
cat("cv.method = 'all_units':     r.cv =", out.all$r.cv, "\n")
cat("cv.method = 'treated_units': r.cv =", out.tr$r.cv, "\n")


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


## ----placebo_ife_plot, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.ife.p, ylab = "Effect of D on Y", main = "Estimated ATT (IFE)",
     cex.text = 0.8, stats = c("placebo.p","equiv.p"))


## ----placebo_mc_plot, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.mc.p, cex.text = 0.8, stats = c("placebo.p","equiv.p"),
     main = "Estimated ATT (MC)")


## ----simdata_ife_loo, eval=TRUE, cache = TRUE, message = FALSE, results = 'hide'----
out.ife.loo <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
  method = "ife", force = "two-way", se = TRUE, parallel = TRUE, cores = 16, nboots = 200, loo = TRUE)
out.mc.loo <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
  method = "mc", force = "two-way", se = TRUE, parallel = TRUE, cores = 16, nboots = 200, loo = TRUE)


## ----pretrend_ife, eval = TRUE, cache = TRUE, fig.width = 6, fig.height = 4.5, warning = FALSE----
plot(out.ife.loo, type = "equiv", ylim = c(-4,4), loo = TRUE,
     cex.legend = 0.6, main = "Testing Pre-Trend (IFEct)", cex.text = 0.8)


## ----pretrend_mc, eval = TRUE, cache = TRUE, fig.width = 6, fig.height = 4.5, warning = FALSE----
plot(out.mc.loo, type = "equiv", ylim = c(-4,4), loo = TRUE,
     cex.legend = 0.6, main = "Testing Pre-Trend (MC)", cex.text = 0.8)


## ----carryover_ife, eval = TRUE, cache = TRUE, message = FALSE, results='hide'----
out.ife.c <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "ife", r = 2, CV = 0,
  parallel = TRUE, cores = 16, se = TRUE,
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))

out.mc.c <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "mc",  lambda = out.mc$lambda.cv,
  CV = 0, parallel = TRUE, cores = 16, se = TRUE,
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))


## ----carryover_ife_plot, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 5----
plot(out.ife.c, type = "exit", ylim = c(-2.5,4.5),
          cex.text = 0.8, main = "Carryover Effects (IFE)")


## ----carryover_mc_plot, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 5----
plot(out.mc.c, type = "exit", ylim = c(-2.5,4.5),
          cex.text = 0.8, main = "Carryover Effects (MC)")


## ----carryover_rm, eval = TRUE, cache = TRUE, message = FALSE, results='hide', fig.width = 6, fig.height = 4.5----
out.ife.rm.test <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "ife", r = 2, CV = 0,
  parallel = TRUE, cores = 16, se = TRUE,  carryover.rm = 3,
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))# remove three periods

plot(out.ife.rm.test, cex.text = 0.8, stats.pos = c(5, 2.5))

