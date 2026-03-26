## ----.common, include = FALSE-------------------------------------------------
source("_common.R")


## ----setup-seed, echo = FALSE-------------------------------------------------
set.seed(1234)


## ----load-packages, message = FALSE, warning = FALSE--------------------------
data(sim_base)
data(sim_gsynth)


## ----panelview-treatment, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5----
library(panelView)
panelview(Y ~ D, data = sim_base, index = c("id","time"),
  axis.lab = "time", xlab = "Time", ylab = "Unit",
  gridOff = TRUE, by.timing = TRUE,
  background = "white", main = "Simulated Data: Treatment Status")


## ----panelview-outcome, fig.width = 6, fig.height = 4.5, warning = FALSE------
panelview(Y ~ D, data = sim_base, index = c("id","time"),
  axis.lab = "time", xlab = "Time", ylab = "Unit",
  theme.bw = TRUE, type = "outcome", 
  main = "Simulated Data: Outcome")


## ----simdata_fect_nose, eval=TRUE, cache = TRUE, message = FALSE, results = 'hide'----
out.fect <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
                 method = "fe", force = "two-way")


## ----fect_plot_nose, fig.width = 6, fig.height = 4.5--------------------------
plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y",
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)


## ----simdata_fect, eval=TRUE, cache = TRUE, message = FALSE, results = 'hide'----
out.fect <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  method = "fe", force = "two-way", se = TRUE,
  parallel = TRUE, cores = 16, nboots = 1000)


## ----fect_plot_nse, fig.width = 6, fig.height = 4.5---------------------------
plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y",
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8, stats = "F.p")


## ----exit_fect, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.fect, type = "exit", main = "Exit Plot (FEct)")


## ----print-fect---------------------------------------------------------------
print(out.fect)


## ----extract-estimates, eval = FALSE------------------------------------------
# out.fect$est.att
# out.fect$est.avg
# out.fect$beta


## ----extract-bootstrap, eval = FALSE------------------------------------------
# out.fect$eff.boot


## ----fect_placebo, eval=TRUE, cache=TRUE, message=FALSE, results='hide', fig.width=6, fig.height=4.5----
out.fect.placebo <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  force = "two-way", method = "fe",
  se = TRUE, nboots = 1000, parallel = TRUE, cores = 16,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.fect.placebo, cex.text = 0.8)


## ----fect_carryover, eval=TRUE, cache=TRUE, message=FALSE, results='hide'-----
out.fect.carry <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  force = "two-way", method = "fe",
  se = TRUE, nboots = 1000, parallel = TRUE, cores = 16,
  carryoverTest = TRUE, carryover.period = c(1, 3))


## ----fect_carryover_plot, eval=TRUE, cache=TRUE, warning=FALSE, fig.width=6, fig.height=5----
plot(out.fect.carry, type = "exit", cex.text = 0.8, main = "Carryover Effects (FEct)")


## ----fect_loo, eval=TRUE, cache = TRUE, message = FALSE-----------------------
out.fect.loo <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  method = "fe", force = "two-way", se = TRUE, loo = TRUE,
  parallel = TRUE, cores = 16, nboots = 1000)


## ----plot-gap-loo, fig.width = 6, fig.height = 4.5----------------------------
plot(out.fect.loo,main = "Estimated ATT (FEct) -- LOO",
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)


## ----cumu_effect, cache = TRUE------------------------------------------------
out <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                        method = "ife", time.component.from = "nevertreated",
                        force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 1000, vartype = 'bootstrap',
                        parallel = TRUE, cores = 16, keep.sims=TRUE)
cumu.out <- effect(out)


## ----cumu_effect_plot, cache = TRUE-------------------------------------------
print(cumu.out)
plot(cumu.out)


## ----cumu_effect_byperiod, cache = TRUE---------------------------------------
effect(out, cumu=FALSE)


## ----cumu_effect_subset, cache = TRUE-----------------------------------------
effect(out, cumu=TRUE, id=c(101,102,103), period=c(1,5))


## ----effect-mc, cache = TRUE--------------------------------------------------
out_mc <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                        method = "mc", force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 1000, vartype = 'bootstrap',
                        parallel = TRUE, cores = 16, keep.sims=TRUE)
plot(effect(out_mc))


## ----effect-jackknife, cache = TRUE-------------------------------------------
out_jack <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                        method = "mc", force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 1000, vartype = 'jackknife',
                        parallel = TRUE, cores = 16, keep.sims=TRUE)
plot(effect(out_jack))


## ----simdata_bal, eval=TRUE, cache = TRUE-------------------------------------
out.bal <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  balance.period = c(-3, 4), force = "two-way", method = "ife",
  CV = FALSE, r = 2, se = TRUE, nboots = 1000, parallel = TRUE, cores = 16)


## ----plot-balanced-att, fig.width = 6, fig.height = 4.5-----------------------
plot(out.bal, main = "Estimated ATT (Balanced Sample)")


## ----plot-balanced-custom, fig.width = 6, fig.height = 4.5--------------------
plot(out.bal, main = "Estimated ATT (Balanced Sample)",
  post.color = "red", count.color = "blue")


## ----simdata_panelview_cohort, fig.width = 6, fig.height = 4.5, warning = FALSE----
panelview(Y ~ D, data = sim_base, index = c("id","time"), by.timing = TRUE,
  axis.lab = "time", xlab = "Time", ylab = "Unit",
  background = "white", main = "Simulated Data: Treatment Status")


## ----get_cohort---------------------------------------------------------------
sim_base.cohort <- get.cohort(data = sim_base,D = 'D',index = c("id","time"))
print(table(sim_base.cohort[,'Cohort']))


## ----get_cohort2--------------------------------------------------------------
sim_base.cohort2 <- get.cohort(data = sim_base,D = 'D',index = c("id","time"),
                               entry.time = list(c(21,27),c(30,33)))
print(table(sim_base.cohort2[,'Cohort']))


## ----simdata_fe_cohort, eval = TRUE, cache = TRUE, message = FALSE, results='hide'----
out.fe.g <- fect(Y ~ D + X1 + X2, data = sim_base.cohort, index = c("id","time"),
          force = "two-way", method = "fe",
          se = TRUE, nboots = 1000, parallel = TRUE, cores = 16, group = 'Cohort')


## ----cohort_plot1, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.fe.g, show.group = "Cohort:22",
          xlim = c(-15, 10), ylim = c(-10, 10))


## ----simdata_w, eval=TRUE, cache = TRUE---------------------------------------
sim_base$Weight <- abs(rnorm(n = dim(sim_base)[1]))
out.w <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  force = "two-way", method = "ife", W = 'Weight',
  CV = FALSE, r = 2, se = TRUE, nboots = 1000, parallel = TRUE, cores = 16)


## ----plot-weighted-att, fig.width = 6, fig.height = 4.5-----------------------
plot(out.w, main = "Estimated Weighted ATT")

