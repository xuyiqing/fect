## ----.common, include = FALSE-------------------------------------------------
source("_common.R")


## ----setup-seed, echo = FALSE-------------------------------------------------
set.seed(1234)


## ----load-packages, warning=FALSE, message=FALSE------------------------------
data(sim_gsynth)
data(turnout)
ls()


## ----head-sim_gsynth----------------------------------------------------------
head(sim_gsynth)


## ----sim-panelview-status, cache = FALSE, fig.height=7, fig.width=7, warning=FALSE----
library(panelView)
panelview(Y ~ D, data = sim_gsynth,  index = c("id","time"), pre.post = TRUE) 


## ----sim-panelview-outcome, cache = FALSE,fig.height=5, fig.width=7-----------
panelview(Y ~ D, data = sim_gsynth,  index = c("id","time"), type = "outcome") 


## ----sim2_onecore, cache = TRUE-----------------------------------------------
system.time(
out <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"), 
            method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 5), 
            se = TRUE, nboots = 1000, vartype = 'parametric', 
            parallel = TRUE, cores = 16))


## ----print-results, eval = FALSE----------------------------------------------
# print(out)
# out$est.att
# out$est.avg
# out$beta


## ----sim2, cache = TRUE, warning = FALSE--------------------------------------
system.time(
out <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"), method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 5), se = TRUE, nboots = 1000, vartype = 'parametric', parallel = TRUE, cores = 16)
)


## ----simJack, cache = TRUE, message = FALSE-----------------------------------
out2 <- fect(Y ~ D + X1 + X2, data = sim_gsynth,  index = c("id","time"),
               method = "gsynth", force = "two-way",
               CV = TRUE, r = c(0, 5), se = TRUE,
               vartype = "jackknife",
               parallel = TRUE, cores = 16)



## ----rcv_dispatcher_gsc, eval = FALSE-----------------------------------------
# fit <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id", "time"),
#             method = "gsynth", force = "two-way",
#             CV = TRUE, r = c(0, 5),
#             cv.method = "rolling",
#             cv.buffer = 1, cv.nobs = 3, k = 20, cv.prop = 0.1,
#             cv.rule = "1se",
#             se = TRUE, vartype = "parametric")
# fit$r.cv     # selected r --- use this for downstream inference


## ----sim_gap1, fig.height=5, fig.width=7--------------------------------------
a <- plot(out) # by default, type = "gap"
print(a)


## ----sim-gap-connected, fig.height=5, fig.width=7-----------------------------
plot(out, connected = TRUE)


## ----sim-gap-line, fig.height=5, fig.width=7----------------------------------
plot(out, connected = TRUE, show.points = FALSE)


## ----sim_gap2, fig.height=5, fig.width=7--------------------------------------
plot(out, type = "gap", ylim = c(-6,12), xlab = "Period", 
     main = "Estimated ATT (Gsynth)")


## ----sim-counterfactual, cache = FALSE, fig.height=5, fig.width=7-------------
plot(out, type = "counterfactual")


## ----sim-counterfactual-all, cache = FALSE,fig.height=5, fig.width=7----------
plot(out, type = "counterfactual", raw = "all")


## ----sim-counterfactual-band, cache = FALSE,fig.height=5, fig.width=7---------
plot(out, type = "counterfactual", raw = "band")


## ----sim_status, cache = FALSE,fig.height=7, fig.width=7----------------------
plot(out, type = "status", yticklabels="0", 
     xticklabels=c("5", "10", "15","20", "25", "30") )


## ----sim_L, cache = TRUE, message = FALSE, results='hide', fig.height=7, fig.width=7----
plot(out, type = "loadings")


## ----sim_F, message = FALSE, results='hide', fig.height=5, fig.width=7--------
plot(out, type = "factors", xlab = "Time")


## ----sim_box, cache = FALSE,fig.height=5, fig.width=8-------------------------
plot(out, type = "box", xlab = "time",
     xticklabels=c("-19", "-15", "-10", "-5","0","5","10") )


## ----sim_box2, eval = FALSE, fig.height=7, fig.width=7------------------------
# plot(out, type = "box", xlim = c(-15, 10),
#      xticklabels=c( "-15", "-10", "-5","0","5","10"))


## ----calendar, cache = FALSE,fig.height=5, fig.width=7------------------------
plot(out,type = "calendar")


## ----sim-equiv, cache = FALSE, fig.height=5, fig.width=7----------------------
plot(out, type = "equiv", ylim = c(-5, 5))


## ----sim-equiv-no-stats, cache = FALSE, fig.height=5, fig.width=7-------------
plot(out, type = "equiv", show.stats =  FALSE)


## ----sim-equiv-reposition, cache = FALSE, fig.height=5, fig.width=7-----------
plot(out, type = "equiv", stats.pos = c(-19, 4.5), ylim = c(-5, 5))


## ----turnout-panelview-status, cache = FALSE, warning=FALSE, fig.height=10, fig.width=7----
panelview(turnout ~ policy_edr, data = turnout, 
          index = c("abb","year"), pre.post = TRUE, 
          by.timing = TRUE) 


## ----turnout-panelview-outcome, cache = FALSE, warning =FALSE, fig.height=5, fig.width=7----
panelview(turnout ~ policy_edr, data = turnout,
          index = c("abb","year"), type = "outcome", 
          main = "EDR Reform and Turnout")


## ----turnout_did, cache = TRUE------------------------------------------------
out0 <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor, 
               data = turnout, index = c("abb","year"), 
               se = TRUE, method = "gsynth",
               r = 0, CV = FALSE, force = "two-way", 
               nboots = 1000, seed = 02139)


## ----turnout_did_gap, fig.height=5, fig.width=7-------------------------------
plot(out0, type = "gap", xlim = c(-15, 5), ylim=c(-15, 10))


## ----turnout_est, cache = TRUE------------------------------------------------
out_turnout <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor, 
              data = turnout,  index = c("abb","year"), 
              se = TRUE, method = "gsynth", vartype = "parametric",
              r = c(0, 5), CV = TRUE, force = "two-way", 
              nboots = 1000, seed = 02139, keep.sims = TRUE)


## ----turnout-implied-weights--------------------------------------------------
dim(out_turnout$wgt.implied)
sort(out_turnout$wgt.implied[,8])


## ----turnout_gap, fig.height=5, fig.width=7-----------------------------------
plot(out_turnout, xlim = c(-10, 5), ylim=c(-10, 10))


## ----turnout-status-plot, fig.height=12, fig.width=7--------------------------
plot(out_turnout, type = "status",xlab = "Year", ylab = "State", main = "Treatment Status", 
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960, 
                   1968, 1976, 1984, 1992, 2000, 2008), xangle=10)


## ----turnout_counterfactual, fig.height=5, fig.width=7------------------------
plot(out_turnout, type = "counterfactual")


## ----turnout_gap2, fig.height=5, fig.width=7----------------------------------
plot(out_turnout, type = "counterfactual", id = "WI", main = "Wisconsin")


## ----turnout_box, fig.height=5, fig.width=7-----------------------------------
plot(out_turnout, type = "box", 
     xticklabels=c("-20", "-15", "-10", "-5","0","5","10"))


## ----turnout_calendar, fig.height=5, fig.width=7------------------------------
plot(out_turnout, type = "calendar", ylim = c(-15,15))


## ----turnout_F, message = FALSE, results = 'hide', fig.height=5, fig.width=7, warning=FALSE----
plot(out_turnout, type = "factors", xlab = "Year")


## ----turnout_L, message = FALSE, results = 'hide', fig.height=7, fig.width=7, warning=FALSE----
plot(out_turnout, type = "loadings")


## ----turnout-cumulative, fig.height=5, fig.width=7----------------------------
cumu.turnout <- estimand(out_turnout, "att.cumu", "event.time")
esplot(cumu.turnout, Period = "event.time",
       Estimate = "estimate", CI.lower = "ci.lo", CI.upper = "ci.hi",
       Count = "n_cells",
       main = "Cumulative Effect of EDR Reform on Turnout",
       ylab = "Cumulative ATT")


## ----create-unbalanced-data---------------------------------------------------
set.seed(123456)
turnout.ub <- turnout[-c(which(turnout$abb=="WY")[1:15], 
                         sample(1:nrow(turnout),50,replace=FALSE)),]


## ----turnout_ub_panelview_miss, cache = FALSE,fig.height=7, fig.width=7-------
panelview(turnout ~ policy_edr + policy_mail_in + policy_motor, 
          data = turnout.ub,  index = c("abb","year"), 
          pre.post = TRUE) 


## ----turnout_ub_est, cache = TRUE, message = FALSE----------------------------
out_ub <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor,
              data = turnout.ub,  index = c("abb","year"),
              se = TRUE, method = "gsynth",
              r = c(0, 5), CV = TRUE, force = "two-way",
              parallel = TRUE, cores = 16, min.T0 = 8,
              nboots = 1000, seed = 02139)


## ----turnout_ub_param, cache = TRUE, message = FALSE--------------------------
out_ub_param <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor,
                     data = turnout.ub, index = c("abb","year"),
                     se = TRUE, method = "gsynth", vartype = "parametric",
                     r = c(0, 5), CV = TRUE, force = "two-way",
                     parallel = TRUE, cores = 16, min.T0 = 8,
                     nboots = 1000, seed = 02139)


## ----turnout_ub_ci_compare----------------------------------------------------
ci_width <- function(out) mean(out$est.att[, "CI.upper"] -
                                out$est.att[, "CI.lower"], na.rm = TRUE)
data.frame(
  vartype  = c("bootstrap", "parametric"),
  CI_width = c(ci_width(out_ub), ci_width(out_ub_param))
)


## ----turnout_ub_panelview_miss2, fig.height=12, fig.width=7-------------------
plot(out_ub, type = "status",
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960, 
                   1968, 1976, 1984, 1992, 2000, 2008),
     xangle=10)


## ----turnout_ub_obs_2, fig.height=7, fig.width=7------------------------------
plot(out_ub, type = "status", xlab = "Year", ylab = "State",
     main = "Treatment Status", id = out_ub$id[out_ub$tr],
     xlim = c(1920,2012), 
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960,
                   1968, 1976, 1984, 1992, 2000, 2008))



## ----turnout_ub_gap, fig.height=5, fig.width=7--------------------------------
plot(out_ub, type = "gap", xlim = c(-10, 5), ylim = c(-10, 15))


## ----turnout_ub_gap_param, fig.height=5, fig.width=7--------------------------
plot(out_ub_param, type = "gap", xlim = c(-10, 5), ylim = c(-10, 15))


## ----bounded-vs-unbounded-fit, eval = TRUE, cache = TRUE, message = FALSE, results = "hide"----
# Unbounded gsynth (default v2.2.x behavior --- standard Xu 2017 estimator)
out.unbounded <- fect(Y ~ D, data = sim_gsynth, index = c("id", "time"),
                      method = "gsynth", r = 2,
                      se = TRUE, vartype = "parametric", nboots = 100,
                      CV = FALSE, parallel = FALSE)

# Bounded gsynth (new in v2.3.0 --- simplex projection)
out.bounded   <- fect(Y ~ D, data = sim_gsynth, index = c("id", "time"),
                      method = "gsynth", r = 2,
                      loading.bound = "simplex",
                      se = TRUE, vartype = "parametric", nboots = 100,
                      CV = FALSE, parallel = FALSE)


## ----bounded-vs-unbounded-att-------------------------------------------------
cat("Unbounded ATT =", round(out.unbounded$att.avg, 3), "\n")
cat("Bounded   ATT =", round(out.bounded$att.avg,   3), "\n")
cat("\ngamma.loading (CV-selected) =", round(out.bounded$gamma.loading, 4), "\n")
cat("\nPer-treated projection residuals (loading.proj.resid):\n")
print(round(out.bounded$loading.proj.resid, 3))


## ----bounded-vs-unbounded-gap, fig.width = 9, fig.height = 4, message = FALSE----
library(gridExtra)
gap.un  <- plot(out.unbounded, type = "gap", main = "Unbounded gsynth")
gap.bd  <- plot(out.bounded,   type = "gap", main = "Bounded gsynth (simplex)")
grid.arrange(gap.un, gap.bd, ncol = 2)


## ----bounded-vs-unbounded-overlap, fig.width = 9, fig.height = 4--------------
ov.un <- plot(out.unbounded, type = "loading.overlap",
              main = "Unbounded loadings")
ov.bd <- plot(out.bounded,   type = "loading.overlap",
              main = "Bounded loadings (projected)")
grid.arrange(ov.un, ov.bd, ncol = 2)


## ----bounded-vs-unbounded-weights, fig.width = 6.5, fig.height = 5, message = FALSE----
library(ggplot2)
W.un <- as.matrix(out.unbounded$wgt.implied)
W.bd <- as.matrix(out.bounded$wgt.implied)
df.w <- data.frame(
  unbounded = as.numeric(W.un),
  bounded   = as.numeric(W.bd)
)

ggplot(df.w, aes(x = unbounded, y = bounded)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0,
              color = "grey60", linetype = "dashed", linewidth = 0.4) +
  geom_point(alpha = 0.4, size = 1.4, color = "#3F6A99") +
  labs(x = "Unbounded weight (Moore-Penrose pseudo-inverse)",
       y = "Bounded weight (simplex-projected)",
       title = "Implicit donor weights: bounded vs unbounded",
       subtitle = "One point per (treated, control) pair; dashed line is y = x") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())


## ----cfe_nt_demo, eval=TRUE, cache=TRUE, message=FALSE, results='hide'--------
out.cfe.nt <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                   method = "cfe", force = "two-way",
                   time.component.from = "nevertreated",
                   Q.type = "linear",
                   se = FALSE, CV = TRUE, r = c(0, 5),
                   parallel = TRUE, cores = 16)


## ----cfe-nt-summary-----------------------------------------------------------
cat("CFE + nevertreated: r.cv =", out.cfe.nt$r.cv,
    ", ATT =", round(out.cfe.nt$att.avg, 3), "\n")


## ----cfe_nt_vs_gsynth, eval=TRUE, cache=TRUE, message=FALSE, results='hide'----
# Model 1: gsynth (pure IFE, r = 2)
out.gsynth.comp <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                        method = "gsynth", force = "two-way",
                        r = 2, se = FALSE, CV = FALSE)

# Model 2: CFE + nevertreated with r = 2 only (equivalent to gsynth)
out.cfe.nt.comp <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                        method = "cfe", force = "two-way",
                        time.component.from = "nevertreated",
                        r = 2, se = FALSE, CV = FALSE)

# Model 3: CFE + nevertreated with r = 2 and linear trend (overspecified)
out.cfe.nt.lin <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                       method = "cfe", force = "two-way",
                       time.component.from = "nevertreated",
                       Q.type = "linear", r = 2, se = FALSE, CV = FALSE)


## ----cfe_nt_mspe, eval=TRUE, cache=TRUE---------------------------------------
mspe.comp <- fect_mspe(list(gsynth_r2 = out.gsynth.comp,
                            CFE_r2 = out.cfe.nt.comp,
                            CFE_linear_r2 = out.cfe.nt.lin), seed = 1234)
print(mspe.comp$summary[, c("Model", "MSPE", "RMSE", "MAD")])

