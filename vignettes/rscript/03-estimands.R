## ----.common, include = FALSE-------------------------------------------------
source("_common.R")


## ----setup-estimand, eval = TRUE, message = FALSE, results = 'hide'-----------
library(dplyr)

set.seed(1)
N <- 60; TT <- 25
df <- expand.grid(id = 1:N, time = 1:TT)
treat_start <- sample(c(NA, 8:18), N, replace = TRUE)
df$D <- ifelse(is.na(treat_start[df$id]) | df$time < treat_start[df$id],
               0, 1)
## Higher intercept (2.0) and smaller noise (sd = 0.1) keep Y safely
## positive throughout, so the v2.4.2+ cell-drop hard-error in
## log.att / aptt does not fire on benign bootstrap noise.
df$Y <- exp(2.0 + 0.05 * df$time + 0.3 * df$D + rnorm(nrow(df), sd = 0.1))

fit <- fect(Y ~ D, data = df, index = c("id", "time"),
            method = "fe", force = "two-way",
            se = TRUE, nboots = 200, parallel = FALSE,
            keep.sims = TRUE)


## ----est-att-event-time, eval = TRUE------------------------------------------
estimand(fit, "att", "event.time") |>
  head(8)


## ----plot-att-event-time, eval = TRUE, fig.width = 6, fig.height = 4.5--------
est <- estimand(fit, "att", "event.time")
esplot(est, Period = "event.time",
       Estimate = "estimate", CI.lower = "ci.lo", CI.upper = "ci.hi",
       Count = "n_cells",
       main = "Per-event-time ATT")


## ----est-cumu-event-time, eval = TRUE-----------------------------------------
estimand(fit, "att.cumu", "event.time") |>
  head(8)


## ----plot-cumu-event-time, eval = TRUE, fig.width = 6, fig.height = 4.5-------
cumu <- estimand(fit, "att.cumu", "event.time")
esplot(cumu, Period = "event.time",
       Estimate = "estimate", CI.lower = "ci.lo", CI.upper = "ci.hi",
       Count = "n_cells",
       main = "Cumulative ATT",
       ylab = "Cumulative Effect")


## ----est-cumu-overall, eval = TRUE--------------------------------------------
estimand(fit, "att.cumu", "overall", window = c(1, 5))


## ----est-aptt, eval = TRUE----------------------------------------------------
estimand(fit, "aptt", "event.time") |>
  head(8)


## ----plot-aptt-event-time, eval = TRUE, fig.width = 6, fig.height = 4.5-------
aptt <- estimand(fit, "aptt", "event.time")
esplot(aptt, Period = "event.time",
       Estimate = "estimate", CI.lower = "ci.lo", CI.upper = "ci.hi",
       Count = "n_cells",
       main = "APTT by event time",
       ylab = "Average Proportional TE on the Treated")


## ----est-log-att, eval = TRUE-------------------------------------------------
estimand(fit, "log.att", "event.time") |>
  head(8)


## ----plot-log-att-event-time, eval = TRUE, fig.width = 6, fig.height = 4.5----
log_att <- estimand(fit, "log.att", "event.time")
esplot(log_att, Period = "event.time",
       Estimate = "estimate", CI.lower = "ci.lo", CI.upper = "ci.hi",
       Count = "n_cells",
       main = "Log-scale ATT",
       ylab = "Mean log(Y) − log(Y0)")


## ----setup-placebo, eval = TRUE, message = FALSE, results = 'hide'------------
## Re-fit on the same DGP with placeboTest = TRUE so that the dispatcher
## can recompute APTT / log-ATT at the placebo cells.
fit_placebo <- fect(Y ~ D, data = df, index = c("id", "time"),
                    method = "fe", force = "two-way",
                    se = TRUE, nboots = 200, parallel = FALSE,
                    keep.sims = TRUE,
                    placeboTest = TRUE, placebo.period = c(-2, 0))


## ----plot-aptt-placebo, eval = TRUE, fig.width = 6, fig.height = 4.5----------
aptt_placebo <- estimand(fit_placebo, "aptt", "event.time",
                         test = "placebo")
esplot(aptt_placebo, Period = "event.time",
       Estimate = "estimate", CI.lower = "ci.lo", CI.upper = "ci.hi",
       Count = "n_cells",
       main = "APTT placebo (pre-treatment)",
       ylab = "Average Proportional TE on the Treated")


## ----plot-log-att-placebo, eval = TRUE, fig.width = 6, fig.height = 4.5-------
log_att_placebo <- estimand(fit_placebo, "log.att", "event.time",
                            test = "placebo")
esplot(log_att_placebo, Period = "event.time",
       Estimate = "estimate", CI.lower = "ci.lo", CI.upper = "ci.hi",
       Count = "n_cells",
       main = "Log-ATT placebo (pre-treatment)",
       ylab = "Mean log(Y) − log(Y0)")


## ----setup-carryover, eval = TRUE, message = FALSE, results = 'hide'----------
## Build a panel with treatment reversals for the carryover demo.
set.seed(2)
df_rev <- df
treat_end <- pmin(treat_start[df_rev$id] + sample(5:10, N, replace = TRUE),
                  TT + 1L)
df_rev$D <- ifelse(is.na(treat_start[df_rev$id]) |
                   df_rev$time < treat_start[df_rev$id] |
                   df_rev$time >= treat_end[df_rev$id],
                   0, 1)
df_rev$Y <- exp(2.0 + 0.05 * df_rev$time + 0.3 * df_rev$D +
                rnorm(nrow(df_rev), sd = 0.1))


## ----panelview-carryover, eval = TRUE, fig.width = 7, fig.height = 4.5--------
panelView::panelview(Y ~ D, data = df_rev, index = c("id", "time"),
                     by.timing = TRUE, axis.lab = "time",
                     main = "Treatment-reversal panel for carryover demo")


## ----fit-carryover, eval = TRUE, message = FALSE, results = 'hide'------------
fit_carry <- fect(Y ~ D, data = df_rev, index = c("id", "time"),
                  method = "fe", force = "two-way",
                  se = TRUE, nboots = 200, parallel = FALSE,
                  keep.sims = TRUE,
                  carryoverTest = TRUE, carryover.period = c(1, 2))


## ----plot-aptt-carryover, eval = TRUE, fig.width = 6, fig.height = 4.5--------
aptt_carry <- estimand(fit_carry, "aptt", "event.time",
                       test = "carryover")
esplot(aptt_carry, Period = "event.time",
       Estimate = "estimate", CI.lower = "ci.lo", CI.upper = "ci.hi",
       Count = "n_cells",
       main = "APTT carryover (post-reversal)",
       ylab = "Average Proportional TE on the Treated")


## ----plot-log-att-carryover, eval = TRUE, fig.width = 6, fig.height = 4.5-----
log_att_carry <- estimand(fit_carry, "log.att", "event.time",
                          test = "carryover")
esplot(log_att_carry, Period = "event.time",
       Estimate = "estimate", CI.lower = "ci.lo", CI.upper = "ci.hi",
       Count = "n_cells",
       main = "Log-ATT carryover (post-reversal)",
       ylab = "Mean log(Y) − log(Y0)")


## ----est-window-overall, eval = TRUE------------------------------------------
estimand(fit, "att", "overall", window = c(1, 5))


## ----imp-out-glance, eval = TRUE----------------------------------------------
po <- imputed_outcomes(fit)
head(po)


## ----imp-out-rep, eval = TRUE-------------------------------------------------
po_rep <- imputed_outcomes(fit, replicates = TRUE)
nrow(po_rep) == nrow(po) * 200    # one row per (cell, replicate)


## ----custom-estimand, eval = TRUE---------------------------------------------
po |>
  group_by(event.time) |>
  summarise(sd_eff = sd(eff, na.rm = TRUE),
            n     = dplyr::n(),
            .groups = "drop") |>
  head(8)


## ----keep-sims-error, eval = FALSE--------------------------------------------
# # No bootstrap/jackknife results available. Choose keep.sims = TRUE in fect().


## ----parametric-att, eval = FALSE---------------------------------------------
# fit_para <- fect(Y ~ D, data = sim_linear, index = c("id", "time"),
#                  method = "ife", force = "two-way",
#                  r = 2, CV = FALSE, se = TRUE, nboots = 200,
#                  keep.sims = TRUE,
#                  vartype = "parametric",
#                  time.component.from = "nevertreated",
#                  parallel = FALSE)
# 
# est <- estimand(fit_para, "att", "event.time")
# head(est)
# #>   event.time     estimate        se      ci.lo     ci.hi n_cells    vartype
# #> 1        -39 -0.019958112 0.2427538 -0.4957469 0.4558307      80 parametric
# #> 2        -38 -0.006695018 0.1693391 -0.3385935 0.3252035      80 parametric
# #> ...

