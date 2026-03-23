##############################
# 07-gsynth.R
# Generated from 07-gsynth.Rmd
##############################
rm(list = ls())
set.seed(1234)

## --- load-packages ---
library(fect)
data(fect)
ls()

## --- head-sim_gsynth ---
head(sim_gsynth)

## --- sim-panelview-status ---
library(panelView)
panelview(Y ~ D, data = sim_gsynth,  index = c("id","time"), pre.post = TRUE)

## --- sim-panelview-outcome ---
panelview(Y ~ D, data = sim_gsynth,  index = c("id","time"), type = "outcome")

## --- sim2_onecore ---
system.time(
out <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
            method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 5),
            se = TRUE, nboots = 1000, vartype = 'parametric',
            parallel = FALSE))

## --- print-results ---
## print(out)
## out$est.att
## out$est.avg
## out$beta

## --- sim2 ---
system.time(
out <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"), method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 5), se = TRUE, nboots = 1000,vartype = 'parametric', parallel = TRUE, cores = 16)
)

## --- simJack ---
out2 <- fect(Y ~ D + X1 + X2, data = sim_gsynth,  index = c("id","time"),
               method = "gsynth", force = "two-way",
               CV = TRUE, r = c(0, 5), se = TRUE,
               vartype = "jackknife",
               parallel = TRUE, cores = 8)

## --- sim_gap1 ---
a <- plot(out) # by default, type = "gap"
print(a)

## --- sim_gap1a ---
plot(out, theme.bw = FALSE)

## --- sim-gap-connected ---
plot(out, connected = TRUE)

## --- sim-gap-line ---
plot(out, connected = TRUE, show.points = FALSE)

## --- sim_gap2 ---
plot(out, type = "gap", ylim = c(-6,12), xlab = "Period",
     main = "Estimated ATT (Gsynth)")

## --- sim-counterfactual ---
plot(out, type = "counterfactual")

## --- sim-counterfactual-all ---
plot(out, type = "counterfactual", raw = "all")

## --- sim-counterfactual-band ---
plot(out, type = "counterfactual", raw = "band")

## --- sim_status ---
plot(out, type = "status", yticklabels="0",
     xticklabels=c("5", "10", "15","20", "25", "30") )

## --- sim_L ---
plot(out, type = "loadings")

## --- sim_F ---
plot(out, type = "factors", xlab = "Time")

## --- sim_box ---
plot(out, type = "box", xlab = "time",
     xticklabels=c("-19", "-15", "-10", "-5","0","5","10") )

## --- sim_box2 ---
## plot(out, type = "box", xlim = c(-15, 10),
##      xticklabels=c( "-15", "-10", "-5","0","5","10"))

## --- calendar ---
plot(out,type = "calendar")

## --- sim-equiv ---
plot(out, type = "equiv", ylim = c(-5, 5))

## --- sim-equiv-no-stats ---
plot(out, type = "equiv", show.stats =  FALSE)

## --- sim-equiv-reposition ---
plot(out, type = "equiv", stats.pos = c(-19, 4.5), ylim = c(-5, 5))

## --- turnout-panelview-status ---
panelview(turnout ~ policy_edr, data = turnout,
          index = c("abb","year"), pre.post = TRUE,
          by.timing = TRUE)

## --- turnout-panelview-outcome ---
panelview(turnout ~ policy_edr, data = turnout,
          index = c("abb","year"), type = "outcome",
          main = "EDR Reform and Turnout",
          by.group = TRUE)

## --- turnout_did ---
out0 <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor,
               data = turnout, index = c("abb","year"),
               se = TRUE, method = "gsynth",
               r = 0, CV = FALSE, force = "two-way",
               nboots = 1000, seed = 02139)

## --- turnout_did_gap ---
plot(out0, type = "gap", xlim = c(-15, 5), ylim=c(-15, 10))

## --- turnout_est ---
out_turnout <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor,
              data = turnout,  index = c("abb","year"),
              se = TRUE, method = "gsynth", vartype = "parametric",
              r = c(0, 5), CV = TRUE, force = "two-way",
              nboots = 1000, seed = 02139, keep.sims = TRUE)

## --- turnout-implied-weights ---
dim(out_turnout$wgt.implied)
sort(out_turnout$wgt.implied[,8])

## --- turnout_gap ---
plot(out_turnout, xlim = c(-10, 5), ylim=c(-15, 10))

## --- turnout-status-plot ---
plot(out_turnout, type = "status",xlab = "Year", ylab = "State", main = "Treatment Status",
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960,
                   1968, 1976, 1984, 1992, 2000, 2008), xangle=10)

## --- turnout_counterfactual ---
plot(out_turnout, type = "counterfactual")

## --- turnout_gap2 ---
plot(out_turnout, type = "counterfactual", id = "WI", main = "Wisconsin")

## --- turnout_box ---
plot(out_turnout, type = "box",
     xticklabels=c("-20", "-15", "-10", "-5","0","5","10"))

## --- turnout_calendar ---
plot(out_turnout, type = "calendar", ylim = c(-15,15))

## --- turnout_F ---
plot(out_turnout, type = "factors", xlab = "Year")

## --- turnout_L ---
plot(out_turnout, type = "loadings")

## --- create-unbalanced-data ---
set.seed(123456)
turnout.ub <- turnout[-c(which(turnout$abb=="WY")[1:15],
                         sample(1:nrow(turnout),50,replace=FALSE)),]

## --- turnout_ub_panelview_miss ---
panelview(turnout ~ policy_edr + policy_mail_in + policy_motor,
          data = turnout.ub,  index = c("abb","year"),
          pre.post = TRUE)

## --- turnout_ub_est ---
out_ub <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor,
              data = turnout.ub,  index = c("abb","year"),
              se = TRUE, method = "gsynth",
              r = c(0, 5), CV = TRUE, force = "two-way",
              parallel = TRUE, min.T0 = 8,
              nboots = 1000, seed = 02139)

## --- turnout_ub_panelview_miss2 ---
plot(out_ub, type = "status",
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960,
                   1968, 1976, 1984, 1992, 2000, 2008),
     xangle=10)

## --- turnout_ub_obs_2 ---
plot(out_ub, type = "status", xlab = "Year", ylab = "State",
     main = "Treatment Status", id = out_ub$id[out_ub$tr],
     xlim = c(1920,2012),
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960,
                   1968, 1976, 1984, 1992, 2000, 2008))

## --- turnout_ub_gap ---
plot(out_ub, type = "gap", ylim = c(-10, 20))

## --- cfe_nt_demo ---
out.cfe.nt <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                   method = "cfe", force = "two-way",
                   time.component.from = "nevertreated",
                   Q.type = "linear",
                   se = FALSE, CV = TRUE, r = c(0, 5),
                   parallel = TRUE)

## --- cfe-nt-summary ---
cat("CFE + nevertreated: r.cv =", out.cfe.nt$r.cv,
    ", ATT =", round(out.cfe.nt$att.avg, 3), "\n")

## --- cfe_nt_vs_gsynth ---
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

## --- cfe_nt_mspe ---
mspe.comp <- fect_mspe(list(gsynth_r2 = out.gsynth.comp,
                            CFE_r2 = out.cfe.nt.comp,
                            CFE_linear_r2 = out.cfe.nt.lin), seed = 1234)
print(mspe.comp$summary[, c("Model", "MSPE", "RMSE", "MAD")])
