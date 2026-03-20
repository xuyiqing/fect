##############################
# 02-fect.R
# Generated from 02-fect.Rmd
##############################
rm(list = ls())
set.seed(1234)

library(fect)
data(fect)

## --- panelview-treatment ---
library(panelView)
panelview(Y ~ D, data = sim_base, index = c("id","time"),
  axis.lab = "time", xlab = "Time", ylab = "Unit",
  gridOff = TRUE, by.timing = TRUE,
  background = "white", main = "Simulated Data: Treatment Status")

## --- panelview-outcome ---
panelview(Y ~ D, data = sim_base, index = c("id","time"),
  axis.lab = "time", xlab = "Time", ylab = "Unit",
  theme.bw = TRUE, type = "outcome",
  main = "Simulated Data: Outcome")

## --- simdata_fect_nose ---
out.fect <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
                 method = "fe", force = "two-way")

## --- fect_plot_nose ---
plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y",
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## --- simdata_fect ---
out.fect <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  method = "fe", force = "two-way", se = TRUE,
  cores = 8, parallel = TRUE, nboots = 1000)

## --- fect_plot_nse ---
plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y",
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8, stats = "F.p")

## --- exit_fect ---
plot(out.fect, type = "exit", main = "Exit Plot (FEct)")

## --- print-fect ---
print(out.fect)

## --- extract-estimates ---
## out.fect$est.att
## out.fect$est.avg
## out.fect$beta

## --- extract-bootstrap ---
## out.fect$eff.boot

## --- fect_placebo ---
out.fect.placebo <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  force = "two-way", method = "fe",
  se = TRUE, cores = 8, nboots = 200, parallel = TRUE,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.fect.placebo, cex.text = 0.8)

## --- fect_carryover ---
out.fect.carry <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  force = "two-way", method = "fe",
  se = TRUE, cores = 8, nboots = 200, parallel = TRUE,
  carryoverTest = TRUE, carryover.period = c(1, 3))

## --- fect_carryover_plot ---
plot(out.fect.carry, type = "exit", cex.text = 0.8, main = "Carryover Effects (FEct)")

## --- fect_loo ---
out.fect.loo <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  method = "fe", force = "two-way", se = TRUE, loo = TRUE,
  cores = 8, parallel = TRUE, nboots = 200)

## --- plot-gap-loo ---
plot(out.fect.loo,main = "Estimated ATT (FEct) -- LOO",
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## --- cumu_effect ---
out <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                        method = "ife", time.component.from = "nevertreated",
                        force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 200, vartype = 'bootstrap',
                        parallel = FALSE, keep.sims=TRUE)
cumu.out <- effect(out)

## --- cumu_effect_plot ---
print(cumu.out)
plot(cumu.out)

## --- cumu_effect_byperiod ---
effect(out, cumu=FALSE)

## --- cumu_effect_subset ---
effect(out, cumu=TRUE, id=c(101,102,103), period=c(1,5))

## --- effect-mc ---
out_mc <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                        method = "mc", force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 200, vartype = 'bootstrap',
                        parallel = FALSE, keep.sims=TRUE)
plot(effect(out_mc))

## --- effect-jackknife ---
out_jack <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                        method = "mc", force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 200, vartype = 'jackknife',
                        parallel = FALSE, keep.sims=TRUE)
plot(effect(out_jack))

## --- simdata_bal ---
out.bal <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  balance.period = c(-3, 4), force = "two-way", method = "ife",
  CV = FALSE, r = 2, se = TRUE, nboots = 200, parallel = TRUE)

## --- plot-balanced-att ---
plot(out.bal, main = "Estimated ATT (Balanced Sample)")

## --- plot-balanced-custom ---
plot(out.bal, main = "Estimated ATT (Balanced Sample)",
  post.color = "red", count.color = "blue")

## --- simdata_panelview_cohort ---
panelview(Y ~ D, data = sim_base, index = c("id","time"), by.timing = TRUE,
  axis.lab = "time", xlab = "Time", ylab = "Unit",
  background = "white", main = "Simulated Data: Treatment Status")

## --- get_cohort ---
sim_base.cohort <- get.cohort(data = sim_base,D = 'D',index = c("id","time"))
print(table(sim_base.cohort[,'Cohort']))

## --- get_cohort2 ---
sim_base.cohort2 <- get.cohort(data = sim_base,D = 'D',index = c("id","time"),
                               entry.time = list(c(21,27),c(30,33)))
print(table(sim_base.cohort2[,'Cohort']))

## --- simdata_fe_cohort ---
out.fe.g <- fect(Y ~ D + X1 + X2, data = sim_base.cohort, index = c("id","time"),
          force = "two-way", method = "fe",
          se = TRUE, nboots = 200, parallel = TRUE, group = 'Cohort')

## --- cohort_plot1 ---
plot(out.fe.g, show.group = "Cohort:22",
          xlim = c(-15, 10), ylim = c(-10, 10))

## --- simdata_w ---
sim_base$Weight <- abs(rnorm(n = dim(sim_base)[1]))
out.w <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  force = "two-way", method = "ife", W = 'Weight',
  CV = FALSE, r = 2, se = TRUE, nboots = 200, parallel = TRUE)

## --- plot-weighted-att ---
plot(out.w, main = "Estimated Weighted ATT")
