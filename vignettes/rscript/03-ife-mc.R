##########################
# Factor-Augmented Methods
##########################
rm(list = ls())
library(fect)
data(fect)

set.seed(1234)

##########################
# From: ./03-ife-mc.Rmd
##########################

##########################
# Interactive fixed effects
##########################

out.ife <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
          force = "two-way", method = "ife", CV = TRUE, r = c(0, 5),
          se = TRUE, cores = 8, nboots = 1000, parallel = TRUE)
print(out.ife)

plot(out.ife, main = "Estimated ATT (IFEct)")

##########################
# Matrix completion
##########################

out.mc <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
          force = "two-way", method = "mc", CV = TRUE,
          se = TRUE, cores = 8, nboots = 1000, parallel = TRUE)
print(out.mc)

plot(out.mc, main = "Estimated ATT (MC)")

##########################
# Cross-validation
##########################

out.cv <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
               method = "ife", CV = TRUE, r = c(0, 5),
               se = FALSE, parallel = TRUE)
cat("Selected r:", out.cv$r.cv, "\n")

## Comparing cv.method choices
out.all <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                method = "ife", CV = TRUE, r = c(0, 5),
                cv.method = "all_units", se = FALSE, parallel = TRUE)

out.tr <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
               method = "ife", CV = TRUE, r = c(0, 5),
               cv.method = "treated_units", se = FALSE, parallel = TRUE)

cat("cv.method = 'all_units':     r.cv =", out.all$r.cv, "\n")
cat("cv.method = 'treated_units': r.cv =", out.tr$r.cv, "\n")

## Scoring criteria
out.mspe <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                 method = "ife", CV = TRUE, r = c(0, 5),
                 criterion = "mspe", se = FALSE, parallel = TRUE)

out.pc <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
               method = "ife", CV = TRUE, r = c(0, 5),
               criterion = "pc", se = FALSE, parallel = TRUE)

cat("criterion = 'mspe': r.cv =", out.mspe$r.cv, "\n")
cat("criterion = 'pc':   r.cv =", out.pc$r.cv, "\n")

##########################
# Diagnostics
##########################

## Placebo tests
out.ife.p <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "ife",  r = 2, CV = 0,
  parallel = TRUE, se = TRUE,
  nboots = 200, placeboTest = TRUE, placebo.period = c(-2, 0))

out.mc.p <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "mc",  lambda = out.mc$lambda.cv,
  CV = 0, parallel = TRUE, se = TRUE,
  nboots = 200, placeboTest = TRUE, placebo.period = c(-2, 0))

plot(out.ife.p, ylab = "Effect of D on Y", main = "Estimated ATT (IFE)",
     cex.text = 0.8, stats = c("placebo.p","equiv.p"))

plot(out.mc.p, cex.text = 0.8, stats = c("placebo.p","equiv.p"),
     main = "Estimated ATT (MC)")

## LOO pre-trend test
# (computationally expensive, set eval=FALSE in Rmd)
# out.ife.loo <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
#   method = "ife", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200, loo = TRUE)
# out.mc.loo <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
#   method = "mc", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200, loo = TRUE)

## Equivalence test
plot(out.ife, type = "equiv", ylim = c(-4,4),
     cex.legend = 0.6, main = "Testing Pre-Trend (IFEct)", cex.text = 0.8)

plot(out.mc, type = "equiv", ylim = c(-4,4),
     cex.legend = 0.6, main = "Testing Pre-Trend (MC)", cex.text = 0.8)

## Tests for (no) carryover effects
out.ife.c <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "ife", r = 2, CV = 0,
  parallel = TRUE, se = TRUE,
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))

out.mc.c <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "mc",  lambda = out.mc$lambda.cv,
  CV = 0, parallel = TRUE, se = TRUE,
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))

plot(out.ife.c, type = "exit", ylim = c(-2.5,4.5),
          cex.text = 0.8, main = "Carryover Effects (IFE)")

plot(out.mc.c, type = "exit", ylim = c(-2.5,4.5),
          cex.text = 0.8, main = "Carryover Effects (MC)")

## Carryover removal
out.ife.rm.test <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "ife", r = 2, CV = 0,
  parallel = TRUE, se = TRUE,  carryover.rm = 3,
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))

plot(out.ife.rm.test, cex.text = 0.8, stats.pos = c(5, 2))

##########################
# Cumulative effects
##########################

out <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                        method = "ife", factors.from = "nevertreated",
                        force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 200, vartype = 'bootstrap',
                        parallel = FALSE, keep.sims=TRUE)
cumu.out <- effect(out)

print(cumu.out)
plot(cumu.out)

effect(out, cumu=FALSE)

effect(out, cumu=TRUE, id=c(101,102,103), period=c(1,5))

## MC cumulative effects
out_mc <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                        method = "mc", force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 200, vartype = 'bootstrap',
                        parallel = FALSE, keep.sims=TRUE)
plot(effect(out_mc))

## Jackknife cumulative effects
out_jack <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id","time"),
                        method = "mc", force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 200, vartype = 'jackknife',
                        parallel = FALSE, keep.sims=TRUE)
plot(effect(out_jack))
