##########################
# Factor-Augmented Approaches
##########################
rm(list = ls())
library(fect)
data(fect)

set.seed(1234)

##########################
# From: ./03-ife-mc.Rmd
##########################

##########################
# Interactive fixed effects (IFE)
##########################

out.ife <- fect(Y ~ D + X1 + X2, data = simdata2, index = c("id","time"),
          force = "two-way", method = "ife", CV = TRUE, r = c(0, 5),
          se = TRUE, cores = 8, nboots = 1000, parallel = TRUE)
print(out.ife)

plot(out.ife, main = "Estimated ATT (IFEct)")

##########################
# Matrix completion (MC)
##########################

out.mc <- fect(Y ~ D + X1 + X2, data = simdata2, index = c("id","time"),
          force = "two-way", method = "mc", CV = TRUE,
          se = TRUE, cores = 8, nboots = 1000, parallel = TRUE)
print(out.mc)

plot(out.mc, main = "Estimated ATT (MC)")

##########################
# Quick diagnostic check
##########################

plot(out.ife, type = "equiv", ylim = c(-4, 4),
     cex.legend = 0.6, main = "Testing Pre-Trend (IFEct)", cex.text = 0.8)

##########################
# Cumulative effects
##########################

out <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"),
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
out_mc <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"),
                        method = "mc", force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 200, vartype = 'bootstrap',
                        parallel = FALSE, keep.sims=TRUE)
plot(effect(out_mc))

## Jackknife cumulative effects
out_jack <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"),
                        method = "mc", force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 200, vartype = 'jackknife',
                        parallel = FALSE, keep.sims=TRUE)
plot(effect(out_jack))
