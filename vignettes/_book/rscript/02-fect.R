##########################
# Install Packages 
##########################
rm(list = ls())
install.packages("fect")
devtools::install_github("xuyiqing/fect")
installed.packages()["fect", "Version"]
devtools::install_github('xuyiqing/panelView')
install_all <- function(packages) {
  installed_pkgs <- installed.packages()[, "Package"]
  for (pkg in packages) {
    if (!pkg %in% installed_pkgs) {
      install.packages(pkg)
    }
  }
}
packages <- c("abind", "doParallel", "doRNG", "fixest", "foreach", "future", 
              "GGally", "ggplot2", "grid", "gridExtra", "Mass", 
              "panelView", "Rcpp")
install_all(packages)
library(fect)
data(fect)
ls()

##########################
# From: ./02-fect.Rmd
##########################

set.seed(1234)
rm(list = ls())

library(fect)
data(fect)
ls()

library(panelView)
panelview(Y ~ D, data = simdata, index = c("id","time"), 
  axis.lab = "time", xlab = "Time", ylab = "Unit", 
  gridOff = TRUE, by.timing = TRUE,
  background = "white", main = "Simulated Data: Treatment Status")

panelview(Y ~ D, data = simdata, index = c("id","time"), 
  axis.lab = "time", xlab = "Time", ylab = "Unit", 
  theme.bw = TRUE, type = "outcome", 
  main = "Simulated Data: Outcome")

out.fect <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                 method = "fe", force = "two-way")

plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y", 
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

out.fect <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"), 
  method = "fe", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200)

plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y", 
  cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8, stats = "F.p")

print(out.fect)

out.fect$est.att
out.fect$est.avg
out.fect$est.beta

out.ife <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"), 
          force = "two-way", method = "ife", CV = TRUE, r = c(0, 5), 
          se = TRUE, nboots = 200, parallel = TRUE)

plot(out.ife, main = "Estimated ATT (IFEct)")

out.mc <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"), 
          force = "two-way", method = "mc", CV = TRUE, 
          se = TRUE, nboots = 200, parallel = TRUE)

plot(out.mc, main = "Estimated ATT (MC)")

out.bal <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"), 
  balance.period = c(-3, 4), force = "two-way", method = "ife", 
  CV = FALSE, r = 2, se = TRUE, nboots = 200, parallel = TRUE)

plot(out.bal, main = "Estimated ATT (Balanced Sample)")

plot(out.bal, main = "Estimated ATT (Balanced Sample)", 
  color = "red", count.color = "blue")

simdata$Weight <- abs(rnorm(n = dim(simdata)[1]))
out.w <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"), 
  force = "two-way", method = "ife", W = 'Weight',
  CV = FALSE, r = 2, se = TRUE, nboots = 200, parallel = TRUE)

plot(out.w, main = "Estimated Weighted ATT")

plot(out.ife, type = "box", xlim = c(-15, 10))

plot(out.ife, type = "calendar", xlim = c(1, 35))

out.fect.p <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", parallel = TRUE, se = TRUE, CV = 0,
  nboots = 200, placeboTest = TRUE, placebo.period = c(-2, 0))

out.ife.p <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "ife",  r = 2, CV = 0,
  parallel = TRUE, se = TRUE,
  nboots = 200, placeboTest = TRUE, placebo.period = c(-2, 0))

out.mc.p <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "mc",  lambda = out.mc$lambda.cv, 
  CV = 0, parallel = TRUE, se = TRUE,
  nboots = 200, placeboTest = TRUE, placebo.period = c(-2, 0))

plot(out.fect.p, cex.text = 0.8, stats = c("placebo.p","equiv.p"), 
     main = "Estimated ATT (TWFE)")

plot(out.ife.p, ylab = "Effect of D on Y", main = "Estimated ATT (IFE)", 
     cex.text = 0.8, stats = c("placebo.p","equiv.p"))

plot(out.mc.p, cex.text = 0.8, stats = c("placebo.p","equiv.p"),
     main = "Estimated ATT (MC)")

plot(out.fect, type = "equiv", ylim = c(-4,4), 
     cex.legend = 0.6, main = "Testing Pre-Trend (FEct)", cex.text = 0.8)

plot(out.ife, type = "equiv", ylim = c(-4,4), 
     cex.legend = 0.6, main = "Testing Pre-Trend (IFEct)", cex.text = 0.8)

plot(out.mc, type = "equiv", ylim = c(-4,4),
     cex.legend = 0.6, main = "Testing Pre-Trend (MC)", cex.text = 0.8)

out.fect.loo <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"), 
  method = "fe", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200, loo = TRUE)
out.ife.loo <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"), 
  method = "ife", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200, loo = TRUE)
out.mc.loo <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"), 
  method = "mc", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200, loo = TRUE)

plot(out.fect.loo, type = "equiv", ylim = c(-4,4), loo = TRUE,
     cex.legend = 0.6, main = "Testing Pre-Trend LOO (FEct)", cex.text = 0.8)

plot(out.ife.loo, type = "equiv", ylim = c(-4,4), loo = TRUE,
     cex.legend = 0.6, main = "Testing Pre-Trend LOO (IFEct)", cex.text = 0.8)

plot(out.mc.loo, type = "equiv", ylim = c(-4,4), loo = TRUE, 
     cex.legend = 0.6, main = "Testing Pre-Trend LOO (MC)", cex.text = 0.8)

plot(out.fect, type = "exit", ylim = c(-2.5,4.5), main = "What Happens after the Treatment Switches Off?")

plot(out.ife, type = "exit", ylim = c(-2.5,4.5), main = "Exit Plot (IFE)")

plot(out.mc, type = "exit", ylim = c(-2.5,4.5), main = "Exit Plot (MC)")

out.fect.c <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", parallel = TRUE, se = TRUE, CV = 0, 
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))

out.ife.c <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "ife", r = 2, CV = 0, 
  parallel = TRUE, se = TRUE,
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))

out.mc.c <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "mc",  lambda = out.mc$lambda.cv, 
  CV = 0, parallel = TRUE, se = TRUE, 
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))

plot(out.fect.c, type = "exit", ylim = c(-2.5,4.5), 
          cex.text = 0.8, main = "Carryover Effects (FE)")

plot(out.ife.c, type = "exit", ylim = c(-2.5,4.5), 
          cex.text = 0.8, main = "Carryover Effects (IFE)")

plot(out.mc.c, type = "exit", ylim = c(-2.5,4.5), 
          cex.text = 0.8, main = "Carryover Effects (MC)")

out.ife.rm.test <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
  force = "two-way", method = "ife", r = 2, CV = 0, 
  parallel = TRUE, se = TRUE,  carryover.rm = 3,
  nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))# remove three periods

plot(out.ife.rm.test, cex.text = 0.8, stats.pos = c(5, 2))

out <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"),
                        method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 200, vartype = 'bootstrap',
                        parallel = FALSE, keep.sims=TRUE)
cumu.out <- effect(out)

print(cumu.out)
plot(cumu.out)

effect(out, cumu=FALSE)

effect(out, cumu=TRUE, id=c(101,102,103), period=c(1,5))

out_mc <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"),
                        method = "mc", force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 200, vartype = 'bootstrap',
                        parallel = FALSE, keep.sims=TRUE)
plot(effect(out_mc))

out_jack <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"),
                        method = "mc", force = "two-way", CV = TRUE, r = c(0, 5),
                        se = TRUE, nboots = 200, vartype = 'jackknife',
                        parallel = FALSE, keep.sims=TRUE)
plot(effect(out_jack))

simdata[,"FE3"] <- sample(c(1,2,3,4,5), size = dim(simdata)[1], replace = TRUE)
out.cfe <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"), 
  method = "cfe", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200,
  sfe = c("FE3"), cfe = list(c("id","time")))
plot(out.cfe)

out.poly <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"), 
  method = "polynomial", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200,
  degree = 2)
plot(out.poly)

plot(out.ife, show.count = FALSE)

plot(out.ife, type = "counterfactual")

plot(out.fect, type = 'status', axis.lab = "time", cex.axis  = 0.6)

plot(out.fect.p, type = 'status', axis.lab = "both", id = c(101:120), cex.axis  = 0.6)

plot(out.fect.c, type = 'status', axis.lab = "off", gridOff = TRUE)

plot(out.ife.rm.test, type = 'status', axis.lab = "off", gridOff = TRUE)

panelview(Y ~ D, data = simdata, index = c("id","time"), by.timing = TRUE,
  axis.lab = "time", xlab = "Time", ylab = "Unit", 
  background = "white", main = "Simulated Data: Treatment Status")

simdata.cohort <- get.cohort(data = simdata,D = 'D',index = c("id","time"))
print(table(simdata.cohort[,'Cohort']))

simdata.cohort2 <- get.cohort(data = simdata,D = 'D',index = c("id","time"),
                               entry.time = list(c(21,27),c(30,33)))
print(table(simdata.cohort2[,'Cohort']))

out.ife.g <- fect(Y ~ D + X1 + X2, data = simdata.cohort, index = c("id","time"), 
          force = "two-way", method = "ife", CV = TRUE, r = c(0, 5), 
          se = TRUE, nboots = 200, parallel = TRUE, group = 'Cohort')
out.ife.g.p <- fect(Y ~ D + X1 + X2, data = simdata.cohort, index = c("id","time"), 
          force = "two-way", method = "ife", CV = FALSE, 
          placeboTest = TRUE, placebo.period = c(-2,0), 
          se = TRUE, nboots = 200, parallel = TRUE, group = 'Cohort')

plot(out.ife.g, show.group = "Cohort:22", 
          xlim = c(-15, 10), ylim = c(-10, 10))
