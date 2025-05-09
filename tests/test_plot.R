
## Fixed Effects Counterfactuals Estimators

### Starting from v2.0.0, **gsynth*** is integrated into **fect**.

## Date: 1/18/2025

## ## Install from Github (development version)
# devtools::install_github('xuyiqing/panelView')
# devtools::install_github('xuyiqing/fect')


###################################################
##  Load data
###################################################
set.seed(1234)
data(fect)
ls()


##################################################
## Look at data
##################################################

# Here, we load simdata, a dataset of 200 units and 35 time periods.
# For further details, see Liu, Wang, Xu (2022).
library(panelView)
panelview(Y ~ D, data = simdata, index = c("id","time"),
          axis.lab = "time", xlab = "Time", ylab = "Unit",
          gridOff = TRUE, by.timing = TRUE,
          background = "white", main = "Simulated Data: Treatment Status")
panelview(Y ~ D, data = simdata, index = c("id","time"),
          axis.lab = "time", xlab = "Time", ylab = "Unit",
          theme.bw = TRUE, type = "outcome", main = "Simulated Data: Outcome")

# We load simgsynth, a dataset with 5 treated units, 45 control units, and 30 time periods.
# Unlike simdata, simgsynth does not have treatment reversal.
# For further details, see Xu(2017).
panelview(Y ~ D, data = simgsynth,  index = c("id","time"),
          pre.post = TRUE, axis.lab = "time")
panelview(Y ~ D, data = simgsynth,  index = c("id","time"), type = "outcome")

# We load turnout to show how to run gsynth for mutiperiod DID.
# For further details, see Xu(2017).
panelview(turnout ~ policy_edr, data = turnout,
          index = c("abb","year"), pre.post = TRUE,
          by.timing = TRUE)
panelview(turnout ~ policy_edr, data = turnout,
          index = c("abb","year"), type = "outcome",
          main = "EDR Reform and Turnout",
          by.group = TRUE)


##################################################
## FEct
##################################################

# Base FEct
out.fect <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                 method = "fe", force = "two-way")
plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

# FEct with uncertainty estimates
out.fect <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                 method = "fe", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200)
plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8, stats = "F.p")
print(out.fect)

# IFEct
out.ife <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                force = "two-way", method = "ife", CV = TRUE, r = c(0, 5),
                se = TRUE, nboots = 200, parallel = TRUE)
plot(out.ife, main = "Estimated ATT (IFEct)")

# MC
out.mc <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
               force = "two-way", method = "mc", CV = TRUE,
               se = TRUE, nboots = 200, parallel = TRUE)
plot(out.mc, main = "Estimated ATT (MC)")

# Balanced Treated Sample
# by specifying `balance.period`
out.bal <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                balance.period = c(-3, 4), force = "two-way", method = "ife",
                CV = FALSE, r = 2, se = TRUE, nboots = 200, parallel = TRUE)
plot(out.bal, main = "Estimated ATT (Balanced Sample)")


# Weighted Treatment Effect
# by specifying `W`
simdata$Weight <- abs(rnorm(n = dim(simdata)[1]))
out.w <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
              force = "two-way", method = "ife", W = 'Weight',
              CV = FALSE, r = 2, se = TRUE, nboots = 200, parallel = TRUE)
plot(out.w, main = "Estimated Weighted ATT")


##################################################
## Heterogeneous treatment effects
##################################################
# Box plots
plot(out.ife, type = "box", xlim = c(-15, 10))

# By calendar time
plot(out.ife, type = "calendar", xlim = c(1, 35))


##################################################
## Diagnostic tests
##################################################

# Placebo tests
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
plot(out.mc.p, vis = "none", cex.text = 0.8, stats = c("placebo.p","equiv.p"),
     main = "Estimated ATT (MC)")

# Equivalence test
plot(out.fect, type = "equiv", ylim = c(-4,4),
     cex.legend = 0.6, main = "Testing Pre-Trend (FEct)", cex.text = 0.8)
plot(out.ife, type = "equiv", ylim = c(-4,4),
     cex.legend = 0.6, main = "Testing Pre-Trend (IFEct)", cex.text = 0.8)
plot(out.mc, type = "equiv", ylim = c(-4,4),
     cex.legend = 0.6, main = "Testing Pre-Trend (MC)", cex.text = 0.8)

# Leave one out (slow)
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

# Exiting treatment and the carry-over test
plot(out.fect, type = "exit", ylim = c(-2.5,4.5), main = "What Happens after the Treatment Switches Off?")
plot(out.ife, type = "exit", ylim = c(-2.5,4.5), main = "Exit Plot (IFE)")
plot(out.mc, type = "exit", ylim = c(-2.5,4.5), main = "Exit Plot (MC)")

## Carry-over effects
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
plot(out.ife.rm.test, cex.text = 0.8, stats.pos = c(7, 2.5))


##################################################
## Complex Fixed Effects
##################################################

# Set additional fixed effect variable
simdata[,"FE3"] <- sample(c(1,2,3,4,5), size = dim(simdata)[1], replace = TRUE)

# Complex fixed effects
out.cfe <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                method = "cfe", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200,
                sfe = c("FE3"), cfe = list(c("id","time")))
plot(out.cfe)

# Polynomial
out.poly <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                 method = "polynomial", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200,
                 degree = 2)
plot(out.poly)


##################################################
## Visualization Options
##################################################

# Remove bars counting #treated units
plot(out.ife, count = FALSE)

# Show estimates using points
plot(out.ife, show.points = TRUE)
plot(out.ife.p, show.points = TRUE)
plot(out.ife.c, type = "exit", show.points = TRUE, stats.pos = c(2.3,4))

# Visualizing treatment status (and obs being used in placebo and carryover effect tests)
plot(out.fect, type = 'status', axis.lab = "time", cex.axis  = 0.6)
plot(out.fect.p, type = 'status', axis.lab = "both", id = c(101:120), cex.axis  = 0.6)
plot(out.fect.c, type = 'status', axis.lab = "off", gridOff = TRUE)
plot(out.ife.rm.test, type = 'status', axis.lab = "off", gridOff = TRUE)


##################################################
## Cohort Effects
##################################################
# Get cohort for the dataset
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
plot(out.ife.g.p, show.group = "Cohort:21", #could not find function "diagtest"
     xlim = c(-15, 10), ylim = c(-10, 10))

##################################################
## Gsynth
##################################################

# Estimation
system.time(
  out <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"),
              method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 5),
              se = TRUE, nboots = 1000, vartype = 'parametric',
              parallel = FALSE))
# Parallel computing
out <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"),
            method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 5),
            se = TRUE, nboots = 1000,vartype = 'parametric', parallel = TRUE, cores = 16)

out2 <- fect(Y ~ D + X1 + X2, data = simgsynth,  index = c("id","time"),
             method = "gsynth", force = "two-way",
             CV = TRUE, r = c(0, 5), se = TRUE,
             vartype = "jackknife",
             parallel = TRUE, cores = 10)


##################################################
## Visualization Options
##################################################

# Block DID
# Gap Plot
a <- plot(out) # by default, type = "gap"
print(a)
plot(out, theme.bw = FALSE)
plot(out, show.points = FALSE)
plot(out, type = "gap", ylim = c(-3,12), xlab = "Period",
     main = "Estimated ATT (FEct)")

plot(out, type = "status", yticklabels="0",
     xticklabels=c("5", "10", "15","20", "25", "30") )

# Box plot
plot(out, type = "box", xlab = "time",
     xticklabels=c("-19", "-15", "-10", "-5","0","5","10") )
plot(out, type = "box", xlim = c(-15, 10),
     xticklabels=c( "-15", "-10", "-5","0","5","10"))

# Calendar plot
plot(out,type = "calendar")

# Loadings plot
plot(out, type = "loadings")

# Factors plot
plot(out, type = "factors", xlab = "Time")

# Counterfactual plot
plot(out, type = "counterfactual")
plot(out, type = "counterfactual", raw = "all")
plot(out, type = "counterfactual", raw = "band")

# Equiv plot
plot(out, type = "equiv", ylim = c(-5, 5))
plot(out, type = "equiv", show.stats =  FALSE)
plot(out, type = "equiv", stats.pos = c(-19, 4.5), ylim = c(-5, 5))


##################################################
### Comparison with Other Estimators
##################################################

# IFEct
out.ife <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"),
                force = "two-way", method = "ife", CV = TRUE, r = c(0, 5),
                se = TRUE, nboots = 200, parallel = TRUE)
plot(out.ife, main = "Estimated ATT (EM)")

# MC
out.mc <- fect(Y ~ D + X1 + X2, data = simgsynth,
               index = c("id","time"),
               force = "two-way", method = "mc", CV = TRUE,
               se = TRUE, nboots = 200, parallel = TRUE)
plot(out.mc, main = "Estimated ATT (MC)")

# Staggered DID
out0 <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor,
             data = turnout, index = c("abb","year"),
             se = TRUE, method = "gsynth",
             r = 0, CV = FALSE, force = "two-way",
             nboots = 1000, seed = 02139)
plot(out0, type = "gap", xlim = c(-15, 15))


# Estimation w/ factors
out_turnout <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor,
                    data = turnout,  index = c("abb","year"),
                    se = TRUE, method = "gsynth",
                    r = c(0, 5), CV = TRUE, force = "two-way",
                    nboots = 1000, seed = 02139)


# Implied weights
dim(out_turnout$wgt.implied)
sort(out_turnout$wgt.implied[,8])

plot(out_turnout, xlim = c(-10, 5), ylim=c(-15,15))
plot(out_turnout, type = "status",xlab = "Year", ylab = "State", main = "Treatment Status",
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960,
                   1968, 1976, 1984, 1992, 2000, 2008), xangle=10)
plot(out_turnout, type = "gap", id = "WI", main = "Wisconsin")
plot(out_turnout, type = "box",
     xticklabels=c("-20", "-15", "-10", "-5","0","5","10"))
plot(out_turnout, type = "calendar", ylim = c(-15,15))
plot(out_turnout, type = "factors", xlab = "Year")
plot(out_turnout, type = "loadings")


# Unbalanced panels
set.seed(123456)
turnout.ub <- turnout[-c(which(turnout$abb=="WY")[1:15],
                         sample(1:nrow(turnout),50,replace=FALSE)),]
panelview(turnout ~ policy_edr + policy_mail_in + policy_motor,
          data = turnout.ub,  index = c("abb","year"),
          pre.post = TRUE)



# Estimation
out_ub <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor,
               data = turnout.ub,  index = c("abb","year"),
               se = TRUE, method = "gsynth",
               r = c(0, 5), CV = TRUE, force = "two-way",
               parallel = TRUE, min.T0 = 8,
               nboots = 1000, seed = 02139)
plot(out_ub, type = "status",
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960,
                   1968, 1976, 1984, 1992, 2000, 2008), xangle=10)
plot(out_ub, type = "status", xlab = "Year", ylab = "State",
     main = "Treatment Status", id = out_ub$id[out_ub$tr],
     xlim = c(1920,2012),
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960,
                   1968, 1976, 1984, 1992, 2000, 2008))
plot(out_ub, type = "gap", ylim = c(-10, 20))


out.mc_ub <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor,
                  min.T0 = 8, data = turnout.ub,
                  index = c("abb","year"), method = "mc",
                  se = TRUE, nboots = 1000, seed = 02139)
plot(out.mc_ub, main = "Estimated ATT (MC)", ylim = c(-10, 20))


out.ife.rm.test <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                        force = "two-way", method = "ife", r = 2, CV = 0,
                        parallel = TRUE, se = TRUE,  carryover.rm = 5,
                        nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))# remove three periods

plot(out.ife.rm.test, cex.text = 0.8, stats.pos = c(7, 2.5))

out.ife.rm.test <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                        force = "two-way", method = "ife", r = 2, CV = 0,
                        parallel = TRUE, se = TRUE,  carryover.rm = 3,
                        nboots = 200, carryoverTest = TRUE, carryover.period = c(1, 3))# remove three periods# remove three periods
