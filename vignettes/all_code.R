# From: ./01-start.Rmd
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

# From: ./02-fect.Rmd
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

#  devtools:: install_github("xuyiqing/paneltools" if not already installed
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

# From: ./03-plots.Rmd
# load libraries and data
library(ggplot2)
library(panelView)
library(fect)
data(fect)
ls()

out <- fect(Y = "general_sharetotal_A_all", 
            D = "cand_A_all", 
            X = c("cand_H_all", "cand_B_all"), 
            index = c("district_final", "cycle"), 
            data = gs2020, method = "fe", 
            force = "two-way", se = TRUE, 
            parallel = TRUE, nboots = 1000)

out.hh <- fect(nat_rate_ord ~ indirect, 
               data = hh2019,
               index = c("bfs","year"),
               method = 'fe', se = TRUE, 
               parallel = TRUE, nboots = 1000,
               keep.sims = TRUE)

plot(out) # the effect co-ethnic mobilization
plot(out.hh) # the effect of indirect democrazy on naturalization rate

plot(out, start0 = TRUE,  # Shift time so treatment begins at 0
     main = "Custom Starting Period")

plot(out, 
     xlim = c(-10, 1),  # only show time periods -8 to 1
     ylim = c(-0.15, 0.30),  # set y-range
     xlab = "Custom Time Axis", # x-axis label
     ylab = "Estimated ATT", # y-axis label
     xangle = 90, # rotate x-axis labels by 90°
     xbreaks = seq(-10, 1, by = 2), 
     gridOff = TRUE,
     # Label x-axis from -12 to 1 with a break of 2
     main = "Axis and Legend Customization")

plot(out, plot.ci = "0.9", 
     main = "90% confidence intervals")

plot(out, 
     ylim = c(-0.15, 0.3), # set yrange
     theme.bw = FALSE,      # Change the color theme
     cex.main = 1.25,     # Scale for the main title
     cex.axis = 1.2,      # Axis tick label size
     cex.lab = 1.2,       # Axis label size
     cex.legend = 1,    # Legend text size
     cex.text = 1.2,        # Annotation text size
     main = "Text and Theme Customization")

plot(out, 
     preset = "vibrant", # Use vibrant colors
     main = "Vibrant Preset Colors: Grumbach and Sahn (2020)")
plot(out.hh, 
     preset = "vibrant", # Use vibrant colors
     main = "Vibrant Preset Colors: Hainmueller and Hangartner (2019)")

plot(out.hh, 
     preset = "vibrant", # Use vibrant colors
     color = "green4", # Color of the estimates and CIs
     main = "Change Estimates' Color: Hainmueller and Hangartner (2019)")

plot(out, 
     preset = "grayscale", # Use grayscale colors
     main = "Grayscale Preset Colors")

plot(out, 
     color = "green4", # color of the estimates and CIs
     connected = TRUE, # Connect the points with lines
     est.lwidth = 1.2, # Makes the lines thicker
     est.pointsize = 3 # Makes the points larger
)

plot(out, 
     connected = TRUE, 
     ci.outline = TRUE,
     main = "The Effect of Coethnic Mobilization") # Outline the confidence interval band
plot(out.hh, 
     preset = "vibrant",
     ci.outline = TRUE,
     main = "The Effect of Indirect Democracy") # Outline the confidence interval band

plot(out, 
     est.lwidth = 1.5, # Makes the confidence intervals thicker
     est.pointsize = 3, # Makes the points larger
     lcolor = c("red","skyblue"),      # Color for horizontal and vertical lines
     lwidth = 2,                # Widths of the horizontal and vertical lines
     main = "Line Customization")

plot(out, 
     count.color = "lightblue", # Color of the histogram bars
     count.outline.color = "darkblue", # Outline color of the histogram bars
     count.alpha = 0.2, # Opacity of the histogram bars
     main = "Count Histogram Customization")

plot(out, type = "counterfactual",
     main = "Grumbach & Sahn (2020): Treated vs. Counterfactuals",
     ylab = "Proportion of Asian Donation",
     legend.pos = "bottom")

plot(out.hh, type = "counterfactual",
     main = "Hainmueller & Hangartner (2019): Treated vs. Counterfactuals",
     ylab = "Naturalization Rate",
     legend.pos = "bottom")

plot(out.hh, type = "counterfactual",
     main = "Hainmueller & Hangartner (2019): Treated vs. Counterfactuals",
     ylab = "Naturalization Rate",
     legend.pos = "bottom",
     ci.outline = TRUE, # Outline the confidence interval band
     color = "red3", # Color for the main line
     counterfactual.color = "green4") # Color for the counterfactual line

plot(out, type = "counterfactual", raw = "all")

plot(out, type = "counterfactual", raw = "band")

plot(out, type = "counterfactual",
     count.color = "black", # Color for the count histogram
     count.alpha = 1, # Opacity for the count histogram
     color = "red", # Color for the main line
     counterfactual.color = "purple", # Color for the counterfactual line
     counterfactual.raw.treated.color = "orange", # Color for the treated units
     counterfactual.linetype = "dotted", # Line type for the counterfactual line
     raw = "all", 
     main = "Counterfactual Plot with Custom Colors")

plot(effect(out.hh), main = "Cumulative Effect of Indirect Democracy",
     ylab = "Cumulative Effect on Naturalization Rate")

# flag units that ever have a 1 to 0 change in d
rev_flag <- tapply(gs2020[["cand_A_all"]],
                   gs2020[["district_final"]],
                   function(x) any(diff(x) < 0))

# units with no reversals
good_units <- names(rev_flag)[!rev_flag]

# subset the desired rows
gs2020_no_reversals <- gs2020[gs2020[["district_final"]] %in% good_units, ]


out_no_reversals <- fect(Y = "general_sharetotal_A_all", 
                         D = "cand_A_all" , 
                         X = c("cand_H_all", "cand_B_all") ,
                         index = c("district_final", "cycle"), 
                         data = gs2020_no_reversals,
                         method = "fe", 
                         force =  "two-way", 
                         se = TRUE, parallel = TRUE,
                         nboots = 100, 
                         keep.sims = TRUE)

plot(effect(out_no_reversals), xlim = c(1, 2))

plot(out, type = "equiv", bound = "equiv", tost.threshold = 0.1, 
     ylim = c(-0.15, 0.15))

plot(out, type = "equiv", bound = "min", ylim = c(-0.15, 0.15))

plot(out, type = "equiv", tost.threshold = 0.1, ylim = c(-0.15, 0.15))

plot(out, type = "equiv",
     ylim = c(-0.25, 0.25),
     stats = c("F.p", "equiv.p"),
     stats.labs = c("F Test P-value", "Equivalence P-value"),
     stats.pos = c(-8, 0.2),   # (x, y) position for the stats text
     show.stats = TRUE,     # Can be switched off to hide all test stats
     main = "Statistical Test Annotations")

out_fe_placebo <- fect(Y = "general_sharetotal_A_all", D = "cand_A_all", X = c("cand_H_all", "cand_B_all"), data = gs2020, 
                       index = c("district_final", "cycle"), force = "two-way",
                       method = "fe", CV = FALSE, parallel = TRUE,
                       se = TRUE, nboots = 1000, placeboTest = TRUE,
                       placebo.period = c(-2, 0))

plot(out_fe_placebo)

plot(out_fe_placebo, placebo.color = "green4")

plot(out_fe_placebo, type = "exit")

out_fe_carryover <- fect(Y = "general_sharetotal_A_all", D = "cand_A_all", X = c("cand_H_all", "cand_B_all"), data = gs2020, 
                       index = c("district_final", "cycle"), force = "two-way",
                         parallel = TRUE, se = TRUE, CV = FALSE,
                         nboots = 1000, carryoverTest = TRUE,
                         carryover.period = c(1, 3))
plot(out_fe_carryover)

plot(out_fe_carryover, type = "status",
     status.treat.color      = "#D55E00",  # Color for treated units
     status.control.color    = "#0072B2",  # Color for control units
     status.carryover.color  = "#CC79A7",  # Color for carryover units
     status.missing.color    = "#009E73",  # Color for missing data
     status.background.color = "#F3EAD2",  # Background color
     main = "Status Plot")

plot(out, type = "box", xlim = c(-12, 3))

plot(out, type = "calendar", main = "The Effect of Coethnic Mobilization")
plot(out.hh, type = "calendar", xlim = c(1995, 2009),
     main = "The Effect of Indirect Democracy")

# From: ./04-gsynth.Rmd
set.seed(1234)
rm(list = ls())

library(fect)
data(fect)
ls()

head(simgsynth)

library(panelView)
panelview(Y ~ D, data = simgsynth,  index = c("id","time"), pre.post = TRUE) 

panelview(Y ~ D, data = simgsynth,  index = c("id","time"), type = "outcome") 

system.time(
out <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"), 
            method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 5), 
            se = TRUE, nboots = 1000, vartype = 'parametric', 
            parallel = FALSE))

print(out)
out$est.att
out$est.avg
out$est.beta

system.time(
out <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"), method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 5), se = TRUE, nboots = 1000,vartype = 'parametric', parallel = TRUE, cores = 16)
)

out2 <- fect(Y ~ D + X1 + X2, data = simgsynth,  index = c("id","time"), 
               method = "gsynth", force = "two-way", 
               CV = TRUE, r = c(0, 5), se = TRUE,
               vartype = "jackknife", 
               parallel = TRUE, cores = 16)


a <- plot(out) # by default, type = "gap"
print(a)

plot(out, theme.bw = FALSE) 

plot(out, connected = TRUE) 

plot(out, connected = TRUE, show.points = FALSE) 

plot(out, type = "gap", ylim = c(-3,12), xlab = "Period", 
     main = "Estimated ATT (FEct)")

plot(out, type = "status", yticklabels="0", 
     xticklabels=c("5", "10", "15","20", "25", "30") )

plot(out, type = "box", xlab = "time",
     xticklabels=c("-19", "-15", "-10", "-5","0","5","10") )

plot(out, type = "box", xlim = c(-15, 10), 
     xticklabels=c( "-15", "-10", "-5","0","5","10"))

plot(out,type = "calendar")

plot(out, type = "loadings")

plot(out, type = "factors", xlab = "Time")

plot(out, type = "counterfactual")

plot(out, type = "counterfactual", raw = "all")

plot(out, type = "counterfactual", raw = "band")

plot(out, type = "equiv", ylim = c(-5, 5))

plot(out, type = "equiv", show.stats =  FALSE)

plot(out, type = "equiv", stats.pos = c(-19, 4.5), ylim = c(-5, 5))

out.ife <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"), 
          force = "two-way", method = "ife", CV = TRUE, r = c(0, 5), 
          se = TRUE, nboots = 200, parallel = TRUE)

plot(out.ife, main = "Estimated ATT (EM)")

out.mc <- fect(Y ~ D + X1 + X2, data = simgsynth, 
               index = c("id","time"), 
                force = "two-way", method = "mc", CV = TRUE, 
                se = TRUE, nboots = 200, parallel = TRUE)


plot(out.mc, main = "Estimated ATT (MC)")

panelview(turnout ~ policy_edr, data = turnout, 
          index = c("abb","year"), pre.post = TRUE, 
          by.timing = TRUE) 

panelview(turnout ~ policy_edr, data = turnout, 
          index = c("abb","year"), type = "outcome", 
          main = "EDR Reform and Turnout", 
          by.group = TRUE)

out0 <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor, 
               data = turnout, index = c("abb","year"), 
               se = TRUE, method = "gsynth",
               r = 0, CV = FALSE, force = "two-way", 
               nboots = 1000, seed = 02139)

plot(out0, type = "gap", xlim = c(-15, 15))

out_turnout <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor, 
              data = turnout,  index = c("abb","year"), 
              se = TRUE, method = "gsynth", vartype = "parametric",
              r = c(0, 5), CV = TRUE, force = "two-way", 
              nboots = 1000, seed = 02139, keep.sims = TRUE)

dim(out_turnout$wgt.implied)
sort(out_turnout$wgt.implied[,8])

plot(out_turnout, xlim = c(-10, 5), ylim=c(-15,15))

plot(out_turnout, type = "status",xlab = "Year", ylab = "State", main = "Treatment Status", 
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960, 
                   1968, 1976, 1984, 1992, 2000, 2008), xangle=10)

plot(out_turnout, type = "counterfactual")

plot(out_turnout, type = "counterfactual", id = "WI", main = "Wisconsin")

plot(out_turnout, type = "box", 
     xticklabels=c("-20", "-15", "-10", "-5","0","5","10"))

plot(out_turnout, type = "calendar", ylim = c(-15,15))

plot(out_turnout, type = "factors", xlab = "Year")

plot(out_turnout, type = "loadings")

set.seed(123456)
turnout.ub <- turnout[-c(which(turnout$abb=="WY")[1:15], 
                         sample(1:nrow(turnout),50,replace=FALSE)),]

panelview(turnout ~ policy_edr + policy_mail_in + policy_motor, 
          data = turnout.ub,  index = c("abb","year"), 
          pre.post = TRUE) 

out_ub <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor, 
              data = turnout.ub,  index = c("abb","year"), 
              se = TRUE, method = "gsynth", 
              r = c(0, 5), CV = TRUE, force = "two-way", 
              parallel = TRUE, min.T0 = 8, 
              nboots = 1000, seed = 02139)

plot(out_ub, type = "status",
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960, 
                   1968, 1976, 1984, 1992, 2000, 2008),
     xangle=10)

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

# From: ./05-panel.Rmd
# install packages from CRAN
packages <- c("dplyr", "fixest", "did", "didimputation", 
              "panelView", "ggplot2", "bacondecomp", "HonestDiD",
              "DIDmultiplegtDYN", "PanelMatch")
install.packages(setdiff(packages, rownames(installed.packages())))  

# install most up-to-date "fect" from Github
if ("fect" %in% rownames(installed.packages()) == FALSE) {
  devtools:: install_github("xuyiqing/fect")
}

# install forked "HonestDiD" package compatible with "fect"
if ("HonestDiDFEct" %in% rownames(installed.packages()) == FALSE) {
  devtools:: install_github("lzy318/HonestDiDFEct")
}

library(dplyr)
library(readstata13)
library(fixest)
library(did)
library(fect)
library(panelView)
library(PanelMatch)
library(ggplot2)
library(bacondecomp)
library(fect)
library(didimputation)
library(doParallel)
library(HonestDiD)
library(HonestDiDFEct)
library(DIDmultiplegtDYN) # may require XQuartz 

data(fect)
data <- hh2019
head(data)

panelview(nat_rate_ord ~ indirect, data = data, index = c("bfs","year"), 
  xlab = "Year", ylab = "Unit", display.all = T,
  gridOff = TRUE, by.timing = TRUE)

panelview(data = data,Y='nat_rate_ord',
          D='indirect',index=c("bfs","year"),
          by.timing = TRUE, display.all = TRUE,
          type = "outcome", by.cohort = TRUE)

# remember to cluster standard errors
model.twfe.0 <- feols(nat_rate_ord~indirect|bfs+year,
                      data=data, cluster = "bfs") 
print(model.twfe.0)

data.complete <- data[which(!is.na(data$nat_rate_ord)),] # bacon requires no missingness in the data
df_bacon <- bacon(nat_rate_ord~indirect,
                  data = data.complete,
                  id_var = "bfs",
                  time_var = "year")
ggplot(df_bacon) +
   aes(x = weight, y = estimate, shape = factor(type), color = factor(type)) +
   labs(x = "Weight", y = "Estimate", shape = "Type", color = 'Type') +
   geom_point()

print(aggregate(df_bacon$estimate * df_bacon$weight, 
                list(df_bacon$type), FUN=sum))

# drop always treated units
df <- as.data.frame(data %>% 
                      group_by(bfs) %>%
                      mutate(treatment_mean = mean(indirect,na.rm = TRUE)))
df.use <- df[which(df$treatment_mean<1),]

# Re-estimate TWFE on this Sub-sample
model.twfe.1 <- feols(nat_rate_ord~indirect|bfs+year,
                      data=df.use, cluster = "bfs")
print(model.twfe.1)

df.use <- get.cohort(df.use, D = "indirect", index=c("bfs","year"), 
                     start0 = TRUE)
head(df.use[,-5],19)

# Dynamic TWFE
df.twfe <- df.use
# drop always treated units
df.twfe$treat <- as.numeric(df.twfe$treatment_mean>0) 
df.twfe[which(is.na(df.twfe$Time_to_Treatment)),'Time_to_Treatment'] <- 0 # can be an arbitrary value
twfe.est <- feols(nat_rate_ord ~ i(Time_to_Treatment, treat, ref = -1)| bfs + year, 
                  data = df.twfe, cluster = "bfs")
twfe.output <- as.matrix(twfe.est$coeftable)
print(round(twfe.output, 3))

twfe.output <- as.data.frame(twfe.output)
twfe.output$Time <- c(c(-18:-2),c(0:17))+1 
p.twfe <- esplot(twfe.output,Period = 'Time',Estimate = 'Estimate',
                               SE = 'Std. Error', xlim = c(-12,10))
p.twfe

twfe.output <- as.data.frame(twfe.est$coeftable)
twfe.output$Time <- c(c(-18:-2),c(0:17)) 
p.twfe <- esplot(twfe.output, Period = 'Time',
                 Estimate = 'Estimate', SE = 'Std. Error', 
                 xlim = c(-12,10),start0 = TRUE)
p.twfe

df.st <- NULL
target.cohorts <- setdiff(unique(df.use$Cohort),"Control")
k <- 1
for(cohort in target.cohorts){
  df.sub <- df.use[which(df.use$Cohort%in%c(cohort,"Control")),]
  df.sub$stack <- k
  df.st <- rbind(df.st,df.sub)
  k <- k + 1
}
df.st$st_unit <- as.numeric(factor(paste0(df.st$stack,'-',df.st$bfs)))
df.st$st_year <- as.numeric(factor(paste0(df.st$stack,'-',df.st$year)))
model.st <- feols(nat_rate_ord~indirect|st_unit+st_year,
                  data=df.st, cluster = "st_unit")

print(model.st)

df.st$treat <- as.numeric(df.st$treatment_mean>0)
df.st[which(is.na(df.st$Time_to_Treatment)),'Time_to_Treatment'] <- 1000 
# note, this "1000" can be arbitrary value

st.est <- feols(nat_rate_ord ~ 
                  i(Time_to_Treatment, treat, ref = -1)| st_unit + 
                  st_year,data = df.st,cluster = "st_unit")

# make plot
st.output <- as.data.frame(st.est$coeftable)
st.output$Time <- c(c(-18:-2),c(0:17))+1 
p.st <- esplot(st.output,Period = 'Time',Estimate = 'Estimate',
                               SE = 'Std. Error', xlim = c(-12,10))
p.st

knitr::include_graphics("fig/fig_iw.png")

df.sa <- df.use
df.sa[which(is.na(df.sa$FirstTreat)),"FirstTreat"] <- 1000 
# above, replace NA with an arbitrary number 

model.sa.1 <- feols(nat_rate_ord~sunab(FirstTreat,year)|bfs+year,
                    data = df.sa, cluster = "bfs")
summary(model.sa.1,agg = "ATT")

sa.output <- as.data.frame(as.matrix(model.sa.1$coeftable))
sa.output$Time <- c(c(-18:-2),c(0:17)) + 1
p.sa <- esplot(sa.output,Period = 'Time',Estimate = 'Estimate',
                             SE = 'Std. Error', xlim = c(-12,10))
p.sa

df.cs <- df.use
df.cs[which(is.na(df.cs$FirstTreat)),"FirstTreat"] <- 0 # replace NA with 0
cs.est.1 <- att_gt(yname = "nat_rate_ord",
                 gname = "FirstTreat",
                 idname = "bfs",
                 tname = "year",
                 xformla = ~1,
                 control_group = "nevertreated",
                 allow_unbalanced_panel = TRUE,
                 data = df.cs,
                 est_method = "reg")
cs.est.att.1 <- aggte(cs.est.1, type = "simple", na.rm=T, bstrap = F)
print(cs.est.att.1)

cs.att.1 <- aggte(cs.est.1, type = "dynamic",
                  bstrap=FALSE, cband=FALSE, na.rm=T) 
print(cs.att.1)

cs.output <- cbind.data.frame(Estimate = cs.att.1$att.egt,
                              SE = cs.att.1$se.egt,
                              time = cs.att.1$egt + 1)
p.cs.1 <- esplot(cs.output,Period = 'time',Estimate = 'Estimate',
                               SE = 'SE', xlim = c(-12,10))
p.cs.1

cs.est.1.u <- att_gt(yname = "nat_rate_ord",
                 gname = "FirstTreat",
                 idname = "bfs",
                 tname = "year",
                 xformla = ~1,
                 control_group = "nevertreated",
                 allow_unbalanced_panel = TRUE,
                 data = df.cs,
                 est_method = "reg", 
                 base_period = "universal")
cs.att.1.u <- aggte(cs.est.1.u, type = "dynamic",
                    bstrap=FALSE, cband=FALSE, na.rm=T) 
cs.output.u <- cbind.data.frame(Estimate = cs.att.1.u$att.egt,
                                SE = cs.att.1.u$se.egt,
                                time = cs.att.1.u$egt + 1)
p.cs.1.u <- esplot(cs.output.u,Period = 'time',Estimate = 'Estimate',
                               SE = 'SE', xlim = c(-12,10))
p.cs.1.u

knitr::include_graphics("fig/fig_cs.png")

cs.est.2 <- att_gt(yname = "nat_rate_ord",
                   gname = "FirstTreat",
                   idname = "bfs",
                   tname = "year",
                   xformla = ~1,
                   control_group = "notyettreated",
                   allow_unbalanced_panel = TRUE,
                   data = df.cs,
                   est_method = "reg")
cs.est.att.2 <- aggte(cs.est.2, type = "simple",na.rm=T, bstrap = F)
print(cs.est.att.2)

cs.est.2.u <- att_gt(yname = "nat_rate_ord", gname = "FirstTreat",
                     idname = "bfs", tname = "year", xformla = ~1,
                     control_group = "notyettreated",
                     allow_unbalanced_panel = TRUE,
                     data = df.cs, est_method = "reg", 
                     base_period = "universal")

cs.att.2.u <- aggte(cs.est.2.u, type = "dynamic",
                    bstrap=FALSE, cband=FALSE, na.rm=T) 

# plot
cs.output.u <- cbind.data.frame(Estimate = cs.att.2.u$att.egt,
                                SE = cs.att.2.u$se.egt,
                                time = cs.att.2.u$egt + 1)
p.cs.2.u <- esplot(cs.output.u,Period = 'time',Estimate = 'Estimate',
                               SE = 'SE', xlim = c(-12,10))
p.cs.2.u

knitr::include_graphics("fig/fig_pm.png")

didm.results <- did_multiplegt_dyn(
      df = df.use,
      outcome = "nat_rate_ord",
      group = "bfs",
      controls = NULL,
      time = "year",
      treatment = "indirect",
      effects = 12, 
      placebo = 9,
      cluster = "bfs",
      graph_off = TRUE
    )
print(didm.results)

T.post <- dim(didm.results$results$Effects)[1]
T.pre <- dim(didm.results$results$Placebos)[1]
didm.vis <- rbind(didm.results$results$Placebos,didm.results$results$Effects)
didm.vis <- as.data.frame(didm.vis)
didm.vis[,'Time'] <- c(c(-1:-(T.pre)),c(1:T.post))
est.dynamic <- didm.vis[,c(9,1,2,3,4)]
colnames(est.dynamic) <- c("T","estimate","se","lb","ub")
p.didm <- esplot(est.dynamic,Period = 'T',Estimate = 'estimate',
                               SE = 'se', xlim = c(-9, 9))
p.didm

df.pm <- df.use
# we need to convert the unit and time indicator to integer
df.pm[,"bfs"] <- as.integer(as.factor(df.pm[,"bfs"]))
df.pm[,"year"] <- as.integer(as.factor(df.pm[,"year"]))
df.pm <- df.pm[,c("bfs","year","nat_rate_ord","indirect")]

# Pre-processes and balances panel data
df.pm <- PanelData(panel.data = df.pm,
                    unit.id = "bfs",
                    time.id = "year",
                    treatment = "indirect",
                    outcome = "nat_rate_ord")

PM.results <- PanelMatch(lag=3, 
                         refinement.method = "none", 
                         panel.data = df.pm, 
                         qoi = "att", 
                         lead = c(0:3), 
                         match.missing = TRUE)

## For pre-treatment dynamic effects
PM.results.placebo <- PanelMatch(lag=3, 
                         refinement.method = "none", 
                         panel.data = df.pm, 
                         qoi = "att", 
                         lead = c(0:3), 
                         match.missing = TRUE,
                         placebo.test = TRUE)


# ATT
PE.results.pool <- PanelEstimate(PM.results, panel.data = df.pm, pooled = TRUE)
summary(PE.results.pool)

# Dynamic Treatment Effects
PE.results <- PanelEstimate(PM.results, panel.data = df.pm)
PE.results.placebo <- placebo_test(PM.results.placebo, panel.data = df.pm, plot = F)

# obtain lead and lag (placebo) estimates
est_lead <- as.vector(PE.results$estimate)
est_lag <- as.vector(PE.results.placebo$estimates)
sd_lead <- apply(PE.results$bootstrapped.estimates,2,sd)
sd_lag <- apply(PE.results.placebo$bootstrapped.estimates,2,sd)
coef <- c(est_lag, 0, est_lead)
sd <- c(sd_lag, 0, sd_lead)
pm.output <- cbind.data.frame(ATT=coef, se=sd, t=c(-2:4))

# plot
p.pm <- esplot(data = pm.output,Period = 't',
               Estimate = 'ATT',SE = 'se')
p.pm

knitr::include_graphics("fig/fig_fect.png")

out.fect <- fect(nat_rate_ord~indirect, data = df, 
                 index = c("bfs","year"),
                 method = 'fe', se = TRUE)
print(out.fect$est.avg)

fect.output <- as.matrix(out.fect$est.att)
head(fect.output)

df.impute <- df.use
df.impute[which(is.na(df.impute$FirstTreat)),"FirstTreat"] <- 0 
# above, replace NA with 0

out.impute <- did_imputation(data = df.impute,
                               yname = "nat_rate_ord",
                               gname = "FirstTreat",
                               tname = "year",
                               idname = "bfs",
                               cluster_var = "bfs")
out.impute

fect.output <- as.data.frame(fect.output)
fect.output$Time <- c(-17:18)
p.fect <- esplot(fect.output,Period = 'Time',Estimate = 'ATT',
                   SE = 'S.E.',CI.lower = "CI.lower", 
                   CI.upper = 'CI.upper',xlim = c(-12,10))
p.fect

model.impute <- did_imputation(data = df.impute,
                               yname = "nat_rate_ord",
                               gname = "FirstTreat",
                               tname = "year",
                               idname = "bfs",
                               cluster_var = "bfs",
                               pretrends = c(-13:-1),
                               horizon = TRUE)
model.impute$term <- as.numeric(model.impute$term)+1 
# above, set 1 as the first post-treatment period

# plot
to_plot <- as.data.frame(model.impute)
esplot(data=to_plot,Period = "term", 
       Estimate = 'estimate', SE = 'std.error',
       xlim = c(-12,10))
out.impute

out.fect.balance <- fect(nat_rate_ord~indirect, data = df, 
                         index = c("bfs","year"),
                         method = 'fe', se = TRUE, 
                         balance.period = c(-2,4))
# att
print(out.fect.balance$est.balance.avg)

# event study plot
fect.balance.output <- as.data.frame(out.fect.balance$est.balance.att)
fect.balance.output$Time <- c(-2:4)
p.fect.balance <- esplot(fect.balance.output,Period = 'Time',
                         Estimate = 'ATT', SE = 'S.E.',
                         CI.lower = "CI.lower", 
                         CI.upper = 'CI.upper')
p.fect.balance

data(fect)
data <- gs2020
data$cycle <- as.integer(as.numeric(data$cycle/2))
head(data)

y <- "general_sharetotal_A_all"
d <- "cand_A_all"
unit <- "district_final"
time <- "cycle"
controls <- c("cand_H_all", "cand_B_all")
index <- c("district_final", "cycle")

panelview(Y=y, D=d, X=controls, index = index, data = data, 
          xlab = "Time Period", ylab = "Unit", gridOff = TRUE, 
          by.timing = TRUE, cex.legend=5, cex.axis= 5, 
          cex.main = 10, cex.lab = 5)


model.twfe <- feols(general_sharetotal_A_all ~ cand_A_all + 
                      cand_H_all + cand_B_all | district_final + cycle,
                    data=data, cluster = "district_final") 
summary(model.twfe)

data_cohort <- get.cohort(data, index = index, D=d,start0 = TRUE)
# Generate a dummy variable treat
data_cohort$treat <- 0
data_cohort[which(data_cohort$Cohort!='Control'),'treat'] <- 1
data_cohort[which(is.na(data_cohort$Time_to_Treatment)), "treat"] <- 0

# remove observations that starts with treated status
remove <- intersect(which(is.na(data_cohort$Time_to_Treatment)),
                    which(data_cohort[,d]==1)) 
if(length(remove)>0){data_cohort <- data_cohort[-remove,]}

# replace missingness in Time_to_Treatment with an arbitrary number
data_cohort[which(is.na(data_cohort$Time_to_Treatment)), "Time_to_Treatment"] <- 999 

twfe.est <- feols(general_sharetotal_A_all ~ 
                    i(Time_to_Treatment, treat, ref = -1) + 
                    cand_H_all +cand_B_all | district_final + cycle,  
                  data = data_cohort, cluster = "district_final")

twfe.output <- as.data.frame(twfe.est$coeftable[c(1:25),])
twfe.output$Time <- c(c(-16:-2),c(0:9)) + 1 

# plot
p.twfe <- esplot(twfe.output,Period = 'Time',Estimate = 'Estimate',
                               SE = 'Std. Error', xlim = c(-15,1))
p.twfe

df.pm <- data_cohort
# we need to convert the unit and time indicator to integer
df.pm[,"district_final"] <- as.integer(as.factor(df.pm[,"district_final"]))
df.pm[,"cycle"] <- as.integer(as.factor(df.pm[,"cycle"]))
df.pm <- df.pm[,c("district_final","cycle","cand_A_all", 
                  "general_sharetotal_A_all")]

# Pre-processes and balances panel data
df.pm <- PanelData(panel.data = df.pm,
                    unit.id = "district_final",
                    time.id = "cycle",
                    treatment = "cand_A_all",
                    outcome = "general_sharetotal_A_all")

PM.results <- PanelMatch(lag=4, 
                         refinement.method = "none", 
                         panel.data = df.pm, 
                         qoi = "att", 
                         lead = 0, 
                         match.missing = TRUE)

## For pre-treatment dynamic effects
PM.results.placebo <- PanelMatch(lag=4, 
                         refinement.method = "none", 
                         panel.data = df.pm, 
                         qoi = "att", 
                         lead = 0, 
                         match.missing = TRUE,
                         placebo.test = TRUE)


PE.results.pool <- PanelEstimate(PM.results, panel.data = df.pm, pooled = TRUE)
summary(PE.results.pool)

# Dynamic Treatment Effects
PE.results <- PanelEstimate(PM.results, panel.data = df.pm)
PE.results.placebo <- placebo_test(PM.results.placebo, panel.data = df.pm,
                                   plot = FALSE)

est_lead <- as.vector(PE.results$estimate)
est_lag <- as.vector(PE.results.placebo$estimates)
sd_lead <- apply(PE.results$bootstrapped.estimates,2,sd)
sd_lag <- apply(PE.results.placebo$bootstrapped.estimates,2,sd)
coef <- c(est_lag, 0, est_lead)
sd <- c(sd_lag, 0, sd_lead)
pm.output <- cbind.data.frame(ATT=coef, se=sd, t=c(-3:1))

# plot
p.pm <- esplot(data = pm.output,Period = 't',
               Estimate = 'ATT',SE = 'se')
p.pm

model.fect <- fect(Y = "general_sharetotal_A_all", D = "cand_A_all", 
                   X= c("cand_H_all", "cand_B_all"), data = data, 
                   method = "fe", index = index, se = TRUE, 
                   parallel = TRUE, seed = 1234, force = "two-way")

print(model.fect$est.avg)

fect.output <- as.data.frame(model.fect$est.att)
fect.output$Time <- c(-15:10)
p.fect <- esplot(fect.output,Period = 'Time',Estimate = 'ATT',
                   SE = 'S.E.',CI.lower = "CI.lower", 
                   CI.upper = 'CI.upper', xlim = c(-15,1))
p.fect

plot(model.fect)

plot(model.fect, type = 'exit')

out.fect.p <- fect(Y = y, X = controls, D = d, data = data, index = index,
                   method = 'fe', se = TRUE, placeboTest = TRUE,
                   placebo.period = c(-2,0))

plot(out.fect.p, proportion = 0.1, stats = "placebo.p")

out.fect.c <- fect(Y = y, X = controls, D = d, data = data, index = index,
                   method = 'fe', se = TRUE, carryoverTest = TRUE, carryover.period = c(1,2))

# plot
plot(out.fect.c,  stats = "carryover.p", ylim = c(-0.15, 0.20))

out.fect.balance <- fect(Y = y, X = controls, D = d, data = data, 
                         index = index, method = 'fe', se = TRUE,
                         balance.period = c(-3,1))

# att
print(out.fect.balance$est.balance.avg)

# event study plot
fect.balance.output <- as.data.frame(out.fect.balance$est.balance.att)
fect.balance.output$Time <- c(-3:1)
p.fect.balance <- esplot(fect.balance.output,Period = 'Time',Estimate = 'ATT',
                   SE = 'S.E.',CI.lower = "CI.lower", 
                   CI.upper = 'CI.upper')
p.fect.balance

res_st <- did_wrapper(
  data   = hh2019,
  Y      = "nat_rate_ord",
  D      = "indirect",
  index  = c("bfs", "year"),
  method = "st",
  se     = "default"
)
print(res_st)

res_st <- did_wrapper(
  data   = hh2019,
  Y      = "nat_rate_ord",
  D      = "indirect",
  index  = c("bfs", "year"),
  method = "st",
  se     = "boot"
)
print(res_st)

esplot(data = res_st, main = "Stacked DID", xlim = c(-12,10))

# From: ./06-sens.Rmd
# install packages from CRAN
packages <- c("dplyr", "panelView", "ggplot2") # Removed HonestDiD, doParallel
install.packages(setdiff(packages, rownames(installed.packages())))  

# install most up-to-date "fect" from Github
if ("fect" %in% rownames(installed.packages()) == FALSE) {
  devtools:: install_github("xuyiqing/fect")
}

# install forked "HonestDiD" package compatible with "fect"
if ("HonestDiDFEct" %in% rownames(installed.packages()) == FALSE) {
  devtools:: install_github("lzy318/HonestDiDFEct") # This is used by fect_sens
}

library(dplyr)
library(fect)
library(panelView)
library(ggplot2)
library(HonestDiDFEct) # Required for fect_sens to work

data(fect)
data <- hh2019
head(data)

out.fect.placebo <- fect(nat_rate_ord~indirect, data = hh2019, 
                         index = c("bfs","year"),
                         method = 'fe', se = TRUE, 
                         placeboTest = TRUE, placebo.period = c(-2,0))

# Define post-treatment periods and sensitivity parameters for fect_sens
T.post <- 10 # Number of post-treatment periods based on original analysis
post_periods_vec <- 1:T.post

# Parameters for Relative Magnitude (RM) restriction
Mbar_vec_avg_rm <- seq(0, 1, by = 0.1)    # For average ATT plot
Mbar_vec_period_rm <- c(0, 0.5)          # For period-by-period ATT plot

# Parameters for Smoothness restriction
M_vec_avg_smooth <- seq(0, 0.25, by = 0.05) # For average ATT plot
M_vec_period_smooth <- c(0, 0.1)           # For period-by-period ATT plot

# Run sensitivity analysis using fect_sens
# This function augments out.fect.placebo with sensitivity results
out.fect.placebo <- fect_sens(
  fect.out      = out.fect.placebo,
  post.periods  = post_periods_vec,
  Mbarvec       = Mbar_vec_avg_rm,
  periodMbarvec = Mbar_vec_period_rm,
  Mvec          = M_vec_avg_smooth,
  periodMvec    = M_vec_period_smooth,
  parallel      = FALSE # Set to TRUE for parallel processing if desired
)

plot(out.fect.placebo,
     type = "sens",
     restrict = "rm",
     main = "Relative Magnitude Restriction")

plot(out.fect.placebo,
    type = "sens_es",
    restrict = "rm",
    main = "ATTs with Robust Confidence Sets (RM)",
    ylab = "Coefficients and 95% CI",
    xlim = c(-12,10), 
    ylim = c(-6,8), 
    show.count = TRUE)

plot(out.fect.placebo,
    type = "sens_es",
    restrict = "rm",
    main = "ATTs with Robust Confidence Sets (RM)",
    ylab = "Coefficients and 95% CI",
    xlim = c(-12,10), 
    ylim = c(-6,8), 
    show.count = TRUE,
    sens.colors = c("blue", "red"))

plot(out.fect.placebo,
    type = "sens",
    restrict = "sm",
    main = "Smoothness Restriction")

plot(out.fect.placebo,
    type = "sens_es",
    restrict = "sm",
    main = "ATTs with Robust Confidence Sets (Smoothness)",
    ylab = "Coefficients and 95% CI",
    xlim = c(-12,10), # Adjusted to match original detailed plot
    ylim = c(-12,15),
    show.count = TRUE)

# From: ./07-cheetsheet.Rmd
