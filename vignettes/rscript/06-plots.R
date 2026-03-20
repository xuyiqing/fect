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
# From: ./06-plots.Rmd
##########################

# load libraries and data
library(ggplot2)
library(panelView)
library(fect)
data(fect)
ls()

##--- Gap Plot ---##

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

plot(out) # the effect of co-ethnic mobilization
plot(out.hh) # the effect of indirect democracy on naturalization rate

# Starting period
plot(out, start0 = TRUE,
     main = "Custom Starting Period")

# Connected estimates
plot(out,
     post.color = "green4",
     connected = TRUE,
     est.lwidth = 1.2,
     est.pointsize = 3)

plot(out,
     connected = TRUE,
     ci.outline = TRUE,
     main = "The Effect of Coethnic Mobilization")
plot(out.hh,
     preset = "vibrant",
     ci.outline = TRUE,
     main = "The Effect of Indirect Democracy")

# Presets
plot(out,
     preset = "vibrant",
     main = "Vibrant Preset Colors: Grumbach and Sahn (2020)")
plot(out.hh,
     preset = "vibrant",
     main = "Vibrant Preset Colors: Hainmueller and Hangartner (2019)")

plot(out,
     preset = "grayscale",
     main = "Grayscale Preset Colors")

# Colors
plot(out.hh,
     preset = "vibrant",
     post.color = "green4",
     main = "Change Estimates' Color: Hainmueller and Hangartner (2019)")

# Confidence intervals
plot(out, plot.ci = "0.9",
     main = "90% confidence intervals")

# Count bars
plot(out,
     count.color = "lightblue",
     count.outline.color = "darkblue",
     count.alpha = 0.2,
     main = "Count Histogram Customization")

# Axis customization
plot(out,
     xlim = c(-10, 1),
     ylim = c(-0.15, 0.30),
     xlab = "Custom Time Axis",
     ylab = "Estimated ATT",
     xangle = 90,
     xbreaks = seq(-10, 1, by = 2),
     gridOff = TRUE,
     main = "Axis and Legend Customization")

# Text sizes
plot(out,
     ylim = c(-0.15, 0.3),
     theme.bw = FALSE,
     cex.main = 1.25,
     cex.axis = 1.2,
     cex.lab = 1.2,
     cex.legend = 1,
     cex.text = 1.2,
     main = "Text and Theme Customization")

# Reference lines
plot(out,
     est.lwidth = 1.5,
     est.pointsize = 3,
     lcolor = c("red","skyblue"),
     lwidth = 2,
     main = "Line Customization")

##--- Pretrend Tests ---##

# Placebo test
out_fe_placebo <- fect(Y = "general_sharetotal_A_all", D = "cand_A_all", X = c("cand_H_all", "cand_B_all"), data = gs2020,
                       index = c("district_final", "cycle"), force = "two-way",
                       method = "fe", CV = FALSE, parallel = TRUE,
                       se = TRUE, nboots = 1000, placeboTest = TRUE,
                       placebo.period = c(-2, 0))

plot(out_fe_placebo)

plot(out_fe_placebo, placebo.color = "green4")

# Equivalence test
plot(out, type = "equiv", bound = "equiv", tost.threshold = 0.1,
     ylim = c(-0.15, 0.15))

plot(out, type = "equiv", bound = "min", ylim = c(-0.15, 0.15))

plot(out, type = "equiv", tost.threshold = 0.1, ylim = c(-0.15, 0.15))

plot(out, type = "equiv",
     ylim = c(-0.25, 0.25),
     stats = c("F.p", "equiv.p"),
     stats.labs = c("F Test P-value", "Equivalence P-value"),
     stats.pos = c(-8, 0.2),
     show.stats = TRUE,
     main = "Statistical Test Annotations")

##--- Carryover Test ---##

plot(out_fe_placebo, type = "exit")

out_fe_carryover <- fect(Y = "general_sharetotal_A_all", D = "cand_A_all", X = c("cand_H_all", "cand_B_all"), data = gs2020,
                       index = c("district_final", "cycle"), force = "two-way",
                         parallel = TRUE, se = TRUE, CV = FALSE,
                         nboots = 1000, carryoverTest = TRUE,
                         carryover.period = c(1, 3))
plot(out_fe_carryover)

##--- Cumulative Effects ---##

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

##--- Effect Heterogeneity ---##

# Box plot
plot(out, type = "box", xlim = c(-12, 3))

# Calendar plot
plot(out, type = "calendar", main = "The Effect of Coethnic Mobilization")
plot(out.hh, type = "calendar", xlim = c(1995, 2009),
     main = "The Effect of Indirect Democracy")

##--- Other Plot Types ---##

# Status plot
plot(out_fe_carryover, type = "status",
     status.treat.color      = "#D55E00",
     status.control.color    = "#0072B2",
     status.carryover.color  = "#CC79A7",
     status.missing.color    = "#009E73",
     status.background.color = "#F3EAD2",
     main = "Status Plot")

# Counterfactual plot
plot(out, type = "counterfactual",
     main = "Grumbach & Sahn (2020): Treated vs. Counterfactuals",
     ylab = "Proportion of Asian Donation",
     legend.pos = "bottom")

plot(out.hh, type = "counterfactual",
     main = "Hainmueller & Hangartner (2019): Treated vs. Counterfactuals",
     ylab = "Naturalization Rate",
     legend.pos = "top")

plot(out.hh, type = "counterfactual",
     main = "Hainmueller & Hangartner (2019): Treated vs. Counterfactuals",
     ylab = "Naturalization Rate",
     legend.pos = "bottom",
     ci.outline = TRUE,
     color = "red3",
     counterfactual.color = "green4")

plot(out, type = "counterfactual", raw = "all")

plot(out, type = "counterfactual", raw = "band")

plot(out, type = "counterfactual",
     count.color = "black",
     count.alpha = 1,
     color = "red",
     counterfactual.color = "purple",
     counterfactual.raw.treated.color = "orange",
     counterfactual.linetype = "dotted",
     raw = "all",
     main = "Counterfactual Plot with Custom Colors")

##--- Standalone esplot() ---##

# Basic usage with data frames
es_data <- data.frame(
  Time = as.numeric(rownames(out$est.att)),
  ATT = out$est.att[, "ATT"],
  CI.lower = out$est.att[, "CI.lower"],
  CI.upper = out$est.att[, "CI.upper"]
)

esplot(es_data, Period = "Time",
       main = "Event Study Plot with esplot()",
       ylab = "Estimated ATT",
       xlab = "Periods Since Treatment")

# Using fect objects directly
esplot(out, main = "Direct from fect object")

# Connected line style
esplot(es_data, Period = "Time",
       connected = TRUE,
       main = "Connected Event Study Plot",
       ylab = "Estimated ATT")

# Highlighting periods
esplot(es_data, Period = "Time",
       highlight.periods = c(-2, -1, 0),
       highlight.colors = c("orange", "orange", "red"),
       main = "Highlighting Key Periods",
       ylab = "Estimated ATT")

# ---- Factors and loadings ----
out_ife <- fect(nat_rate_ord ~ indirect,
                data = hh2019,
                index = c("bfs", "year"),
                method = "ife", r = 2,
                se = TRUE, parallel = TRUE, nboots = 200)

plot(out_ife, type = "factors", main = "Estimated Latent Factors")
plot(out_ife, type = "factors", include.FE = FALSE,
     main = "Factors without Fixed Effects")
plot(out_ife, type = "loadings", main = "Factor Loadings")
