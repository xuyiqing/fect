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
# From: ./03-plots.Rmd
##########################

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

