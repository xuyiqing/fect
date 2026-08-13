## ----.common, include = FALSE-------------------------------------------------
source("_common.R")


## ----load, message=FALSE------------------------------------------------------
# load libraries and data
library(ggplot2)
library(panelView)
data(gs2020)
data(hh2019)
ls()


## ----est, cache = TRUE--------------------------------------------------------
out <- fect(Y = "general_sharetotal_A_all",
            D = "cand_A_all",
            X = c("cand_H_all", "cand_B_all"),
            index = c("district_final", "cycle"),
            data = gs2020, method = "fe",
            force = "two-way", se = TRUE,
            parallel = TRUE, cores = 16, nboots = 1000)

out.hh <- fect(nat_rate_ord ~ indirect,
               data = hh2019,
               index = c("bfs","year"),
               method = 'fe', se = TRUE,
               parallel = TRUE, cores = 16, nboots = 1000,
               keep.sims = TRUE)


## ----plot-gap-default---------------------------------------------------------
plot(out) # the effect of co-ethnic mobilization
plot(out.hh) # the effect of indirect democracy on naturalization rate


## ----begin-post-customization-------------------------------------------------
plot(out, start0 = TRUE,
     main = "Custom Starting Period")


## ----connected-estimates------------------------------------------------------
plot(out,
     post.color = "green4",
     connected = TRUE,
     est.lwidth = 1.2,
     est.pointsize = 3)


## ----ci-outline---------------------------------------------------------------
plot(out,
     connected = TRUE,
     ci.outline = TRUE,
     main = "The Effect of Coethnic Mobilization")
plot(out.hh,
     preset = "vibrant",
     ci.outline = TRUE,
     main = "The Effect of Indirect Democracy")


## ----preset-vibrant-----------------------------------------------------------
plot(out,
     preset = "vibrant",
     main = "Vibrant Preset Colors: Grumbach and Sahn (2020)")
plot(out.hh,
     preset = "vibrant",
     main = "Vibrant Preset Colors: Hainmueller and Hangartner (2019)")


## ----preset-grayscale---------------------------------------------------------
plot(out,
     preset = "grayscale",
     main = "Grayscale Preset Colors")


## ----preset-vibrant2----------------------------------------------------------
plot(out.hh,
     preset = "vibrant",
     post.color = "green4",
     main = "Change Estimates' Color: Hainmueller and Hangartner (2019)")


## ----ci-raw-customization-----------------------------------------------------
plot(out, plot.ci = "0.9",
     main = "90% confidence intervals")


## ----count-histogram-customization--------------------------------------------
plot(out,
     count.color = "lightblue",
     count.outline.color = "darkblue",
     count.alpha = 0.2,
     main = "Count Histogram Customization")


## ----axis-legend-customization------------------------------------------------
plot(out,
     xlim = c(-10, 1),
     ylim = c(-0.15, 0.30),
     xlab = "Custom Time Axis",
     ylab = "Estimated ATT",
     xangle = 90,
     xbreaks = seq(-10, 1, by = 2),
     gridOff = TRUE,
     main = "Axis and Legend Customization")


## ----text-customization-------------------------------------------------------
plot(out,
     ylim = c(-0.15, 0.3),
     theme.bw = FALSE,
     cex.main = 1.25,
     cex.axis = 1.2,
     cex.lab = 1.2,
     cex.legend = 1,
     cex.text = 1.2,
     main = "Text and Theme Customization")


## ----title-bold-centered------------------------------------------------------
plot(out, main = "Bold Centered Title") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## ----title-multi--------------------------------------------------------------
plot(out, main = "Centered, Bold, Larger") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))


## ----legacy-style-default-----------------------------------------------------
plot(out, legacy.style = TRUE)


## ----line-bound-customization-------------------------------------------------
plot(out,
     est.lwidth = 1.5,
     est.pointsize = 3,
     lcolor = c("red","skyblue"),
     lwidth = 2,
     main = "Line Customization")


## ----counterfactual-----------------------------------------------------------
plot(out, type = "counterfactual",
     main = "Grumbach & Sahn (2020): Treated vs. Counterfactuals",
     ylab = "Proportion of Asian Donation",
     legend.pos = "bottom")


## ----counterfactual-hh--------------------------------------------------------
plot(out.hh, type = "counterfactual",
     main = "Hainmueller & Hangartner (2019): Treated vs. Counterfactuals",
     ylab = "Naturalization Rate",
     legend.pos = "top")


## ----counterfactual-colors----------------------------------------------------
plot(out.hh, type = "counterfactual",
     main = "Hainmueller & Hangartner (2019): Treated vs. Counterfactuals",
     ylab = "Naturalization Rate",
     legend.pos = "bottom",
     ci.outline = TRUE,
     color = "red3",
     counterfactual.color = "green4")


## ----counterfactual-rawall----------------------------------------------------
plot(out, type = "counterfactual", raw = "all")


## ----counterfactual-rawband---------------------------------------------------
plot(out, type = "counterfactual", raw = "band")


## ----counterfactual-colors2---------------------------------------------------
plot(out, type = "counterfactual",
     count.color = "black",
     count.alpha = 1,
     color = "red",
     counterfactual.color = "purple",
     counterfactual.raw.treated.color = "orange",
     counterfactual.linetype = "dotted",
     raw = "all",
     main = "Counterfactual Plot with Custom Colors")


## ----placebo_fit, cache = TRUE------------------------------------------------
out_fe_placebo <- fect(Y = "general_sharetotal_A_all", D = "cand_A_all", X = c("cand_H_all", "cand_B_all"), data = gs2020,
                       index = c("district_final", "cycle"), force = "two-way",
                       method = "fe", CV = FALSE, parallel = TRUE, cores = 16,
                       se = TRUE, nboots = 1000, placeboTest = TRUE,
                       placebo.period = c(-2, 0))


## ----placebo------------------------------------------------------------------
plot(out_fe_placebo)


## ----plot-placebo-connected---------------------------------------------------
plot(out_fe_placebo, connected = TRUE, preset = "grayscale",
     main = "Placebo Test with Connected Estimates")


## ----plot-placebo-color-------------------------------------------------------
plot(out_fe_placebo, placebo.color = "green4")


## ----plot-placebo-fill--------------------------------------------------------
plot(out_fe_placebo, highlight.fill = TRUE,
     main = "Placebo test with background rectangle")


## ----plot-equiv-bound---------------------------------------------------------
plot(out, type = "equiv", bound = "equiv", tost.threshold = 0.1,
     ylim = c(-0.15, 0.15))


## ----plot-equiv-min-----------------------------------------------------------
plot(out, type = "equiv", bound = "min", ylim = c(-0.15, 0.15))


## ----plot-equiv-both----------------------------------------------------------
plot(out, type = "equiv", tost.threshold = 0.1, ylim = c(-0.15, 0.15))


## ----stats-customization------------------------------------------------------
plot(out, type = "equiv",
     ylim = c(-0.25, 0.25),
     stats = c("F.p", "equiv.p"),
     stats.labs = c("F Test P-value", "Equivalence P-value"),
     stats.pos = c(-8, 0.2),
     show.stats = TRUE,
     main = "Statistical Test Annotations")


## ----plot-exit-default--------------------------------------------------------
plot(out_fe_placebo, type = "exit")


## ----carryover_fit, cache = TRUE----------------------------------------------
out_fe_carryover <- fect(Y = "general_sharetotal_A_all", D = "cand_A_all", X = c("cand_H_all", "cand_B_all"), data = gs2020,
                       index = c("district_final", "cycle"), force = "two-way",
                         parallel = TRUE, cores = 16, se = TRUE, CV = FALSE,
                         nboots = 1000, carryoverTest = TRUE,
                         carryover.period = c(1, 3))


## ----carryover----------------------------------------------------------------
plot(out_fe_carryover)


## ----plot-cumulative-hh-------------------------------------------------------
cumu.hh <- estimand(out.hh, "att.cumu", "event.time")
esplot(cumu.hh, Period = "event.time",
       Estimate = "estimate", CI.lower = "ci.lo", CI.upper = "ci.hi",
       main = "Cumulative Effect of Indirect Democracy",
       ylab = "Cumulative Effect on Naturalization Rate")


## ----subset-no-reversals------------------------------------------------------
# flag units that ever have a 1 to 0 change in d
rev_flag <- tapply(gs2020[["cand_A_all"]],
                   gs2020[["district_final"]],
                   function(x) any(diff(x) < 0))

# units with no reversals
good_units <- names(rev_flag)[!rev_flag]

# subset the desired rows
gs2020_no_reversals <- gs2020[gs2020[["district_final"]] %in% good_units, ]



## ----no-reversals-est, cache = TRUE-------------------------------------------
out_no_reversals <- fect(Y = "general_sharetotal_A_all",
                         D = "cand_A_all" ,
                         X = c("cand_H_all", "cand_B_all") ,
                         index = c("district_final", "cycle"),
                         data = gs2020_no_reversals,
                         method = "fe",
                         force =  "two-way",
                         se = TRUE, parallel = TRUE, cores = 16,
                         nboots = 1000,
                         keep.sims = TRUE)


## ----cumulative-effects-------------------------------------------------------
cumu.gs <- estimand(out_no_reversals, "att.cumu", "event.time")
esplot(cumu.gs, Period = "event.time",
       Estimate = "estimate", CI.lower = "ci.lo", CI.upper = "ci.hi",
       xlim = c(1, 2))


## ----plot-box-hte-------------------------------------------------------------
plot(out, type = "box", xlim = c(-12, 3))


## ----plot-calendar-hte--------------------------------------------------------
plot(out, type = "calendar", main = "The Effect of Coethnic Mobilization")
plot(out.hh, type = "calendar", xlim = c(1995, 2009),
     main = "The Effect of Indirect Democracy")


## ----plot-hte-covariate-------------------------------------------------------
plot(out, type = "hte", covariate = "cand_B_all",
     main = "HTE by Black Candidate Presence",
     xlab = "Black Candidate Indicator",
     ylab = "Effect on Asian Donation Share")


## ----plot-hte-discrete-ch6----------------------------------------------------
plot(out, type = "hte", covariate = "cand_H_all",
     covariate.labels = c("No Hispanic Candidate", "Hispanic Candidate"),
     main = "HTE by Hispanic Candidate Presence",
     ylab = "Effect on Asian Donation Share")


## ----status-------------------------------------------------------------------
plot(out_fe_carryover, type = "status",
     status.treat.color      = "#D55E00",
     status.control.color    = "#0072B2",
     status.carryover.color  = "#CC79A7",
     status.missing.color    = "#009E73",
     status.background.color = "#F3EAD2",
     main = "Status Plot")


## ----est-ife, cache = TRUE----------------------------------------------------
out_ife <- fect(nat_rate_ord ~ indirect,
                data = hh2019,
                index = c("bfs", "year"),
                method = "ife", r = 2,
                se = TRUE, parallel = TRUE, cores = 16, nboots = 1000)


## ----plot-factors-------------------------------------------------------------
plot(out_ife, type = "factors", main = "Estimated Latent Factors")


## ----plot-factors-nofe--------------------------------------------------------
plot(out_ife, type = "factors", include.FE = FALSE,
     main = "Factors without Fixed Effects")


## ----plot-loadings------------------------------------------------------------
plot(out_ife, type = "loadings", main = "Factor Loadings")


## ----plot-loading-overlap, fig.width = 6, fig.height = 5----------------------
plot(out_ife, type = "loading.overlap")


## ----loading-overlap-r1-fit, cache = TRUE, message = FALSE, warning = FALSE----
out_ife_r1 <- fect(nat_rate_ord ~ indirect, data = hh2019,
                   index = c("bfs", "year"), method = "ife", r = 1,
                   se = TRUE, parallel = TRUE, cores = 16, nboots = 500)


## ----plot-loading-overlap-r1, fig.width = 6, fig.height = 5, message = FALSE, warning = FALSE----
plot(out_ife_r1, type = "loading.overlap")


## ----esplot-basic, fig.width = 6, fig.height = 4.5----------------------------
# Create example data from a fect result
es_data <- data.frame(Time = as.numeric(rownames(out$est.att)),
  ATT = out$est.att[, "ATT"],
  CI.lower = out$est.att[, "CI.lower"],
  CI.upper = out$est.att[, "CI.upper"]
)

esplot(es_data, Period = "Time",
       main = "Event Study Plot with esplot()",
       ylab = "Estimated ATT",
       xlab = "Periods Since Treatment",
       xlim = c(-15, 5), ylim = c(-0.3, 0.7))


## ----esplot-fect-object, fig.width = 6, fig.height = 4.5----------------------
esplot(out, main = "Direct from fect object")


## ----esplot-connected, fig.width = 6, fig.height = 4.5------------------------
esplot(es_data, Period = "Time",
       connected = TRUE,
       main = "Connected Event Study Plot",
       ylab = "Estimated ATT",
       xlim = c(-15, 5), ylim = c(-0.3, 0.7))


## ----esplot-highlight, fig.width = 6, fig.height = 4.5------------------------
esplot(es_data, Period = "Time",
       highlight.periods = c(-2, -1, 0),
       highlight.colors = c("orange", "orange", "red"),
       main = "Highlighting Key Periods",
       ylab = "Estimated ATT",
       xlim = c(-15, 5), ylim = c(-0.3, 0.7))

