\name{plot.fect}
\alias{plot.fect}
\title{Plotting}
\description{Visualizes estimation results of the matrix completion method.}
\usage{\method{plot}{fect}(x, type = NULL, loo = "FALSE", 
            highlight = NULL, plot.ci = NULL, show.points = NULL,
            show.group = NULL, bound = NULL, vis = NULL,
            count = TRUE, proportion = 0.3, pre.periods = NULL,
            f.threshold = NULL, tost.threshold = NULL,
            effect.bound.ratio = FALSE, 
            stats = NULL, stats.labs = NULL, 
            main = NULL, 
            xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, 
            gridOff = FALSE, legendOff = FALSE, legend.pos = NULL, legend.nrow = NULL,
            legend.labs = NULL, stats.pos = NULL, theme.bw = TRUE,
            nfactors = NULL, include.FE = TRUE,
            id = NULL,
            cex.main = NULL, cex.main.sub = NULL, cex.axis = NULL,
            cex.lab = NULL, cex.legend = NULL, cex.text = NULL,
            axis.adjust = FALSE, axis.lab = "both", axis.lab.gap = c(0, 0),
            start0 = FALSE,
            return.test = FALSE,
            balance = NULL, weight = NULL, \dots)
}
\arguments{
<<<<<<< Updated upstream
  \item{x}{a \code{\link{fect}} object.}
  \item{type}{a string specifying the type of the plot. 
    \code{"gap"} plots the estimated period-wise ATT (dynamic treatment effects), 
    \code{"equiv"} shows the visualization of the equivalence test, 
    \code{"exit"} plots the estimated period-wise switch-off effects,
    \code{"status"} shows the treatment status of all observations, 
    \code{"factors"} plots the estimated factor and time fixed effects,
    \code{"loadings"} plots the estimated factor loadings and unit fixed effects,
    \code{"calendar"} plots the estimated treatment effects for each calendar period,
    and \code{"box"} visualizes the estimated individualistic treatment effects of observations.}
  \item{loo}{a logical flag indicating whether to use the leave-one-period-out pre-treatment effects for the visualization and tests.}  
  \item{highlight}{a logical flag indicating whether to highlight the periods for the carryover and placebo test.}
  \item{plot.ci}{a string specifying the confidence interval. Choose from: "0.9", "0.95", or "none".}
  \item{show.points}{a logical flag indicating whether to represent treatment effects by points or point-ranges. Default to TRUE.}
  \item{show.group}{a string indicating the group to be visualized.}
  \item{bound}{a string that specifies the bounds to be plotted for equivalence test 
    for pre-treatment fit checking. Choose from: "both", "equiv", "min" and "none".}
  \item{vis}{A string specifying whether to plot the dots for placebo plots.}
  \item{count}{a logical flag controlling whether to show the count of each term 
    for gap plot.}
  \item{proportion}{a positive value specifying periods at which observations equal to 
    or greater than the proporation of the largest number of observations at a certain 
    period. Default to 0.3.}
  \item{pre.periods}{a vector specifying the range of pre-treatment period used for goodness-of-fit test. If left blank, all 
    pre-treatment periods specified by \code{proportion} will be used. Ignore if \code{se = FALSE}.}
  \item{f.threshold}{a numeric value specifying the threshold for the F-statistic in the equivalent test.
    Ignore if \code{se = FALSE}.  Deafult is \code{f.threshold = 0.5}.}
  \item{tost.threshold}{a numeric value specifying the threshold for the two-one-sided t-test.
    If \code{alpha=0.05}, TOST checks whether the 90% confidence intervals for the estimated period-wise ATT exceeds this threshold.
    The default value is 0.36 times the standard deviation of the outcome variable after two-way fixed effects are partialed out.}
  \item{effect.bound.ratio}{a logical value specifiying whether to annotate 
    the ratio of estimated average treatment effects / minimun bound.}
  \item{stats}{a string that specifices what statistics to be shown. For \code{"gap"} plot, choose from \code{c("none", "F.p", "equiv.p", "F.equiv.p",)}, for 
  \code{"placebo"} plot, choose from \code{c("none", "placebo.p", "equiv.placebo.p")}, and for "carryover" plot, choose from \code{c("none", "carryover.p")}.}
  \item{stats.labs}{a string specifying the label for the statistics specified by \code{stats}.}
  \item{main}{a string that controls the title of the plot. If not
    supplied, no title will be shown.}
  \item{xlim}{a two-element numeric vector specifying the range of x-axis. When 
    class of time varible is string, must specify not original value but a counting number e.g. \code{xlim=c(1,30)}.}
  \item{ylim}{a two-element numeric vector specifying the range of y-axis.}
  \item{xlab}{a string indicating the label of the x-axis.}
  \item{ylab}{a string indicating the label of the y-axis.}
  \item{gridOff}{a logical flag indicating whether to remove the grid lines for the status plot.}
  \item{legendOff}{a logical flag controlling whether to show the legend.}
  \item{legend.pos}{a string specifying the position of legend. If left blank, 
    legend will be setted at the bottom.}
  \item{legend.nrow}{an integer speficying rows of legend.}
  \item{legend.labs}{a string vector for user-defined legends.}
  \item{stats.pos}{a numeric vector of length 2 specifying the postion for 
    labels of test statistic value.}
  \item{theme.bw}{a logical flag specifying whether to use the black and 
    white theme.}
  \item{nfactors}{a integer controlling the number of factors to be shown when \code{type='factors'}.}
  \item{include.FE}{a loical flag indicating whether to keep the fixed effects when \code{type='factors'} or \code{type='loadings'}.}
  \item{id}{a string vector specifying a sub-group of units that are to be 
    plotted for treatment status(\code{type = "status"}).}
  \item{cex.main}{a numeric value (pt) specifying the fontsize of the main title.}
  \item{cex.main.sub}{a numeric value (pt) specifying the fontsize of the subtitles.}
  \item{cex.axis}{a numeric value (pt) specifying the fontsize of the texts on the axes.}
  \item{cex.lab}{a numeric value (pt) specifying the fontsize of the axis titles.}
  \item{cex.legend}{a numeric value (pt) specifying the fontsize of the legend.}
  \item{cex.text}{a numeric value (pt) specifying the fontsize of the annotations.}
  \item{axis.adjust}{a logic flag indicating whether to adjust labels on x-axis. 
    Useful when class of time variable is string and data magnitude is large.}
  \item{axis.lab}{a string indicating whether labels on the x- and y-axis will be shown. There are four options: \code{"both"} (default): labels on both axes will be shown; \code{"unit"}: only labels on y-axis will be shown; \code{"time"}: only labels on the x-axis will be shown; "none": no labels will be shown.}
  \item{axis.lab.gap}{a numeric vector setting the gaps between labels on the x- or y-axis for \code{"missing"} plot. Default is \code{axis.lab.gap = c(0, 0)}, which means that all labels will be shown. Useful for datasets with large N or T.}
  \item{start0}{a logical flag indicating whether to index the start of the treatment as period 0 rather than period 1. Default to FALSE.}  
  \item{return.test}{a logical flag indicating whether to return the results of statistical tests.}
  \item{balance}{a logical flag indicating whether to plot the dynamic treatment effects for the balance sample.}
  \item{weight}{a logical flag indicating whether to plot the weighted dynamic treatment effects.}
  \item{\dots}{other argv.}
}
\details{
  \code{plot.fect} visualizes the estimation results obtained by \code{fect}.
}
\value{
  \item{p}{a ggplot2 object saving the graph.}
  \item{test.out}{a list storing the results of statistical tests if \code{return.test=TRUE}.}
=======
  \item{x}{A fitted \code{\link{fect}} object.}
  \item{type}{Plot type: \code{"gap"}, \code{"equiv"}, \code{"status"}, \code{"exit"},
    \code{"factors"}, \code{"loadings"}, \code{"calendar"}, \code{"box"}, or \code{"counterfactual"}.}
  \item{loo}{Logical; if \code{TRUE}, use leave-one-out for equivalence tests.}
  \item{highlight}{Logical; if \code{TRUE}, highlight certain periods (e.g., placebo).}
  \item{plot.ci}{What confidence interval to show: e.g., \code{"0.9"}, \code{"0.95"}, or \code{"none"}.}
  \item{show.points}{Logical; if \code{TRUE}, show point estimates on the plot.}
  \item{show.group}{Optional; specify a subgroup to plot.}
  \item{bound}{Which bounds to display in bounding/equivalence tests: \code{"none"}, \code{"min"}, \code{"equiv"}, \code{"both"}.}
  \item{vis}{For placebo/carryover tests, how lines are displayed: \code{"connected"} or \code{"none"}.}
  \item{count}{Logical; if \code{TRUE}, shows a histogram of the number of observations.}
  \item{proportion}{Restrict time points to those with at least this fraction of the max units.}
  \item{pre.periods}{Specific pre-treatment periods used for bounding.}
  \item{f.threshold}{F-test threshold in equivalence checks.}
  \item{tost.threshold}{Threshold for TOST-based checks.}
  \item{effect.bound.ratio}{If \code{TRUE}, show ratio of ATT to the minimal bound.}
  \item{stats}{Which stats to display, e.g. \code{"F.p"}, \code{"equiv.p"}, etc.}
  \item{stats.labs}{Custom labels for the displayed stats.}
  \item{raw}{For \code{"counterfactual"} plot: \code{"none"}, \code{"band"}, or \code{"all"}.}
  \item{main}{Main plot title.}
  \item{xlim}{x-axis limits.}
  \item{ylim}{y-axis limits.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{xangle}{Angle in degrees for x-axis labels.}
  \item{yangle}{Angle in degrees for y-axis labels.}
  \item{xbreaks}{Major x-axis tick locations.}
  \item{ybreaks}{Major y-axis tick locations.}
  \item{xticklabels}{Labels for x-axis ticks.}
  \item{yticklabels}{Labels for y-axis ticks.}
  \item{gridOff}{Logical; if \code{TRUE}, removes grid lines.}
  \item{legendOff}{Logical; if \code{TRUE}, hides the legend.}
  \item{legend.pos}{Legend position.}
  \item{legend.nrow}{Number of rows in legend.}
  \item{legend.labs}{Custom legend labels.}
  \item{stats.pos}{Position (x,y) for test statistics text.}
  \item{show.stats}{Logical; if \code{TRUE}, shows test statistic text on plot.}
  \item{theme.bw}{Logical; if \code{TRUE}, uses a black-and-white ggplot2 theme.}
  \item{nfactors}{Number of factors to plot for \code{"factors"} or \code{"loadings"}.}
  \item{include.FE}{Whether to include fixed effects as "factors".}
  \item{id}{Vector specifying a subset of units to plot.}
  \item{cex.main}{Size of main title.}
  \item{cex.main.sub}{Size of a possible subtitle.}
  \item{cex.axis}{Size of axis tick labels.}
  \item{cex.lab}{Size of axis labels.}
  \item{cex.legend}{Size of legend text.}
  \item{cex.text}{Size of text annotations.}
  \item{axis.adjust}{If \code{TRUE}, attempts to adjust axis labels for large data.}
  \item{axis.lab}{Controls axis labeling in \code{"status"} plots. Could be \code{"both"}, \code{"time"}, \code{"unit"}, or \code{"off"}.}
  \item{axis.lab.gap}{Controls spacing for axis labels in \code{"status"} plots.}
  \item{shade.post}{If \code{TRUE}, shades post-treatment periods in some plot types.}
  \item{start0}{If \code{TRUE}, re-index first post-treatment period as 0.}
  \item{return.test}{If \code{TRUE}, returns test statistics as well as the plot.}
  \item{balance}{If \code{TRUE}, uses "balance sample" results.}
  \item{weight}{If \code{TRUE}, uses weighting from the \code{fect} object.}
  \item{line.color}{Color for main lines.}
  \item{line.width}{Numeric vector specifying line widths.}
  \item{lcolor}{Color for reference lines.}
  \item{lwidth}{Width of reference lines.}
  \item{...}{Additional graphical parameters.}
}
\details{
  \code{plot.fect} can produce various diagnostic or substantive plots for objects from \code{\link{fect}}. Depending on \code{type}, it shows treatment effects, equivalence tests, factor/loadings, or raw vs. counterfactual outcomes.
}
\value{
  \item{p}{A \code{ggplot2} object with the requested plot.}
  \item{test.out}{If \code{return.test=TRUE}, includes test statistics (placebo, F-test, equivalence p-values, etc.).}
>>>>>>> Stashed changes
}
\author{
  Licheng Liu; Ye Wang; Yiqing Xu, Ziyi Liu
}
\references{
  Jushan Bai. 2009. "Panel Data Models with Interactive Fixed Effects." Econometrica.
  
  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models." Political 
  Analysis. 
  
  Athey, Susan, et al. 2021 "Matrix completion methods for causal panel data models." Journal of the American Statistical Association. 
  
  Licheng Liu, et al. 2022. "A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional 
  Data." American Journal of Political Science. 

  For more details about the matrix completion method, see \url{https://github.com/susanathey/MCPanel}. 
}
\seealso{
  \code{\link{fect}} and \code{\link{print.fect}}
}
\examples{
library(fect)
data(fect)
out <- fect(Y ~ D + X1 + X2, data = simdata1, 
            index = c("id","time"), force = "two-way",
            CV = TRUE, r = c(0, 5), se = 0, parallel = FALSE) 
plot(out)
plot(out,type='status')
plot(out,show.points = FALSE)
}



