\encoding{UTF-8}
\name{plot.fect}
\alias{plot.fect}
\title{Plot Method for 'fect' Objects}
\description{
  Visualizes results from a \code{\link{fect}} estimation.
}
\usage{
\method{plot}{fect}(
    x,
    type = NULL,
    loo = FALSE,
    highlight = NULL,
    plot.ci = NULL,
    show.points = NULL,
    show.group = NULL,
    bound = NULL,
    vis = NULL,
    count = TRUE,
    proportion = 0.3,
    pre.periods = NULL,
    f.threshold = NULL,
    tost.threshold = NULL,
    effect.bound.ratio = FALSE,
    stats = NULL,
    stats.labs = NULL,
    raw = "none",
    main = NULL,
    xlim = NULL,
    ylim = NULL,
    xlab = NULL,
    ylab = NULL,
    xangle = NULL,
    yangle = NULL,
    xbreaks = NULL,
    ybreaks = NULL,
    xticklabels = NULL,
    yticklabels = NULL,
    gridOff = FALSE,
    legendOff = FALSE,
    legend.pos = NULL,
    legend.nrow = NULL,
    legend.labs = NULL,
    stats.pos = NULL,
    show.stats = TRUE,
    theme.bw = TRUE,
    nfactors = NULL,
    include.FE = TRUE,
    id = NULL,
    cex.main = NULL,
    cex.main.sub = NULL,
    cex.axis = NULL,
    cex.lab = NULL,
    cex.legend = NULL,
    cex.text = NULL,
    axis.adjust = FALSE,
    axis.lab = "both",
    axis.lab.gap = c(0, 0),
    shade.post = FALSE,
    start0 = FALSE,
    return.test = FALSE,
    balance = NULL,
    weight = NULL,
    line.color = NULL,
    line.width = c(1.2, 0.5),
    lcolor = NULL,
    lwidth = NULL,
    ...
)
}
\arguments{
  \item{x}{A fitted \code{\link{fect}} object.}
  \item{type}{Plot type: \code{"gap"}, \code{"equiv"}, \code{"status"}, \code{"exit"},
    \code{"factors"}, \code{"loadings"}, \code{"calendar"}, \code{"box"}, or \code{"counterfactual"}.}
  \item{loo}{Logical; if \code{TRUE}, use leave-one-out for equivalence tests.}
  \item{highlight}{Logical; if \code{TRUE}, highlight certain periods (e.g., placebo).}
  \item{plot.ci}{Specifies the confidence interval to be plotted. Options include:
    \code{"0.9"}, \code{"0.95"}, or \code{"none"}. Default depends on the context.}
  \item{show.points}{Logical; if \code{TRUE}, shows point estimates on the plot.}
  \item{show.group}{Optional; specify a subgroup to plot.}
  \item{bound}{Which bounds to display in bounding/equivalence tests:
    \code{"none"}, \code{"min"}, \code{"equiv"}, or \code{"both"}.}
  \item{vis}{Controls line display for placebo/carryover tests: \code{"connected"} or \code{"none"}.}
  \item{count}{Logical; if \code{TRUE}, shows a histogram of observation counts.}
  \item{proportion}{Restrict plotted time points to those with at least this fraction of the max units.}
  \item{pre.periods}{Specifies pre-treatment periods used for bounding.}
  \item{f.threshold}{F-test threshold in equivalence checks.}
  \item{tost.threshold}{Threshold for TOST-based checks.}
  \item{effect.bound.ratio}{If \code{TRUE}, shows the ratio of ATT to the minimal bound.}
  \item{stats}{Specifies which statistics to display: e.g., \code{"F.p"}, \code{"equiv.p"}.}
  \item{stats.labs}{Custom labels for the displayed statistics.}
  \item{raw}{For \code{"counterfactual"} plot: \code{"none"}, \code{"band"}, or \code{"all"}.}
  \item{main}{Main plot title.}
  \item{xlim}{Specifies x-axis limits.}
  \item{ylim}{Specifies y-axis limits.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{xangle}{Angle (in degrees) for x-axis labels.}
  \item{yangle}{Angle (in degrees) for y-axis labels.}
  \item{xbreaks}{Positions of major ticks on the x-axis.}
  \item{ybreaks}{Positions of major ticks on the y-axis.}
  \item{xticklabels}{Labels for the x-axis ticks.}
  \item{yticklabels}{Labels for the y-axis ticks.}
  \item{gridOff}{Logical; if \code{TRUE}, removes grid lines.}
  \item{legendOff}{Logical; if \code{TRUE}, hides the legend.}
  \item{legend.pos}{Position of the legend.}
  \item{legend.nrow}{Number of rows in the legend.}
  \item{legend.labs}{Custom legend labels.}
  \item{stats.pos}{(x,y) coordinates for displaying text of test statistics.}
  \item{show.stats}{Logical; if \code{TRUE}, displays test statistics on the plot.}
  \item{theme.bw}{Logical; if \code{TRUE}, uses a black-and-white ggplot2 theme.}
  \item{nfactors}{Number of factors to plot when \code{type} is \code{"factors"} or \code{"loadings"}.}
  \item{include.FE}{Logical; if \code{TRUE}, includes fixed effects as additional "factors".}
  \item{id}{When specified, plots only the chosen unit(s).}
  \item{cex.main}{Font size of the main title.}
  \item{cex.main.sub}{Font size of a possible subtitle.}
  \item{cex.axis}{Font size of axis tick labels.}
  \item{cex.lab}{Font size of axis labels.}
  \item{cex.legend}{Font size of legend text.}
  \item{cex.text}{Font size for text annotations.}
  \item{axis.adjust}{If \code{TRUE}, attempts to adjust axis labels for large data.}
  \item{axis.lab}{Controls axis labeling in \code{"status"} plots. Could be \code{"both"},
    \code{"time"}, \code{"unit"}, or \code{"off"}.}
  \item{axis.lab.gap}{Gap spacing for axis labels in \code{"status"} plots.}
  \item{shade.post}{If \code{TRUE}, shades post-treatment periods in certain plots.}
  \item{start0}{If \code{TRUE}, re-indexes the first post-treatment period as 0.}
  \item{return.test}{If \code{return.test=TRUE}, returns test statistics along with the plot.}
  \item{balance}{Logical; if \code{TRUE}, uses a "balance sample" from the \code{fect} object.}
  \item{weight}{Logical; if \code{TRUE}, applies weighting from the \code{fect} object.}
  \item{line.color}{Color for main plot lines.}
  \item{line.width}{A numeric vector specifying line widths for different line types.}
  \item{lcolor}{Color for reference lines (e.g., 0 line, vertical lines).}
  \item{lwidth}{Line width for reference lines.}
  \item{...}{Additional graphical parameters passed to internal plotting routines.}
}
\details{
  \code{plot.fect} generates various visualizations for objects estimated via \code{\link{fect}}. Depending on 
  the selected \code{type}, it can show treatment-effect dynamics, equivalence test outcomes, factor/loadings 
  for factor models, or raw vs. counterfactual trajectories. 
}
\value{
  \item{p}{A \code{ggplot2} object containing the requested plot.}
  \item{test.out}{If \code{return.test=TRUE}, includes relevant test statistics (F-tests, equivalence p-values, etc.).}
}
\author{
  Licheng Liu, Ye Wang, Yiqing Xu, Ziyi Liu
}
\references{
  Xu, Y. (2017). Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models.
  \emph{Political Analysis, 25}(1), 57–76.
}
\examples{
library(fect)
data(fect)
out.fect <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"), 
  method = "fe", force = "two-way", se = TRUE, parallel = TRUE, nboots = 200)
plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y")
}
