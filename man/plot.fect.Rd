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
    restrict = "rm",
    loo = FALSE,
    highlight = NULL,
    plot.ci = NULL,
    show.points = TRUE,
    show.group = NULL,
    bound = NULL,
    show.count = TRUE,
    proportion = 0.3,
    pre.periods = NULL,
    f.threshold = NULL,
    tost.threshold = NULL,
    effect.bound.ratio = FALSE,
    stats = NULL,
    stats.labs = NULL,
    raw = "none",
    vis = NULL,
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
    gridOff = NULL,
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
    lcolor = NULL,
    lwidth = NULL,
    line.color = NULL,
    line.width = NULL,
    theme = NULL,
    connected = NULL,
    ci.outline = FALSE,
    color =  NULL,
    est.lwidth = 0.8,
    est.pointsize = 3,
    count.color =  NULL,
    placebo.color =  NULL,
    carryover.color =  NULL,
    carryover.rm.color = NULL,
    sens.original.color = NULL,
    sens.colors =  NULL,
    counterfactual.color =  NULL,
    counterfactual.raw.controls.color =  NULL,
    counterfactual.raw.treated.color =  NULL,
    counterfactual.linetype =  NULL,
    box.control =  NULL,
    box.treat =  NULL,
    calendar.color =  NULL,
    calendar.lcolor =  NULL,
    equiv.color =  NULL,
    status.treat.color =  NULL,
    status.control.color = NULL,
    status.missing.color =  NULL,
    status.removed.color =  NULL,
    status.placebo.color = NULL,
    status.carryover.color =  NULL,
    status.carryover.rm.color =  NULL,
    status.balanced.post.color = NULL,
    status.balanced.pre.color = NULL,
    status.background.color = NULL,
    ...
)
}
\arguments{
  \item{x}{A fitted \code{\link{fect}} object.}
  \item{type}{Plot type. Options include: \code{"gap"}, \code{"equiv"}, \code{"status"}, \code{"exit"},
    \code{"factors"}, \code{"loadings"}, \code{"calendar"}, \code{"box"}, \code{"counterfactual"},
    \code{"sens"} (sensitivity analysis plot), \code{"sens_es"} (event-study style sensitivity plot),
    or \code{"cumul"} (cumulative effect plot). Default is context-dependent.}
  \item{restrict}{For sensitivity plots (\code{type = "sens"} or \code{type = "sens_es"}). Specifies the type of restriction: \code{"rm"} (restriction on M) or \code{"sm"} (smoothness restriction). Default is \code{"rm"}.}
  \item{loo}{Logical; if \code{TRUE} (default is \code{FALSE}), use leave-one-out for equivalence tests.}
  \item{highlight}{Logical or \code{NULL} (default). If \code{TRUE}, highlight certain periods (e.g., placebo or carryover periods in relevant tests). If \code{NULL}, defaults to \code{TRUE} if a placebo or carryover test is being plotted, \code{FALSE} otherwise.}
  \item{plot.ci}{Specifies the confidence interval to be plotted. Options include:
    \code{"0.9"}, \code{"0.95"}, or \code{"none"}. Default is \code{NULL}, which is context-dependent (e.g., "0.95" for gap plots, "0.9" for equivalence plots, "none" if no CIs available).}
  \item{show.points}{Logical; if \code{TRUE} (default), shows point estimates on event study style plots.}
  \item{show.group}{Optional character string or \code{NULL} (default); specify a subgroup (defined in the \code{fect} call) to plot.}
  \item{bound}{Which bounds to display in bounding/equivalence tests:
    \code{"none"} (default for most plots), \code{"min"}, \code{"equiv"}, or \code{"both"} (default for \code{type="equiv"}).}
  \item{show.count}{Logical; if \code{TRUE} (default), shows a histogram of observation counts in relevant plots (e.g., gap, calendar, box plots).}
  \item{proportion}{Numeric (0 to 1); restrict plotted time points in event study plots to those with at least this fraction of the maximum number of treated units. Default is \code{0.3}.}
  \item{pre.periods}{Optional numeric vector; specifies pre-treatment periods used for bounding and equivalence tests. Default is \code{NULL}, using periods from the \code{fect} object.}
  \item{f.threshold}{Numeric or \code{NULL} (default); F-test threshold for equivalence checks. If \code{NULL}, uses value from the \code{fect} object.}
  \item{tost.threshold}{Numeric or \code{NULL} (default); Threshold for TOST-based equivalence checks (e.g., for pre-trend, placebo, carryover). If \code{NULL}, uses value from the \code{fect} object.}
  \item{effect.bound.ratio}{Logical; if \code{TRUE} (default is \code{FALSE}), shows the ratio of ATT to the minimal bound in equivalence plots.}
  \item{stats}{Character vector or \code{NULL} (default); specifies which test statistics to display on the plot (e.g., \code{"F.p"}, \code{"equiv.p"}, \code{"placebo.p"}). Default is context-dependent.}
  \item{stats.labs}{Character vector or \code{NULL} (default); custom labels for the displayed statistics. Must match the length of \code{stats}.}
  \item{raw}{For \code{type = "counterfactual"} plot: \code{"none"} (default), \code{"band"} (shows 5-95\% quantile band for raw control/treated data), or \code{"all"} (shows all raw control/treated unit trajectories).}
  \item{vis}{(Deprecated) Formerly controlled line display for placebo/carryover tests. This is now handled by \code{highlight} and internal logic.}
  \item{main}{Character string or \code{NULL} (default); main plot title. If \code{NULL}, a default title is used based on plot type.}
  \item{xlim}{Numeric vector of length 2 or \code{NULL} (default); specifies x-axis limits.}
  \item{ylim}{Numeric vector of length 2 or \code{NULL} (default); specifies y-axis limits.}
  \item{xlab}{Character string or \code{NULL} (default); x-axis label. If \code{NULL}, a default label is used.}
  \item{ylab}{Character string or \code{NULL} (default); y-axis label. If \code{NULL}, a default label is used.}
  \item{xangle}{Numeric or \code{NULL} (default); angle (in degrees) for x-axis tick labels.}
  \item{yangle}{Numeric or \code{NULL} (default); angle (in degrees) for y-axis tick labels.}
  \item{xbreaks}{Numeric vector or \code{NULL} (default); positions of major ticks on the x-axis.}
  \item{ybreaks}{Numeric vector or \code{NULL} (default); positions of major ticks on the y-axis.}
  \item{xticklabels}{Character vector or \code{NULL} (default); labels for the x-axis ticks. Used in status plots.}
  \item{yticklabels}{Character vector or \code{NULL} (default); labels for the y-axis ticks. Used in status plots.}
  \item{gridOff}{Logical or \code{NULL} (default). If \code{TRUE}, removes major and minor grid lines. If \code{NULL}, grid is off for most plots but on for \code{type = "status"} (tile outlines).}
  \item{legendOff}{Logical; if \code{TRUE} (default is \code{FALSE}), hides the legend.}
  \item{legend.pos}{Character string or \code{NULL} (default); position of the legend (e.g., "bottom", "top", "none").}
  \item{legend.nrow}{Integer or \code{NULL} (default); number of rows in the legend.}
  \item{legend.labs}{Character vector or \code{NULL} (default); custom legend labels.}
  \item{stats.pos}{Numeric vector of length 2 or \code{NULL} (default); (x,y) coordinates for displaying text of test statistics.}
  \item{show.stats}{Logical; if \code{TRUE} (default), displays test statistics on the plot if \code{stats} is not "none".}
  \item{theme.bw}{Logical; if \code{TRUE} (default), uses a black-and-white ggplot2 theme (\code{theme_bw()}).}
  \item{nfactors}{Integer or \code{NULL} (default); number of factors to plot when \code{type} is \code{"factors"} or \code{"loadings"}. Defaults to min(r.cv, 4).}
  \item{include.FE}{Logical; if \code{TRUE} (default), includes fixed effects as additional "factors" if applicable.}
  \item{id}{Character vector or \code{NULL} (default); when specified, plots only the chosen unit(s). For counterfactual plots, only one ID is supported.}
  \item{cex.main}{Numeric or \code{NULL} (default); font size multiplier for the main title. Base size is 16pt.}
  \item{cex.main.sub}{Numeric or \code{NULL} (default); font size multiplier for a possible subtitle. Base size is 16pt.}
  \item{cex.axis}{Numeric or \code{NULL} (default); font size multiplier for axis tick labels. Base size is 15pt.}
  \item{cex.lab}{Numeric or \code{NULL} (default); font size multiplier for axis labels. Base size is 15pt.}
  \item{cex.legend}{Numeric or \code{NULL} (default); font size multiplier for legend text. Base size is 15pt.}
  \item{cex.text}{Numeric or \code{NULL} (default); font size multiplier for text annotations (e.g., stats values). Base size is 5pt.}
  \item{axis.adjust}{Logical; if \code{TRUE} (default is \code{FALSE}), attempts to adjust axis labels for large data (e.g., rotates x-axis labels).}
  \item{axis.lab}{Controls axis labeling in \code{type="status"} plots. Options: \code{"both"} (default), \code{"time"}, \code{"unit"}, or \code{"off"}.}
  \item{axis.lab.gap}{Numeric vector of length 1 or 2; gap spacing for axis labels in \code{type="status"} plots. Default is \code{c(0, 0)}.}
  \item{shade.post}{Logical; if \code{FALSE} (default). If \code{TRUE}, shades post-treatment periods in certain plots (e.g., counterfactual plot where it becomes \code{TRUE} if initially \code{NULL} and not explicitly set).}
  \item{start0}{Logical; if \code{TRUE} (default is \code{FALSE}), re-indexes the first post-treatment period as 0 in event-time plots.}
  \item{return.test}{Logical; if \code{TRUE} (default is \code{FALSE}), returns a list containing the plot object and relevant test statistics.}
  \item{balance}{Logical or \code{NULL} (default); if \code{TRUE}, uses a "balance sample" from the \code{fect} object if available.}
  \item{weight}{Logical or \code{NULL} (default); if \code{TRUE}, applies weighting from the \code{fect} object if available.}
  \item{lcolor}{Character vector of length 1 or 2, or \code{NULL} (default); color(s) for reference lines (e.g., horizontal y=0 line, vertical pre/post treatment line). Default is theme-dependent.}
  \item{lwidth}{Numeric vector of length 1 or 2, or \code{NULL} (default); line width(s) for reference lines. Default is theme-dependent (e.g., 1 or 1.5).}
  \item{line.color}{(Deprecated) Use \code{color} for main estimate lines or \code{lcolor} for reference lines instead.}
  \item{line.width}{(Deprecated) Use \code{est.lwidth} for main estimate lines or \code{lwidth} for reference lines instead.}
  \item{theme}{A character string specifying a color theme: \code{NULL} (default), \code{"vibrant"}, or \code{"grayscale"}. This sets a predefined group of colors for many plot elements. Individual color arguments below can override theme settings.}
  \item{connected}{Logical or \code{NULL} (default). For event study type plots (e.g., \code{type="gap"}), if \code{TRUE}, connects points with lines; if \code{FALSE}, shows point-ranges. If \code{NULL}, behavior is theme-dependent (default \code{FALSE} for base theme, \code{TRUE} for "vibrant").}
  \item{ci.outline}{Logical; if \code{FALSE} (default). For event study plots where CIs are shown as ribbons, if \code{TRUE}, adds an outline to the CI ribbon.}
  \item{color}{Character string or \code{NULL} (default). Color for main plot lines/points (e.g., ATT line). Default is theme-dependent (e.g., "black").}
  \item{est.lwidth}{Numeric; line width for the main estimate line in event study plots. Default is \code{0.8}.}
  \item{est.pointsize}{Numeric; size of points for the main estimate in event study plots. Default is \code{3}.}
  \item{count.color}{Color for count bars/histogram. Default is \code{NULL} (theme-dependent).}
  \item{placebo.color}{Color for placebo test elements. Default is \code{NULL} (theme-dependent).}
  \item{carryover.color}{Color for carryover test elements. Default is \code{NULL} (theme-dependent).}
  \item{carryover.rm.color}{Color for removed carryover elements. Default is \code{NULL} (theme-dependent).}
  \item{sens.original.color}{Color for original estimate in sensitivity plots. Default is \code{NULL} (theme-dependent).}
  \item{sens.colors}{Vector of colors for different M values in sensitivity plots. Default is \code{NULL} (theme-dependent).}
  \item{counterfactual.color}{Color for counterfactual line. Default is \code{NULL} (theme-dependent).}
  \item{counterfactual.raw.controls.color}{Color for raw control lines in counterfactual plot. Default is \code{NULL} (theme-dependent).}
  \item{counterfactual.raw.treated.color}{Color for raw treated lines in counterfactual plot. Default is \code{NULL} (theme-dependent).}
  \item{counterfactual.linetype}{Linetype for counterfactual line. Default is \code{NULL} (theme-dependent, e.g., "dashed").}
  \item{box.control}{Fill color for control group boxplots. Default is \code{NULL} (theme-dependent).}
  \item{box.treat}{Fill color for treated group boxplots. Default is \code{NULL} (theme-dependent).}
  \item{calendar.color}{Color for fitted line in calendar time plots. Default is \code{NULL} (theme-dependent).}
  \item{calendar.lcolor}{Color for average ATT line in calendar time plots. Default is \code{NULL} (theme-dependent).}
  \item{equiv.color}{Color for equivalence bounds. Default is \code{NULL} (theme-dependent).}
  \item{status.treat.color}{Color for treated status in status plots. Default is \code{NULL} (theme-dependent).}
  \item{status.control.color}{Color for control status in status plots. Default is \code{NULL} (theme-dependent).}
  \item{status.missing.color}{Color for missing status in status plots. Default is \code{NULL} (theme-dependent).}
  \item{status.removed.color}{Color for removed status in status plots. Default is \code{NULL} (theme-dependent).}
  \item{status.placebo.color}{Color for placebo status in status plots. Default is \code{NULL} (theme-dependent).}
  \item{status.carryover.color}{Color for carryover status in status plots. Default is \code{NULL} (theme-dependent).}
  \item{status.carryover.rm.color}{Color for removed carryover status in status plots. Default is \code{NULL} (theme-dependent).}
  \item{status.balanced.post.color}{Color for balanced post-treatment status. Default is \code{NULL} (theme-dependent).}
  \item{status.balanced.pre.color}{Color for balanced pre-treatment status. Default is \code{NULL} (theme-dependent).}
  \item{status.background.color}{Background color for status plots. Default is \code{NULL} (theme-dependent, e.g., "white").}
  \item{...}{Additional graphical parameters passed to internal plotting routines, primarily those accepted by \code{esplot} for event study style plots.}
}
\details{
  \code{plot.fect} generates various visualizations for objects estimated via \code{\link{fect}}. Depending on
  the selected \code{type}, it can show treatment-effect dynamics (\code{"gap"}, \code{"exit"}, \code{"cumul"}),
  equivalence test outcomes (\code{"equiv"}), factor/loadings for factor models (\code{"factors"}, \code{"loadings"}),
  raw vs. counterfactual trajectories (\code{"counterfactual"}), treatment status maps (\code{"status"}),
  effects by calendar time (\code{"calendar"}), distribution of individual effects (\code{"box"}),
  or sensitivity analysis results (\code{"sens"}, \code{"sens_es"}).
  The \code{theme} argument allows for quick customization of plot aesthetics using predefined color schemes.
  Many individual color parameters can be set to further customize the appearance.
}
\value{
  A \code{ggplot2} object representing the plot. If \code{return.test=TRUE}, a list containing:
  \item{p}{The \code{ggplot2} object.}
  \item{test.out}{Relevant test statistics (F-tests, equivalence p-values, etc.), if applicable.}
}
\author{
  Licheng Liu, Ye Wang, Yiqing Xu, Ziyi Liu
}
\references{
  Xu, Y. (2017). Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models.
  \emph{Political Analysis, 25}(1), 57–76.

  Liu, L., Wang, Y., & Xu, Y. (2022). A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional Data.
  \emph{American Journal of Political Science, 66}(1), 220-237. (Provides context for \code{esplot} which is used internally)

  Rambachan, A., & Roth, J. (2023). A More Credible Approach to Parallel Trends.
  \emph{Review of Economic Studies, 90}(5), 2555-2591. (Provides context for sensitivity analysis plots)
}
\examples{
library(fect)
# For CRAN checks, use a small number of bootstraps
# In practice, use a larger number (e.g., nboots = 200 or more)
if(requireNamespace("ggplot2") && requireNamespace("ggrepel")) {
  data(simdata)
  # Estimate with fixed effects method
  out.fect <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
    method = "fe", force = "two-way", se = TRUE, parallel = FALSE, nboots = 5) # nboots low for example

  # Default gap plot
  plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y")

  # Counterfactual plot for the first treated unit
  # Need to know the ID of a treated unit. Let's find one.
  treated_ids <- unique(simdata$id[simdata$D == 1])
  if (length(treated_ids) > 0) {
    plot(out.fect, type = "counterfactual", id = treated_ids[1],
         main = paste("Counterfactual for Unit", treated_ids[1]))
  }

  # Status plot
  plot(out.fect, type = "status")

  # Cumulative effect plot (if est.eff is available from fect call)
  # This example might not have it by default, but showing how to call
  # out.fect.cumul <- fect(Y ~ D, data = simdata, index = c("id","time"), method = "fe",
  #                        cumulative = TRUE, se = TRUE, parallel = FALSE, nboots = 5)
  # if (exists("out.fect.cumul")) {
  #   plot(out.fect.cumul, type = "cumul", main = "Cumulative ATT")
  # }
}
}
