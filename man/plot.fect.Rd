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
    return.data = FALSE,
    balance = NULL,
    weight = NULL,
    lcolor = NULL,
    lwidth = NULL,
    ltype = NULL,
    line.color = NULL,
    line.width = NULL,
    count = NULL,
    preset = NULL,
    connected = NULL,
    ci.outline = FALSE,
    color =  NULL,
    est.lwidth = NULL,
    est.pointsize = NULL,
    count.color =  NULL,
    count.alpha = NULL,
    count.outline.color = NULL,
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
    calendar.color = NULL,
    calendar.lcolor = NULL,
    calendar.cicolor = NULL,
    heterogeneous.color = NULL,
    heterogeneous.lcolor = NULL,
    heterogeneous.cicolor = NULL,
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
    covariate = NULL,
    covariate.labels = NULL,
    ...
)
}
\arguments{
  \item{x}{A fitted \code{\link{fect}} object.}
  \item{type}{Plot type. Options include: \code{"gap"}, \code{"equiv"} (equivalence test), \code{"status"} (treatment status), \code{"exit"} (exiting treatment effects),
    \code{"factors"}, \code{"loadings"}, \code{"calendar"} (ATT by calendar time), \code{"box"} (individual effects distribution), \code{"counterfactual"},
    \code{"sens"} (sensitivity analysis plot, e.g., for Rambachan & Roth bounds), \code{"sens_es"} (event-study style sensitivity plot),
    or \code{"cumul"} (cumulative effect plot). Default is context-dependent (e.g., "gap", or "equiv" if \code{loo=TRUE}).}
  \item{restrict}{For sensitivity plots (\code{type = "sens"} or \code{type = "sens_es"}). Specifies the type of restriction: \code{"rm"} (restriction on M, for relative magnitude bounds) or \code{"sm"} (smoothness restriction, for bounds based on second differences). Default is \code{"rm"}.}
  \item{loo}{Logical; if \code{TRUE} (default is \code{FALSE}), use leave-one-out estimates for pre-treatment periods in equivalence tests or gap plots.}
  \item{highlight}{Logical or \code{NULL} (default). If \code{TRUE}, highlights specific periods such as placebo or carryover periods in relevant tests. If \code{NULL}, defaults to \code{TRUE} if a placebo or carryover test is being plotted, and \code{FALSE} otherwise.}
  \item{plot.ci}{Specifies the confidence interval to be plotted. Options include:
    \code{"0.9"}, \code{"0.95"}, or \code{"none"}. Default is \code{NULL}, which is context-dependent (e.g., "0.95" for gap plots, "0.9" for equivalence plots, "none" if no CIs are available in the \code{fect} object).}
  \item{show.points}{Logical; if \code{TRUE} (default), shows point estimates on event study style plots (e.g., gap, equiv, exit plots).}
  \item{show.group}{Optional character string or \code{NULL} (default); if specified, plots results for a particular subgroup defined in the \code{fect} call.}
  \item{bound}{Which bounds to display in bounding/equivalence tests:
    \code{"none"} (default for most plots), \code{"min"} (minimal effect bound), \code{"equiv"} (equivalence range), or \code{"both"} (default for \code{type="equiv"}).}
  \item{show.count}{Logical; if \code{TRUE} (default), shows a histogram of observation counts (number of treated units) in relevant plots like gap, calendar, or box plots.}
  \item{proportion}{Numeric (0 to 1); for event study plots, restricts plotted time points to those where the number of treated units is at least this fraction of the maximum number of treated units observed across all periods. Default is \code{0.3}.}
  \item{pre.periods}{Optional numeric vector; specifies pre-treatment periods to be used for bounding and equivalence tests. Default is \code{NULL}, which uses the pre-treatment periods defined in the \code{fect} object.}
  \item{f.threshold}{Numeric or \code{NULL} (default); F-test threshold for equivalence checks. If \code{NULL}, uses the value from the \code{fect} object.}
  \item{tost.threshold}{Numeric or \code{NULL} (default); threshold for TOST-based equivalence checks (e.g., for pre-trend, placebo, or carryover tests). If \code{NULL}, uses the value from the \code{fect} object.}
  \item{effect.bound.ratio}{Logical; if \code{TRUE} (default is \code{FALSE}), shows the ratio of ATT to the minimal bound in equivalence plots.}
  \item{stats}{Character vector or \code{NULL} (default); specifies which test statistics to display on the plot (e.g., \code{"F.p"}, \code{"equiv.p"}, \code{"placebo.p"}). Default is context-dependent based on plot type and tests performed.}
  \item{stats.labs}{Character vector or \code{NULL} (default); custom labels for the displayed statistics. Must match the length of \code{stats}.}
  \item{raw}{For \code{type = "counterfactual"} plot: \code{"none"} (default; shows average treated and counterfactual lines), \code{"band"} (shows 5-95\% quantile band for raw control/treated data), or \code{"all"} (shows all raw control/treated unit trajectories).}
  \item{vis}{(Deprecated) Formerly controlled line display for placebo/carryover tests. This is now handled by \code{highlight} and internal logic.}
  \item{main}{Character string or \code{NULL} (default); main plot title. If \code{NULL}, a default title is used based on the plot type.}
  \item{xlim}{Numeric vector of length 2 or \code{NULL} (default); specifies x-axis limits.}
  \item{ylim}{Numeric vector of length 2 or \code{NULL} (default); specifies y-axis limits.}
  \item{xlab}{Character string or \code{NULL} (default); x-axis label. If \code{NULL}, a default label is used based on the plot type.}
  \item{ylab}{Character string or \code{NULL} (default); y-axis label. If \code{NULL}, a default label is used based on the plot type.}
  \item{xangle}{Numeric or \code{NULL} (default); angle (in degrees) for x-axis tick labels.}
  \item{yangle}{Numeric or \code{NULL} (default); angle (in degrees) for y-axis tick labels.}
  \item{xbreaks}{Numeric vector or \code{NULL} (default); positions of major ticks on the x-axis.}
  \item{ybreaks}{Numeric vector or \code{NULL} (default); positions of major ticks on the y-axis.}
  \item{xticklabels}{Character vector or \code{NULL} (default); custom labels for the x-axis ticks. Primarily used in status plots.}
  \item{yticklabels}{Character vector or \code{NULL} (default); custom labels for the y-axis ticks. Primarily used in status plots.}
  \item{gridOff}{Logical or \code{NULL} (default). If \code{TRUE}, removes major and minor grid lines. If \code{NULL} (which defaults to \code{FALSE} internally for this argument), the grid is off for most plots but on for \code{type = "status"} (where it forms tile outlines).}
  \item{legendOff}{Logical; if \code{TRUE} (default is \code{FALSE}), hides the legend.}
  \item{legend.pos}{Character string or \code{NULL} (default); position of the legend (e.g., "bottom", "top", "right", "left", "none", or c(x,y) coordinates). Default is typically "bottom".}
  \item{legend.nrow}{Integer or \code{NULL} (default); number of rows in the legend.}
  \item{legend.labs}{Character vector or \code{NULL} (default); custom legend labels.}
  \item{stats.pos}{Numeric vector of length 2 or \code{NULL} (default); (x,y) coordinates for displaying the text of test statistics on the plot.}
  \item{show.stats}{Logical; if \code{TRUE} (default), displays test statistics on the plot if \code{stats} is not "none".}
  \item{theme.bw}{Logical; if \code{TRUE} (default), uses a black-and-white ggplot2 theme (\code{theme_bw()}). If \code{FALSE} and no \code{preset} is specified, uses \code{theme_grey()}.}
  \item{nfactors}{Integer or \code{NULL} (default); number of factors to plot when \code{type} is \code{"factors"} or \code{"loadings"}. Defaults to \code{min(r.cv, 4)}, where \code{r.cv} is the number of estimated factors.}
  \item{include.FE}{Logical; if \code{TRUE} (default), includes fixed effects as additional "factors" in factor/loading plots if applicable (i.e., if unit or time fixed effects were estimated as part of the factor model).}
  \item{id}{Character vector or \code{NULL} (default); when specified, plots only the chosen unit(s). For counterfactual plots (\code{type="counterfactual"}), only one ID is typically supported for individual unit display; otherwise, averages are shown.}
  \item{cex.main}{Numeric or \code{NULL} (default); font size multiplier for the main plot title. Base size is 16pt.}
  \item{cex.main.sub}{Numeric or \code{NULL} (default); font size multiplier for a possible subtitle. Base size is 16pt.}
  \item{cex.axis}{Numeric or \code{NULL} (default); font size multiplier for axis tick labels. Base size is 15pt.}
  \item{cex.lab}{Numeric or \code{NULL} (default); font size multiplier for axis labels. Base size is 15pt.}
  \item{cex.legend}{Numeric or \code{NULL} (default); font size multiplier for legend text. Base size is 15pt.}
  \item{cex.text}{Numeric or \code{NULL} (default); font size multiplier for text annotations on the plot (e.g., statistics values). Base size is 5pt.}
  \item{axis.adjust}{Logical; if \code{TRUE} (default is \code{FALSE}), attempts to adjust axis labels for plots with many time points (e.g., by rotating x-axis labels).}
  \item{axis.lab}{Controls axis labeling in \code{type="status"} plots. Options: \code{"both"} (default), \code{"time"}, \code{"unit"}, or \code{"off"}.}
  \item{axis.lab.gap}{Numeric vector of length 1 or 2; gap spacing for axis labels in \code{type="status"} plots (sets frequency of labels). Default is \code{c(0, 0)} (label every tick).}
  \item{shade.post}{Logical; if \code{FALSE} (default). If \code{TRUE}, shades post-treatment periods in certain plots. For counterfactual plots, if \code{NULL} (the default for the argument is \code{FALSE}), it effectively becomes \code{TRUE} internally unless explicitly set to \code{FALSE}.}
  \item{start0}{Logical; if \code{TRUE} (default is \code{FALSE}), re-indexes the first post-treatment period as 0 in event-time plots (gap, equiv, exit). The period before treatment becomes -1.}
  \item{return.test}{Logical; if \code{TRUE} (default is \code{FALSE}), returns a list containing the plot object and relevant test statistics instead of just the plot.}
  \item{return.data}{Logical; if \code{TRUE} (default is \code{FALSE}), returns a data frame with the data used to construct the plot instead of returning a plot object. If \code{return.test=TRUE} as well, returns a list with \code{data} and \code{test.out}.}
  \item{balance}{Logical or \code{NULL} (default); if \code{TRUE}, uses a "balance sample" (if available from the \code{fect} object, e.g., from setting \code{balanced.sample = TRUE} in \code{fect}) for plotting.}
  \item{weight}{Logical or \code{NULL} (default); if \code{TRUE}, applies weighting (if available from the \code{fect} object, e.g., from providing weights in \code{fect}) to the estimates being plotted.}
  \item{lcolor}{Character vector of length 1 or 2, or \code{NULL} (default). If \code{NULL}, uses a default color (e.g., \code{"#AAAAAA70"} if \code{theme.bw=TRUE}, else \code{"white"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{lwidth}{Numeric vector of length 1 or 2, or \code{NULL} (default). If \code{NULL}, uses a default width (e.g., \code{1.5} if \code{theme.bw=TRUE}, else \code{2} for the base style). The \code{preset} argument may define a different width. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{ltype}{Character vector of length 1 or 2, or \code{NULL} (default). If \code{NULL}, uses a default linetype (e.g., \code{c("solid", "solid")} for the base style). The \code{preset} argument may define a different linetype (e.g., \code{preset="vibrant"} uses \code{c("solid", "dashed")}). Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{line.color}{(Deprecated) Use \code{color} for main estimate lines or \code{lcolor} for reference lines instead.}
  \item{line.width}{(Deprecated) Use \code{est.lwidth} for main estimate lines or \code{lwidth} for reference lines instead.}
  \item{count}{(Deprecated) Use \code{show.count} to control visibility and \code{count.color}, \code{count.alpha}, \code{count.outline.color} for appearance.}
  \item{preset}{A character string specifying a color and style preset: \code{NULL} (default, uses a base style, often with \code{theme_bw} if \code{theme.bw=TRUE}), \code{"vibrant"}, or \code{"grayscale"}. This sets a group of defaults for many plot elements. Individual appearance arguments (like \code{color}, \code{lcolor}, etc.) can override settings from a \code{preset} or the base style.}
  \item{connected}{Logical or \code{NULL} (default). For event study type plots (e.g., \code{type="gap"}), if \code{TRUE}, connects points with lines; if \code{FALSE}, shows point-ranges. If \code{NULL}, defaults to \code{FALSE} for the base style, and \code{TRUE} if \code{preset = "vibrant"}.}
  \item{ci.outline}{Logical; if \code{FALSE} (default). For event study plots where CIs are shown as ribbons, if \code{TRUE}, adds an outline to the CI ribbon.}
  \item{color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for main plot lines/points (e.g., \code{"black"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{est.lwidth}{Numeric or \code{NULL}; line width for the main estimate line in event study plots. If \code{NULL} (default), the width is determined based on \code{connected} and \code{show.points} (e.g., `0.6` for point-ranges with the base style, `0.7` for connected points with markers with the base style). The \code{preset} argument may influence this default. Specifying a numeric value directly overrides any other setting.}
  \item{est.pointsize}{Numeric or \code{NULL}; size of points for the main estimate in event study plots. If \code{NULL} (default), the size is determined based on \code{connected} and \code{show.points} (e.g., `2` for point-ranges with the base style, `1.2` for connected points with markers with the base style). The \code{preset} argument may influence this default. Specifying a numeric value directly overrides any other setting.}
  \item{count.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for count bars (e.g., \code{"grey70"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{count.alpha}{Numeric (0-1) or \code{NULL} (default). If \code{NULL}, uses a default alpha for count bars (e.g., \code{0.4} for the base style). The \code{preset} argument may define a different alpha. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{count.outline.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default outline color for count bars (e.g., \code{"grey69"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{placebo.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for placebo elements (e.g., \code{"blue"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{carryover.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for carryover elements (e.g., \code{"red"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{carryover.rm.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for removed carryover elements (e.g., \code{"blue"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{sens.original.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for original estimates in sensitivity plots (e.g., \code{"darkblue"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{sens.colors}{Vector of colors or \code{NULL} (default). If \code{NULL}, uses a default palette for sensitivity bounds (e.g., \code{c("#218C23","#FF34B4","#FF521B","#2B59C3")} for the base style). The \code{preset} argument may define a different palette. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{counterfactual.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for counterfactual lines (e.g., \code{"steelblue"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{counterfactual.raw.controls.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for raw control lines (e.g., \code{"#4682B420"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{counterfactual.raw.treated.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for raw treated lines (e.g., \code{"#77777750"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{counterfactual.linetype}{Character string or \code{NULL} (default). If \code{NULL}, uses a default linetype for counterfactual lines (e.g., \code{"longdash"} for the base style). The \code{preset} argument may define a different linetype. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{box.control}{Character string or \code{NULL} (default). If \code{NULL}, uses a default fill color for control boxplots (e.g., \code{"skyblue"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{box.treat}{Character string or \code{NULL} (default). If \code{NULL}, uses a default fill color for treated boxplots (e.g., \code{"pink"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{calendar.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for fitted lines in calendar plots (e.g., \code{"skyblue"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{calendar.lcolor}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for average ATT lines in calendar plots (e.g., \code{"red"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{calendar.cicolor}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for CI in calendar plots (e.g., \code{"red"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{heterogeneous.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for fitted lines in heterogeneous plots (e.g., \code{"red"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{heterogeneous.lcolor}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for average ATT lines in heterogeneous plots (e.g., \code{"red"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{heterogeneous.cicolor}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for CI in heterogeneous plots (e.g., \code{"red"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{equiv.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for equivalence bounds (e.g., \code{"red"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{status.treat.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for treated status (e.g., \code{"#06266F"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{status.control.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for control status (e.g., \code{"#B0C4DE"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{status.missing.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for missing status (e.g., \code{"#FFFFFF"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{status.removed.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for removed status (e.g., \code{"#A9A9A9"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{status.placebo.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for placebo status (e.g., \code{"#66C2A5"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{status.carryover.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for carryover status (e.g., \code{"#E78AC3"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{status.carryover.rm.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for removed carryover status (e.g., \code{"#ffc425"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{status.balanced.post.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for balanced post-treatment status (e.g., \code{"#00852B"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{status.balanced.pre.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default color for balanced pre-treatment status (e.g., \code{"#A5CA18"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{status.background.color}{Character string or \code{NULL} (default). If \code{NULL}, uses a default background color for status plots (e.g., \code{"gray90"} for the base style). The \code{preset} argument may define a different color. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{covariate}{Character string or \code{NULL} (default). If \code{NULL}, uses a default covariate for heterogeneous plots (e.g., \code{"X"} for the base style). The \code{preset} argument may define a different covariate. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{covariate.labels}{Character vector or \code{NULL} (default). If \code{NULL}, uses a default covariate labels for heterogeneous plots (e.g., \code{"X1"} for the base style). The \code{preset} argument may define a different covariate labels. Specifying this argument directly overrides any \code{preset} or base style setting.}
  \item{...}{Additional graphical parameters passed to internal plotting routines, primarily those accepted by \code{esplot} for event study style plots (gap, equiv, exit, sens_es, cumul).}
}
\details{
  \code{plot.fect} generates various visualizations for objects estimated via \code{\link{fect}}. Depending on
  the selected \code{type}, it can show treatment-effect dynamics (\code{"gap"}, \code{"exit"}, \code{"cumul"}),
  equivalence test outcomes (\code{"equiv"}), factor/loadings for factor models (\code{"factors"}, \code{"loadings"}),
  raw vs. counterfactual trajectories (\code{"counterfactual"}), treatment status (\code{"status"}),
  effects by calendar time (\code{"calendar"}), distribution of individual effects (\code{"box"}),
  or sensitivity analysis results (\code{"sens"}, \code{"sens_es"}).
  The \code{preset} argument allows for quick customization of plot aesthetics using predefined color schemes.
  Many individual color and style parameters can be set to further customize the appearance, and these will override any settings from a \code{preset} or the base style.
  The function heavily relies on \code{ggplot2} for plotting and \code{esplot} for event-study style visualizations.
}
\value{
  By default, a \code{ggplot2} object representing the plot.

  If \code{return.data=TRUE}, returns a list containing:
  \item{p}{The \code{ggplot2} object.}
  \item{data}{A named list of data frames containing the data used to construct each layer of the plot (no extra all-NA columns are added).}
  If \code{return.test=TRUE} as well, the list also contains:
  \item{test.out}{A list or data frame containing relevant test statistics (such as F-tests, equivalence p-values, placebo test p-values, etc.), if applicable to the plot type and options chosen.}

  If \code{return.test=TRUE} and \code{return.data=FALSE}, returns a list containing \code{p} and \code{test.out}.
}
\author{
  Licheng Liu, Ye Wang, Yiqing Xu, Ziyi Liu, Rivka Lipkovitz
}
\references{
  Xu, Y. (2017). Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models.
  \emph{Political Analysis, 25}(1), 57–76.

  Liu, L., Wang, Y., & Xu, Y. (2022). A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional Data.
  \emph{American Journal of Political Science, 66}(1), 220-237. (Provides context for \code{esplot} which is used internally for many plot types)

  Rambachan, A., & Roth, J. (2023). A More Credible Approach to Parallel Trends.
  \emph{Review of Economic Studies, 90}(5), 2555-2591. (Provides context for sensitivity analysis plots like \code{type="sens"} and \code{type="sens_es"})
}
\examples{
library(fect)
# For CRAN checks, use a small number of bootstraps
# In practice, use a larger number (e.g., nboots = 200 or more)
if(requireNamespace("ggplot2") && requireNamespace("ggrepel")) {
  data(simdata)
  # Estimate with fixed effects method
  out.fect <- fect(
    Y ~ D + X1 + X2,
    data = simdata,
    index = c("id","time"),
    method = "fe",
    force = "two-way",
    se = TRUE,
    parallel = FALSE,
    nboots = 5 # nboots low for example
  )

  # Default gap plot
  plot(out.fect, main = "Estimated ATT (FEct)", ylab = "Effect of D on Y")

  # Gap plot with vibrant preset and custom line color
  # plot(out.fect, preset = "vibrant", color = "darkgreen",
  #      main = "Estimated ATT (Vibrant Preset, Custom Line)")

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

  # Example for sensitivity plot (requires IFE/GSYNTH method and sensitivity analysis)
  # \donttest{
  #   out.ife <- fect(Y ~ D, data = simdata, index = c("id","time"),
  #                   method = "ife", se = TRUE, r = 2,
  #                   sensitivity.analysis = TRUE, sensitivity.plot = FALSE, # run analysis
  #                   parallel = FALSE, nboots = 5) # nboots low for example
  #   if (!is.null(out.ife$sensitivity.rm)) {
  #     plot(out.ife, type = "sens", restrict = "rm",
  #          main = "Sensitivity Analysis (Relative Magnitude)")
  #     plot(out.ife, type = "sens_es", restrict = "rm",
  #          main = "Event-Study Sensitivity (Relative Magnitude)")
  #   }
  # }
}
}
