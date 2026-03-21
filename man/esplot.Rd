\name{esplot}
\alias{esplot}
\title{Event Study Visualization}
\description{
  Visualize dynamic treatment effects and create an event study plot.
  This function offers flexibility in displaying estimates, confidence intervals,
  and various annotations. It can handle data directly or from `did_wrapper` objects,
  calculate confidence intervals from standard errors if needed, and allows for
  connected (line/ribbon) or point-range style plots.
}
\usage{
esplot(data, Period = NULL, Estimate = "ATT", SE = NULL,
       CI.lower = "CI.lower", CI.upper = "CI.upper", Count = NULL,
       proportion = 0.3, est.lwidth = NULL, est.pointsize = NULL,
       show.points = FALSE, fill.gap = TRUE, start0 = FALSE,
       only.pre = FALSE, only.post = FALSE, show.count = NULL,
       stats = NULL, stats.labs = NULL, highlight.periods = NULL,
       highlight.colors = NULL, lcolor = NULL, lwidth = NULL,
       ltype = c("solid", "solid"), connected = FALSE, ci.outline = FALSE,
       main = NULL, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
       gridOff = FALSE, stats.pos = NULL, theme.bw = TRUE,
       cex.main = NULL, cex.axis = NULL, cex.lab = NULL,
       cex.text = NULL, axis.adjust = FALSE, color = "#000000",
       count.color = "gray70", count.alpha = 0.4,
       count.outline.color = "grey69")
}
\arguments{
  \item{data}{The input data for the event study plot. Can be a \code{data.frame} or an object of class \code{did_wrapper} (in which case, the \code{est.att} component will be used).}
  \item{Period}{The name of the column in \code{data} representing the relative time period. If \code{NULL} (the default), the function will attempt to automatically identify the period column from common names (e.g., 'time', 'Time', 'period', 'Period', 'event.time', 'event_time', 'rel_time') or use numeric rownames if available.}
  \item{Estimate}{The name of the column in \code{data} containing the point estimates (e.g., Average Treatment Effect on the Treated). Default is \code{"ATT"}.}
  \item{SE}{The name of the column in \code{data} containing the standard errors. If columns for \code{CI.lower} and \code{CI.upper} (as specified by their respective arguments, defaulting to \code{"CI.lower"} and \code{"CI.upper"}) are not found in \code{data}, and an \code{SE} column is provided and found, 95\% confidence intervals will be calculated using \code{Estimate +/- 1.96 * SE}. Default is \code{NULL}.}
  \item{CI.lower}{The name of the column in \code{data} for the lower bound of the confidence interval. Default is \code{"CI.lower"}. If this column is not found, it may be calculated from \code{Estimate} and \code{SE} (if \code{SE} is provided and found in \code{data}).}
  \item{CI.upper}{The name of the column in \code{data} for the upper bound of the confidence interval. Default is \code{"CI.upper"}. If this column is not found, it may be calculated from \code{Estimate} and \code{SE} (if \code{SE} is provided and found in \code{data}).}
  \item{Count}{Optional. The name of the column in \code{data} indicating a count measure (e.g., number of observations) for each time period. Used if \code{show.count = TRUE} or for \code{xlim} determination based on \code{proportion}. Default is \code{NULL}.}
  \item{proportion}{Numeric, between 0 and 1. If \code{Count} is specified and \code{xlim} is not, this proportion is used to determine the default x-axis limits. Periods where the count is below \code{proportion * max(Count)} might be excluded from the default view. Default is \code{0.3}.}
  \item{est.lwidth}{Numeric. The line width for the estimate line (if \code{connected = TRUE}) or the main vertical line of the point-range (if \code{connected = FALSE}). Default is \code{NULL}, which resolves to \code{0.6} if \code{connected = FALSE} (the default for \code{connected}), \code{1.2} if \code{connected = TRUE} and \code{show.points = FALSE} (the default for \code{show.points}), and \code{0.7} if \code{connected = TRUE} and \code{show.points = TRUE}.}
  \item{est.pointsize}{Numeric. The size of the points. If \code{connected = TRUE} and \code{show.points = TRUE}, this is the size of points at integer time periods. If \code{connected = FALSE}, this controls the size of the central point in the point-range (via \code{fatten} aesthetic). Default is \code{NULL}, which resolves to \code{2} if \code{connected = FALSE} (the default for \code{connected}), \code{3} if \code{connected = TRUE} and \code{show.points = FALSE} (the default for \code{show.points}), and \code{1.2} if \code{connected = TRUE} and \code{show.points = TRUE}.}
  \item{show.points}{Logical. If \code{connected = TRUE}, whether to display points at integer time periods on top of the line and ribbon. Default is \code{FALSE}.}
  \item{fill.gap}{Logical. If \code{connected = FALSE}, whether to fill gaps in the sequence of time periods with an estimate and confidence interval of 0. This is useful when some integer time periods are missing from the input data. Default is \code{TRUE}.}
  \item{start0}{Logical. If \code{TRUE}, the vertical line separating pre- and post-treatment periods is drawn at \code{x = -0.5}, implying period 0 is the first post-treatment period. If \code{FALSE} (default), the line is at \code{x = 0.5}, implying period 0 is the last pre-treatment period.}
  \item{only.pre}{Logical. If \code{TRUE}, the plot will only display pre-treatment periods. The vertical separator line will be omitted. Default is \code{FALSE}.}
  \item{only.post}{Logical. If \code{TRUE}, the plot will only display post-treatment periods. The vertical separator line will be omitted. Default is \code{FALSE}.}
  \item{show.count}{Logical or \code{NULL}. Whether to display a bar plot of the values from the \code{Count} column at the bottom of the main plot. If \code{NULL} (default), it's treated as \code{FALSE}.}
  \item{stats}{Optional. A numeric vector of statistics (e.g., p-values) to be printed on the plot.}
  \item{stats.labs}{Optional. A character vector of labels corresponding to the \code{stats} values. Must be the same length as \code{stats}.}
  \item{highlight.periods}{Optional. A numeric vector of time periods to highlight with different colors. For \code{connected = TRUE}, these define intervals from \code{period - 0.5} to \code{period + 0.5}. For \code{connected = FALSE}, individual points at these periods are highlighted.}
  \item{highlight.colors}{Optional. A character vector of colors corresponding to \code{highlight.periods}. If \code{NULL} and \code{highlight.periods} is provided, default rainbow colors are used. Must be the same length as \code{highlight.periods}.}
  \item{lcolor}{Optional. Color(s) for the reference lines. Can be a single color (applied to both horizontal y=0 line and vertical pre/post separator line) or a vector of two colors (first for horizontal, second for vertical). If \code{NULL}, defaults to \code{"#aaaaaa"} if \code{theme.bw = TRUE}, otherwise \code{"white"}. Default is \code{NULL}.}
  \item{lwidth}{Optional. Line width(s) for the reference lines. Can be a single width or a vector of two widths (similar to \code{lcolor}). If \code{NULL}, defaults to \code{1.5} if \code{theme.bw = TRUE}, otherwise \code{2}. Default is \code{NULL}.}
  \item{ltype}{Optional. Linetype(s) for the reference lines. Can be a single linetype (applied to both horizontal y=0 line and vertical pre/post separator line) or a vector of two linetypes (first for horizontal, second for vertical). Default is \code{c("solid", "solid")}.}
  \item{connected}{Logical. If \code{TRUE}, estimates and confidence intervals are plotted as a connected line with a ribbon. This involves interpolating values between observed time points (at 0.5 steps by default) to create a smoother appearance. If \code{FALSE} (default), \code{geom_pointrange} is used, showing individual estimates and their CIs as points with ranges for each observed time period.}
  \item{ci.outline}{Logical. If \code{connected = TRUE}, whether to draw an outline around the confidence interval ribbon. The outline color is a slightly darker version of the fill color. Default is \code{FALSE}.}
  \item{main}{Optional. The main title for the plot. If \code{NULL} (default), a default title "Estimated Dynamic Treatment Effects" is used. If an empty string \code{""} is provided, no title is displayed.}
  \item{xlim}{Optional. A numeric vector of length 2 specifying the x-axis limits (\code{c(min, max)}). If \code{NULL} (default), limits are determined automatically based on the data range, potentially filtered by \code{proportion} if \code{Count} is used.}
  \item{ylim}{Optional. A numeric vector of length 2 specifying the y-axis limits (\code{c(min, max)}). If \code{NULL} (default), limits are determined automatically to encompass all estimates and confidence intervals, with potential expansion if \code{show.count = TRUE}.}
  \item{xlab}{Optional. The label for the x-axis. If \code{NULL} (default), "Time Relative to Treatment" is used. If an empty string \code{""} is provided, no label is displayed.}
  \item{ylab}{Optional. The label for the y-axis. If \code{NULL} (default), "Effect on Y" is used. If an empty string \code{""} is provided, no label is displayed.}
  \item{gridOff}{Logical. Whether to turn off major and minor grid lines. Default is \code{TRUE}.}
  \item{stats.pos}{Optional. A numeric vector of length 2 (\code{c(x, y)}) specifying the coordinates for the top-left position of the \code{stats} text block. If \code{NULL} (default), the position is automatically determined.}
  \item{theme.bw}{Logical. Whether to use \code{ggplot2::theme_bw()}. Default is \code{TRUE}.}
  \item{cex.main}{Optional. Numeric scaling factor for the plot title font size. The base size used by ggplot is 16. Default is \code{NULL} (uses base size 16).}
  \item{cex.axis}{Optional. Numeric scaling factor for the axis tick mark labels font size. The base size used by ggplot is 15. Default is \code{NULL} (uses base size 15).}
  \item{cex.lab}{Optional. Numeric scaling factor for the axis title (x and y labels) font size. The base size used by ggplot is 15. Default is \code{NULL} (uses base size 15).}
  \item{cex.text}{Optional. Numeric scaling factor for annotated text elements (e.g., \code{stats} text, count label). The base size used by ggplot for annotation is 5. Default is \code{NULL} (uses base size 5).}
  \item{axis.adjust}{Logical. If \code{TRUE}, x-axis tick labels are rotated 45 degrees for better readability with many labels. Default is \code{FALSE}.}
  \item{color}{Character. The primary color for plotting estimates, points, lines, and confidence interval fills/lines (unless overridden by \code{highlight.colors} for specific periods). Default is \code{"#000000"} (black).}
  \item{count.color}{Character. The fill color for the bars if \code{show.count = TRUE}. Default is \code{"gray70"}.}
  \item{count.alpha}{Numeric. Alpha transparency for the count bars if \code{show.count = TRUE}. Default is \code{0.4}.}
  \item{count.outline.color}{Character. The color for the outline of count bars if \code{show.count = TRUE}. Default is \code{"grey69"}.}
}
\value{
  \item{p}{A \code{ggplot} object representing the event study plot.}
}
\author{
  Licheng Liu, Yiqing Xu, Ziyi Liu, Zhongyu Yin, Rivka Lipkovitz
}
\examples{
# Basic example with simulated data
set.seed(123)
event_data <- data.frame(
  time = -5:5,
  ATT = cumsum(rnorm(11, 0, 0.2)) + c(rep(0,5), 0, 0.5, 1, 1.2, 1.5, 1.3),
  SE = runif(11, 0.1, 0.3)
)
event_data$CI.lower <- event_data$ATT - 1.96 * event_data$SE
event_data$CI.upper <- event_data$ATT + 1.96 * event_data$SE
event_data$count <- sample(50:150, 11, replace = TRUE)
event_data$count[event_data$time == -5 | event_data$time == 5] <- 20 # for proportion demo

# Default plot (point-range)
esplot(event_data, Period = "time", Estimate = "ATT",
       CI.lower = "CI.lower", CI.upper = "CI.upper")

# # Connected plot with ribbon
# esplot(event_data, Period = "time", Estimate = "ATT",
#        CI.lower = "CI.lower", CI.upper = "CI.upper",
#        connected = TRUE, show.points = TRUE)

# # Connected plot using SE for CI calculation
# event_data_no_ci <- event_data[, c("time", "ATT", "SE", "count")]
# esplot(event_data_no_ci, Period = "time", Estimate = "ATT", SE = "SE",
#        connected = TRUE, ci.outline = TRUE, color = "blue")

# # Show count bars and stats
# esplot(event_data, Period = "time", Estimate = "ATT",
#        CI.lower = "CI.lower", CI.upper = "CI.upper", Count = "count",
#        show.count = TRUE, stats = c(0.03, 0.12), stats.labs = c("P-val Pre", "P-val Post"),
#        main = "Event Study with Counts and Stats", proportion = 0.2)

# # Highlight specific periods (connected)
# esplot(event_data, Period = "time", Estimate = "ATT", SE = "SE",
#        connected = TRUE, highlight.periods = c(-1, 2),
#        highlight.colors = c("orange", "green"),
#        main = "Highlighted Periods (Connected)")

# # Highlight specific periods (point-range)
# esplot(event_data, Period = "time", Estimate = "ATT", SE = "SE",
#        connected = FALSE, highlight.periods = c(-1, 2),
#        highlight.colors = c("orange", "green"),
#        main = "Highlighted Periods (Point-Range)")

# # Only post-treatment period, custom labels
# esplot(event_data, Period = "time", Estimate = "ATT", SE = "SE",
#        only.post = TRUE, xlab = "Years Post-Intervention", ylab = "Impact Metric",
#        start0 = TRUE, color = "darkred", est.lwidth = 1.5)

# Using did_wrapper object (conceptual example, requires `did` package and setup)
# if (requireNamespace("did", quietly = TRUE)) {
#   # Assume `did_out` is an output from `did::att_gt` or similar
#   # and `did_wrapper_obj` is created, e.g.,
#   # did_wrapper_obj <- list(est.att = event_data) # Simplified for example
#   # class(did_wrapper_obj) <- "did_wrapper"
#   # esplot(did_wrapper_obj) # Would use defaults: Period="time", Estimate="ATT"
# }

}
