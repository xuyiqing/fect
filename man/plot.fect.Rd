\name{plot.fect}
\alias{plot.fect}
\title{Plotting}
\description{Visualizes estimation results of the matrix completion method.}
\usage{\method{plot}{fect}(x, type = "raw", gap.type = "on", 
            bound = "both", vis = "connected",
            count = TRUE, proportion = NULL,
            equiv.bound = NULL, equiv.bound.size = 1,
            effect.bound.ratio = FALSE,
            xlim = NULL, ylim = NULL,
            xlab = NULL, ylab = NULL, 
            label.size = 5, label.pos = NULL, legend.pos = NULL,
            theme.bw = FALSE,
            legendOff = FALSE, gridOff = FALSE, main = NULL, id = NULL, 
            cex.main = 15, cex.main.sub = 12, cex.axis = 8,
            cex.lab = 12, cex.legend = 12,
            axis.adjust = FALSE, axis.lab = "both", axis.lab.gap = c(0, 0), \dots)
}
\arguments{
  \item{x}{a \code{\link{fect}} object.}
  \item{type}{a string that specifies the type of the plot. Must be one
    of the following: "gap" (plotting the treatment effect); "raw" (plotting the raw data); "missing" (plotting status of each unit at each time point).}
  \item{gap.type}{a string that specifies the type of the gap plot. Must be one
    of the following: "on" (plotting the switch-on treatment effect); "off" (plotting the switch-off treatment effect). default is \code{gap.type = "on"}.}
  \item{bound}{a string that specifies the bounds to be plotted for equivalence test 
    for pre-treatment fit checking. Choose from: "both", "equiv", "min" and "none".}
  \item{vis}{For placebo plot. Whether to plot the dots.}
  \item{count}{a logical flag controlling whether to show the count of each term 
    for gap plot.}
  \item{proportion}{a positive value ...}
  \item{equiv.bound}{User defined equivalence margin.}
  \item{equiv.bound.size}{a scale parameter for estimated ATT that is used as 
  the equivalence margin.}
  \item{effect.bound.ratio}{a logical value specifiying whether to annotate 
  the ratio of estimated average treatment effects / minimun bound.}
  \item{xlim}{a two-element numeric vector specifying the range of x-axis. When 
    class of time varible is string, must specify not original value but a counting number e.g. \code{xlim=c(1,30)}.}
  \item{ylim}{a two-element numeric vector specifying the range of y-axis.}
  \item{xlab}{a string indicating the label of the x-axis.}
  \item{ylab}{a string indicating the label of the y-axis.}
  \item{label.size}{...}
  \item{label.pos}{...}
  \item{legend.pos}{...}
  \item{theme.bw}{a logical flag indicating whether black and white theme will be used.}
  \item{legendOff}{a logical flag controlling whether to show the legend.}
  \item{gridOff}{...}
  \item{main}{a string that controls the title of the plot. If not
    supplied, no title will be shown.}
  \item{id}{a vector specifying units to be plotted.}
  \item{cex.main}{a numeric value (pt) specifying the fontsize of the main title.}
  \item{cex.main.sub}{a numeric value (pt) specifying the fontsize of the subtitles.}
  \item{cex.axis}{a numeric value (pt) specifying the fontsize of the texts on the axes.}
  \item{cex.lab}{a numeric value (pt) specifying the fontsize of the axis titles.}
  \item{cex.legend}{a numeric value (pt) specifying the fontsize of the legend.}
  \item{axis.adjust}{a logic flag indicating whether to adjust labels on x-axis. 
    Useful when class of time variable is string and data magnitude is large.}
  \item{axis.lab}{a string indicating whether labels on the x- and y-axis will be shown. There are four options: \code{"both"} (default): labels on both axes will be shown; \code{"unit"}: only labels on y-axis will be shown; \code{"time"}: only labels on the x-axis will be shown; "none": no labels will be shown.}
  \item{axis.lab.gap}{a numeric vector setting the gaps between labels on the x- or y-axis for \code{"missing"} plot. Default is \code{axis.lab.gap = c(0, 0)}, which means that all labels will be shown. Useful for datasets with large N or T.}
  \item{\dots}{other argv.}
}
\details{
  \code{plot.fect} visualizes the raw data used by the matrix completion  method.
}
\author{
  Licheng Liu; Ye Wang; Yiqing Xu
}
\references{
  Jushan Bai. 2009. "Panel Data Models with Interactive Fixed
  Effects." Econometrica 77:1229--1279.

  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference
  with Interactive Fixed Effects Models." Political Analysis, Vol. 25, 
  Iss. 1, January 2017, pp. 57-76. Available at: \url{https://doi.org/10.1017/pan.2016.2}.

  Athey S, Bayati M, Doudchenko N, et al. Matrix completion methods for causal panel data models[J]. arXiv preprint arXiv:1710.10251, 2017.

  For more details about the matrix completion method, see \url{https://github.com/susanathey/MCPanel}. 
}
\seealso{
  \code{\link{fect}} and \code{\link{print.fect}}
}



