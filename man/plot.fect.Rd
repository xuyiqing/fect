\name{plot.fect}
\alias{plot.fect}
\title{Plotting}
\description{Visualizes estimation results of the matrix completion method.}
\usage{\method{plot}{fect}(x, type = "gap", switch.on = TRUE,
            bound = NULL, vis = "connected",
            count = TRUE, proportion = 0.3,
            equiv.bound = NULL, equiv.bound.size = 1,
            effect.bound.ratio = FALSE, p.value = TRUE, main = NULL, 
            xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, 
            legendOff = FALSE, legend.pos = NULL, legend.nrow = NULL,
            legend.labs = NULL, text.pos = NULL, theme.bw = TRUE,
            gridOff = FALSE, id = NULL,
            cex.main = NULL, cex.main.sub = NULL, cex.axis = NULL,
            cex.lab = NULL, cex.legend = NULL, cex.text = NULL,
            axis.adjust = FALSE, axis.lab = "both", axis.lab.gap = c(0, 0), \dots)
}
\arguments{
  \item{x}{a \code{\link{fect}} object.}
  \item{type}{a string that specifies the type of the plot. Must be one
    of the following: "gap" (plotting the treatment effect); "raw" (plotting the raw data); "missing" (plotting status of each unit at each time point).}
  \item{switch.on}{a logical flag that specifies the type of the gap plot. Must be one
    of the following: \code{TRUE} (plotting the switch-on treatment effect); 
    \code{FALSE} (plotting the switch-off treatment effect). default is 
    \code{switch.on = TRUE}.}
  \item{bound}{a string that specifies the bounds to be plotted for equivalence test 
    for pre-treatment fit checking. Choose from: "both", "equiv", "min" and "none".}
  \item{vis}{A string specifying whether to plot the dots for placebo plots.}
  \item{count}{a logical flag controlling whether to show the count of each term 
    for gap plot.}
  \item{proportion}{a positive value specifying periods at which observations equal to 
    or greater than the proporation of the largest number of observations at a certain 
    period.}
  \item{equiv.bound}{User defined equivalence margin.}
  \item{equiv.bound.size}{a scale parameter for estimated ATT that is used as 
  the equivalence margin.}
  \item{effect.bound.ratio}{a logical value specifiying whether to annotate 
  the ratio of estimated average treatment effects / minimun bound.}
  \item{p.value}{a logical value specifying whether to label p value for placebo test.}
  \item{main}{a string that controls the title of the plot. If not
    supplied, no title will be shown.}
  \item{xlim}{a two-element numeric vector specifying the range of x-axis. When 
    class of time varible is string, must specify not original value but a counting number e.g. \code{xlim=c(1,30)}.}
  \item{ylim}{a two-element numeric vector specifying the range of y-axis.}
  \item{xlab}{a string indicating the label of the x-axis.}
  \item{ylab}{a string indicating the label of the y-axis.}
  \item{legendOff}{a logical flag controlling whether to show the legend.}
  \item{legend.pos}{a string specifying the position of legend. If left blank, 
    legend will be setted at the bottom.}
  \item{legend.nrow}{an integer speficying rows of legend.}
  \item{legend.labs}{a string vector for user-defined legends.}
  \item{text.pos}{a numeric vector of length 2 specifying the postion for 
    labels of equivalence bounds.}
  \item{theme.bw}{a logical flag specifying whether to use the black and 
    white theme.}
  \item{gridOff}{a logical flag specifying whether to move grid lines.}
  \item{id}{a string vector specifying a sub-group of units that are to be 
    plotted for treatment status(\code{type = "missing"}).}
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
  \item{\dots}{other argv.}
}
\details{
  \code{plot.fect} visualizes the estimation results obtained by \code{fect}.
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



