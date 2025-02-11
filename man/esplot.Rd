
\name{esplot}
\alias{esplot}
\title{Event Study Visualization}
\description{Visualize the dynamic treatment effects and create an event study plot.}
\usage{
    esplot(data, Period, Estimate, SE, CI.lower = NULL,
        CI.upper = NULL, Count = NULL, fill.gap = TRUE,
        start0 = FALSE, show.count = NULL, stats = NULL,
        stats.labs = NULL, highlight.periods = NULL,
        highlight.colors = NULL, main = NULL,
        xlim = NULL, ylim = NULL,
        xlab = NULL, ylab = NULL,
        gridOff = FALSE, stats.pos = NULL,
        theme.bw = TRUE, cex.main = NULL,
        cex.axis = NULL, cex.lab = NULL,
        cex.text = NULL, axis.adjust = FALSE)
}
\arguments{
    \item{data}{The input data for the event study plot.}
    \item{Period}{The variable name of the relative period.}
    \item{Estimate}{The variable name of the estimated coefficient.}
    \item{SE}{The variable name of the standard errors. Omitted if \code{CI.lower} and \code{CI.upper} are specified.}
    \item{CI.lower}{Optional. The variable name of the lower confidence interval.}
    \item{CI.upper}{Optional. The variable name of the upper confidence interval.}
    \item{Count}{Optional. The variable name indicating the number of treated observations for each time period.}
    \item{fill.gap}{Logical. Whether to use 0 values to fill gaps in dynamic treatment effects. Default is \code{TRUE}.}
    \item{start0}{Logical. Whether period 0 should be considered as the first post-treatment period. Default is \code{FALSE}.}
    \item{show.count}{Logical. Whether to display the bar plot of the number of observations at the bottom of the plot.}
    \item{stats}{Optional. A list of statistics to be printed on the plot.}
    \item{stats.labs}{Optional. The names corresponding to \code{stats}.}
    \item{highlight.periods}{Optional. A list of periods to highlight.}
    \item{highlight.colors}{Optional. A list of colors for highlighting periods.}
    \item{main}{Optional. The title of the plot.}
    \item{xlim}{Optional. The limits of the x-axis.}
    \item{ylim}{Optional. The limits of the y-axis.}
    \item{xlab}{Optional. The label for the x-axis.}
    \item{ylab}{Optional. The label for the y-axis.}
    \item{gridOff}{Logical. Whether to turn off grid lines. Default is \code{FALSE}.}
    \item{stats.pos}{Optional. The position of the printed statistics.}
    \item{theme.bw}{Logical. Whether to use a black-and-white theme. Default is \code{TRUE}.}
    \item{cex.main}{Optional. The font size for the plot title.}
    \item{cex.axis}{Optional. The font size for the axis labels.}
    \item{cex.lab}{Optional. The font size for the axis tick labels.}
    \item{cex.text}{Optional. The font size for other text elements in the plot.}
    \item{axis.adjust}{Logical. Whether to adjust the axis position. Default is \code{FALSE}.}
}
\value{
\item{p}{A ggplot object of the event study plot.}
}
\author{
Licheng Liu, Yiqing Xu, Ziyi Liu, Zhongyu Yin
}

