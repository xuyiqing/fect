\name{att.cumu}
\alias{att.cumu}
\title{Calculate Cumulative Treatment Effects}
\description{Calculate cumulative treatment effects}
\usage{att.cumu(x, period = NULL, weighted = TRUE, alpha = 0.05, type = "on", plot = FALSE)} 
\arguments{
  \item{x}{a \code{\link{fect}} object.}
  \item{period}{a two-element numeric vector specifying the range of term during which treatment effects are to be accumulated. 
  e.g. \code{time.on.lim = c(-1,1)}.}
  \item{weighted}{a logical flag specifying whether to calculate weigthed cumulative treatment effects based on counts at each period. Default is 
  \code{weighted = TRUE}.}
  \item{alpha}{significant levels.}
  \item{type}{a string that specifies the type. Must be one of the following: "on" (switch-on treatment effect); "off" (switch-off treatment effect). Default 
  is \code{type = "on"}.}
  \item{plot}{A logical flag indicating whether to plot cumulative effects.}
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
  \code{\link{fect}} and \code{\link{plot.fect}}
}


