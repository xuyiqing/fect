\name{att.cumu}
\alias{att.cumu}
\title{Calculate Cumulative Treatment Effects}
\description{Calculate cumulative treatment effects}
\usage{att.cumu(x, period = NULL, weighted = TRUE, alpha = 0.05, type = "on", plot = FALSE)} 
\arguments{
  \item{x}{a \code{\link{fect}} object.}
  \item{period}{a two-element numeric vector specifying the range of term during which treatment effects are to be accumulated. 
  e.g. \code{period = c(-1,1)}.}
  \item{weighted}{a logical flag specifying whether to calculate weigthed cumulative treatment effects based on counts at each period. Default is 
  \code{weighted = TRUE}.}
  \item{alpha}{a numerical value that specfies significant level.}
  \item{type}{a string that specifies the type. Must be one of the following: "on" (switch-on treatment effect); "off" (switch-off treatment effect). Default 
  is \code{type = "on"}.}
  \item{plot}{A logical flag indicating whether to plot cumulative effects. 
  Default is \code{plot = FALSE}.}
}
\author{
  Licheng Liu; Ye Wang; Yiqing Xu 
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
  \code{\link{fect}} and \code{\link{plot.fect}}
}


