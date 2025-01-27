\encoding{UTF-8}
\name{att.cumu}
\alias{att.cumu}
\title{Calculate Cumulative Treatment Effects}
\description{
  Calculate cumulative treatment effects based on the results of a \code{\link{fect}} object.
}
\usage{
att.cumu(x, period = NULL, weighted = TRUE, alpha = 0.05, type = "on", plot = FALSE)
}
\arguments{
  \item{x}{A \code{\link{fect}} object.}
  \item{period}{A two-element numeric vector specifying the range of terms during which treatment effects are to be accumulated, e.g., \code{period = c(-1, 1)}.}
  \item{weighted}{A logical flag specifying whether to calculate weighted cumulative treatment effects based on counts at each period. Default is \code{TRUE}.}
  \item{alpha}{A numerical value specifying the significance level. Default is \code{0.05}.}
  \item{type}{A string that specifies the type of effect to calculate. Must be one of the following: \code{"on"} (switch-on treatment effect) or \code{"off"} (switch-off treatment effect). Default is \code{"on"}.}
  \item{plot}{A logical flag indicating whether to plot cumulative effects. Default is \code{FALSE}.}
}
\author{
  Licheng Liu, Ye Wang, and Yiqing Xu
}
\references{
  Athey, S., Bayati, M., Doudchenko, N., Imbens, G., and Khosravi, K. (2021).
  Matrix completion methods for causal panel data models.
  \emph{Journal of the American Statistical Association}, 116(536), 1716-1730.

  Bai, J. (2009).
  Panel data models with interactive fixed effects.
  \emph{Econometrica}, 77(4), 1229-1279.

  Liu, L., Wang, Y., and Xu, Y. (2022).
  A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional Data.
  \emph{American Journal of Political Science}, 68(1), 160-176.

  Xu, Y. (2017).
  Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models.
  \emph{Political Analysis}, 25(1), 57-76.
}
\seealso{
  \code{\link{fect}}, \code{\link{plot.fect}}
}