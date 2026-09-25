\encoding{UTF-8}
\name{att.cumu}
\alias{att.cumu}
\title{Calculate Cumulative Treatment Effects}
\description{
  Calculate cumulative treatment effects based on the results of a \code{\link{fect}} object.
}
\usage{
att.cumu(x, period = NULL, weighted = FALSE, alpha = 0.05, type = "on", plot = FALSE)
}
\arguments{
  \item{x}{A \code{\link{fect}} object.}
  \item{period}{A two-element numeric vector specifying the range of terms during which treatment effects are to be accumulated, e.g., \code{period = c(-1, 1)}. Default \code{NULL}: from event time 1 to the last event time.}
  \item{weighted}{A logical flag choosing how the per-period ATTs are combined. Each per-period ATT, \eqn{ATT_j}, is the average effect over the treated cells at event time \eqn{j} (the values in \code{x$att}). \code{FALSE} (the default) returns their running sum: the cumulative ATT from event time \eqn{s} to \eqn{k} is \eqn{\sum_{j=s}^{k} ATT_j}{ATT_s + ... + ATT_k}. \code{TRUE} returns the count-weighted version: \eqn{L} times the average effect over all treated cells in the \eqn{L} event times from \eqn{s} to \eqn{k}, that is \eqn{L \sum_j n_j ATT_j / \sum_j n_j}{L * sum(n_j * ATT_j) / sum(n_j)}, where \eqn{n_j} is the number of treated cells at event time \eqn{j}. This was \code{fect}'s cumulative ATT before version 2.4.6. The two agree when every event time has the same number of treated cells.}
  \item{alpha}{A numerical value specifying the significance level. Default is \code{0.05}.}
  \item{type}{A string that specifies the type of effect to calculate. Must be one of the following: \code{"on"} (switch-on treatment effect) or \code{"off"} (switch-off treatment effect). Default is \code{"on"}.}
  \item{plot}{A logical flag indicating whether to plot cumulative effects. Default is \code{FALSE}.}
}
\value{
  A numeric matrix with one row per accumulation window, from \code{period[1]}
  to each subsequent period up to \code{period[2]}, with columns \code{start},
  \code{end}, and \code{catt}: the cumulative ATT over that window, defined by
  \code{weighted} (by default, the sum of the per-period ATTs from \code{start}
  to \code{end}). The first row is the ATT of the first period. When the fit
  has standard errors (\code{se = TRUE}), four further columns give the
  bootstrap standard error (\code{S.E.}), the confidence interval at level
  \code{alpha} (\code{CI.lower}, \code{CI.upper}), and the \code{p.value};
  without standard errors only the first three columns are returned. With
  \code{plot = TRUE} the cumulative effects are also plotted.
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