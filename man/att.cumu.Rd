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
  \item{weighted}{A logical flag that chooses how the per-period ATTs are combined. The per-period ATT at event time \eqn{j}, \eqn{\mathrm{ATT}_j}{ATT_j}, is the average effect over the treated cells at that event time (the values in \code{x$att}, or \code{x$att.off} with \code{type = "off"}), and \eqn{n_j} is the number of those cells. Take a window from event time \eqn{s} to event time \eqn{k}, with \eqn{L = k - s + 1} event times. \code{FALSE} (the default) returns the running sum \eqn{\sum_{j=s}^{k} \mathrm{ATT}_j}{ATT_s + ... + ATT_k}. \code{TRUE} returns the count-weighted version \eqn{L \sum_{j=s}^{k} n_j \mathrm{ATT}_j / \sum_{j=s}^{k} n_j}{L * sum(n_j * ATT_j) / sum(n_j)}, that is, \eqn{L} times the average effect over all treated cells in the window. Before version 2.4.7, \code{weighted} had no effect and \code{att.cumu()} always returned the count-weighted version. The two agree when every event time in the window has the same number of treated cells.}
  \item{alpha}{A numerical value specifying the significance level. Default is \code{0.05}.}
  \item{type}{A string that specifies the type of effect to calculate. Must be one of the following: \code{"on"} (switch-on treatment effect) or \code{"off"} (switch-off treatment effect). Default is \code{"on"}.}
  \item{plot}{A logical flag indicating whether to plot cumulative effects. Default is \code{FALSE}.}
}
\value{
  A numeric matrix with one row per accumulation window, from \code{period[1]}
  to each subsequent period up to \code{period[2]}, with columns \code{start},
  \code{end}, and \code{catt}: the cumulative ATT over that window, as
  defined by \code{weighted} (by default, the sum of the per-period ATTs
  from \code{start} to \code{end}). When the \code{fect} object carries
  uncertainty estimates, four further columns give the bootstrap standard
  error (\code{S.E.}), the confidence interval at level \code{alpha}
  (\code{CI.lower}, \code{CI.upper}), and the \code{p.value}; without
  standard errors (\code{se = FALSE}) only the first three columns are
  returned. Each replicate's cumulative ATT is built the same way as
  \code{catt}. With \code{weighted = FALSE}, a replicate that has no treated
  cell at one of the window's event times is left out of that row, as in
  \code{\link{effect}}, so the two functions use the same replicates. The
  interval and p-value depend on the fit's variance type
  (\code{x$vartype}). For bootstrap fits they come from the quantiles of
  the replicate estimates and from the share of replicates on each side
  of zero. For parametric fits they are \code{catt} plus or minus the
  normal critical value times \code{S.E.} and the two-sided normal
  p-value, because parametric replicates are simulated under no treatment
  effect and centred near zero. The first row repeats the one-period
  estimate from \code{x$est.att}. With \code{plot = TRUE} the cumulative
  effects are also plotted.
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