\name{fect_sens}
\alias{fect_sens}
\title{Sensitivity Analysis for fect Objects under Relative Magnitude and Smoothness Restrictions}
\description{
Conducts sensitivity analyses on \code{fect} model objects under relative magnitude (RM) and smoothness (C-LF) assumptions, producing robust confidence intervals and parameter sets.
}
\usage{
fect_sens(
  fect.out,
  post.periods = NA,
  l_vec = NA,
  Mbarvec = seq(0, 1, by = 0.1),
  Mvec = seq(0, 0.25, 0.05),
  periodMbarvec = c(0, 0.5),
  periodMvec = c(0, 0.1),
  parallel = FALSE,
  cores = NULL
)
}
\arguments{
  \item{fect.out}{A fitted \code{fect} object with uncertainty estimates (\code{se = TRUE}) and a placebo test (\code{placeboTest = TRUE} with a range \code{placebo.period = c(a, b)}): the placebo periods are the pre-treatment estimates against which violations are measured.}
  \item{post.periods}{Vector of post-treatment event times for the analysis (1 is the first treated period). Default \code{NA}: every event time from 1 to the last one in \code{fect.out$est.att}.}
  \item{l_vec}{Optional weighting vector, one weight per element of \code{post.periods}, for averaging the ATT across post-treatment periods. Default weights by the number of treated units in each period.}
  \item{Mbarvec}{Values of \code{Mbar} for the relative-magnitude analysis of the average effect. \code{Mbar} is relative: post-treatment violations may be up to \code{Mbar} times the largest pre-treatment violation. \code{NULL} skips this part.}
  \item{Mvec}{Values of \code{M} for the smoothness analysis of the average effect, in the units of the outcome: how much the slope of the violation may change between consecutive periods. \code{NULL} skips this part.}
  \item{periodMbarvec}{Values of \code{Mbar} for period-specific RM sensitivity analysis.}
  \item{periodMvec}{Values of \code{M} for period-specific smoothness sensitivity analysis.}
  \item{parallel}{Logical; if \code{TRUE}, uses parallel computation where supported. Default is \code{FALSE}.}
  \item{cores}{Optional integer; number of workers when \code{parallel = TRUE}.}
}
\details{
This function:
\enumerate{
  \item Extracts ATT estimates and variance-covariance matrices from \code{fect.out}.
  \item Constructs weighted averages across post-treatment periods.
  \item Computes robust confidence sets under RM and C-LF assumptions.
  \item Optionally computes robust bounds for each post-treatment period separately.
}
Robust sets are computed using functions from the \pkg{HonestDiDFEct} package.
}
\value{
An updated \code{fect} object including:
\item{sensitivity.rm}{Relative-magnitude results, a list: \code{results} (robust confidence sets for the average effect, one row per \code{Mbar}), \code{original} (the original confidence set for the average effect) and \code{periods} (robust sets for each post-treatment period, one row per period and \code{Mbar}). Parts not requested are \code{NULL}.}
\item{sensitivity.smooth}{Smoothness results, with the same three parts, by \code{M}.}
}
\author{
Rivka Lipkovitz
}
\examples{
\donttest{
## HonestDiDFEct is not on CRAN: remotes::install_github("lzy318/HonestDiDFEct")
if (requireNamespace("HonestDiDFEct", quietly = TRUE)) {
  library(fect)
  data(simdata)
  ## The sensitivity analysis uses pre-treatment placebo estimates, so fit
  ## with placeboTest = TRUE (three placebo periods here).
  out <- fect(Y ~ D, data = simdata, index = c("id", "time"), method = "fe",
              se = TRUE, nboots = 50, parallel = FALSE,
              placeboTest = TRUE, placebo.period = c(-2, 0))
  if (!is.matrix(out$att.vcov) && is.matrix(out$att.boot)) {
    out$att.vcov <- cov(t(out$att.boot), use = "pairwise.complete.obs")
  }
  out_sens <- fect_sens(
    fect.out = out,
    post.periods = 1:4,
    Mbarvec = seq(0, 1, by = 0.2),
    Mvec = c(0, 0.05, 0.1),
    periodMbarvec = c(0, 0.5),
    periodMvec = c(0, 0.1)
  )
  names(out_sens$sensitivity.rm)
  names(out_sens$sensitivity.smooth)
}
}
}
