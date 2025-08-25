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
  \item{fect.out}{A fitted \code{fect} object.}
  \item{post.periods}{Vector of post-treatment periods for sensitivity analysis. Default uses all available post-treatment periods.}
  \item{l_vec}{Optional weighting vector for averaging ATT across post-treatment periods. Default weights by treated-unit counts.}
  \item{Mbarvec}{Values of \code{Mbar} for overall RM-based sensitivity analysis.}
  \item{Mvec}{Values of \code{M} for overall smoothness-based sensitivity analysis.}
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
\item{sensitivity.rm}{Relative Magnitude (RM) sensitivity results for average and period-by-period ATT.}
\item{sensitivity.smooth}{Smoothness (C-LF) sensitivity results for average and period-by-period ATT.}
}
\author{
Rivka Lipkovitz
}
\examples{
\dontrun{
out <- fect(Y ~ D, data = df, index = c("unit", "time"), method = "fe", se = TRUE)
out_sens <- fect_sens(
  fect.out = out,
  post.periods = c(1, 2, 3, 4),
  Mbarvec = seq(0, 1, by = 0.2),
  Mvec = c(0, 0.05, 0.1),
  periodMbarvec = c(0, 0.5),
  periodMvec = c(0, 0.1)
)
names(out_sens$sensitivity.rm)
names(out_sens$sensitivity.smooth)
}
}
