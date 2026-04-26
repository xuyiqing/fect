\name{fect_mspe}
\alias{fect_mspe}
\title{Mean Squared Prediction Error Evaluation for fect Objects}
\description{
Evaluates the prediction accuracy of one or more \code{fect} model fits by hiding a subset of control observations and comparing counterfactual predictions against actual values.
}
\usage{
fect_mspe(
  out.fect,
  seed = NULL,
  cv.method = "rolling",
  cv.nobs = 3,
  cv.donut = 1,
  cv.buffer = 1,
  cv.prop = 0.2,
  min.T0 = 5,
  k = 5,
  criterion = "mspe",
  W = NULL,
  norm.para = NULL,
  proportion = 0
)
}
\arguments{
  \item{out.fect}{A fitted \code{fect} object or a named list of fitted \code{fect} objects to compare.}
  \item{seed}{Optional integer; random seed for reproducibility.}
  \item{cv.method}{Character; cross-validation masking strategy. One of \code{"rolling"} (default), \code{"block"} (random scattered anchors with contiguous-block masking), or \code{"loo"}. Legacy aliases \code{"all_units"} (= \code{"block"}) and \code{"treated_units"} are still accepted but emit a deprecation message.}
  \item{cv.nobs}{Integer; number of observations to mask per unit per fold. Default is 3.}
  \item{cv.donut}{Integer; number of periods around treatment onset to exclude from masking. Used by block CV. Default is 1.}
  \item{cv.buffer}{Integer; number of past-side cells masked from training (but not scored) before each rolling-window holdout. Used only by \code{cv.method = "rolling"}. Default is 1.}
  \item{cv.prop}{Numeric; for block CV, proportion of observations to mask per round; for rolling CV, fraction of eligible units sampled per fold. Default is 0.2.}
  \item{min.T0}{Integer; minimum number of pre-treatment periods required. Default is 5.}
  \item{k}{Integer; number of cross-validation folds. Default is 5.}
  \item{criterion}{Character; scoring criterion. One of \code{"mspe"}, \code{"wmspe"}, \code{"gmspe"}, \code{"wgmspe"}, \code{"mad"}, \code{"moment"}, \code{"gmoment"}. Default is \code{"mspe"}.}
  \item{W}{Optional TT x N observation weight matrix. Default is \code{NULL}.}
  \item{norm.para}{Optional normalization vector. Default is \code{NULL}.}
  \item{proportion}{Numeric; proportion cutoff for count.T.cv (same as fect_cv). Default is 0.}
}
\value{
A list containing:
\item{summary}{A data frame with mean scores across replications for each model, including MSPE, WMSPE, GMSPE, WGMSPE, MAD, Moment, GMoment, RMSE, and Bias.}
\item{records}{A data frame with per-replication results including Rep, Model, Hidden_N, and all score columns.}
\item{fits}{The refitted \code{fect} objects from the last replication.}
\item{criterion}{The scoring criterion used.}
\item{scores}{The scores from the last replication.}
}
\examples{
\dontrun{
out <- fect(Y ~ D, data = df, index = c("unit", "time"), method = "ife", r = 2)
mspe <- fect_mspe(out.fect = out, hide_n = 10, seed = 42, n_rep = 5)
mspe$summary
}
}
