\name{fect_mspe}
\alias{fect_mspe}
\title{Mean Squared Prediction Error Evaluation for fect Objects}
\description{
Evaluates the prediction accuracy of one or more \code{fect} model fits by hiding a subset of untreated observations, refitting, and comparing the model's predictions there with the actual values.
}
\usage{
fect_mspe(
  out.fect,
  seed = NULL,
  cv.method = "rolling",
  cv.nobs = 3,
  cv.donut = 1,
  cv.buffer = 1,
  cv.prop = 0.1,
  min.T0 = 5,
  k = 20,
  criterion = "mspe",
  W = NULL,
  norm.para = NULL,
  proportion = 0
)
}
\arguments{
  \item{out.fect}{A fitted \code{fect} object or a named list of fitted \code{fect} objects to compare.}
  \item{seed}{Optional integer; random seed for reproducibility.}
  \item{cv.method}{Character; cross-validation masking strategy. One of \code{"rolling"} (default) or \code{"block"} (random scattered anchors with contiguous-block masking). Legacy aliases \code{"all_units"} (= \code{"block"}) and \code{"treated_units"} are still accepted but emit a deprecation message.}
  \item{cv.nobs}{Integer; number of observations to mask per unit per fold. Default is 3.}
  \item{cv.donut}{Integer; number of periods around treatment onset to exclude from masking. Used by block CV. Default is 1.}
  \item{cv.buffer}{Integer; number of past-side cells masked from training (but not scored) before each rolling-window holdout. Used only by \code{cv.method = "rolling"}. Default is 1.}
  \item{cv.prop}{Numeric; for block CV, proportion of observations to mask per round; for rolling CV, fraction of eligible units sampled per fold. Default is 0.1.}
  \item{min.T0}{Integer; minimum number of pre-treatment periods required. Default is 5.}
  \item{k}{Integer; number of cross-validation folds. Default is 20.}
  \item{criterion}{Character; scoring criterion. One of \code{"mspe"}, \code{"wmspe"}, \code{"gmspe"}, \code{"wgmspe"}, \code{"mad"}, \code{"moment"}, \code{"gmoment"}. Default is \code{"mspe"}.}
  \item{W}{Optional TT x N observation weight matrix. Default is \code{NULL}.}
  \item{norm.para}{Optional normalization vector. Default is \code{NULL}.}
  \item{proportion}{Numeric; proportion cutoff for count.T.cv (same as fect_cv). Default is 0.}
}
\details{
Each fit is refitted with the function that made it (\code{fect()}, or \code{gsynth()} for a fit made by \pkg{gsynth}), with its own arguments, the number of factors it used (\code{r.cv}), \code{CV = FALSE} and \code{se = FALSE}, on the data with the held-out cells hidden. The data and the arguments named in the fit's call are looked up where \code{fect_mspe()} is called, and for a fit made with a formula also where the formula was written; for the data, a value of that name that is not a data frame (such as the function \code{df}) is skipped. A \pkg{gsynth} fit uses the data stored in it. Call \code{fect_mspe()} where those objects are visible.

Each held-out cell is scored with the refit's full prediction there: the intercept, the unit and time effects, the cell's covariates times their coefficients, the factors and, for \code{method = "cfe"}, the other CFE terms. (Before version 2.4.7 the scores left out the covariates, and for never-treated fits also the intercept and the fixed effects.) A held-out cell whose covariates are missing in the data is not scored. For CFE fits with extra fixed effects in \code{index}, the refit cannot place a held-out cell in its group, so the prediction there uses the effect of the group of hidden cells instead of the cell's own group. A never-treated CFE fit with \code{Z}, \code{Q} or extra fixed effects in \code{index} does not store those terms for treated units, so it is not scored at held-out cells of treated units (such cells are held out only when a not-yet-treated fit comes first in the list).
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
\donttest{
library(fect)
data(simdata)
out <- fect(Y ~ D, data = simdata, index = c("id", "time"), method = "ife",
            r = 2, CV = FALSE, se = FALSE, parallel = FALSE)
mspe <- fect_mspe(out.fect = out, seed = 42, k = 5)
mspe$summary
}
}
