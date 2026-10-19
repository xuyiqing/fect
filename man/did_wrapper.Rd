\name{did_wrapper}
\alias{did_wrapper}
\title{A Multi-Method Difference-in-Differences Estimator}
\description{
Implements multiple difference-in-differences (DiD) estimators under a unified interface, supporting overall ATT and event-study estimates across staggered adoption settings.
}
\usage{
did_wrapper(
  data,
  Y,
  D,
  X = NULL,
  index,
  method = c("twfe", "st", "iw", "cs_never", "cs_notyet", "didm"),
  se = c("default", "boot", "bootstrap", "jackknife"),
  nboots = 200,
  parallel = TRUE,
  core = NULL,
  time_to_treat_var = "Time_to_Treatment",
  treat_indicator = "treat",
  csdid.base_period = "universal",
  didm.effects = NA,
  didm.placebo = NA
)
}
\arguments{
  \item{data}{Input data frame.}
  \item{Y}{Outcome variable name (string).}
  \item{D}{Treatment indicator variable name (string).}
  \item{X}{Optional character vector of covariate column names. Used by every method except \code{"st"}, which ignores it.}
  \item{index}{Character vector of unit and time variable names, e.g., \code{c("id", "time")}.}
  \item{method}{DiD method: \code{"twfe"} (two-way fixed effects, the default), \code{"st"} (stacked DiD), \code{"iw"} (Sun and Abraham's interaction-weighted estimator, via \code{fixest::sunab}), \code{"cs_never"} or \code{"cs_notyet"} (Callaway and Sant'Anna with never-treated or not-yet-treated controls; needs \pkg{did} 2.5.1 or later), or \code{"didm"} (de Chaisemartin and D'Haultfoeuille, via \pkg{DIDmultiplegtDYN}).}
  \item{se}{Standard error method. \code{"default"} (the default) uses each estimator's own standard errors: clustered by unit for \code{"twfe"}, \code{"st"} and \code{"iw"}, and the defaults of \pkg{did} and \pkg{DIDmultiplegtDYN} for the other methods, with a 95\% interval. \code{"boot"} (or \code{"bootstrap"}) is a cluster bootstrap that draws units with replacement (a unit drawn twice enters as two units), with \code{nboots} draws, a basic (reflected) 95\% interval, and a p-value from the share of draws on the other side of zero. \code{"jackknife"} drops one unit at a time; its standard error is \eqn{\sqrt{(n-1)/n \sum_i (\hat\theta_{(i)} - \bar\theta)^2}} and the interval and p-value are normal. Before version 2.4.8, the bootstrap kept a unit drawn several times only once and the jackknife omitted the \eqn{(n-1)/n} scaling, so both gave standard errors that were too small.}
  \item{nboots}{Number of bootstrap replications for \code{se = "boot"}. For \code{se = "jackknife"} it is set to the number of units.}
  \item{parallel}{Logical; use parallel computation for bootstrapping.}
  \item{core}{Number of CPU cores to use if \code{parallel = TRUE}.}
  \item{time_to_treat_var}{Name of the event-time column created by \code{\link{get.cohort}}; change it only if your data already use that name for something else.}
  \item{treat_indicator}{Name of the ever-treated indicator column used by the \code{"twfe"} and \code{"st"} event studies; created if absent.}
  \item{csdid.base_period}{Base period for the Callaway--Sant'Anna estimators, passed to \code{did::att_gt(base_period = )}: \code{"universal"} (default; every event time compared with the period just before treatment) or \code{"varying"}.}
  \item{didm.effects}{Number of post-treatment effects to estimate with \code{method = "didm"}; required for that method.}
  \item{didm.placebo}{Number of placebo (pre-treatment) estimates with \code{method = "didm"}; required for that method.}
}
\details{
This function:
\enumerate{
  \item Drops always-treated units.
  \item Constructs event time and cohort variables.
  \item Computes ATT using the specified DiD method.
  \item Constructs event-study ATT curves.
  \item Optionally estimates bootstrap or jackknife standard errors.
}
Supported methods include TWFE, stacked DiD, interaction-weighted DiD, Callaway–Sant'Anna estimators, placebo DiD, and DIDmultiplegt.
}
\value{
A list of class \code{"did_wrapper"}:
\item{est.avg}{A one-row data frame with columns \code{ATT.avg}, \code{S.E.}, \code{CI.lower}, \code{CI.upper} and \code{p.value} (95\% interval).}
\item{est.att}{A data frame of event-study estimates with columns \code{ATT}, \code{S.E.}, \code{CI.lower}, \code{CI.upper}, \code{p.value} and \code{count} (number of treated units). Row names are event times, with 1 the first treated period and 0 the last period before treatment, as in \code{fect}. Can be passed to \code{\link{esplot}}.}
}
\author{
Rivka Lipkovitz
}
\examples{
\donttest{
if (requireNamespace("fixest", quietly = TRUE)) {
  library(fect)
  data(simdata)
  result_twfe <- did_wrapper(data = simdata, Y = "Y", D = "D",
                             index = c("id", "time"), method = "twfe")
  result_twfe$est.avg
  head(result_twfe$est.att)
}
}
}
