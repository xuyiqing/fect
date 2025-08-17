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
  \item{X}{Optional covariate vector for adjustment.}
  \item{index}{Character vector of unit and time variable names, e.g., \code{c("id", "time")}.}
  \item{method}{DiD method: \code{"twfe"}, \code{"st"}, \code{"iw"}, \code{"cs_never"}, \code{"cs_notyet"}, or \code{"didm"}.}
  \item{se}{Standard error method: \code{"default"}, \code{"boot"}, \code{"bootstrap"}, or \code{"jackknife"}.}
  \item{nboots}{Number of bootstrap replications (if applicable).}
  \item{parallel}{Logical; use parallel computation for bootstrapping.}
  \item{core}{Number of CPU cores to use if \code{parallel = TRUE}.}
  \item{time_to_treat_var}{Name of time-to-treatment variable; used internally.}
  \item{treat_indicator}{Name of treatment cohort indicator variable; used internally.}
  \item{csdid.base_period}{Baseline period choice for Callaway–Sant'Anna estimators.}
  \item{didm.effects}{Effects vector for \code{didm}; required if \code{method = "didm"}.}
  \item{didm.placebo}{Placebo vector for \code{didm}; required if \code{method = "didm"}.}
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
\item{est.avg}{Data frame with overall ATT, standard error, confidence interval, and p-value.}
\item{est.att}{Event-study ATT estimates by relative period, including standard errors and confidence intervals.}
}
\author{
Rivka Lipkovitz
}
\examples{
\dontrun{
result_twfe <- did_wrapper(
  data   = df,
  Y      = "outcome",
  D      = "treat",
  index  = c("id", "time"),
  method = "twfe"
)
result_twfe$est.avg
result_twfe$est.att
}
}
