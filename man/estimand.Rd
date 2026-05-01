\name{estimand}
\alias{estimand}
\title{Post-hoc Estimand Dispatcher}
\description{
  Computes a post-hoc estimand from a \code{fect} fit, with bootstrap or
  jackknife uncertainty. The \code{type} argument selects from a closed
  enum of mathematically-defined estimands; the \code{by} argument
  selects the grouping axis. The accessor \code{\link{imputed_outcomes}}
  is the underlying long-form data source for any estimand the
  dispatcher does not ship natively.
}

\usage{
estimand(
  fit,
  type   = c("att", "att.cumu", "aptt", "log.att"),
  by     = c("event.time", "cohort", "calendar.time", "overall"),
  cells       = NULL,
  weights     = NULL,
  window      = NULL,
  direction   = c("on", "off"),
  vartype     = c("bootstrap", "jackknife", "parametric", "none"),
  conf.level  = 0.95,
  ci.method   = c("basic", "percentile")
)
}

\arguments{
  \item{fit}{A \code{fect} object.}

  \item{type}{Estimand type. \code{"att"} returns per-cell mean
    treatment effect aggregated per group; \code{"att.cumu"} returns
    cumulative ATT through each event time (replaces \code{\link{effect}});
    \code{"aptt"} returns the average proportional treatment effect on
    the treated (Chen and Roth 2024 QJE); \code{"log.att"} returns the
    mean log-scale treatment effect. \code{"aptt"} and \code{"log.att"}
    require \code{keep.sims = TRUE} at fit time for SE/CI.}

  \item{by}{Grouping axis. One of \code{"event.time"} (default),
    \code{"cohort"}, \code{"calendar.time"}, \code{"overall"}, or any
    column name resolvable in the fit's panel data.}

  \item{cells}{Optional filter on which treated cells to include.
    Accepts \code{NULL} (default), a logical vector, or a one-sided
    formula. See \code{\link{imputed_outcomes}}.}

  \item{weights}{Aggregation-weight handling. \code{NULL} (default) uses
    \code{fit$W.agg} if the fit was built with \code{W} or \code{W.agg};
    otherwise uniform.}

  \item{window}{Optional event-time window \code{c(L, R)}; convenience
    sugar for \code{cells = ~ event.time >= L & event.time <= R}.}

  \item{direction}{Either \code{"on"} (default) or \code{"off"}; see
    \code{\link{imputed_outcomes}}.}

  \item{vartype}{\code{"bootstrap"} (default), \code{"jackknife"},
    \code{"parametric"}, or \code{"none"}. Selects which variance
    method to source replicates from. The output \code{vartype}
    column reports the method actually used at fit time (read from
    \code{fit$vartype}), which may differ from this argument value
    if the fit was produced with a different setting --- the argument
    is informational and does not re-aggregate replicates.}

  \item{conf.level}{Two-sided confidence level. Defaults to 0.95.}

  \item{ci.method}{\code{"basic"} (reflected; matches fect's existing
    \code{est.att} convention; default) or \code{"percentile"}.}
}

\value{
  A data frame with columns \code{<by_key>}, \code{estimate}, \code{se},
  \code{ci.lo}, \code{ci.hi}, \code{n_cells}, and \code{vartype}.
  Always tidy regardless of \code{type} or \code{by}.
}

\section{Numerical equality with existing slots}{
  \code{estimand(fit, "att", "event.time")} returns \code{estimate},
  \code{se}, \code{ci.lo}, \code{ci.hi} byte-identical to columns
  \code{ATT}, \code{S.E.}, \code{CI.lower}, \code{CI.upper} of
  \code{fit$est.att}, when default arguments are used. This invariant
  is asserted by package tests.
}

\seealso{
  \code{\link{imputed_outcomes}} for the underlying long-form accessor;
  \code{\link{fect}} for the fitting interface.
}

\examples{
\dontrun{
library(fect)
fit <- fect(Y ~ D, data = simdata, index = c("id", "time"),
            method = "fe", force = "two-way",
            se = TRUE, nboots = 200, parallel = FALSE)

## Default: per-event-time ATT (matches fit$est.att numerically).
est <- estimand(fit, "att", "event.time")
head(est)
}
}

\references{
  Chen, Jiafeng, and Jonathan Roth. 2024. "Logs with Zeros? Some Problems
  and Solutions." \emph{Quarterly Journal of Economics} 139 (2):
  891--936.
}
