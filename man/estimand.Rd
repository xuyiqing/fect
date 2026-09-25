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
  test        = c("none", "placebo", "carryover"),
  cells       = NULL,
  weights     = NULL,
  window      = NULL,
  direction   = c("on", "off"),
  vartype     = c("bootstrap", "jackknife", "parametric", "none"),
  conf.level  = 0.95,
  ci.method   = NULL
)
}

\arguments{
  \item{fit}{A \code{fect} object.}

  \item{type}{Estimand type. \code{"att"} returns per-cell mean
    treatment effect aggregated per group; \code{"att.cumu"} returns
    cumulative ATT through each event time, the running sum of the
    per-period ATTs (replaces \code{\link{effect}});
    \code{"aptt"} returns the average proportional treatment effect on
    the treated (Chen and Roth 2024 QJE); \code{"log.att"} returns the
    mean log-scale treatment effect. SEs and CIs need the replicate
    draws saved with \code{keep.sims = TRUE} at fit time, except for
    \code{"att"} with \code{by = "event.time"} and default arguments and
    \code{"att.cumu"} with \code{by = "overall"}, which use the
    pre-aggregated draws that every \code{se = TRUE} fit keeps. Without
    them, the other calls stop.}

  \item{by}{Grouping axis. One of \code{"event.time"} (default),
    \code{"cohort"}, \code{"calendar.time"}, \code{"overall"}, or any
    column name resolvable in the fit's panel data.}

  \item{test}{Selects which subset of cells to aggregate over.
    \code{"none"} (default) uses the standard treated post-treatment
    cells. \code{"placebo"} restricts aggregation to the
    pre-treatment cells in \code{fit$placebo.period} that were
    masked-and-imputed during the placebo fit; e.g.\
    \code{estimand(fit, "aptt", "event.time", test = "placebo")}
    returns a per-event-time placebo APTT series. Requires
    \code{placeboTest = TRUE} at fit time and forces
    \code{direction = "on"}. \code{"carryover"} is the analogous
    extension for early post-reversal cells in
    \code{fit$carryover.period}; requires \code{carryoverTest = TRUE}
    at fit time, a panel with reversals, and forces
    \code{direction = "off"}. Both are incompatible with
    \code{type = "att.cumu"}.}

  \item{cells}{Optional filter on which treated cells to include.
    Accepts \code{NULL} (default), a logical vector, or a one-sided
    formula. See \code{\link{imputed_outcomes}}.}

  \item{weights}{Aggregation-weight handling. \code{NULL} (default) uses
    \code{fit$W.agg} if the fit was built with \code{W} or \code{W.agg};
    otherwise uniform.}

  \item{window}{Optional event-time window \code{c(L, R)}; convenience
    sugar for \code{cells = ~ event.time >= L & event.time <= R}. It
    works only with \code{by = "overall"}; the event-time series stop
    when it is given.}

  \item{direction}{Either \code{"on"} (default) or \code{"off"}; see
    \code{\link{imputed_outcomes}}.}

  \item{vartype}{\code{"bootstrap"} (default), \code{"jackknife"},
    \code{"parametric"}, or \code{"none"}. Selects which variance
    method to source replicates from. The output \code{vartype}
    column reports the method actually used at fit time (read from
    \code{fit$vartype}); this argument is informational and does not
    re-aggregate replicates. Wild-bootstrap fits (\code{vartype =
    "wild"} at fit time) are consumed transparently and reported as
    \code{"wild"} in the output column.}

  \item{conf.level}{Two-sided confidence level. Defaults to 0.95. Not
    used for \code{"att.cumu"}; the per-event-time ATT without a test
    accepts only 0.95.}

  \item{ci.method}{One of \code{"basic"}, \code{"percentile"},
    \code{"bc"} (bias-corrected percentile), \code{"bca"}
    (bias-corrected accelerated; Efron 1987 in full), or
    \code{"normal"} (Wald: \eqn{\hat\theta \pm z \cdot SE}). Default
    is \code{NULL}, which triggers a per-type default:
    \code{"att"} -> \code{"normal"} (matches what \code{fit$est.att}
    already uses internally),
    \code{"aptt"} -> \code{"bca"} and \code{"log.att"} ->
    \code{"bca"} (ratio / log estimators benefit from joint
    bias and acceleration corrections when the bootstrap distribution
    is skewed; \code{"bc"} can degenerate at the boundary when the
    point estimate falls outside the bootstrap support). Pass an
    explicit value to override. \code{"bca"} requires bootstrap- or
    wild-bootstrap-sourced replicates. For \code{"att.cumu"} the
    interval comes from \code{\link{att.cumu}} (\code{by = "overall"})
    or \code{\link{effect}} (\code{by = "event.time"}), whatever
    \code{ci.method} is: bootstrap percentiles for bootstrap fits, and
    for parametric fits the estimate plus or minus a critical value
    (normal for \code{"overall"}, t for \code{"event.time"}) times
    the SE. \code{NULL} still resolves to \code{"basic"} there, so
    \code{estimand()} warns about tail quantiles when \code{nboots} is
    below 1000, and a jackknife fit needs \code{ci.method = "normal"}.
    For \code{"att"} with \code{by = "event.time"} and no test, only
    the default is available (it returns \code{fit$est.att}); other
    values stop, so use \code{by = "overall"}.}
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

  \code{estimand(fit, "att", "overall")$se} equals the SE in
  \code{fit$est.avg} for unweighted bootstrap and parametric fits,
  whatever the
  order of the treated units in the panel: every replicate is read
  through its own units (\code{fit$colnames.boot}).
}

\seealso{
  \code{\link{imputed_outcomes}} for the underlying long-form accessor;
  \code{\link{fect}} for the fitting interface.
}

\examples{
\donttest{
library(fect)
data(simdata)
fit <- fect(Y ~ D, data = simdata, index = c("id", "time"),
            method = "fe", force = "two-way",
            se = TRUE, nboots = 50, parallel = FALSE)

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
