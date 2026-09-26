\name{imputed_outcomes}
\alias{imputed_outcomes}
\title{Imputed Potential Outcomes Accessor}
\description{
  Returns the cell-level imputed potential-outcome surface from a fitted
  \code{fect} object as a long-form data frame. This is the rock-bottom
  accessor for any post-hoc estimand: read off the columns and aggregate
  however you like. Used internally by \code{\link{estimand}} for the
  shipped estimand types; exposed publicly so users can compose custom
  estimands the package does not ship.
}

\usage{
imputed_outcomes(
  fit,
  cells = NULL,
  replicates = FALSE,
  direction = c("on", "off")
)
}

\arguments{
  \item{fit}{A \code{fect} object (output of \code{\link{fect}}).}

  \item{cells}{Optional filter on which treated cells to return. Accepts
    \code{NULL} (default; all treated cells), a logical vector aligned
    with the unfiltered row order, or a one-sided formula evaluated
    against the long-form data
    (e.g., \code{~ event.time \%in\% 1:5 & !id \%in\% bad_ids}).}

  \item{replicates}{Logical. \code{FALSE} (default) returns one row per
    treated cell with the point-estimate \code{Y0_hat}. \code{TRUE}
    returns the selected treated cells of every bootstrap or parametric
    replicate, with a \code{replicate} column; requires the fit to have
    been built with \code{keep.sims = TRUE}, and is not available for
    jackknife fits. Replicate \eqn{b} is a refit on a resampled panel
    whose units are \code{fit$colnames.boot[[b]]}, so its rows are the
    copies of the selected cells that it contains: under the case
    bootstrap a unit drawn twice appears twice and an undrawn unit not at
    all, so the number of rows varies by replicate; under the parametric
    bootstrap every treated cell appears once per replicate.}

  \item{direction}{Either \code{"on"} (default; event time relative to
    treatment onset) or \code{"off"} (relative to treatment exit; only
    meaningful on reversal panels).}
}

\value{
  A data frame with one row per treated cell, or, with
  \code{replicates = TRUE}, one row per copy of a treated cell in each
  replicate (rows in replicate order). Columns:
  \describe{
    \item{\code{id}}{Unit identifier from \code{fit$id}.}
    \item{\code{time}}{Calendar time from \code{fit$rawtime}.}
    \item{\code{event.time}}{Event time at this cell (relative to onset
      or exit per \code{direction}).}
    \item{\code{cohort}}{First-treatment calendar time for the unit
      (\code{direction = "on"}) or first-exit calendar time
      (\code{direction = "off"}).}
    \item{\code{treated}}{Logical; always \code{TRUE} in returned rows.}
    \item{\code{Y_obs}}{Observed outcome at this cell.}
    \item{\code{Y0_hat}}{Imputed counterfactual outcome
      (\code{Y_obs - eff}).}
    \item{\code{eff}}{Cell-level score: contribution to the ATT estimator.
      For imputation estimators this equals \code{Y_obs - Y0_hat}; for
      doubly-robust estimators it equals
      \code{(Y_obs - Y0_hat) + eff_debias}.}
    \item{\code{eff_debias}}{Debias correction; 0 for plain imputation
      estimators, populated for DR estimators.}
    \item{\code{W.agg}}{Aggregation weight at this cell: the weight given
      by \code{W} or \code{W.agg} (stored on the fit as
      \code{fit$W.agg}); 1 if the fit was built without \code{W} or
      \code{W.agg}, or with \code{W.est} only.}
    \item{\code{replicate}}{(only when \code{replicates = TRUE}) the
      replicate number, 1 to the number of replicates kept. Within a
      replicate, \code{eff} and \code{Y0_hat} are that replicate's values,
      and \code{Y_obs} and \code{W.agg} those of the source cell. Averaging
      \code{eff} within a replicate (unweighted fit, no \code{cells}
      filter) gives that replicate's \code{att.avg.boot} draw.}
  }
}

\section{Memory cost}{
  With \code{replicates = TRUE} the returned data frame has about
  \code{n_treated_cells * nboots} rows (exactly that many for the
  parametric bootstrap; under the case bootstrap the count varies by
  replicate). For typical panels this is manageable; for large panels
  (\eqn{TT \times N \ge 50{,}000} and \code{nboots} \eqn{\ge 500})
  consider filtering via \code{cells} before expansion.
}

\section{Slot contract}{
  This function reads from \code{fit$Y.dat}, \code{fit$D.dat},
  \code{fit$T.on} (or \code{T.off}), \code{fit$eff},
  \code{fit$eff.boot} and \code{fit$colnames.boot} (when
  \code{replicates = TRUE}),
  \code{fit$eff_debias} (when populated), \code{fit$W.agg} (fits
  weighted with \code{W} or \code{W.agg}),
  \code{fit$id}, \code{fit$rawtime}. These slots have a documented
  stable shape; see the package design notes.
}

\seealso{
  \code{\link{estimand}} for the typed dispatcher built on top of this
  accessor; \code{\link{fect}} for the fitting interface.
}

\examples{
\donttest{
library(fect)
data(simdata)
fit <- fect(Y ~ D, data = simdata, index = c("id", "time"),
            method = "fe", force = "two-way",
            se = TRUE, nboots = 50, parallel = FALSE,
            keep.sims = TRUE)

## Point-estimate long form: one row per treated cell.
po <- imputed_outcomes(fit)
head(po)

## Filter to first 5 event times.
po5 <- imputed_outcomes(fit, cells = ~ event.time \%in\% 1:5)

## Bootstrap replicate expansion (requires keep.sims = TRUE): the rows of
## replicate b are the copies of the treated cells that replicate b drew.
po_b <- imputed_outcomes(fit, replicates = TRUE)
head(table(po_b$replicate))    # rows per replicate vary

## Averaging eff within a replicate gives its overall ATT draw.
head(tapply(po_b$eff, po_b$replicate, mean))
head(c(fit$att.avg.boot))
}
}
