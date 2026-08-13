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
    expands by bootstrap/jackknife replicate; requires the fit to have
    been built with \code{keep.sims = TRUE}.}

  \item{direction}{Either \code{"on"} (default; event time relative to
    treatment onset) or \code{"off"} (relative to treatment exit; only
    meaningful on reversal panels).}
}

\value{
  A data frame with one row per treated cell or per (treated cell,
  replicate). Columns:
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
    \item{\code{W.agg}}{Aggregation weight at this cell; 1 if the fit
      was built without \code{W} or \code{W.agg}.}
    \item{\code{replicate}}{(only when \code{replicates = TRUE})
      1..\code{nboots} for bootstrap, or the dropped-unit index for
      jackknife.}
  }
}

\section{Memory cost}{
  With \code{replicates = TRUE} the returned data frame has
  \code{n_treated_cells * nboots} rows. For typical panels this is
  manageable; for large panels (\eqn{TT \times N \ge 50{,}000} and
  \code{nboots} \eqn{\ge 500}) consider filtering via \code{cells}
  before expansion.
}

\section{Slot contract}{
  This function reads from \code{fit$Y.dat}, \code{fit$D.dat},
  \code{fit$T.on} (or \code{T.off}), \code{fit$eff},
  \code{fit$eff.boot} (when \code{replicates = TRUE}),
  \code{fit$eff_debias} (when populated), \code{fit$W.agg},
  \code{fit$id}, \code{fit$rawtime}. These slots have a documented
  stable shape; see the package design notes.
}

\seealso{
  \code{\link{estimand}} for the typed dispatcher built on top of this
  accessor; \code{\link{fect}} for the fitting interface.
}

\examples{
\dontrun{
library(fect)
fit <- fect(Y ~ D, data = simdata, index = c("id", "time"),
            method = "fe", force = "two-way",
            se = TRUE, nboots = 200, parallel = FALSE,
            keep.sims = TRUE)

## Point-estimate long form: one row per treated cell.
po <- imputed_outcomes(fit)
head(po)

## Filter to first 5 event times.
po5 <- imputed_outcomes(fit, cells = ~ event.time \%in\% 1:5)

## Bootstrap replicate expansion (requires keep.sims = TRUE).
po_b <- imputed_outcomes(fit, replicates = TRUE)
nrow(po_b) == nrow(po) * 200    # one row per (cell, replicate)
}
}
