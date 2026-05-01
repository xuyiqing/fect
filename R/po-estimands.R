###############################################################################
## Post-hoc estimands API
##
## Public surface (v2.4.0):
##   imputed_outcomes(fit, cells, replicates, direction)
##   estimand(fit, type, by, cells, weights, direction, vartype, ...)
##
## This file owns the slot contract documented in
##   statsclaw-workspace/fect/ref/po-estimands-contract.md
## See section 3 (Slot contract on the fit object) for the invariants enforced by
## .validate_po_contract().
###############################################################################


## ---------------------------------------------------------------------------
## Soft-deprecation helper for v2.4.0
## ---------------------------------------------------------------------------

## Internal: emit a one-time-per-session message that the named function
## is soft-deprecated, pointing users at the unified estimand() API.
## Suppressed during internal calls from estimand() itself.
.fect_deprecation_message_once <- function(fn) {
    if (isTRUE(getOption("fect.suppress_estimand_deprecation", FALSE))) {
        return(invisible(NULL))
    }
    key <- paste0("fect.deprecated_warned_", fn)
    if (isTRUE(getOption(key, FALSE))) {
        return(invisible(NULL))
    }
    message(
        sprintf(
            "`%s()` is soft-deprecated as of fect 2.4.0; prefer the unified ",
            fn
        ),
        "`estimand(fit, \"att.cumu\", ...)` API for the same numerics. ",
        "See `?estimand` and the \"Alternative estimands\" vignette section. ",
        "This function will continue to work; removal not before v3.0.0."
    )
    do.call(options, setNames(list(TRUE), key))
    invisible(NULL)
}


## Internal: evaluate `expr` with the deprecation messages on effect()
## and att.cumu() temporarily suppressed. Used by .compute_att_cumu_*()
## so that calling estimand("att.cumu", ...) does not surface the
## deprecation message that fires only on user-facing calls.
.with_no_deprecation <- function(expr) {
    old <- getOption("fect.suppress_estimand_deprecation", FALSE)
    options(fect.suppress_estimand_deprecation = TRUE)
    on.exit(options(fect.suppress_estimand_deprecation = old), add = TRUE)
    force(expr)
}


## ---------------------------------------------------------------------------
## Slot contract
## ---------------------------------------------------------------------------

## Internal: validate that `fit` satisfies the slot contract.
## Called at the top of every public-API entry point. Raises an error
## that points at the contract on violation.
.validate_po_contract <- function(fit, require_replicates = FALSE) {
    if (!inherits(fit, "fect") && !is.list(fit)) {
        stop("`fit` must be a fect object (output of fect()).", call. = FALSE)
    }

    required <- c("Y.dat", "D.dat", "I.dat", "T.on",
                  "eff", "id", "rawtime")
    missing_slots <- setdiff(required, names(fit))
    if (length(missing_slots) > 0) {
        stop(
            "fit object missing required slots: ",
            paste(missing_slots, collapse = ", "),
            ". The new post-hoc estimand API requires fects produced by ",
            "v2.3.x or later. Refit with the current package version.",
            call. = FALSE
        )
    }

    if (!is.matrix(fit$Y.dat)) {
        stop("Slot contract: fit$Y.dat must be a matrix.", call. = FALSE)
    }
    TT <- nrow(fit$Y.dat)
    N  <- ncol(fit$Y.dat)

    for (slot in c("D.dat", "I.dat", "T.on", "eff")) {
        v <- fit[[slot]]
        if (!is.matrix(v) || !identical(dim(v), c(TT, N))) {
            stop(
                "Slot contract: fit$", slot,
                " must be a TT x N matrix matching fit$Y.dat (",
                TT, " x ", N, ").",
                call. = FALSE
            )
        }
    }

    if (length(fit$id) != N) {
        stop(
            "Slot contract: length(fit$id) (", length(fit$id),
            ") must equal ncol(fit$Y.dat) (", N, ").",
            call. = FALSE
        )
    }
    if (length(fit$rawtime) != TT) {
        stop(
            "Slot contract: length(fit$rawtime) (", length(fit$rawtime),
            ") must equal nrow(fit$Y.dat) (", TT, ").",
            call. = FALSE
        )
    }

    ## eff is the cell-level score per the contract:
    ##   eff[t, i] = Y_obs[t, i] - Y_hat[t, i]
    ## where Y_hat is the prediction from the control-fit. At treated
    ## cells this is the cell-level treatment effect estimate; at
    ## non-treated cells it is the model residual (typically small but
    ## not exactly zero). All ATT-style estimands aggregate eff over
    ## treated cells only.

    ## eff.boot, if populated, must match TT x N x nboots.
    if (!is.null(fit$eff.boot)) {
        eb <- fit$eff.boot
        if (!is.array(eb) || length(dim(eb)) != 3) {
            stop(
                "Slot contract: fit$eff.boot must be a 3D array.",
                call. = FALSE
            )
        }
        if (!identical(dim(eb)[1:2], c(TT, N))) {
            stop(
                "Slot contract: fit$eff.boot first two dimensions (",
                dim(eb)[1], " x ", dim(eb)[2],
                ") must match TT x N (", TT, " x ", N, ").",
                call. = FALSE
            )
        }
    }

    ## eff_debias, if populated, matches eff shape.
    if (!is.null(fit$eff_debias)) {
        eb <- fit$eff_debias
        if (!is.matrix(eb) || !identical(dim(eb), c(TT, N))) {
            stop(
                "Slot contract: fit$eff_debias must be a TT x N matrix ",
                "(or NULL).",
                call. = FALSE
            )
        }
    }

    if (require_replicates && is.null(fit$eff.boot)) {
        stop(
            "No bootstrap/jackknife results available. ",
            "Choose keep.sims = TRUE in fect().",
            call. = FALSE
        )
    }

    invisible(TRUE)
}


## ---------------------------------------------------------------------------
## Public: imputed_outcomes()
## ---------------------------------------------------------------------------

#' Imputed potential outcomes accessor
#'
#' Returns the cell-level imputed potential-outcome surface from a fect fit
#' as a long-form data frame. This is the rock-bottom accessor for any
#' post-hoc estimand: read off the columns and aggregate however you like.
#'
#' @param fit A fect object (output of \code{\link{fect}}).
#' @param cells Optional filter on which treated cells to return. Accepts
#'   \code{NULL} (default; all treated cells), a logical vector aligned
#'   with the unfiltered row order, or a one-sided rlang formula
#'   evaluated against the long-form data (e.g.,
#'   \code{~ event.time \%in\% 1:5 & !id \%in\% bad_ids}).
#' @param replicates Logical. \code{FALSE} (default) returns one row per
#'   treated cell with the point-estimate \code{Y0_hat}. \code{TRUE}
#'   expands by bootstrap/jackknife replicate; requires the fit to have
#'   been built with \code{keep.sims = TRUE}.
#' @param direction Either \code{"on"} (default; event time relative to
#'   treatment onset) or \code{"off"} (relative to treatment exit; only
#'   meaningful on reversal panels).
#'
#' @return A data frame with one row per (treated cell) or per
#'   (treated cell, replicate). Columns:
#'   \describe{
#'     \item{\code{id}}{unit identifier from \code{fit$id}.}
#'     \item{\code{time}}{calendar time from \code{fit$rawtime}.}
#'     \item{\code{event.time}}{event time at this cell (relative to onset
#'       or exit per \code{direction}).}
#'     \item{\code{cohort}}{first-treatment calendar time for the unit
#'       (\code{direction = "on"}) or first-exit calendar time
#'       (\code{direction = "off"}).}
#'     \item{\code{treated}}{logical; always \code{TRUE} in the returned
#'       rows (untreated cells are filtered out).}
#'     \item{\code{Y_obs}}{observed outcome at this cell.}
#'     \item{\code{Y0_hat}}{imputed counterfactual outcome
#'       (\code{Y_obs - eff}).}
#'     \item{\code{eff}}{cell-level score: contribution to the ATT
#'       estimator. For imputation estimators
#'       \code{eff = Y_obs - Y0_hat}; for doubly-robust estimators
#'       \code{eff = (Y_obs - Y0_hat) + eff_debias}.}
#'     \item{\code{eff_debias}}{debias correction; 0 for plain imputation
#'       estimators, populated for DR estimators.}
#'     \item{\code{W.agg}}{aggregation weight at this cell; 1 if the fit
#'       was built without \code{W} or \code{W.agg}.}
#'     \item{\code{replicate}}{(only when \code{replicates = TRUE})
#'       1..\code{nboots} for bootstrap, or the dropped-unit index for
#'       jackknife.}
#'   }
#'
#' @section Memory cost:
#'   With \code{replicates = TRUE} the returned data frame has
#'   \code{n_treated_cells * nboots} rows. For typical panels this is
#'   manageable; for large panels (\eqn{TT \times N \ge 50{,}000} and
#'   \eqn{nboots \ge 500}) consider filtering via \code{cells} before
#'   expansion.
#'
#' @seealso \code{\link{estimand}} for the typed dispatcher built on top
#'   of this accessor; \code{\link{fect}} for the fitting interface.
#'
#' @examples
#' \dontrun{
#' library(fect)
#' fit <- fect(Y ~ D, data = simdata, index = c("id", "time"),
#'             method = "fe", force = "two-way",
#'             se = TRUE, nboots = 200, parallel = FALSE,
#'             keep.sims = TRUE)
#'
#' ## Point-estimate long form: one row per treated cell.
#' po <- imputed_outcomes(fit)
#' head(po)
#'
#' ## Filter to first 5 event times.
#' po5 <- imputed_outcomes(fit, cells = ~ event.time \%in\% 1:5)
#'
#' ## Bootstrap replicate expansion (requires keep.sims = TRUE).
#' po_b <- imputed_outcomes(fit, replicates = TRUE)
#' nrow(po_b) == nrow(po) * 200    # one row per (cell, replicate)
#' }
#'
#' @export
imputed_outcomes <- function(fit,
                              cells = NULL,
                              replicates = FALSE,
                              direction = c("on", "off")) {

    direction <- match.arg(direction)
    .validate_po_contract(fit, require_replicates = replicates)

    Y      <- fit$Y.dat
    D      <- fit$D.dat
    eff    <- fit$eff
    eff_db <- if (is.null(fit$eff_debias)) {
                  matrix(0, nrow = nrow(eff), ncol = ncol(eff))
              } else fit$eff_debias

    TT <- nrow(Y)
    N  <- ncol(Y)

    Tev <- if (direction == "on") fit$T.on else fit$T.off
    if (is.null(Tev)) {
        stop(
            "direction = \"off\" requested, but fit$T.off is NULL ",
            "(this fit has no treatment exit). Use direction = \"on\".",
            call. = FALSE
        )
    }

    ## Treated-cell mask: cells where D = 1 (direction = "on") or where the
    ## unit has been treated at least once and is now post-exit (direction
    ## = "off"). Use D.dat == 1 for "on"; T.off non-NA flags exit-relative
    ## cells for "off".
    if (direction == "on") {
        treated_mask <- !is.na(D) & D == 1
    } else {
        treated_mask <- !is.na(Tev)
    }

    ## Cell coordinates of treated cells.
    idx <- which(treated_mask, arr.ind = TRUE)
    if (nrow(idx) == 0) {
        stop("No treated cells found for direction = \"", direction, "\".",
             call. = FALSE)
    }
    t_idx <- idx[, "row"]
    i_idx <- idx[, "col"]

    ## Cohort: calendar time of first treatment ("on") or first exit
    ## ("off") for each unit. NA for units that are never treated /
    ## never exit.
    cohort_per_unit <- if (direction == "on") {
        apply(D == 1 & !is.na(D), 2, function(d) {
            pos <- which(d)
            if (length(pos) == 0) NA_real_ else fit$rawtime[min(pos)]
        })
    } else {
        apply(Tev == 1 & !is.na(Tev), 2, function(d) {
            pos <- which(d)
            if (length(pos) == 0) NA_real_ else fit$rawtime[min(pos)]
        })
    }

    W_mat <- if (is.null(fit$W.agg)) {
                 matrix(1, nrow = TT, ncol = N)
             } else fit$W.agg

    df_point <- data.frame(
        id          = fit$id[i_idx],
        time        = fit$rawtime[t_idx],
        event.time  = Tev[cbind(t_idx, i_idx)],
        cohort      = cohort_per_unit[i_idx],
        treated     = TRUE,
        Y_obs       = Y[cbind(t_idx, i_idx)],
        Y0_hat      = Y[cbind(t_idx, i_idx)] - eff[cbind(t_idx, i_idx)],
        eff         = eff[cbind(t_idx, i_idx)],
        eff_debias  = eff_db[cbind(t_idx, i_idx)],
        W.agg       = W_mat[cbind(t_idx, i_idx)],
        stringsAsFactors = FALSE
    )

    ## Apply cells filter (point-form data) before optional replicate
    ## expansion to keep memory footprint down.
    if (!is.null(cells)) {
        df_point <- .apply_cells_filter(df_point, cells)
    }

    if (!replicates) {
        rownames(df_point) <- NULL
        return(df_point)
    }

    ## Replicate expansion: one row per (cell, b) where b indexes
    ## bootstrap or jackknife replicate. Y_obs and W.agg are constant
    ## across replicates; Y0_hat and eff vary per replicate.
    eb <- fit$eff.boot
    nboots <- dim(eb)[3]

    n_cells <- nrow(df_point)

    ## Re-derive cell coordinates against the (potentially-filtered)
    ## df_point so we map back into the original eff.boot array.
    t_idx_f <- match(df_point$time, fit$rawtime)
    i_idx_f <- match(df_point$id,   fit$id)

    df_rep <- df_point[rep(seq_len(n_cells), each = nboots), , drop = FALSE]
    df_rep$replicate <- rep(seq_len(nboots), times = n_cells)

    eff_rep <- vapply(
        seq_len(n_cells),
        function(k) eb[t_idx_f[k], i_idx_f[k], ],
        numeric(nboots)
    )
    ## eff_rep is nboots x n_cells; flatten in (cell-major, replicate-minor)
    ## order to match df_rep's row order.
    df_rep$eff    <- as.vector(eff_rep)
    df_rep$Y0_hat <- df_rep$Y_obs - df_rep$eff

    rownames(df_rep) <- NULL
    df_rep
}


## ---------------------------------------------------------------------------
## Public: estimand()
## ---------------------------------------------------------------------------

#' Post-hoc estimand dispatcher
#'
#' Computes a post-hoc estimand from a \code{fect} fit, with bootstrap or
#' jackknife uncertainty. The \code{type} argument selects from a closed
#' enum of mathematically-defined estimands; the \code{by} argument
#' selects the grouping axis. The accessor \code{\link{imputed_outcomes}}
#' is the underlying long-form data source for any estimand the
#' dispatcher does not ship natively.
#'
#' @param fit A fect object (output of \code{\link{fect}}).
#' @param type Estimand type. One of:
#'   \describe{
#'     \item{\code{"att"}}{Per-cell mean treatment effect, aggregated
#'       per group: \eqn{\mathrm{ATT}_g = \mathrm{mean}_{(t,i)\in g, D=1}(Y_{ti} - \widehat Y_{ti}(0))}.}
#'     \item{\code{"att.cumu"}}{Cumulative ATT through each event time.
#'       Replaces \code{\link{effect}} for the unified API.}
#'     \item{\code{"aptt"}}{Average proportional treatment effect on the
#'       treated (Chen & Roth 2024 QJE):
#'       \eqn{\mathrm{APTT}_g = \mathrm{mean}_g(Y - \widehat Y(0)) / \mathrm{mean}_g(\widehat Y(0))}.
#'       Requires \code{keep.sims = TRUE} at fit time.}
#'     \item{\code{"log.att"}}{Mean log-scale treatment effect:
#'       \eqn{\mathrm{logATT}_g = \mathrm{mean}_g(\log Y - \log \widehat Y(0))}.
#'       Requires \code{keep.sims = TRUE}.}
#'   }
#' @param by Grouping axis. One of \code{"event.time"} (default;
#'   per-event-time series), \code{"cohort"}, \code{"calendar.time"},
#'   \code{"overall"} (one row), or any column name resolvable in the
#'   fit's panel data.
#' @param test Selects which subset of cells to aggregate over (v2.4.3+).
#'   \code{"none"} (default) uses standard treated post-treatment
#'   cells. \code{"placebo"} restricts to pre-treatment cells in
#'   \code{fit$placebo.period} that were masked-and-imputed during the
#'   placebo fit; produces e.g.\ a per-event-time placebo APTT series
#'   for credibility checks. Requires \code{placeboTest = TRUE} at fit
#'   time; forces \code{direction = "on"}. \code{"carryover"} is the
#'   analogous extension on the reversal side; requires
#'   \code{carryoverTest = TRUE} + a panel with reversals; forces
#'   \code{direction = "off"}. Both incompatible with
#'   \code{type = "att.cumu"}.
#' @param cells Optional filter on which treated cells to include.
#'   Accepts \code{NULL} (default; all treated cells), a logical vector,
#'   or a one-sided formula. See \code{\link{imputed_outcomes}}.
#' @param weights Aggregation-weight handling. \code{NULL} (default) uses
#'   \code{fit$W.agg} if the fit was built with \code{W} or \code{W.agg};
#'   otherwise uniform.
#' @param window Optional event-time window \code{c(L, R)}; convenience
#'   sugar for \code{cells = ~ event.time >= L & event.time <= R}.
#' @param direction Either \code{"on"} (default) or \code{"off"}; see
#'   \code{\link{imputed_outcomes}}.
#' @param vartype \code{"bootstrap"} (default), \code{"jackknife"},
#'   \code{"parametric"}, or \code{"none"}. Selects which variance method
#'   to source replicates from. The output \code{vartype} column reports
#'   the method actually used at fit time (read from \code{fit$vartype}),
#'   which may differ from this argument value if the fit was produced
#'   with a different setting --- the argument is informational and does
#'   not re-aggregate replicates.
#' @param conf.level Two-sided confidence level. Defaults to 0.95.
#' @param ci.method One of \code{"basic"} (reflected),
#'   \code{"percentile"} (raw bootstrap quantiles), \code{"bc"}
#'   (bias-corrected percentile; Efron 1987 minus the acceleration),
#'   or \code{"normal"} (Wald: \eqn{\hat\theta \pm z \cdot SE}).
#'   Default is \code{NULL}, which triggers a per-type default:
#'   \code{"att"} -> \code{"normal"} (matches what \code{fit$est.att}
#'   already uses), \code{"att.cumu"} -> \code{"percentile"} (matches
#'   what \code{att.cumu()} does internally), \code{"aptt"} ->
#'   \code{"bc"} and \code{"log.att"} -> \code{"bc"} (ratio / log
#'   estimators benefit from bias correction when the bootstrap
#'   distribution is skewed). Pass an explicit value to override.
#'
#' @return A data frame with columns \code{<by_key>}, \code{estimate},
#'   \code{se}, \code{ci.lo}, \code{ci.hi}, \code{n_cells}, and
#'   \code{vartype}. Always tidy regardless of \code{type} or \code{by}.
#'
#' @section Numerical equality with existing slots:
#'   \code{estimand(fit, "att", "event.time")} returns
#'   \code{estimate}, \code{se}, \code{ci.lo}, \code{ci.hi}
#'   byte-identical to columns \code{ATT}, \code{S.E.}, \code{CI.lower},
#'   \code{CI.upper} of \code{fit$est.att}, when default arguments are
#'   used. This invariant is asserted by package tests.
#'
#' @seealso \code{\link{imputed_outcomes}} for the underlying long-form
#'   accessor; \code{\link{fect}} for the fitting interface.
#'
#' @examples
#' \dontrun{
#' library(fect)
#' fit <- fect(Y ~ D, data = simdata, index = c("id", "time"),
#'             method = "fe", force = "two-way",
#'             se = TRUE, nboots = 200, parallel = FALSE)
#'
#' ## Default: per-event-time ATT (matches fit$est.att numerically).
#' est <- estimand(fit, "att", "event.time")
#' head(est)
#' }
#'
#' @export
estimand <- function(fit,
                     type   = c("att", "att.cumu", "aptt", "log.att"),
                     by     = c("event.time", "cohort", "calendar.time",
                                "overall"),
                     test        = c("none", "placebo", "carryover"),
                     cells       = NULL,
                     weights     = NULL,
                     window      = NULL,
                     direction   = c("on", "off"),
                     vartype     = c("bootstrap", "jackknife", "parametric", "none"),
                     conf.level  = 0.95,
                     ci.method   = NULL) {

    type      <- match.arg(type)
    test      <- match.arg(test)
    vartype   <- match.arg(vartype)

    ## test = "placebo" / "carryover": auto-set direction to the semantic
    ## default. Users should not need to remember the pairing.
    if (test == "placebo") {
        direction <- "on"
    } else if (test == "carryover") {
        direction <- "off"
    } else {
        direction <- match.arg(direction)
    }

    ## Cumulative semantics are undefined for placebo / carryover cells.
    if (test != "none" && type == "att.cumu") {
        stop("estimand(type = \"att.cumu\") is incompatible with ",
             "test = \"", test, "\". Cumulative effects are defined ",
             "relative to treatment onset; placebo and carryover cells ",
             "do not have a meaningful cumulative anchor.",
             call. = FALSE)
    }

    ## ci.method = NULL triggers per-type defaults (v2.4.2+).
    ## See statsclaw-workspace/fect/ref/v242-vartype-cimethod-design.md.
    if (is.null(ci.method)) {
        ci.method <- switch(type,
            "att"      = "normal",      ## matches what fit$est.att uses (Wald: theta +- z*SE)
            "att.cumu" = "percentile",  ## matches att.cumu() internals
            "aptt"     = "bca",         ## ratio: bootstrap-bias + skew -> BCa (Efron 1987)
            "log.att"  = "bca"          ## log: same rationale
        )
    }
    ci.method <- match.arg(ci.method,
                           c("basic", "percentile", "bc", "bca", "normal"))

    by_canon <- c("event.time", "cohort", "calendar.time", "overall")
    if (length(by) > 1L) {
        by <- match.arg(by, by_canon)
    } else {
        if (!is.character(by) || length(by) != 1L) {
            stop("`by` must be a single string.", call. = FALSE)
        }
        if (!(by %in% by_canon)) {
            ## Future: resolve as user column on fit's panel data.
            stop("user-column `by` not yet supported in this release; ",
                 "use one of: ",
                 paste(shQuote(by_canon), collapse = ", "),
                 ".", call. = FALSE)
        }
    }

    if (!is.null(window)) {
        if (!is.null(cells)) {
            stop("Pass either `cells` or `window`, not both. ",
                 "`window = c(L, R)` is sugar for ",
                 "`cells = ~ event.time >= L & event.time <= R`.",
                 call. = FALSE)
        }
        if (!is.numeric(window) || length(window) != 2L) {
            stop("`window` must be a numeric vector of length 2: c(L, R).",
                 call. = FALSE)
        }
        L <- window[1]; R <- window[2]
        cells <- stats::as.formula(
            sprintf("~ event.time >= %s & event.time <= %s", L, R)
        )
    }

    .validate_po_contract(fit)

    if (type == "att") {
        return(.estimand_att(fit, by, cells, weights, direction,
                             vartype, conf.level, ci.method, test))
    }
    if (type == "att.cumu") {
        return(.estimand_att_cumu(fit, by, cells, weights, direction,
                                  vartype, conf.level, ci.method, window))
    }
    if (type == "aptt") {
        return(.estimand_aptt(fit, by, cells, weights, direction,
                              vartype, conf.level, ci.method, test))
    }
    if (type == "log.att") {
        return(.estimand_log_att(fit, by, cells, weights, direction,
                                 vartype, conf.level, ci.method, test))
    }

    stop("type = \"", type, "\" is part of the v2.4.0 API surface but ",
         "is not yet implemented at this commit. Stay tuned.",
         call. = FALSE)
}


## Internal: type = "att" dispatcher.
.estimand_att <- function(fit, by, cells, weights, direction,
                          vartype, conf.level, ci.method, test = "none") {

    ## Fast path: by = "event.time" + default args + direction = "on".
    ## Reads directly from fit$est.att for byte-equality with the existing
    ## fect output surface. Triggered when:
    ##   - by = "event.time"
    ##   - cells = NULL (no filter)
    ##   - weights = NULL (use fit's W.agg if any)
    ##   - direction = "on"
    ##   - conf.level = 0.95 (fit$est.att uses 95% by default)
    ##   - ci.method = "normal" (fit$est.att uses Wald: theta +- z*SE)
    ##                          --- this is the v2.4.2 default for type="att"
    is_fast_path <- by == "event.time" &&
                    is.null(cells) &&
                    is.null(weights) &&
                    direction == "on" &&
                    abs(conf.level - 0.95) < 1e-12 &&
                    ci.method == "normal" &&
                    test == "none"

    if (is_fast_path) {
        return(.estimand_att_fast_event_time(fit))
    }

    if (by == "event.time" && test != "none") {
        return(.estimand_att_event_time(fit, cells, weights, direction,
                                        vartype, conf.level, ci.method,
                                        test))
    }

    if (by == "overall") {
        return(.estimand_att_overall(fit, cells, weights, direction,
                                     vartype, conf.level, ci.method,
                                     test))
    }

    stop("estimand(type = \"att\") with by = \"", by, "\" is part of ",
         "the v2.4.0 API surface but is not yet implemented at this ",
         "commit. Default form estimand(fit, \"att\", \"event.time\") ",
         "and estimand(fit, \"att\", \"overall\", window = ...) work.",
         call. = FALSE)
}


## Per-event-time ATT slow path. Used when test = "placebo" /
## "carryover" forces per-cell aggregation.
.estimand_att_event_time <- function(fit, cells, weights, direction,
                                     vartype, conf.level, ci.method,
                                     test) {

    if (!is.null(weights)) {
        stop("estimand(\"att\", \"event.time\", test = \"", test, "\") ",
             "with non-default weights is not yet supported.",
             call. = FALSE)
    }
    if (!is.null(cells)) {
        stop("estimand(\"att\", \"event.time\", test = \"", test, "\") ",
             "with `cells` filter is not yet supported.",
             call. = FALSE)
    }

    mask_info <- .test_cell_mask(fit, test, direction)
    base_mask <- mask_info$mask
    Tev       <- mask_info$Tev

    ets <- sort(unique(Tev[base_mask]))
    if (length(ets) == 0) {
        stop("No cells satisfy test = \"", test, "\".", call. = FALSE)
    }

    nboots <- if (is.null(fit$eff.boot)) 0L else dim(fit$eff.boot)[3]

    estimate <- numeric(length(ets))
    se_vec   <- rep(NA_real_, length(ets))
    ci_lo    <- rep(NA_real_, length(ets))
    ci_hi    <- rep(NA_real_, length(ets))
    n_cells  <- integer(length(ets))

    for (k in seq_along(ets)) {
        et <- ets[k]
        cell_mask <- base_mask & Tev == et
        n_cells[k] <- sum(cell_mask)

        eff_t <- fit$eff[cell_mask]
        estimate[k] <- mean(eff_t, na.rm = TRUE)

        if (nboots > 0L && vartype != "none") {
            eff_boot_cells <- apply(fit$eff.boot, 3,
                                    function(eb) eb[cell_mask])
            if (!is.matrix(eff_boot_cells)) {
                eff_boot_cells <- matrix(eff_boot_cells,
                                         nrow = sum(cell_mask))
            }
            att_b <- colMeans(eff_boot_cells, na.rm = TRUE)

            jack_v <- if (ci.method == "bca") {
                .cell_jackknife("att", eff = eff_t)
            } else NULL

            ci <- .compute_ci(estimate[k], att_b, ci.method, conf.level,
                              jack = jack_v)
            se_vec[k] <- ci$se
            ci_lo[k]  <- ci$ci.lo
            ci_hi[k]  <- ci$ci.hi
        }
    }

    used_vartype <- if (vartype == "none") "none"
                    else if (is.null(fit$vartype)) "bootstrap"
                    else fit$vartype

    data.frame(
        event.time = ets,
        estimate   = estimate,
        se         = se_vec,
        ci.lo      = ci_lo,
        ci.hi      = ci_hi,
        n_cells    = n_cells,
        vartype    = used_vartype,
        stringsAsFactors = FALSE
    )
}


## Compute overall ATT (single scalar) over treated cells, optionally
## filtered. Reads from fit$eff and fit$D.dat directly; bootstrap from
## fit$eff.boot if available, else delegates to the pre-aggregated
## fit$att.avg.boot when no cells filter is active.
.estimand_att_overall <- function(fit, cells, weights, direction,
                                  vartype, conf.level, ci.method,
                                  test = "none") {

    if (!is.null(weights)) {
        stop("estimand(\"att\", \"overall\") with non-default weights ",
             "is not yet supported in v2.4.0.", call. = FALSE)
    }
    if (test != "none" && !is.null(cells)) {
        stop("estimand(\"att\", \"overall\") with both test = \"", test,
             "\" and `cells` is not supported. The test = ... argument ",
             "already filters cells to the placebo / carryover window.",
             call. = FALSE)
    }

    mask_info <- .test_cell_mask(fit, test, direction)
    treated_mask <- mask_info$mask
    Tev          <- mask_info$Tev

    ## Apply cells filter at the event-time / id level, not via long-form
    ## conversion (faster). Build a per-cell mask matching shape(eff).
    if (!is.null(cells)) {
        po <- imputed_outcomes(fit, direction = direction)
        po_filtered <- .apply_cells_filter(po, cells)
        ## Reconstruct the cell-level mask from the filtered po (id, time).
        keep_idx <- match(
            paste(po_filtered$id, po_filtered$time, sep = "::"),
            paste(fit$id[col(fit$D.dat)], fit$rawtime[row(fit$D.dat)],
                  sep = "::")
        )
        keep_mask <- matrix(FALSE, nrow = nrow(fit$D.dat),
                            ncol = ncol(fit$D.dat))
        keep_mask[keep_idx] <- TRUE
        cell_mask <- treated_mask & keep_mask
    } else {
        cell_mask <- treated_mask
    }

    n_cells <- sum(cell_mask)
    if (n_cells == 0L) {
        stop("No treated cells satisfy the requested filter.",
             call. = FALSE)
    }

    estimate <- mean(fit$eff[cell_mask], na.rm = TRUE)

    se_val <- NA_real_
    ci_lo  <- NA_real_
    ci_hi  <- NA_real_

    if (vartype != "none") {
        if (is.null(fit$eff.boot)) {
            stop("estimand(\"att\", \"overall\") with non-default ",
                 "filter requires keep.sims = TRUE in fect() so the ",
                 "per-cell bootstrap surface is available.",
                 call. = FALSE)
        }
        nboots <- dim(fit$eff.boot)[3]
        att_b  <- vapply(seq_len(nboots), function(b) {
            mean(fit$eff.boot[, , b][cell_mask], na.rm = TRUE)
        }, numeric(1))

        jack_v <- if (ci.method == "bca") {
            .cell_jackknife("att", eff = fit$eff[cell_mask])
        } else NULL

        ci <- .compute_ci(estimate, att_b, ci.method, conf.level,
                          jack = jack_v)
        se_val <- ci$se
        ci_lo  <- ci$ci.lo
        ci_hi  <- ci$ci.hi
    }

    used_vartype <- if (vartype == "none") "none"
                    else if (is.null(fit$vartype)) "bootstrap"
                    else fit$vartype

    data.frame(
        estimate = estimate,
        se       = se_val,
        ci.lo    = ci_lo,
        ci.hi    = ci_hi,
        n_cells  = as.integer(n_cells),
        vartype  = used_vartype,
        stringsAsFactors = FALSE
    )
}


## Fast path implementation: byte-equality with fit$est.att. No
## recomputation; just reshapes the existing matrix into the tidy
## return schema.
.estimand_att_fast_event_time <- function(fit) {

    if (is.null(fit$est.att)) {
        stop("fit$est.att is missing; refit with se = TRUE for SE/CI ",
             "or use vartype = \"none\".",
             call. = FALSE)
    }

    ## Detect vartype actually used at fit time. fect populates a
    ## `vartype` element on the fit object; default to "bootstrap" if
    ## absent.
    used_vartype <- if (is.null(fit$vartype)) "bootstrap" else fit$vartype

    out <- data.frame(
        event.time = as.numeric(rownames(fit$est.att)),
        estimate   = unname(fit$est.att[, "ATT"]),
        se         = unname(fit$est.att[, "S.E."]),
        ci.lo      = unname(fit$est.att[, "CI.lower"]),
        ci.hi      = unname(fit$est.att[, "CI.upper"]),
        n_cells    = unname(fit$est.att[, "count"]),
        vartype    = used_vartype,
        stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out
}


## Internal: type = "att.cumu" dispatcher.
.estimand_att_cumu <- function(fit, by, cells, weights, direction,
                               vartype, conf.level, ci.method, window) {

    if (direction != "on") {
        stop("estimand(\"att.cumu\") with direction = \"off\" is not ",
             "supported (cumulative effects are defined relative to ",
             "treatment onset).", call. = FALSE)
    }
    if (!is.null(weights)) {
        stop("estimand(\"att.cumu\") with non-default weights is not ",
             "yet supported in v2.4.0.", call. = FALSE)
    }
    if (!is.null(cells) && by != "overall") {
        stop("estimand(\"att.cumu\") with `cells` is supported only when ",
             "by = \"overall\". For event-time series, use ",
             "by = \"event.time\" with no cells filter.",
             call. = FALSE)
    }

    ## Implementation strategy: derive the cumulative bootstrap
    ## distribution from fit$att.boot (already pre-aggregated per event
    ## time) for the by = "event.time" and by = "overall" paths. These
    ## paths do not require keep.sims = TRUE. Cohort / calendar-time /
    ## user-column paths (deferred) need per-cell eff.boot.

    if (by == "event.time") {
        ## Event-time cumulative series: matches effect(fit)$catt.
        return(.compute_att_cumu_event_time(fit, conf.level, ci.method,
                                            vartype))
    }

    if (by == "overall") {
        ## Single scalar: matches the final row of att.cumu(fit, period).
        if (is.null(window)) {
            ## Default to the maximum post-treatment window.
            available_post <- fit$time[fit$time >= 1]
            if (length(available_post) == 0) {
                stop("No post-treatment event times in fit$time.",
                     call. = FALSE)
            }
            window <- c(min(available_post), max(available_post))
        }
        return(.compute_att_cumu_overall(fit, window, conf.level,
                                         ci.method, vartype))
    }

    stop("estimand(\"att.cumu\") with by = \"", by, "\" is not yet ",
         "supported in v2.4.0. Use by = \"event.time\" or ",
         "by = \"overall\".", call. = FALSE)
}


## Compute event-time cumulative ATT series. Delegates to the existing
## effect() function for the actual numerics (which use a per-cell mean
## approach, not a simple cumsum of fit$att). Reshapes the matrix into
## the tidy estimand return schema. Numerical equality with effect() is
## by construction.
.compute_att_cumu_event_time <- function(fit, conf.level, ci.method,
                                          vartype) {

    if (isTRUE(fit$hasRevs == 1)) {
        stop("Cumulative effects are not well-defined for panels with ",
             "treatment reversals.", call. = FALSE)
    }

    if (is.null(fit$eff.boot)) {
        stop("No bootstrap/jackknife results available. ",
             "Choose keep.sims = TRUE in fect().",
             call. = FALSE)
    }

    out_eff <- .with_no_deprecation(effect(fit, cumu = TRUE, plot = FALSE))
    if (is.null(out_eff) || is.null(out_eff$effect.est.att)) {
        stop("Cumulative effects unavailable for this fit.",
             call. = FALSE)
    }

    M <- out_eff$effect.est.att
    counts <- fit$est.att[match(rownames(M), rownames(fit$est.att)), "count"]
    data.frame(
        event.time = as.numeric(rownames(M)),
        estimate   = unname(M[, "ATT"]),
        se         = unname(M[, "S.E."]),
        ci.lo      = unname(M[, "CI.lower"]),
        ci.hi      = unname(M[, "CI.upper"]),
        n_cells    = as.integer(unname(counts)),
        vartype    = if (is.null(fit$vartype)) "bootstrap" else fit$vartype,
        stringsAsFactors = FALSE
    )
}


## Compute single overall cumulative ATT in a window. Delegates to the
## existing att.cumu() function for the canonical math (count-weighted
## across event times, percentile bootstrap CI). Returns the final row
## reshaped into the tidy estimand schema.
.compute_att_cumu_overall <- function(fit, window, conf.level, ci.method,
                                       vartype) {

    if (isTRUE(fit$hasRevs == 1)) {
        stop("Cumulative effects are not well-defined for panels with ",
             "treatment reversals.", call. = FALSE)
    }

    out_acc <- .with_no_deprecation(att.cumu(fit, period = window,
                                              plot = FALSE))
    final <- out_acc[nrow(out_acc), ]

    ## att.cumu() column name is "catt" when SE is present;
    ## "cumulative att" when no SE. Pick whichever is there.
    catt_col <- if ("catt" %in% names(final)) "catt" else "cumulative att"

    has_se <- "S.E." %in% names(final)

    et <- as.numeric(rownames(fit$est.att))
    in_window <- !is.na(et) & et >= window[1] & et <= window[2]
    n_cells_window <- sum(fit$est.att[in_window, "count"], na.rm = TRUE)

    data.frame(
        estimate = unname(final[catt_col]),
        se       = if (has_se) unname(final["S.E."])     else NA_real_,
        ci.lo    = if (has_se) unname(final["CI.lower"]) else NA_real_,
        ci.hi    = if (has_se) unname(final["CI.upper"]) else NA_real_,
        n_cells  = as.integer(n_cells_window),
        vartype  = if (is.null(fit$vartype)) "bootstrap" else fit$vartype,
        stringsAsFactors = FALSE
    )
}


## Internal: type = "aptt" dispatcher.
##
## APTT (Chen & Roth 2024 QJE):
##   APTT_g = mean_g(Y_obs - Y0_hat) / mean_g(Y0_hat)
## per group g (event time, cohort, etc.). Numerator and denominator are
## both computed per replicate and the ratio is taken inside each
## replicate, so the bootstrap distribution is the distribution of
## ratios, not the ratio of mean distributions.
.estimand_aptt <- function(fit, by, cells, weights, direction,
                           vartype, conf.level, ci.method, test = "none") {

    if (!is.null(weights)) {
        stop("estimand(\"aptt\") with non-default weights is not yet ",
             "supported in v2.4.0.", call. = FALSE)
    }
    if (!is.null(cells)) {
        stop("estimand(\"aptt\") with `cells` filter is not yet ",
             "supported in v2.4.0.", call. = FALSE)
    }
    if (is.null(fit$eff.boot) && vartype != "none") {
        stop("No bootstrap/jackknife results available. ",
             "Choose keep.sims = TRUE in fect().",
             call. = FALSE)
    }

    if (by == "event.time") {
        return(.compute_aptt_event_time(fit, conf.level, ci.method,
                                        vartype, direction, test))
    }

    stop("estimand(\"aptt\") with by = \"", by, "\" is not yet ",
         "supported in v2.4.0. Use by = \"event.time\".",
         call. = FALSE)
}


## Compute per-event-time APTT with bootstrap CI.
.compute_aptt_event_time <- function(fit, conf.level, ci.method, vartype,
                                      direction, test = "none") {

    mask_info <- .test_cell_mask(fit, test, direction)
    treated_mask <- mask_info$mask
    Tev          <- mask_info$Tev

    ets <- sort(unique(Tev[treated_mask]))
    if (length(ets) == 0) {
        stop("No cells satisfy test = \"", test, "\".", call. = FALSE)
    }

    nboots <- if (is.null(fit$eff.boot)) 0L else dim(fit$eff.boot)[3]

    estimate <- numeric(length(ets))
    se_vec   <- rep(NA_real_, length(ets))
    ci_lo    <- rep(NA_real_, length(ets))
    ci_hi    <- rep(NA_real_, length(ets))
    n_cells  <- integer(length(ets))

    alpha <- 1 - conf.level
    probs <- c(alpha / 2, 1 - alpha / 2)

    for (k in seq_along(ets)) {
        et <- ets[k]
        cell_mask <- treated_mask & Tev == et
        n_cells[k] <- sum(cell_mask)

        ## Point estimate.
        eff_t <- fit$eff[cell_mask]
        Y_t   <- fit$Y.dat[cell_mask]
        Y0_t  <- Y_t - eff_t
        num   <- mean(eff_t, na.rm = TRUE)
        den   <- mean(Y0_t, na.rm = TRUE)
        estimate[k] <- num / den

        ## Bootstrap distribution per replicate.
        if (nboots > 0L && vartype != "none") {
            ## eff_boot_cells: rows = cells in this group, cols = replicates.
            eff_boot_cells <- apply(fit$eff.boot, 3, function(eb) eb[cell_mask])
            if (!is.matrix(eff_boot_cells)) {
                eff_boot_cells <- matrix(eff_boot_cells,
                                         nrow = sum(cell_mask))
            }
            Y0_boot <- Y_t - eff_boot_cells

            ## Hard-error on cell-drop pathology (v2.4.2+).
            ## See .compute_log_att_event_time for the full rationale.
            ## For APTT specifically: when E(Y0_b) crosses zero in any
            ## replicate, the denominator blows up (or flips sign),
            ## producing wildly unstable per-replicate APTT values.
            mean_Y0_per_rep <- colMeans(Y0_boot, na.rm = TRUE)
            n_bad_reps <- sum(abs(mean_Y0_per_rep) < 1e-10 |
                              is.na(mean_Y0_per_rep))
            if (n_bad_reps > 0L) {
                pct_bad <- 100 * n_bad_reps / length(mean_Y0_per_rep)
                stop(sprintf(
                    "APTT bootstrap is unreliable at event time %s\n  (E(Y0_hat) at the point = %.4f, but %d of %d bootstrap replicates\n  (%.1f%%) have E(Y0_b) ~ 0, blowing up the APTT denominator and\n  producing wildly unstable per-replicate ratios).\n\n  Options:\n    1. Filter out cells where E(Y0_hat) is small relative to E(Y):\n         estimand(fit, \"aptt\", \"event.time\",\n                  cells = ~ abs(Y0_hat) > <threshold>)\n    2. Transform the outcome to keep Y0_hat away from zero\n    3. Use a different estimand: estimand(\"att\", ...) does not have\n       this denominator instability",
                    as.character(et), den, n_bad_reps,
                    length(mean_Y0_per_rep), pct_bad
                ), call. = FALSE)
            }

            aptt_b  <- colMeans(eff_boot_cells, na.rm = TRUE) /
                       mean_Y0_per_rep

            ## Cell-level jackknife for the BCa acceleration parameter.
            ## Only computed when bca is requested (cheap; no model refits).
            jack_v <- if (ci.method == "bca") {
                .cell_jackknife("aptt", eff = eff_t, Y0 = Y0_t)
            } else NULL

            ci <- .compute_ci(estimate[k], aptt_b, ci.method, conf.level,
                              jack = jack_v)
            se_vec[k] <- ci$se
            ci_lo[k]  <- ci$ci.lo
            ci_hi[k]  <- ci$ci.hi
        }
    }

    used_vartype <- if (vartype == "none") "none"
                    else if (is.null(fit$vartype)) "bootstrap"
                    else fit$vartype

    data.frame(
        event.time = ets,
        estimate   = estimate,
        se         = se_vec,
        ci.lo      = ci_lo,
        ci.hi      = ci_hi,
        n_cells    = n_cells,
        vartype    = used_vartype,
        stringsAsFactors = FALSE
    )
}


## Internal: type = "log.att" dispatcher.
##
## logATT_g = mean_g(log(Y_obs) - log(Y0_hat)) over treated cells.
## Hard-stops if any cell included in the aggregation has Y_obs <= 0 or
## Y0_hat <= 0 (log undefined). Caller must pre-transform the outcome
## (e.g. log(Y + c) inside fect) so that all imputed and observed
## outcomes are strictly positive before requesting log.att.
.estimand_log_att <- function(fit, by, cells, weights, direction,
                              vartype, conf.level, ci.method,
                              test = "none") {

    if (!is.null(weights)) {
        stop("estimand(\"log.att\") with non-default weights is not ",
             "yet supported in v2.4.0.", call. = FALSE)
    }
    if (!is.null(cells)) {
        stop("estimand(\"log.att\") with `cells` filter is not yet ",
             "supported in v2.4.0.", call. = FALSE)
    }
    if (is.null(fit$eff.boot) && vartype != "none") {
        stop("No bootstrap/jackknife results available. ",
             "Choose keep.sims = TRUE in fect().",
             call. = FALSE)
    }

    ## Hard-stop on Y <= 0 / Y0_hat <= 0 at the point-estimate level.
    ## log.att is mathematically undefined for non-positive imputed or
    ## observed outcomes; silently dropping cells contaminates both the
    ## point estimate and the bootstrap distribution. Force the caller
    ## to pre-transform (typical fix: log(Y + c) for some c > 0 chosen
    ## so all Y0_hat > 0; refit; then re-call estimand).
    mask_info  <- .test_cell_mask(fit, test, direction)
    cells_mask <- mask_info$mask
    Y_chk  <- fit$Y.dat[cells_mask]
    Y0_chk <- Y_chk - fit$eff[cells_mask]
    n_bad_Y  <- sum(!is.na(Y_chk)  & Y_chk  <= 0)
    n_bad_Y0 <- sum(!is.na(Y0_chk) & Y0_chk <= 0)
    if (n_bad_Y + n_bad_Y0 > 0L) {
        min_Y  <- if (n_bad_Y  > 0L) min(Y_chk[!is.na(Y_chk)],   na.rm = TRUE) else NA_real_
        min_Y0 <- if (n_bad_Y0 > 0L) min(Y0_chk[!is.na(Y0_chk)], na.rm = TRUE) else NA_real_
        stop(sprintf(
            "log.att requires Y > 0 and Y0_hat > 0 in all treated cells.\n  Found %d cell(s) with Y <= 0 (min Y = %s) and %d cell(s) with Y0_hat <= 0 (min Y0_hat = %s).\n  log(Y) and log(Y0_hat) are undefined; silent dropping would bias the point estimate.\n\n  Fix: refit fect on a strictly-positive outcome, e.g.\n      data$Y_pos <- data$Y + (abs(min(data$Y)) + 1)\n      fit <- fect(Y_pos ~ D, ...)\n      estimand(fit, \"log.att\", \"event.time\")\n  Then back out the original-scale interpretation as needed.",
            n_bad_Y,  if (is.na(min_Y))  "NA" else sprintf("%.4f", min_Y),
            n_bad_Y0, if (is.na(min_Y0)) "NA" else sprintf("%.4f", min_Y0)
        ), call. = FALSE)
    }

    if (by == "event.time") {
        return(.compute_log_att_event_time(fit, conf.level, ci.method,
                                           vartype, direction, test))
    }

    stop("estimand(\"log.att\") with by = \"", by, "\" is not yet ",
         "supported in v2.4.0. Use by = \"event.time\".",
         call. = FALSE)
}


## Compute per-event-time log-ATT. Drops cells where either Y_obs or
## Y0_hat is non-positive (would give -Inf or NaN under log).
.compute_log_att_event_time <- function(fit, conf.level, ci.method,
                                         vartype, direction,
                                         test = "none") {

    mask_info <- .test_cell_mask(fit, test, direction)
    treated_mask <- mask_info$mask
    Tev          <- mask_info$Tev

    ets <- sort(unique(Tev[treated_mask]))
    if (length(ets) == 0) {
        stop("No cells satisfy test = \"", test, "\".", call. = FALSE)
    }

    nboots <- if (is.null(fit$eff.boot)) 0L else dim(fit$eff.boot)[3]

    estimate <- numeric(length(ets))
    se_vec   <- rep(NA_real_, length(ets))
    ci_lo    <- rep(NA_real_, length(ets))
    ci_hi    <- rep(NA_real_, length(ets))
    n_cells  <- integer(length(ets))

    alpha <- 1 - conf.level
    probs <- c(alpha / 2, 1 - alpha / 2)

    for (k in seq_along(ets)) {
        et <- ets[k]
        cell_mask <- treated_mask & Tev == et

        eff_t <- fit$eff[cell_mask]
        Y_t   <- fit$Y.dat[cell_mask]
        Y0_t  <- Y_t - eff_t

        ## Caller-level hard-stop in .estimand_log_att already guarantees
        ## all cells have Y > 0 and Y0_hat > 0; only NA filtering needed.
        ok <- !is.na(Y_t) & !is.na(Y0_t)
        n_cells[k] <- sum(ok)

        if (sum(ok) == 0L) {
            estimate[k] <- NA_real_
            next
        }

        log_diff <- log(Y_t[ok]) - log(Y0_t[ok])
        estimate[k] <- mean(log_diff, na.rm = TRUE)

        if (nboots > 0L && vartype != "none") {
            eff_boot_cells <- apply(fit$eff.boot, 3,
                                    function(eb) eb[cell_mask])
            if (!is.matrix(eff_boot_cells)) {
                eff_boot_cells <- matrix(eff_boot_cells,
                                         nrow = sum(cell_mask))
            }
            eff_ok    <- eff_boot_cells[ok, , drop = FALSE]
            Y_ok      <- Y_t[ok]
            Y0_b_ok   <- Y_ok - eff_ok

            ## Hard-error on cell-drop pathology (v2.4.2+).
            ##
            ## When a cell used in the point estimate has Y0_b <= 0 in
            ## a non-trivial fraction of bootstrap replicates, log(Y0_b)
            ## returns NaN and colMeans(..., na.rm = TRUE) silently
            ## averages over fewer cells in that replicate, breaking
            ## the basic bootstrap principle and contaminating the
            ## bootstrap distribution.
            ##
            ## Threshold: trigger when the WORST cell has Y0_b <= 0 in
            ## > 5% of replicates. Sub-threshold cells are tolerated
            ## (small dropping is benign at the bootstrap-distribution
            ## scale; >5% indicates a genuinely unstable cell that
            ## needs filtering or a different estimand).
            n_reps_per_cell <- rowSums(Y0_b_ok <= 0, na.rm = TRUE)
            n_total_reps    <- ncol(Y0_b_ok)
            drop_frac       <- n_reps_per_cell / n_total_reps
            worst_idx       <- which.max(drop_frac)
            worst_frac      <- drop_frac[worst_idx]
            if (length(worst_frac) && worst_frac > 0.05) {
                worst_Y0 <- Y0_t[ok][worst_idx]
                stop(sprintf(
                    "log-ATT bootstrap is unreliable at event time %s.\n  The worst cell has Y0_hat = %.4f but %d of %d bootstrap replicates\n  (%.1f%%) have Y0_b <= 0 for it, so log(Y0_b) is undefined and the\n  per-replicate average silently drops the cell. This contaminates the\n  bootstrap distribution and yields meaningless inference.\n\n  Options:\n    1. Filter out unstable cells:\n         estimand(fit, \"log.att\", \"event.time\",\n                  cells = ~ Y0_hat > <threshold>)\n    2. Transform the outcome before fect: log(Y + c) for some c > 0\n    3. Use a different estimand: estimand(\"att\", ...) does not have\n       this pathology",
                    as.character(et), worst_Y0,
                    n_reps_per_cell[worst_idx], n_total_reps,
                    100 * worst_frac
                ), call. = FALSE)
            }

            log_Y0_b  <- log(Y0_b_ok)
            log_Y     <- log(Y_ok)
            log_diff_b <- log_Y - log_Y0_b
            logatt_b  <- colMeans(log_diff_b, na.rm = TRUE)

            ## Cell-level jackknife on the per-cell log-diff vector.
            jack_v <- if (ci.method == "bca") {
                .cell_jackknife("log.att", log_diff = log_diff)
            } else NULL

            ci <- .compute_ci(estimate[k], logatt_b, ci.method, conf.level,
                              jack = jack_v)
            se_vec[k] <- ci$se
            ci_lo[k]  <- ci$ci.lo
            ci_hi[k]  <- ci$ci.hi
        }
    }

    used_vartype <- if (vartype == "none") "none"
                    else if (is.null(fit$vartype)) "bootstrap"
                    else fit$vartype

    data.frame(
        event.time = ets,
        estimate   = estimate,
        se         = se_vec,
        ci.lo      = ci_lo,
        ci.hi      = ci_hi,
        n_cells    = n_cells,
        vartype    = used_vartype,
        stringsAsFactors = FALSE
    )
}


## ---------------------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------------------

## Compute SE and (lo, hi) confidence-interval bounds from a bootstrap
## distribution `boot` (a numeric vector of replicate-level estimates),
## a point `estimate`, and a chosen `ci.method`.
##
## Supports five methods:
##   - "basic":     ci = (2*est - q_hi, 2*est - q_lo) [reflected]
##   - "percentile": ci = (q_lo, q_hi)
##   - "bc":        bias-corrected percentile (z0 only, no acceleration)
##                  ci = (q_{Phi(2*z0 + z_alpha/2)},
##                        q_{Phi(2*z0 + z_{1-alpha/2})})
##                  where z0 = Phi^-1(mean(boot < est)).
##   - "bca":       bias-corrected accelerated (Efron 1987 full BCa).
##                  Requires `jack` (a numeric vector of leave-one-out
##                  point estimates over the cells in the aggregation
##                  group) to compute the acceleration parameter.
##                  Cutoffs: a_lo = Phi(z0 + (z0+z_alpha/2) /
##                                          (1 - a*(z0+z_alpha/2)))
##                  Handles bootstrap-bias + bootstrap-skew jointly;
##                  default for ratio (aptt) and log (log.att) estimands
##                  where the bootstrap distribution is inherently
##                  skewed and bc alone degenerates at the boundary.
##   - "normal":    ci = est +/- z_{1-alpha/2} * SE   [Wald, symmetric]
##
## `jack` is the per-cell leave-one-out vector of within-group point
## estimates; required for "bca", ignored otherwise. The acceleration
## is a = sum((mean(jack) - jack)^3) / (6 * (sum(...^2))^1.5).
##
## Returns a list with elements `se`, `ci.lo`, `ci.hi`.
.compute_ci <- function(estimate, boot, ci.method, conf.level, jack = NULL) {
    alpha <- 1 - conf.level
    se    <- stats::sd(boot, na.rm = TRUE)

    if (ci.method == "normal") {
        z <- stats::qnorm(1 - alpha / 2)
        return(list(se    = se,
                    ci.lo = estimate - z * se,
                    ci.hi = estimate + z * se))
    }

    probs <- c(alpha / 2, 1 - alpha / 2)
    qs    <- stats::quantile(boot, probs = probs, na.rm = TRUE)

    if (ci.method == "percentile") {
        return(list(se = se, ci.lo = unname(qs[1]), ci.hi = unname(qs[2])))
    }
    if (ci.method == "basic") {
        return(list(se    = se,
                    ci.lo = 2 * estimate - unname(qs[2]),
                    ci.hi = 2 * estimate - unname(qs[1])))
    }
    if (ci.method == "bc") {
        valid <- !is.na(boot)
        if (sum(valid) == 0L) {
            return(list(se = NA_real_, ci.lo = NA_real_, ci.hi = NA_real_))
        }
        ## z0 = bias correction = Phi^-1(P(boot < estimate))
        p_below <- mean(boot[valid] < estimate)
        ## Clamp to avoid +/-Inf at the boundaries
        p_below <- pmin(pmax(p_below, 1e-6), 1 - 1e-6)
        z0      <- stats::qnorm(p_below)
        z_lo    <- stats::qnorm(alpha / 2)
        z_hi    <- stats::qnorm(1 - alpha / 2)
        a_lo    <- stats::pnorm(2 * z0 + z_lo)
        a_hi    <- stats::pnorm(2 * z0 + z_hi)
        bc_qs   <- stats::quantile(boot, probs = c(a_lo, a_hi), na.rm = TRUE)
        return(list(se = se, ci.lo = unname(bc_qs[1]), ci.hi = unname(bc_qs[2])))
    }
    if (ci.method == "bca") {
        if (is.null(jack)) {
            stop("ci.method = \"bca\" requires the per-cell jackknife ",
                 "vector via the `jack` argument; the caller must compute ",
                 "leave-one-out within-group estimates and pass them.",
                 call. = FALSE)
        }
        valid <- !is.na(boot)
        jack_valid <- !is.na(jack)
        if (sum(valid) == 0L || sum(jack_valid) < 2L) {
            return(list(se = NA_real_, ci.lo = NA_real_, ci.hi = NA_real_))
        }
        ## z0: bias correction
        p_below <- mean(boot[valid] < estimate)
        p_below <- pmin(pmax(p_below, 1e-6), 1 - 1e-6)
        z0      <- stats::qnorm(p_below)
        ## a: acceleration via cell-level jackknife
        jack_v   <- jack[jack_valid]
        jack_bar <- mean(jack_v)
        dev      <- jack_bar - jack_v
        num      <- sum(dev^3)
        den      <- 6 * (sum(dev^2))^1.5
        a        <- if (den > 1e-12) num / den else 0
        ## BCa cutoffs: handles z0 -> +/- inf via the (1 - a*z) denominator
        z_lo  <- stats::qnorm(alpha / 2)
        z_hi  <- stats::qnorm(1 - alpha / 2)
        adjust <- function(z_q) {
            denom <- 1 - a * (z0 + z_q)
            ## Guard against a*(z0+z_q) -> 1 (denom -> 0); fall back to bc.
            if (abs(denom) < 1e-8) {
                return(stats::pnorm(2 * z0 + z_q))
            }
            stats::pnorm(z0 + (z0 + z_q) / denom)
        }
        a_lo  <- adjust(z_lo)
        a_hi  <- adjust(z_hi)
        bca_qs <- stats::quantile(boot, probs = c(a_lo, a_hi), na.rm = TRUE)
        return(list(se = se, ci.lo = unname(bca_qs[1]), ci.hi = unname(bca_qs[2])))
    }
    stop("Unknown ci.method = \"", ci.method, "\".", call. = FALSE)
}


## Compute the per-cell jackknife vector for a within-group functional T.
## For aptt: T(eff, Y0) = mean(eff) / mean(Y0). leave-one-out:
##   theta_jack[i] = mean(eff[-i]) / mean(Y0[-i])
## For log.att: T(Y, Y0) = mean(log(Y) - log(Y0)). leave-one-out:
##   theta_jack[i] = mean(log(Y[-i]) - log(Y0[-i]))
## For att (level): T = mean(eff). theta_jack[i] = mean(eff[-i]).
##
## Used for the BCa acceleration parameter without requiring model refits
## (the influence is computed at the aggregation step, holding the model
## fixed). This is the standard practice for BCa when leave-one-unit-out
## refits are too expensive.
.cell_jackknife <- function(type, ...) {
    args <- list(...)
    if (type == "aptt") {
        eff <- args$eff; Y0 <- args$Y0
        n   <- length(eff)
        if (n < 2L) return(rep(NA_real_, n))
        sum_eff <- sum(eff, na.rm = TRUE)
        sum_Y0  <- sum(Y0,  na.rm = TRUE)
        ## leave-one-out means: (sum - eff_i) / (n - 1)
        num <- (sum_eff - eff) / (n - 1)
        den <- (sum_Y0  - Y0)  / (n - 1)
        return(num / den)
    }
    if (type == "log.att") {
        ld <- args$log_diff
        n  <- length(ld)
        if (n < 2L) return(rep(NA_real_, n))
        ## leave-one-out mean of log_diff
        sum_ld <- sum(ld, na.rm = TRUE)
        return((sum_ld - ld) / (n - 1))
    }
    if (type == "att") {
        eff <- args$eff
        n   <- length(eff)
        if (n < 2L) return(rep(NA_real_, n))
        sum_eff <- sum(eff, na.rm = TRUE)
        return((sum_eff - eff) / (n - 1))
    }
    stop("Unknown jackknife type = \"", type, "\".", call. = FALSE)
}


## Build the cell-level base mask for a given test (none / placebo /
## carryover) and direction. Returns a list with `mask` (logical
## matrix matching shape(fit$D.dat)) and `Tev` (the relevant event-
## time matrix from fit$T.on or fit$T.off).
##
## test = "none":      treated post-treatment cells (the default ATT
##                     surface): D.dat == 1 with non-NA Tev.
## test = "placebo":   pre-treatment cells masked during the placebo
##                     fit, identified by Tev within fit$placebo.period.
##                     Requires fit$placeboTest == TRUE.
## test = "carryover": early post-reversal cells masked during the
##                     carryover fit, identified by Tev within
##                     fit$carryover.period (Tev = T.off). Requires
##                     fit$carryoverTest == TRUE and hasRevs.
##
## v2.4.3+ (closes issue #131, ajunquera).
.test_cell_mask <- function(fit, test, direction) {

    Tev <- if (direction == "on") fit$T.on else fit$T.off
    if (is.null(Tev)) {
        stop("direction = \"", direction, "\" requested, but fit$T.",
             direction, " is NULL.", call. = FALSE)
    }

    if (test == "none") {
        return(list(
            mask = !is.na(fit$D.dat) & fit$D.dat == 1 & !is.na(Tev),
            Tev  = Tev
        ))
    }

    if (test == "placebo") {
        if (!isTRUE(as.logical(fit$placeboTest)) ||
            is.null(fit$placebo.period)) {
            stop("test = \"placebo\" requires the fit to have been run ",
                 "with placeboTest = TRUE. The placebo estimand is only ",
                 "meaningful when the placebo cells were masked from the ",
                 "fit (out-of-sample predictions); a standard fit's ",
                 "pre-treatment residuals are in-sample and would not be ",
                 "an honest credibility check. Refit with: ",
                 "fect(..., placeboTest = TRUE, placebo.period = c(L, R)).",
                 call. = FALSE)
        }
        pp <- fit$placebo.period
        if (length(pp) == 1L) pp <- c(pp, pp)
        return(list(
            mask = !is.na(Tev) & Tev >= pp[1] & Tev <= pp[2],
            Tev  = Tev
        ))
    }

    if (test == "carryover") {
        if (!isTRUE(as.logical(fit$carryoverTest)) ||
            is.null(fit$carryover.period)) {
            stop("test = \"carryover\" requires the fit to have been run ",
                 "with carryoverTest = TRUE. The carryover estimand is ",
                 "only meaningful when the early post-reversal cells were ",
                 "masked from the fit (out-of-sample predictions). Refit ",
                 "with: fect(..., carryoverTest = TRUE, ",
                 "carryover.period = c(L, R)).",
                 call. = FALSE)
        }
        if (!isTRUE(fit$hasRevs == 1)) {
            stop("test = \"carryover\" requires a panel with treatment ",
                 "reversals (fit$hasRevs == 1).", call. = FALSE)
        }
        cp <- fit$carryover.period
        if (length(cp) == 1L) cp <- c(cp, cp)
        return(list(
            mask = !is.na(Tev) & Tev >= cp[1] & Tev <= cp[2],
            Tev  = Tev
        ))
    }

    stop("Unknown test = \"", test, "\".", call. = FALSE)
}


## Apply a `cells` filter (NULL, logical, or one-sided formula/function)
## against a long-form data frame. Returns the filtered data frame.
.apply_cells_filter <- function(df, cells) {
    if (is.null(cells)) return(df)

    if (is.logical(cells)) {
        if (length(cells) != nrow(df)) {
            stop(
                "cells: logical filter has length ", length(cells),
                " but data has ", nrow(df), " rows.",
                call. = FALSE
            )
        }
        return(df[cells & !is.na(cells), , drop = FALSE])
    }

    if (inherits(cells, "formula")) {
        if (length(cells) != 2L) {
            stop("cells: formula must be one-sided (e.g., `~ event.time > 0`).",
                 call. = FALSE)
        }
        mask <- eval(cells[[2]], envir = df, enclos = parent.frame())
        if (!is.logical(mask) || length(mask) != nrow(df)) {
            stop("cells: formula must evaluate to a logical vector of ",
                 "length nrow(data).", call. = FALSE)
        }
        return(df[mask & !is.na(mask), , drop = FALSE])
    }

    if (is.function(cells)) {
        mask <- cells(df)
        if (!is.logical(mask) || length(mask) != nrow(df)) {
            stop("cells: function must return a logical vector of ",
                 "length nrow(data).", call. = FALSE)
        }
        return(df[mask & !is.na(mask), , drop = FALSE])
    }

    stop("cells must be NULL, a logical vector, a one-sided formula, ",
         "or a function.", call. = FALSE)
}
