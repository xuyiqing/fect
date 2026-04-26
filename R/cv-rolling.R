###################################################################
## Standard rolling-window cross-validation for rank selection
##
## Standalone user-facing wrapper. For each of `k` folds:
##   - samples `cv.prop` fraction of eligible units (controls + treated
##     pre-treatment); only sampled units carry a mask. Other units
##     stay fully observed and contribute training data at every
##     period --- preserves factor identification at the masked tails.
##   - per sampled unit, samples a random anchor time `t*` and:
##       - holds out `cv.nobs` observations starting at `t*` (scored for MSPE)
##       - drops `cv.buffer` observations immediately before `t*` from training
##         (buffer, attenuates AR-leakage from training to holdout)
##       - drops ALL observations from `t* + cv.nobs` onward from training
##         (the rolling-window step --- training cannot see the future of
##         the held-out block)
## For treated units, eligible times are the pre-treatment observations
## only (post-treatment is the imputation target, never masked).
## Refits fect with the masked panel for each candidate `r`, scores
## MSPE at the held-out positions via `fit$Y.ct.full`, and aggregates
## across folds. Applies the requested `cv.rule` to pick `r.cv`.
##
## Closes the forward-leakage channel that the existing
## `cv.method = "all_units" / "treated_units"` (random contiguous-block
## masking) leaves open at `cv.donut = 0 / 1` under serially correlated
## residuals: under rolling-window CV, the train/test boundary on the
## future side is closed by construction.
##
## Supports method = "ife" (IFE-EM, time.component.from = "notyettreated"),
## method = "gsynth" (GSC, time.component.from = "nevertreated"), and
## method = "cfe" (Complex Fixed Effects, time.component.from =
## "notyettreated"). All three paths populate Y.ct.full at masked
## positions, so MSPE scoring works uniformly. CFE-specific arguments
## (Z, gamma, Q, Q.type, kappa, extra index columns) are forwarded to
## the inner fect() call via `...`. CFE rolling CV picks `r` only;
## the user holds non-`r` CFE components fixed at their spec.
###################################################################

#' Rolling-window cross-validation for rank selection
#'
#' Picks the number of factors `r` for an interactive-fixed-effects
#' model via standard rolling-window cross-validation. For each of `k`
#' folds, a fraction `cv.prop` of eligible units (controls plus
#' treated pre-treatment) is sampled; only sampled units carry a mask.
#' For each sampled unit, a random anchor time `t*` is drawn and the
#' fold's training set excludes `t* - cv.buffer, ..., t* - 1` (buffer),
#' `t*, ..., t* + cv.nobs - 1` (the held-out, scored block), and
#' `t* + cv.nobs, ..., end_of_eligible(t)` (rolling-window future drop;
#' for treated units, `end_of_eligible` is the cell strictly before
#' treatment onset). MSPE is scored at the held-out block only and
#' averaged across folds.
#'
#' Per-fold unit sampling is required: masking every eligible unit at
#' the same time leaves no donor data at the masked time points and
#' breaks factor identification. Sampling `cv.prop` of units per fold
#' keeps unsampled units fully observed at all periods.
#'
#' This is the standard time-series CV design (cf. `forecast::tsCV`,
#' `tidymodels::sliding_window`, `caret::createTimeSlices`) adapted to
#' panel data: each sampled unit gets its own anchor per fold, drawn
#' uniformly from valid positions.
#'
#' @param formula A model formula, e.g. `Y ~ D + X1 + X2`.
#' @param data A long-format data frame.
#' @param index Character vector identifying the panel structure. For
#'   `method = "ife"` and `method = "gsynth"`, length 2:
#'   `c("unit", "time")`. For `method = "cfe"`, length >= 2: the first
#'   two entries are `c("unit", "time")` and any additional entries are
#'   extra grouping fixed-effect columns forwarded to the inner
#'   `fect()` call's `index =` argument.
#' @param method Estimator. One of `"ife"` (IFE-EM, internally
#'   `time.component.from = "notyettreated"`), `"gsynth"` (GSC,
#'   internally `time.component.from = "nevertreated"`), or `"cfe"`
#'   (Complex Fixed Effects, internally
#'   `time.component.from = "notyettreated"`). All three paths populate
#'   `Y.ct.full` at masked positions, so MSPE scoring works uniformly.
#'   For CFE, rolling CV picks `r` only; CFE-specific arguments
#'   (`Z`, `gamma`, `Q`, `Q.type`, `kappa`, extra index columns) are
#'   forwarded via `...` and held fixed at their user-supplied values.
#' @param r.max Largest candidate rank to evaluate. CV is run over
#'   `0:r.max`.
#' @param cv.nobs Length of the held-out (scored) block per unit per
#'   fold. Default 3.
#' @param cv.buffer Number of observations immediately BEFORE the
#'   held-out block to drop from training (the past-side buffer that
#'   attenuates AR-leakage). Default 1. Analogous to `cv.donut` in the
#'   existing `cv.method = "all_units"` / `"treated_units"` strategies,
#'   but applied only on the past side: the future side is dropped by
#'   construction.
#' @param k Number of folds. Each fold draws a fresh sample of units
#'   and a fresh set of per-unit anchors; the per-r MSPE is averaged
#'   across folds and the SE used by the `"1se"` rule reflects
#'   fold-to-fold variability. Default 20 (matches the default for
#'   the existing CV strategies).
#' @param cv.prop Fraction of eligible units sampled per fold. Only
#'   sampled units receive a mask in that fold; the rest stay fully
#'   observed and contribute training data at every period. Default
#'   0.1 (paired with `k = 20` for ~2x coverage of every eligible
#'   unit across folds). Across `k` folds, every eligible unit lands
#'   in the holdout roughly `k * cv.prop` times in expectation. Must
#'   satisfy `0 < cv.prop <= 1`. On small panels (n_eligible < 30)
#'   consider raising further, since per-fold MSPE precision scales
#'   with `cv.prop * n_eligible * cv.nobs`.
#' @param cv.rule Rule for picking `r` from the MSPE curve: `"1se"`
#'   (default), `"min"`, or `"1pct"`.
#' @param min.T0 Minimum observations required strictly before the
#'   anchor. Sets the lower bound on valid anchor positions. Default 5.
#' @param force One of `"none"`, `"unit"`, `"time"`, `"two-way"`. Default
#'   `"unit"`.
#' @param seed Optional integer base seed; per-fold seeds derive from
#'   `seed + fold_id` for reproducibility. Default `NULL` (use the
#'   ambient RNG).
#' @param verbose If TRUE (default), print per-fold per-r MSPE.
#' @param ... Additional arguments forwarded to `fect()`. For
#'   `method = "cfe"`, the user holds CFE structural arguments (`Z`,
#'   `gamma`, `Q`, `Q.type`, `Q.bspline.degree`, `kappa`, etc.) fixed
#'   at their spec via `...`; rolling CV varies only `r`.
#'
#' @return List with components:
#'   - `r.cv`: chosen rank.
#'   - `cv.rule`: rule applied.
#'   - `mspe`: data.frame of per-r MSPE (averaged across folds), SE
#'     across folds, and held-out cell counts.
#'   - `mspe.per.fold`: r-by-k matrix of per-fold MSPE.
#'   - `k`, `cv.nobs`, `cv.buffer`, `cv.prop`: parameters used.
#'   - `n.units.masked`: distinct units that contributed to at least
#'     one fold's holdout.
#'
#' @examples
#' \dontrun{
#'   library(fect)
#'   data(simdata)
#'   res <- r.cv.rolling(Y ~ D, data = simdata, index = c("id", "time"),
#'                       method = "ife", r.max = 5,
#'                       cv.nobs = 3, cv.buffer = 1, k = 20)
#'   res$r.cv
#'   ## then use the chosen r in a CV-disabled fit:
#'   fit <- fect(Y ~ D, data = simdata, index = c("id", "time"),
#'               method = "ife", time.component.from = "notyettreated",
#'               CV = FALSE, r = res$r.cv, se = TRUE)
#' }
#'
#' @export
r.cv.rolling <- function(formula,
                          data,
                          index,
                          method = c("ife", "gsynth", "cfe"),
                          r.max = 5L,
                          cv.nobs = 3L,
                          cv.buffer = 1L,
                          k = 20L,
                          cv.prop = 0.1,
                          cv.rule = c("1se", "min", "1pct"),
                          min.T0 = 5L,
                          force = "unit",
                          seed = NULL,
                          verbose = TRUE,
                          ...) {
    cv.rule <- match.arg(cv.rule)
    method  <- match.arg(method)
    fect_method <- method
    ## CFE on the notyettreated path also populates Y.ct.full at masked
    ## cells, so MSPE scoring works uniformly. Only "gsynth" uses
    ## the nevertreated factor estimation sample.
    fect_tcf <- if (identical(method, "gsynth")) "nevertreated" else "notyettreated"

    r.max     <- as.integer(r.max);     if (r.max     < 0L) stop("r.max must be >= 0.")
    cv.nobs   <- as.integer(cv.nobs);   if (cv.nobs   < 1L) stop("cv.nobs must be >= 1.")
    cv.buffer <- as.integer(cv.buffer); if (cv.buffer < 0L) stop("cv.buffer must be >= 0.")
    k         <- as.integer(k);         if (k         < 1L) stop("k must be >= 1.")
    min.T0    <- as.integer(min.T0);    if (min.T0    < 1L) stop("min.T0 must be >= 1.")
    cv.prop   <- as.numeric(cv.prop)
    if (!is.finite(cv.prop) || cv.prop <= 0 || cv.prop > 1) {
        stop("cv.prop must satisfy 0 < cv.prop <= 1.")
    }

    ## CFE may carry extra index columns (additional grouping FEs); IFE
    ## and gsynth still require length-2 (unit, time).
    if (length(index) < 2L) {
        stop("index must be a character vector of length >= 2 ",
             "(c(unit, time) for ife/gsynth; c(unit, time, ...) for cfe).")
    }
    if (!all(index %in% colnames(data))) {
        stop("index columns not found in data: ",
             paste(setdiff(index, colnames(data)), collapse = ", "))
    }
    Yname <- as.character(formula[[2L]])
    if (!Yname %in% colnames(data)) {
        stop("Outcome '", Yname, "' not found in data.")
    }

    ## Identify treated units and their treatment-onset times.
    ## Treated units' eligible holdout window is the pre-treatment
    ## observed times only (post-treatment is the imputation target,
    ## never masked).
    Dname <- NULL
    rhs_terms <- attr(stats::terms(formula), "term.labels")
    if (length(rhs_terms) >= 1L) Dname <- rhs_terms[1L]
    treat_onset <- list()  # named: unit -> first treated time
    if (!is.null(Dname) && Dname %in% colnames(data)) {
        for (u in unique(as.character(data[[index[1L]]]))) {
            rows <- which(as.character(data[[index[1L]]]) == u)
            d_u <- data[[Dname]][rows]
            t_u <- data[[index[2L]]][rows]
            treated_rows <- which(d_u >= 1)
            if (length(treated_rows) > 0L) {
                treat_onset[[u]] <- min(t_u[treated_rows], na.rm = TRUE)
            }
        }
    }

    ## Build per-unit eligible-time vectors:
    ##   - controls: all observed times
    ##   - treated:  observed times strictly before treatment onset
    elig_obs_times <- by(data, data[[index[1L]]], function(d) {
        u <- as.character(d[[index[1L]]][1L])
        t_all <- sort(d[[index[2L]]])
        if (!is.null(treat_onset[[u]])) t_all[t_all < treat_onset[[u]]]
        else                            t_all
    })
    ## Eligible: at least min.T0 + cv.nobs observations.
    eligible <- vapply(elig_obs_times, length,
                       integer(1)) >= (min.T0 + cv.nobs)
    elig_obs_times <- elig_obs_times[eligible]
    if (length(elig_obs_times) == 0L) {
        stop("r.cv.rolling: no eligible units have enough pre-treatment ",
             "observations (need >= min.T0 + cv.nobs = ",
             min.T0 + cv.nobs, ").")
    }

    n_eligible_units <- length(elig_obs_times)
    n_sample_per_fold <- max(1L, as.integer(round(cv.prop * n_eligible_units)))
    if (isTRUE(verbose)) {
        message(sprintf(
            "r.cv.rolling: %d eligible units (controls + treated pre-treatment); k = %d folds; cv.nobs = %d, cv.buffer = %d, cv.prop = %.3g (-> %d units sampled per fold).",
            n_eligible_units, k, cv.nobs, cv.buffer, cv.prop, n_sample_per_fold
        ))
    }

    ## --------------------------------------------------------------------
    ## Per-fold per-r CV loop.
    ## --------------------------------------------------------------------
    r.grid <- 0L:r.max
    fold_mspe <- matrix(NA_real_, nrow = length(r.grid), ncol = k,
                        dimnames = list(paste0("r", r.grid),
                                        paste0("fold", seq_len(k))))
    fold_n    <- matrix(0L, nrow = length(r.grid), ncol = k,
                        dimnames = dimnames(fold_mspe))
    units_seen <- character(0)
    data_keys  <- paste(data[[index[1L]]], data[[index[2L]]], sep = "_")

    for (fold_id in seq_len(k)) {
        if (!is.null(seed)) set.seed(as.integer(seed) + fold_id)

        ## Sample cv.prop fraction of eligible units to mask in this fold.
        ## Unsampled units stay fully observed and contribute training
        ## data at every period (preserves factor identification at the
        ## masked time points).
        sampled_idx <- sort(sample.int(n_eligible_units, n_sample_per_fold))
        sampled_units <- names(elig_obs_times)[sampled_idx]

        ## Sample an anchor for each sampled unit.
        anchor_records <- lapply(sampled_units, function(u) {
            obs_t <- elig_obs_times[[u]]
            n_obs <- length(obs_t)
            ## Anchor index in obs_t such that:
            ##   - >= min.T0 obs before  -> a_idx >= min.T0 + 1
            ##   - >= cv.nobs obs from a -> a_idx <= n_obs - cv.nobs + 1
            ## For treated units, obs_t already excludes post-treatment
            ## cells, so the holdout block is guaranteed to stay in the
            ## pre-treatment window.
            valid <- seq.int(min.T0 + 1L, n_obs - cv.nobs + 1L)
            if (length(valid) == 0L) return(NULL)
            a_idx <- valid[sample.int(length(valid), 1L)]
            holdout_t <- obs_t[a_idx:(a_idx + cv.nobs - 1L)]
            buf_t <- if (cv.buffer > 0L) {
                bs <- max(1L, a_idx - cv.buffer)
                obs_t[bs:(a_idx - 1L)]
            } else integer(0)
            drop_t <- if ((a_idx + cv.nobs) <= n_obs) {
                obs_t[(a_idx + cv.nobs):n_obs]
            } else integer(0)
            list(unit = u,
                 holdout_t = holdout_t,
                 buf_t = buf_t,
                 drop_t = drop_t)
        })
        anchor_records <- anchor_records[!vapply(anchor_records, is.null,
                                                 logical(1))]
        if (length(anchor_records) == 0L) {
            if (isTRUE(verbose)) message(sprintf("  fold %d: no valid anchors; skipped.", fold_id))
            next
        }

        ## Build mask key sets.
        scored_keys <- unlist(lapply(anchor_records, function(r) {
            paste(r$unit, r$holdout_t, sep = "_")
        }))
        buffer_keys <- unlist(lapply(anchor_records, function(r) {
            if (length(r$buf_t) > 0L) paste(r$unit, r$buf_t, sep = "_") else character(0)
        }))
        dropped_keys <- unlist(lapply(anchor_records, function(r) {
            if (length(r$drop_t) > 0L) paste(r$unit, r$drop_t, sep = "_") else character(0)
        }))
        all_mask_keys <- c(scored_keys, buffer_keys, dropped_keys)
        units_seen <- union(units_seen, vapply(anchor_records,
                                               function(r) r$unit,
                                               character(1)))

        ## Apply mask: set Y to NA at all masked cells.
        data_masked <- data
        data_masked[[Yname]][data_keys %in% all_mask_keys] <- NA_real_

        ## Score positions and observed Y.
        score_rows <- which(data_keys %in% scored_keys)
        Y_obs_score <- data[[Yname]][score_rows]

        ## Per-r fits for this fold.
        for (idx in seq_along(r.grid)) {
            r_try <- r.grid[idx]
            fit <- tryCatch(
                suppressWarnings(suppressMessages(fect(
                    formula = formula, data = data_masked, index = index,
                    method = fect_method, time.component.from = fect_tcf,
                    force = force, CV = FALSE, r = r_try,
                    min.T0 = min.T0, se = FALSE, na.rm = FALSE, ...))),
                error = function(e) {
                    if (isTRUE(verbose)) message(sprintf(
                        "  fold %d, r=%d: fit error: %s",
                        fold_id, r_try, conditionMessage(e)))
                    NULL
                })
            if (is.null(fit) || is.null(fit$Y.ct.full)) next

            ## Match by raw labels: fit$id and fit$rawtime are the unit/time
            ## labels in fect's coordinate system (integer for ife/gsynth, can
            ## be character/integer for cfe). Match against the same column
            ## values from `data` and accept whichever type comes back.
            unit_ord <- fit$id
            year_ord <- fit$rawtime
            unit_ord_int <- suppressWarnings(as.integer(unit_ord))
            year_ord_int <- suppressWarnings(as.integer(year_ord))
            unit_match <- if (any(is.na(unit_ord_int))) unit_ord else unit_ord_int
            year_match <- if (any(is.na(year_ord_int))) year_ord else year_ord_int
            Y_pred <- numeric(length(score_rows))
            for (kk in seq_along(score_rows)) {
                rr <- score_rows[kk]
                ui <- match(data[[index[1L]]][rr], unit_match)
                ti <- match(data[[index[2L]]][rr], year_match)
                Y_pred[kk] <- if (!is.na(ui) && !is.na(ti)) {
                    fit$Y.ct.full[ti, ui]
                } else NA_real_
            }
            e2 <- (Y_obs_score - Y_pred)^2
            e2 <- e2[is.finite(e2)]
            if (length(e2) == 0L) next
            fold_mspe[idx, fold_id] <- mean(e2)
            fold_n[idx, fold_id]    <- length(e2)
            if (isTRUE(verbose)) {
                message(sprintf("  fold %d, r=%d: MSPE=%.4g (n=%d)",
                                fold_id, r_try, fold_mspe[idx, fold_id],
                                length(e2)))
            }
        }
    }

    ## Aggregate per-r across folds.
    mspe_per_r <- rowMeans(fold_mspe, na.rm = TRUE)
    n_folds_used <- rowSums(!is.na(fold_mspe))
    se_per_r <- vapply(seq_along(r.grid), function(i) {
        x <- fold_mspe[i, ]
        x <- x[is.finite(x)]
        if (length(x) <= 1L) NA_real_
        else stats::sd(x) / sqrt(length(x))
    }, numeric(1))
    n_per_r <- rowSums(fold_n)

    chosen_idx <- .fect_apply_cv_rule(mspe_per_r, ses = se_per_r,
                                      rule = cv.rule)
    if (is.na(chosen_idx)) {
        stop("r.cv.rolling: all candidate fits failed across all folds.")
    }
    r_cv <- r.grid[chosen_idx]
    if (isTRUE(verbose)) {
        message(sprintf(
            "r.cv.rolling: r.cv = %d (rule = %s; aggregated over %d folds).",
            r_cv, cv.rule, k
        ))
    }

    list(r.cv          = as.integer(r_cv),
         cv.rule       = cv.rule,
         mspe          = data.frame(r = r.grid,
                                    mspe = mspe_per_r,
                                    se = se_per_r,
                                    n_holdout = n_per_r,
                                    n_folds_used = n_folds_used),
         mspe.per.fold = fold_mspe,
         k             = k,
         cv.nobs       = cv.nobs,
         cv.buffer     = cv.buffer,
         cv.prop       = cv.prop,
         n.units.masked = length(units_seen))
}
