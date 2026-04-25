###################################################################
## Standard rolling-window cross-validation for rank selection
##
## Standalone user-facing wrapper. For each of `k` folds and each
## eligible control unit, samples a random anchor time `t*` and:
##   - holds out `cv.nobs` observations starting at `t*` (scored for MSPE)
##   - drops `cv.buffer` observations immediately before `t*` from training
##     (buffer, attenuates AR-leakage from training to holdout)
##   - drops ALL observations from `t* + cv.nobs` onward from training
##     (the rolling-window step --- training cannot see the future of
##     the held-out block)
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
## Supports method = "ife" (IFE-EM, time.component.from = "notyettreated")
## and method = "gsynth" (GSC, time.component.from = "nevertreated"). Both
## paths populate Y.ct.full at masked control positions, so MSPE scoring
## works uniformly.
###################################################################

#' Rolling-window cross-validation for rank selection
#'
#' Picks the number of factors `r` for an interactive-fixed-effects
#' model via standard rolling-window cross-validation. For each of `k`
#' folds and each eligible control unit, a random anchor time `t*` is
#' sampled. The fold's training set excludes
#' `t* - cv.buffer, ..., t* - 1` (buffer), `t*, ..., t* + cv.nobs - 1`
#' (the held-out, scored block), and `t* + cv.nobs, ..., last(t)`
#' (rolling-window future drop). MSPE is scored at the held-out block
#' only and averaged across folds.
#'
#' This is the standard time-series CV design (cf. `forecast::tsCV`,
#' `tidymodels::sliding_window`, `caret::createTimeSlices`) adapted to
#' panel data: each unit gets its own anchor per fold, drawn uniformly
#' from valid positions.
#'
#' @param formula A model formula, e.g. `Y ~ D + X1 + X2`.
#' @param data A long-format data frame.
#' @param index Length-2 character vector: `c("unit", "time")`.
#' @param method Estimator. One of `"ife"` (IFE-EM, internally
#'   `time.component.from = "notyettreated"`) or `"gsynth"` (GSC,
#'   internally `time.component.from = "nevertreated"`). Both paths
#'   populate `Y.ct.full` at masked control positions, so the CV
#'   behavior is identical across the two estimators.
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
#' @param k Number of folds. Each fold draws a fresh set of per-unit
#'   anchors; the per-r MSPE is averaged across folds and the SE used
#'   by the `"1se"` rule reflects fold-to-fold variability. Default 10
#'   (matches the default for the existing CV strategies).
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
#' @param ... Additional arguments forwarded to `fect()` (e.g. `covs`).
#'
#' @return List with components:
#'   - `r.cv`: chosen rank.
#'   - `cv.rule`: rule applied.
#'   - `mspe`: data.frame of per-r MSPE (averaged across folds), SE
#'     across folds, and held-out cell counts.
#'   - `mspe.per.fold`: r-by-k matrix of per-fold MSPE.
#'   - `k`, `cv.nobs`, `cv.buffer`: parameters used.
#'   - `n.units.masked`: distinct units that contributed to at least
#'     one fold's holdout.
#'
#' @examples
#' \dontrun{
#'   library(fect)
#'   data(simdata)
#'   res <- r.cv.rolling(Y ~ D, data = simdata, index = c("id", "time"),
#'                       method = "ife", r.max = 5,
#'                       cv.nobs = 3, cv.buffer = 1, k = 10)
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
                          method = c("ife", "gsynth"),
                          r.max = 5L,
                          cv.nobs = 3L,
                          cv.buffer = 1L,
                          k = 10L,
                          cv.rule = c("1se", "min", "1pct"),
                          min.T0 = 5L,
                          force = "unit",
                          seed = NULL,
                          verbose = TRUE,
                          ...) {
    cv.rule <- match.arg(cv.rule)
    method  <- match.arg(method)
    fect_method <- "ife"
    fect_tcf <- if (identical(method, "gsynth")) "nevertreated" else "notyettreated"

    r.max     <- as.integer(r.max);     if (r.max     < 0L) stop("r.max must be >= 0.")
    cv.nobs   <- as.integer(cv.nobs);   if (cv.nobs   < 1L) stop("cv.nobs must be >= 1.")
    cv.buffer <- as.integer(cv.buffer); if (cv.buffer < 0L) stop("cv.buffer must be >= 0.")
    k         <- as.integer(k);         if (k         < 1L) stop("k must be >= 1.")
    min.T0    <- as.integer(min.T0);    if (min.T0    < 1L) stop("min.T0 must be >= 1.")

    if (length(index) != 2L) stop("index must be a length-2 character vector.")
    if (!all(index %in% colnames(data))) {
        stop("index columns not found in data: ",
             paste(setdiff(index, colnames(data)), collapse = ", "))
    }
    Yname <- as.character(formula[[2L]])
    if (!Yname %in% colnames(data)) {
        stop("Outcome '", Yname, "' not found in data.")
    }

    ## Identify treated units (excluded from rolling holdout sampling).
    Dname <- NULL
    rhs_terms <- attr(stats::terms(formula), "term.labels")
    if (length(rhs_terms) >= 1L) Dname <- rhs_terms[1L]
    treated_units <- character(0)
    if (!is.null(Dname) && Dname %in% colnames(data)) {
        agg <- tapply(data[[Dname]], data[[index[1L]]], max, na.rm = TRUE)
        treated_units <- as.character(names(agg)[agg >= 1])
    }

    ## Build per-unit observed-time vectors for control units.
    ctrl_obs_times <- by(data, data[[index[1L]]], function(d) {
        sort(d[[index[2L]]])
    })
    ctrl_obs_times <- ctrl_obs_times[
        !names(ctrl_obs_times) %in% treated_units
    ]
    ## Eligible: at least min.T0 + cv.nobs observations.
    eligible <- vapply(ctrl_obs_times, length,
                       integer(1)) >= (min.T0 + cv.nobs)
    ctrl_obs_times <- ctrl_obs_times[eligible]
    if (length(ctrl_obs_times) == 0L) {
        stop("r.cv.rolling: no control units have enough observations ",
             "(need >= min.T0 + cv.nobs = ", min.T0 + cv.nobs, ").")
    }

    if (isTRUE(verbose)) {
        message(sprintf(
            "r.cv.rolling: %d eligible control units; k = %d folds; cv.nobs = %d, cv.buffer = %d.",
            length(ctrl_obs_times), k, cv.nobs, cv.buffer
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

        ## Sample an anchor for each eligible unit.
        anchor_records <- lapply(names(ctrl_obs_times), function(u) {
            obs_t <- ctrl_obs_times[[u]]
            n_obs <- length(obs_t)
            ## Anchor index in obs_t such that:
            ##   - >= min.T0 obs before  -> a_idx >= min.T0 + 1
            ##   - >= cv.nobs obs from a -> a_idx <= n_obs - cv.nobs + 1
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

            unit_ord <- as.integer(fit$id)
            year_ord <- as.integer(fit$rawtime)
            Y_pred <- numeric(length(score_rows))
            for (kk in seq_along(score_rows)) {
                rr <- score_rows[kk]
                ui <- match(data[[index[1L]]][rr], unit_ord)
                ti <- match(data[[index[2L]]][rr], year_ord)
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
         n.units.masked = length(units_seen))
}
