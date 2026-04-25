###################################################################
## Rolling (forward-only) cross-validation for rank selection
##
## Standalone user-facing wrapper. Does its own CV loop in user space:
## for each candidate r in 0..r.max, masks the LAST cv.nobs observations
## of each unit (forward-only), refits fect with the masked panel, and
## scores MSPE at the held-out positions using fit$Y.ct.full. Applies the
## requested cv.rule.
##
## Closes the forward-leakage channel that the random-anchor masking in
## the standard cv.method = "all_units" / "treated_units" leaves open at
## cv.donut = 0 / 1 under serially correlated residuals.
##
## Currently supports method = "ife" with time.component.from = "notyettreated"
## (the only fect path that populates Y.ct.full at masked control positions).
## A migration message is emitted for unsupported configurations.
###################################################################

#' Rolling (forward-only) CV for rank selection
#'
#' Picks the number of factors `r` for an interactive-fixed-effects model
#' via deterministic forward-only cross-validation. For each unit, the
#' LAST `cv.nobs` observations are held out; training uses everything else
#' (the unit's earlier observations + all other units' full data). Unlike
#' `fect`'s default `cv.method = "all_units"` (random contiguous-block
#' masking), this design closes the forward-leakage channel that AR-
#' correlated residuals exploit at `cv.donut = 0` / `1`, and tends to
#' recover much smaller `r` on panels with serially correlated errors.
#'
#' @param formula A model formula, e.g. `Y ~ D + X1 + X2`.
#' @param data A long-format data frame.
#' @param index Length-2 character vector: `c("unit", "time")`.
#' @param method Estimator. Currently only `"ife"` (with
#'   `time.component.from = "notyettreated"`) is supported. Other paths
#'   error with a migration message because they don't populate
#'   `Y.ct.full` at masked control positions.
#' @param r.max Largest candidate rank to evaluate. CV is run over
#'   `0:r.max`.
#' @param cv.nobs Number of trailing observations per unit to hold out
#'   (default 3).
#' @param cv.rule Rule for picking `r` from the MSPE curve: `"1se"`
#'   (default), `"min"`, or `"1pct"`.
#' @param min.T0 Minimum pre-treatment observations per unit. Units with
#'   fewer than `min.T0 + cv.nobs` total observations are not masked
#'   (mirrors fect's existing safety check).
#' @param force One of `"none"`, `"unit"`, `"time"`, `"two-way"`. Default
#'   `"unit"`.
#' @param verbose If TRUE (default), print per-r MSPE.
#' @param ... Additional arguments forwarded to `fect()` (e.g. `covs`).
#'
#' @return List with components: `r.cv` (chosen rank), `cv.rule` (rule
#'   used), `mspe` (data.frame of per-r MSPE + SE + n_holdout),
#'   `n.units.masked` (number of units contributing to the held-out set).
#'
#' @examples
#' \dontrun{
#'   library(fect)
#'   data(simdata1)
#'   res <- r.cv.rolling(Y ~ D, data = simdata1, index = c("id", "time"),
#'                       method = "ife", r.max = 5, cv.nobs = 3)
#'   res$r.cv
#'   ## then use the chosen r in a CV-disabled fit:
#'   fit <- fect(Y ~ D, data = simdata1, index = c("id", "time"),
#'               method = "ife", time.component.from = "notyettreated",
#'               CV = FALSE, r = res$r.cv, se = TRUE)
#' }
#'
#' @export
r.cv.rolling <- function(formula,
                          data,
                          index,
                          method = "ife",
                          r.max = 5L,
                          cv.nobs = 3L,
                          cv.rule = c("1se", "min", "1pct"),
                          min.T0 = 5L,
                          force = "unit",
                          verbose = TRUE,
                          ...) {
    cv.rule <- match.arg(cv.rule)
    if (!identical(method, "ife")) {
        stop("r.cv.rolling: only method = 'ife' (with time.component.from = ",
             "'notyettreated') is supported in this version. ",
             "method = '", method, "' is not yet supported because that path ",
             "does not populate Y.ct.full at masked control positions.")
    }
    r.max <- as.integer(r.max)
    if (r.max < 0L) stop("r.max must be >= 0.")
    cv.nobs <- as.integer(cv.nobs)
    if (cv.nobs < 1L) stop("cv.nobs must be >= 1.")

    ## Build wide-form panel for masking. Ordering: rows = sorted unique
    ## time, cols = sorted unique unit.
    if (length(index) != 2L) stop("index must be a length-2 character vector.")
    if (!all(index %in% colnames(data))) {
        stop("index columns not found in data: ",
             paste(setdiff(index, colnames(data)), collapse = ", "))
    }
    Yname <- as.character(formula[[2L]])
    if (!Yname %in% colnames(data)) {
        stop("Outcome '", Yname, "' not found in data.")
    }

    units <- sort(unique(data[[index[1L]]]))
    times <- sort(unique(data[[index[2L]]]))
    TT <- length(times); N <- length(units)
    II <- matrix(0L, TT, N)
    for (k in seq_len(nrow(data))) {
        ti <- match(data[[index[2L]]][k], times)
        ui <- match(data[[index[1L]]][k], units)
        if (!is.na(ti) && !is.na(ui)) II[ti, ui] <- 1L
    }

    ## Identify treated units to exclude from the rolling mask. We mask
    ## controls only because masking treated obs intersects with the
    ## treatment indicator and confounds the rank selection.
    Dname <- NULL
    rhs_terms <- attr(stats::terms(formula), "term.labels")
    if (length(rhs_terms) >= 1L) Dname <- rhs_terms[1L]
    treated_units <- character(0)
    if (!is.null(Dname) && Dname %in% colnames(data)) {
        agg <- tapply(data[[Dname]], data[[index[1L]]], max, na.rm = TRUE)
        treated_units <- as.character(names(agg)[agg >= 1])
    }
    is_control <- !(as.character(units) %in% treated_units)

    ## Build the rolling mask on controls only.
    II_for_mask <- II
    II_for_mask[, !is_control] <- 0L
    cv.id <- cv.sample.rolling(II_for_mask, NULL,
                                cv.count = cv.nobs,
                                cv.treat = FALSE,
                                min.T0 = min.T0)$cv.id
    if (length(cv.id) == 0L) {
        stop("r.cv.rolling: no control units have enough observations ",
             "(need > min.T0 + cv.nobs = ", min.T0 + cv.nobs, ").")
    }

    ## Translate cv.id (linear in vec(II)) back to (unit, time) keys for
    ## masking the long-form data frame.
    j_idx <- (cv.id - 1L) %/% TT + 1L
    t_idx <- (cv.id - 1L) %%  TT + 1L
    ho_units <- units[j_idx]
    ho_times <- times[t_idx]
    ho_keys  <- paste(ho_units, ho_times, sep = "_")
    n_units_masked <- length(unique(j_idx))

    if (isTRUE(verbose)) {
        message(sprintf("r.cv.rolling: %d held-out positions across %d units (cv.nobs = %d).",
                        length(cv.id), n_units_masked, cv.nobs))
    }

    data_keys <- paste(data[[index[1L]]], data[[index[2L]]], sep = "_")
    mask_rows <- which(data_keys %in% ho_keys)
    Y_obs <- data[[Yname]][mask_rows]

    ## Per-r CV loop.
    r.grid   <- 0L:r.max
    mspe_vec <- rep(NA_real_, length(r.grid))
    se_vec   <- rep(NA_real_, length(r.grid))
    n_vec    <- rep(0L,        length(r.grid))

    for (idx in seq_along(r.grid)) {
        r_try <- r.grid[idx]
        data_masked <- data
        data_masked[[Yname]][mask_rows] <- NA_real_
        fit <- tryCatch(
            suppressWarnings(suppressMessages(fect(
                formula = formula, data = data_masked, index = index,
                method = "ife", time.component.from = "notyettreated",
                force = force, CV = FALSE, r = r_try,
                min.T0 = min.T0, se = FALSE, na.rm = FALSE, ...))),
            error = function(e) {
                if (isTRUE(verbose)) message(sprintf("  r=%d: fit error: %s", r_try, conditionMessage(e)))
                NULL
            })
        if (is.null(fit) || is.null(fit$Y.ct.full)) next

        unit_ord <- as.integer(fit$id)
        year_ord <- as.integer(fit$rawtime)
        Y_pred <- numeric(length(mask_rows))
        for (k in seq_along(mask_rows)) {
            rr <- mask_rows[k]
            ui <- match(data[[index[1L]]][rr], unit_ord)
            ti <- match(data[[index[2L]]][rr], year_ord)
            Y_pred[k] <- if (!is.na(ui) && !is.na(ti)) fit$Y.ct.full[ti, ui] else NA_real_
        }
        e2 <- (Y_obs - Y_pred)^2
        e2 <- e2[is.finite(e2)]
        if (length(e2) == 0L) next
        mspe_vec[idx] <- mean(e2)
        se_vec[idx]   <- stats::sd(e2) / sqrt(length(e2))
        n_vec[idx]    <- length(e2)
        if (isTRUE(verbose)) {
            message(sprintf("  r=%d: MSPE=%.4g (n=%d)", r_try, mspe_vec[idx], length(e2)))
        }
    }

    chosen_idx <- .fect_apply_cv_rule(mspe_vec, ses = se_vec, rule = cv.rule)
    if (is.na(chosen_idx)) {
        stop("r.cv.rolling: all candidate fits failed.")
    }
    r_cv <- r.grid[chosen_idx]
    if (isTRUE(verbose)) {
        message(sprintf("r.cv.rolling: r.cv = %d (rule = %s).", r_cv, cv.rule))
    }

    list(r.cv          = as.integer(r_cv),
         cv.rule       = cv.rule,
         mspe          = data.frame(r = r.grid, mspe = mspe_vec,
                                    se = se_vec, n_holdout = n_vec),
         n.units.masked = n_units_masked)
}
