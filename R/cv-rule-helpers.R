###################################################################
## CV rule selection helpers (added v2.3.0)
##
## The classic fect rule is a 1% relative-tolerance check applied
## sequentially during the (r) loop:
##     update r.cv if (best - current) > 0.01 * best
##
## v2.3.0 introduces the 1-SE rule (Breiman, Friedman, Olshen & Stone 1984;
## Hastie, Tibshirani & Friedman 2009 §7.10): pick the smallest r whose mean
## CV criterion is within one fold-SE of the best. Default `cv.rule = "1se"`.
## The 1% rule is preserved as `cv.rule = "1pct"` for byte-identical
## reproducibility of pre-2.3.0 fits. A pure-min rule is also available as
## `cv.rule = "min"`.
##
## The fold SE is computed from per-fold criterion values:
##     SE_r = sd(score_per_fold_r) / sqrt(K)
##
## Per-fold scores are obtained by calling .score_residuals() on each fold's
## residuals separately (instead of pooling residuals first), then taking
## column-wise sd / sqrt(K).
###################################################################


## Compute per-fold criterion vectors and pooled criterion vector for a list
## of fold-result lists (each list item containing $resid, $time_idx, optional
## $obs_w). Returns:
##   $pooled : named vector of criteria computed from pooled residuals (legacy)
##   $folds  : K x ncriteria matrix of per-fold criterion values
##   $mean   : named vector of mean-across-folds criterion values
##   $se     : named vector of (sd / sqrt(K)) per criterion
##
## When K == 1, $se is set to 0 for all criteria (no fold dispersion).
.fect_cv_aggregate_folds <- function(fold_list,
                                     count.T.cv,
                                     use_weight = 0L,
                                     norm.para = NULL) {
    ## Pooled (legacy) score
    all_resid    <- unlist(lapply(fold_list, `[[`, "resid"))
    all_time_idx <- unlist(lapply(fold_list, `[[`, "time_idx"))
    all_obs_w    <- if (use_weight == 1L)
                        unlist(lapply(fold_list, `[[`, "obs_w"))
                    else NULL

    pooled <- .score_residuals(
        resid         = all_resid,
        obs_weights   = all_obs_w,
        time_index    = all_time_idx,
        count_weights = count.T.cv,
        norm.para     = norm.para
    )

    ## Per-fold scores
    per_fold <- lapply(fold_list, function(fl) {
        .score_residuals(
            resid         = fl$resid,
            obs_weights   = if (use_weight == 1L) fl$obs_w else NULL,
            time_index    = fl$time_idx,
            count_weights = count.T.cv,
            norm.para     = norm.para
        )
    })
    score_mat <- do.call(rbind, per_fold)
    K <- nrow(score_mat)

    means <- if (K >= 1L) colMeans(score_mat, na.rm = TRUE) else pooled
    ses <- if (K > 1L) {
        apply(score_mat, 2L, function(x) {
            x_finite <- x[is.finite(x)]
            if (length(x_finite) <= 1L) return(0.0)
            stats::sd(x_finite) / sqrt(length(x_finite))
        })
    } else {
        rep(0.0, length(pooled))
    }
    names(ses) <- names(pooled)

    list(pooled = pooled, folds = score_mat, mean = means, se = ses)
}


## Apply a CV rule to select the row index of CV.out matching the chosen
## hyper-parameter (r or lambda).
##
## means: numeric vector indexed by row (one per candidate hyper-param value).
##        NAs and Inf are treated as missing (excluded from selection).
## ses:   numeric vector of SE per row, OR NULL (used only when rule == "1se").
## rule:  one of "1se" (default), "min", "1pct".
##
## Returns the row index (integer). Uses the smallest index among ties
## consistent with bias-toward-parsimony.
.fect_apply_cv_rule <- function(means, ses = NULL, rule = c("1se", "min", "1pct")) {
    rule <- match.arg(rule)
    valid <- is.finite(means)
    if (!any(valid)) return(NA_integer_)

    if (rule == "min") {
        m_min <- min(means[valid])
        return(min(which(valid & means == m_min)))
    }

    if (rule == "1pct") {
        m_min <- min(means[valid])
        threshold <- m_min * 1.01
        return(min(which(valid & means <= threshold)))
    }

    ## "1se"
    i_min <- min(which(valid & means == min(means[valid])))
    se_min <- if (!is.null(ses) && is.finite(ses[i_min])) ses[i_min] else 0
    threshold <- means[i_min] + se_min
    min(which(valid & means <= threshold))
}


## Validate user-supplied cv.rule argument and return a single canonical string.
.fect_validate_cv_rule <- function(cv.rule) {
    if (is.null(cv.rule)) return("1se")
    if (!is.character(cv.rule) || length(cv.rule) != 1L) {
        stop("'cv.rule' must be a single character string: ",
             "'1se' (default), 'min', or '1pct'.")
    }
    if (!cv.rule %in% c("1se", "min", "1pct")) {
        stop("'cv.rule' must be one of '1se', 'min', '1pct'. Got: '",
             cv.rule, "'.")
    }
    cv.rule
}
