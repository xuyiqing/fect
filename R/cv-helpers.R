###################################################################
## CV Helper Functions — single canonical implementation of one
## (r, fold) or (lambda, fold) scoring iteration for each call site.
##
## All functions are internal (prefix ".") and not exported.
## Used by cv.R (notyettreated path) and fect_nevertreated.R (nevertreated path).
###################################################################

## Threshold (Nco * TT) above which CV parallelism auto-enables for each method.
## Used when parallel = TRUE (auto mode); overridden when parallel = "cv".
## CFE value is a placeholder — empirical tuning is deferred to Phase 3.
.CV_PARALLEL_THRESH <- list(ife = 20000L, mc = 20000L, cfe = 60000L)


## Build a PSOCK cluster that inherits the parent's `.libPaths()` so that
## workers can find `fect` (and other user-installed packages) regardless
## of how the parent R session was launched (e.g. Quarto render, RStudio,
## bare Rscript). The default `future::plan(future::multisession)` does
## not propagate libPaths, which causes "no package called 'fect'" worker
## init errors in non-default R session contexts. Use this helper anywhere
## the package would otherwise call `future::plan(future::multisession,
## workers = cores)`. Returns the cluster object; caller is responsible
## for calling `future::plan(future::cluster, workers = cl)` and the
## subsequent on.exit cleanup.
.fect_make_future_cluster <- function(cores) {
    cl <- parallelly::makeClusterPSOCK(
        workers      = cores,
        rscript_libs = .libPaths(),
        autoStop     = TRUE
    )
    ## Pre-load packages on workers with messages + warnings suppressed
    ## (v2.4.2+). Without this, each worker fires "package was built under
    ## R version X.Y.Z" warnings from the user's local R / package skew on
    ## first use of mvtnorm / future / etc., once per worker. These are
    ## informational, not actionable, but they clutter user-facing output
    ## (especially under parametric bootstrap where mvtnorm::rmvnorm
    ## fires per worker).
    parallel::clusterEvalQ(cl, {
        suppressPackageStartupMessages({
            suppressWarnings({
                ## Pre-load packages used inside parallel paths.
                requireNamespace("mvtnorm",       quietly = TRUE)
                requireNamespace("future",        quietly = TRUE)
                requireNamespace("future.apply",  quietly = TRUE)
                requireNamespace("doParallel",    quietly = TRUE)
                requireNamespace("foreach",       quietly = TRUE)
            })
        })
        ## Suppress "built under R version" warnings for any subsequent
        ## library() / requireNamespace() inside this worker.
        options(warn.conflicts = FALSE)
    })
    cl
}


## Run `expr` with "package was built under R version" warnings
## suppressed (they fire from parallel workers loading packages
## compiled against a slightly newer R; informational, not actionable).
## Other warnings pass through normally. Use to wrap parallel-execution
## entry points (future_lapply / foreach calls) so user-facing output
## isn't cluttered. Added v2.4.2.
.fect_with_quiet_pkg_warnings <- function(expr) {
    withCallingHandlers(
        expr,
        warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("was built under R version", msg, fixed = TRUE)) {
                invokeRestart("muffleWarning")
            }
        }
    )
}


## Normalize the cv.method argument across all CV entry points
## (fect_cv, fect_binary_cv, fect_nevertreated, fect_mspe).
##
## Adds `"block"` as a forward-looking alias for `"all_units"` (added in
## v2.3.0; the names will be unified in v2.4.0 alongside a new `cv.units`
## parameter). Emits a one-time deprecation `message()` when the user
## passes the legacy `"all_units"` or `"treated_units"` strings, then
## maps `"block"` -> `"all_units"` so all internal code paths stay
## unchanged in this release.
##
## Returns the canonical (internal) cv.method string. `loo` and
## `rolling` pass through unchanged.
.fect_normalize_cv_method <- function(cv.method,
                                      allowed = c("rolling", "block",
                                                  "all_units",
                                                  "treated_units",
                                                  "loo")) {
    cv.method <- match.arg(cv.method, allowed)
    if (cv.method == "all_units") {
        message('Note: cv.method = "all_units" is being deprecated; ',
                'use cv.method = "block" instead. ',
                'Both names produce identical behavior in this release.')
    }
    if (cv.method == "block") {
        cv.method <- "all_units"  # internal dispatch alias
    }
    if (cv.method == "treated_units") {
        message('Note: cv.method = "treated_units" will be deprecated ',
                'when the cv.units parameter is added (v2.4.0). ',
                'The legacy name still works for now.')
    }
    cv.method
}


## --------------------------------------------------------------------------
## 2a. IFE all_units variant (notyettreated context)
##     Called from cv.R IFE loop (both serial and parallel branches)
## --------------------------------------------------------------------------
.fect_cv_score_one_ife_all <- function(
    ii,            # fold index (integer, 1..k)
    r,             # candidate rank (integer)
    YY,            # TT×N outcome matrix (masked-zero version)
    Y0CV,          # TT×N×k array of initial Y0 per fold
    X,             # TT×N×p covariate array
    II,            # TT×N observation indicator
    W.use,         # TT×N weight matrix or scalar 0
    WW,            # TT×N weight matrix (unmasked) or NULL
    beta0CV,       # p×1×k or 1×0×k array of initial beta0 per fold
    rmCV,          # list of k integer vectors (masked indices)
    estCV,         # list of k integer vectors (scored indices)
    T.on,          # TT×N matrix of event-time indicators
    force,         # integer: force FE type
    cv_tol,        # numeric tolerance for CV fits
    max.iteration, # integer
    use_weight     # integer: 0 or 1
) {
    II.cv <- II
    II.cv[rmCV[[ii]]] <- 0
    YY.cv <- YY
    YY.cv[rmCV[[ii]]] <- 0
    if (use_weight) {
        W.use2 <- W.use
        W.use2[rmCV[[ii]]] <- 0
    } else {
        W.use2 <- as.matrix(0)
    }
    est.cv.fit <- inter_fe_ub(
        YY.cv, as.matrix(Y0CV[, , ii]), X, II.cv,
        W.use2, as.matrix(beta0CV[, , ii]),
        r, force, cv_tol, max.iteration
    )$fit
    resid_ii <- YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]]
    idx_ii <- as.character(T.on[estCV[[ii]]])
    idx_ii[which(is.na(idx_ii))] <- "Control"
    obs_w_ii <- if (use_weight == 1) WW[estCV[[ii]]] else NULL
    list(resid = resid_ii, time_idx = idx_ii, obs_w = obs_w_ii)
}


## --------------------------------------------------------------------------
## 2b. MC all_units variant (notyettreated context)
##     Helper defined for symmetry; parallel branch wired in Phase 2.
## --------------------------------------------------------------------------
.fect_cv_score_one_mc_all <- function(
    ii,
    lambda_i,
    YY,
    Y0CV,
    X,
    II,
    W.use,
    WW,
    beta0CV,
    rmCV,
    estCV,
    T.on,
    force,
    cv_tol,
    max.iteration,
    use_weight
) {
    II.cv <- II; II.cv[rmCV[[ii]]] <- 0
    YY.cv <- YY; YY.cv[rmCV[[ii]]] <- 0
    W.use2 <- if (use_weight) {
        W2 <- W.use; W2[rmCV[[ii]]] <- 0; W2
    } else {
        as.matrix(0)
    }
    est.cv.fit <- inter_fe_mc(
        YY.cv, as.matrix(Y0CV[, , ii]), X, II.cv,
        W.use2, as.matrix(beta0CV[, , ii]),
        1L, lambda_i, force, cv_tol, max.iteration
    )$fit
    resid_ii <- YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]]
    idx_ii <- as.character(T.on[estCV[[ii]]]); idx_ii[is.na(idx_ii)] <- "Control"
    obs_w_ii <- if (use_weight == 1) WW[estCV[[ii]]] else NULL
    list(resid = resid_ii, time_idx = idx_ii, obs_w = obs_w_ii)
}


## --------------------------------------------------------------------------
## 2c. nevertreated IFE all_units variant
##     Replaces the inline body in fect_nevertreated.R IFE all_units foreach block
## --------------------------------------------------------------------------
.fect_cv_score_one_ife_nt_all <- function(
    ii,
    YY.co,
    Y0CV.co,
    X.co,
    II.co,
    W.use,
    W,
    beta0CV.co,
    rmCV,
    estCV,
    r,
    force,
    cv_tol,
    max.iteration
) {
    II.co.cv <- II.co; II.co.cv[rmCV[[ii]]] <- 0
    YY.co.cv <- YY.co; YY.co.cv[rmCV[[ii]]] <- 0
    W.cv <- if (!is.null(W)) {
        W2 <- W.use; W2[rmCV[[ii]]] <- 0; W2
    } else {
        as.matrix(0)
    }
    est.cv.fit <- inter_fe_ub(
        YY.co.cv, as.matrix(Y0CV.co[, , ii]), X.co, II.co.cv,
        W.cv, as.matrix(beta0CV.co[, , ii]),
        r, force, cv_tol, max.iteration
    )$fit
    resid_ii <- YY.co[estCV[[ii]]] - est.cv.fit[estCV[[ii]]]
    list(
        resid    = resid_ii,
        time_idx = rep("Control", length(resid_ii)),
        obs_w    = if (!is.null(W)) W.use[estCV[[ii]]] else NULL
    )
}


## --------------------------------------------------------------------------
## 2d. nevertreated IFE treated_units variant
##     Uses precomputed F.hat, U.tr, pre from the outer r-loop.
## --------------------------------------------------------------------------
.fect_cv_score_one_ife_nt_tr <- function(
    ii,
    U.tr,
    F.hat,
    pre,
    r,
    force,
    rmCV.tr,
    estCV.tr,
    W.tr,
    T.on,
    tr,
    TT,
    Ntr
) {
    pre.cv <- pre; pre.cv[rmCV.tr[[ii]]] <- 0
    if (r == 0) {
        if (force %in% c(1, 3)) {
            alpha.cv <- rep(0, Ntr)
            for (j in 1:Ntr) {
                pre_t <- which(pre.cv[, j] == 1)
                if (length(pre_t) > 0) alpha.cv[j] <- mean(U.tr[pre_t, j])
            }
            fitted.cv <- matrix(alpha.cv, TT, Ntr, byrow = TRUE)
        } else {
            fitted.cv <- matrix(0, TT, Ntr)
        }
    } else {
        fitted.cv <- matrix(0, TT, Ntr)
        for (j in 1:Ntr) {
            pre_t <- which(pre.cv[, j] == 1)
            if (length(pre_t) >= ncol(F.hat)) {
                F.pre <- as.matrix(F.hat[pre_t, , drop = FALSE])
                lam <- try(
                    solve(t(F.pre) %*% F.pre) %*% t(F.pre) %*% U.tr[pre_t, j],
                    silent = TRUE
                )
                if (!"try-error" %in% class(lam)) fitted.cv[, j] <- F.hat %*% lam
            }
        }
    }
    resid_ii <- U.tr[estCV.tr[[ii]]] - fitted.cv[estCV.tr[[ii]]]
    obs_w_ii <- if (!is.null(W.tr)) W.tr[estCV.tr[[ii]]] else NULL
    idx_ii <- NULL
    if (!is.null(T.on)) {
        T.on.tr <- T.on[, tr, drop = FALSE]
        idx_ii <- as.character(T.on.tr[estCV.tr[[ii]]])
        idx_ii[is.na(idx_ii)] <- "Control"
    }
    list(resid = resid_ii, time_idx = idx_ii, obs_w = obs_w_ii)
}


## --------------------------------------------------------------------------
## 2e. nevertreated CFE all_units variant
##     Same structure as 2c but calls complex_fe_ub with full CFE argument list
## --------------------------------------------------------------------------
.fect_cv_score_one_cfe_nt_all <- function(
    ii,
    YY.co,
    Y0CV.co,
    X.co,
    II.co,
    W.use,
    W,
    beta0CV.co,
    X.extra.FE.co.B,
    X.Z.co,
    X.Q.co,
    X.gamma.co,
    X.kappa.co,
    Zgamma.id,
    kappaQ.id,
    rmCV,
    estCV,
    r,
    force,
    cv_tol,
    max.iteration
) {
    II.co.cv <- II.co; II.co.cv[rmCV[[ii]]] <- 0
    YY.co.cv <- YY.co; YY.co.cv[rmCV[[ii]]] <- 0
    W.cv <- if (!is.null(W)) {
        W2 <- W.use; W2[rmCV[[ii]]] <- 0; W2
    } else {
        as.matrix(0)
    }
    est.cv.co <- complex_fe_ub(
        YY.co.cv, as.matrix(Y0CV.co[, , ii]), X.co,
        X.extra.FE.co.B, X.Z.co, X.Q.co, X.gamma.co, X.kappa.co,
        Zgamma.id, kappaQ.id,
        II.co.cv, W.cv, as.matrix(beta0CV.co[, , ii]),
        r, force = force, cv_tol, max.iteration
    )
    resid_ii <- YY.co[estCV[[ii]]] - est.cv.co$fit[estCV[[ii]]]
    list(
        resid    = resid_ii,
        time_idx = rep("Control", length(resid_ii)),
        obs_w    = if (!is.null(W)) W.use[estCV[[ii]]] else NULL
    )
}


## --------------------------------------------------------------------------
## 2f. nevertreated CFE treated_units variant
##     Same structure as 2d: F.hat, U.tr, pre passed in as arguments.
##     Body is identical to IFE treated_units — factor loadings computation
##     is independent of whether the control fit used IFE or CFE.
## --------------------------------------------------------------------------
.fect_cv_score_one_cfe_nt_tr <- function(
    ii,
    U.tr,
    F.hat,
    pre,
    r,
    force,
    rmCV.tr,
    estCV.tr,
    W.tr,
    T.on,
    tr,
    TT,
    Ntr
) {
    pre.cv <- pre; pre.cv[rmCV.tr[[ii]]] <- 0
    if (r == 0) {
        if (force %in% c(1, 3)) {
            alpha.cv <- rep(0, Ntr)
            for (j in 1:Ntr) {
                pre_t <- which(pre.cv[, j] == 1)
                if (length(pre_t) > 0) alpha.cv[j] <- mean(U.tr[pre_t, j])
            }
            fitted.cv <- matrix(alpha.cv, TT, Ntr, byrow = TRUE)
        } else {
            fitted.cv <- matrix(0, TT, Ntr)
        }
    } else {
        fitted.cv <- matrix(0, TT, Ntr)
        for (j in 1:Ntr) {
            pre_t <- which(pre.cv[, j] == 1)
            if (length(pre_t) >= ncol(F.hat)) {
                F.pre <- as.matrix(F.hat[pre_t, , drop = FALSE])
                lam <- try(
                    solve(t(F.pre) %*% F.pre) %*% t(F.pre) %*% U.tr[pre_t, j],
                    silent = TRUE
                )
                if (!"try-error" %in% class(lam)) fitted.cv[, j] <- F.hat %*% lam
            }
        }
    }
    resid_ii <- U.tr[estCV.tr[[ii]]] - fitted.cv[estCV.tr[[ii]]]
    obs_w_ii <- if (!is.null(W.tr)) W.tr[estCV.tr[[ii]]] else NULL
    idx_ii <- NULL
    if (!is.null(T.on)) {
        T.on.tr <- T.on[, tr, drop = FALSE]
        idx_ii <- as.character(T.on.tr[estCV.tr[[ii]]])
        idx_ii[is.na(idx_ii)] <- "Control"
    }
    list(resid = resid_ii, time_idx = idx_ii, obs_w = obs_w_ii)
}


## --------------------------------------------------------------------------
## 3. Rolling-window mask builder
##    Shared helper for cv.method = "rolling" across fect_cv,
##    fect_nevertreated, fect_binary_cv, and fect_mspe.
##
##    Returns a list of length `k`, each element a list with the same
##    shape as cv.sample()'s output:
##        list(cv.id  = integer,   ## flat column-major indices to mask
##             est.id = integer)   ## flat column-major indices to score
##
##    Per fold (uses set.seed(seed + fold_id) when `seed` is supplied):
##      1. Sample max(1L, round(cv.prop * n_eligible)) units from the
##         eligible pool. Eligible = controls (all observed times) +
##         treated (observed times strictly before treatment onset).
##         A unit is eligible only if it has at least min.T0 + cv.nobs
##         observations in its eligible window.
##      2. For each sampled unit, draw a random anchor index a_idx in
##         seq.int(min.T0 + 1L, n_obs - cv.nobs + 1L). The mask
##         partitions cells into:
##           - holdout: obs_t[a_idx:(a_idx + cv.nobs - 1L)]
##                      -> added to BOTH cv.id (masked) and est.id (scored)
##           - past-side buffer (if cv.buffer > 0):
##             obs_t[max(1, a_idx - cv.buffer):(a_idx - 1L)]
##                      -> added to cv.id only (masked but not scored)
##           - drop-future: obs_t[(a_idx + cv.nobs):n_obs] (if any)
##                      -> added to cv.id only (rolling-window step)
##
##    (t, j) is mapped to vec(II) position via (j - 1L) * TT + t.
##
##    II : TT x N indicator matrix (II[t, j] = 1 means unit j observed at t)
##    D  : TT x N treatment indicator (used to derive treatment-onset times)
## --------------------------------------------------------------------------
.build_cv_mask_rolling <- function(II, D, k, cv.nobs, cv.buffer, cv.prop,
                                    min.T0, seed = NULL) {
    II <- as.matrix(II)
    D  <- as.matrix(D)
    if (!all(dim(II) == dim(D))) {
        stop(".build_cv_mask_rolling: II and D must have the same shape.")
    }
    TT <- nrow(II)
    N  <- ncol(II)
    k         <- as.integer(k);         if (k         < 1L) stop("k must be >= 1.")
    cv.nobs   <- as.integer(cv.nobs);   if (cv.nobs   < 1L) stop("cv.nobs must be >= 1.")
    cv.buffer <- as.integer(cv.buffer); if (cv.buffer < 0L) stop("cv.buffer must be >= 0.")
    min.T0    <- as.integer(min.T0);    if (min.T0    < 1L) stop("min.T0 must be >= 1.")
    if (!is.finite(cv.prop) || cv.prop <= 0 || cv.prop > 1) {
        stop(".build_cv_mask_rolling: cv.prop must satisfy 0 < cv.prop <= 1.")
    }

    ## Per-unit treatment-onset time (NA when never treated).
    onset <- rep(NA_integer_, N)
    for (j in seq_len(N)) {
        treated_t <- which(D[, j] >= 1)
        if (length(treated_t) > 0L) onset[j] <- min(treated_t)
    }

    ## Per-unit eligible time vector:
    ##   - controls: all observed t (II == 1)
    ##   - treated:  observed t strictly before onset
    elig_times <- vector("list", N)
    for (j in seq_len(N)) {
        obs_t <- which(II[, j] == 1L)
        if (!is.na(onset[j])) obs_t <- obs_t[obs_t < onset[j]]
        elig_times[[j]] <- obs_t
    }

    ## Filter to units with enough eligible observations.
    elig_lengths <- vapply(elig_times, length, integer(1))
    eligible_units <- which(elig_lengths >= (min.T0 + cv.nobs))
    n_eligible <- length(eligible_units)

    if (n_eligible == 0L) {
        stop(".build_cv_mask_rolling: no eligible units have enough ",
             "observations (need >= min.T0 + cv.nobs = ",
             min.T0 + cv.nobs, ").")
    }

    n_sample_per_fold <- max(1L, as.integer(round(cv.prop * n_eligible)))

    folds <- vector("list", k)
    for (fold_id in seq_len(k)) {
        if (!is.null(seed)) set.seed(as.integer(seed) + fold_id)

        sampled <- if (n_sample_per_fold >= n_eligible) {
            eligible_units
        } else {
            eligible_units[sort(sample.int(n_eligible, n_sample_per_fold))]
        }

        cv_id_acc  <- integer(0)
        est_id_acc <- integer(0)

        for (j in sampled) {
            obs_t <- elig_times[[j]]
            n_obs <- length(obs_t)
            valid <- seq.int(min.T0 + 1L, n_obs - cv.nobs + 1L)
            if (length(valid) == 0L) next
            a_idx <- valid[sample.int(length(valid), 1L)]

            holdout_t <- obs_t[a_idx:(a_idx + cv.nobs - 1L)]
            buf_t <- if (cv.buffer > 0L) {
                bs <- max(1L, a_idx - cv.buffer)
                obs_t[bs:(a_idx - 1L)]
            } else integer(0)
            drop_t <- if ((a_idx + cv.nobs) <= n_obs) {
                obs_t[(a_idx + cv.nobs):n_obs]
            } else integer(0)

            base <- (j - 1L) * TT
            holdout_idx <- base + as.integer(holdout_t)
            buf_idx     <- if (length(buf_t)  > 0L) base + as.integer(buf_t)  else integer(0)
            drop_idx    <- if (length(drop_t) > 0L) base + as.integer(drop_t) else integer(0)

            cv_id_acc  <- c(cv_id_acc, holdout_idx, buf_idx, drop_idx)
            est_id_acc <- c(est_id_acc, holdout_idx)
        }

        folds[[fold_id]] <- list(
            cv.id  = sort(unique(cv_id_acc)),
            est.id = sort(unique(est_id_acc))
        )
    }

    folds
}
