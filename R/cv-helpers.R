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
    est.cv.fit <- fect:::inter_fe_ub(
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
    est.cv.fit <- fect:::inter_fe_mc(
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
    est.cv.fit <- fect:::inter_fe_ub(
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
    est.cv.co <- fect:::complex_fe_ub(
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
