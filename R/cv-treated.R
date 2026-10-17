## ======================================================================
## The treated-unit rank rule: cv.method = "pretreatment" (method = "gsynth")
## ======================================================================
##
## Rank choice from the treated units' own pretreatment outcomes (Xu 2017,
## Algorithm 1, extended to staggered adoption and unbalanced panels). For
## each candidate rank r the controls' gsynth fit gives the grand mean,
## the time effects, the covariate coefficients and the factors. Each treated
## unit's untreated outcome is mu + xi_t + X_t' b + alpha_i + F_t' lambda_i;
## the rule refits the unit's intercept (when the model has unit effects) and
## loadings on part of its pretreatment periods and predicts the rest:
##
##   "loo" (default): hold out each pretreatment period in turn and train on
##       the unit's other pretreatment periods;
##   "fwd" (a check for long pretreatment periods): score period c + 1 from
##       periods 1..c, for the last K cutoffs, where the shortest training set
##       has one more period than the largest rank needs.
##
## A unit's loss at rank r is the mean squared prediction error over its scored
## periods; the score of rank r is the sum of the units' losses weighted by
## each unit's share of the ATT (its observed post-treatment cells). The
## candidate set and the scored units are fixed before any score is computed,
## and every rank is scored on the same cells. The chosen rank is the smallest
## minimizer.
##
## The arithmetic is that of the prototype `treated_rank()` (gsynth-revisit,
## R/rank_rule.R); the functions here take matrices, not fits, so they can be
## tested on their own.

## the value of cv.method this file implements (one constant, so a rename is
## one line)
.FECT_CV_PRETREATMENT <- "pretreatment"

## share of the ATT weight the units scored at the largest rank must hold
.FECT_CV_TREATED_MIN_SHARE <- 0.8

## Scored tasks of one unit: a list of (train, score) period indices. `p` is
## the unit's pretreatment periods in time order, `k` the number of regressors
## at the largest candidate rank.
.fect_cv_treated_unit_tasks <- function(p, k, window) {
    n <- length(p)
    if (window == "loo") {
        if (n - 1 < k) return(list())
        lapply(seq_len(n), function(s) list(train = p[-s], score = p[s]))
    } else {
        K <- n - k - 1  # shortest training set has k + 1 periods
        if (K < 1) return(list())
        lapply(n - seq_len(K), function(c0) {
            list(train = p[seq_len(c0)], score = p[c0 + 1])
        })
    }
}

## The pieces the rule needs for one rank: the treated units' outcome net of
## the offset (R, TT x Ntr), the pre- and post-treatment indicators (TT x Ntr),
## and the regressors Z = [1 if has_a, F] (TT x k).
.fect_cv_treated_parts <- function(R, pre, post, F, has_a, ids) {
    TT <- nrow(R)
    Z <- cbind(if (has_a) rep(1, TT), F)
    if (is.null(Z)) Z <- matrix(0, TT, 0)
    list(ids = ids, R = R, pre = pre, post = post, Z = Z, has_a = has_a,
         r = ncol(F))
}

## Scores of both windows at every rank in `parts` (a list named by rank, the
## ranks in increasing order). Returns the score table, both picks and the
## scored-unit counts, with the names of the prototype's treated_rank().
## `eligible` (optional, named by rank or aligned with `parts`): ranks that
## cannot be picked (the control fit did not converge); their scores are still
## computed and reported. NULL (the default) leaves the plain rule and its
## returned list exactly as they were.
.fect_cv_treated_rank <- function(parts, has_a,
                                  min.share = .FECT_CV_TREATED_MIN_SHARE,
                                  eligible = NULL) {
    ranks <- sort(as.integer(names(parts)))
    if (length(ranks) == 0L) return(NULL)
    parts <- parts[as.character(ranks)]
    ids <- parts[[1]]$ids
    P1 <- parts[[1]]
    n_pre <- colSums(P1$pre)
    w_att <- colSums(P1$post)
    w_att <- w_att / sum(w_att)
    ## Rbar: the largest rank whose leave-one-out tasks exist for units
    ## holding at least min.share of the ATT weight
    ok_at <- function(rb) n_pre - 1 >= rb + has_a
    cand <- ranks[vapply(ranks, function(rb) {
        sum(w_att[ok_at(rb)]) >= min.share
    }, TRUE)]
    if (length(cand) == 0L) {
        stop("cv.method = \"pretreatment\": too few pretreatment periods to ",
             "score any rank; supply r or use another cv.method.",
             call. = FALSE)
    }
    rbar <- max(cand)
    ranks <- ranks[ranks <= rbar]
    elig <- NULL
    if (!is.null(eligible)) {
        elig <- if (is.null(names(eligible))) {
            as.logical(eligible)[seq_along(ranks)]
        } else {
            as.logical(eligible[as.character(ranks)])
        }
        elig[is.na(elig)] <- TRUE
    }
    k <- rbar + has_a
    pre_idx <- lapply(seq_along(ids), function(j) which(P1$pre[, j]))

    score <- function(window) {
        tasks <- lapply(seq_along(ids), function(j) {
            .fect_cv_treated_unit_tasks(pre_idx[[j]], k, window)
        })
        use <- which(vapply(tasks, length, 1L) > 0L)
        if (length(use) == 0L) {
            return(list(s = rep(NA_real_, length(ranks)), n = 0L, share = 0))
        }
        wu <- w_att[use]
        wu <- wu / sum(wu)
        s <- vapply(ranks, function(r) {
            P <- parts[[as.character(r)]]
            tot <- 0
            for (u in seq_along(use)) {
                j <- use[u]
                e <- vapply(tasks[[j]], function(t) {
                    if (ncol(P$Z) == 0) return(P$R[t$score, j])
                    Zt <- P$Z[t$train, , drop = FALSE]
                    ## rank deficient: the task is unusable
                    if (qr(Zt)$rank < ncol(Zt)) return(NA_real_)
                    P$R[t$score, j] -
                        drop(P$Z[t$score, , drop = FALSE] %*%
                                 qr.solve(Zt, P$R[t$train, j]))
                }, 0)
                if (anyNA(e)) return(Inf)  # never chosen; reported
                tot <- tot + wu[u] * mean(e^2)
            }
            tot
        }, 0)
        list(s = s, n = length(use), share = sum(w_att[use]))
    }
    lo <- score("loo")
    fw <- score("fwd")
    pick <- function(s) {
        if (all(!is.finite(s))) return(NA_integer_)
        if (!is.null(elig)) {
            ok <- is.finite(s) & elig
            if (any(ok)) return(ranks[which.min(ifelse(ok, s, Inf))])
        }
        ranks[which.min(s)]
    }
    ## every rank with a finite score is ineligible: the unrestricted pick
    fell_back <- function(s) {
        !is.null(elig) && any(is.finite(s)) && !any(is.finite(s) & elig)
    }
    res <- list(
        table = data.frame(r = ranks, loo = lo$s, fwd = fw$s),
        r_loo = pick(lo$s), r_fwd = pick(fw$s), rbar = rbar,
        n_tr = length(ids), n_loo = lo$n, n_fwd = fw$n,
        share_loo = lo$share, share_fwd = fw$share,
        note = sprintf(
            paste0("ranks 0-%d; leave one out scores %d of %d treated units ",
                   "(%.0f%% of the ATT weight), forward %d (%.0f%%)"),
            rbar, lo$n, length(ids), 100 * lo$share, fw$n, 100 * fw$share)
    )
    if (!is.null(elig)) {
        res$eligible <- elig
        res$fallback <- if (identical(getOption("fect.cv.treated.window", "loo"),
                                      "fwd")) fell_back(fw$s) else fell_back(lo$s)
    }
    res
}

## The rule inside fect_nevertreated(): fit the controls at each candidate
## rank at the user's tolerance, build the pieces of the treated units, run the
## rule once. `fit.rank(r)` returns the controls' fit at rank r (the solver
## wrapper of fect_nevertreated), `Y.tr` and `X.tr` are the treated units'
## outcome (TT x Ntr) and covariates (TT x Ntr x p). Returns the per-rank fits,
## the rule's output and the window that sets the pick.
.fect_cv_treated_run <- function(ranks, fit.rank, Y.tr, X.tr, p, force,
                                 pre, post, ids, exclude.nonconverged = FALSE,
                                 max.iteration = NULL) {
    TT <- nrow(Y.tr)
    Ntr <- ncol(Y.tr)
    has_a <- force %in% c(1, 3)
    has_xi <- force %in% c(2, 3)
    fits <- vector("list", length(ranks))
    parts <- vector("list", length(ranks))
    for (i in seq_along(ranks)) {
        est.co <- fit.rank(ranks[i])
        fits[[i]] <- est.co
        ## the offset fect subtracts before it fits the treated unit
        U <- Y.tr
        if (p > 0) {
            beta <- est.co$beta
            if (!is.null(est.co$validX) && est.co$validX == 0) {
                beta <- matrix(0, p, 1)
            } else {
                beta[is.nan(beta)] <- 0
            }
            for (j in seq_len(p)) {
                U <- U - matrix(X.tr[, , j], TT, Ntr) * beta[j]
            }
        }
        U <- U - matrix(est.co$mu, TT, Ntr)
        if (has_xi) U <- U - matrix(c(est.co$xi), TT, Ntr)
        F <- if (ranks[i] > 0) {
            as.matrix(est.co$factor)
        } else {
            matrix(0, TT, 0)
        }
        parts[[i]] <- .fect_cv_treated_parts(U, pre, post, F, has_a, ids)
    }
    names(parts) <- as.character(ranks)
    ## a rank whose control fit stopped at max.iteration cannot be picked
    ## (only the serial = "ar1" alternation asks for this)
    elig <- NULL
    if (isTRUE(exclude.nonconverged) && !is.null(max.iteration)) {
        elig <- vapply(fits, function(f) {
            !(length(f$niter) == 1L && is.finite(f$niter) &&
                  f$niter >= max.iteration)
        }, TRUE)
        names(elig) <- as.character(ranks)
    }
    res <- .fect_cv_treated_rank(parts, has_a, eligible = elig)
    window <- getOption("fect.cv.treated.window", "loo")
    if (!window %in% c("loo", "fwd")) {
        stop("option fect.cv.treated.window must be \"loo\" or \"fwd\".",
             call. = FALSE)
    }
    list(fits = fits, rule = res, window = window)
}
