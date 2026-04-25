## -----------------------------------------------------------------------------
## loading_bound.R
##
## Entropy-regularized simplex projection of treated factor loadings for GSC.
##
## Solves, per treated unit i:
##
##   minimize   (1/gamma) * sum_j w_j * log(w_j / q_j)
##            + || u_pre - F_pre %*% t(Lambda_co) %*% w ||_2^2
##   over       w in Delta_{Nco}  (simplex)
##
## The resulting loading is lambda_hat = t(Lambda_co) %*% w, which lies in the
## convex hull of control loadings by construction. Under the factor model this
## yields Y_hat(0) = F %*% lambda_hat = F %*% t(Lambda_co) %*% w, a convex
## combination of factor-implied control outcomes.
##
## See statsclaw-workspace/fect/runs/REQ-bounded-loadings/spec.md for the full
## design and the theoretical justification.
## -----------------------------------------------------------------------------

.bounded_loading_obj <- function(theta, u_pre, F_pre, Lambda_co, gamma, q) {
    ## Softmax reparameterization: w_j = exp(theta_j) / sum_k exp(theta_k).
    th   <- theta - max(theta)            # numerical stability
    ew   <- exp(th)
    w    <- ew / sum(ew)
    w_safe <- pmax(w, 1e-300)
    lam  <- crossprod(Lambda_co, w)       # r-vector
    res  <- u_pre - F_pre %*% lam
    ent  <- sum(w * (log(w_safe) - log(q)))
    val  <- (1 / gamma) * ent + sum(res^2)
    if (!is.finite(val)) val <- 1e30      # L-BFGS-B needs finite values
    val
}

.bounded_loading_grad <- function(theta, u_pre, F_pre, Lambda_co, gamma, q) {
    th   <- theta - max(theta)
    ew   <- exp(th)
    w    <- ew / sum(ew)
    w_safe <- pmax(w, 1e-300)
    lam  <- crossprod(Lambda_co, w)
    res  <- u_pre - F_pre %*% lam

    ## d/dw of fit term: -2 * Lambda_co %*% t(F_pre) %*% res   (Nco-vector)
    g_fit_w <- -2 * as.numeric(Lambda_co %*% crossprod(F_pre, res))

    ## d/dw of entropy term: (1/gamma) * (log(w/q) + 1)
    g_ent_w <- (1 / gamma) * (log(w_safe / q) + 1)

    g_w <- g_fit_w + g_ent_w

    ## Chain rule for softmax: g_theta = w * (g_w - sum(w * g_w))
    g <- w * (g_w - sum(w * g_w))
    g[!is.finite(g)] <- 0
    g
}

.mirror_descent_bounded_loading <- function(
    u_pre, F_pre, Lambda_co, gamma, q,
    tol = 1e-8, max_iter = 500L
) {
    Nco <- nrow(Lambda_co)
    w   <- rep(1 / Nco, Nco)
    obj_new <- Inf
    for (iter in seq_len(max_iter)) {
        lam  <- crossprod(Lambda_co, w)
        res  <- u_pre - F_pre %*% lam
        g_w  <- -2 * as.numeric(Lambda_co %*% crossprod(F_pre, res)) +
                 (1 / gamma) * (log(pmax(w, 1e-300) / q) + 1)
        step <- 1 / (iter + 1)             # diminishing step
        ## Clip gradient to prevent exp() overflow on the log-simplex step
        g_w <- pmax(pmin(g_w, 1e6), -1e6)
        lw   <- log(pmax(w, 1e-300)) - step * g_w
        lw   <- lw - max(lw)
        ew   <- exp(lw)
        ew_sum <- sum(ew)
        if (!is.finite(ew_sum) || ew_sum <= 0) {
            ## Numerical failure; fall back to uniform and terminate
            w <- rep(1 / Nco, Nco)
            break
        }
        w_new <- ew / ew_sum
        if (any(!is.finite(w_new))) {
            w <- rep(1 / Nco, Nco)
            break
        }

        obj_new <- sum(res^2) +
            (1 / gamma) * sum(w_new * (log(pmax(w_new, 1e-300)) - log(q)))
        delta <- max(abs(w_new - w))
        if (isTRUE(is.finite(delta) && delta < tol)) {
            w <- w_new
            break
        }
        w <- w_new
    }
    list(w = w, n_iter = iter, obj = obj_new)
}

.solve_bounded_loading <- function(
    u_pre,
    F_pre,
    Lambda_co,
    gamma,
    q        = NULL,
    tol      = 1e-8,
    max_iter = 200L,
    fallback = TRUE
) {
    u_pre     <- as.numeric(u_pre)
    F_pre     <- as.matrix(F_pre)
    Lambda_co <- as.matrix(Lambda_co)
    Nco       <- nrow(Lambda_co)
    stopifnot(length(u_pre) == nrow(F_pre))
    stopifnot(ncol(F_pre) == ncol(Lambda_co))

    if (Nco < 2L) {
        stop("bounded-loading projection requires at least 2 control units; got Nco = ",
             Nco, ".", call. = FALSE)
    }
    if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma <= 0) {
        stop("gamma must be a single finite positive numeric; got ", gamma, ".",
             call. = FALSE)
    }
    if (is.null(q)) {
        q <- rep(1 / Nco, Nco)
    } else {
        q <- as.numeric(q)
        if (length(q) != Nco || any(q <= 0) || abs(sum(q) - 1) > 1e-8) {
            stop("q must be a positive Nco-vector summing to 1.", call. = FALSE)
        }
    }

    ## Degenerate: all-zero Lambda_co -- return uniform w
    if (max(abs(Lambda_co)) < 1e-12) {
        w <- q
        return(list(
            w          = w,
            lambda_hat = as.numeric(crossprod(Lambda_co, w)),
            obj        = sum(u_pre^2),
            converged  = TRUE,
            n_iter     = 0L,
            method     = "degenerate-return-uniform"
        ))
    }

    if (Nco <= ncol(Lambda_co)) {
        warning(
            "Nco (", Nco, ") <= r (", ncol(Lambda_co),
            "); conv(Lambda_co) has empty interior in R^r and the projection may ",
            "be highly constrained.",
            call. = FALSE
        )
    }

    theta0 <- rep(0, Nco)
    opt <- tryCatch(
        stats::optim(
            par     = theta0,
            fn      = .bounded_loading_obj,
            gr      = .bounded_loading_grad,
            method  = "L-BFGS-B",
            u_pre   = u_pre,
            F_pre   = F_pre,
            Lambda_co = Lambda_co,
            gamma   = gamma,
            q       = q,
            control = list(maxit = max_iter, factr = max(tol / .Machine$double.eps, 1))
        ),
        error = function(e) list(convergence = 99L, message = conditionMessage(e))
    )

    converged_lbfgs <- isTRUE(opt$convergence == 0L)

    if (converged_lbfgs) {
        th   <- opt$par - max(opt$par)
        ew   <- exp(th)
        w    <- ew / sum(ew)
        lam  <- as.numeric(crossprod(Lambda_co, w))
        return(list(
            w          = w,
            lambda_hat = lam,
            obj        = opt$value,
            converged  = TRUE,
            n_iter     = if (!is.null(opt$counts)) unname(opt$counts[1]) else NA_integer_,
            method     = "lbfgs"
        ))
    }

    if (!fallback) {
        return(list(
            w          = rep(1 / Nco, Nco),
            lambda_hat = as.numeric(crossprod(Lambda_co, rep(1 / Nco, Nco))),
            obj        = NA_real_,
            converged  = FALSE,
            n_iter     = NA_integer_,
            method     = "lbfgs-failed"
        ))
    }

    md <- .mirror_descent_bounded_loading(
        u_pre     = u_pre,
        F_pre     = F_pre,
        Lambda_co = Lambda_co,
        gamma     = gamma,
        q         = q,
        tol       = tol,
        max_iter  = 500L
    )
    list(
        w          = md$w,
        lambda_hat = as.numeric(crossprod(Lambda_co, md$w)),
        obj        = md$obj,
        converged  = TRUE,
        n_iter     = md$n_iter,
        method     = "mirror-descent"
    )
}

## -----------------------------------------------------------------------------
## .cv_gamma_loading
##
## Selects gamma by 5-fold CV over pre-treatment periods. Folds are constructed
## across pre-treatment time indices (shared across treated units); per-fold
## score is the mean squared prediction error of the bounded projection on the
## held-out pre-treatment observations.
## -----------------------------------------------------------------------------

.cv_gamma_loading <- function(
    U_tr_pre,        # T0 x Ntr matrix OR list of length Ntr of unit-specific vectors
    F_hat_pre,       # T0 x r matrix (balanced case) OR list per unit (unbalanced)
    Lambda_co,       # Nco x r matrix
    gamma_grid,      # positive numeric vector
    cv_k        = 5L,
    tol         = 1e-8,
    max_iter    = 200L
) {
    Lambda_co <- as.matrix(Lambda_co)
    stopifnot(all(gamma_grid > 0))

    is_balanced <- is.matrix(U_tr_pre) && is.matrix(F_hat_pre)

    if (is_balanced) {
        T0  <- nrow(U_tr_pre)
        Ntr <- ncol(U_tr_pre)
        fold_id <- sample(rep_len(seq_len(cv_k), T0))
    } else {
        Ntr <- length(U_tr_pre)
        fold_id <- NULL   # per-unit folding for unbalanced
    }

    cv_mspe <- numeric(length(gamma_grid))
    for (gi in seq_along(gamma_grid)) {
        gamma_i <- gamma_grid[gi]
        sq_err  <- 0
        n_obs   <- 0L

        if (is_balanced) {
            for (ii in seq_len(cv_k)) {
                tr_idx <- which(fold_id != ii)
                te_idx <- which(fold_id == ii)
                if (length(tr_idx) < ncol(F_hat_pre) || length(te_idx) == 0L) next
                F_tr <- F_hat_pre[tr_idx, , drop = FALSE]
                F_te <- F_hat_pre[te_idx, , drop = FALSE]
                for (i in seq_len(Ntr)) {
                    u_tr <- U_tr_pre[tr_idx, i]
                    u_te <- U_tr_pre[te_idx, i]
                    sol <- .solve_bounded_loading(
                        u_pre     = u_tr,
                        F_pre     = F_tr,
                        Lambda_co = Lambda_co,
                        gamma     = gamma_i,
                        tol       = tol,
                        max_iter  = max_iter
                    )
                    pred <- as.numeric(F_te %*% sol$lambda_hat)
                    sq_err <- sq_err + sum((u_te - pred)^2)
                    n_obs  <- n_obs + length(u_te)
                }
            }
        } else {
            for (i in seq_len(Ntr)) {
                u_i <- as.numeric(U_tr_pre[[i]])
                F_i <- as.matrix(F_hat_pre[[i]])
                T0i <- length(u_i)
                if (T0i < cv_k + ncol(F_i)) next
                fold_i <- sample(rep_len(seq_len(cv_k), T0i))
                for (ii in seq_len(cv_k)) {
                    tr_idx <- which(fold_i != ii)
                    te_idx <- which(fold_i == ii)
                    if (length(tr_idx) < ncol(F_i) || length(te_idx) == 0L) next
                    sol <- .solve_bounded_loading(
                        u_pre     = u_i[tr_idx],
                        F_pre     = F_i[tr_idx, , drop = FALSE],
                        Lambda_co = Lambda_co,
                        gamma     = gamma_i,
                        tol       = tol,
                        max_iter  = max_iter
                    )
                    pred <- as.numeric(F_i[te_idx, , drop = FALSE] %*% sol$lambda_hat)
                    sq_err <- sq_err + sum((u_i[te_idx] - pred)^2)
                    n_obs  <- n_obs + length(te_idx)
                }
            }
        }

        cv_mspe[gi] <- if (n_obs > 0L) sq_err / n_obs else Inf
    }

    best <- which.min(cv_mspe)
    list(
        gamma_cv   = gamma_grid[best],
        cv_mspe    = cv_mspe,
        gamma_grid = gamma_grid
    )
}

## Default gamma grid: 9 log-spaced points over [0.01, 100]
.default_gamma_grid <- function() 10 ^ seq(-2, 2, length.out = 9L)
