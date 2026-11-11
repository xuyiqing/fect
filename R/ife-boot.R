## Parametric bootstrap for IFE-EM (fect 2.5.0, experimental): Algorithm
## "A parametric bootstrap for IFE-EM" of the gsynth-revisit paper
## (Supplementary Material), with its extension to covariates and missing
## cells. fect_boot() calls these helpers for method = "ife",
## time.component.from = "notyettreated", vartype = "parametric".
##
##   1. The point fit (IFE-EM with r factors on the untreated cells) gives the
##      fitted untreated outcome of every cell, X beta included.
##   2. Split the never-treated controls at random into G = min(Nco, 30)
##      groups (one control per group when Nco <= 30). For each group, refit
##      IFE-EM with r_e factors without the group's controls, r_e = max(r - 1,
##      0) when Nco >= 30, else r. Each control in the group gets a noise
##      vector: its outcome net of the refit's grand mean, time effects and
##      X beta, regressed on an intercept and the refit's factors over its
##      observed periods; the residual times sqrt(n / (n - r_e - 1)), n its
##      number of observed periods. A control observed in r_e + 1 or fewer
##      periods gets none.
##   3. For b = 1..B: every unit, controls and treated alike, draws a noise
##      vector with replacement; a cell observed in the data whose drawn
##      vector is missing there takes the noise of another control drawn at
##      random among those observed in that period. The bootstrap panel is the
##      fitted untreated outcome plus the noise, so the effect is zero in it;
##      it has the data's missing cells. IFE-EM is refitted with r factors and
##      the same mask, and its average effect is a draw of the estimation
##      error under no effect.
##
## Random numbers are drawn in the order of the replication code
## (sim/ife_boot_x.R: ife_noise_group_x() then ife_pb_x()): the group split,
## then for each draw the N noise vectors and then the patches, so the same
## seed gives the same groups and draws. All draws are made before the refits,
## so they do not depend on parallel computing.

## Noise rank of the algorithm.
.ife_noise_rank <- function(r, n.co) if (n.co >= 30) max(r - 1L, 0L) else as.integer(r)

## Leave-group-out noise vectors: a T x Nco matrix (NA where a control is not
## observed or has no vector), on the scale of Y.
.ife_noise_group <- function(Y, X, D, I, II, T.on, co, r.noise, G, force, tol,
                             max.iteration) {
    TT <- nrow(Y)
    Nco <- length(co)
    p <- if (is.null(X) || length(dim(X)) < 3) 0L else dim(X)[3]
    grp <- if (G >= Nco) seq_len(Nco) else sample(rep_len(seq_len(G), Nco))
    E <- matrix(NA_real_, TT, Nco)
    failed <- character(0)
    for (g in unique(grp)) {
        out.id <- co[grp == g]
        keep <- setdiff(seq_len(ncol(Y)), out.id)
        f <- tryCatch(
            fect_fe(Y = Y[, keep, drop = FALSE],
                    X = if (p > 0) X[, keep, , drop = FALSE] else NULL,
                    D = D[, keep, drop = FALSE], W = NULL,
                    I = I[, keep, drop = FALSE], II = II[, keep, drop = FALSE],
                    T.on = T.on[, keep, drop = FALSE], r.cv = r.noise,
                    force = force, hasRevs = 0, tol = tol,
                    max.iteration = max.iteration, boot = TRUE),
            error = function(e) e)
        if (inherits(f, "error")) {
            failed <- c(failed, conditionMessage(f))
            next
        }
        Z <- cbind(rep(1, TT), if (r.noise > 0) as.matrix(f$factor))
        base <- rep(as.numeric(f$mu), TT)
        if (!is.null(f$xi)) base <- base + as.vector(f$xi)
        beta <- if (p > 0) as.vector(f$beta) else numeric(0)
        beta[!is.finite(beta)] <- 0
        for (k in which(grp == g)) {
            j <- co[k]
            R <- Y[, j] - base
            if (p > 0) R <- R - drop(matrix(X[, j, , drop = FALSE], TT, p) %*% beta)
            ok <- I[, j] == 1
            n <- sum(ok)
            if (n <= ncol(Z)) next
            Zo <- Z[ok, , drop = FALSE]
            E[ok, k] <- (R[ok] - drop(Zo %*% qr.solve(Zo, R[ok]))) * sqrt(n / (n - ncol(Z)))
        }
    }
    list(E = E, groups = grp, failed = failed)
}

## The B draws: each unit's noise vector (column of the pool) and the patches
## (cell, pool column) for observed cells whose drawn vector is missing.
.ife_draws <- function(pool, I, B) {
    TT <- nrow(I)
    N <- ncol(I)
    patch <- any(I == 0) || anyNA(pool)
    obs <- I == 1
    lapply(seq_len(B), function(b) {
        idx <- sample.int(ncol(pool), N, replace = TRUE)
        pt <- NULL
        if (patch) {
            E <- pool[, idx, drop = FALSE]
            need <- which(obs & is.na(E), arr.ind = TRUE)
            if (nrow(need)) {
                pt <- matrix(0L, nrow(need), 2L)
                for (q in seq_len(nrow(need))) {
                    t <- need[q, 1]
                    av <- which(is.finite(pool[t, ]))
                    if (!length(av))
                        stop("IFE-EM bootstrap: no control has a noise value in period ", t,
                             call. = FALSE)
                    pt[q, ] <- c(t + (need[q, 2] - 1L) * TT, av[sample.int(length(av), 1)])
                }
            }
        }
        list(idx = idx, pt = pt)
    })
}

## The noise matrix of one draw (0 in unobserved cells).
.ife_noise_matrix <- function(pool, d, I) {
    E <- pool[, d$idx, drop = FALSE]
    if (!is.null(d$pt)) {
        TT <- nrow(pool)
        rows <- (d$pt[, 1] - 1L) %% TT + 1L
        E[d$pt[, 1]] <- pool[cbind(rows, d$pt[, 2])]
    }
    E[I == 0] <- 0
    E
}
