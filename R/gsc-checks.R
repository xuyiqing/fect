## Overlap checks for the generalized synthetic control (GSC) estimator
## (fect 2.5.0), following the procedures of the gsynth-revisit paper
## (Supplementary Material, "The two overlap checks"):
##
##   overlap_check()  the check that needs no rank: before any fit, can each
##                    treated unit's demeaned pretreatment path be rebuilt as a
##                    weighted average of the never-treated units' paths, as
##                    well as the never-treated units rebuild one another?
##   loading_check()  the check at the chosen rank: are the treated units'
##                    estimated loadings inside the controls' hull, allowing
##                    for the noise in their own estimates? Optionally drops
##                    the flagged units, chooses the rank again and refits.
##   extrapolation    the share of treated post-treatment counterfactuals
##                    outside the range of the observed outcomes.
##
## The reference implementation is the paper's replication code
## (R/overlap_free.R and R/rank_rule.R). The 2.5.0 run log in workbench
## (statsclaw-workspace/fect/runs/) records the comparison on the paper's cases.

## ---------------------------------------------------------------------------
## Least squares on the simplex
## ---------------------------------------------------------------------------

## min over w >= 0, sum(w) = 1 of |y - X w|^2, for an m x J matrix X.
## Dense quadratic program (quadprog) with a small ridge for numerical
## stability when J <= dense.max (the replication code's solver, kept so that
## the two agree to rounding); above that, an exact active-set method that
## solves the same problem on a growing set of columns: the solution of a
## least-squares problem on the simplex has at most m + 1 positive weights,
## so the working set stays small. Returns the weights.
.gsc_simplex_ls <- function(y, X, dense.max = 200L, ridge = 1e-10) {
    y <- as.vector(y)
    X <- as.matrix(X)
    J <- ncol(X)
    if (J == 1L) return(1)
    qp <- function(Xs, eps0) {
        D <- crossprod(Xs)
        s <- max(mean(diag(D)), 1e-12)
        k <- ncol(Xs)
        for (eps in eps0 * c(1, 1e2, 1e4)) {
            w <- tryCatch(
                quadprog::solve.QP(D + diag(eps * s, k), drop(crossprod(Xs, y)),
                                   cbind(rep(1, k), diag(k)), c(1, rep(0, k)),
                                   meq = 1)$solution,
                error = function(e) NULL)
            if (!is.null(w)) return(w)
        }
        ## badly scaled data (outcomes in the tens of thousands): the same
        ## problem divided by a constant has the same solution
        sc <- max(abs(Xs), abs(y), 1e-12)
        if (sc != 1) {
            D <- crossprod(Xs / sc)
            s <- max(mean(diag(D)), 1e-12)
            for (eps in eps0 * c(1, 1e2, 1e4)) {
                w <- tryCatch(
                    quadprog::solve.QP(D + diag(eps * s, k), drop(crossprod(Xs / sc, y / sc)),
                                       cbind(rep(1, k), diag(k)), c(1, rep(0, k)),
                                       meq = 1)$solution,
                    error = function(e) NULL)
                if (!is.null(w)) return(w)
            }
        }
        stop("least squares on the simplex failed", call. = FALSE)
    }
    if (J <= dense.max) return(qp(X, ridge))
    ## active set: start from the closest column, add the columns that
    ## violate the optimality conditions most, until none does
    scale <- max(mean(colSums(X^2)), 1e-12)
    S <- which.min(colSums((X - y)^2))
    for (it in seq_len(10L * J)) {
        wS <- if (length(S) == 1L) 1 else qp(X[, S, drop = FALSE], ridge)
        res <- drop(X[, S, drop = FALSE] %*% wS) - y
        g <- drop(crossprod(X, res))              # half the gradient
        lam <- sum(wS * g[S])                      # the multiplier of sum(w) = 1
        viol <- g - lam                            # optimal when all >= 0
        viol[S] <- pmax(viol[S], 0)
        if (min(viol) >= -1e-12 * scale) break
        add <- setdiff(order(viol)[seq_len(min(5L, J))], S)
        add <- add[viol[add] < -1e-12 * scale]
        keep <- S[wS > 0]
        if (!length(keep)) keep <- S
        S <- c(keep, add)
    }
    w <- numeric(J)
    w[S] <- wS
    w
}

## RMSE of y (length k) rebuilt from the columns of X (k x J) after both are
## demeaned over the window.
.gsc_sc_rmse <- function(y, X) {
    y <- y - mean(y)
    X <- sweep(X, 2, colMeans(X))
    if (ncol(X) == 1L) return(sqrt(mean((y - X[, 1])^2)))
    w <- .gsc_simplex_ls(y, X)
    sqrt(mean((y - drop(X %*% w))^2))
}

## Run `expr` with the RNG seeded, restoring the caller's RNG state.
.gsc_with_seed <- function(seed, expr) {
    old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
    on.exit({
        if (is.null(old)) {
            if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
                rm(".Random.seed", envir = globalenv())
        } else assign(".Random.seed", old, envir = globalenv())
    })
    set.seed(seed)
    expr
}

## ---------------------------------------------------------------------------
## The check that needs no rank
## ---------------------------------------------------------------------------

.of_na <- function(reason, n.periods = NA_integer_, n.donors = NA_integer_)
    list(status = "not assessable", reason = reason, error = NA_real_,
         yardstick = NA_real_, median = NA_real_, score = NA_real_, flag = NA,
         pct = NA_real_, n.periods = n.periods, n.donors = n.donors,
         e.co = numeric(0))

.of_judge <- function(err, e.co, Y, rows, donors, level, collapse) {
    med <- stats::median(e.co)
    sdy <- stats::sd(as.vector(Y[rows, donors]))
    if (!is.finite(med) || med < collapse * max(sdy, 1e-12))
        return(.of_na("yardstick collapsed", length(rows), length(donors)))
    q <- stats::quantile(e.co, level, names = FALSE)
    list(status = "ok", reason = "", error = err, yardstick = q, median = med,
         score = err / med, flag = err > q, pct = mean(e.co < err),
         n.periods = length(rows), n.donors = length(donors), e.co = e.co)
}

.of_donors <- function(Y, rows, co) co[colSums(is.na(Y[rows, co, drop = FALSE])) == 0]

.of_unit <- function(Y, i, pre, co, cache, pars) {
    rows <- pre[!is.na(Y[pre, i])]
    if (length(rows) < pars$min.periods)
        return(.of_na("too few periods", length(rows), NA))
    donors <- .of_donors(Y, rows, co)
    if (length(donors) < pars$min.donors)
        return(.of_na("too few donors", length(rows), length(donors)))
    ## reuse the yardstick of a window already seen (same periods and donors);
    ## the cache is a list searched in full, since a key built from the
    ## indices can exceed R's limit on the length of a name
    e.co <- NULL
    for (ce in cache$entries)
        if (identical(ce$rows, rows) && identical(ce$donors, donors)) {
            e.co <- ce$e.co
            break
        }
    if (is.null(e.co)) {
        held <- donors
        if (length(held) > pars$max.heldout)
            held <- sort(.gsc_with_seed(pars$seed, sample(donors, pars$max.heldout)))
        e.co <- vapply(held, function(j)
            .gsc_sc_rmse(Y[rows, j], Y[rows, setdiff(donors, j), drop = FALSE]), 0)
        cache$entries[[length(cache$entries) + 1L]] <- list(rows = rows, donors = donors, e.co = e.co)
    }
    .of_judge(.gsc_sc_rmse(Y[rows, i], Y[rows, donors, drop = FALSE]), e.co,
              Y, rows, donors, pars$level, pars$collapse)
}

.of_cohort <- function(Y, grp, pre, co, seed, pars) {
    obs <- !is.na(Y[pre, grp, drop = FALSE])
    rows <- pre[rowMeans(obs) >= 0.5]
    if (length(rows) < pars$min.periods)
        return(.of_na("too few periods", length(rows), NA))
    donors <- .of_donors(Y, rows, co)
    if (length(donors) < pars$min.donors + 1L)
        return(.of_na("too few donors", length(rows), length(donors)))
    y <- rowMeans(Y[rows, grp, drop = FALSE], na.rm = TRUE)
    m <- min(length(grp), floor(length(donors) / 2))
    e.co <- .gsc_with_seed(seed, vapply(seq_len(pars$groups), function(g) {
        s <- sample(donors, m)
        .gsc_sc_rmse(rowMeans(Y[rows, s, drop = FALSE]),
                     Y[rows, setdiff(donors, s), drop = FALSE])
    }, 0))
    out <- .of_judge(.gsc_sc_rmse(y, Y[rows, donors, drop = FALSE]), e.co,
                     Y, rows, donors, pars$level, pars$collapse)
    out$group.size <- length(grp)
    out$yard.size <- m
    out
}

## Wide outcome and treatment matrices (time x unit) from long data, with the
## rows that miss any variable of the formula or the index dropped.
.gsc_wide <- function(formula, data, index) {
    if (!inherits(formula, "formula"))
        stop("\"formula\" must be a formula such as Y ~ D + X1.", call. = FALSE)
    if (!is.data.frame(data)) stop("\"data\" must be a data frame.", call. = FALSE)
    if (length(index) != 2L || !all(index %in% names(data)))
        stop("\"index\" must name the unit and time variables of \"data\".", call. = FALSE)
    vars <- all.vars(formula)
    if (length(vars) < 2L)
        stop("The formula needs an outcome and a treatment: Y ~ D (+ covariates).", call. = FALSE)
    miss <- setdiff(vars, names(data))
    if (length(miss))
        stop("variable(s) not in the data: ", paste(miss, collapse = ", "), call. = FALSE)
    yname <- vars[1]; dname <- vars[2]
    d <- as.data.frame(data)[, unique(c(vars, index))]
    d <- d[stats::complete.cases(d), , drop = FALSE]
    uid <- as.character(d[[index[1]]])
    ids <- sort(unique(uid))
    tt <- sort(unique(d[[index[2]]]))
    if (anyDuplicated(paste(uid, d[[index[2]]], sep = "\r")))
        stop("The data have more than one row for some unit and time.", call. = FALSE)
    pos <- cbind(match(d[[index[2]]], tt), match(uid, ids))
    wide <- function(v) {
        m <- matrix(NA_real_, length(tt), length(ids), dimnames = list(as.character(tt), ids))
        m[pos] <- as.numeric(d[[v]])
        m
    }
    D <- wide(dname)
    if (any(!D %in% c(0, 1, NA)))
        stop("The treatment must be 0 or 1.", call. = FALSE)
    list(Y = wide(yname), D = D, ids = ids, time = tt, yname = yname, dname = dname)
}

#' Overlap check that needs no rank
#'
#' Checks, before any model is fitted, whether each treated unit's
#' pretreatment outcomes can be rebuilt from the never-treated units as well
#' as the never-treated units rebuild one another. This is the first overlap
#' check of the gsynth-revisit paper (Supplementary Material, "The check that
#' needs no rank").
#'
#' For a treated unit, take its observed pretreatment periods and the
#' never-treated units observed in all of them (its donors). Subtract from
#' each path its mean over these periods. The unit's error is the root mean
#' squared error of rebuilding its demeaned path as a weighted average of the
#' donors' demeaned paths, with weights that are non-negative and sum to one.
#' The yardstick is the same error for each donor rebuilt from the other
#' donors (at most \code{max.heldout} donors, a random subset above that).
#' The unit's percentile \code{pct} is the share of donors rebuilt with a
#' smaller error, and the unit is flagged when its error exceeds the
#' \code{level} quantile of the donors' errors. A treated unit has one donor
#' more than a held-out donor, so the flag is conservative.
#'
#' A unit is "not assessable", rather than passing, when it has fewer than
#' four pretreatment periods or five donors, or when the donors rebuild one
#' another almost exactly (median donor error below \code{1e-8} times their
#' outcome standard deviation).
#'
#' The cohort table does the same for the average path of each adoption
#' cohort (the periods in which at least half of the cohort is observed),
#' with a yardstick of \code{groups} random groups of donors of the cohort's
#' size (at most half the donors), each group's average rebuilt from the
#' donors outside it.
#'
#' The check uses no rank, no factor model and no post-treatment outcome. It
#' runs on the outcome as given: covariates in the formula are used only to
#' drop the rows that miss them, as \code{fect()} does. Like every
#' pretreatment check, it cannot see a factor that is flat before adoption.
#'
#' @param formula a formula \code{Y ~ D} or \code{Y ~ D + X1 + ...}; the
#'   treatment must be absorbing (once 1, always 1).
#' @param data a data frame in long form.
#' @param index the names of the unit and time variables.
#' @param level the quantile of the donors' errors above which a unit is
#'   flagged. The paper uses 0.90.
#' @param max.heldout the largest number of donors in the yardstick.
#' @param groups the number of random donor groups in a cohort's yardstick.
#' @param seed seed for the random choices (the held-out donors when there
#'   are more than \code{max.heldout}, and the donor groups); the caller's
#'   random-number state is restored afterwards.
#' @return An object of class \code{"fect_overlap"}: a list with
#'   \code{units} (one row per treated unit: \code{unit}, \code{onset} (first
#'   treated period), \code{weight} (the unit's share of the treated
#'   post-treatment cells, its weight in the average effect), \code{status}
#'   (\code{"ok"} or \code{"not assessable"}), \code{reason}, \code{error},
#'   \code{yardstick} (the \code{level} quantile of the donors' errors),
#'   \code{median} (their median), \code{score} (error over median),
#'   \code{pct}, \code{flag}, \code{n.periods}, \code{n.donors});
#'   \code{cohorts} (the same by adoption cohort, with \code{n.units});
#'   \code{yard} (each unit's donor errors); \code{summary}; \code{level}.
#' @seealso \code{\link{loading_check}}, \code{\link{fect_rank}}
#' @examples
#' data(simgsynth)
#' oc <- overlap_check(Y ~ D, data = simgsynth, index = c("id", "time"))
#' oc
#' @export
overlap_check <- function(formula, data, index, level = 0.90,
                          max.heldout = 200L, groups = 100L, seed = 1L) {
    if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
        stop("\"level\" must be a number between 0 and 1.", call. = FALSE)
    if (!is.numeric(max.heldout) || length(max.heldout) != 1L || max.heldout < 2)
        stop("\"max.heldout\" must be a number of at least 2.", call. = FALSE)
    if (!is.numeric(groups) || length(groups) != 1L || groups < 1)
        stop("\"groups\" must be a positive number.", call. = FALSE)
    w <- .gsc_wide(formula, data, index)
    out <- .overlap_panel(w$Y, w$D, ids = w$ids, level = level,
                          max.heldout = as.integer(max.heldout),
                          groups = as.integer(groups), seed = seed)
    out$units$onset <- w$time[out$units$onset]
    out$cohorts$onset <- w$time[out$cohorts$onset]
    out$call <- match.call()
    out
}

## The check on wide matrices: Y and D are T x N (NA allowed), units in
## columns. `onset` in the returned tables is the row of the first treated
## period. Internal: overlap_check() and the reference tests call it.
.overlap_panel <- function(Y, D, co = which(colSums(D, na.rm = TRUE) == 0),
                           ids = colnames(Y), level = 0.90, max.heldout = 200L,
                           groups = 100L, seed = 1L) {
    if (is.null(ids)) ids <- as.character(seq_len(ncol(Y)))
    pars <- list(level = level, max.heldout = max.heldout, groups = groups,
                 seed = seed, min.periods = 4L, min.donors = 5L, collapse = 1e-8)
    tr <- setdiff(which(colSums(D, na.rm = TRUE) > 0), co)
    if (!length(tr)) stop("No treated unit.", call. = FALSE)
    if (length(co) < pars$min.donors)
        warning("Fewer than five never-treated units: no treated unit is assessable.",
                call. = FALSE)
    rev <- vapply(tr, function(i) {
        d <- D[!is.na(D[, i]), i]
        any(diff(d) < 0)
    }, TRUE)
    if (any(rev))
        stop("The treatment switches off for some units; the check needs an ",
             "absorbing treatment.", call. = FALSE)
    onset <- vapply(tr, function(i) min(which(D[, i] == 1)), 0L)
    post <- vapply(tr, function(i) sum(D[, i] == 1 & !is.na(Y[, i]), na.rm = TRUE), 0)
    wt <- post / sum(post)
    cache <- new.env()
    cache$entries <- list()
    ux <- lapply(seq_along(tr), function(k)
        .of_unit(Y, tr[k], seq_len(onset[k] - 1L), co, cache, pars))
    cols <- c("status", "reason", "error", "yardstick", "median", "score",
              "flag", "pct", "n.periods", "n.donors")
    units <- do.call(rbind, lapply(seq_along(tr), function(k)
        data.frame(unit = ids[tr[k]], onset = onset[k], weight = wt[k],
                   ux[[k]][cols], stringsAsFactors = FALSE)))
    yard <- stats::setNames(lapply(ux, `[[`, "e.co"), ids[tr])
    cohorts <- do.call(rbind, lapply(sort(unique(onset)), function(o) {
        g <- tr[onset == o]
        x <- .of_cohort(Y, g, seq_len(o - 1L), co, seed = seed - 1L + o, pars)
        data.frame(onset = o, n.units = length(g), weight = sum(wt[onset == o]),
                   x[cols], stringsAsFactors = FALSE)
    }))
    rownames(units) <- rownames(cohorts) <- NULL
    ok.u <- units$status == "ok"
    ok.c <- cohorts$status == "ok"
    summary <- data.frame(
        n.tr = length(tr), n.co = length(co), n.assessable = sum(ok.u),
        n.flagged = sum(units$flag %in% TRUE),
        weight.assessable = sum(units$weight[ok.u]),
        weight.flagged = sum(units$weight[units$flag %in% TRUE]),
        cohorts = nrow(cohorts), cohort.weight.assessable = sum(cohorts$weight[ok.c]),
        cohort.weight.flagged = sum(cohorts$weight[cohorts$flag %in% TRUE]))
    structure(list(units = units, cohorts = cohorts, yard = yard,
                   summary = summary, level = level),
              class = "fect_overlap")
}

#' @rdname overlap_check
#' @param x an object returned by \code{overlap_check()}.
#' @param ... ignored.
#' @export
print.fect_overlap <- function(x, ...) {
    s <- x$summary
    cat("Overlap check without a rank (level ", format(x$level), ")\n", sep = "")
    cat(sprintf("Treated units: %d; never-treated donors: %d\n", s$n.tr, s$n.co))
    cat(sprintf("Assessable: %d of %d (%.0f%% of the effect's weight)\n",
                s$n.assessable, s$n.tr, 100 * s$weight.assessable))
    cat(sprintf("Flagged: %d (%.0f%% of the effect's weight)\n",
                s$n.flagged, 100 * s$weight.flagged))
    na <- x$units[x$units$status != "ok", , drop = FALSE]
    if (nrow(na)) {
        tab <- table(na$reason)
        cat("Not assessable: ", paste(sprintf("%s (%d)", names(tab), tab), collapse = ", "),
            "\n", sep = "")
    }
    fl <- x$units[x$units$flag %in% TRUE, , drop = FALSE]
    if (nrow(fl)) {
        fl <- fl[order(-fl$pct, -fl$score), c("unit", "onset", "weight", "error", "yardstick", "pct")]
        cat("\nFlagged units:\n")
        print(utils::head(fl, 20), row.names = FALSE, digits = 3)
        if (nrow(fl) > 20) cat("... and", nrow(fl) - 20, "more\n")
    }
    cat("\nCohorts: ", sum(x$cohorts$status == "ok"), " of ", nrow(x$cohorts),
        " assessable, ", sum(x$cohorts$flag %in% TRUE), " flagged\n", sep = "")
    invisible(x)
}

#' @rdname overlap_check
#' @export
plot.fect_overlap <- function(x, ...) {
    u <- x$units[x$units$status == "ok", , drop = FALSE]
    if (!nrow(u)) stop("No assessable treated unit to plot.", call. = FALSE)
    u <- u[order(u$pct, u$score), ]
    u$unit <- factor(u$unit, levels = u$unit)
    u$Flag <- ifelse(u$flag, "flagged", "not flagged")
    ggplot(u, aes(x = .data$pct, y = .data$unit)) +
        geom_vline(xintercept = x$level, linetype = "dashed", colour = "grey40") +
        geom_point(aes(size = .data$weight, colour = .data$Flag)) +
        scale_colour_manual(values = c(flagged = "#a8003b", `not flagged` = "#004a8c"),
                            name = NULL) +
        scale_size_continuous(name = "Weight", range = c(1, 4),
                              guide = if (length(unique(u$weight)) > 1) "legend" else "none") +
        scale_x_continuous(limits = c(0, 1)) +
        labs(x = "Percentile among the donors' errors", y = NULL,
             title = "Overlap check without a rank",
             subtitle = sprintf("%d of %d assessable units flagged at level %s",
                                sum(u$flag), nrow(u), format(x$level))) +
        theme_bw() +
        theme(axis.text.y = if (nrow(u) > 40) element_blank() else element_text())
}

## ---------------------------------------------------------------------------
## The check at the chosen rank
## ---------------------------------------------------------------------------

## The parts of a GSC fit that the rank rule and the loading check need, for
## the fit's treated units: the outcome net of the grand mean, time effects
## and covariates (R), the regressors of each unit's loading regression (Z:
## an intercept with unit effects, and the factors), and the pretreatment and
## post-treatment cells. Checks that a regression on all pretreatment periods
## reproduces fect's own unit effects and loadings.
.gsc_treated_parts <- function(fit) {
    tr <- fit$tr
    T <- nrow(fit$Y.dat)
    r <- if (is.null(fit$factor)) 0L else ncol(as.matrix(fit$factor))
    has.a <- fit$force %in% c(1, 3)
    F <- if (r > 0) as.matrix(fit$factor) else matrix(0, T, 0)
    A <- if (has.a) as.vector(fit$alpha.tr) else rep(0, length(tr))
    L <- if (r > 0) as.matrix(fit$lambda.tr) else matrix(0, length(tr), 0)
    M <- fit$Y.ct[, tr, drop = FALSE] - outer(rep(1, T), A) - F %*% t(L)
    R <- fit$Y.dat[, tr, drop = FALSE] - M
    obs <- fit$I.dat[, tr, drop = FALSE] == 1
    pre <- fit$D.dat[, tr, drop = FALSE] == 0 & obs
    post <- fit$D.dat[, tr, drop = FALSE] == 1 & obs
    Z <- cbind(if (has.a) 1, F)
    if (ncol(Z) > 0) {
        dev <- vapply(seq_along(tr), function(j) {
            p <- which(pre[, j])
            if (length(p) < ncol(Z)) return(0)
            b <- qr.solve(Z[p, , drop = FALSE], R[p, j])
            max(abs(b - c(if (has.a) A[j], L[j, ])))
        }, 0)
        scale <- max(1, stats::sd(fit$Y.dat[, tr, drop = FALSE][obs]), na.rm = TRUE)
        if (max(dev) > 1e-6 * scale)
            stop(sprintf("rebuilt loadings differ from the fit's (max %.2g) at r = %d",
                         max(dev), r), call. = FALSE)
    }
    list(ids = as.character(fit$id[tr]), R = R, pre = pre, post = post, Z = Z,
         has.a = has.a, r = r)
}

## Loadings of the controls from the same regression as the treated units':
## each control's outcome net of the grand mean, time effects and covariates
## on [1, F] over the periods in which at least half the treated units are
## untreated.
.gsc_control_loadings <- function(fit) {
    r <- if (is.null(fit$factor)) 0L else ncol(as.matrix(fit$factor))
    if (r == 0) return(NULL)
    co <- fit$co; tr <- fit$tr
    T <- nrow(fit$Y.dat)
    has.a <- fit$force %in% c(1, 3)
    F <- as.matrix(fit$factor)
    res <- as.matrix(fit$res)
    if (ncol(res) != length(co)) res <- res[, co, drop = FALSE]
    fitted <- fit$Y.dat[, co, drop = FALSE] - res
    A <- if (has.a) as.vector(fit$alpha.co) else rep(0, length(co))
    M <- fitted - outer(rep(1, T), A) - F %*% t(as.matrix(fit$lambda.co))
    R <- fit$Y.dat[, co, drop = FALSE] - M
    win <- rowMeans(fit$D.dat[, tr, drop = FALSE] == 0) >= 0.5
    Z <- cbind(if (has.a) 1, F)
    obs <- fit$I.dat[, co, drop = FALSE] == 1
    L <- t(vapply(seq_along(co), function(j) {
        p <- which(win & obs[, j])
        if (length(p) < ncol(Z) + 1) return(rep(NA_real_, r))
        b <- qr.solve(Z[p, , drop = FALSE], R[p, j])
        utils::tail(b, r)
    }, numeric(r)))
    if (r == 1) L <- matrix(L, ncol = 1)
    rownames(L) <- as.character(fit$id[co])
    L[stats::complete.cases(L), , drop = FALSE]
}

## Is x a convex combination of the rows of Lco? Exact least squares on the
## simplex; inside when the distance is zero up to rounding.
.gsc_in_hull <- function(Lco, x) {
    w <- .gsc_simplex_ls(x, t(Lco), dense.max = 0L)
    gap <- x - drop(crossprod(Lco, w))
    max(abs(gap)) <= 1e-7 * max(1, abs(Lco), abs(x))
}

## Squared distance from l to the controls' hull in the metric V:
## min over the simplex of (l - Lco'w)' V^{-1} (l - Lco'w).
.gsc_maha_hull <- function(Lco, l, V) {
    Vi <- solve(V)
    n <- nrow(Lco)
    if (n <= 400) {
        D <- Lco %*% Vi %*% t(Lco)
        D <- D + diag(1e-8 * max(1, mean(diag(D))), n)
        w <- quadprog::solve.QP(D, Lco %*% Vi %*% l, cbind(rep(1, n), diag(n)),
                                c(1, rep(0, n)), meq = 1)$solution
    } else {
        C <- chol(Vi)                     # Vi = C'C
        w <- .gsc_simplex_ls(drop(C %*% l), C %*% t(Lco), dense.max = 0L)
    }
    g <- l - drop(crossprod(Lco, w))
    drop(t(g) %*% Vi %*% g)
}

## Share of the treated post-treatment counterfactuals (D = 1) outside the
## range of the observed outcomes, overall and by treated unit.
.gsc_extrapolation <- function(fit) {
    tr <- fit$tr
    post <- fit$D.dat[, tr, drop = FALSE] == 1
    ct <- fit$Y.ct[, tr, drop = FALSE]
    rg <- range(fit$Y.dat[fit$I.dat == 1], na.rm = TRUE)
    out <- (ct < rg[1] | ct > rg[2]) & post
    out[!post] <- NA
    by.unit <- colMeans(out, na.rm = TRUE)
    list(share = mean(ct[post] < rg[1] | ct[post] > rg[2], na.rm = TRUE),
         by.unit = data.frame(unit = as.character(fit$id[tr]), share = unname(by.unit),
                              stringsAsFactors = FALSE),
         range = rg)
}

## Is this fit a GSC fit (factors from the never-treated units only)?
.gsc_is_gsc <- function(fit) {
    m <- fit$method
    tc <- fit$time.component.from
    identical(m, "gsynth") || (identical(m, "ife") && identical(tc, "nevertreated"))
}

## The loading check on a fit (no trimming): table, loadings and window.
.loading_table <- function(fit, level = 0.95) {
    Lco <- .gsc_control_loadings(fit)
    if (is.null(Lco)) return(NULL)
    Ltr <- as.matrix(fit$lambda.tr)
    ids <- as.character(fit$id[fit$tr])
    P <- .gsc_treated_parts(fit)
    r <- ncol(Ltr)
    inside <- vapply(seq_len(nrow(Ltr)), function(i) .gsc_in_hull(Lco, Ltr[i, ]), TRUE)
    d2 <- vapply(seq_len(nrow(Ltr)), function(i) {
        if (inside[i]) return(0)
        p <- which(P$pre[, i])
        Z <- P$Z[p, , drop = FALSE]
        k <- ncol(Z)
        if (length(p) <= k) return(NA_real_)
        e <- P$R[p, i] - drop(Z %*% qr.solve(Z, P$R[p, i]))
        V <- (sum(e^2) / (length(p) - k)) * solve(crossprod(Z))
        if (P$has.a) V <- V[-1, -1, drop = FALSE]
        .gsc_maha_hull(Lco, Ltr[i, ], V)
    }, 0)
    cutoff <- stats::qchisq(level, r)
    flagged <- !inside & is.finite(d2) & d2 > cutoff
    list(table = data.frame(unit = ids, in.hull = inside, d2 = d2, flagged = flagged,
                            stringsAsFactors = FALSE),
         Lco = Lco, Ltr = Ltr, cutoff = cutoff,
         window = sum(rowMeans(fit$D.dat[, fit$tr, drop = FALSE] == 0) >= 0.5))
}

#' Loading check at the chosen rank
#'
#' Checks whether each treated unit's estimated factor loadings lie inside
#' the region the controls' loadings enclose (their convex hull), allowing
#' for the noise in the treated unit's own estimate. This is the second
#' overlap check of the gsynth-revisit paper (Supplementary Material, "The
#' loading check at the chosen rank"). With \code{trim = TRUE}, the flagged
#' units are dropped, the rank is chosen again on the units kept (forward
#' rule, \code{\link{fect_rank}}), and the trimmed estimate is returned
#' beside the full one.
#'
#' Both groups' loadings come from the same regression: each unit's outcome,
#' net of the grand mean, time effects and covariates, on the estimated
#' factors (and an intercept when the model has unit effects). A treated
#' unit uses its own pretreatment periods (these are the fit's loadings); a
#' control uses the periods in which at least half the treated units are
#' untreated. A treated unit is inside when its loadings are a convex
#' combination of the controls' (computed by least squares on the simplex:
#' distance zero). Otherwise its squared distance to the hull is measured in
#' the metric of its loadings' sampling variance,
#' \eqn{\hat\sigma_i^2 (Z_i'Z_i)^{-1}}, with \eqn{\hat\sigma_i^2} the residual
#' variance of its loading regression with a degrees-of-freedom correction
#' (the intercept's row and column removed). It is flagged when the distance
#' exceeds the \code{level} quantile of a chi-square distribution with r
#' degrees of freedom. The cutoff treats the controls' loadings as fixed and
#' the treated unit's errors as uncorrelated over time, so it is a screen,
#' not a test.
#'
#' The output also gives the extrapolation share: the share of the treated
#' post-treatment counterfactuals that lie outside the range of the observed
#' outcomes.
#'
#' Trimming changes the estimand: the trimmed estimate is the effect for the
#' treated units the controls can represent. Report it beside the full
#' estimate, not in its place.
#'
#' @param fit a fit of \code{fect()} with \code{method = "gsynth"} (or
#'   \code{method = "ife"} with \code{time.component.from = "nevertreated"})
#'   and at least one factor.
#' @param level the chi-square quantile for the flag.
#' @param trim if \code{TRUE}, drop the flagged units, choose the rank again
#'   on the units kept and refit with the settings of \code{fit}.
#' @param data the data of \code{fit}; needed only for \code{trim = TRUE}
#'   when the data cannot be found from the call that made \code{fit}.
#' @param r the candidate ranks for choosing the rank again.
#' @return An object of class \code{"fect_loading"}: a list with
#'   \code{table} (one row per treated unit: \code{unit}, \code{in.hull},
#'   \code{d2} (0 inside the hull), \code{flagged}), \code{r}, \code{level},
#'   \code{cutoff}, \code{n.flagged}, \code{extrapolation} (\code{share},
#'   \code{by.unit}, \code{range}), \code{Lco} and \code{Ltr} (the controls'
#'   and treated units' loadings), \code{window} (the controls' loading
#'   periods), and with \code{trim = TRUE} \code{trim}: \code{status},
#'   \code{out} (the units dropped), \code{r} (the rank chosen again),
#'   \code{rank} (its \code{fect_rank} object), \code{est} (the trimmed
#'   estimate: \code{ATT}, \code{S.E.}, \code{CI.lower}, \code{CI.upper},
#'   \code{p.value}, \code{n.tr}), \code{full} (the same for \code{fit}),
#'   \code{n.flagged.after} and \code{fit} (the trimmed fit).
#' @seealso \code{\link{overlap_check}}, \code{\link{fect_rank}}
#' @examples
#' data(simgsynth)
#' fit <- fect(Y ~ D, data = simgsynth, index = c("id", "time"),
#'             method = "gsynth", r = 2, se = FALSE)
#' loading_check(fit)
#' @export
loading_check <- function(fit, level = 0.95, trim = FALSE, data = NULL, r = 0:5) {
    if (!inherits(fit, "fect"))
        stop("\"fit\" must be a fit of fect().", call. = FALSE)
    if (!.gsc_is_gsc(fit))
        stop("loading_check() needs a GSC fit: method = \"gsynth\" (or method = \"ife\" with ",
             "time.component.from = \"nevertreated\").", call. = FALSE)
    if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
        stop("\"level\" must be a number between 0 and 1.", call. = FALSE)
    if (!is.null(fit$hasRevs) && isTRUE(fit$hasRevs == 1))
        stop("loading_check() needs an absorbing treatment.", call. = FALSE)
    rk <- if (is.null(fit$factor)) 0L else ncol(as.matrix(fit$factor))
    if (rk == 0)
        stop("The fit has no factors (r = 0): there are no loadings to check.", call. = FALSE)
    h <- .loading_table(fit, level)
    out <- list(table = h$table, r = rk, level = level, cutoff = h$cutoff,
                n.flagged = sum(h$table$flagged), extrapolation = .gsc_extrapolation(fit),
                Lco = h$Lco, Ltr = h$Ltr, window = h$window)
    if (isTRUE(trim)) out$trim <- .loading_trim(fit, h, data, r, level, parent.frame())
    structure(out, class = "fect_loading")
}

.est_row <- function(f) {
    a <- f$est.avg
    if (is.null(a)) {
        return(data.frame(ATT = unname(f$att.avg), S.E. = NA_real_, CI.lower = NA_real_,
                          CI.upper = NA_real_, p.value = NA_real_, n.tr = f$Ntr))
    }
    data.frame(ATT = a[1, "ATT.avg"], S.E. = a[1, "S.E."], CI.lower = a[1, "CI.lower"],
               CI.upper = a[1, "CI.upper"], p.value = a[1, "p.value"], n.tr = f$Ntr,
               row.names = NULL)
}

## The data and the call of a fit, for refits.
.fit_call_data <- function(fit, data, env) {
    cl <- fit$call
    if (is.null(cl)) stop("The fit has no call; refit it with fect().", call. = FALSE)
    if (is.null(data)) {
        data <- tryCatch(eval(cl$data, env), error = function(e) NULL)
        if (!is.data.frame(data))
            stop("Cannot find the data of the fit; pass them as \"data\".", call. = FALSE)
    }
    if (is.null(fit$index)) stop("The fit has no index.", call. = FALSE)
    list(call = cl, data = data)
}

.loading_trim <- function(fit, h, data, r, level, env) {
    out.ids <- h$table$unit[h$table$flagged]
    base <- list(status = "", r.full = ncol(h$Ltr), n.tr = nrow(h$table),
                 n.out = length(out.ids), out = out.ids, full = .est_row(fit))
    if (!length(out.ids)) {
        base$status <- "no treated unit significantly outside the hull"
        return(base)
    }
    if (length(out.ids) == nrow(h$table)) {
        base$status <- "every treated unit is significantly outside the hull"
        return(base)
    }
    cd <- .fit_call_data(fit, data, env)
    unit.var <- fit$index[1]
    d <- cd$data[!as.character(cd$data[[unit.var]]) %in% out.ids, , drop = FALSE]
    rk <- .rank_from_call(cd$call, d, r, env)
    r2 <- if (is.na(rk$r.fwd)) base$r.full else as.integer(rk$r.fwd)
    base$r <- r2
    base$rank <- rk
    cl <- cd$call
    cl[[1]] <- fect                       # the call may name fect.formula(), not exported
    cl$data <- d
    cl$r <- r2
    cl$CV <- FALSE
    cl$cv.method <- NULL
    f <- tryCatch(eval(cl, env), error = function(e) e)
    if (inherits(f, "error")) {
        base$status <- paste("trimmed fit failed:", conditionMessage(f))
        return(base)
    }
    base$n.flagged.after <- if (r2 >= 1) sum(.loading_table(f, level)$table$flagged) else NA_integer_
    base$n.outside.exact <- sum(!h$table$in.hull)
    base$est <- .est_row(f)
    base$fit <- f
    base$status <- "ok"
    base
}

#' @rdname loading_check
#' @param x an object returned by \code{loading_check()}.
#' @param ... ignored.
#' @export
print.fect_loading <- function(x, ...) {
    t <- x$table
    cat(sprintf("Loading check at r = %d (chi-square level %s, cutoff %.2f)\n",
                x$r, format(x$level), x$cutoff))
    cat(sprintf("Treated units: %d; inside the controls' hull: %d; flagged: %d\n",
                nrow(t), sum(t$in.hull), x$n.flagged))
    cat(sprintf("Extrapolation: %.1f%% of treated post-treatment counterfactuals outside the observed outcome range\n",
                100 * x$extrapolation$share))
    if (x$n.flagged > 0) {
        fl <- t[t$flagged, c("unit", "d2")]
        fl <- fl[order(-fl$d2), ]
        cat("\nFlagged units:\n")
        print(utils::head(fl, 20), row.names = FALSE, digits = 3)
        if (nrow(fl) > 20) cat("... and", nrow(fl) - 20, "more\n")
    }
    if (!is.null(x$trim)) {
        tr <- x$trim
        cat("\nTrimmed: ", tr$status, "\n", sep = "")
        if (identical(tr$status, "ok")) {
            e <- rbind(full = tr$full, trimmed = tr$est)
            e$r <- c(tr$r.full, tr$r)
            print(e, digits = 3)
            cat("The trimmed estimate is the effect for the treated units the controls can",
                "represent; report it beside the full one.\n")
        }
    }
    invisible(x)
}

#' @rdname loading_check
#' @export
plot.fect_loading <- function(x, ...) {
    co <- x$Lco
    tr <- x$Ltr
    lab <- x$table$unit
    fl <- x$table$flagged
    st <- ifelse(fl, "flagged", ifelse(x$table$in.hull, "inside", "outside, not flagged"))
    cols <- c(flagged = "#a8003b", inside = "#004a8c", `outside, not flagged` = "grey40")
    if (ncol(co) == 1) {
        dt <- data.frame(l = tr[, 1], st = st, lab = lab)
        p <- ggplot() +
            geom_jitter(aes(x = co[, 1], y = 0), width = 0, height = 0.05,
                        colour = "grey60", size = 1.5) +
            geom_point(data = dt, aes(x = .data$l, y = 0.2, colour = .data$st),
                       shape = 17, size = 2.6) +
            geom_text(data = dt[fl, , drop = FALSE], aes(x = .data$l, y = 0.2, label = .data$lab),
                      vjust = -1, size = 3) +
            scale_colour_manual(values = cols, name = "Treated") +
            scale_y_continuous(breaks = c(0, 0.2), labels = c("controls", "treated")) +
            labs(x = "Loading", y = NULL)
    } else {
        hh <- grDevices::chull(co[, 1], co[, 2])
        dt <- data.frame(l1 = tr[, 1], l2 = tr[, 2], st = st, lab = lab)
        p <- ggplot() +
            geom_polygon(aes(x = co[hh, 1], y = co[hh, 2]), fill = "grey80", alpha = 0.4,
                         colour = "grey60") +
            geom_point(aes(x = co[, 1], y = co[, 2]), colour = "grey50", size = 1.5) +
            geom_point(data = dt, aes(x = .data$l1, y = .data$l2, colour = .data$st),
                       shape = 17, size = 2.6) +
            geom_text(data = dt[fl, , drop = FALSE],
                      aes(x = .data$l1, y = .data$l2, label = .data$lab), vjust = -1, size = 3) +
            scale_colour_manual(values = cols, name = "Treated") +
            labs(x = "Loading on factor 1", y = "Loading on factor 2")
    }
    p + labs(title = sprintf("Loadings at r = %d", x$r),
             subtitle = paste0("Grey: controls (loadings from their first ", x$window,
                               " periods).\nFlagged: outside the hull by more than the noise allows",
                               if (ncol(co) > 2) "; first two factors shown" else "")) +
        theme_bw()
}
