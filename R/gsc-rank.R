## Rank choice from the treated units' pretreatment outcomes (fect 2.5.0):
## the forward rule and the leave-one-out rule of the gsynth-revisit paper
## (Algorithm "Rank choice from the treated units' pretreatment outcomes" and
## its extensions to staggered adoption and missing cells). The reference
## implementation is the paper's replication code (R/rank_rule.R,
## treated_rank() and window_minimum()).
##
## For each candidate rank r, GSC is fitted once. The controls' fit (grand
## mean, time effects, covariate coefficients, factors) is the fit's. Each
## treated unit's untreated outcome is m + alpha_i + F lambda_i; the rule
## refits the unit's intercept (with unit effects) and loadings on part of
## its pretreatment periods and predicts the rest:
##   forward:       score period c + 1 from periods 1..c, for the last cutoffs,
##                  where the shortest training set has one more period than
##                  the largest rank needs;
##   leave one out: hold out each pretreatment period in turn and train on the
##                  unit's other pretreatment periods (Xu 2017, Algorithm 1).
## The score is the squared prediction error averaged over a unit's scored
## periods, averaged over units with weights equal to each unit's share of
## the treated post-treatment cells (its weight in the average effect). The
## chosen rank is the smallest minimizer. The candidate ranks of a window go
## up to the largest rank (at most the largest rank fitted) at which the
## units that can be scored carry at least 80% of that weight; units with too
## few pretreatment periods are left out of the score, not the estimate.

.RANK_MIN_SHARE <- 0.8
.RANK_MIN_PREDICTIONS <- 8L

## Scored tasks of one unit: list of (train, score) period indices. `p` is the
## unit's pretreatment periods in time order; k = number of regressors at the
## window's largest rank.
.rank_unit_tasks <- function(p, k, window) {
    n <- length(p)
    if (window == "loo") {
        if (n - 1 < k) return(list())
        lapply(seq_len(n), function(s) list(train = p[-s], score = p[s]))
    } else {
        K <- n - k - 1
        if (K < 1) return(list())
        lapply(n - seq_len(K), function(c0) list(train = p[seq_len(c0)], score = p[c0 + 1]))
    }
}

## Scores of both windows at every rank in `fits` (a list of GSC fits named
## by rank, consecutive from 0). Returns the score table, both picks and the
## counts of scored units.
.rank_scores <- function(fits) {
    fits <- fits[!vapply(fits, is.null, TRUE)]
    ranks <- sort(as.integer(names(fits)))
    if (!length(ranks)) return(NULL)
    parts <- lapply(fits[as.character(ranks)], function(f)
        tryCatch(.gsc_treated_parts(f), error = function(e) e))
    bad <- vapply(parts, inherits, TRUE, "error")
    dropped <- ""
    if (any(bad)) {
        k <- which(bad)[1]
        dropped <- sprintf("; r >= %d left out (%s)", ranks[k], conditionMessage(parts[[k]]))
        if (k == 1) stop(conditionMessage(parts[[k]]), call. = FALSE)
        ranks <- ranks[seq_len(k - 1)]
        parts <- parts[seq_len(k - 1)]
    }
    ids <- Reduce(intersect, lapply(parts, `[[`, "ids"))
    has.a <- parts[[1]]$has.a
    jj0 <- match(ids, parts[[1]]$ids)
    w.att <- colSums(parts[[1]]$post[, jj0, drop = FALSE])
    w.att <- w.att / sum(w.att)
    pre.j <- lapply(jj0, function(j) which(parts[[1]]$pre[, j]))
    ok.at <- function(rb, window)
        vapply(pre.j, function(p) length(.rank_unit_tasks(p, rb + has.a, window)) > 0, TRUE)
    rbar.of <- function(window) {
        cand <- ranks[vapply(ranks, function(rb)
            sum(w.att[ok.at(rb, window)]) >= .RANK_MIN_SHARE, TRUE)]
        if (length(cand)) max(cand) else NA_integer_
    }
    rbars <- c(loo = rbar.of("loo"), fwd = rbar.of("fwd"))
    if (all(is.na(rbars)))
        return(list(table = NULL, r.loo = NA_integer_, r.fwd = NA_integer_,
                    note = "no rank can be scored: too few pretreatment periods"))
    rbar <- max(rbars, na.rm = TRUE)
    ranks <- ranks[ranks <= rbar]
    score <- function(window) {
        rb <- rbars[[window]]
        if (is.na(rb)) return(list(s = rep(NA_real_, length(ranks)), n = 0L, share = 0))
        k <- rb + has.a
        use <- ids[ok.at(rb, window)]
        wu <- w.att[match(use, ids)]
        wu <- wu / sum(wu)
        s <- vapply(ranks, function(r) {
            if (r > rb) return(NA_real_)
            P <- parts[[as.character(r)]]
            tot <- 0
            for (j in seq_along(use)) {
                jj <- match(use[j], P$ids)
                p <- which(P$pre[, jj])
                tk <- .rank_unit_tasks(p, k, window)
                e <- vapply(tk, function(t) {
                    if (ncol(P$Z) == 0) return(P$R[t$score, jj])
                    Zt <- P$Z[t$train, , drop = FALSE]
                    if (qr(Zt)$rank < ncol(Zt)) return(NA_real_)
                    P$R[t$score, jj] - drop(P$Z[t$score, , drop = FALSE] %*%
                                                qr.solve(Zt, P$R[t$train, jj]))
                }, 0)
                if (anyNA(e)) return(Inf)        # rank deficient: never chosen; reported
                tot <- tot + wu[j] * mean(e^2)
            }
            tot
        }, 0)
        list(s = s, n = length(use), share = sum(w.att[match(use, ids)]))
    }
    lo <- score("loo")
    fw <- score("fwd")
    pick <- function(s) if (all(!is.finite(s))) NA_integer_ else ranks[which.min(s)]
    list(table = data.frame(r = ranks, fwd = fw$s, loo = lo$s),
         r.fwd = pick(fw$s), r.loo = pick(lo$s), rbar = rbar,
         rbar.fwd = rbars[["fwd"]], rbar.loo = rbars[["loo"]],
         n.tr = length(ids), n.fwd = fw$n, n.loo = lo$n,
         share.fwd = fw$share, share.loo = lo$share,
         note = paste0(sprintf(paste0("forward: ranks 0-%s, scores %d of %d treated units (%.0f%% of the ",
                                      "effect's weight); leave one out: ranks 0-%s, %d (%.0f%%)"),
                               rbars[["fwd"]], fw$n, length(ids), 100 * fw$share,
                               rbars[["loo"]], lo$n, 100 * lo$share), dropped))
}

## The forward rule's minimum window: at least 8 one-step forward predictions
## per treated unit at rank `rmax` (15 pretreatment periods with unit effects
## and rmax = 5), for treated units carrying at least 80% of the effect's
## weight. Any fit of the case will do: the windows do not depend on the rank.
.rank_window <- function(fit, rmax = 5L) {
    P <- .gsc_treated_parts(fit)
    n.pre <- colSums(P$pre)
    w <- colSums(P$post)
    w <- w / sum(w)
    preds <- pmax(n.pre - (rmax + P$has.a) - 1, 0)
    share <- sum(w[preds >= .RANK_MIN_PREDICTIONS])
    list(n.pre.min = min(n.pre), n.pre.med = stats::median(n.pre),
         preds.med = stats::median(preds), share8 = share,
         all8 = all(preds >= .RANK_MIN_PREDICTIONS), meets = share >= .RANK_MIN_SHARE)
}

## Candidate ranks from fect's `r` argument: a vector of ranks, or a range
## c(lower, upper) as in fect(CV = TRUE).
.rank_candidates <- function(r) {
    if (is.null(r)) return(0:5)
    if (!is.numeric(r) || any(!is.finite(r)) || any(r < 0) || any(r != round(r)))
        stop("\"r\" must be non-negative whole numbers.", call. = FALSE)
    r <- as.integer(r)
    if (length(r) == 1L) return(0:r)
    if (length(r) == 2L) return(seq(min(r), max(r)))
    sort(unique(r))
}

## Rank choice from point fits made by `fit_at(r)`, at every candidate rank.
## A rank fect cannot fit is left out with every rank above it (the
## candidates must be consecutive from the lowest).
.rank_from_fitter <- function(fit_at, ranks, keep.fits = FALSE) {
    if (ranks[1] != 0L || any(diff(ranks) != 1L))
        stop("The candidate ranks must run from 0 without gaps, e.g. r = 0:5.", call. = FALSE)
    fits <- vector("list", length(ranks))
    names(fits) <- ranks
    failed <- integer(0)
    why <- character(0)
    for (k in seq_along(ranks)) {
        f <- tryCatch(fit_at(ranks[k]), error = function(e) e)
        if (inherits(f, "error")) {
            failed <- ranks[k]
            why <- conditionMessage(f)
            break
        }
        fits[[k]] <- f
    }
    if (length(failed)) fits <- fits[seq_len(match(failed, ranks) - 1L)]
    out <- tryCatch(.rank_scores(fits), error = function(e)
        list(table = NULL, r.loo = NA_integer_, r.fwd = NA_integer_,
             note = paste("rank rule failed:", conditionMessage(e))))
    if (is.null(out)) out <- list(table = NULL, r.loo = NA_integer_, r.fwd = NA_integer_,
                                  note = "no rank could be fitted")
    if (length(failed)) {
        out$note <- sprintf("%s; fect could not fit r = %d (%s)", out$note, failed, why)
        out$failed <- data.frame(r = failed, reason = why, stringsAsFactors = FALSE)
    }
    out$window <- if (length(fits) && !is.null(fits[[1]]))
        tryCatch(.rank_window(fits[[1]], max(ranks)), error = function(e) NULL)
    ## estimates and extrapolation at every rank fitted
    if (length(fits)) {
        est <- do.call(rbind, lapply(names(fits), function(k) {
            f <- fits[[k]]
            if (is.null(f)) return(NULL)
            e <- .est_row(f)
            e$extrapolation <- tryCatch(.gsc_extrapolation(f)$share, error = function(e) NA_real_)
            cbind(r = as.integer(k), e)
        }))
        out$estimates <- est
        if (!is.null(out$table)) {
            out$table <- merge(out$table, est, by = "r", all.x = TRUE, sort = TRUE)
        }
    }
    if (keep.fits) out$fits <- fits
    out
}

## The rank rule on the data of a call (loading_check(trim = TRUE)): point
## fits with the call's settings, no standard errors.
.rank_from_call <- function(cl, data, r, env) {
    ranks <- .rank_candidates(r)
    fit_at <- function(k) {
        c2 <- cl
        c2[[1]] <- fect                   # the call may name fect.formula(), not exported
        c2$data <- data
        c2$r <- k
        c2$CV <- FALSE
        c2$se <- FALSE
        c2$parallel <- FALSE
        c2$cv.method <- NULL
        suppressWarnings(suppressMessages(eval(c2, env)))
    }
    structure(.rank_from_fitter(fit_at, ranks), class = "fect_rank")
}

#' Rank choice from the treated units' pretreatment outcomes
#'
#' Chooses the number of factors for the generalized synthetic control
#' (GSC) estimator by how well each candidate rank predicts the treated
#' units' own pretreatment outcomes: the forward rule and the leave-one-out
#' rule of the gsynth-revisit paper (Algorithm "Rank choice from the treated
#' units' pretreatment outcomes"). GSC (\code{method = "gsynth"}) is fitted
#' once at every candidate rank.
#'
#' For each rank, the controls' fit (grand mean, time effects, covariate
#' coefficients and factors) is the fit's. Each treated unit's intercept
#' (with unit effects) and loadings are refitted on part of its
#' pretreatment periods and used to predict another:
#' \describe{
#'   \item{forward}{predicts period c + 1 from periods 1 to c, for the last
#'     cutoffs; the shortest training set has one more period than the
#'     largest candidate rank needs. Trains only on earlier periods, as the
#'     post-treatment prediction does.}
#'   \item{leave one out}{holds out each pretreatment period in turn and
#'     trains on the others (Xu 2017, Algorithm 1). Under serially
#'     correlated errors the later periods carry information about the
#'     held-out error, so it leans toward more factors.}
#' }
#' Each score averages the squared prediction errors over a unit's scored
#' periods, then over units, weighting each unit by its share of the treated
#' post-treatment cells (its weight in the average effect); with staggered
#' adoption the tasks are in each unit's own event time. The chosen rank is
#' the smallest minimizer. A window's candidate ranks go up to the largest
#' rank at which the units with enough pretreatment periods to be scored
#' carry at least 80\% of that weight; the other units are left out of the
#' score (not of the estimate), and \code{note} says how many were scored. A
#' rank at which some unit's training regression is rank deficient gets an
#' infinite score.
#'
#' \code{window} reports whether the forward window is long enough: the
#' paper asks for at least 8 one-step forward predictions per treated unit
#' at the largest candidate rank (15 pretreatment periods with unit effects
#' and ranks up to 5), for units carrying at least 80\% of the weight.
#'
#' The interval reported at the chosen rank does not account for the
#' choice. The rule does not establish a flat pretrend or recover a
#' structural rank.
#'
#' @param formula,data,index as in \code{\link{fect}}.
#' @param r the candidate ranks, from 0 without gaps (default 0:5); a pair
#'   \code{c(0, 5)} is read as a range, as in \code{fect(CV = TRUE)}.
#' @param se if \code{TRUE}, every fit has standard errors (by default the
#'   parametric bootstrap, Algorithm 2 of Xu 2017), so the table gives the
#'   estimate and its interval at every rank. The choice itself needs none.
#' @param keep.fits if \code{TRUE}, return the fits.
#' @param ... further arguments to \code{fect()} (for example \code{force},
#'   \code{min.T0}, \code{nboots}, \code{vartype}, \code{seed}). \code{CV} and
#'   \code{method} are set by \code{fect_rank()}.
#' @return An object of class \code{"fect_rank"}: a list with \code{table}
#'   (one row per rank: \code{r}, the scores \code{fwd} and \code{loo},
#'   \code{ATT}, \code{S.E.}, \code{CI.lower}, \code{CI.upper},
#'   \code{p.value}, \code{n.tr} and \code{extrapolation}, the share of
#'   treated post-treatment counterfactuals outside the observed outcome
#'   range), \code{r.fwd} and \code{r.loo} (the two picks), \code{rbar},
#'   \code{rbar.fwd} and \code{rbar.loo} (the largest candidate rank of each
#'   window), \code{n.tr}, \code{n.fwd}, \code{n.loo}, \code{share.fwd},
#'   \code{share.loo} (scored units and their weight), \code{window}
#'   (\code{meets}: the forward window is long enough; \code{share8},
#'   \code{all8}, \code{n.pre.min}, \code{n.pre.med}, \code{preds.med}),
#'   \code{note}, and \code{failed} when fect could not fit a rank.
#' @seealso \code{\link{loading_check}}, \code{\link{overlap_check}};
#'   \code{fect(..., CV = TRUE, cv.method = "forward")} fits at the forward
#'   rule's rank directly.
#' @examples
#' data(simgsynth)
#' rk <- fect_rank(Y ~ D, data = simgsynth, index = c("id", "time"), r = 0:3)
#' rk
#' @export
fect_rank <- function(formula = NULL, data, index, r = 0:5, se = FALSE,
                      keep.fits = FALSE, ...) {
    dots <- list(...)
    for (a in c("CV", "method", "cv.method"))
        if (a %in% names(dots))
            stop("fect_rank() sets \"", a, "\" itself; leave it out.", call. = FALSE)
    if (!is.null(dots$time.component.from) && !identical(dots$time.component.from, "nevertreated"))
        stop("fect_rank() fits GSC, whose factors come from the never-treated units.", call. = FALSE)
    dots$time.component.from <- NULL
    ranks <- .rank_candidates(r)
    fit_at <- function(k) {
        args <- c(list(formula = formula, data = data, index = index, method = "gsynth",
                       CV = FALSE, r = k, se = se), dots)
        suppressWarnings(suppressMessages(do.call(fect, args)))
    }
    out <- .rank_from_fitter(fit_at, ranks, keep.fits = keep.fits)
    out$call <- match.call()
    structure(out, class = "fect_rank")
}

#' @rdname fect_rank
#' @param x an object returned by \code{fect_rank()}.
#' @export
print.fect_rank <- function(x, ...) {
    cat("Rank choice from the treated units' pretreatment outcomes\n")
    cat("Forward rule: r = ", if (is.na(x$r.fwd)) "none" else x$r.fwd,
        "; leave one out: r = ", if (is.na(x$r.loo)) "none" else x$r.loo, "\n", sep = "")
    if (!is.null(x$note)) cat(x$note, "\n", sep = "")
    if (!is.null(x$window)) {
        w <- x$window
        cat(sprintf(paste0("Forward window: %s (units with at least 8 forward predictions at the ",
                           "largest rank carry %.0f%% of the weight; median pretreatment periods %s)\n"),
                    if (isTRUE(w$meets)) "long enough" else "SHORT, read the forward pick with care",
                    100 * w$share8, format(w$n.pre.med)))
    }
    if (!is.null(x$table)) {
        cat("\n")
        tb <- x$table
        print(tb, row.names = FALSE, digits = 4)
    }
    invisible(x)
}

## fect(CV = TRUE, cv.method = "forward"): choose the rank with the forward
## rule, then fit at that rank. `cl` is the fect.default call; `env` the
## frame it was made in. Point fits use the call's settings without standard
## errors.
.fect_forward_cv <- function(cl, env, args) {
    method <- args$method
    tcf <- args$time.component.from
    if (!(identical(method, "gsynth") || (identical(method, "ife") && identical(tcf, "nevertreated"))))
        stop("cv.method = \"forward\" chooses the rank of GSC: use method = \"gsynth\" ",
             "(or method = \"ife\" with time.component.from = \"nevertreated\").", call. = FALSE)
    if (!is.null(args$lambda) || isTRUE(args$binary))
        stop("cv.method = \"forward\" is for a continuous outcome and a factor model.", call. = FALSE)
    ranks <- .rank_candidates(args$r)
    cl[[1]] <- fect.default
    fit_at <- function(k) {
        c2 <- cl
        c2$r <- k
        c2$CV <- FALSE
        c2$se <- FALSE
        c2$parallel <- FALSE
        c2$cv.method <- "rolling"
        suppressWarnings(suppressMessages(eval(c2, env)))
    }
    rk <- structure(.rank_from_fitter(fit_at, ranks), class = "fect_rank")
    if (is.na(rk$r.fwd))
        stop("cv.method = \"forward\" could not choose a rank: ", rk$note, call. = FALSE)
    message("Forward rule: r = ", rk$r.fwd, " (", rk$note, ")")
    if (!is.null(rk$window) && !isTRUE(rk$window$meets))
        message("The forward window is short: fewer than 80% of the effect's weight has 8 ",
                "forward predictions at the largest rank.")
    c2 <- cl
    c2$r <- rk$r.fwd
    c2$CV <- FALSE
    c2$cv.method <- "rolling"
    out <- eval(c2, env)
    out$r.cv <- rk$r.fwd
    out$CV <- TRUE
    out$cv.method <- "forward"
    out$rank.forward <- rk
    out
}
