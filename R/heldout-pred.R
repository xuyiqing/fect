## The model's prediction at held-out cells, for fect_mspe() and
## r.cv.rolling().
##
## Both scorers hide a cell by setting its outcome to NA and refitting with
## fect(). fect() drops a row whose outcome is missing and stores 0 for that
## cell's covariates, so the refit's Y.ct.full there leaves out X times beta.
## On the never-treated path (method = "gsynth", or "ife" or "cfe" with
## time.component.from = "nevertreated") Y.ct.full is not the fitted value
## at all: its control columns hold only the factor part F lambda', and a
## treated unit's hidden cell holds 0. Before 2.4.7 the scorers read
## Y.ct.full there.
##
## .fect_heldout_pred() returns base + X beta, where X holds the cells'
## covariates from the data the cells were hidden in (newdata) and base is
## the refit's fitted value at the cell with the covariates at 0:
##   - never-treated path, control unit: fit$est$fit, the fitted values of
##     the model estimated on the never-treated units (intercept, unit and
##     time effects, factors and, for cfe, its other terms);
##   - never-treated path, treated unit: mu + alpha_i + xi_t + F_t lambda_i,
##     with the unit's own loadings and unit effect (NA for a cfe fit with
##     Z, Q or extra fixed effects: those terms of a treated unit are not
##     stored in the fit);
##   - otherwise: fit$Y.ct.full, the fitted value at the missing cell.
## The fit is not changed.
##
## fit: a fect() output (the refit). rr, cc: the cells' row (time) and
## column (unit) positions in the fit's TT x N matrices, no NA. newdata: a
## data frame with one row per cell that has the covariates fit$X.
.fect_heldout_pred <- function(fit, rr, cc, newdata) {
    yct <- if (!is.null(fit$Y.ct.full)) fit$Y.ct.full else fit$Y.ct
    pred <- yct[cbind(rr, cc)]
    nt.path <- identical(fit$time.component.from, "nevertreated") &&
        !identical(fit$method, "fe") && !is.null(fit$co)
    if (nt.path) {
        j <- match(cc, fit$co)
        in.co <- !is.na(j)
        if (any(in.co)) {
            est.fit <- fit$est$fit
            if (is.matrix(est.fit) && ncol(est.fit) == length(fit$co) &&
                nrow(est.fit) == nrow(yct)) {
                pred[in.co] <- est.fit[cbind(rr[in.co], j[in.co])]
            } else {
                ## Not reached by the scorers: a hidden control cell makes
                ## the never-treated panel unbalanced, and fect() then
                ## stores the fitted values in est$fit.
                pred[in.co] <- NA_real_
            }
        }
        k <- match(cc, fit$tr)
        in.tr <- !is.na(k)
        cfe.terms <- identical(fit$method, "cfe") &&
            (length(fit$gamma) > 0 || length(fit$kappa) > 0 ||
             length(fit$index) > 2)
        if (any(in.tr) && cfe.terms) {
            ## A treated unit's Z, Q and extra fixed-effect terms are not
            ## stored in the fit, so its prediction cannot be rebuilt:
            ## the cell is not scored. (fect_mspe() hides treated cells
            ## of a never-treated fit only when a not-yet-treated fit
            ## comes first in its list.)
            pred[in.tr] <- NA_real_
        } else if (any(in.tr)) {
            rt <- rr[in.tr]
            kt <- k[in.tr]
            b <- rep(fit$mu, length(rt))
            if (!is.null(fit$alpha.tr)) {
                b <- b + as.matrix(fit$alpha.tr)[kt, 1]
            }
            if (!is.null(fit$xi)) {
                b <- b + as.matrix(fit$xi)[rt, 1]
            }
            if (!is.null(fit$factor) && !is.null(fit$lambda.tr) &&
                NCOL(fit$factor) > 0) {
                b <- b + rowSums(as.matrix(fit$factor)[rt, , drop = FALSE] *
                                 as.matrix(fit$lambda.tr)[kt, , drop = FALSE])
            }
            pred[in.tr] <- b
        }
    }
    pred + .fect_heldout_xbeta(fit, newdata)
}

## X times beta at the held-out cells. A coefficient fect() could not
## estimate (NA) adds nothing, as in the refit's own fitted values.
.fect_heldout_xbeta <- function(fit, newdata) {
    xn <- fit$X
    if (length(xn) == 0 || is.null(fit$beta)) {
        return(0)
    }
    b <- as.numeric(fit$beta)
    if (length(b) != length(xn)) {
        stop("Internal error: fit$beta has ", length(b), " entries for ",
             length(xn), " covariates.", call. = FALSE)
    }
    b[!is.finite(b)] <- 0
    ## column by column with [[ ]], which reads a data.frame, a tibble and
    ## a data.table the same way
    xb <- numeric(nrow(newdata))
    for (j in seq_along(xn)) {
        xb <- xb + as.numeric(newdata[[xn[j]]]) * b[j]
    }
    xb
}
