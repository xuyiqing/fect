## -----------------------------------------------------------------------------
## conformal.R
##
## Cross-sectional conformal inference for separated, symmetric counterfactual
## estimators (vartype = "conformal").  The treated counterfactual prediction
## error is calibrated against the leave-one-control-out prediction errors of the
## donors and read off as a RANK, rather than resampled and normal-wrapped
## (parametric bootstrap) or quantile-resampled (nonparametric bootstrap).
##
## Design: statsclaw-workspace/fect/runs/REQ-conformal/spec.md
## Decisions (2026-06-08): jackknife+ default (>= 1 - 2 alpha); scores
## "studentized" (default), "ratio" (Abadie post/pre RMSPE), "meanabs", "rmse";
## staggered adoption via valid_controls(predictive = "notyettreated").
##
## This file holds the estimator-AGNOSTIC core: scores, the conformal p-value,
## and interval inversion.  The calibration loop (leave-one-control-out via
## impute_Y0) and the vartype dispatch live in boot.R / default.R and feed this
## core a control-by-period residual matrix plus the treated gap path.
## -----------------------------------------------------------------------------

## ---- 1. nonconformity score -------------------------------------------------
## e.post : numeric, a unit's post-treatment residual (prediction-error) path
## e.pre  : numeric, the same unit's pre-treatment residual path (for normaliz.)
## type   : "studentized" | "ratio" | "meanabs" | "rmse"
## Returns a scalar.  Validity (exchangeability) is identical across types; the
## choice affects interval width / adaptivity only.
.conformal_score <- function(e.post, e.pre, type = "studentized") {
  e.post <- e.post[is.finite(e.post)]
  if (length(e.post) == 0L) return(NA_real_)
  rms <- function(x) sqrt(mean(x^2))
  if (type == "meanabs") {
    return(abs(mean(e.post)))
  } else if (type == "rmse") {
    return(rms(e.post))
  } else if (type == "ratio") {
    ## Abadie post/pre RMSPE ratio.
    ep <- e.pre[is.finite(e.pre)]
    denom <- if (length(ep) > 0L) rms(ep) else NA_real_
    if (is.na(denom) || denom <= 0) return(NA_real_)
    return(rms(e.post) / denom)
  } else if (type == "studentized") {
    ep <- e.pre[is.finite(e.pre)]
    denom <- if (length(ep) > 1L) stats::sd(ep) else NA_real_
    if (is.na(denom) || denom <= 0) return(NA_real_)
    return(mean(abs(e.post)) / denom)
  }
  stop("conformal: unknown score type '", type, "'.")
}

## ---- 2. conformal p-value (weighted-capable) --------------------------------
## s.tr : scalar treated score; s.co : numeric vector of control scores.
## w.co : optional non-negative control weights (weighted conformal); NULL = 1.
## w.tr : treated self-weight (default 1).  Returns p in (0, 1].
.conformal_pval <- function(s.tr, s.co, w.co = NULL, w.tr = 1) {
  ok <- is.finite(s.co)
  s.co <- s.co[ok]
  if (is.null(w.co)) w.co <- rep(1, length(s.co)) else w.co <- w.co[ok]
  num <- w.tr + sum(w.co[s.co >= s.tr])
  den <- w.tr + sum(w.co)
  num / den
}

## effective sample size for weighted conformal
.conformal_neff <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0L) return(0)
  (sum(w)^2) / sum(w^2)
}

## ---- 3. interval inversion --------------------------------------------------
## Closed form for the average effect under the "meanabs" score.
## g.tr : scalar treated mean post-period gap; g.co : control mean gaps.
## Returns c(lower, upper) for att.avg.  Width = +/- the rank-corrected
## (1 - alpha) quantile of |g.co|.
.conformal_ci_meanabs <- function(g.tr, g.co, alpha = 0.05, w.co = NULL) {
  ag <- abs(g.co[is.finite(g.co)])
  n  <- length(ag)
  k  <- ceiling((1 - alpha) * (n + 1))      # conformal quantile rank
  if (k > n) {                              # alpha < 1/(n+1): unbounded
    return(c(-Inf, Inf))
  }
  q <- sort(ag)[k]
  c(g.tr - q, g.tr + q)
}

## Grid inversion for ratio / studentized / rmse (refit-free: control scores are
## tau-independent; only the treated score is recomputed on the grid).
## gap.tr.post / gap.tr.pre : treated gap paths.  s.co : control scores.
## Returns c(lower, upper) = range of accepted constant effects tau.
.conformal_ci_grid <- function(gap.tr.post, gap.tr.pre, s.co, type, alpha = 0.05,
                               w.co = NULL, w.tr = 1, n.grid = 401L, span = NULL) {
  center <- mean(gap.tr.post[is.finite(gap.tr.post)])
  if (is.null(span)) {
    sd.co <- stats::sd(s.co[is.finite(s.co)])
    span  <- max(abs(gap.tr.post), na.rm = TRUE) + 6 * (if (is.finite(sd.co)) sd.co else 1)
  }
  grid <- seq(center - span, center + span, length.out = n.grid)
  accept <- vapply(grid, function(tau) {
    s.tr <- .conformal_score(gap.tr.post - tau, gap.tr.pre, type)
    if (!is.finite(s.tr)) return(FALSE)
    .conformal_pval(s.tr, s.co, w.co = w.co, w.tr = w.tr) > alpha
  }, logical(1))
  if (!any(accept)) return(c(NA_real_, NA_real_))
  rng <- range(grid[accept])
  ## flag if the acceptance set hit a grid edge (interval may be unbounded)
  if (accept[1] || accept[n.grid]) {
    attr(rng, "edge") <- TRUE
  }
  rng
}
