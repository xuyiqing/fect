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

## Outcome-scale denominator a given score divides by, used to size the grid in
## the SAME units as tau (the grid is over outcome-scale constant effects, so the
## span must be outcome-scaled, not score-scaled).
.conformal_denom <- function(e.pre, type) {
  ep <- e.pre[is.finite(e.pre)]
  if (type == "ratio")       return(if (length(ep) > 0L) sqrt(mean(ep^2)) else NA_real_)
  if (type == "studentized") return(if (length(ep) > 1L) stats::sd(ep)      else NA_real_)
  1  # rmse / meanabs: numerator already in outcome units
}

## Grid inversion for ratio / studentized / rmse (refit-free: control scores are
## tau-independent; only the treated score is recomputed on the grid).
## gap.tr.post / gap.tr.pre : treated gap paths.  s.co : control scores.
## Returns c(lower, upper) = range of accepted constant effects tau.  Two special
## returns, both carried as attributes so the caller can react:
##   attr "empty" = "rejected_all_tau" : acceptance set is genuinely empty (no
##     constant effect is consistent at level alpha; common with per-unit
##     studentization when the treated pre-period sd is small) -> c(NA, NA).
##   attr "edge"  = TRUE               : acceptance reached the (expanded) grid
##     edge; the interval is effectively unbounded on that side.
.conformal_ci_grid <- function(gap.tr.post, gap.tr.pre, s.co, type, alpha = 0.05,
                               w.co = NULL, w.tr = 1, n.grid = 401L, span = NULL) {
  gp <- gap.tr.post[is.finite(gap.tr.post)]
  center <- mean(gp)
  ## span in OUTCOME units: beyond |tau - center| ~ (max control score) * denom
  ## the treated score must exceed every control score, so acceptance is
  ## impossible -- this bounds where the grid needs to look.
  denom <- .conformal_denom(gap.tr.pre, type)
  if (!is.finite(denom) || denom <= 0) denom <- 1
  s.hi <- suppressWarnings(max(s.co[is.finite(s.co)]))
  if (!is.finite(s.hi)) s.hi <- 1
  if (is.null(span)) span <- max(abs(gp)) + (s.hi + 1) * denom

  accept_on <- function(span) {
    grid <- seq(center - span, center + span, length.out = n.grid)
    keep <- vapply(grid, function(tau) {
      s.tr <- .conformal_score(gap.tr.post - tau, gap.tr.pre, type)
      if (!is.finite(s.tr)) return(FALSE)
      .conformal_pval(s.tr, s.co, w.co = w.co, w.tr = w.tr) > alpha
    }, logical(1))
    list(grid = grid, keep = keep)
  }
  a <- accept_on(span)
  ## adaptive expansion: if acceptance touches an edge the true interval may run
  ## further (guards the span heuristic against unusual gap/score scales).
  tries <- 0L
  while (any(a$keep) && (a$keep[1L] || a$keep[length(a$keep)]) && tries < 3L) {
    span <- span * 2; a <- accept_on(span); tries <- tries + 1L
  }
  if (!any(a$keep)) {
    out <- c(NA_real_, NA_real_)
    attr(out, "empty") <- "rejected_all_tau"
    return(out)
  }
  ## an accepted grid endpoint means the interval runs past it -> unbounded side.
  lo <- min(a$grid[a$keep]); hi <- max(a$grid[a$keep])
  if (a$keep[1L])               lo <- -Inf
  if (a$keep[length(a$keep)])   hi <-  Inf
  rng <- c(lo, hi)
  if (is.infinite(lo) || is.infinite(hi)) attr(rng, "edge") <- TRUE
  rng
}

## ---- 4. calibration: deterministic leave-one-control-out --------------------
## Mirror of boot.R::draw.error but HELD-OUT (not resampled): each valid control
## is predicted once from the other controls, giving its out-of-fold residual
## path. Estimator-agnostic — built-ins via impute_Y0(method); a custom separated
## learner via `conformal.fit`. Inputs are the preprocessed matrices available
## inside fect_boot(): Y, D, X, I, II, T.on (all TT x N; X is TT x N x p or NULL),
## r.cv (selected rank), and eff (TT x N point-fit gap matrix).
##
## Returns: list(att, ci, p.value, score.tr, score.co, n.calib, n.eff, score, form).
## Scope (verified): block design (common onset), method gsynth/ife, nevertreated
## (separated) calibration. Staggered uses the union post-window (per-cohort
## windows are a TODO); pooled (notyettreated) fits fall back to nevertreated with
## a warning. MC is not yet routed through impute_Y0 (upstream stop()).
conformal_calibrate <- function(Y, D, X = NULL, I, II, T.on, r.cv, eff,
                                method = "gsynth", predictive = "nevertreated",
                                force = 3L, hasRevs = 0L, tol = 1e-5,
                                max.iteration = 1000L, norm.para = NULL,
                                score = "studentized", alpha = 0.05,
                                conformal.full = FALSE, conformal.weight = "none",
                                conformal.fit = NULL) {

  TT <- nrow(Y); N <- ncol(Y)
  sum.D <- colSums(D)
  id.tr <- which(sum.D > 0); id.co <- which(sum.D == 0)
  Ntr <- length(id.tr); Nco <- length(id.co)
  if (Nco < 2L) stop("conformal: need at least 2 controls.")

  ## --- separation guard: conformal calibration is controls-only (nevertreated)
  if (!identical(predictive, "nevertreated")) {
    warning("vartype = \"conformal\" requires a separated fit; ",
            "calibrating on controls only (predictive = \"nevertreated\").")
    predictive <- "nevertreated"
  }

  ## --- valid controls (enough pre/post obs); reuse fect's screen
  valid.co <- valid_controls(list(D = D, I = I, r.cv = r.cv), method, predictive, force)
  valid.co <- intersect(valid.co, id.co)
  if (length(valid.co) < 2L) {
    stop("conformal: fewer than 2 valid controls after screening.")
  }

  ## --- post / pre period windows (block: common onset; staggered: union)
  post.idx <- which(rowSums(D[, id.tr, drop = FALSE]) > 0)   # any treated on
  pre.idx  <- setdiff(seq_len(TT), post.idx)
  if (length(post.idx) == 0L || length(pre.idx) == 0L) {
    stop("conformal: could not identify pre/post windows from D.")
  }
  if (Ntr > 1L) {
    onsets <- vapply(id.tr, function(i) min(which(D[, i] == 1)), integer(1))
    if (length(unique(onsets)) > 1L) {
      warning("conformal: staggered onsets detected; using the union post-window. ",
              "Per-cohort windows are not yet implemented.")
    }
  }

  ## --- one held-out fit per valid control, reusing the parametric refit path
  d.pattern <- D[, id.tr[1]]                  # a treated D column (defines onset)
  sub3 <- function(A, idx) if (is.null(A)) NULL else A[, idx, , drop = FALSE]

  ## fake-treated column = held-out control j's DATA carrying a treated unit's
  ## TIMING (D, T.on, II must agree with the assigned onset; only Y/I are j's).
  loo_gap <- function(j) {
    co.rest <- setdiff(valid.co, j)
    if (!is.null(conformal.fit)) {
      y0 <- conformal.fit(Y = Y, X = X, time = seq_len(TT),
                          control.ids = co.rest, target.id = j, T0 = max(pre.idx))
      return(Y[, j] - y0)
    }
    Y.ps   <- cbind(Y[, j],           Y[, co.rest, drop = FALSE])
    D.ps   <- cbind(d.pattern,        D[, co.rest, drop = FALSE])
    Ton.ps <- cbind(T.on[, id.tr[1]], T.on[, co.rest, drop = FALSE])
    I.ps   <- cbind(I[, j],           I[, co.rest, drop = FALSE])
    II.ps  <- cbind(I[, j] * (d.pattern == 0), II[, co.rest, drop = FALSE])
    X.ps   <- if (is.null(X)) NULL else
              array(c(X[, j, , drop = FALSE], X[, co.rest, , drop = FALSE]),
                    dim = c(TT, 1L + length(co.rest), dim(X)[3]))
    synth <- try(impute_Y0(
      method = method, predictive = predictive,
      Y = Y.ps, X = X.ps, D = D.ps, W = NULL, I = I.ps, II = II.ps,
      T.on = Ton.ps, tuning = r.cv, boot = 1,
      force = force, hasRevs = hasRevs, tol = tol,
      max.iteration = max.iteration, norm.para = norm.para
    ), silent = TRUE)
    if (inherits(synth, "try-error") || !("eff" %in% names(synth))) return(rep(NA_real_, TT))
    g <- as.matrix(synth$eff.tr)[, 1]
    if (!is.null(norm.para)) g <- g / norm.para[1]
    g
  }

  ## --- control scores (tau-independent) and mean post gaps
  obs <- (I == 1)
  S.co <- rep(NA_real_, length(valid.co))
  g.co <- rep(NA_real_, length(valid.co))
  for (m in seq_along(valid.co)) {
    j  <- valid.co[m]
    gj <- loo_gap(j)
    ep <- gj[post.idx][obs[post.idx, j]]
    eq <- gj[pre.idx ][obs[pre.idx,  j]]
    S.co[m] <- .conformal_score(ep, eq, score)
    g.co[m] <- mean(ep, na.rm = TRUE)
  }
  keep <- is.finite(S.co)
  S.co <- S.co[keep]; g.co <- g.co[keep]
  Ncal <- length(S.co)

  ## --- treated statistic (averaged over treated units, post window)
  eff.tr.mat <- eff[, id.tr, drop = FALSE]
  tr.post.path <- rowMeans(eff.tr.mat[post.idx, , drop = FALSE], na.rm = TRUE)
  tr.pre.path  <- rowMeans(eff.tr.mat[pre.idx,  , drop = FALSE], na.rm = TRUE)
  g.tr <- mean(tr.post.path, na.rm = TRUE)

  ## --- weights (overlap) — placeholder until loading_bound wiring (Phase 3)
  w.co <- NULL
  if (identical(conformal.weight, "overlap")) {
    warning("conformal.weight = \"overlap\" not yet wired; using unweighted.")
  }

  ## --- interval
  if (score == "meanabs") {
    ci <- .conformal_ci_meanabs(g.tr, g.co, alpha = alpha, w.co = w.co)
  } else {
    ci <- .conformal_ci_grid(tr.post.path, tr.pre.path, S.co, score,
                             alpha = alpha, w.co = w.co)
  }
  ## status: "ok" | "unbounded" (alpha too small for Ncal) | "empty" (a
  ## constant-effect score that rejects every tau -- correct conformal output,
  ## not a failure; warn so callers/users know to read it as "no constant effect
  ## consistent at level alpha", or switch to score = "meanabs").
  status <- "ok"
  if (identical(attr(ci, "empty"), "rejected_all_tau")) {
    status <- "empty"
    warning("conformal: no constant effect is consistent at level ", alpha,
            " under score = \"", score, "\" (the treated path is rejected at ",
            "every tau). Reported interval is empty; score = \"meanabs\" gives ",
            "a level interval for the average effect that cannot be empty.")
  } else if (any(is.infinite(ci))) {
    status <- "unbounded"
  }
  S.tr <- .conformal_score(tr.post.path, tr.pre.path, score)
  p.value <- .conformal_pval(S.tr, S.co, w.co = w.co)

  list(att = g.tr, ci = c(ci[1], ci[2]), p.value = p.value,
       score.tr = S.tr, score.co = S.co, status = status,
       n.calib = Ncal, n.eff = if (is.null(w.co)) Ncal else .conformal_neff(w.co),
       score = score, form = if (conformal.full) "full" else "jackknife+")
}
