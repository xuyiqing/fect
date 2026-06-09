## -----------------------------------------------------------------------------
## conformal.R
##
## Cross-sectional conformal inference for separated, symmetric counterfactual
## estimators (vartype = "conformal").  The treated unit's post-treatment
## prediction error is calibrated against the leave-one-control-out prediction
## errors of the donors and read off as a RANK, rather than resampled and
## normal-wrapped (parametric bootstrap) or quantile-resampled (nonparametric
## bootstrap).
##
## Design: statsclaw-workspace/fect/runs/REQ-conformal/spec.md
## Decision (2026-06-08): FAMILY A ONLY.  Every interval is a level statistic
##
##     S_i(tau) = | m_i - tau | / scale_i ,   interval  m_tr +/- scale_tr * Q ,
##
## where m_i is a location of unit i's post-period gaps (center: mean | median |
## per-horizon) and scale_i is a per-unit scale from its pre-period gaps
## (scale: none | sd | rmspe | mad | diff | model-se).  Q is the conformal
## (1 - alpha) quantile of the control scores |m_j| / scale_j.  This family is
## never empty (the numerator vanishes at tau = m_tr) and needs no grid inversion.
## The path / constant-effect scores (rmse / studentized-path / ratio), which
## could return an empty acceptance set, were removed because the empty case is
## confusing to report.
##
## This file is the estimator-AGNOSTIC core: centers, scales, the conformal
## quantile / p-value, and the closed-form interval.  The calibration loop
## (leave-one-control-out via impute_Y0) and the vartype dispatch live in
## boot.R / default.R and feed this core the per-unit gap paths.
## -----------------------------------------------------------------------------

## ---- 1. location (center) of a unit's post-period gap path -------------------
## "per-horizon" is not a scalar; conformal_calibrate handles it by calling this
## per post period, so here we cover the scalar centers only.
.conformal_center <- function(e.post, type = "mean") {
  e <- e.post[is.finite(e.post)]
  if (length(e) == 0L) return(NA_real_)
  if (type == "median") return(stats::median(e))
  mean(e)
}

## ---- 2. per-unit scale from a unit's pre-period gap path ---------------------
## A common (unit-invariant) scale cancels in the ranking, so only per-unit
## variation matters; that is why "none" (meanabs) and a per-unit scale differ.
## "model-se" is not computable from e.pre alone -- conformal_calibrate supplies
## it per unit -- so this returns NA to signal the caller must override.
.conformal_scale <- function(e.pre, type = "none") {
  if (type == "none") return(1)
  if (type == "model-se") return(NA_real_)        # supplied externally
  e <- e.pre[is.finite(e.pre)]
  if (type == "sd")    return(if (length(e) > 1L) stats::sd(e)          else NA_real_)
  if (type == "rmspe") return(if (length(e) > 0L) sqrt(mean(e^2))       else NA_real_)
  if (type == "mad")   return(if (length(e) > 1L) stats::mad(e)         else NA_real_)
  if (type == "diff")  return(if (length(e) > 1L) stats::sd(diff(e))    else NA_real_)
  stop("conformal: unknown scale type '", type, "'.")
}

## ---- 3. conformal quantile and p-value (weighted-capable) --------------------
## Weighted (1 - alpha) quantile of the control scores, with the treated unit held
## as a +Inf atom of weight w.tr (Tibshirani et al. 2019).  Unweighted, this is
## the ceil((1 - alpha)(n + 1))-th smallest control score.  Returns Inf when the
## level is unreachable (alpha below the resolution floor 1 / (n + 1)).
.conformal_quantile <- function(s.co, alpha, w.co = NULL, w.tr = 1) {
  ok <- is.finite(s.co); s <- s.co[ok]
  n <- length(s); if (n == 0L) return(Inf)
  w <- if (is.null(w.co)) rep(1, n) else w.co[ok]
  ord <- order(s); s <- s[ord]; w <- w[ord]
  cum <- cumsum(w) / (sum(w) + w.tr)
  idx <- which(cum >= 1 - alpha)
  if (length(idx) == 0L) return(Inf)      # mass reached only at the +Inf atom
  s[idx[1L]]
}

## p in (0, 1]: rank of the treated score among the controls.
.conformal_pval <- function(s.tr, s.co, w.co = NULL, w.tr = 1) {
  ok <- is.finite(s.co); s.co <- s.co[ok]
  if (is.null(w.co)) w.co <- rep(1, length(s.co)) else w.co <- w.co[ok]
  (w.tr + sum(w.co[s.co >= s.tr])) / (w.tr + sum(w.co))
}

## effective sample size for weighted conformal
.conformal_neff <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0L) return(0)
  (sum(w)^2) / sum(w^2)
}

## ---- 4. closed-form level interval ------------------------------------------
## m.tr, scale.tr : treated center and per-unit scale.
## m.co, scale.co : control centers and per-unit scales (vectors).
## Control score S_j = |m_j| / scale_j (no effect under H0); treated
## S_tr(tau) = |m.tr - tau| / scale.tr.  Acceptance {tau : S_tr(tau) <= Q} =
## m.tr +/- scale.tr * Q.  NEVER empty; unbounded c(-Inf, Inf) only when alpha is
## below the resolution floor (too few controls for the requested level).
.conformal_ci_level <- function(m.tr, scale.tr, m.co, scale.co, alpha = 0.05, w.co = NULL) {
  ok <- is.finite(m.co) & is.finite(scale.co) & scale.co > 0
  s.co <- abs(m.co[ok]) / scale.co[ok]
  w <- if (is.null(w.co)) NULL else w.co[ok]
  Q <- .conformal_quantile(s.co, alpha, w.co = w)
  if (!is.finite(Q) || !is.finite(scale.tr) || scale.tr <= 0) return(c(-Inf, Inf))
  half <- scale.tr * Q
  c(m.tr - half, m.tr + half)
}

## ---- 5. calibration: deterministic leave-one-control-out --------------------
## Mirror of boot.R::draw.error but HELD-OUT (not resampled): each valid control
## is predicted once from the other controls, giving its out-of-fold residual
## path. Estimator-agnostic -- built-ins via impute_Y0(method); a custom separated
## learner via `conformal.fit`. Inputs are the preprocessed matrices available
## inside fect_boot(): Y, D, X, I, II, T.on (all TT x N; X is TT x N x p or NULL),
## r.cv (selected rank), and eff (TT x N point-fit gap matrix).
##
## Returns: list(att, ci, p.value, score.tr, score.co, status, n.calib, n.eff,
##               scale, center, weight, form).
## Scope (verified): block design (common onset), method gsynth/ife, nevertreated
## (separated) calibration. Staggered uses the union post-window (per-cohort
## windows are a TODO, Phase 4); pooled (notyettreated) fits fall back to
## nevertreated with a warning. MC is not yet routed through impute_Y0.
conformal_calibrate <- function(Y, D, X = NULL, I, II, T.on, r.cv, eff,
                                method = "gsynth", predictive = "nevertreated",
                                force = 3L, hasRevs = 0L, tol = 1e-5,
                                max.iteration = 1000L, norm.para = NULL,
                                scale = "none", center = "mean", weight = "cell",
                                band.type = "pointwise", cutoff = "per-period",
                                alpha = 0.05, conformal.fit = NULL) {

  TT <- nrow(Y); N <- ncol(Y)
  sum.D <- colSums(D)
  id.tr <- which(sum.D > 0); id.co <- which(sum.D == 0)
  Ntr <- length(id.tr); Nco <- length(id.co)
  if (Nco < 2L) stop("conformal: need at least 2 controls.")

  ## model-se needs the estimator's per-unit prediction SE (not wired yet).
  if (identical(scale, "model-se")) {
    stop("conformal: scale = \"model-se\" is not yet implemented; use one of ",
         "\"none\", \"sd\", \"rmspe\", \"mad\", \"diff\".", call. = FALSE)
  }

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

  ## --- one held-out fit per valid control, reusing the parametric refit path.
  ## fake-treated column = held-out control j's DATA carrying a treated unit's
  ## TIMING (D, T.on, II must agree with the assigned onset; only Y/I are j's).
  d.pattern <- D[, id.tr[1]]                  # a treated D column (defines onset)
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

  ## --- full per-control LOO gap paths (calendar-indexed), and a per-unit scale
  ## from each control's pre-period gaps. Keeping the whole path lets us serve
  ## both the scalar average interval AND the per-period band from one pass.
  obs <- (I == 1)
  G.co  <- matrix(NA_real_, length(valid.co), TT)   # control x calendar gap
  sc.co <- rep(NA_real_, length(valid.co))
  for (mi in seq_along(valid.co)) {
    j  <- valid.co[mi]
    gj <- loo_gap(j)
    gj[!obs[, j]] <- NA_real_
    G.co[mi, ] <- gj
    sc.co[mi]  <- .conformal_scale(gj[pre.idx], scale)
  }

  ## --- treated aggregate across treated units, combined per `weight`:
  ##   cell      = weight each unit by its post-obs count (per-treated-cell ATT),
  ##   unit      = each treated unit counts equally,
  ##   precision = inverse pre-period variance (down-weights noisy units).
  ## The scalar center is the weighted mean of the per-unit post centers; the
  ## per-period band uses the matching weighted treated trajectory. They coincide
  ## for a single treated unit and (mean center) for a balanced panel.
  ## Multi-treated NOTE: the aggregate is ranked against single-control gaps; a
  ## placebo-AVERAGE calibration (ranking against averages of Ntr controls) is a
  ## refinement deferred to the simulation phase.
  trc  <- eff[, id.tr, drop = FALSE]                              # TT x Ntr
  m.i  <- apply(trc[post.idx, , drop = FALSE], 2L, function(z) .conformal_center(z, center))
  n.i  <- colSums(!is.na(trc[post.idx, , drop = FALSE]))
  ## precision weight uses each unit's pre-period VARIANCE directly (decoupled
  ## from the `scale` knob, which only normalizes the score), so it down-weights
  ## noisy units even under the default scale = "none".
  var.i <- apply(trc[pre.idx, , drop = FALSE], 2L, function(z) {
    z <- z[is.finite(z)]; if (length(z) > 1L) stats::var(z) else NA_real_
  })
  wi   <- switch(weight,
                 cell      = n.i,
                 unit      = rep(1, Ntr),
                 precision = ifelse(is.finite(var.i) & var.i > 0, 1 / var.i, 0))
  wi[!is.finite(wi)] <- 0
  if (!any(wi > 0)) wi <- rep(1, Ntr)
  Wt    <- matrix(wi, TT, Ntr, byrow = TRUE); Wt[is.na(trc)] <- 0
  denom <- rowSums(Wt)
  tr.path <- ifelse(denom > 0, rowSums(trc * Wt, na.rm = TRUE) / denom, NA_real_)  # calendar
  sc.tr   <- .conformal_scale(tr.path[pre.idx], scale)

  w.co <- NULL   # overlap density-ratio weighting deferred to a later phase

  ## --- scalar average effect: weighted mean of per-unit post centers, ranked vs
  ## the per-control post centers.
  m.co  <- apply(G.co[, post.idx, drop = FALSE], 1L,
                 function(z) .conformal_center(z, center))
  keep  <- is.finite(m.co) & is.finite(sc.co) & sc.co > 0
  okm   <- is.finite(m.i) & wi > 0
  m.tr  <- if (any(okm)) sum(wi[okm] * m.i[okm]) / sum(wi[okm]) else
           .conformal_center(tr.path[post.idx], center)
  ci    <- .conformal_ci_level(m.tr, sc.tr, m.co[keep], sc.co[keep],
                               alpha = alpha, w.co = w.co)
  Ncal  <- sum(keep)
  status <- if (any(is.infinite(ci))) "unbounded" else "ok"
  s.tr     <- if (is.finite(sc.tr) && sc.tr > 0) abs(m.tr) / sc.tr else NA_real_
  s.co.vec <- abs(m.co[keep]) / sc.co[keep]
  p.value  <- .conformal_pval(s.tr, s.co.vec, w.co = w.co)

  ## --- per-period band(s) (calendar-indexed). Standardized control gaps
  ## g.tilde[j,t] = |G.co[j,t]| / sc.co[j]; the treated deviation at period t is
  ## |tr.path[t] - tau| / sc.tr. Three constructions:
  ##   pointwise + per-period (default): Q_t = conformal quantile of g.tilde[,t].
  ##   pointwise + pooled: one Q over all post-window g.tilde values.
  ##   simultaneous: one multiplier c = quantile of the per-control MAX over the
  ##     post window (sup-t / uniform band over the post path); pre-period rows
  ##     keep the per-period placebo band.
  ## The pre-period rows always show the per-period placebo band (for pre-trends).
  ok.co <- which(is.finite(sc.co) & sc.co > 0)
  Gv   <- G.co[ok.co, , drop = FALSE]
  scv  <- sc.co[ok.co]
  Gtil <- abs(Gv) / scv                      # ncal x TT, row j divided by scv[j]

  Qpool <- NULL
  if (identical(cutoff, "pooled")) {
    pv <- as.vector(Gtil[, post.idx, drop = FALSE]); pv <- pv[is.finite(pv)]
    Qpool <- .conformal_quantile(pv, alpha)
  }
  Mj    <- apply(Gtil[, post.idx, drop = FALSE], 1L,
                 function(z) { z <- z[is.finite(z)]; if (length(z)) max(z) else NA_real_ })
  c.sim <- .conformal_quantile(Mj[is.finite(Mj)], alpha)

  mk_band <- function(mode) {
    b <- matrix(NA_real_, TT, 4L,
                dimnames = list(NULL, c("eff", "CI.lower", "CI.upper", "p.value")))
    for (t in seq_len(TT)) {
      gt <- Gtil[, t]; okt <- is.finite(gt)
      if (!is.finite(tr.path[t]) || !is.finite(sc.tr) || sc.tr <= 0 || sum(okt) < 2L) next
      Qt <- if (mode == "simultaneous" && t %in% post.idx) c.sim
            else if (identical(cutoff, "pooled")) Qpool
            else .conformal_quantile(gt[okt], alpha)
      half <- if (is.finite(Qt)) sc.tr * Qt else Inf
      b[t, ] <- c(tr.path[t], tr.path[t] - half, tr.path[t] + half,
                  .conformal_pval(abs(tr.path[t]) / sc.tr, gt[okt]))
    }
    b
  }
  band     <- mk_band(band.type)                 # populates est.att
  band.sim <- if (identical(band.type, "simultaneous")) band else mk_band("simultaneous")

  list(att = m.tr, ci = c(ci[1], ci[2]), p.value = p.value,
       score.tr = s.tr, score.co = s.co.vec, status = status,
       n.calib = Ncal, n.eff = if (is.null(w.co)) Ncal else .conformal_neff(w.co),
       scale = scale, center = center, weight = weight,
       band.type = band.type, cutoff = cutoff, form = "jackknife+",
       band = band, band.sim = band.sim, post.idx = post.idx, pre.idx = pre.idx)
}
