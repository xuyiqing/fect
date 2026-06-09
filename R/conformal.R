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

  ## --- treatment-incidence structure (handles block AND staggered uniformly).
  ## post.idx is the UNION post window (any treated unit on); nt[t] is the number
  ## of treated units in post at calendar t; onsets are per-unit first-treated
  ## periods. For a control's leave-one-out fit, the union window is masked so its
  ## gaps are held-out everywhere a treated unit is post (block: a single onset).
  Dtr        <- D[, id.tr, drop = FALSE]                    # TT x Ntr, 1 = post
  nt         <- rowSums(Dtr)                                 # treated post count / period
  union.post <- nt > 0
  post.idx   <- which(union.post)
  pre.idx    <- setdiff(seq_len(TT), post.idx)
  if (length(post.idx) == 0L || length(pre.idx) == 0L) {
    stop("conformal: could not identify pre/post windows from D.")
  }
  onsets    <- apply(Dtr, 2L, function(d) { w <- which(d == 1); if (length(w)) min(w) else NA_integer_ })
  staggered <- length(unique(onsets[is.finite(onsets)])) > 1L

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
    ## mask the UNION post window so the control's gaps are held-out wherever any
    ## treated unit is post (block: a single onset; staggered: the union).
    II.ps  <- cbind(I[, j] * as.numeric(!union.post), II[, co.rest, drop = FALSE])
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

  ## --- per-treated-unit summaries, each over the unit's OWN post / pre window
  ## (block: every unit shares the union window; staggered: per-unit windows).
  ## `weight`: cell = per-treated-cell ATT (n_i); unit = equal per unit;
  ## precision = inverse pre-period variance (decoupled from the score `scale`).
  trc        <- eff[, id.tr, drop = FALSE]                       # TT x Ntr
  obstr      <- (I[, id.tr, drop = FALSE] == 1)
  postmask.i <- (Dtr == 1)                                       # own post cells
  premask.i  <- (Dtr == 0) & obstr                              # own observed pre cells
  m.i  <- vapply(seq_len(Ntr), function(i)
            .conformal_center(trc[postmask.i[, i], i], center), numeric(1))
  n.i  <- colSums(postmask.i)
  var.i <- vapply(seq_len(Ntr), function(i) {
            z <- trc[premask.i[, i], i]; z <- z[is.finite(z)]
            if (length(z) > 1L) stats::var(z) else NA_real_ }, numeric(1))
  wi   <- switch(weight, cell = n.i, unit = rep(1, Ntr),
                 precision = ifelse(is.finite(var.i) & var.i > 0, 1 / var.i, 0))
  wi[!is.finite(wi)] <- 0
  if (!any(wi > 0)) wi <- rep(1, Ntr)

  ## treated trajectory (calendar): weighted mean over the units that are POST at
  ## t (effect window) or, where none are, over the units that are PRE at t
  ## (placebo window). For block this is the per-period weighted mean throughout.
  wsum  <- function(e, w) { ok <- is.finite(e) & is.finite(w) & w > 0
                            if (any(ok)) sum(e[ok] * w[ok]) / sum(w[ok]) else NA_real_ }
  Wt    <- matrix(wi, TT, Ntr, byrow = TRUE)
  postW <- Wt * postmask.i; preW <- Wt * premask.i
  tr.path <- vapply(seq_len(TT), function(t)
               wsum(trc[t, ], if (nt[t] > 0) postW[t, ] else preW[t, ]), numeric(1))
  sc.tr   <- .conformal_scale(tr.path[pre.idx], scale)

  w.co <- NULL   # overlap density-ratio weighting deferred to a later phase

  ## --- scalar average effect: weighted mean of per-unit post centers, ranked
  ## against per-control placebo centers (each control's gap aggregated over the
  ## union post window the SAME way the treated ATT is, i.e. nt-weighted = the
  ## per-treated-cell ATT; for block nt is constant so this is a plain mean).
  ntp  <- nt[post.idx]
  m.co <- apply(G.co, 1L, function(g) wsum(g[post.idx], ntp))
  keep <- is.finite(m.co) & is.finite(sc.co) & sc.co > 0
  okm  <- is.finite(m.i) & wi > 0
  m.tr <- if (any(okm)) sum(wi[okm] * m.i[okm]) / sum(wi[okm]) else
          .conformal_center(tr.path[post.idx], center)
  ci   <- .conformal_ci_level(m.tr, sc.tr, m.co[keep], sc.co[keep], alpha = alpha, w.co = w.co)
  Ncal <- sum(keep)
  status <- if (any(is.infinite(ci))) "unbounded" else "ok"
  s.tr     <- if (is.finite(sc.tr) && sc.tr > 0) abs(m.tr) / sc.tr else NA_real_
  s.co.vec <- abs(m.co[keep]) / sc.co[keep]
  p.value  <- .conformal_pval(s.tr, s.co.vec, w.co = w.co)

  ## --- standardized control gaps g.tilde[j,t] = |G.co[j,t]| / sc.co[j], the
  ## per-control MAX over the post window (for the sup-t multiplier), and the
  ## pooled post-window values (for cutoff = "pooled"). The band builders take a
  ## level `a` so the outer (alpha) and inner (2*alpha) bands reuse one pass.
  ok.co  <- which(is.finite(sc.co) & sc.co > 0)
  Gv     <- G.co[ok.co, , drop = FALSE]; scv <- sc.co[ok.co]
  Gtil   <- abs(Gv) / scv                      # ncal x TT
  Mj     <- apply(Gtil[, post.idx, drop = FALSE], 1L,
                  function(z) { z <- z[is.finite(z)]; if (length(z)) max(z) else NA_real_ })
  Mj     <- Mj[is.finite(Mj)]
  pvpool <- if (identical(cutoff, "pooled")) {
              v <- as.vector(Gtil[, post.idx, drop = FALSE]); v[is.finite(v)]
            } else NULL
  qt_post <- function(gp, a) {                 # post-window cutoff at level a
    if (!is.null(pvpool)) .conformal_quantile(pvpool, a) else .conformal_quantile(gp, a)
  }

  ## --- CALENDAR band (per calendar period) -> est.eff.calendar.
  mk_cal <- function(mode, a) {
    c.sim <- .conformal_quantile(Mj, a)
    b <- matrix(NA_real_, TT, 4L, dimnames = list(NULL, c("eff", "CI.lower", "CI.upper", "p.value")))
    for (t in seq_len(TT)) {
      gt <- Gtil[, t]; okt <- is.finite(gt)
      if (!is.finite(tr.path[t]) || !is.finite(sc.tr) || sc.tr <= 0 || sum(okt) < 2L) next
      Qt <- if (mode == "simultaneous" && union.post[t]) c.sim
            else if (union.post[t]) qt_post(gt[okt], a)
            else .conformal_quantile(gt[okt], a)
      half <- if (is.finite(Qt)) sc.tr * Qt else Inf
      b[t, ] <- c(tr.path[t], tr.path[t] - half, tr.path[t] + half,
                  .conformal_pval(abs(tr.path[t]) / sc.tr, gt[okt]))
    }
    b
  }
  band     <- mk_cal(band.type, alpha)
  band.sim <- if (identical(band.type, "simultaneous")) band else mk_cal("simultaneous", alpha)

  ## --- EVENT-TIME band (per relative period) -> est.att. For block this is the
  ## calendar band re-indexed; for staggered it aggregates cohorts at each relative
  ## time, pooling the control gaps at the matching calendar cells (cross-sectional
  ## exchangeability across units, with stationarity across cohort onsets).
  rel.mat <- T.on[, id.tr, drop = FALSE]
  etimes  <- sort(unique(rel.mat[obstr]))
  mk_et <- function(mode, a) {
    c.sim <- .conformal_quantile(Mj, a)
    b <- matrix(NA_real_, length(etimes), 4L,
                dimnames = list(as.character(etimes), c("eff", "CI.lower", "CI.upper", "p.value")))
    for (k in seq_along(etimes)) {
      cells <- which(rel.mat == etimes[k] & obstr, arr.ind = TRUE)   # (t, unit) pairs
      if (!nrow(cells)) next
      tg  <- wsum(trc[cells], wi[cells[, 2L]])
      ts  <- unique(cells[, 1L])                                     # calendar periods at this rel time
      ispost <- any(union.post[ts])
      gp  <- as.vector(Gtil[, ts, drop = FALSE]); gp <- gp[is.finite(gp)]
      if (!is.finite(tg) || !is.finite(sc.tr) || sc.tr <= 0 || length(gp) < 2L) next
      Qk <- if (mode == "simultaneous" && ispost) c.sim
            else if (ispost) qt_post(gp, a)
            else .conformal_quantile(gp, a)
      half <- if (is.finite(Qk)) sc.tr * Qk else Inf
      b[k, ] <- c(tg, tg - half, tg + half, .conformal_pval(abs(tg) / sc.tr, gp))
    }
    b
  }
  a.inner       <- min(2 * alpha, 0.5)                 # inner (1 - 2*alpha) band
  band.et       <- mk_et(band.type, alpha)
  band.et.inner <- mk_et(band.type, a.inner)
  band.et.sim   <- if (identical(band.type, "simultaneous")) band.et else mk_et("simultaneous", alpha)

  list(att = m.tr, ci = c(ci[1], ci[2]), p.value = p.value,
       score.tr = s.tr, score.co = s.co.vec, status = status,
       n.calib = Ncal, n.eff = if (is.null(w.co)) Ncal else .conformal_neff(w.co),
       scale = scale, center = center, weight = weight,
       band.type = band.type, cutoff = cutoff, staggered = staggered, form = "jackknife+",
       band = band, band.sim = band.sim, band.et = band.et,
       band.et.inner = band.et.inner, band.et.sim = band.et.sim,
       post.idx = post.idx, pre.idx = pre.idx)
}
