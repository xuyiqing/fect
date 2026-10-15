###############################################################################
## vartype = "conformal": pooled information-matched Jackknife+ (v2)
###############################################################################
##
## Implements the pooled information-matched Jackknife+ interval of Ziyi Liu,
## "Unit-Level Conformal Inference for Panel Counterfactuals" (write-up v2):
## Section 3.4 (pooled query and background information), Section 4.3
## (Jackknife+ for a pooled estimator), the calibration q_gamma of eq.
## (jplus-calibration-main), Section 4.5 (date-specific coverage; q_gamma = 0
## gives the whole real line) and the appendix algorithm "Pooled
## information-matched Jackknife+ interval". The p-value is the interval-dual
## p-value of eq. (placebo-jplus-pvalue-main) evaluated at c = 0.
##
## For one treated unit (unit 0) and one post-treatment date t, the never-
## treated controls C form the calibration pool and N = |C| + 1. In fold j:
##   * the background B_j = C \ {j} keeps its complete panels;
##   * unit 0 and control j are queries: each exposes only its pre-treatment
##     history (control j: the treated unit's pre-treatment window), and
##     their post-treatment outcomes are hidden;
##   * the two queries are predicted SEPARATELY against B_j. gsynth: one
##     factor estimation on B_j, then each query's loading projection. ife
##     (EM, not-yet-treated estimator) and mc: one component fit per query
##     (query + B_j);
##   * tuning is fixed at the original fit's values (r, lambda, force,
##     covariates, tol, max.iteration); there is no cross-validation inside
##     folds.
## With R_j0(t) = |Y_jt - M_j0(t)|, a_j = M_0j(t) - R_j0(t) and
## b_j = M_0j(t) + R_j0(t), the interval for Y_0t(0) is [a_(q), b_(N - q)] and
## the effect interval is [Y_0t - b_(N - q), Y_0t - a_(q)]. The point estimate
## is the original fit's (all controls).
##
## EXTENSION BEYOND v2 (Yiqing, 2026-10-14; needs Ziyi's sign-off): intervals
## for aggregated effects (est.att by event time, est.avg, est.avg.unit and the
## calendar-time ATT). In each fold the aggregated treated query is the
## weighted average of the treated units' fold predictions (each treated unit
## predicted with its own pre-period against the same B_j). The control's
## error uses the same weights over the same cells: for a cell (i, t) the
## control is predicted over treated unit i's window, so with staggered
## adoption the control's error at event time l is the cohort-size weighted
## average of its signed errors over the cohorts' calendar windows. The
## absolute value of that weighted signed error is R_j, and the endpoint
## construction is the same as above. Comparing an average of several treated
## units with a single control's error is conservative when many units are
## treated (the control's error does not average down).
###############################################################################


## Post-fit functions that read bootstrap draws or a variance matrix have no
## conformal version; stop with a message rather than a downstream error.
.fect_stop_if_conformal <- function(fit, what) {
  if (is.list(fit) && identical(fit[["vartype"]], "conformal")) {
    stop(what, " is not available for a conformal fit (vartype = ",
         "\"conformal\"): it needs bootstrap draws or a variance matrix, ",
         "which the Jackknife+ does not produce. Refit with vartype = ",
         "\"bootstrap\" or \"parametric\".", call. = FALSE)
  }
  invisible(NULL)
}


## q_gamma of v2 eq. (jplus-calibration-main):
## max{q in 0..floor(N/2) : (2q - 1)_+ / N <= gamma}. q = 0 means the counting
## bound certifies no finite interval.
.conformal_q <- function(N, gamma) {
  if (!is.finite(N) || N < 1) return(0L)
  qs <- 0:floor(N / 2)
  ok <- pmax(2 * qs - 1, 0) <= gamma * N + 1e-10
  as.integer(max(qs[ok]))
}

## Jackknife+ interval and interval-dual p-value from the fold quantities.
## M0: the treated query's fold predictions (length J); R: the controls'
## absolute fold errors (length J); y: the observed treated outcome; gamma:
## miscoverage. Returns the prediction interval for Y(0), the effect interval,
## the p-value at c = 0, N and q.
.conformal_interval <- function(M0, R, y, gamma) {
  ok <- is.finite(M0) & is.finite(R)
  M0 <- M0[ok]
  R <- R[ok]
  J <- length(M0)
  N <- J + 1L
  q <- if (J >= 1L) .conformal_q(N, gamma) else 0L
  out <- c(pred.lower = -Inf, pred.upper = Inf,
           CI.lower = -Inf, CI.upper = Inf,
           p.value = 1, N = N, q = q)
  if (J < 1L || !is.finite(y)) {
    out[c("p.value")] <- NA_real_
    if (!is.finite(y)) out[c("CI.lower", "CI.upper")] <- NA_real_
    return(out)
  }
  a <- M0 - R
  b <- M0 + R
  ## interval-dual p-value at c = 0 (v2 eq. placebo-jplus-pvalue-main):
  ## strict comparisons, the larger of the two tail counts
  m <- max(sum(y < a), sum(y > b))
  out["p.value"] <- min(1, (2 * (N - m) - 1) / N)
  if (q >= 1L) {
    a.s <- sort(a)
    b.s <- sort(b)
    lo <- a.s[q]
    hi <- b.s[N - q]
    out[c("pred.lower", "pred.upper")] <- c(lo, hi)
    out[c("CI.lower", "CI.upper")] <- c(y - hi, y - lo)
  }
  out
}


## Build the column of a pseudo-treated copy of control j for query window k:
## the control exposes only the treated window's observed pre-treatment dates
## and hides every date from the window's first post-treatment date on.
.conformal_copy <- function(j, win, Y, X, I, TT) {
  tt <- seq_len(TT)
  D.c <- as.numeric(tt >= win$start)
  I.c <- I[, j]
  I.c[tt < win$start] <- as.numeric(tt[tt < win$start] %in% win$pre)
  II.c <- I.c * (1 - D.c)
  Y.c <- Y[, j]
  Y.c[I.c == 0] <- 0
  T.on.c <- tt + win$offset
  T.on.c[I.c == 0] <- NA
  X.c <- if (!is.null(X) && dim(X)[3] > 0) X[, j, , drop = FALSE] else NULL
  list(Y = Y.c, D = D.c, I = I.c, II = II.c, T.on = T.on.c, X = X.c)
}

## Bind sub-panel columns: background controls (complete) plus query columns.
.conformal_panel <- function(cols, copies, Y, X, D, I, II, T.on) {
  Ys <- Y[, cols, drop = FALSE]
  Ds <- D[, cols, drop = FALSE]
  Is <- I[, cols, drop = FALSE]
  IIs <- II[, cols, drop = FALSE]
  Ts <- T.on[, cols, drop = FALSE]
  p <- if (is.null(X)) 0 else dim(X)[3]
  Xs <- if (p > 0) X[, cols, , drop = FALSE] else NULL
  for (cp in copies) {
    Ys <- cbind(Ys, cp$Y)
    Ds <- cbind(Ds, cp$D)
    Is <- cbind(Is, cp$I)
    IIs <- cbind(IIs, cp$II)
    Ts <- cbind(Ts, cp$T.on)
    if (p > 0) Xs <- abind::abind(Xs, cp$X, along = 2)
  }
  if (p == 0) Xs <- array(0, dim = c(nrow(Ys), ncol(Ys), 0))
  list(Y = Ys, D = Ds, I = Is, II = IIs, T.on = Ts, X = Xs)
}

## One component fit with the original fit's tuning; returns Y.ct (TT x n)
## or NULL when the fit fails.
.conformal_component_fit <- function(pan, method, r, lambda, validF, force,
                                     tol, max.iteration) {
  fit <- try(suppressMessages(suppressWarnings(
    if (method == "gsynth") {
      fect_nevertreated(
        Y = pan$Y, X = pan$X, D = pan$D, W = NULL, I = pan$I, II = pan$II,
        T.on = pan$T.on, CV = 0, r = r, binary = FALSE, force = force,
        hasRevs = 0, tol = tol, max.iteration = max.iteration, boot = 1,
        norm.para = NULL, parallel = FALSE
      )
    } else if (method == "ife") {
      fect_fe(
        Y = pan$Y, X = pan$X, D = pan$D, W = NULL, I = pan$I, II = pan$II,
        T.on = pan$T.on, r.cv = r, binary = FALSE, force = force,
        hasRevs = 0, tol = tol, max.iteration = max.iteration, boot = 1,
        norm.para = NULL
      )
    } else {
      fect_mc(
        Y = pan$Y, X = pan$X, D = pan$D, W = NULL, I = pan$I, II = pan$II,
        T.on = pan$T.on, lambda.cv = lambda, force = force, hasF = validF,
        hasRevs = 0, tol = tol, max.iteration = max.iteration, boot = 1,
        norm.para = NULL
      )
    }
  )), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit$Y.ct)) return(NULL)
  as.matrix(fit$Y.ct)
}

## Fold j: background B_j = co \ {j}. Returns the treated queries' predictions
## (TT x Ntr) and control j's predictions over each query window (TT x K).
.conformal_fold <- function(jj, ctx) {
  co <- ctx$co
  tr <- ctx$tr
  TT <- ctx$TT
  j <- co[jj]
  Bj <- co[-jj]
  nB <- length(Bj)
  K <- length(ctx$wins)
  ks <- which(ctx$elig[jj, ])
  M0 <- matrix(NA_real_, TT, length(tr))
  Mj <- matrix(NA_real_, TT, K)
  copies <- lapply(ks, function(k) {
    .conformal_copy(j, ctx$wins[[k]], ctx$Y, ctx$X, ctx$I, TT)
  })
  fit_cols <- function(cols, cps) {
    pan <- .conformal_panel(cols, cps, ctx$Y, ctx$X, ctx$D, ctx$I, ctx$II,
                            ctx$T.on)
    .conformal_component_fit(pan, ctx$method, ctx$r, ctx$lambda, ctx$validF,
                             ctx$force, ctx$tol, ctx$max.iteration)
  }
  n.fail <- 0L
  if (ctx$method == "gsynth") {
    ## the factors come from B_j alone and every query is projected on its
    ## own pre-treatment history, so one call gives all the separate
    ## component predictions (no query supplies outcomes to another's fit)
    yct <- fit_cols(c(Bj, tr), copies)
    if (is.null(yct)) {
      n.fail <- 1L
    } else {
      M0[, ] <- yct[, nB + seq_along(tr)]
      if (length(ks) > 0) {
        Mj[, ks] <- yct[, nB + length(tr) + seq_along(ks)]
      }
    }
  } else {
    ## EM (ife) and matrix completion: a query's cells enter the low-rank
    ## fit, so each query gets its own component fit with B_j
    for (ii in seq_along(tr)) {
      yct <- fit_cols(c(Bj, tr[ii]), list())
      if (is.null(yct)) {
        n.fail <- n.fail + 1L
      } else {
        M0[, ii] <- yct[, nB + 1L]
      }
    }
    for (kk in seq_along(ks)) {
      yct <- fit_cols(Bj, copies[kk])
      if (is.null(yct)) {
        n.fail <- n.fail + 1L
      } else {
        Mj[, ks[kk]] <- yct[, nB + 1L]
      }
    }
  }
  list(M0 = M0, Mj = Mj, n.fail = n.fail)
}


## Main entry, called from fect_boot() after the point fit `out`.
conformal_jplus <- function(Y, X, D, I, II, T.on, out, method, force, tol,
                            max.iteration, alpha, norm.para = NULL,
                            parallel = FALSE, keep.sims = FALSE) {
  TT <- nrow(Y)
  N.all <- ncol(Y)
  D[is.na(D)] <- 0
  scale <- if (is.null(norm.para)) 1 else norm.para[1]
  if (!is.null(X) && length(dim(X)) == 3 && dim(X)[3] == 0) X <- NULL

  tr <- which(colSums(D * I) > 0)
  co <- which(colSums(D) == 0)
  Ntr <- length(tr)
  Nco <- length(co)

  r <- if (!is.null(out$r.cv)) as.numeric(out$r.cv[1]) else 0
  ## lambda.cv is on the outcome's scale; the component fits run on the
  ## (possibly normalized) Y without norm.para, so divide as fect_mc() does
  lambda <- if (is.null(out$lambda.cv)) NULL else out$lambda.cv / scale
  validF <- if (!is.null(out$validF)) out$validF else 1

  ## query windows: one per distinct (observed pre-treatment dates, event
  ## time offset) among the treated units
  win.key <- character(Ntr)
  wins <- list()
  for (ii in seq_len(Ntr)) {
    i <- tr[ii]
    pre <- which(D[, i] == 0 & I[, i] == 1)
    obs.on <- which(!is.na(T.on[, i]))
    offset <- T.on[obs.on[1], i] - obs.on[1]
    key <- paste0(paste(pre, collapse = ","), "|", offset)
    win.key[ii] <- key
    if (is.null(wins[[key]])) {
      wins[[key]] <- list(pre = pre, start = max(pre) + 1L, offset = offset)
    }
  }
  key.idx <- match(win.key, names(wins))
  K <- length(wins)

  ## eligibility: control j can take window k's query role only if it is
  ## observed at every one of the window's pre-treatment dates
  elig <- matrix(FALSE, Nco, K)
  for (k in seq_len(K)) {
    elig[, k] <- colSums(I[wins[[k]]$pre, co, drop = FALSE] == 1) ==
      length(wins[[k]]$pre)
  }

  ctx <- list(Y = Y, X = X, D = D, I = I, II = II, T.on = T.on, co = co,
              tr = tr, TT = TT, wins = wins, elig = elig, method = method,
              r = r, lambda = lambda, validF = validF, force = force,
              tol = tol, max.iteration = max.iteration)

  run_fold <- function(jj) .conformal_fold(jj, ctx)
  use.future <- isTRUE(parallel) && Nco > 1 &&
    !inherits(future::plan(), "sequential")
  folds <- if (use.future) {
    future.apply::future_lapply(seq_len(Nco), run_fold, future.seed = TRUE)
  } else {
    lapply(seq_len(Nco), run_fold)
  }
  n.fail <- sum(vapply(folds, function(f) f$n.fail, integer(1)))

  ## treated post-treatment cells and their fold quantities
  cell <- which(D[, tr, drop = FALSE] == 1 & I[, tr, drop = FALSE] == 1,
                arr.ind = TRUE)
  colnames(cell) <- c("t", "i")
  cell <- cell[order(cell[, "i"], cell[, "t"]), , drop = FALSE]
  n.cell <- nrow(cell)
  cM <- matrix(NA_real_, n.cell, Nco)  # treated fold predictions M_0j
  cE <- matrix(NA_real_, n.cell, Nco)  # control signed errors Y_jt - M_j0
  for (jj in seq_len(Nco)) {
    f <- folds[[jj]]
    cM[, jj] <- f$M0[cell]
    j <- co[jj]
    kk <- key.idx[cell[, "i"]]
    mj <- f$Mj[cbind(cell[, "t"], kk)]
    obs.j <- I[cell[, "t"], j] == 1
    e <- Y[cell[, "t"], j] - mj
    e[!obs.j] <- NA
    cE[, jj] <- e
  }
  y.cell <- Y[cbind(cell[, "t"], tr[cell[, "i"]])]
  et.cell <- T.on[cbind(cell[, "t"], tr[cell[, "i"]])]
  eff.cell <- out$eff[cbind(cell[, "t"], tr[cell[, "i"]])]

  ## a linear aggregate over cells with weights w (sum 1): the fold
  ## prediction and the control's error both use the same weights
  agg <- function(sel, w, gamma) {
    if (length(sel) == 0) {
      return(c(pred.lower = NA, pred.upper = NA, CI.lower = NA,
               CI.upper = NA, p.value = NA, N = NA, q = NA))
    }
    Mbar <- colSums(cM[sel, , drop = FALSE] * w)
    Ebar <- colSums(cE[sel, , drop = FALSE] * w)
    ybar <- sum(y.cell[sel] * w)
    v <- .conformal_interval(Mbar, abs(Ebar), ybar, gamma)
    v[c("pred.lower", "pred.upper", "CI.lower", "CI.upper")] <-
      v[c("pred.lower", "pred.upper", "CI.lower", "CI.upper")] * scale
    v
  }

  ## per-unit, per-date intervals (v2, date-specific)
  unit.tab <- matrix(NA_real_, n.cell, 7)
  colnames(unit.tab) <- c("pred.lower", "pred.upper", "CI.lower", "CI.upper",
                          "p.value", "N", "q")
  for (s in seq_len(n.cell)) {
    unit.tab[s, ] <- agg(s, 1, alpha)
  }
  unit.df <- data.frame(
    unit = tr[cell[, "i"]], time = cell[, "t"], event.time = et.cell,
    Y = y.cell * scale, Y.ct = (y.cell * scale) - eff.cell, eff = eff.cell,
    CI.lower = unit.tab[, "CI.lower"], CI.upper = unit.tab[, "CI.upper"],
    p.value = unit.tab[, "p.value"], N = unit.tab[, "N"],
    q = unit.tab[, "q"], stringsAsFactors = FALSE
  )

  ## event-time ATT (extension beyond v2)
  time.on <- out$time
  att.mat <- matrix(NA_real_, length(time.on), 6,
                    dimnames = list(time.on, c("ATT", "S.E.", "CI.lower",
                                               "CI.upper", "p.value",
                                               "count")))
  att.mat[, "ATT"] <- out$att
  att.mat[, "count"] <- out$count
  att90.mat <- att.mat
  att.info <- matrix(NA_real_, length(time.on), 3,
                     dimnames = list(time.on, c("event.time", "N", "q")))
  att.info[, "event.time"] <- time.on
  for (l in seq_along(time.on)) {
    if (time.on[l] < 1) next
    sel <- which(et.cell == time.on[l])
    if (length(sel) == 0) next
    w <- rep(1 / length(sel), length(sel))
    v <- agg(sel, w, alpha)
    att.mat[l, c("CI.lower", "CI.upper", "p.value")] <-
      v[c("CI.lower", "CI.upper", "p.value")]
    att.info[l, c("N", "q")] <- v[c("N", "q")]
    v90 <- agg(sel, w, 2 * alpha)
    att90.mat[l, c("CI.lower", "CI.upper", "p.value")] <-
      v90[c("CI.lower", "CI.upper", "p.value")]
  }

  ## overall averages (extension beyond v2): cell-weighted and unit-weighted
  all.sel <- seq_len(n.cell)
  v.avg <- agg(all.sel, rep(1 / n.cell, n.cell), alpha)
  n.per.unit <- tabulate(cell[, "i"], nbins = Ntr)
  w.unit <- 1 / (sum(n.per.unit > 0) * n.per.unit[cell[, "i"]])
  v.avg.unit <- agg(all.sel, w.unit, alpha)
  est.avg <- t(as.matrix(c(out$att.avg, NA_real_,
                           v.avg[c("CI.lower", "CI.upper", "p.value")])))
  colnames(est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")
  est.avg.unit <- t(as.matrix(c(out$att.avg.unit, NA_real_,
                                v.avg.unit[c("CI.lower", "CI.upper",
                                             "p.value")])))
  colnames(est.avg.unit) <- c("ATT.avg.unit", "S.E.", "CI.lower", "CI.upper",
                              "p.value")

  ## calendar-time ATT (extension beyond v2)
  cal.eff <- as.numeric(out$eff.calendar)
  cal.N <- if (!is.null(out$N.calendar)) as.numeric(out$N.calendar) else
    rep(NA_real_, TT)
  est.cal <- cbind(cal.eff, NA_real_, NA_real_, NA_real_, NA_real_, cal.N)
  colnames(est.cal) <- c("ATT-calendar", "S.E.", "CI.lower", "CI.upper",
                         "p.value", "count")
  for (s in seq_len(TT)) {
    sel <- which(cell[, "t"] == s)
    if (length(sel) == 0 || !is.finite(cal.eff[s])) next
    v <- agg(sel, rep(1 / length(sel), length(sel)), alpha)
    est.cal[s, c("CI.lower", "CI.upper", "p.value")] <-
      v[c("CI.lower", "CI.upper", "p.value")]
  }
  cal.fit <- if (!is.null(out$eff.calendar.fit)) {
    as.numeric(out$eff.calendar.fit)
  } else {
    rep(NA_real_, TT)
  }
  est.cal.fit <- cbind(cal.fit, NA_real_, NA_real_, NA_real_, NA_real_, cal.N)
  colnames(est.cal.fit) <- colnames(est.cal)

  att.bound <- att90.mat[, c("CI.lower", "CI.upper"), drop = FALSE]
  rownames(att.bound) <- time.on

  N.used <- unit.df$N
  q.used <- unit.df$q
  status <- if (all(is.na(q.used))) {
    "no intervals"
  } else if (all(q.used[!is.na(q.used)] == 0)) {
    "infinite: too few controls for a finite interval at this level"
  } else if (any(q.used[!is.na(q.used)] == 0)) {
    "partly infinite: some intervals are the whole real line"
  } else {
    "ok"
  }

  conformal <- list(
    procedure = "pooled information-matched Jackknife+ (Liu, v2)",
    method = method,
    level = 1 - alpha,
    gamma = alpha,
    status = status,
    r = r,
    lambda = if (method == "mc") out$lambda.cv else NULL,
    n.controls = Nco,
    n.controls.ineligible = sum(rowSums(elig) == 0),
    n.failed.fits = n.fail,
    N.range = if (length(N.used)) range(N.used, na.rm = TRUE) else NA,
    unit = unit.df,
    att = data.frame(att.info, row.names = NULL),
    avg = data.frame(type = c("ATT.avg", "ATT.avg.unit"),
                     N = c(v.avg["N"], v.avg.unit["N"]),
                     q = c(v.avg["q"], v.avg.unit["q"]),
                     row.names = NULL),
    extension = paste0(
      "Intervals for est.att, est.avg, est.avg.unit and the calendar-time ",
      "ATT aggregate the treated fold predictions and the controls' signed ",
      "errors with the same cell weights; this is an extension beyond the ",
      "v2 write-up (which covers one treated unit and one date) and awaits ",
      "Ziyi Liu's sign-off. It is conservative when many units are treated."
    )
  )
  if (isTRUE(keep.sims)) {
    conformal$folds <- list(
      cell = data.frame(unit = tr[cell[, "i"]], time = cell[, "t"],
                        event.time = et.cell),
      control = co,
      M0 = cM * scale,
      E = cE * scale
    )
  }

  est.beta <- NULL
  if (!is.null(X) && !is.null(out$beta)) {
    b <- as.numeric(out$beta)
    est.beta <- cbind(b, NA_real_, NA_real_, NA_real_, NA_real_)
  }

  list(
    est.avg = est.avg,
    est.avg.unit = est.avg.unit,
    est.att = att.mat,
    est.att90 = att90.mat,
    att.bound = att.bound,
    est.eff.calendar = est.cal,
    est.eff.calendar.fit = est.cal.fit,
    est.beta = est.beta,
    vartype = "conformal",
    conformal = conformal
  )
}
