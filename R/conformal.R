###############################################################################
## vartype = "conformal": the Jackknife+ for full-path estimators
###############################################################################
##
## Implements the Jackknife+ for full-path estimators of Ziyi Liu and Yiqing
## Xu, "Unit-Level Conformal Inference for Panel Counterfactuals" (2026; the
## draft of 2026-11-05 is the reference for the section names below): the
## section "Held-out units and what they hide" (a full-path estimator hides
## each held-out unit's whole post-treatment path), the section "The
## Jackknife+" (fold j holds out unit 0 and control j, each predicted in a
## separate fit against the other controls; the calibration q_gamma; coverage
## 1 - (2q - 1)_+/N; q_gamma = 0 gives the whole real line) and the SM
## algorithm "Full-path Jackknife+ interval". GSC, IFE-EM and matrix
## completion fit the whole panel, so they are full-path estimators there.
## The p-value is the interval-dual p-value of the SM section on pre-treatment
## diagnostics, evaluated at c = 0. Three more of the paper's procedures are
## options, each in its own section below: the pre-treatment diagnostics
## (conformal.diag), the block Jackknife+ for average effects
## (conformal.block) and the relabeling test (conformal.relabel). Full
## conformal inference (for target-date estimators) is not implemented.
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
## AGGREGATED TARGETS (Yiqing, 2026-10-14; beyond the paper): intervals
## for aggregated effects (est.att by event time, est.avg, est.avg.unit and the
## calendar-time ATT). In each fold the aggregated treated query is the
## weighted average of the treated units' fold predictions (each treated unit
## predicted with its own pre-period against the same B_j). The control's
## error uses the same weights over the same cells: for a cell (i, t) the
## control is predicted over treated unit i's window, so with staggered
## adoption the control's error at event time l is the cohort-size weighted
## average of its signed errors over the cohorts' calendar windows. The
## absolute value of that weighted signed error is R_j, and the endpoint
## construction is the same as above.
## Coverage for these aggregated targets is NOT established by the unit-level
## theorem, which is about one treated unit and one pseudo-treated control's
## error. The paper (section "Multiple treated units") says so: an average of
## K prediction errors need not be exchangeable with a single control's error,
## and its block Jackknife+ (controls in blocks of K) is the construction with
## a guarantee (conformal.block = TRUE replaces these intervals with it).
## These intervals are not wide enough in general: with nine controls and 100 treated
## units whose outcomes are iid Bernoulli(0.01), a zero predictor and no
## effect, gamma = 0.1 gives N = 10 and q = 1 and exact coverage
## 1 - 0.99^9 * (1 - 0.99^100) = 0.42, not 0.90. They need separate
## assumptions and validation.
##
## Failed folds (rule D). The calibration pool of a treated cell (or of an
## aggregate) is fixed by the observation pattern alone, before any fit is
## looked at: a control is in the pool iff it is eligible for the treated
## unit's window and observed at the date (for an aggregate: in the pool of
## every aggregated cell). Controls outside the pool are dropped and N is the
## pool size plus one. A control inside the pool whose fold fit failed (a
## component fit returned NULL, or its prediction or error is non-finite) is
## NOT dropped: it keeps its slot with score R_j = +Inf, i.e. a_j = -Inf and
## b_j = +Inf, and N and q are unchanged. The interval is therefore widened or
## unchanged, never shortened. At least q failed slots give the whole real
## line. Pool exclusions (n.excluded) and failed folds (n.failed) are reported
## as separate counts, and one warning is issued when any fit failed.
##
## What the guarantee needs. All of it is an assumption about the data and the
## fits; fect checks none of it.
##  (1) Exchangeable scores. The unit-level coverage bound 1 - (2q - 1)/N holds
##      if the all-pairs score table is exchangeable given the information used.
##  (2) Rule D is valid relative to a finite reference completion. The theorem
##      is stated for a finite-valued rule. Take a reference procedure with the
##      same pool, N and q that equals the fold quantities of every successful
##      fold, supplies a finite, unit-symmetric fallback on every failed fold,
##      and whose score table satisfies (1). Replacing a failed lower endpoint
##      by -Inf can only lower each lower order statistic, and replacing a
##      failed upper endpoint by +Inf can only raise each upper one. So the
##      interval here contains the reference interval on every data set, covers
##      at least as often, and has a p-value at least as large. fect does not
##      build or check the reference completion.
##  (3) Pool by observation pattern. The pool is a function of which cells are
##      observed. The guarantee is conditional on that pattern and assumes that,
##      given it, the prediction errors of the pool controls and of the treated
##      unit are exchangeable: missingness carries no information about the size
##      of the errors beyond the pattern. If missingness depends on outcomes, or
##      on unobserved unit traits that also drive the error, coverage is not
##      established.
##  (4) Tuning. r and lambda are frozen from the original fit, which uses the
##      treated units' pre-treatment cells and the controls' full panels. The
##      paper allows tuning fixed in advance, or selected once and reused in
##      every refit if the score condition holds given the selected values
##      (SM proposition on tuning); it recommends fixing it in advance and
##      warns that tuning on the treated unit's fit can invalidate coverage.
##      With CV = TRUE, fect selects r or lambda on that fit, so the guarantee
##      rests on the condition given the selected value; CV = FALSE with a
##      given r or lambda avoids the issue. fect_boot() prints a message
##      when CV chose the value.
##  A fold fit that stops at max.iteration without converging is used as it is;
##  it is not counted as failed.
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


## Calibration q_gamma(N), the counting bound of the paper's Jackknife+
## in one exact convention for floating point:
##   q_gamma(N) = max{ q in 0..floor(N/2) : q = 0 or fl((2q - 1) / N) <= gamma }.
## 2q - 1 and N are exact integers in double arithmetic, `/` is the single
## correctly rounded IEEE division, and `gamma` is the double the caller passed,
## compared as it is: there is no tolerance and no rounding of `gamma`.
## fl((2q - 1) / N) is nondecreasing in q, so the admissible q form a prefix
## and the maximum is well defined. q = 0 means no finite interval. The
## p-value in .conformal_interval() uses the same division and the same
## comparison to gamma, so the two cannot disagree.
##  * If gamma equals a cutoff (2q - 1)/N as a rational number (0.1 at N = 10,
##    0.05 at N = 20, 0.29 at N = 100), both sides are the same double and that
##    q is admitted, as the exact definition requires.
##  * A gamma below a cutoff, however slightly (0.1 - 1e-12 at N = 10), does
##    not admit it: q = 0.
##  * The product form gamma * N is not used: 0.29 * 100 is 28.999999999999996.
##  * Duality is exact for 0 < gamma < 1: p <= gamma iff m >= N - q iff 0 lies
##    outside the effect interval. For gamma >= 1 (reached only as 2 * alpha
##    with alpha >= 0.5) q = floor(N/2) and no duality is claimed.
##  * 2 * alpha (est.att90, att.bound) is exact in floating point.
##  * An alpha that is computed rather than typed (1 - 0.9) is below 0.1 and is
##    treated as such.
## N < 1, non-finite N or non-finite gamma give 0.
.conformal_q <- function(N, gamma) {
  if (!is.finite(N) || N < 1 || !is.finite(gamma)) return(0L)
  qs <- 0:floor(N / 2)
  ok <- qs == 0L | (2 * qs - 1) / N <= gamma
  as.integer(max(qs[ok]))
}

## Jackknife+ interval and interval-dual p-value from the fold quantities.
## M0: the treated query's fold predictions (length J); R: the controls'
## absolute fold errors (length J); y: the observed treated outcome; gamma:
## miscoverage. One entry per control IN THE POOL: the caller has already
## removed controls outside the pool, and this function never drops a control.
## failed: NULL or a logical vector of length J marking pool controls whose
## fold has no usable quantity (rule D, see the file header). Such a slot is
## kept with a = -Inf and b = +Inf (score R = +Inf), so N = J + 1 and q do not
## change; M0[j] and R[j] are ignored there and may be NA. Entries with
## non-finite M0 or R are treated as failed too. The p-value uses STRICT
## comparisons (y < a, y > b): a failed slot (a = -Inf, b = +Inf) is never a
## lower- or upper-tail violation, so the p-value stays dual to the interval.
## Returns the prediction interval for Y(0), the effect interval, the p-value
## at c = 0, N, q and n.failed (the number of failed slots).
.conformal_interval <- function(M0, R, y, gamma, failed = NULL) {
  J <- length(M0)
  if (length(R) != J) {
    stop("`M0` and `R` must have the same length (", J, " vs ", length(R),
         ").", call. = FALSE)
  }
  if (is.null(failed)) {
    failed <- rep(FALSE, J)
  } else if (length(failed) != J) {
    stop("`failed` must have the same length as `M0` (", length(failed),
         " vs ", J, ").", call. = FALSE)
  }
  failed <- as.logical(failed) | !is.finite(M0) | !is.finite(R)
  N <- J + 1L
  q <- if (J >= 1L) .conformal_q(N, gamma) else 0L
  out <- c(pred.lower = -Inf, pred.upper = Inf,
           CI.lower = -Inf, CI.upper = Inf,
           p.value = 1, N = N, q = q, n.failed = sum(failed))
  if (J < 1L || !is.finite(y)) {
    out[c("p.value")] <- NA_real_
    if (!is.finite(y)) out[c("CI.lower", "CI.upper")] <- NA_real_
    return(out)
  }
  a <- M0 - R
  b <- M0 + R
  a[failed] <- -Inf
  b[failed] <- Inf
  ## interval-dual p-value at c = 0 (paper, SM pre-treatment diagnostics):
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
## Failure masks (rule D): fail0 (length Ntr) is TRUE for a treated query
## whose component fit returned NULL; failj (length K) is TRUE for a window in
## which control j was queried (window in ks) and whose component fit returned
## NULL. A window outside ks was never attempted: it is out of the pool, not
## failed. gsynth makes one call per fold, so one NULL fails every treated
## query and every window in ks. M0[, ii] / Mj[, k] stay NA for failed or
## unattempted entries.
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
  fail0 <- logical(length(tr))
  failj <- logical(K)
  if (ctx$method == "gsynth") {
    ## the factors come from B_j alone and every query is projected on its
    ## own pre-treatment history, so one call gives all the separate
    ## component predictions (no query supplies outcomes to another's fit)
    yct <- fit_cols(c(Bj, tr), copies)
    if (is.null(yct)) {
      n.fail <- 1L
      fail0[] <- TRUE
      failj[ks] <- TRUE
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
        fail0[ii] <- TRUE
      } else {
        M0[, ii] <- yct[, nB + 1L]
      }
    }
    for (kk in seq_along(ks)) {
      yct <- fit_cols(Bj, copies[kk])
      if (is.null(yct)) {
        n.fail <- n.fail + 1L
        failj[ks[kk]] <- TRUE
      } else {
        Mj[, ks[kk]] <- yct[, nB + 1L]
      }
    }
  }
  list(M0 = M0, Mj = Mj, n.fail = n.fail, fail0 = fail0, failj = failj)
}


## Event time of each date for one unit: 1 at its first treated date, NA
## where the unit is unobserved or never treated.
.conformal_t_on <- function(d, i.obs) {
  TT <- length(d)
  if (!any(d == 1)) return(rep(NA_real_, TT))
  out <- seq_len(TT) - min(which(d == 1)) + 1
  out[i.obs == 0] <- NA
  out
}

## Run a list of jobs, in parallel when a future plan is set.
.conformal_lapply <- function(n, fun, parallel) {
  use.future <- isTRUE(parallel) && n > 1 &&
    !inherits(future::plan(), "sequential")
  if (use.future) {
    future.apply::future_lapply(seq_len(n), fun, future.seed = TRUE)
  } else {
    lapply(seq_len(n), fun)
  }
}

## The fold engine of the Jackknife+: windows, eligibility, the J fold fits,
## and for every treated post-treatment cell the treated fold predictions
## (cM), the controls' signed errors (cE), pool membership and failures.
## `tun` holds method, r, lambda, validF, force, tol and max.iteration.
.conformal_engine <- function(Y, X, D, I, II, T.on, tun, parallel = FALSE) {
  TT <- nrow(Y)
  tr <- which(colSums(D * I) > 0)
  co <- which(colSums(D) == 0)
  Ntr <- length(tr)
  Nco <- length(co)

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

  ctx <- c(list(Y = Y, X = X, D = D, I = I, II = II, T.on = T.on, co = co,
                tr = tr, TT = TT, wins = wins, elig = elig), tun)
  folds <- .conformal_lapply(Nco, function(jj) .conformal_fold(jj, ctx),
                             parallel)
  n.fail <- sum(vapply(folds, function(f) f$n.fail, integer(1)))

  ## treated post-treatment cells and their fold quantities
  cell <- which(D[, tr, drop = FALSE] == 1 & I[, tr, drop = FALSE] == 1,
                arr.ind = TRUE)
  colnames(cell) <- c("t", "i")
  cell <- cell[order(cell[, "i"], cell[, "t"]), , drop = FALSE]
  n.cell <- nrow(cell)
  cM <- matrix(NA_real_, n.cell, Nco)  # treated fold predictions M_0j
  cE <- matrix(NA_real_, n.cell, Nco)  # control signed errors Y_jt - M_j0
  ## rule D: pool membership comes from the observation pattern only; a
  ## failed slot is in the pool and has no usable fold quantity
  F0 <- matrix(unlist(lapply(folds, function(f) f$fail0)), nrow = Ntr)
  Fj <- matrix(unlist(lapply(folds, function(f) f$failj)), nrow = K)
  ii.cell <- cell[, "i"]
  kk <- key.idx[ii.cell]
  inpool <- matrix(FALSE, n.cell, Nco)
  failed <- matrix(FALSE, n.cell, Nco)
  for (jj in seq_len(Nco)) {
    f <- folds[[jj]]
    cM[, jj] <- f$M0[cell]
    j <- co[jj]
    mj <- f$Mj[cbind(cell[, "t"], kk)]
    obs.j <- I[cell[, "t"], j] == 1
    e <- Y[cell[, "t"], j] - mj
    e[!obs.j] <- NA
    cE[, jj] <- e
    inpool[, jj] <- elig[jj, kk] & obs.j
    failed[, jj] <- inpool[, jj] &
      (F0[ii.cell, jj] | Fj[kk, jj] | !is.finite(cM[, jj]) |
         !is.finite(cE[, jj]))
  }
  list(ctx = ctx, tr = tr, co = co, wins = wins, key.idx = key.idx,
       elig = elig, cell = cell, cM = cM, cE = cE, inpool = inpool,
       failed = failed, F0 = F0, Fj = Fj, n.fail = n.fail,
       y.cell = Y[cbind(cell[, "t"], tr[cell[, "i"]])],
       et.cell = T.on[cbind(cell[, "t"], tr[cell[, "i"]])])
}

## Interval for a linear aggregate over cells with weights w (sum 1): the
## fold prediction and the control's error both use the same weights.
## (rule D) a control is in the aggregate's pool iff it is in the pool of
## every aggregated cell (out of pool wins over failed); among pool controls
## a failure on any aggregated cell makes the slot a = -Inf, b = +Inf.
.conformal_agg_names <- c("pred.lower", "pred.upper", "CI.lower", "CI.upper",
                          "p.value", "N", "q", "n.excluded", "n.failed")
.conformal_make_agg <- function(cM, cE, inpool, failed, y.cell, scale) {
  Nco <- ncol(cM)
  function(sel, w, gamma) {
    if (length(sel) == 0) {
      v <- rep(NA_real_, length(.conformal_agg_names))
      names(v) <- .conformal_agg_names
      return(v)
    }
    pool <- colSums(inpool[sel, , drop = FALSE]) == length(sel)
    fl <- colSums(failed[sel, pool, drop = FALSE]) > 0
    Mbar <- colSums(cM[sel, pool, drop = FALSE] * w)
    Ebar <- colSums(cE[sel, pool, drop = FALSE] * w)
    ybar <- sum(y.cell[sel] * w)
    v <- .conformal_interval(Mbar, abs(Ebar), ybar, gamma, failed = fl)
    v[c("pred.lower", "pred.upper", "CI.lower", "CI.upper")] <-
      v[c("pred.lower", "pred.upper", "CI.lower", "CI.upper")] * scale
    c(v[c("pred.lower", "pred.upper", "CI.lower", "CI.upper", "p.value",
          "N", "q")],
      n.excluded = Nco - sum(pool), n.failed = unname(v["n.failed"]))
  }
}

## Intervals for the aggregated targets: the event-time ATT (at alpha and
## 2 alpha), the overall averages and the calendar-time ATT.
.conformal_targets <- function(agg, cell, et.cell, alpha, out, TT, Ntr) {
  time.on <- out$time
  att.mat <- matrix(NA_real_, length(time.on), 6,
                    dimnames = list(time.on, c("ATT", "S.E.", "CI.lower",
                                               "CI.upper", "p.value",
                                               "count")))
  att.mat[, "ATT"] <- out$att
  att.mat[, "count"] <- out$count
  att90.mat <- att.mat
  att.info <- matrix(NA_real_, length(time.on), 5,
                     dimnames = list(time.on, c("event.time", "N", "q",
                                                "n.excluded", "n.failed")))
  att.info[, "event.time"] <- time.on
  for (l in seq_along(time.on)) {
    if (time.on[l] < 1) next
    sel <- which(et.cell == time.on[l])
    if (length(sel) == 0) next
    w <- rep(1 / length(sel), length(sel))
    v <- agg(sel, w, alpha)
    att.mat[l, c("CI.lower", "CI.upper", "p.value")] <-
      v[c("CI.lower", "CI.upper", "p.value")]
    att.info[l, c("N", "q", "n.excluded", "n.failed")] <-
      v[c("N", "q", "n.excluded", "n.failed")]
    v90 <- agg(sel, w, 2 * alpha)
    att90.mat[l, c("CI.lower", "CI.upper", "p.value")] <-
      v90[c("CI.lower", "CI.upper", "p.value")]
  }

  ## overall averages: cell-weighted and unit-weighted
  n.cell <- nrow(cell)
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

  ## calendar-time ATT
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
  list(att.mat = att.mat, att90.mat = att90.mat, att.info = att.info,
       v.avg = v.avg, v.avg.unit = v.avg.unit, est.avg = est.avg,
       est.avg.unit = est.avg.unit, est.cal = est.cal)
}


###############################################################################
## Block Jackknife+ (Liu and Xu 2026, sections "Multiple treated units" and
## "Staggered adoption"; SM algorithms "Block Jackknife+ interval for the
## average effect" and "Block Jackknife+ intervals by relative period")
###############################################################################
## The K treated units form one block. The never-treated controls observed at
## every date are split at random into G = floor(J/K) blocks of K; control
## k of block g is matched with treated unit k and takes its window (its
## pre-treatment dates and adoption date). The other controls (the leftover
## ones and those with missing outcomes) are training units in every fit.
## Fold g hides the treated units' and block g's outcomes from their dates on
## and fits ONCE on the remaining cells; that one fit predicts every hidden
## cell. For a target (a set of treated cells with weights), the fold's
## treated prediction M_g is the weighted average of the treated units'
## predictions and R_g is the absolute weighted average of the matched
## controls' errors over the same cells (control k's cell at the same date as
## treated unit k's cell). The interval then follows the single-unit rule with
## N = G + 1. Coverage is at least 1 - (2q - 1)_+/(G + 1) for the average
## effect, if the treated block and the control blocks are exchangeable given
## the information used, the partition and matching use no outcomes, and the
## tuning is fixed in advance. Effects may differ across treated units.

## Draw the partition: a list of G vectors of K column positions; element k
## of each block is matched with treated unit k.
.conformal_block_draw <- function(co, I, Ntr) {
  pool <- co[colSums(I[, co, drop = FALSE] == 1) == nrow(I)]
  G <- floor(length(pool) / Ntr)
  if (G < 1) {
    stop("conformal.block = TRUE needs at least as many never-treated ",
         "controls observed at every date as treated units (", Ntr,
         "); this panel has ", length(pool), ".", call. = FALSE)
  }
  perm <- pool[sample.int(length(pool))]
  lapply(seq_len(G), function(g) perm[(g - 1) * Ntr + seq_len(Ntr)])
}

## The G block folds for the cells of an engine `eng`. `blocks` are column
## positions in eng's panel. Returns cM, cE, inpool and failed with one
## column per block, in the format .conformal_make_agg() reads.
.conformal_block_engine <- function(eng, blocks, parallel = FALSE) {
  ctx <- eng$ctx
  tr <- eng$tr
  Ntr <- length(tr)
  TT <- ctx$TT
  G <- length(blocks)
  run_fold <- function(g) {
    Bg <- blocks[[g]]
    bg <- setdiff(ctx$co, Bg)
    copies <- lapply(seq_len(Ntr), function(k) {
      .conformal_copy(Bg[k], ctx$wins[[eng$key.idx[k]]], ctx$Y, ctx$X,
                      ctx$I, TT)
    })
    pan <- .conformal_panel(c(bg, tr), copies, ctx$Y, ctx$X, ctx$D, ctx$I,
                            ctx$II, ctx$T.on)
    yct <- .conformal_component_fit(pan, ctx$method, ctx$r, ctx$lambda,
                                    ctx$validF, ctx$force, ctx$tol,
                                    ctx$max.iteration)
    if (is.null(yct)) return(NULL)
    nb <- length(bg)
    list(T = yct[, nb + seq_len(Ntr), drop = FALSE],
         C = yct[, nb + Ntr + seq_len(Ntr), drop = FALSE])
  }
  fits <- .conformal_lapply(G, run_fold, parallel)
  cell <- eng$cell
  n.cell <- nrow(cell)
  cM <- matrix(NA_real_, n.cell, G)
  cE <- matrix(NA_real_, n.cell, G)
  failed <- matrix(TRUE, n.cell, G)
  tt <- cell[, "t"]
  ii <- cell[, "i"]
  for (g in seq_len(G)) {
    f <- fits[[g]]
    if (is.null(f)) next
    cM[, g] <- f$T[cbind(tt, ii)]
    cE[, g] <- ctx$Y[cbind(tt, blocks[[g]][ii])] - f$C[cbind(tt, ii)]
    failed[, g] <- !is.finite(cM[, g]) | !is.finite(cE[, g])
  }
  list(cM = cM, cE = cE, inpool = matrix(TRUE, n.cell, G), failed = failed,
       n.fail = sum(vapply(fits, is.null, logical(1))))
}


###############################################################################
## Pre-treatment diagnostics (Liu and Xu 2026, section "Pre-treatment
## diagnostics"; SM algorithm "Rolling pre-treatment diagnostics")
###############################################################################
## For the m dates before each treated unit's adoption (event times 1 - m to
## 0) and each such date s, the diagnostic treats s as if it were the adoption
## date: unit 0 and the held-out control show only their outcomes before s,
## their outcomes from s on are hidden, and the training controls keep their
## complete panels. The point estimate is the prediction error at s from a fit
## with all controls (the pseudo-effect, zero under no anticipation); the
## interval and the p-value are the Jackknife+ ones at s. So the diagnostic at
## event time e is the procedure above run on a panel in which each treated
## unit's adoption is moved to its date at e. A treated unit enters the
## diagnostic at e only if it is observed at that date and before it.

## Point predictions of the treated units of an engine's panel against all
## its controls (gsynth: one call; ife and mc: one fit per treated unit).
.conformal_point <- function(eng) {
  ctx <- eng$ctx
  tr <- eng$tr
  co <- eng$co
  out <- matrix(NA_real_, ctx$TT, length(tr))
  fit_cols <- function(cols) {
    pan <- .conformal_panel(cols, list(), ctx$Y, ctx$X, ctx$D, ctx$I, ctx$II,
                            ctx$T.on)
    .conformal_component_fit(pan, ctx$method, ctx$r, ctx$lambda, ctx$validF,
                             ctx$force, ctx$tol, ctx$max.iteration)
  }
  if (ctx$method == "gsynth") {
    yct <- fit_cols(c(co, tr))
    if (!is.null(yct)) out[, ] <- yct[, length(co) + seq_along(tr)]
  } else {
    for (ii in seq_along(tr)) {
      yct <- fit_cols(c(co, tr[ii]))
      if (!is.null(yct)) out[, ii] <- yct[, length(co) + 1L]
    }
  }
  out
}

.conformal_diag <- function(Y, X, D, I, T.on, m, eng0, tun, alpha, scale,
                            blocks = NULL, parallel = FALSE) {
  TT <- nrow(Y)
  tt <- seq_len(TT)
  tr <- eng0$tr
  co <- eng0$co
  A <- vapply(tr, function(i) min(which(D[, i] == 1)), numeric(1))
  unit.rows <- list()
  att.rows <- list()
  n.fail <- 0L
  for (e in (1 - m):0) {
    s <- A + e - 1
    ok <- vapply(seq_along(tr), function(ii) {
      s[ii] >= 2 && I[s[ii], tr[ii]] == 1 && any(I[seq_len(s[ii] - 1), tr[ii]] == 1)
    }, logical(1))
    if (!any(ok)) next
    cols <- c(co, tr[ok])
    De <- D[, cols, drop = FALSE]
    Ie <- I[, cols, drop = FALSE]
    Te <- T.on[, cols, drop = FALSE]
    pos <- length(co) + seq_len(sum(ok))
    for (k in seq_along(pos)) {
      De[, pos[k]] <- as.numeric(tt >= s[ok][k])
      Te[, pos[k]] <- .conformal_t_on(De[, pos[k]], Ie[, pos[k]])
    }
    IIe <- Ie * (1 - De)
    Xe <- if (is.null(X)) NULL else X[, cols, , drop = FALSE]
    eng <- .conformal_engine(Y[, cols, drop = FALSE], Xe, De, Ie, IIe, Te,
                             tun, parallel)
    n.fail <- n.fail + eng$n.fail
    s.e <- s[ok]
    sel <- which(eng$cell[, "t"] == s.e[eng$cell[, "i"]])
    pt <- .conformal_point(eng)
    agg <- .conformal_make_agg(eng$cM, eng$cE, eng$inpool, eng$failed,
                               eng$y.cell, scale)
    est <- (eng$y.cell[sel] -
              pt[cbind(eng$cell[sel, "t"], eng$cell[sel, "i"])]) * scale
    tab <- t(vapply(sel, function(x) agg(x, 1, alpha),
                    numeric(length(.conformal_agg_names))))
    colnames(tab) <- .conformal_agg_names
    unit.rows[[length(unit.rows) + 1L]] <- data.frame(
      unit = cols[eng$tr[eng$cell[sel, "i"]]], time = eng$cell[sel, "t"],
      event.time = e, Y = eng$y.cell[sel] * scale,
      Y.ct = eng$y.cell[sel] * scale - est, eff = est,
      CI.lower = tab[, "CI.lower"], CI.upper = tab[, "CI.upper"],
      p.value = tab[, "p.value"], N = tab[, "N"], q = tab[, "q"],
      n.excluded = tab[, "n.excluded"], n.failed = tab[, "n.failed"],
      stringsAsFactors = FALSE
    )
    ## the average over the treated units at e (one unit: its own interval)
    w <- rep(1 / length(sel), length(sel))
    if (!is.null(blocks)) {
      bl <- lapply(blocks, function(b) match(b[ok], cols))
      beng <- .conformal_block_engine(eng, bl, parallel)
      n.fail <- n.fail + beng$n.fail
      agg <- .conformal_make_agg(beng$cM, beng$cE, beng$inpool, beng$failed,
                                 eng$y.cell, scale)
    }
    v <- agg(sel, w, alpha)
    v90 <- agg(sel, w, 2 * alpha)
    att.rows[[length(att.rows) + 1L]] <- data.frame(
      event.time = e, ATT = mean(est), CI.lower = v[["CI.lower"]],
      CI.upper = v[["CI.upper"]], p.value = v[["p.value"]],
      CI.lower.90 = v90[["CI.lower"]], CI.upper.90 = v90[["CI.upper"]],
      p.value.90 = v90[["p.value"]], count = length(sel), N = v[["N"]],
      q = v[["q"]], n.excluded = v[["n.excluded"]],
      n.failed = v[["n.failed"]]
    )
  }
  unit.df <- if (length(unit.rows)) do.call(rbind, unit.rows) else NULL
  if (!is.null(unit.df)) {
    unit.df <- unit.df[order(unit.df$unit, unit.df$time), , drop = FALSE]
    rownames(unit.df) <- NULL
  }
  list(m = m, unit = unit.df,
       att = if (length(att.rows)) do.call(rbind, att.rows) else NULL,
       n.failed.fits = n.fail)
}


###############################################################################
## Relabeling test (Liu and Xu 2026, sections "Multiple treated units" and
## "Staggered adoption"; SM algorithm "Relabeling test, Monte Carlo p-value")
###############################################################################
## Tests the sharp null of no effect in any treated cell. Under it every
## observed outcome is untreated, so the observed panel is complete. Each draw
## gives the treated units' adoption dates to K units drawn at random from all
## units (treated units included; draw b assigns treated unit k's date to unit
## rho_b(k)), refits the point estimator on the observed panel with that
## pattern (the other units are controls) and records the statistic: the
## absolute average estimated effect over the pattern's treated cells, and by
## event time the absolute average effect at that event time. The p-value is
## (1 + #{T_b >= T_obs}) / (B + 1), exact under relabeling exchangeability:
## the actual treated units, absent treatment, are as easy to predict as any
## other K units. With random assignment it is Fisher's randomization test. A
## draw whose fit fails counts as at least as extreme as the observed value.
.conformal_relabel <- function(Y, X, D, I, tun, B, scale, time.post,
                               parallel = FALSE) {
  TT <- nrow(Y)
  N <- ncol(Y)
  tt <- seq_len(TT)
  tr <- which(colSums(D * I) > 0)
  K <- length(tr)
  A <- vapply(tr, function(i) min(which(D[, i] == 1)), numeric(1))
  stat <- function(units) {
    Db <- matrix(0, TT, N)
    Tb <- matrix(NA_real_, TT, N)
    for (k in seq_len(K)) {
      Db[, units[k]] <- as.numeric(tt >= A[k])
      Tb[, units[k]] <- .conformal_t_on(Db[, units[k]], I[, units[k]])
    }
    pan <- .conformal_panel(seq_len(N), list(), Y, X, Db, I, I * (1 - Db),
                            Tb)
    yct <- .conformal_component_fit(pan, tun$method, tun$r, tun$lambda,
                                    tun$validF, tun$force, tun$tol,
                                    tun$max.iteration)
    if (is.null(yct)) return(NULL)
    on <- Db == 1 & I == 1
    eff <- (Y[on] - yct[on]) * scale
    et <- Tb[on]
    by <- vapply(time.post, function(l) {
      x <- eff[et == l]
      if (length(x)) mean(x) else NA_real_
    }, numeric(1))
    list(avg = mean(eff), by = by)
  }
  draws <- matrix(replicate(B, sample.int(N, K)), nrow = K)
  obs <- stat(tr)
  sims <- .conformal_lapply(B, function(b) stat(draws[, b]), parallel)
  n.failed <- sum(vapply(sims, is.null, logical(1)))
  avg.b <- vapply(sims, function(x) if (is.null(x)) NA_real_ else x$avg,
                  numeric(1))
  paths <- vapply(sims, function(x) {
    if (is.null(x)) rep(NA_real_, length(time.post)) else x$by
  }, numeric(length(time.post)))
  paths <- matrix(paths, nrow = length(time.post),
                  dimnames = list(time.post, NULL))
  ## a failed or undefined draw counts as at least as extreme
  pv <- function(tb, tobs) {
    if (!is.finite(tobs)) return(NA_real_)
    (1 + sum(!is.finite(tb) | tb >= tobs)) / (B + 1)
  }
  if (is.null(obs)) {
    warning("vartype = \"conformal\": the fit with the observed treatment ",
            "pattern failed, so the relabeling test has no p-value.",
            call. = FALSE)
    obs <- list(avg = NA_real_, by = rep(NA_real_, length(time.post)))
  }
  if (n.failed > 0) {
    warning("vartype = \"conformal\": ", n.failed, " of ", B, " relabeling ",
            "fits failed; each counts as at least as extreme as the observed ",
            "statistic, so the p-values are conservative.", call. = FALSE)
  }
  att <- data.frame(
    event.time = time.post, ATT = obs$by,
    p.value = vapply(seq_along(time.post), function(l) {
      pv(abs(paths[l, ]), abs(obs$by[l]))
    }, numeric(1)),
    row.names = NULL
  )
  list(B = B, statistic = "absolute average effect over the treated cells",
       ATT.avg = obs$avg, p.value = pv(abs(avg.b), abs(obs$avg)),
       n.failed = n.failed, att = att, draws.avg = avg.b, paths = paths,
       units = draws)
}


## Caption for the gap plot of a conformal fit: what the intervals of the
## averages rest on, and where the pre-treatment intervals come from.
.conformal_plot_note <- function(cf) {
  K <- if (is.data.frame(cf$unit)) length(unique(cf$unit$unit)) else 0L
  out <- character(0)
  if (identical(cf$aggregate, "block") && !is.null(cf$block)) {
    out <- c(out, paste0("Block Jackknife+ intervals: ", cf$block$G,
                         " control blocks of ", cf$block$K, ", N = ",
                         cf$block$N, "."))
  } else if (K > 1L) {
    out <- c(out, paste0("Intervals for averages over ", K, " treated units ",
                         "have no coverage guarantee;\nconformal.block = TRUE ",
                         "gives block Jackknife+ intervals."))
  }
  if (!is.null(cf$diag) && !is.null(cf$diag$att)) {
    out <- c(out, paste0("Pre-treatment: diagnostic intervals at the last ",
                         cf$diag$m, " dates, each predicted from earlier ",
                         "dates."))
  }
  if (length(out) == 0L) return(NULL)
  paste(out, collapse = "\n")
}


## Main entry, called from fect_boot() after the point fit `out`.
## diag: number m of pre-treatment diagnostic dates (0: none); block: use the
## block Jackknife+ for the aggregated targets; relabel: number B of
## relabeling draws (0: no relabeling test).
conformal_jplus <- function(Y, X, D, I, II, T.on, out, method, force, tol,
                            max.iteration, alpha, norm.para = NULL,
                            parallel = FALSE, keep.sims = FALSE,
                            diag = 0, block = FALSE, relabel = 0) {
  TT <- nrow(Y)
  D[is.na(D)] <- 0
  scale <- if (is.null(norm.para)) 1 else norm.para[1]
  if (!is.null(X) && length(dim(X)) == 3 && dim(X)[3] == 0) X <- NULL

  co.all <- which(colSums(D) == 0)
  if (length(co.all) < 2) {
    stop("vartype = \"conformal\" needs at least two never-treated units ",
         "(the calibration pool); this panel has ", length(co.all), ".",
         call. = FALSE)
  }

  r <- if (!is.null(out$r.cv)) as.numeric(out$r.cv[1]) else 0
  ## lambda.cv is on the outcome's scale; the component fits run on the
  ## (possibly normalized) Y without norm.para, so divide as fect_mc() does
  lambda <- if (is.null(out$lambda.cv)) NULL else out$lambda.cv / scale
  validF <- if (!is.null(out$validF)) out$validF else 1
  tun <- list(method = method, r = r, lambda = lambda, validF = validF,
              force = force, tol = tol, max.iteration = max.iteration)

  eng <- .conformal_engine(Y, X, D, I, II, T.on, tun, parallel)
  tr <- eng$tr
  co <- eng$co
  Ntr <- length(tr)
  Nco <- length(co)
  n.fail <- eng$n.fail
  cell <- eng$cell
  n.cell <- nrow(cell)
  cM <- eng$cM
  cE <- eng$cE
  failed <- eng$failed
  F0 <- eng$F0
  Fj <- eng$Fj
  elig <- eng$elig
  y.cell <- eng$y.cell
  et.cell <- eng$et.cell
  eff.cell <- out$eff[cbind(cell[, "t"], tr[cell[, "i"]])]
  agg <- .conformal_make_agg(cM, cE, eng$inpool, failed, y.cell, scale)

  ## per-unit, per-date intervals
  unit.tab <- matrix(NA_real_, n.cell, length(.conformal_agg_names))
  colnames(unit.tab) <- .conformal_agg_names
  for (s in seq_len(n.cell)) {
    unit.tab[s, ] <- agg(s, 1, alpha)
  }
  unit.df <- data.frame(
    unit = tr[cell[, "i"]], time = cell[, "t"], event.time = et.cell,
    Y = y.cell * scale, Y.ct = (y.cell * scale) - eff.cell, eff = eff.cell,
    CI.lower = unit.tab[, "CI.lower"], CI.upper = unit.tab[, "CI.upper"],
    p.value = unit.tab[, "p.value"], N = unit.tab[, "N"],
    q = unit.tab[, "q"], n.excluded = unit.tab[, "n.excluded"],
    n.failed = unit.tab[, "n.failed"], stringsAsFactors = FALSE
  )

  ## aggregated targets: the pooled extension, or the block Jackknife+
  blocks <- NULL
  block.info <- NULL
  if (isTRUE(block)) {
    blocks <- .conformal_block_draw(co, I, Ntr)
    beng <- .conformal_block_engine(eng, blocks, parallel)
    agg.t <- .conformal_make_agg(beng$cM, beng$cE, beng$inpool, beng$failed,
                                 y.cell, scale)
    G <- length(blocks)
    used <- unlist(blocks)
    block.info <- list(
      G = G, K = Ntr, N = G + 1L, q = .conformal_q(G + 1L, alpha),
      partition = data.frame(block = rep(seq_len(G), each = Ntr),
                             unit = used, matched.to = rep(tr, G),
                             row.names = NULL),
      n.training = Nco - length(used),
      n.failed.fits = beng$n.fail
    )
  } else {
    agg.t <- agg
  }
  tg <- .conformal_targets(agg.t, cell, et.cell, alpha, out, TT, Ntr)
  att.mat <- tg$att.mat
  att90.mat <- tg$att90.mat
  att.info <- tg$att.info
  v.avg <- tg$v.avg
  v.avg.unit <- tg$v.avg.unit
  est.avg <- tg$est.avg
  est.avg.unit <- tg$est.avg.unit
  est.cal <- tg$est.cal

  ## pre-treatment diagnostics: fill the rows of event times 1 - m to 0
  diag.out <- NULL
  if (diag > 0) {
    diag.out <- .conformal_diag(Y, X, D, I, T.on, diag, eng, tun, alpha,
                                scale, blocks = blocks, parallel = parallel)
    da <- diag.out$att
    if (!is.null(da)) {
      rr <- match(da$event.time, out$time)
      ok <- !is.na(rr)
      rr <- rr[ok]
      da <- da[ok, , drop = FALSE]
      att.mat[rr, "ATT"] <- att90.mat[rr, "ATT"] <- da$ATT
      att.mat[rr, c("CI.lower", "CI.upper", "p.value")] <-
        as.matrix(da[, c("CI.lower", "CI.upper", "p.value")])
      att90.mat[rr, c("CI.lower", "CI.upper", "p.value")] <-
        as.matrix(da[, c("CI.lower.90", "CI.upper.90", "p.value.90")])
      att.mat[rr, "count"] <- att90.mat[rr, "count"] <- da$count
      att.info[rr, c("N", "q", "n.excluded", "n.failed")] <-
        as.matrix(da[, c("N", "q", "n.excluded", "n.failed")])
    }
  }

  cal.fit <- if (!is.null(out$eff.calendar.fit)) {
    as.numeric(out$eff.calendar.fit)
  } else {
    rep(NA_real_, TT)
  }
  cal.N <- est.cal[, "count"]
  est.cal.fit <- cbind(cal.fit, NA_real_, NA_real_, NA_real_, NA_real_, cal.N)
  colnames(est.cal.fit) <- colnames(est.cal)

  att.bound <- att90.mat[, c("CI.lower", "CI.upper"), drop = FALSE]
  rownames(att.bound) <- out$time

  ## relabeling test of the sharp null of no effect
  relabel.out <- NULL
  if (relabel > 0) {
    time.post <- sort(unique(et.cell))
    relabel.out <- .conformal_relabel(Y, X, D, I, tun, relabel, scale,
                                      time.post, parallel)
  }

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

  ## rule D bookkeeping: failures are reported apart from pool exclusions.
  ## n.failed.folds counts controls whose fold had a component fit that
  ## returned nothing (F0, Fj) OR an in-pool slot with a non-finite prediction
  ## or error (failed), so the note and the warning never say 0 while an
  ## interval is widened.
  n.failed.folds <- sum(((colSums(F0) + colSums(Fj)) > 0) |
                          (colSums(failed) > 0))
  n.intervals.affected <- sum(unit.df$n.failed > 0, na.rm = TRUE)
  any.failed <- n.fail > 0 || any(failed)
  if (any.failed) {
    status <- paste0(status, "; ", n.failed.folds, " failed fold(s), affected ",
                     "intervals widened")
    warning("vartype = \"conformal\": ", n.failed.folds, " fold(s) failed (",
            n.fail, " fit(s) returned no prediction). Failed folds are kept ",
            "in the calibration pool with infinite endpoints, so the ",
            "affected intervals are widened, not dropped, and N and q are ",
            "not reduced. ", n.intervals.affected, " reported unit-date ",
            "interval(s) are affected; see fit$conformal$unit$n.failed.",
            call. = FALSE)
  }
  n.fail.extra <- (if (is.null(block.info)) 0L else block.info$n.failed.fits) +
    (if (is.null(diag.out)) 0L else diag.out$n.failed.fits)
  if (n.fail.extra > 0) {
    warning("vartype = \"conformal\": ", n.fail.extra, " fit(s) of the ",
            "block Jackknife+ or the pre-treatment diagnostics returned no ",
            "prediction. Their folds keep infinite endpoints, so the ",
            "affected intervals are widened.", call. = FALSE)
  }

  conformal <- list(
    procedure = "Jackknife+ for full-path estimators (Liu and Xu 2026)",
    method = method,
    level = 1 - alpha,
    gamma = alpha,
    status = status,
    r = r,
    lambda = if (method == "mc") out$lambda.cv else NULL,
    ## n.controls: never-treated controls (the maximum pool);
    ## n.controls.ineligible: controls usable for no treated unit's window
    ## (a count of controls, not of cells); n.failed.fits: component fits that
    ## returned NULL; n.failed.folds: controls whose fold had a failed fit or
    ## an in-pool slot with a non-finite prediction or error (union);
    ## n.intervals.affected: rows of `unit` with n.failed > 0
    n.controls = Nco,
    n.controls.ineligible = sum(rowSums(elig) == 0),
    n.failed.fits = n.fail,
    n.failed.folds = n.failed.folds,
    n.intervals.affected = n.intervals.affected,
    N.range = if (length(N.used)) range(N.used, na.rm = TRUE) else NA,
    unit = unit.df,
    att = data.frame(att.info, row.names = NULL),
    avg = data.frame(type = c("ATT.avg", "ATT.avg.unit"),
                     N = c(v.avg["N"], v.avg.unit["N"]),
                     q = c(v.avg["q"], v.avg.unit["q"]),
                     n.excluded = c(v.avg["n.excluded"],
                                    v.avg.unit["n.excluded"]),
                     n.failed = c(v.avg["n.failed"], v.avg.unit["n.failed"]),
                     row.names = NULL),
    aggregate = if (isTRUE(block)) "block" else "pooled",
    block = block.info,
    diag = diag.out,
    relabel = relabel.out,
    extension = if (isTRUE(block)) {
      paste0(
        "Intervals for est.att, est.avg, est.avg.unit and the calendar-time ",
        "ATT come from the block Jackknife+ of Liu and Xu (2026): the treated ",
        "units form one block, the controls observed at every date are split ",
        "at random into blocks of the same size, and N = G + 1. Coverage of ",
        "the average effect holds if the treated block and the control blocks ",
        "are exchangeable and the tuning is fixed in advance."
      )
    } else {
      paste0(
        "Intervals for est.att, est.avg, est.avg.unit and the calendar-time ",
        "ATT aggregate the treated fold predictions and the controls' signed ",
        "errors with the same cell weights. This goes beyond Liu and Xu ",
        "(2026), whose guarantee is for one treated unit and one date. ",
        "Coverage for these aggregated targets is not established by the ",
        "unit-level theorem: an average of treated units' errors need not be ",
        "exchangeable with a single control's error. These intervals require ",
        "separate assumptions and validation; the paper's block Jackknife+ ",
        "(conformal.block = TRUE) has a guarantee for an average effect."
      )
    }
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
