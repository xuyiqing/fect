###############################################################################
## Full conformal inference and the Jackknife+ for target-date estimators
## (Liu and Xu 2026, sections "Full conformal inference for target-date
## estimators" and "The Jackknife+"; SM "Full conformal coverage").
###############################################################################
## A target-date estimator predicts the treated unit's untreated outcome at a
## date t from its own pre-treatment outcomes and the controls' pre-treatment
## outcomes and outcomes at t, and from no other post-treatment date.
## Synthetic control and the difference-in-differences (two-way fixed effects)
## imputation with one treated unit are examples. Such an estimator is written
## here as a RULE:
##
##   rule(query.pre, background.pre, background.target) -> prediction
##
## query.pre: the held-out unit's outcomes at the pre-treatment dates (length
## T0); background.pre: the training units' outcomes at those dates (T0 x M);
## background.target: the training units' outcomes at the target date (length
## M). The rule never sees the held-out unit's outcome at the target date.
##
## Full conformal inference (paper, Sec. 2.3). Candidate effect c: complete
## the treated unit's target outcome as Y0t - c. Hold out each unit i alone,
## hide its target outcome and predict it from all the others, the treated
## unit included with its completed outcome. Score S_i(c) = |error_i(c)|; for
## the treated unit S_0(c) = |tau_hat - c|. p(c) = #{i : S_i(c) >= S_0(c)}/N
## and the set is {c : p(c) > gamma}; its coverage is at least
## 1 - floor(gamma N)/N when the scores are exchangeable at the true effect.
##
## Exact inversion. When the rule's prediction is affine in the completed
## outcome (DID, simplex synthetic control: the weights are fitted on
## pre-treatment outcomes only), error_i(c) = a_i + b_i c, with a_0 = tau_hat
## and b_0 = -1. Then S_i(c) >= S_0(c) iff
##   ((a_i - a_0) + (b_i - b_0) c) * ((a_i + a_0) + (b_i + b_0) c) >= 0,
## a product of two linear functions, so the set is found exactly from their
## roots by a sweep: no grid. When every control's slope satisfies
## |b_i| <= 1 (the built-in rules: DID and the time mean put weight 1/(N - 1)
## on the completed outcome, the unit mean 0, simplex synthetic control a
## weight in [0, 1]), each control's acceptance set is an interval containing
## tau_hat, so the set is a single interval; a user rule with a larger weight
## can give a union of intervals. The built-in rules are affine by
## construction. A user rule is checked numerically at five candidates far on
## both sides of c = 0, and the ends of its exact set are confirmed by
## refitting just inside and just outside them; a rule that fails either
## check falls back to a grid of candidates (each refits every control).
##
## Average over dates (one treated unit). If every unit's slope b_i is the same
## at every target date (true for DID and for synthetic control, whose weights
## do not depend on the target date), the error averaged over the dates is
## abar_i + b_i cbar under any candidate path with average cbar, so the set of
## average effects is the same exact inversion with abar in place of a.
##
## Jackknife+ (paper, Sec. 2.4): hold out the treated unit together with each
## control j, predict both from the other controls; the interval follows
## .conformal_interval() with N = J + 1.
##
## The comparisons are exact (no tolerance), as in the reference code of the
## paper's simulations; ties count in favor of the candidate.
###############################################################################


## ---- built-in rules --------------------------------------------------------

## Difference in differences: the query's pre-treatment mean plus the
## background's average change from its pre-treatment means to the target
## date. With one treated unit and controls observed at every date, this is
## fect's two-way fixed-effects imputation (paper, SM lemma on saturated
## event-date TWFE).
.fc_rule_did <- function(query.pre, background.pre, background.target) {
  mean(query.pre) + mean(background.target - colMeans(background.pre))
}

## force = "unit": the query's pre-treatment mean.
.fc_rule_unit <- function(query.pre, background.pre, background.target) {
  mean(query.pre)
}

## force = "time": the background's mean at the target date.
.fc_rule_time <- function(query.pre, background.pre, background.target) {
  mean(background.target)
}

## Synthetic control with simplex weights fitted on pre-treatment outcomes,
## with a small ridge term eta * ||w||^2 so that the weights are unique. The
## defaults reproduce the reference implementation of Liu and Xu (2026).
.fc_sc_weights <- function(y, X, eta = 1e-5, tol = 1e-8) {
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("rule = \"sc\" needs the quadprog package: ",
         "install.packages(\"quadprog\").", call. = FALSE)
  }
  m <- ncol(X)
  if (m == 1L) return(1)
  gram <- crossprod(X) + eta * diag(m)
  sol <- quadprog::solve.QP(Dmat = 2 * gram,
                            dvec = 2 * drop(crossprod(X, y)),
                            Amat = cbind(rep(1, m), diag(m)),
                            bvec = c(1, rep(0, m)), meq = 1L)
  w <- as.numeric(sol$solution)
  w[abs(w) < tol] <- 0
  w / sum(w)
}

.fc_make_rule_sc <- function(eta = 1e-5) {
  force(eta)
  function(query.pre, background.pre, background.target) {
    w <- .fc_sc_weights(query.pre, background.pre, eta = eta)
    sum(w * background.target)
  }
}

.fc_resolve_rule <- function(rule, sc.eta = 1e-5) {
  if (is.function(rule)) {
    if (length(formals(rule)) < 3L && !("..." %in% names(formals(rule)))) {
      stop("A rule function must take three arguments: query.pre, ",
           "background.pre and background.target.", call. = FALSE)
    }
    return(list(fun = rule, name = "user rule", affine = FALSE))
  }
  if (!(is.character(rule) && length(rule) == 1L)) {
    stop("rule must be \"did\", \"sc\" or a function.", call. = FALSE)
  }
  ## the built-in rules are affine in the completed outcome by construction
  ## (DID and the time mean put weight 1/(N - 1) on it, the unit mean 0, and
  ## synthetic control a simplex weight fitted on pre-treatment outcomes)
  switch(rule,
         did = list(fun = .fc_rule_did, name = "did", affine = TRUE),
         sc = list(fun = .fc_make_rule_sc(sc.eta),
                   name = paste0("sc (ridge-simplex, eta = ", format(sc.eta),
                                 ")"), affine = TRUE),
         unit = list(fun = .fc_rule_unit, name = "unit", affine = TRUE),
         time = list(fun = .fc_rule_time, name = "time", affine = TRUE),
         stop("rule must be \"did\", \"sc\" or a function; got \"", rule,
              "\".", call. = FALSE))
}

## One call of the rule, checked: a single finite number or NA.
.fc_call <- function(fun, qp, bp, bt) {
  v <- tryCatch(fun(qp, bp, bt), error = function(e) NA_real_)
  v <- suppressWarnings(as.numeric(v))
  if (length(v) != 1L || !is.finite(v)) NA_real_ else v
}


## ---- studentized score: the weights ------------------------------------------

## Root mean square.
.fc_rms <- function(x) sqrt(mean(x^2))

## In-sample root mean squared residual of the query unit y (its history
## outcomes) when the training units are the columns of bp (history outcomes,
## T0 x M): the rule is called at each history date s with the training units'
## outcomes at s in the place of background.target, and sigma-hat is the root
## mean square of y - fit over the T0 dates (divisor T0). NA if any call fails.
.fc_rms_insample <- function(y, bp, fun) {
  fit <- vapply(seq_len(nrow(bp)),
                function(s) .fc_call(fun, y, bp, bp[s, ]), numeric(1))
  if (any(!is.finite(fit))) return(NA_real_)
  .fc_rms(y - fit)
}

## sigma-hat of every unit of a full conformal fit: unit j is held out and the
## training units are all the others (for the treated unit, the controls). It
## uses history dates only, so it depends on neither the target date's
## outcomes nor the candidate effect. Ypre: T0 x N history outcomes.
.fc_sigma <- function(Ypre, fun) {
  N <- ncol(Ypre)
  vapply(seq_len(N), function(j) {
    .fc_rms_insample(Ypre[, j], Ypre[, -j, drop = FALSE], fun)
  }, numeric(1))
}

## Guard on the weights. A non-finite sigma-hat is a failed fit and is left to
## the caller. A sigma-hat that is zero or tiny relative to the root mean
## square of the unit's own history outcomes (an exact pre-treatment fit) has
## no finite weight: stop. `fold`: for the Jackknife+, the control held out in
## the fold. `where`: the label of a diagnostic date, when the fit is a
## pre-treatment diagnostic (its history is the dates before that date).
.fc_check_sigma <- function(sig, yrms, labels, fold = NULL, where = NULL) {
  bad <- is.finite(sig) & sig <= 1e-8 * yrms
  if (!any(bad)) return(invisible(TRUE))
  nm <- labels[bad]
  who <- paste(utils::head(nm, 5L), collapse = ", ")
  if (length(nm) > 5L) who <- paste0(who, ", ...")
  stop("score = \"studentized\" needs a positive pre-treatment root mean ",
       "squared residual for every unit, but the fit of unit(s) ", who,
       if (!is.null(fold)) paste0(" (in the Jackknife+ fold that holds out ",
                                  "control ", fold, ")"),
       " reproduces its pre-treatment outcomes exactly",
       if (!is.null(where)) {
         paste0(" at diagnostic date ", where, " (the history is the dates ",
                "before it); use a smaller diag. Or use score = ")
       } else ". Use score = ",
       "\"absolute\", or drop that unit from the data before the call.",
       call. = FALSE)
}


## ---- exact inversion of affine scores --------------------------------------

## Sign of the linear function intercept + slope * c as c -> -Inf.
.fc_sign_minf <- function(intercept, slope) {
  if (slope > 0) return(-1L)
  if (slope < 0) return(1L)
  if (intercept > 0) return(1L)
  if (intercept < 0) return(-1L)
  0L
}

## Roots of the comparison S_i(c) = S_0(c) for every control i: a table of
## (root, unit) with one row per root of each linear factor. Controls whose
## comparison is an identity (a factor identically 0) give no roots and
## compare as TRUE everywhere.
.fc_root_table <- function(a, b, i0) {
  a0 <- a[i0]
  b0 <- b[i0]
  rr <- numeric(0)
  uu <- integer(0)
  for (j in setdiff(seq_along(a), i0)) {
    f <- rbind(c(a[j] - a0, b[j] - b0), c(a[j] + a0, b[j] + b0))
    if (any(f[, 1] == 0 & f[, 2] == 0)) next
    for (k in 1:2) {
      if (f[k, 2] != 0) {
        root <- -f[k, 1] / f[k, 2]
        if (is.finite(root)) {
          rr <- c(rr, root)
          uu <- c(uu, j)
        }
      }
    }
  }
  o <- order(rr, uu)
  data.frame(root = rr[o], unit = uu[o])
}

## Comparisons TRUE/FALSE (S_j(c) >= S_0(c)) on every open interval between
## consecutive roots and at every root, by a sweep from -Inf.
.fc_sweep <- function(a, b, i0, rt = NULL) {
  if (is.null(rt)) rt <- .fc_root_table(a, b, i0)
  a0 <- a[i0]
  b0 <- b[i0]
  ctrl <- setdiff(seq_along(a), i0)
  cur <- vapply(ctrl, function(j) {
    s1 <- .fc_sign_minf(a[j] - a0, b[j] - b0)
    s2 <- .fc_sign_minf(a[j] + a0, b[j] + b0)
    s1 == 0L || s2 == 0L || s1 == s2
  }, logical(1))
  names(cur) <- ctrl
  roots <- sort(unique(rt$root))
  int.cmp <- vector("list", length(roots) + 1L)
  root.cmp <- vector("list", length(roots))
  int.cmp[[1L]] <- cur
  for (k in seq_along(roots)) {
    rows <- rt[rt$root == roots[k], , drop = FALSE]
    at <- cur
    at[as.character(unique(rows$unit))] <- TRUE
    root.cmp[[k]] <- at
    cnt <- table(rows$unit)
    odd <- names(cnt)[cnt %% 2L == 1L]
    if (length(odd)) cur[odd] <- !cur[odd]
    int.cmp[[k + 1L]] <- cur
  }
  list(roots = roots, int.cmp = int.cmp, root.cmp = root.cmp)
}

## p-value at one candidate c for affine scores.
.fc_affine_p <- function(a, b, i0, c0) {
  sc <- abs(a + b * c0)
  sc[i0] <- abs(a[i0] + b[i0] * c0)
  mean(sc >= sc[i0])
}

## Merge sorted atoms (lower, upper, lower.closed, upper.closed) into
## disjoint components.
.fc_merge <- function(atoms) {
  if (!nrow(atoms)) return(atoms)
  atoms <- atoms[order(atoms$lower, -atoms$lower.closed), , drop = FALSE]
  out <- atoms[1L, , drop = FALSE]
  for (k in seq_len(nrow(atoms))[-1L]) {
    x <- atoms[k, , drop = FALSE]
    m <- nrow(out)
    touch <- x$lower < out$upper[m] ||
      (x$lower == out$upper[m] && (x$lower.closed || out$upper.closed[m]))
    if (touch) {
      if (x$upper > out$upper[m]) {
        out$upper[m] <- x$upper
        out$upper.closed[m] <- x$upper.closed
      } else if (x$upper == out$upper[m]) {
        out$upper.closed[m] <- out$upper.closed[m] || x$upper.closed
      }
    } else {
      out <- rbind(out, x)
    }
  }
  rownames(out) <- NULL
  out
}

## Exact full conformal set for affine errors a + b c (unit i0 is the treated
## unit, with a[i0] = tau_hat and b[i0] = -1). Returns the components, the
## hull, whether the set is disconnected, and the p-value at c = 0.
.fc_invert_affine <- function(a, b, i0, gamma) {
  N <- length(a)
  empty <- data.frame(lower = numeric(0), upper = numeric(0),
                      lower.closed = logical(0), upper.closed = logical(0))
  rt <- .fc_root_table(a, b, i0)
  sw <- .fc_sweep(a, b, i0, rt)
  roots <- sw$roots
  bd <- c(-Inf, roots, Inf)
  int.ok <- vapply(sw$int.cmp, function(cm) (1 + sum(cm)) / N > gamma,
                   logical(1))
  root.ok <- vapply(sw$root.cmp, function(cm) (1 + sum(cm)) / N > gamma,
                    logical(1))
  atoms <- empty
  for (k in seq_along(int.ok)) {
    if (int.ok[k]) {
      atoms <- rbind(atoms, data.frame(lower = bd[k], upper = bd[k + 1L],
                                       lower.closed = FALSE,
                                       upper.closed = FALSE))
    }
  }
  for (k in seq_along(root.ok)) {
    if (root.ok[k]) {
      atoms <- rbind(atoms, data.frame(lower = roots[k], upper = roots[k],
                                       lower.closed = TRUE,
                                       upper.closed = TRUE))
    }
  }
  comp <- .fc_merge(atoms)
  .fc_summary(comp, p0 = .fc_affine_p(a, b, i0, 0), N = N, gamma = gamma,
              method = "exact")
}

.fc_summary <- function(comp, p0, N, gamma, method,
                        lower.truncated = FALSE, upper.truncated = FALSE) {
  n <- nrow(comp)
  list(components = comp,
       CI.lower = if (n) min(comp$lower) else NA_real_,
       CI.upper = if (n) max(comp$upper) else NA_real_,
       disconnected = n > 1L,
       empty = n == 0L,
       p.value = p0, N = N,
       level = 1 - floor(gamma * N) / N,
       method = method,
       lower.truncated = lower.truncated,
       upper.truncated = upper.truncated)
}


## ---- one treated unit, one target date --------------------------------------

## Ypre: T0 x N pre-treatment outcomes; ytg: length-N outcomes at the target
## date; i0: the treated column; the others are the pool of controls.
## Returns the prediction coefficients: a (error at c = 0), b (slope), whether
## the rule was affine, and the treated prediction.
.fc_coefs <- function(Ypre, ytg, i0, fun, known.affine = FALSE) {
  N <- ncol(Ypre)
  ctrl <- setdiff(seq_len(N), i0)
  h <- stats::sd(ytg)
  if (!is.finite(h) || h == 0) h <- 1
  pred0 <- .fc_call(fun, Ypre[, i0], Ypre[, ctrl, drop = FALSE], ytg[ctrl])
  a <- b <- rep(NA_real_, N)
  a[i0] <- ytg[i0] - pred0
  b[i0] <- -1
  affine <- TRUE
  for (j in ctrl) {
    bg <- setdiff(seq_len(N), j)
    pos <- match(i0, bg)
    bt <- ytg[bg]
    bp <- Ypre[, bg, drop = FALSE]
    p.0 <- .fc_call(fun, Ypre[, j], bp, bt)
    bt[pos] <- ytg[i0] - h
    p.h <- .fc_call(fun, Ypre[, j], bp, bt)
    a[j] <- ytg[j] - p.0
    b[j] <- (p.0 - p.h) / h
    if (!all(is.finite(c(p.0, p.h)))) {
      affine <- FALSE
      next
    }
    ## affinity: the prediction must continue the line at candidates far on
    ## both sides (relative tolerance, so the check does not depend on the
    ## outcome's units); a rule that is affine only near c = 0 fails here or
    ## at the endpoint check in .fc_one_unit()
    if (!known.affine) {
      for (cc in c(-h, 10 * h, -10 * h, 100 * h, -100 * h)) {
        bt[pos] <- ytg[i0] - cc
        p.c <- .fc_call(fun, Ypre[, j], bp, bt)
        lin <- p.0 - b[j] * cc
        if (!is.finite(p.c) ||
            abs(p.c - lin) > 1e-8 * max(abs(p.c), abs(lin), abs(p.0), h)) {
          affine <- FALSE
          break
        }
      }
    }
  }
  list(a = a, b = b, affine = affine, pred = pred0,
       failed = !is.finite(pred0) || any(!is.finite(a)))
}

## Grid inversion for a rule that is not affine: p(c) at each candidate,
## each needing a refit of every control.
## Full conformal p-value of one candidate c by refitting every control.
## `w`: NULL for the absolute score; for the studentized score the weights
## 1/sigma-hat of the units (length N, they do not depend on c0).
.fc_pdirect <- function(Ypre, ytg, i0, fun, c0, tau, w = NULL) {
  N <- ncol(Ypre)
  yc <- ytg
  yc[i0] <- ytg[i0] - c0
  s <- numeric(N)
  s[i0] <- if (is.null(w)) abs(tau - c0) else w[i0] * abs(tau - c0)
  for (j in setdiff(seq_len(N), i0)) {
    bg <- setdiff(seq_len(N), j)
    pj <- .fc_call(fun, Ypre[, j], Ypre[, bg, drop = FALSE], yc[bg])
    s[j] <- if (!is.finite(pj)) Inf else if (is.null(w)) abs(ytg[j] - pj) else
      w[j] * abs(ytg[j] - pj)
  }
  mean(s >= s[i0])
}

.fc_grid <- function(Ypre, ytg, i0, fun, gamma, grid, tau, w = NULL) {
  N <- ncol(Ypre)
  pval <- function(c0) .fc_pdirect(Ypre, ytg, i0, fun, c0, tau, w)
  pv <- vapply(grid, pval, numeric(1))
  acc <- pv > gamma
  comp <- data.frame(lower = numeric(0), upper = numeric(0),
                     lower.closed = logical(0), upper.closed = logical(0))
  if (any(acc)) {
    r <- rle(acc)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    for (k in which(r$values)) {
      comp <- rbind(comp, data.frame(lower = grid[starts[k]],
                                     upper = grid[ends[k]],
                                     lower.closed = TRUE,
                                     upper.closed = TRUE))
    }
  }
  .fc_summary(comp, p0 = pval(0), N = N, gamma = gamma, method = "grid",
              lower.truncated = acc[1L],
              upper.truncated = acc[length(acc)])
}

## Jackknife+ for one treated unit and date: hold out the treated unit with
## each control j, predict both from the other controls. With `studentize`,
## each of the two held-out units gets its own weight from the same fold fit
## (its root mean squared in-sample residual over the history dates), and the
## fold's radius, in the treated unit's units, is E_j * s0_j / sj_j.
.fc_jplus <- function(Ypre, ytg, i0, fun, gamma, studentize = FALSE,
                      labels = NULL, where = NULL) {
  N <- ncol(Ypre)
  ctrl <- setdiff(seq_len(N), i0)
  J <- length(ctrl)
  M0 <- R <- rep(NA_real_, J)
  if (is.null(labels)) labels <- as.character(seq_len(N))
  for (k in seq_len(J)) {
    j <- ctrl[k]
    bg <- setdiff(ctrl, j)
    if (!length(bg)) next
    bp <- Ypre[, bg, drop = FALSE]
    M0[k] <- .fc_call(fun, Ypre[, i0], bp, ytg[bg])
    mj <- .fc_call(fun, Ypre[, j], bp, ytg[bg])
    R[k] <- abs(ytg[j] - mj)
    if (studentize) {
      s0 <- .fc_rms_insample(Ypre[, i0], bp, fun)
      sj <- .fc_rms_insample(Ypre[, j], bp, fun)
      .fc_check_sigma(c(s0, sj),
                      c(.fc_rms(Ypre[, i0]), .fc_rms(Ypre[, j])),
                      labels[c(i0, j)], fold = labels[j], where = where)
      R[k] <- if (is.finite(s0) && is.finite(sj)) R[k] * s0 / sj else NA_real_
    }
  }
  v <- .conformal_interval(M0, R, ytg[i0], gamma)
  list(CI.lower = unname(v["CI.lower"]), CI.upper = unname(v["CI.upper"]),
       p.value = unname(v["p.value"]), N = unname(v["N"]),
       q = unname(v["q"]), n.failed = unname(v["n.failed"]),
       level = unname(1 - max(2 * v["q"] - 1, 0) / v["N"]))
}

## Default grid for a rule that is not affine: centred on the estimate and
## wide enough to hold every control's error at c = tau_hat several times.
.fc_default_grid <- function(tau, cf, n = 401L, wt = NULL) {
  w <- if (is.null(wt)) max(abs(cf$a), na.rm = TRUE) else
    max(abs(cf$a) * wt / wt[1L], na.rm = TRUE)
  if (!is.finite(w) || w == 0) w <- 1
  seq(tau - 4 * w, tau + 4 * w, length.out = n)
}


## ---- panel driver -------------------------------------------------------------

## Wide matrices from long data: Y and D, periods x units.
.fc_wide <- function(data, yname, dname, index) {
  u <- data[[index[1]]]
  t <- data[[index[2]]]
  if (anyDuplicated(paste(u, t, sep = "\r"))) {
    stop("index does not identify the rows: some unit-period pairs repeat.",
         call. = FALSE)
  }
  units <- sort(unique(u))
  times <- sort(unique(t))
  Y <- D <- matrix(NA_real_, length(times), length(units),
                   dimnames = list(as.character(times), as.character(units)))
  pos <- cbind(match(t, times), match(u, units))
  Y[pos] <- as.numeric(data[[yname]])
  D[pos] <- as.numeric(data[[dname]])
  list(Y = Y, D = D, units = units, times = times)
}

## Treated units, their adoption dates, and the never-treated pool.
.fc_design <- function(D) {
  D0 <- D
  D0[is.na(D0)] <- 0
  if (any(!D0 %in% c(0, 1))) {
    stop("The treatment must be coded 0/1.", call. = FALSE)
  }
  tr <- which(colSums(D0) > 0)
  if (!length(tr)) stop("No treated unit in the data.", call. = FALSE)
  co <- which(colSums(D0) == 0)
  T.on <- integer(length(tr))
  for (k in seq_along(tr)) {
    d <- D[, tr[k]]
    obs <- which(!is.na(d))
    first <- obs[which(d[obs] == 1)[1L]]
    if (any(d[obs][obs >= first] == 0)) {
      stop("Full conformal inference needs a treatment that stays on once it ",
           "starts; unit ", colnames(D)[tr[k]], " switches it off.",
           call. = FALSE)
    }
    T.on[k] <- first
  }
  list(tr = tr, co = co, T.on = T.on)
}

## Refit check of an exact set's finite ends: just inside an end the direct
## p-value must exceed gamma, just outside it must not.
.fc_verify_ends <- function(v, Ypre, ytg, fun, gamma, tau, w = NULL) {
  comp <- v$components
  if (!nrow(comp)) return(TRUE)
  h <- stats::sd(ytg)
  if (!is.finite(h) || h == 0) h <- 1
  ends <- c(comp$lower, comp$upper)
  side <- c(rep(1, nrow(comp)), rep(-1, nrow(comp)))   # +1: inside is above
  for (k in seq_along(ends)) {
    e <- ends[k]
    if (!is.finite(e)) next
    eps <- 1e-6 * max(abs(e), h)
    p.in <- .fc_pdirect(Ypre, ytg, 1L, fun, e + side[k] * eps, tau, w)
    p.out <- .fc_pdirect(Ypre, ytg, 1L, fun, e - side[k] * eps, tau, w)
    if (!(p.in > gamma) || (p.out > gamma)) return(FALSE)
  }
  TRUE
}

## One treated unit k: full conformal sets (or Jackknife+ intervals) at each
## target date, and the full conformal set for the average over the post-
## treatment dates. `targets` are row positions; `history` the rows used as
## pre-treatment dates for each target (the dates before adoption for post
## dates; the dates before s for a diagnostic date s).
.fc_one_unit <- function(Y, i, co, targets, history, fun, procedure, gamma,
                         grid, known.affine = FALSE, score = "absolute",
                         diag.lab = NULL) {
  stud <- identical(score, "studentized")
  sigt <- rep(NA_real_, length(targets))
  rows <- vector("list", length(targets))
  coefs <- vector("list", length(targets))
  pools <- vector("list", length(targets))
  sigmas <- vector("list", length(targets))
  sigma.cols <- vector("list", length(targets))
  lab <- if (is.null(colnames(Y))) as.character(seq_len(ncol(Y))) else
    colnames(Y)
  for (s in seq_along(targets)) {
    tt <- targets[s]
    hh <- history[[s]]
    hh <- hh[!is.na(Y[hh, i])]
    ## pool: never-treated units observed at every history date and at tt
    ok <- vapply(co, function(j) all(!is.na(Y[c(hh, tt), j])), logical(1))
    pool <- co[ok]
    pools[[s]] <- pool
    base <- list(time.row = tt, n.history = length(hh), n.pool = length(pool))
    if (is.na(Y[tt, i]) || length(hh) < 1L || length(pool) < 1L) {
      rows[[s]] <- c(base, list(eff = NA_real_, CI.lower = NA_real_,
                                CI.upper = NA_real_, p.value = NA_real_,
                                N = length(pool) + 1, level = NA_real_,
                                disconnected = NA, method = NA_character_,
                                components = list(NULL)))
      next
    }
    cols <- c(i, pool)
    Ypre <- Y[hh, cols, drop = FALSE]
    ytg <- Y[tt, cols]
    if (procedure == "jackknife+") {
      v <- .fc_jplus(Ypre, ytg, 1L, fun, gamma, studentize = stud,
                     labels = lab[cols], where = diag.lab[s])
      if (stud) {
        ## display divisor only (never a weight): the treated unit's
        ## in-sample residual when it is predicted from all the controls,
        ## the same call as the full conformal treated-unit sigma-hat
        s0 <- .fc_rms_insample(Ypre[, 1L], Ypre[, -1L, drop = FALSE], fun)
        sigt[s] <- if (is.finite(s0) && s0 > 1e-8 * .fc_rms(Ypre[, 1L])) s0 else
          NA_real_
      }
      pred <- .fc_call(fun, Ypre[, 1L], Ypre[, -1L, drop = FALSE], ytg[-1L])
      rows[[s]] <- c(base, list(eff = ytg[1L] - pred, CI.lower = v$CI.lower,
                                CI.upper = v$CI.upper, p.value = v$p.value,
                                N = v$N, q = v$q, level = unname(v$level),
                                disconnected = FALSE, method = "jackknife+",
                                components = list(NULL)))
      next
    }
    cf <- .fc_coefs(Ypre, ytg, 1L, fun, known.affine)
    tau <- cf$a[1L]
    if (cf$failed) {
      rows[[s]] <- c(base, list(eff = tau, CI.lower = NA_real_,
                                CI.upper = NA_real_, p.value = NA_real_,
                                N = length(cols), level = NA_real_,
                                disconnected = NA, method = "failed",
                                components = list(NULL)))
      next
    }
    ## studentized score: the weights 1 / sigma-hat, from the history dates
    ## only (no outcome at the target date, no candidate effect)
    wt <- NULL
    if (stud) {
      sg <- .fc_sigma(Ypre, fun)
      sigmas[[s]] <- stats::setNames(sg, lab[cols])
      sigma.cols[[s]] <- cols
      sigt[s] <- sg[1L]
      if (any(!is.finite(sg))) {
        rows[[s]] <- c(base, list(eff = tau, CI.lower = NA_real_,
                                  CI.upper = NA_real_, p.value = NA_real_,
                                  N = length(cols), level = NA_real_,
                                  disconnected = NA, method = "failed",
                                  components = list(NULL)))
        next
      }
      .fc_check_sigma(sg, apply(Ypre, 2L, .fc_rms), lab[cols],
                      where = diag.lab[s])
      wt <- 1 / sg
    }
    v <- NULL
    if (cf$affine) {
      v <- if (stud) .fc_invert_affine(cf$a * wt, cf$b * wt, 1L, gamma) else
        .fc_invert_affine(cf$a, cf$b, 1L, gamma)
      ## a user rule that passed the affinity check: confirm each finite end
      ## of the set by refitting just inside and just outside it; if the
      ## refits disagree, the rule is not affine there and the grid is used
      if (!known.affine &&
          !.fc_verify_ends(v, Ypre, ytg, fun, gamma, tau, wt)) {
        cf$affine <- FALSE
        v <- NULL
      }
    }
    coefs[[s]] <- cf
    if (is.null(v)) {
      g <- if (is.null(grid)) .fc_default_grid(tau, cf, wt = wt) else
        tau + grid
      v <- .fc_grid(Ypre, ytg, 1L, fun, gamma, g, tau, wt)
    }
    rows[[s]] <- c(base, list(eff = tau, CI.lower = v$CI.lower,
                              CI.upper = v$CI.upper, p.value = v$p.value,
                              N = v$N, level = v$level,
                              disconnected = v$disconnected,
                              method = v$method,
                              lower.truncated = v$lower.truncated,
                              upper.truncated = v$upper.truncated,
                              components = list(v$components)))
  }
  list(rows = rows, coefs = coefs, pools = pools, sigma = sigmas,
       sigma.cols = sigma.cols, sigma.treated = sigt)
}

## Average over the post-treatment dates for one treated unit, from the
## per-date affine coefficients. Needs the same pool at every date and slopes
## that do not change across dates.
.fc_avg <- function(res, post.idx, gamma, score = "absolute") {
  cf <- res$coefs[post.idx]
  pools <- res$pools[post.idx]
  no <- function(why) list(eff = NA_real_, CI.lower = NA_real_,
                           CI.upper = NA_real_, p.value = NA_real_,
                           N = NA_real_, level = NA_real_,
                           disconnected = NA, note = why,
                           components = NULL)
  if (identical(score, "studentized")) {
    return(no("not defined for the studentized score"))
  }
  if (!length(cf) || any(vapply(cf, is.null, logical(1)))) {
    return(no("no full conformal set at some date"))
  }
  if (any(vapply(cf, function(x) x$failed || !x$affine, logical(1)))) {
    return(no("the rule is not affine or a fit failed"))
  }
  if (length(unique(vapply(pools, paste, character(1), collapse = ","))) > 1L) {
    return(no("the pool of controls differs across dates"))
  }
  B <- vapply(cf, `[[`, numeric(length(cf[[1L]]$b)), "b")
  B <- matrix(B, ncol = length(cf))
  if (any(abs(B - B[, 1L]) > 1e-10 * pmax(1, abs(B[, 1L])))) {
    return(no("the slopes differ across dates"))
  }
  A <- matrix(vapply(cf, `[[`, numeric(length(cf[[1L]]$a)), "a"),
              ncol = length(cf))
  v <- .fc_invert_affine(rowMeans(A), B[, 1L], 1L, gamma)
  list(eff = mean(A[1L, ]), CI.lower = v$CI.lower, CI.upper = v$CI.upper,
       p.value = v$p.value, N = v$N, level = v$level,
       disconnected = v$disconnected, note = "",
       components = v$components)
}

## Run the procedure for every treated unit. Returns the unit-date table, the
## per-unit averages, the diagnostics and the components.
.fc_panel <- function(Y, D, fun, procedure, gamma, diag = 0L, grid = NULL,
                      units = colnames(Y), times = rownames(Y),
                      known.affine = FALSE, score = "absolute",
                      time.names = NULL) {
  stud <- identical(score, "studentized")
  if (length(time.names) != nrow(Y)) time.names <- times
  ds <- .fc_design(D)
  tr <- ds$tr
  co <- ds$co
  if (length(co) < 1L) {
    stop("Full conformal inference needs never-treated units as controls.",
         call. = FALSE)
  }
  TT <- nrow(Y)
  tab <- list()
  dtab <- list()
  stab <- list()
  dstab <- list()
  sttab <- list()
  dsttab <- list()
  avg <- list()
  comps <- list()
  dcomps <- list()
  skipped <- character(0)
  for (k in seq_along(tr)) {
    i <- tr[k]
    on <- ds$T.on[k]
    pre <- seq_len(on - 1L)
    post <- on:TT
    post <- post[!is.na(Y[post, i])]
    if (length(pre) < 1L || all(is.na(Y[pre, i])) || !length(post)) {
      skipped <- c(skipped, as.character(units[i]))
      next
    }
    res <- .fc_one_unit(Y, i, co, post, rep(list(pre), length(post)), fun,
                        procedure, gamma, grid, known.affine, score)
    for (s in seq_along(post)) {
      r <- res$rows[[s]]
      if (stud && procedure == "full" && !is.null(res$sigma[[s]])) {
        stab[[length(stab) + 1L]] <- data.frame(
          treated = units[i], time = times[post[s]],
          unit = units[res$sigma.cols[[s]]],
          sigma = unname(res$sigma[[s]]), stringsAsFactors = FALSE)
      }
      if (stud) {
        sttab[[length(sttab) + 1L]] <- data.frame(
          treated = units[i], time = times[post[s]],
          sigma = res$sigma.treated[s], stringsAsFactors = FALSE)
      }
      comps[[length(comps) + 1L]] <- r$components[[1L]]
      tab[[length(tab) + 1L]] <- data.frame(
        unit = units[i], time = times[post[s]], event.time = post[s] - on + 1L,
        Y = Y[post[s], i], eff = r$eff, CI.lower = r$CI.lower,
        CI.upper = r$CI.upper, p.value = r$p.value, N = r$N,
        level = r$level, disconnected = r$disconnected,
        method = r$method, n.history = r$n.history,
        lower.truncated = isTRUE(r$lower.truncated),
        upper.truncated = isTRUE(r$upper.truncated),
        stringsAsFactors = FALSE)
    }
    if (procedure == "full") {
      ## the average over dates needs one pool at every date: when the pools
      ## differ, use the controls observed at every date the unit uses
      pk <- vapply(res$pools, paste, character(1), collapse = ",")
      res.avg <- res
      if (!stud && length(unique(pk)) > 1L) {
        hh <- pre[!is.na(Y[pre, i])]
        common <- co[vapply(co, function(j) all(!is.na(Y[c(hh, post), j])),
                            logical(1))]
        res.avg <- .fc_one_unit(Y, i, common, post,
                                rep(list(pre), length(post)), fun, procedure,
                                gamma, grid, known.affine)
      }
      a <- .fc_avg(res.avg, seq_along(post), gamma, score)
      avg[[length(avg) + 1L]] <- data.frame(
        unit = units[i], n.dates = length(post), eff = a$eff,
        CI.lower = a$CI.lower, CI.upper = a$CI.upper, p.value = a$p.value,
        N = a$N, level = a$level, disconnected = a$disconnected,
        note = a$note, stringsAsFactors = FALSE)
    }
    if (diag > 0) {
      ## diagnostic dates: the last `diag` pre-treatment dates; each uses the
      ## dates before it as its history
      ds.rows <- pre[pre > 1L]
      ds.rows <- utils::tail(ds.rows, diag)
      ds.rows <- ds.rows[!is.na(Y[ds.rows, i])]
      if (length(ds.rows)) {
        dres <- .fc_one_unit(Y, i, co, ds.rows,
                             lapply(ds.rows, function(s) seq_len(s - 1L)),
                             fun, procedure, gamma, grid, known.affine, score,
                             diag.lab = as.character(time.names[ds.rows]))
        for (s in seq_along(ds.rows)) {
          r <- dres$rows[[s]]
          if (stud && procedure == "full" && !is.null(dres$sigma[[s]])) {
            dstab[[length(dstab) + 1L]] <- data.frame(
              treated = units[i], time = times[ds.rows[s]],
              unit = units[dres$sigma.cols[[s]]],
              sigma = unname(dres$sigma[[s]]), stringsAsFactors = FALSE)
          }
          if (stud) {
            dsttab[[length(dsttab) + 1L]] <- data.frame(
              treated = units[i], time = times[ds.rows[s]],
              sigma = dres$sigma.treated[s], stringsAsFactors = FALSE)
          }
          dcomps[[length(dcomps) + 1L]] <- r$components[[1L]]
          dtab[[length(dtab) + 1L]] <- data.frame(
            unit = units[i], time = times[ds.rows[s]],
            event.time = ds.rows[s] - on + 1L, Y = Y[ds.rows[s], i],
            eff = r$eff, CI.lower = r$CI.lower, CI.upper = r$CI.upper,
            p.value = r$p.value, N = r$N, level = r$level,
            disconnected = r$disconnected, method = r$method,
            n.history = r$n.history,
            lower.truncated = isTRUE(r$lower.truncated),
            upper.truncated = isTRUE(r$upper.truncated),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  bind <- function(x) if (length(x)) do.call(rbind, x) else NULL
  out <- list(unit = bind(tab), avg = bind(avg), diag = bind(dtab),
              components = comps, diag.components = dcomps,
              n.treated = length(tr) - length(skipped), skipped = skipped,
              n.controls = length(co))
  if (stud && procedure == "full") out$sigma <- bind(c(stab, dstab))
  if (stud) out$sigma.treated <- bind(c(sttab, dsttab))
  out
}


## ---- user-facing function ------------------------------------------------------

#' Conformal inference for synthetic control and other target-date estimators
#'
#' Full conformal confidence sets or Jackknife+ intervals, as the user
#' chooses, for the effect on each treated unit at each post-treatment date.
#' The estimator is synthetic control by default. Any estimator that predicts
#' a treated unit's untreated outcome at a date from pre-treatment outcomes and
#' the controls' outcomes at that date (a target-date estimator) can be used
#' instead, such as difference in differences or a function the user writes.
#' The procedures are those of Liu and Xu (2026).
#'
#' @param formula A formula \code{Y ~ D} with the outcome and the 0/1
#'   treatment. Covariates are not supported.
#' @param data A data frame in long format.
#' @param index A character vector of length two: the unit and time
#'   variables.
#' @param rule The estimator: \code{"sc"} (the default: synthetic control
#'   with simplex weights fitted on the pre-treatment outcomes, with a ridge
#'   term \code{sc.eta}; needs the \pkg{quadprog} package), \code{"did"}
#'   (difference in differences: the held-out unit's pre-treatment mean plus
#'   the training units' average change from their pre-treatment means to the
#'   target date), or a function \code{function(query.pre, background.pre,
#'   background.target)} that returns the held-out unit's predicted outcome
#'   at the target date. \code{query.pre} holds the held-out unit's outcomes
#'   at the pre-treatment dates, \code{background.pre} the training units'
#'   outcomes at those dates (a matrix with one column per unit), and
#'   \code{background.target} the training units' outcomes at the target
#'   date. The function never sees the held-out unit's outcome at the target
#'   date.
#' @param procedure The inference procedure: \code{"full"} (full conformal
#'   inference, the default) or \code{"jackknife+"}. Both have a
#'   finite-sample guarantee under the same conditions. Full conformal
#'   inference is what Liu and Xu (2026) recommend for target-date
#'   estimators: its sets cover at their stated level, while Jackknife+
#'   intervals often cover more than stated and are then wider. For a rule
#'   whose prediction is not affine in the completed outcome, full conformal
#'   inference falls back to a grid and refits every control at each
#'   candidate; the Jackknife+ needs a fixed number of fits and no grid. See
#'   Details.
#' @param alpha Miscoverage level, \code{0 < alpha < 1}.
#' @param diag Number of pre-treatment dates at which to run the procedure as
#'   a placebo, each with the dates before it as its history (0 for none).
#' @param grid For a rule whose prediction is not affine in the completed
#'   outcome: candidate effects, as offsets from the estimate, at which the
#'   full conformal p-value is evaluated. By default 401 points spanning four
#'   times the largest control error on each side. Ignored for affine rules,
#'   whose sets are computed exactly.
#' @param sc.eta Ridge term of \code{rule = "sc"}.
#' @param score The score that ranks the units. \code{"absolute"} (the default) is the
#'   absolute prediction error. \code{"studentized"} is the absolute prediction
#'   error divided by the unit's pre-treatment root mean squared residual, taken
#'   from the fit that holds that unit out (see Details). The studentized score
#'   gives no set for the effect averaged over the dates: \code{avg} holds
#'   \code{NA} with a note. With \code{diag} above 0, the weights at a placebo
#'   date come from the dates before it. Choose the score before looking at the
#'   post-treatment outcomes, not by which gives the narrower set.
#' @param x,y For \code{plot}: \code{x} is a fit of class \code{"conformal_sc"}.
#'   \code{y} is an optional second fit of that class. \code{plot(fc, jp)} draws a
#'   full conformal fit and a Jackknife+ fit together, in either order; the two
#'   must differ in the procedure alone, and the call stops otherwise.
#' @param y.axis For \code{plot}: \code{"outcome"} (the default) draws the estimate and
#'   the set ends in the outcome's units. \code{"t"} divides them by the treated
#'   unit's pre-treatment root mean squared residual and needs
#'   \code{score = "studentized"}. It is not a t-statistic; see "The t axis" in
#'   Details.
#' @param xlab,ylab,main For \code{plot}: axis labels and title.
#' @param ... For \code{plot}: not used by the default plot. With
#'   \code{y.axis = "t"} or a second fit, any argument given here is named in a
#'   warning.
#'
#' @details
#' Each treated unit is analysed separately, with the never-treated units as
#' its controls. For a target date, the pool of controls holds the
#' never-treated units observed at every pre-treatment date of the treated
#' unit and at the target date; \eqn{N} is the pool size plus one.
#'
#' \strong{Full conformal inference.} For a candidate effect \eqn{c}, the
#' treated unit's untreated outcome at the target date is set to its observed
#' outcome minus \eqn{c}. Each unit is held out in turn and predicted by the
#' rule from all the other units, the treated unit included with its
#' completed outcome. With the default score, a unit's score is its absolute
#' prediction error and the treated unit's score is \eqn{|\hat\tau - c|}; the
#' studentized score divides each of these by the unit's own residual size
#' (see below). The p-value of \eqn{c} is the
#' share of units whose score is at least the treated unit's, and the set
#' keeps the candidates with p-value above \code{alpha}. It covers the effect
#' with probability at least \eqn{1 - \lfloor \alpha N \rfloor / N} if the
#' scores are exchangeable across units at the true effect: the treated unit
#' is no harder to predict than a control, and the rule and its tuning are
#' the same for every unit and fixed in advance. When the rule's prediction
#' is affine in the completed outcome, as for \code{"did"} and \code{"sc"},
#' the set is computed exactly. With the absolute score, the set for
#' \code{"did"} and \code{"sc"} is a single interval, because every control's
#' weight on the treated unit lies between 0 and 1. A user rule with a larger
#' weight can give a union of intervals (\code{disconnected = TRUE}, the parts
#' in \code{components}). With the studentized score, a set for \code{"did"} or
#' \code{"sc"} can also be a union of intervals, some of them unbounded, for
#' example when a control's pre-treatment residual is much smaller than the
#' treated unit's, so that the control's weight is large; \code{CI.lower} and
#' \code{CI.upper} then give the hull of the set, \code{disconnected} is
#' \code{TRUE}, and \code{conformal_sc} gives a message. A
#' user rule is treated as affine only if its predictions lie on a line at
#' five candidates far on both sides of zero and refits just inside and just
#' outside the ends of the set agree with it; otherwise the set is evaluated
#' on \code{grid}, and \code{lower.truncated} or \code{upper.truncated} flag
#' a set that reaches the end of the grid (\code{conformal_sc} then gives a
#' message).
#'
#' For each treated unit, the set for the effect averaged over its
#' post-treatment dates (\code{avg}) uses the errors averaged over the dates.
#' It is computed when the rule is affine with the same slope at every date
#' (true for \code{"did"} and \code{"sc"}), with the controls observed at
#' every date the unit uses. Its guarantee needs the units' absolute
#' prediction errors averaged over the dates, rather than the errors at each
#' date, to be exchangeable. It is not computed with the studentized score,
#' which is defined for each date but not for an average of errors.
#'
#' \strong{Jackknife+.} The treated unit is held out together with each
#' control, and both are predicted from the other controls. The interval
#' covers the effect with probability at least \eqn{1 - (2q - 1)/N}, where
#' \eqn{q} is the largest integer with \eqn{(2q - 1)/N \le} \code{alpha}.
#'
#' The p-value reported for each date is for the null of no effect at that
#' date.
#'
#' \strong{The studentized score.} With \code{score = "studentized"}, unit
#' \eqn{i}'s score is its absolute prediction error divided by
#' \eqn{\hat\sigma_i}, the root mean squared residual of its own fit to the
#' pre-treatment dates. To find it, the rule is applied to unit \eqn{i} at each
#' of its pre-treatment dates, with the other units in the fit as the training
#' units and their outcomes at that date in the place of the target date's.
#' \eqn{\hat\sigma_i} is the square root of the mean of the squared differences
#' between unit \eqn{i}'s outcomes and these predictions (the divisor is the
#' number of dates; there is no correction for degrees of freedom). The treated
#' unit's score at the candidate \eqn{c} is \eqn{|\hat\tau - c|/\hat\sigma_0}.
#' The weight \eqn{1/\hat\sigma_i} uses no outcome at the target date and does
#' not depend on \eqn{c}, so for an affine rule the set is still computed
#' exactly. In the Jackknife+, fold \eqn{j} predicts the treated unit and
#' control \eqn{j} from the other controls and gives each its own
#' \eqn{\hat\sigma} from that fold's fit. Control \eqn{j}'s score is its
#' absolute error divided by its \eqn{\hat\sigma_j}, so in the treated unit's
#' units the fold's radius is that error times the treated unit's
#' \eqn{\hat\sigma} divided by \eqn{\hat\sigma_j}. The interval is then formed
#' from the folds as usual. For a user rule, \eqn{\hat\sigma_i} comes from
#' calling the rule at each pre-treatment date with that date's training
#' outcomes as \code{background.target}. This is the unit's in-sample fit only
#' if the rule predicts from \code{background.target} in the same way at every
#' date, which is what the target-date interface asks. A call that returns a
#' value that is not finite makes the fit fail: full conformal inference then
#' gives no set for that date (\code{method} is \code{"failed"}), and a failed
#' Jackknife+ fold widens the interval.
#'
#' \strong{What the guarantee needs.} Both scores carry the same guarantee,
#' \eqn{1 - \lfloor \alpha N \rfloor / N} for full conformal inference and
#' \eqn{1 - (2q - 1)/N} for the Jackknife+, when the units' scores are
#' exchangeable. With the studentized score, the condition concerns each unit's
#' pre-treatment fit together with its post-treatment error: for full
#' conformal inference, the studentized scores of all \eqn{N} units at the true
#' effect must be exchangeable, and for the Jackknife+, the table of
#' studentized scores over all pairs of held-out units must be. Exchangeable
#' complete units with one common rule, fixed in advance, imply it, because
#' the same rule computes each unit's weight from its own fit. Liu and Xu
#' (2026) state this for their bands over several dates, of which a single date
#' is the case here. Choose the score before the post-treatment outcomes are
#' seen. The absolute score suits panels where prediction errors have similar
#' variances across units; when the units' pre-treatment residuals differ
#' widely, the paper recommends the studentized score. The residuals depend on
#' the pre-treatment outcomes only, and the \code{sigma} element of a
#' studentized full conformal fit reports them, so their spread can be looked
#' at before any set. Choosing the score after seeing which one gives the
#' narrower set, or which one excludes zero, voids the guarantee.
#'
#' \strong{Placebo dates.} With \code{diag} above 0, each placebo date runs the
#' same procedure on the dates before it, and with the studentized score each
#' unit's \eqn{\hat\sigma} is computed from those dates only. The paper states
#' its pre-treatment diagnostics for the absolute score. Applying them to the
#' studentized score is an extension in this package, not a result stated in the
#' paper: it uses the same argument, with each unit's \eqn{\hat\sigma}
#' computed from the history before the placebo date.
#'
#' \strong{Errors.} The studentized score needs a positive \eqn{\hat\sigma}
#' for every unit in the fit. If a unit's fit reproduces its pre-treatment
#' outcomes exactly (\eqn{\hat\sigma} is zero, or below \eqn{10^{-8}} times the
#' root mean square of that unit's own pre-treatment outcomes), the function
#' stops and names the units, and at a placebo date the date. An exact fit
#' gives an infinite weight, and leaving the unit out because of its fit would
#' make the pool depend on the data, which the guarantee does not allow. Use
#' \code{score = "absolute"}, a smaller \code{diag}, or remove the unit from the
#' data before the call, a choice made before the results are seen. For
#' example, \code{rule = "did"} with one pre-treatment date fits that date
#' exactly. The \code{plot} method stops for \code{y.axis = "t"} on a fit with
#' the absolute score, and for a pair of fits that differ in more than the
#' procedure; the message names what differs (the rule, \code{alpha}, the
#' score, the diagnostic dates or the data).
#'
#' \strong{The plot.} \code{plot(out)} draws, for each treated unit, the
#' estimate and the set at each post-treatment date, and at each placebo date
#' if \code{diag} is above 0 (to the left of the vertical line). A set that is a
#' union of intervals is drawn by its hull. \code{plot(fc, jp)} draws a full
#' conformal fit and a Jackknife+ fit at each date, the full conformal one in
#' the darker colour. \code{plot(fc, jp, y.axis = "t")} does the same on the
#' scale below.
#'
#' \strong{The t axis.} \code{y.axis = "t"} divides the estimate and the set
#' ends at each date by the treated unit's pre-treatment root mean squared
#' residual, so a value of 3 is three times that residual. For a full conformal
#' fit the divisor is the treated unit's \eqn{\hat\sigma} in that date's fit.
#' For a Jackknife+ fit it is the residual when the treated unit is predicted
#' from all the controls, a choice made for the figure, since each fold uses
#' its own; in a pair, both procedures are divided by this same number. A
#' placebo date is divided by its own residual, computed from the dates before
#' it. The division is done after the set is found, so the sets, and whether
#' they contain zero, do not change. The axis is not a standard error. Its
#' values are not t-statistics in the usual sense, and no critical value such
#' as 1.96 applies. Dates with no positive divisor are not drawn, and a
#' message names them. The default y label is \dQuote{Effect / pre-treatment
#' residual size}. The ggplot's \code{$data} holds the plotted values, and
#' with \code{y.axis = "t"} the divisor in \code{sigma}; a pair adds a
#' \code{procedure} column. For a \code{fect} fit, see \code{\link{plot.fect}}.
#'
#' @return An object of class \code{"conformal_sc"}, a list with
#'   \item{unit}{One row per treated unit and post-treatment date:
#'     \code{eff} (estimate), \code{CI.lower}, \code{CI.upper} (the set's
#'     hull), \code{p.value}, \code{N}, \code{level} (the guaranteed
#'     coverage), \code{disconnected}, \code{method} (\code{"exact"},
#'     \code{"grid"} or \code{"jackknife+"}), \code{n.history}, and
#'     \code{lower.truncated} and \code{upper.truncated} (a grid-based set
#'     that reaches the end of the grid).}
#'   \item{avg}{For \code{procedure = "full"}: one row per treated unit with
#'     the set for the effect averaged over its post-treatment dates. With
#'     \code{score = "studentized"} the set is not computed: the row holds
#'     \code{NA} and a \code{note}, "not defined for the studentized score".}
#'   \item{diag}{The pre-treatment placebo rows, if \code{diag > 0}.}
#'   \item{components}{A list with the components of each set in
#'     \code{unit}, in the same order.}
#'   \item{n.treated, n.controls, skipped}{The treated units analysed, the
#'     never-treated controls, and the treated units left out because they
#'     have no pre-treatment or post-treatment outcome.}
#'   \item{procedure, rule, alpha, score, call}{The settings.}
#'   \item{sigma, sigma.treated}{With \code{score = "studentized"}: \code{sigma.treated} has
#'     one row per treated unit and date of \code{unit} and \code{diag} (columns
#'     \code{treated}, \code{time}, \code{sigma}), the treated unit's
#'     pre-treatment root mean squared residual that \code{plot(y.axis = "t")}
#'     divides by. \code{sigma} (full conformal only; columns \code{treated},
#'     \code{time}, \code{unit}, \code{sigma}) has the residual of every unit in
#'     each date's fit, the post-treatment dates first and then the placebo dates.
#'     Use \code{x[["sigma"]]}, not \code{x$sigma}: on a Jackknife+ fit,
#'     \code{x$sigma} matches \code{sigma.treated}, because \code{$} accepts a
#'     unique prefix.}
#'
#' @references Liu, Ziyi and Yiqing Xu. 2026. "Unit-Level Conformal Inference
#'   for Panel Counterfactuals." arXiv:2611.19209.
#'   \doi{10.48550/arXiv.2611.19209}
#'
#' @seealso \code{\link{fect}} with \code{method = "fe"} and
#'   \code{vartype = "conformal"}, which uses \code{rule = "did"}.
#'
#' @examples
#' data(sim_gsynth)
#' d1 <- subset(sim_gsynth, !(id %in% 102:105))
#' ## difference in differences, full conformal sets
#' out <- conformal_sc(Y ~ D, data = d1, index = c("id", "time"),
#'                     rule = "did", alpha = 0.1)
#' out
#' head(out$unit)
#'
#' ## studentized score, both procedures, and the two together on the t axis,
#' ## with two placebo dates (the weights at a placebo date use the dates
#' ## before it)
#' fc <- conformal_sc(Y ~ D, data = d1, index = c("id", "time"), rule = "did",
#'                    alpha = 0.1, diag = 2, score = "studentized")
#' jp <- conformal_sc(Y ~ D, data = d1, index = c("id", "time"), rule = "did",
#'                    procedure = "jackknife+", alpha = 0.1, diag = 2,
#'                    score = "studentized")
#' range(fc[["sigma"]]$sigma)    ## spread of the pre-treatment residuals
#' head(fc$sigma.treated)        ## what y.axis = "t" divides by
#' plot(fc, jp, y.axis = "t")
#'
#' \donttest{
#' if (requireNamespace("quadprog", quietly = TRUE)) {
#'   ## synthetic control (the default rule): full conformal sets, then
#'   ## Jackknife+ intervals
#'   out.sc <- conformal_sc(Y ~ D, data = d1, index = c("id", "time"),
#'                          alpha = 0.1)
#'   out.sc$avg
#'   out.jp <- conformal_sc(Y ~ D, data = d1, index = c("id", "time"),
#'                          procedure = "jackknife+", alpha = 0.1)
#'   head(out.jp$unit)
#' }
#' ## a user rule: the mean of the training units at the target date, shifted
#' ## by the held-out unit's pre-treatment gap to them
#' my_rule <- function(query.pre, background.pre, background.target) {
#'   mean(background.target) + mean(query.pre - rowMeans(background.pre))
#' }
#' out.my <- conformal_sc(Y ~ D, data = d1, index = c("id", "time"),
#'                        rule = my_rule)
#' }
#' @export
conformal_sc <- function(formula, data, index, rule = "sc",
                         procedure = c("full", "jackknife+"), alpha = 0.1,
                         diag = 0, grid = NULL, sc.eta = 1e-5,
                         score = c("absolute", "studentized")) {
  cl <- match.call()
  procedure <- match.arg(procedure)
  score <- match.arg(score)
  if (!(is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
        alpha > 0 && alpha < 1)) {
    stop("alpha must be a number with 0 < alpha < 1.", call. = FALSE)
  }
  if (!(is.numeric(diag) && length(diag) == 1L && is.finite(diag) &&
        diag >= 0 && diag == round(diag))) {
    stop("diag must be a whole number of 0 or more.", call. = FALSE)
  }
  if (!is.null(grid) && !(is.numeric(grid) && all(is.finite(grid)) &&
                          length(grid) >= 2L)) {
    stop("grid must be a numeric vector of at least two finite offsets.",
         call. = FALSE)
  }
  if (!is.null(grid)) grid <- sort(unique(grid))
  if (!(is.character(index) && length(index) == 2L &&
        all(index %in% names(data)))) {
    stop("index must name the unit and time variables of data.",
         call. = FALSE)
  }
  vars <- all.vars(formula)
  if (length(vars) != 2L) {
    stop("formula must be Y ~ D; covariates are not supported.",
         call. = FALSE)
  }
  if (!all(vars %in% names(data))) {
    stop("formula variables not found in data: ",
         paste(setdiff(vars, names(data)), collapse = ", "), ".",
         call. = FALSE)
  }
  data <- as.data.frame(data)
  if (anyNA(data[[index[1]]]) || anyNA(data[[index[2]]])) {
    stop("The index variables must not have missing values.", call. = FALSE)
  }
  ru <- .fc_resolve_rule(rule, sc.eta)
  w <- .fc_wide(data, vars[1], vars[2], index)
  res <- .fc_panel(w$Y, w$D, ru$fun, procedure, alpha, diag = as.integer(diag),
                   grid = grid, units = w$units, times = w$times,
                   known.affine = isTRUE(ru$affine), score = score)
  if (length(res$skipped)) {
    message("Treated unit(s) without a pre-treatment or post-treatment ",
            "outcome are left out: ", paste(res$skipped, collapse = ", "), ".")
  }
  tr.flag <- c(res$unit$lower.truncated, res$unit$upper.truncated,
               res$diag$lower.truncated, res$diag$upper.truncated)
  if (any(tr.flag %in% TRUE)) {
    message("Some grid-based sets reach the end of the grid ",
            "(lower.truncated or upper.truncated is TRUE): the true set extends ",
            "beyond it. Widen `grid`.")
  }
  if (!is.null(res$unit) && procedure == "full" &&
      any(res$unit$disconnected %in% TRUE)) {
    message("Some full conformal sets are unions of intervals; CI.lower and ",
            "CI.upper give their hull. See $components.")
  }
  structure(c(res, list(procedure = procedure, rule = ru$name,
                        alpha = alpha, score = score, call = cl)),
            class = "conformal_sc")
}

#' @export
print.conformal_sc <- function(x, digits = 3, ...) {
  cat("Call:\n")
  print(x$call)
  lab <- if (x$procedure == "full") "Full conformal inference" else
    "Jackknife+"
  cat("\n", lab, " (Liu and Xu 2026), rule: ", x$rule, ", alpha = ",
      x$alpha, "\n", sep = "")
  if (identical(x$score, "studentized")) {
    cat("Score: studentized (absolute prediction error divided by the ",
        "unit's pre-treatment root mean squared residual)\n", sep = "")
  }
  cat("Treated units: ", x$n.treated, "; never-treated controls: ",
      x$n.controls, "\n", sep = "")
  if (!is.null(x$unit) && nrow(x$unit)) {
    lv <- range(x$unit$level, na.rm = TRUE)
    if (all(is.finite(lv))) {
      cat("Guaranteed coverage of each date's set: ",
          if (lv[1] == lv[2]) format(round(lv[1], 3)) else
            paste(format(round(lv, 3)), collapse = " to "), "\n", sep = "")
    }
    cat("\nBy unit and date:\n")
    u <- x$unit[, c("unit", "time", "event.time", "eff", "CI.lower",
                    "CI.upper", "p.value", "N")]
    print(utils::head(u, 12), digits = digits, row.names = FALSE)
    if (nrow(u) > 12) cat("... ", nrow(u) - 12, " more rows in $unit\n",
                          sep = "")
  }
  if (!is.null(x$avg) && nrow(x$avg) && identical(x$score, "studentized")) {
    cat("\nAverage over the post-treatment dates: not computed with the ",
        "studentized score\n", sep = "")
  } else if (!is.null(x$avg) && nrow(x$avg)) {
    cat("\nAverage over the post-treatment dates, by unit:\n")
    print(x$avg[, c("unit", "n.dates", "eff", "CI.lower", "CI.upper",
                    "p.value", "N")], digits = digits, row.names = FALSE)
  }
  if (!is.null(x$diag) && nrow(x$diag)) {
    nx <- sum(x$diag$p.value <= x$alpha, na.rm = TRUE)
    cat("\nPre-treatment placebo dates: ", nx, " of ", nrow(x$diag),
        " sets exclude zero\n", sep = "")
  }
  invisible(x)
}

## ---- event-study plot: outcome units or the studentized (t) scale ------------

## The columns drawn, and their rows: post-treatment dates, then diagnostic
## dates. `f` is a conformal_sc fit or a list with the same pieces (unit,
## diag, sigma.treated, procedure, score).
.fc_plot_cols <- c("unit", "time", "event.time", "eff", "CI.lower", "CI.upper")

.fc_plot_rows <- function(f) {
  u <- f$unit[, .fc_plot_cols, drop = FALSE]
  n.diag <- 0L
  if (!is.null(f$diag) && nrow(f$diag)) {
    u <- rbind(u, f$diag[, .fc_plot_cols, drop = FALSE])
    n.diag <- nrow(f$diag)
  }
  rownames(u) <- NULL
  list(rows = u, n.diag = n.diag)
}

## Same-length numeric vectors equal to a relative tolerance, NA in the same
## places.
.fc_same_num <- function(a, b, tol) {
  if (length(a) != length(b) || !identical(is.na(a), is.na(b))) return(FALSE)
  k <- !is.na(a)
  all(abs(a[k] - b[k]) <= tol * pmax(1, abs(a[k])))
}

## The divisor of each row, from the fit's sigma.treated; rows without a
## positive finite divisor are dropped (one message names them).
.fc_plot_divide <- function(rows, st, quiet = FALSE) {
  key <- paste(rows$unit, rows$time, sep = "\r")
  rows$sigma <- st$sigma[match(key, paste(st$treated, st$time, sep = "\r"))]
  bad <- !(is.finite(rows$sigma) & rows$sigma > 0)
  if (any(bad)) {
    if (!quiet) {
      nm <- if (length(unique(rows$unit)) > 1L) {
        paste0(rows$unit[bad], ": ", rows$time[bad])
      } else as.character(rows$time[bad])
      message("No divisor (a positive pre-treatment residual) at: ",
              paste(nm, collapse = ", "), "; not drawn.")
    }
    rows <- rows[!bad, , drop = FALSE]
  }
  for (v in c("eff", "CI.lower", "CI.upper")) rows[[v]] <- rows[[v]] / rows$sigma
  rows
}

## Errors for y.axis = "t" on a fit: needs the studentized score and the
## stored divisor.
.fc_plot_need_t <- function(f, what = "score = \"studentized\"") {
  if (!identical(f$score, "studentized")) {
    stop("y.axis = \"t\" needs a fit with the studentized score (", what,
         ").", call. = FALSE)
  }
  if (is.null(f$sigma.treated)) {
    stop("The fit stores no pre-treatment residual for the treated unit; ",
         "refit with the current version.", call. = FALSE)
  }
  invisible(TRUE)
}

## The two fits of plot(x, y) must differ in the procedure alone.
.fc_pair_check <- function(x, y) {
  if (!inherits(y, "conformal_sc")) {
    stop("The second argument of plot(x, y) must be a conformal_sc fit.",
         call. = FALSE)
  }
  if (is.null(y$unit) || !nrow(y$unit)) stop("Nothing to plot.", call. = FALSE)
  pr <- c(x$procedure, y$procedure)
  if (!setequal(pr, c("full", "jackknife+"))) {
    stop("plot(x, y) needs one full conformal fit and one Jackknife+ fit.",
         call. = FALSE)
  }
  sc <- function(f) if (is.null(f$score)) "absolute" else f$score
  diff <- character(0)
  if (!identical(x$rule, y$rule)) diff <- c(diff, "rule")
  if (!isTRUE(abs(x$alpha - y$alpha) < 1e-12)) diff <- c(diff, "alpha")
  if (!identical(sc(x), sc(y))) diff <- c(diff, "score")
  dg <- function(f) if (is.null(f$diag)) character(0) else
    paste(f$diag$unit, f$diag$time)
  if (!identical(dg(x), dg(y))) diff <- c(diff, "diagnostic dates")
  key <- function(f) paste(f$unit$unit, f$unit$time, f$unit$event.time)
  if (!identical(key(x), key(y)) ||
      !.fc_same_num(x$unit$Y, y$unit$Y, 1e-10) ||
      (identical(x$rule, y$rule) &&
       !.fc_same_num(x$unit$eff, y$unit$eff, 1e-10))) diff <- c(diff, "data")
  if (length(diff)) {
    stop("The two fits differ in ", paste(diff, collapse = ", "),
         "; plot them together only when they differ in procedure alone.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.fc_plot_xlab <- "Periods since treatment"
.fc_plot_ylab_t <- "Effect / pre-treatment residual size"

## Caption of the studentized (t) plot and of paired plots: at most two short
## sentences, wrapped to fit a figure 6.5 inches wide.
.fc_plot_caption <- function(tmode, paired, procedure, has.diag, differs, hull) {
  cap <- character(0)
  if (tmode) {
    cap <- paste0("Effects and sets are divided by the treated unit's ",
                  "pre-treatment root mean squared residual.")
  }
  cl <- character(0)
  if (has.diag) cl <- c(cl, "dates before the line are placebo dates")
  if (tmode && (paired || identical(procedure, "jackknife+"))) {
    cl <- c(cl, "the divisor uses all controls")
  }
  if (hull) cl <- c(cl, "split sets are drawn as their hull")
  if (length(cl)) {
    cl <- paste(cl, collapse = "; ")
    cap <- c(cap, paste0(toupper(substr(cl, 1L, 1L)), substring(cl, 2L), "."))
  }
  if (length(cap)) {
    paste(strwrap(paste(cap, collapse = " "), width = 90), collapse = "\n")
  } else NULL
}

## One drawing function for a single fit and for a pair. `a` is a fit (the
## full conformal one when `b` is given), `b` the Jackknife+ fit or NULL.
## Returns the ggplot; its $data holds the plotted columns.
.fc_plot_events <- function(a, b = NULL, y.axis = "outcome", xlab = NULL,
                            ylab = NULL, main = NULL, xlim = NULL,
                            ylim = NULL) {
  tmode <- identical(y.axis, "t")
  paired <- !is.null(b)
  if (!tmode && !paired) {
    ## the default: dev's plot, layer for layer
    u <- a$unit
    if (!is.null(a$diag)) u <- rbind(u, a$diag)
    d <- u
    cap <- NULL
  } else {
    ra <- .fc_plot_rows(a)
    ra.rows <- ra$rows
    if (tmode) ra.rows <- .fc_plot_divide(ra.rows, a$sigma.treated)
    if (paired) {
      rb <- .fc_plot_rows(b)$rows
      if (tmode) rb <- .fc_plot_divide(rb, a$sigma.treated, quiet = TRUE)
      ra.rows$procedure <- "Full conformal"
      rb$procedure <- "Jackknife+"
      ## the Jackknife+ rows come first, so the full conformal set is
      ## drawn on top
      d <- rbind(rb, ra.rows)
      d$procedure <- factor(d$procedure,
                            levels = c("Full conformal", "Jackknife+"))
      nm <- c("unit", "time", "event.time", "procedure", "eff", "CI.lower",
              "CI.upper", if (tmode) "sigma")
    } else {
      d <- ra.rows
      nm <- c(.fc_plot_cols, if (tmode) "sigma")
    }
    d <- d[, nm, drop = FALSE]
    rownames(d) <- NULL
    if (!nrow(d)) stop("Nothing to plot.", call. = FALSE)
    differs <- FALSE
    if (tmode) {
      post <- ra.rows[ra.rows$event.time >= 1, , drop = FALSE]
      differs <- any(vapply(split(post$sigma, post$unit), function(v) {
        length(v) > 1L && diff(range(v)) > 1e-10 * max(abs(v))
      }, logical(1)))
    }
    full <- if (paired) a else if (identical(a$procedure, "full")) a else NULL
    hull <- !is.null(full) &&
      any(c(full$unit$disconnected, full$diag$disconnected) %in% TRUE)
    cap <- .fc_plot_caption(tmode, paired, a$procedure, ra$n.diag > 0L,
                            differs, hull)
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$event.time,
                                       y = .data$eff)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_vline(xintercept = 0.5, colour = "grey80")
  if (paired) {
    p <- p +
      ggplot2::geom_pointrange(
        ggplot2::aes(ymin = .data$CI.lower, ymax = .data$CI.upper,
                     colour = .data$procedure, group = .data$procedure),
        size = 0.2, position = ggplot2::position_dodge(width = 0.5)) +
      ggplot2::scale_colour_manual(
        values = c("Full conformal" = "#004a8c", "Jackknife+" = "#8fb4d9"),
        breaks = c("Full conformal", "Jackknife+"), name = NULL)
  } else {
    p <- p + ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$CI.lower,
                                                   ymax = .data$CI.upper),
                                      size = 0.2)
  }
  p <- p + ggplot2::facet_wrap(~ unit) +
    ggplot2::labs(x = if (is.null(xlab)) .fc_plot_xlab else xlab,
                  y = if (!is.null(ylab)) ylab else if (tmode)
                    .fc_plot_ylab_t else "Effect") +
    ggplot2::theme_bw()
  if (paired) p <- p + ggplot2::theme(legend.position = "bottom")
  if (!is.null(main)) p <- p + ggplot2::labs(title = main)
  if (!is.null(cap)) {
    ## long caption lines run from the left edge of the figure
    p <- p + ggplot2::labs(caption = cap) +
      ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0),
                     plot.caption.position = "plot")
  }
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }
  p
}

## plot.fect(y.axis = "t"): the same figure as plot.conformal_sc(), built from
## fit$conformal. `cl` is the match.call() of plot.fect: any argument that is
## not used here is named in one warning.
.fc_plot_fect_t <- function(x, type, cl, main, xlab, ylab, xlim, ylim,
                            return.data) {
  cf <- x[["conformal"]]
  if (!identical(x[["vartype"]], "conformal") || is.null(cf)) {
    stop("y.axis = \"t\" needs a fit with vartype = \"conformal\".",
         call. = FALSE)
  }
  if (!identical(cf$method, "fe")) {
    stop("y.axis = \"t\" needs a fit with method = \"fe\" and vartype = ",
         "\"conformal\" (full conformal inference); gsynth, ife and mc fits ",
         "have no studentized score.", call. = FALSE)
  }
  if (!(is.null(type) || identical(type, "gap") || identical(type, "es"))) {
    stop("y.axis = \"t\" is available for type = \"gap\" only.", call. = FALSE)
  }
  f <- list(unit = cf$unit, diag = cf$diag$unit,
            sigma.treated = cf$sigma.treated, procedure = "full",
            score = cf$score)
  .fc_plot_need_t(f, "conformal.score = \"studentized\"")
  used <- c("x", "type", "y.axis", "main", "xlab", "ylab", "xlim", "ylim",
            "return.data")
  extra <- setdiff(names(as.list(cl))[-1L], used)
  nun <- sum(!nzchar(names(as.list(cl))[-1L]))
  extra <- c(extra[nzchar(extra)], rep("(unnamed)", nun))
  if (length(extra)) {
    warning(paste(extra, collapse = ", "), ": not used with y.axis = \"t\".",
            call. = FALSE)
  }
  p <- .fc_plot_events(f, NULL, "t", xlab, ylab, main, xlim, ylim)
  if (isTRUE(return.data)) return(p$data)
  p
}

#' @rdname conformal_sc
#' @export
plot.conformal_sc <- function(x, y = NULL, y.axis = c("outcome", "t"),
                              xlab = NULL, ylab = NULL, main = NULL, ...) {
  y.axis <- match.arg(y.axis)
  if (is.null(x$unit) || !nrow(x$unit)) {
    stop("Nothing to plot.", call. = FALSE)
  }
  new <- identical(y.axis, "t") || !is.null(y)
  if (new && ...length()) {
    dn <- names(list(...))
    if (is.null(dn)) dn <- rep("", ...length())
    dn[!nzchar(dn)] <- "(unnamed)"
    warning(paste(dn, collapse = ", "),
            ": not used with y.axis = \"t\" or a second fit.", call. = FALSE)
  }
  b <- NULL
  if (!is.null(y)) {
    .fc_pair_check(x, y)
    if (!identical(x$procedure, "full")) {
      b <- x
      x <- y
    } else {
      b <- y
    }
  }
  if (identical(y.axis, "t")) {
    .fc_plot_need_t(x)
    if (!is.null(b)) .fc_plot_need_t(b)
    if (!is.null(b) &&
        !.fc_same_num(x$sigma.treated$sigma, b$sigma.treated$sigma, 1e-12)) {
      stop("The two fits differ in data; plot them together only when they ",
           "differ in procedure alone.", call. = FALSE)
    }
  }
  .fc_plot_events(x, b, y.axis, xlab, ylab, main)
}


## ---- fect(method = "fe", vartype = "conformal") -------------------------------

## Called from fect_boot() after the point fit `out`. Y, D, I: TT x N
## matrices (Y possibly normalized by norm.para); force: 1 unit, 2 time,
## 3 two-way. The rule is the target-date form of the fixed-effects
## imputation: for force = 3 the difference in differences between the
## treated unit's pre-treatment mean and the controls' change to the target
## date. With one treated unit and controls observed at every date it equals
## fect's estimate; otherwise the two can differ, and the intervals are built
## around the rule's estimate (column eff; fect's estimate is in eff.fit).
conformal_full_fe <- function(Y, D, I, T.on, out, force, alpha,
                              norm.para = NULL, diag = 0,
                              score = "absolute", time.lab = NULL) {
  scale <- if (is.null(norm.para)) 1 else norm.para[1]
  TT <- nrow(Y)
  N <- ncol(Y)
  Yr <- Y * scale
  Yr[I == 0] <- NA
  Dm <- D
  Dm[I == 0] <- NA
  rule <- switch(as.character(force),
                 "1" = .fc_rule_unit, "2" = .fc_rule_time,
                 "3" = .fc_rule_did,
                 stop("vartype = \"conformal\" with method = \"fe\" needs ",
                      "force = \"unit\", \"time\" or \"two-way\".",
                      call. = FALSE))
  rule.name <- switch(as.character(force), "1" = "unit mean",
                      "2" = "time mean", "3" = "difference in differences")
  ## `time.lab`: the user's time labels; they only name dates in messages
  ## (the table's time column stays the position)
  res <- .fc_panel(Yr, Dm, rule, "full", alpha, diag = diag,
                   units = seq_len(N), times = seq_len(TT),
                   time.names = time.lab,
                   known.affine = TRUE, score = score)
  res90 <- .fc_panel(Yr, Dm, rule, "full", 2 * alpha, diag = diag,
                     units = seq_len(N), times = seq_len(TT),
                     time.names = time.lab,
                     known.affine = TRUE, score = score)
  u <- res$unit
  if (length(res$skipped)) {
    message("vartype = \"conformal\": treated unit(s) without a ",
            "pre-treatment or post-treatment outcome get no set.")
  }
  if (is.null(u) || !nrow(u)) {
    stop("vartype = \"conformal\": no treated unit has a pre-treatment date ",
         "and a post-treatment outcome.", call. = FALSE)
  }
  tr.units <- sort(unique(u$unit))
  Ntr <- length(tr.units)
  eff.fit <- out$eff[cbind(u$time, u$unit)]
  dev <- abs(u$eff - eff.fit)
  if (any(dev > 1e-6 * (1 + abs(eff.fit)), na.rm = TRUE)) {
    message("vartype = \"conformal\" with method = \"fe\": the sets are built ",
            "around the target-date (difference-in-differences) estimate ",
            "with the never-treated units as controls, which differs from ",
            "fect's estimate in ", sum(dev > 1e-6 * (1 + abs(eff.fit)),
                                       na.rm = TRUE),
            " cell(s): other treated units or not-yet-treated cells enter ",
            "fect's fit, or controls with missing outcomes at the dates used ",
            "are dropped from the pool. Both estimates are in ",
            "fit$conformal$unit (eff, eff.fit); with one treated unit, the ",
            "ATT rows of est.att and est.avg hold the target-date estimate.")
  }
  unit.df <- data.frame(
    unit = u$unit, time = u$time, event.time = u$event.time, Y = u$Y,
    Y.ct = u$Y - u$eff, eff = u$eff, eff.fit = eff.fit,
    CI.lower = u$CI.lower, CI.upper = u$CI.upper, p.value = u$p.value,
    N = u$N, level = u$level, disconnected = u$disconnected,
    stringsAsFactors = FALSE)

  ## event-time and overall rows: with one treated unit they are that unit's
  ## sets; with several, the averages have no full conformal set here
  time.on <- out$time
  mk <- function(nm) {
    m <- matrix(NA_real_, length(time.on), 6,
                dimnames = list(time.on, c(nm, "S.E.", "CI.lower", "CI.upper",
                                           "p.value", "count")))
    m[, nm] <- out$att
    m[, "count"] <- out$count
    m
  }
  att.mat <- mk("ATT")
  att90.mat <- mk("ATT")
  est.avg <- t(as.matrix(c(out$att.avg, NA_real_, NA_real_, NA_real_,
                           NA_real_)))
  colnames(est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")
  est.avg.unit <- t(as.matrix(c(out$att.avg.unit, NA_real_, NA_real_,
                                NA_real_, NA_real_)))
  colnames(est.avg.unit) <- c("ATT.avg.unit", "S.E.", "CI.lower", "CI.upper",
                              "p.value")
  cal.eff <- as.numeric(out$eff.calendar)
  cal.N <- if (!is.null(out$N.calendar)) as.numeric(out$N.calendar) else
    rep(NA_real_, TT)
  est.cal <- cbind(cal.eff, NA_real_, NA_real_, NA_real_, NA_real_, cal.N)
  colnames(est.cal) <- c("ATT-calendar", "S.E.", "CI.lower", "CI.upper",
                         "p.value", "count")
  diag.att <- NULL
  if (Ntr == 1L) {
    u90 <- res90$unit
    for (l in seq_along(time.on)) {
      k <- which(u$event.time == time.on[l])
      if (length(k) != 1L) next
      att.mat[l, "ATT"] <- att90.mat[l, "ATT"] <- u$eff[k]
      att.mat[l, c("CI.lower", "CI.upper", "p.value")] <-
        unlist(u[k, c("CI.lower", "CI.upper", "p.value")])
      k90 <- which(u90$event.time == time.on[l])
      att90.mat[l, c("CI.lower", "CI.upper", "p.value")] <-
        unlist(u90[k90, c("CI.lower", "CI.upper", "p.value")])
    }
    for (k in seq_len(nrow(u))) {
      est.cal[u$time[k], c("CI.lower", "CI.upper", "p.value")] <-
        unlist(u[k, c("CI.lower", "CI.upper", "p.value")])
    }
    a <- res$avg
    if (!is.null(a) && nrow(a) == 1L && is.finite(a$eff)) {
      est.avg[1, ] <- c(a$eff, NA_real_, a$CI.lower, a$CI.upper, a$p.value)
      est.avg.unit[1, ] <- est.avg[1, ]
    }
    if (!is.null(res$diag)) {
      dg <- res$diag
      dg90 <- res90$diag
      for (k in seq_len(nrow(dg))) {
        l <- match(dg$event.time[k], time.on)
        if (is.na(l)) next
        att.mat[l, "ATT"] <- att90.mat[l, "ATT"] <- dg$eff[k]
        att.mat[l, c("CI.lower", "CI.upper", "p.value")] <-
          unlist(dg[k, c("CI.lower", "CI.upper", "p.value")])
        att90.mat[l, c("CI.lower", "CI.upper", "p.value")] <-
          unlist(dg90[k, c("CI.lower", "CI.upper", "p.value")])
      }
      diag.att <- dg
    }
  }
  est.cal.fit <- cbind(if (!is.null(out$eff.calendar.fit))
    as.numeric(out$eff.calendar.fit) else rep(NA_real_, TT),
    NA_real_, NA_real_, NA_real_, NA_real_, cal.N)
  colnames(est.cal.fit) <- colnames(est.cal)
  att.bound <- att90.mat[, c("CI.lower", "CI.upper"), drop = FALSE]
  rownames(att.bound) <- out$time

  lv <- unit.df$level
  status <- if (all(!is.finite(lv))) {
    "no sets"
  } else if (any(is.infinite(unit.df$CI.lower) | is.infinite(unit.df$CI.upper),
                 na.rm = TRUE)) {
    "partly infinite: too few controls for a finite set at some dates"
  } else {
    "ok"
  }
  if (any(unit.df$disconnected %in% TRUE)) {
    status <- paste0(status, "; ", sum(unit.df$disconnected %in% TRUE),
                     " set(s) are unions of intervals (hull reported; see ",
                     "fit$conformal$components)")
  }
  diag.unit <- NULL
  if (!is.null(res$diag)) {
    dg <- res$diag
    diag.unit <- data.frame(unit = dg$unit, time = dg$time,
                            event.time = dg$event.time, eff = dg$eff,
                            CI.lower = dg$CI.lower, CI.upper = dg$CI.upper,
                            p.value = dg$p.value, N = dg$N,
                            stringsAsFactors = FALSE)
  }
  conformal <- list(
    procedure = "Full conformal inference for target-date estimators (Liu and Xu 2026)",
    procedure.type = "full",
    method = "fe",
    rule = rule.name,
    score = score,
    level = suppressWarnings(min(lv[is.finite(lv)])),
    gamma = alpha,
    status = status,
    n.controls = res$n.controls,
    n.treated = Ntr,
    unit = unit.df,
    avg = res$avg,
    components = res$components,
    aggregate = "none",
    diag = if (diag > 0) list(m = diag, unit = diag.unit, att = diag.att,
                              components = res$diag.components) else NULL,
    extension = if (identical(score, "studentized")) {
      paste0("Studentized score: full conformal sets are reported for each ",
             "treated unit and date (fit$conformal$unit). The average over ",
             "the post-treatment dates has no set with this score, so ",
             "est.avg carries no interval.",
             if (Ntr == 1L) paste0(" With one treated unit, est.att and the ",
                                   "calendar-time ATT carry its sets at each ",
                                   "date.") else paste0(" Averages over ",
                                   "several treated units have no set here."))
    } else if (Ntr == 1L) {
      paste0("With one treated unit, est.att and the calendar-time ATT ",
             "carry its full conformal sets at each date, and est.avg the ",
             "set for its effect averaged over the post-treatment dates.")
    } else {
      paste0("Full conformal sets are reported for each treated unit and ",
             "date (fit$conformal$unit) and for each unit's average over ",
             "its dates (fit$conformal$avg). Averages over several treated ",
             "units have no full conformal set here, so est.att and est.avg ",
             "carry no interval.")
    }
  )
  if (!is.null(res$sigma.treated)) {
    conformal$sigma.treated <- res$sigma.treated
  }
  list(
    est.avg = est.avg,
    est.avg.unit = est.avg.unit,
    est.att = att.mat,
    est.att90 = att90.mat,
    att.bound = att.bound,
    est.eff.calendar = est.cal,
    est.eff.calendar.fit = est.cal.fit,
    est.beta = NULL,
    vartype = "conformal",
    conformal = conformal
  )
}
