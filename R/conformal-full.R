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
## roots by a sweep: no grid. It can be a union of intervals. Affinity is
## checked numerically at a third candidate; a rule that is not affine falls
## back to a grid of candidates (each refits every control).
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
    return(list(fun = rule, name = "user rule"))
  }
  if (!(is.character(rule) && length(rule) == 1L)) {
    stop("rule must be \"did\", \"sc\" or a function.", call. = FALSE)
  }
  switch(rule,
         did = list(fun = .fc_rule_did, name = "did"),
         sc = list(fun = .fc_make_rule_sc(sc.eta),
                   name = paste0("sc (ridge-simplex, eta = ", format(sc.eta),
                                 ")")),
         unit = list(fun = .fc_rule_unit, name = "unit"),
         time = list(fun = .fc_rule_time, name = "time"),
         stop("rule must be \"did\", \"sc\" or a function; got \"", rule,
              "\".", call. = FALSE))
}

## One call of the rule, checked: a single finite number or NA.
.fc_call <- function(fun, qp, bp, bt) {
  v <- tryCatch(fun(qp, bp, bt), error = function(e) NA_real_)
  v <- suppressWarnings(as.numeric(v))
  if (length(v) != 1L || !is.finite(v)) NA_real_ else v
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
.fc_coefs <- function(Ypre, ytg, i0, fun) {
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
    bt[pos] <- ytg[i0] + h
    p.m <- .fc_call(fun, Ypre[, j], bp, bt)
    a[j] <- ytg[j] - p.0
    b[j] <- (p.0 - p.h) / h
    ## affinity: the prediction at -h must continue the line
    if (!all(is.finite(c(p.0, p.h, p.m))) ||
        abs((p.m - p.0) - (p.0 - p.h)) >
        1e-8 * max(1, abs(p.0), abs(p.h), abs(p.m))) {
      affine <- FALSE
    }
  }
  list(a = a, b = b, affine = affine, pred = pred0,
       failed = !is.finite(pred0) || any(!is.finite(a)))
}

## Grid inversion for a rule that is not affine: p(c) at each candidate,
## each needing a refit of every control.
.fc_grid <- function(Ypre, ytg, i0, fun, gamma, grid, tau) {
  N <- ncol(Ypre)
  ctrl <- setdiff(seq_len(N), i0)
  pval <- function(c0) {
    yc <- ytg
    yc[i0] <- ytg[i0] - c0
    s <- numeric(N)
    s[i0] <- abs(tau - c0)
    for (j in ctrl) {
      bg <- setdiff(seq_len(N), j)
      pj <- .fc_call(fun, Ypre[, j], Ypre[, bg, drop = FALSE], yc[bg])
      s[j] <- if (is.finite(pj)) abs(ytg[j] - pj) else Inf
    }
    mean(s >= s[i0])
  }
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
## each control j, predict both from the other controls.
.fc_jplus <- function(Ypre, ytg, i0, fun, gamma) {
  N <- ncol(Ypre)
  ctrl <- setdiff(seq_len(N), i0)
  J <- length(ctrl)
  M0 <- R <- rep(NA_real_, J)
  for (k in seq_len(J)) {
    j <- ctrl[k]
    bg <- setdiff(ctrl, j)
    if (!length(bg)) next
    bp <- Ypre[, bg, drop = FALSE]
    M0[k] <- .fc_call(fun, Ypre[, i0], bp, ytg[bg])
    mj <- .fc_call(fun, Ypre[, j], bp, ytg[bg])
    R[k] <- abs(ytg[j] - mj)
  }
  v <- .conformal_interval(M0, R, ytg[i0], gamma)
  list(CI.lower = unname(v["CI.lower"]), CI.upper = unname(v["CI.upper"]),
       p.value = unname(v["p.value"]), N = unname(v["N"]),
       q = unname(v["q"]), n.failed = unname(v["n.failed"]),
       level = if (v["q"] >= 1) 1 - (2 * v["q"] - 1) / v["N"] else 0)
}

## Default grid for a rule that is not affine: centred on the estimate and
## wide enough to hold every control's error at c = tau_hat several times.
.fc_default_grid <- function(tau, cf, n = 401L) {
  w <- max(abs(cf$a), na.rm = TRUE)
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

## One treated unit k: full conformal sets (or Jackknife+ intervals) at each
## target date, and the full conformal set for the average over the post-
## treatment dates. `targets` are row positions; `history` the rows used as
## pre-treatment dates for each target (the dates before adoption for post
## dates; the dates before s for a diagnostic date s).
.fc_one_unit <- function(Y, i, co, targets, history, fun, procedure, gamma,
                         grid) {
  rows <- vector("list", length(targets))
  coefs <- vector("list", length(targets))
  pools <- vector("list", length(targets))
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
      v <- .fc_jplus(Ypre, ytg, 1L, fun, gamma)
      pred <- .fc_call(fun, Ypre[, 1L], Ypre[, -1L, drop = FALSE], ytg[-1L])
      rows[[s]] <- c(base, list(eff = ytg[1L] - pred, CI.lower = v$CI.lower,
                                CI.upper = v$CI.upper, p.value = v$p.value,
                                N = v$N, q = v$q, level = unname(v$level),
                                disconnected = FALSE, method = "jackknife+",
                                components = list(NULL)))
      next
    }
    cf <- .fc_coefs(Ypre, ytg, 1L, fun)
    coefs[[s]] <- cf
    tau <- cf$a[1L]
    if (cf$failed) {
      rows[[s]] <- c(base, list(eff = tau, CI.lower = NA_real_,
                                CI.upper = NA_real_, p.value = NA_real_,
                                N = length(cols), level = NA_real_,
                                disconnected = NA, method = "failed",
                                components = list(NULL)))
      next
    }
    v <- if (cf$affine) {
      .fc_invert_affine(cf$a, cf$b, 1L, gamma)
    } else {
      g <- if (is.null(grid)) .fc_default_grid(tau, cf) else tau + grid
      .fc_grid(Ypre, ytg, 1L, fun, gamma, g, tau)
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
  list(rows = rows, coefs = coefs, pools = pools)
}

## Average over the post-treatment dates for one treated unit, from the
## per-date affine coefficients. Needs the same pool at every date and slopes
## that do not change across dates.
.fc_avg <- function(res, post.idx, gamma) {
  cf <- res$coefs[post.idx]
  pools <- res$pools[post.idx]
  no <- function(why) list(eff = NA_real_, CI.lower = NA_real_,
                           CI.upper = NA_real_, p.value = NA_real_,
                           N = NA_real_, level = NA_real_,
                           disconnected = NA, note = why,
                           components = NULL)
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
                      units = colnames(Y), times = rownames(Y)) {
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
  avg <- list()
  comps <- list()
  dcomps <- list()
  for (k in seq_along(tr)) {
    i <- tr[k]
    on <- ds$T.on[k]
    pre <- seq_len(on - 1L)
    post <- on:TT
    post <- post[!is.na(Y[post, i])]
    if (length(pre) < 1L) next
    res <- .fc_one_unit(Y, i, co, post, rep(list(pre), length(post)), fun,
                        procedure, gamma, grid)
    for (s in seq_along(post)) {
      r <- res$rows[[s]]
      comps[[length(comps) + 1L]] <- r$components[[1L]]
      tab[[length(tab) + 1L]] <- data.frame(
        unit = units[i], time = times[post[s]], event.time = post[s] - on + 1L,
        Y = Y[post[s], i], eff = r$eff, CI.lower = r$CI.lower,
        CI.upper = r$CI.upper, p.value = r$p.value, N = r$N,
        level = r$level, disconnected = r$disconnected,
        method = r$method, n.history = r$n.history,
        stringsAsFactors = FALSE)
    }
    if (procedure == "full") {
      a <- .fc_avg(res, seq_along(post), gamma)
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
                             fun, procedure, gamma, grid)
        for (s in seq_along(ds.rows)) {
          r <- dres$rows[[s]]
          dcomps[[length(dcomps) + 1L]] <- r$components[[1L]]
          dtab[[length(dtab) + 1L]] <- data.frame(
            unit = units[i], time = times[ds.rows[s]],
            event.time = ds.rows[s] - on + 1L, Y = Y[ds.rows[s], i],
            eff = r$eff, CI.lower = r$CI.lower, CI.upper = r$CI.upper,
            p.value = r$p.value, N = r$N, level = r$level,
            disconnected = r$disconnected, method = r$method,
            n.history = r$n.history, stringsAsFactors = FALSE)
        }
      }
    }
  }
  bind <- function(x) if (length(x)) do.call(rbind, x) else NULL
  list(unit = bind(tab), avg = bind(avg), diag = bind(dtab),
       components = comps, diag.components = dcomps,
       n.treated = length(tr), n.controls = length(co))
}


## ---- user-facing function ------------------------------------------------------

#' Conformal inference for target-date estimators
#'
#' Full conformal confidence sets, or Jackknife+ intervals, for the effect on
#' each treated unit at each post-treatment date, for any estimator that
#' predicts a treated unit's untreated outcome at a date from pre-treatment
#' outcomes and the controls' outcomes at that date. Synthetic control and the
#' difference-in-differences (two-way fixed effects) imputation with one
#' treated unit are such target-date estimators. The procedures are those of
#' Liu and Xu (2026).
#'
#' @param formula A formula \code{Y ~ D} with the outcome and the 0/1
#'   treatment. Covariates are not supported.
#' @param data A data frame in long format.
#' @param index A character vector of length two: the unit and time
#'   variables.
#' @param rule The estimator: \code{"did"} (difference in differences: the
#'   held-out unit's pre-treatment mean plus the training units' average
#'   change from their pre-treatment means to the target date), \code{"sc"}
#'   (synthetic control with simplex weights fitted on the pre-treatment
#'   outcomes, with a ridge term \code{sc.eta}; needs the \pkg{quadprog}
#'   package), or a function \code{function(query.pre, background.pre,
#'   background.target)} that returns the held-out unit's predicted outcome
#'   at the target date. \code{query.pre} holds the held-out unit's outcomes
#'   at the pre-treatment dates, \code{background.pre} the training units'
#'   outcomes at those dates (a matrix with one column per unit), and
#'   \code{background.target} the training units' outcomes at the target
#'   date. The function never sees the held-out unit's outcome at the target
#'   date.
#' @param procedure \code{"full"} (full conformal inference, the default) or
#'   \code{"jackknife+"}.
#' @param alpha Miscoverage level, \code{0 < alpha < 1}.
#' @param diag Number of pre-treatment dates at which to run the procedure as
#'   a placebo, each with the dates before it as its history (0 for none).
#' @param grid For a rule whose prediction is not affine in the completed
#'   outcome: candidate effects, as offsets from the estimate, at which the
#'   full conformal p-value is evaluated. By default 401 points spanning four
#'   times the largest control error on each side. Ignored for affine rules,
#'   whose sets are computed exactly.
#' @param sc.eta Ridge term of \code{rule = "sc"}.
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
#' completed outcome. A unit's score is its absolute prediction error; the
#' treated unit's score is \eqn{|\hat\tau - c|}. The p-value of \eqn{c} is the
#' share of units whose score is at least the treated unit's, and the set
#' keeps the candidates with p-value above \code{alpha}. It covers the effect
#' with probability at least \eqn{1 - \lfloor \alpha N \rfloor / N} if the
#' scores are exchangeable across units at the true effect: the treated unit
#' is no harder to predict than a control, and the rule and its tuning are
#' the same for every unit and fixed in advance. When the rule's prediction
#' is affine in the completed outcome, as for \code{"did"} and \code{"sc"},
#' the set is computed exactly; it can be a union of intervals
#' (\code{disconnected = TRUE}, the parts in \code{components}). Otherwise it
#' is evaluated on \code{grid}, and \code{lower.truncated} or
#' \code{upper.truncated} flag a set that reaches the end of the grid.
#'
#' For each treated unit, the set for the effect averaged over its
#' post-treatment dates (\code{avg}) uses the errors averaged over the dates.
#' It is computed when the rule is affine with the same slope at every date
#' (true for \code{"did"} and \code{"sc"}) and the pool is the same at every
#' date.
#'
#' \strong{Jackknife+.} The treated unit is held out together with each
#' control, and both are predicted from the other controls. The interval
#' covers the effect with probability at least \eqn{1 - (2q - 1)/N}, where
#' \eqn{q} is the largest integer with \eqn{(2q - 1)/N \le} \code{alpha}.
#'
#' The p-value reported for each date is for the null of no effect at that
#' date.
#'
#' @return An object of class \code{"conformal_td"}, a list with
#'   \item{unit}{One row per treated unit and post-treatment date:
#'     \code{eff} (estimate), \code{CI.lower}, \code{CI.upper} (the set's
#'     hull), \code{p.value}, \code{N}, \code{level} (the guaranteed
#'     coverage), \code{disconnected}, \code{method} (\code{"exact"},
#'     \code{"grid"} or \code{"jackknife+"}) and \code{n.history}.}
#'   \item{avg}{For \code{procedure = "full"}: one row per treated unit with
#'     the set for the effect averaged over its post-treatment dates.}
#'   \item{diag}{The pre-treatment placebo rows, if \code{diag > 0}.}
#'   \item{components}{A list with the components of each set in
#'     \code{unit}, in the same order.}
#'   \item{procedure, rule, alpha, call}{The settings.}
#'
#' @references Liu, Licheng and Yiqing Xu. 2026. "Conformal Inference for
#'   Panel Counterfactual Estimators." Working paper.
#'
#' @seealso \code{\link{fect}} with \code{method = "fe"} and
#'   \code{vartype = "conformal"}, which uses \code{rule = "did"}.
#'
#' @examples
#' data(sim_gsynth)
#' d1 <- subset(sim_gsynth, !(id %in% 102:105))
#' out <- conformal_td(Y ~ D, data = d1, index = c("id", "time"),
#'                     rule = "did", alpha = 0.1)
#' out
#' head(out$unit)
#'
#' \donttest{
#' if (requireNamespace("quadprog", quietly = TRUE)) {
#'   out.sc <- conformal_td(Y ~ D, data = d1, index = c("id", "time"),
#'                          rule = "sc", alpha = 0.1)
#'   out.sc$avg
#' }
#' ## a user rule: the mean of the training units at the target date, shifted
#' ## by the held-out unit's pre-treatment gap to them
#' my_rule <- function(query.pre, background.pre, background.target) {
#'   mean(background.target) + mean(query.pre - rowMeans(background.pre))
#' }
#' out.my <- conformal_td(Y ~ D, data = d1, index = c("id", "time"),
#'                        rule = my_rule)
#' }
#' @export
conformal_td <- function(formula, data, index, rule = "did",
                         procedure = c("full", "jackknife+"), alpha = 0.1,
                         diag = 0, grid = NULL, sc.eta = 1e-5) {
  cl <- match.call()
  procedure <- match.arg(procedure)
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
                   grid = grid, units = w$units, times = w$times)
  if (!is.null(res$unit) && procedure == "full" &&
      any(res$unit$disconnected %in% TRUE)) {
    message("Some full conformal sets are unions of intervals; CI.lower and ",
            "CI.upper give their hull. See $components.")
  }
  structure(c(res, list(procedure = procedure, rule = ru$name,
                        alpha = alpha, call = cl)),
            class = "conformal_td")
}

#' @export
print.conformal_td <- function(x, digits = 3, ...) {
  cat("Call:\n")
  print(x$call)
  lab <- if (x$procedure == "full") "Full conformal inference" else
    "Jackknife+"
  cat("\n", lab, " (Liu and Xu 2026), rule: ", x$rule, ", alpha = ",
      x$alpha, "\n", sep = "")
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
  if (!is.null(x$avg) && nrow(x$avg)) {
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

#' @export
plot.conformal_td <- function(x, ...) {
  if (is.null(x$unit) || !nrow(x$unit)) {
    stop("Nothing to plot.", call. = FALSE)
  }
  u <- x$unit
  if (!is.null(x$diag)) u <- rbind(u, x$diag)
  ggplot2::ggplot(u, ggplot2::aes(x = .data$event.time, y = .data$eff)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_vline(xintercept = 0.5, colour = "grey80") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$CI.lower,
                                          ymax = .data$CI.upper),
                             size = 0.2) +
    ggplot2::facet_wrap(~ unit) +
    ggplot2::labs(x = "Periods since treatment", y = "Effect") +
    ggplot2::theme_bw()
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
                              norm.para = NULL, diag = 0) {
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
  res <- .fc_panel(Yr, Dm, rule, "full", alpha, diag = diag,
                   units = seq_len(N), times = seq_len(TT))
  res90 <- .fc_panel(Yr, Dm, rule, "full", 2 * alpha, diag = diag,
                     units = seq_len(N), times = seq_len(TT))
  u <- res$unit
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
            " cell(s) (other treated units or not-yet-treated cells enter ",
            "fect's fit). Both are in fit$conformal$unit (eff, eff.fit).")
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
    extension = if (Ntr == 1L) {
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
