## Training floor for the rolling cross-validation (issue #167).
##
## Every held-out unit keeps at least floor = max(min.T0, 2 * (r.max + 1))
## observed periods before its buffer and held-out block, in r.cv.rolling()
## and in the shared mask builder .build_cv_mask_rolling() (fect(CV = TRUE),
## gsynth, the binary CV, fect_mspe). No unit is dropped by the inner fits,
## so no held-out cell goes unscored. When no unit can meet the floor it is
## lowered to what the data allow (never below min.T0) with a message. The
## fold standard errors the rule used come back as CV.out.se, and the
## default cv.rule is "min" (the anchor study in
## Research_Hub/Packages/fect/triage/2026-09-29-cv-anchor-study/).

skip_on_cran()

.sg <- function() {
  e <- new.env()
  data("sim_gsynth", package = "fect", envir = e)
  e$sim_gsynth
}

## sim_gsynth cut to 14 periods, with the five treated units treated at
## periods 13 and 14 (they keep 12 pre-treatment periods; controls have 14).
.short_panel <- function() {
  d <- .sg()
  treated <- unique(d$id[d$D == 1])
  d <- d[d$time <= 14, ]
  d$D <- as.integer(d$id %in% treated & d$time >= 13)
  d
}

## --------------------------------------------------------------------------
## (b) The study's "floor" fold loop, ported from study.R (run_cv() under
## rule = "floor"): same eligibility, per-fold unit sampling, anchors,
## masking, inner fect() call, scoring and aggregation. study.R fits the
## inner model with min.T0 = FLOOR; the package keeps the user's min.T0
## (5). Under the floor no unit has fewer than FLOOR training periods, so
## the two agree cell for cell.
## --------------------------------------------------------------------------
.study_floor_cv <- function(data, r.max = 3L, k = 20L, cv.prop = 0.1,
                            cv.nobs = 3L, cv.buffer = 1L, seed = 1L,
                            inner.min.T0 = 2L * (r.max + 1L)) {
  FLOOR <- 2L * (r.max + 1L)
  R_GRID <- 0L:r.max
  onset <- tapply(data$time[data$D >= 1], as.character(data$id[data$D >= 1]), min)
  elig <- by(data, data$id, function(d) {
    u <- as.character(d$id[1L]); t_all <- sort(d$time)
    if (u %in% names(onset)) t_all[t_all < onset[[u]]] else t_all
  })
  elig <- elig[vapply(elig, length, integer(1)) >= FLOOR + cv.buffer + cv.nobs]
  n_el <- length(elig)
  n_samp <- max(1L, as.integer(round(cv.prop * n_el)))
  keys <- paste(data$id, data$time, sep = "_")
  fold_mspe <- matrix(NA_real_, length(R_GRID), k)
  fold_n <- matrix(0L, length(R_GRID), k)
  fold_expected <- integer(k)
  for (fold in seq_len(k)) {
    set.seed(seed + fold)
    units <- names(elig)[sort(sample.int(n_el, n_samp))]
    rec <- lapply(units, function(u) {
      obs <- elig[[u]]; n <- length(obs)
      lo <- FLOOR + cv.buffer + 1L; hi <- n - cv.nobs + 1L
      if (hi < lo) return(NULL)
      valid <- seq.int(lo, hi)
      a <- valid[sample.int(length(valid), 1L)]
      list(unit = u,
           hold = obs[a:(a + cv.nobs - 1L)],
           buf  = if (cv.buffer > 0L) obs[max(1L, a - cv.buffer):(a - 1L)] else integer(0),
           drop = if ((a + cv.nobs) <= n) obs[(a + cv.nobs):n] else integer(0))
    })
    rec <- rec[!vapply(rec, is.null, logical(1))]
    if (length(rec) == 0L) next
    scored <- unlist(lapply(rec, function(r) paste(r$unit, r$hold, sep = "_")))
    masked <- c(scored,
                unlist(lapply(rec, function(r) if (length(r$buf)) paste(r$unit, r$buf, sep = "_") else character(0))),
                unlist(lapply(rec, function(r) if (length(r$drop)) paste(r$unit, r$drop, sep = "_") else character(0))))
    dm <- data
    dm$Y[keys %in% masked] <- NA_real_
    rows <- which(keys %in% scored)
    Yobs <- data$Y[rows]
    fold_expected[fold] <- length(rows)
    for (idx in seq_along(R_GRID)) {
      fit <- tryCatch(suppressWarnings(suppressMessages(fect(
        Y ~ D + X1 + X2, data = dm, index = c("id", "time"),
        method = "gsynth", time.component.from = "nevertreated",
        force = "two-way", CV = FALSE, r = R_GRID[idx], min.T0 = inner.min.T0,
        se = FALSE, na.rm = FALSE, parallel = FALSE))),
        error = function(e) NULL)
      if (is.null(fit) || is.null(fit$Y.ct.full)) next
      ui <- match(data$id[rows], fit$id)
      ti <- match(data$time[rows], fit$rawtime)
      ok <- !is.na(ui) & !is.na(ti)
      pred <- rep(NA_real_, length(rows))
      if (any(ok)) {
        pred[ok] <- fect:::.fect_heldout_pred(fit, ti[ok], ui[ok], data[rows[ok], , drop = FALSE])
      }
      e2 <- (Yobs - pred)^2
      e2 <- e2[is.finite(e2)]
      if (length(e2) == 0L) next
      fold_mspe[idx, fold] <- mean(e2)
      fold_n[idx, fold] <- length(e2)
    }
  }
  list(fold_mspe = fold_mspe, fold_n = fold_n, n_expected = sum(fold_expected),
       mspe = rowMeans(fold_mspe, na.rm = TRUE), n_eligible = n_el)
}

test_that("(a) r.cv.rolling() on sim_gsynth at the defaults picks r = 2 with every held-out cell scored", {
  d <- .sg()
  res <- suppressMessages(r.cv.rolling(
    Y ~ D + X1 + X2, data = d, index = c("id", "time"),
    method = "gsynth", r.max = 3, force = "two-way", seed = 1,
    verbose = FALSE, parallel = FALSE))
  ## floor = max(5, 2 * (3 + 1)) = 8 training periods per held-out unit
  expect_equal(res$floor, 8L)
  ## before the floor (dev 0e100bd) this picked r = 0 with 12 of 300 cells unscored
  expect_equal(res$r.cv, 2L)
  expect_true("n_unscored" %in% names(res$mspe))
  expect_true(all(res$mspe$n_unscored == 0L))
  expect_true(all(res$mspe$n_holdout == res$mspe$n_holdout[1L]))
  expect_equal(res$mspe$n_holdout[1L], 20L * 5L * 3L)   # k * units per fold * cv.nobs
  expect_true(all(res$mspe$n_folds_used == 20L))
  expect_true(all(is.finite(res$mspe$se)))
})

test_that("(b) r.cv.rolling() reproduces the study's floor loop fold for fold", {
  d <- .sg()
  ref <- suppressMessages(r.cv.rolling(
    Y ~ D + X1 + X2, data = d, index = c("id", "time"),
    method = "gsynth", r.max = 3, force = "two-way", seed = 1,
    verbose = FALSE, parallel = FALSE))
  mine <- .study_floor_cv(d, r.max = 3L, seed = 1L)
  expect_equal(unname(ref$mspe.per.fold), mine$fold_mspe, tolerance = 1e-8)
  expect_equal(unname(ref$mspe$mspe), unname(mine$mspe), tolerance = 1e-8)
  expect_equal(sum(ref$mspe$n_holdout[1L]), mine$n_expected)
  expect_equal(unname(rowSums(mine$fold_n)), rep(mine$n_expected, 4L))
})

test_that("(c) the mask builder keeps the floor and fect(CV = TRUE) returns CV.out.se", {
  d <- .sg()
  ids <- sort(unique(d$id)); tt <- sort(unique(d$time))
  TT <- length(tt); N <- length(ids)
  D <- matrix(0, TT, N)
  D[cbind(match(d$time, tt), match(d$id, ids))] <- d$D
  II <- matrix(1L, TT, N)
  II[D == 1] <- 0L                           # as fect builds II: untreated observed cells
  folds <- fect:::.build_cv_mask_rolling(
    II = II, D = D, k = 20L, cv.nobs = 3L, cv.buffer = 1L, cv.prop = 0.1,
    min.T0 = 5L, r.max = 3L, seed = 1L)
  expect_length(folds, 20L)
  expect_equal(attr(folds, "floor"), 8L)
  n_pre <- colSums(II)                       # observed periods before treatment
  for (f in folds) {
    II.cv <- II
    II.cv[f$cv.id] <- 0L
    masked_units <- unique((f$cv.id - 1L) %/% TT + 1L)
    train <- colSums(II.cv)[masked_units]
    expect_true(all(train >= 8L))
    expect_true(all(train <= n_pre[masked_units] - 1L - 3L))
    expect_equal(length(f$est.id), 3L * length(masked_units))
  }

  ## fect(CV = TRUE, method = "gsynth"): the nevertreated path
  msgs <- character(0)
  fit <- withCallingHandlers(
    suppressWarnings(fect(Y ~ D + X1 + X2, data = d, index = c("id", "time"),
                          method = "gsynth", force = "two-way", CV = TRUE,
                          r = c(0, 3), se = FALSE, parallel = FALSE, seed = 1)),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_false(any(grepl("dropped|Removed automatically", msgs)))
  expect_equal(fit$N, N)
  expect_true(!is.null(fit$CV.out.se))
  expect_equal(nrow(fit$CV.out.se), nrow(fit$CV.out))
  expect_equal(fit$CV.out.se[, "r"], fit$CV.out[, "r"])
  expect_true(all(c("MSPE", "WMSPE", "GMSPE", "WGMSPE", "MAD", "Moment", "GMoment") %in% colnames(fit$CV.out.se)))
  expect_true(all(is.finite(fit$CV.out.se[, "MSPE"])))
  expect_true(all(fit$CV.out.se[, "MSPE"] > 0))
  expect_true(all(is.finite(fit$CV.out[, "MSPE"])))

  ## fect_cv(): the ife and mc tables carry their own SE tables
  e <- new.env(); data("simdata", package = "fect", envir = e)
  fit_ife <- suppressMessages(suppressWarnings(fect(
    Y ~ D + X1 + X2, data = e$simdata, index = c("id", "time"),
    method = "ife", CV = TRUE, r = c(0, 2), k = 5, se = FALSE,
    parallel = FALSE, seed = 1)))
  expect_true(!is.null(fit_ife$CV.out.ife.se))
  expect_equal(nrow(fit_ife$CV.out.ife.se), nrow(fit_ife$CV.out.ife))
  expect_equal(fit_ife$CV.out.ife.se[, "r"], fit_ife$CV.out.ife[, "r"])
  expect_true(all(is.finite(fit_ife$CV.out.ife.se[, "MSPE"])))
  ## and, for the method fitted, as CV.out and CV.out.se (exact names: the
  ## .se tables make a partial match of fit$CV.out ambiguous)
  expect_identical(fit_ife$CV.out, fit_ife$CV.out.ife)
  expect_identical(fit_ife$CV.out.se, fit_ife$CV.out.ife.se)
  fit_mc <- suppressMessages(suppressWarnings(fect(
    Y ~ D + X1 + X2, data = e$simdata, index = c("id", "time"),
    method = "mc", CV = TRUE, k = 5, se = FALSE, parallel = FALSE, seed = 1)))
  expect_identical(fit_mc$CV.out, fit_mc$CV.out.mc)
  expect_identical(fit_mc$CV.out.se, fit_mc$CV.out.mc.se)
  expect_equal(nrow(fit_mc$CV.out.se), nrow(fit_mc$CV.out))
  expect_true("lambda.norm" %in% colnames(fit_mc$CV.out.se))
})

test_that("(d) a panel too short for the floor lowers it with a message", {
  d <- .short_panel()
  ids <- sort(unique(d$id)); tt <- sort(unique(d$time))
  TT <- length(tt); N <- length(ids)
  D <- matrix(0, TT, N)
  D[cbind(match(d$time, tt), match(d$id, ids))] <- d$D
  II <- matrix(1L, TT, N)
  II[D == 1] <- 0L
  ## r.max = 5: floor 12, a unit needs 12 + 1 + 3 = 16 periods; the
  ## controls have 14, so the floor drops to 14 - 1 - 3 = 10, which can
  ## judge ranks up to 10 / 2 - 1 = 4.
  expect_message(
    folds <- fect:::.build_cv_mask_rolling(
      II = II, D = D, k = 3L, cv.nobs = 3L, cv.buffer = 1L, cv.prop = 0.1,
      min.T0 = 5L, r.max = 5L, seed = 1L),
    regexp = "floor of 10 .*ranks up to about 4")
  expect_equal(attr(folds, "floor"), 10L)
  ## the five treated units (12 pre-periods) cannot meet 10 + 1 + 3
  expect_equal(attr(folds, "n.eligible"), 45L)
  for (f in folds) {
    II.cv <- II
    II.cv[f$cv.id] <- 0L
    masked_units <- unique((f$cv.id - 1L) %/% TT + 1L)
    expect_true(all(colSums(II.cv)[masked_units] >= 10L))
  }
  ## the floor is met at r.max = 3 (floor 8, 12 periods needed): no message
  expect_silent(
    folds3 <- fect:::.build_cv_mask_rolling(
      II = II, D = D, k = 3L, cv.nobs = 3L, cv.buffer = 1L, cv.prop = 0.1,
      min.T0 = 5L, r.max = 3L, seed = 1L))
  expect_equal(attr(folds3, "floor"), 8L)
  ## never below min.T0: with min.T0 = 12 and 14 periods no unit is eligible
  expect_error(
    fect:::.build_cv_mask_rolling(
      II = II, D = D, k = 3L, cv.nobs = 3L, cv.buffer = 1L, cv.prop = 0.1,
      min.T0 = 12L, r.max = 5L, seed = 1L),
    regexp = "no eligible units")

  ## the same message from r.cv.rolling()
  expect_message(
    res <- r.cv.rolling(Y ~ D + X1 + X2, data = d, index = c("id", "time"),
                        method = "gsynth", r.max = 5, k = 2, force = "two-way",
                        seed = 1, verbose = FALSE, parallel = FALSE),
    regexp = "floor of 10 .*ranks up to about 4")
  expect_equal(res$floor, 10L)
  expect_equal(res$n.eligible, 45L)
  expect_true(all(res$mspe$n_unscored == 0L))
})

## sim_gsynth cut to 14 periods, with the five treated units treated from
## period 12 (they keep 11 pre-treatment periods; controls have 14).
.short_panel_11 <- function() {
  d <- .sg()
  treated <- unique(d$id[d$D == 1])
  d <- d[d$time <= 14, ]
  d$D <- as.integer(d$id %in% treated & d$time >= 12)
  d
}

test_that("(f) a unit one period short of floor + cv.buffer + cv.nobs is not eligible", {
  ## r.max = 3: floor 8, so a unit needs 8 + 1 + 3 = 12 periods. The
  ## treated units have 11, one short, so only the 45 controls are
  ## eligible. An eligibility rule of min.T0 + cv.buffer + cv.nobs (= 9)
  ## would admit them, in the mask builder and in r.cv.rolling() alike.
  d <- .short_panel_11()
  ids <- sort(unique(d$id)); tt <- sort(unique(d$time))
  TT <- length(tt); N <- length(ids)
  D <- matrix(0, TT, N)
  D[cbind(match(d$time, tt), match(d$id, ids))] <- d$D
  II <- matrix(1L, TT, N)
  II[D == 1] <- 0L
  treated_cols <- which(colSums(D) > 0)
  expect_length(treated_cols, 5L)
  expect_true(all(colSums(II)[treated_cols] == 11L))
  folds <- fect:::.build_cv_mask_rolling(
    II = II, D = D, k = 3L, cv.nobs = 3L, cv.buffer = 1L, cv.prop = 1,
    min.T0 = 5L, r.max = 3L, seed = 1L)
  expect_equal(attr(folds, "floor"), 8L)
  expect_equal(attr(folds, "n.eligible"), 45L)
  for (f in folds) {
    masked_units <- unique((f$cv.id - 1L) %/% TT + 1L)
    expect_false(any(treated_cols %in% masked_units))
    II.cv <- II
    II.cv[f$cv.id] <- 0L
    expect_true(all(colSums(II.cv)[masked_units] >= 8L))
  }
  res <- suppressMessages(r.cv.rolling(
    Y ~ D + X1 + X2, data = d, index = c("id", "time"),
    method = "gsynth", r.max = 3, k = 2, force = "two-way", seed = 1,
    verbose = FALSE, parallel = FALSE))
  expect_equal(res$floor, 8L)
  expect_equal(res$n.eligible, 45L)
  expect_true(all(res$mspe$n_unscored == 0L))
})

test_that("(e) the default cv.rule is \"min\" and \"1se\" still works", {
  expect_identical(formals(fect)$cv.rule, "min")
  expect_identical(eval(formals(r.cv.rolling)$cv.rule)[1L], "min")
  expect_identical(fect:::.fect_validate_cv_rule(NULL), "min")
  d <- .sg()
  fit_min <- suppressMessages(suppressWarnings(fect(
    Y ~ D + X1 + X2, data = d, index = c("id", "time"),
    method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 2), k = 5,
    se = FALSE, parallel = FALSE, seed = 1)))
  fit_1se <- suppressMessages(suppressWarnings(fect(
    Y ~ D + X1 + X2, data = d, index = c("id", "time"),
    method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 2), k = 5,
    cv.rule = "1se", se = FALSE, parallel = FALSE, seed = 1)))
  ## same folds, same table; the rules read it differently
  expect_equal(fit_min$CV.out[, "MSPE"], fit_1se$CV.out[, "MSPE"], tolerance = 1e-10)
  m <- fit_min$CV.out[, "MSPE"]
  expect_equal(as.integer(fit_min$r.cv), as.integer(fit_min$CV.out[which.min(m), "r"]))
  i_min <- which.min(m)
  expect_equal(as.integer(fit_1se$r.cv),
               as.integer(fit_1se$CV.out[min(which(m <= m[i_min] + fit_1se$CV.out.se[i_min, "MSPE"])), "r"]))
})
