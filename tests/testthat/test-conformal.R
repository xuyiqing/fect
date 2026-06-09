## Tests for the Family-A conformal core (R/conformal.R) and the full
## leave-one-control-out calibration (conformal_calibrate, via impute_Y0).
## Every Family-A interval is a level statistic |m - tau| / scale: closed-form and
## never empty. Coverage is validated through the calibration itself; the fect()
## output-slot wiring (vartype = "conformal") gets its print/plot tests once that
## lands.

test_that("center: mean and median", {
  set.seed(1); e <- rnorm(8)
  expect_equal(fect:::.conformal_center(e, "mean"), mean(e))
  expect_equal(fect:::.conformal_center(e, "median"), stats::median(e))
  expect_true(is.na(fect:::.conformal_center(numeric(0), "mean")))
})

test_that("scale: none / sd / rmspe / mad / diff (+ model-se is external)", {
  set.seed(2); e <- rnorm(20)
  expect_equal(fect:::.conformal_scale(e, "none"), 1)
  expect_equal(fect:::.conformal_scale(e, "sd"), stats::sd(e))
  expect_equal(fect:::.conformal_scale(e, "rmspe"), sqrt(mean(e^2)))
  expect_equal(fect:::.conformal_scale(e, "mad"), stats::mad(e))
  expect_equal(fect:::.conformal_scale(e, "diff"), stats::sd(diff(e)))
  expect_true(is.na(fect:::.conformal_scale(e, "model-se")))   # supplied externally
})

test_that("conformal quantile reduces to the order statistic and guards resolution", {
  set.seed(3); s <- abs(rnorm(20))
  for (a in c(0.05, 0.10, 0.20)) {
    k <- ceiling((1 - a) * (length(s) + 1))
    expected <- if (k > length(s)) Inf else sort(s)[k]
    expect_equal(fect:::.conformal_quantile(s, a), expected)
  }
  ## alpha below 1/(n+1) -> Inf (cannot reach the level)
  expect_true(is.infinite(fect:::.conformal_quantile(s, 1 / 30)))  # 1/30 < 1/21
})

test_that("conformal p-value lies in (0,1] and respects weights", {
  s.co <- rnorm(20)
  p <- fect:::.conformal_pval(0, s.co); expect_gt(p, 0); expect_lte(p, 1)
  expect_equal(fect:::.conformal_pval(100, s.co), 1 / (length(s.co) + 1))
  w <- rep(1, 20); w[which.max(s.co)] <- 0
  expect_true(is.finite(fect:::.conformal_pval(0, s.co, w.co = w)))
})

test_that("level interval brackets the center, scales with the treated scale, never empty", {
  set.seed(4)
  m.co <- rnorm(40); sc.co <- rep(1, 40)
  ci1 <- fect:::.conformal_ci_level(2, 1, m.co, sc.co, alpha = 0.10)
  expect_true(all(is.finite(ci1)))
  expect_lt(ci1[1], 2); expect_gt(ci1[2], 2)            # brackets the center
  ci2 <- fect:::.conformal_ci_level(2, 2, m.co, sc.co, alpha = 0.10)
  expect_equal(diff(ci2), 2 * diff(ci1))                # half-width = scale.tr * Q
  ## below the resolution floor -> unbounded, never empty
  ci3 <- fect:::.conformal_ci_level(2, 1, rnorm(15), rep(1, 15), alpha = 1 / 30)
  expect_true(all(is.infinite(ci3)))
})

test_that("effective sample size is sensible", {
  expect_equal(fect:::.conformal_neff(rep(1, 10)), 10)
  expect_lt(fect:::.conformal_neff(c(rep(1, 9), 100)), 10)
})

## ---- full leave-one-control-out calibration (coverage) ----------------------
## Block DGP, known factor structure, zero true effect. Regression guards: never
## empty, never unbounded (Ncal above the floor), coverage near nominal. Precise
## calibration is the coverage study's job.

.conf_dgp <- function(seed, N = 25, T = 20, T0 = 15, r = 2) {
  set.seed(seed)
  Fm <- matrix(rnorm(T * r), T, r); L <- matrix(rnorm(N * r), N, r)
  D  <- matrix(0, T, N); D[(T0 + 1):T, 1] <- 1
  list(mu = Fm %*% t(L), D = D, N = N, T = T)
}

.conf_cover <- function(seed, scale, reps = 80, alpha = 0.10) {
  g <- .conf_dgp(seed); cov <- logical(reps); bad <- 0L
  for (b in seq_len(reps)) {
    Y <- g$mu + matrix(rnorm(g$N * g$T), g$T, g$N)
    dat <- data.frame(id = rep(1:g$N, each = g$T), time = rep(1:g$T, g$N),
                      Y = as.vector(Y), D = as.vector(g$D))
    fit <- suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
              method = "gsynth", force = 3, CV = FALSE, r = 2, se = FALSE,
              parallel = FALSE))
    cc <- conformal_calibrate(Y = Y, D = g$D, I = fit$I, II = fit$II, T.on = fit$T.on,
              r.cv = fit$r.cv, eff = fit$eff, method = "gsynth",
              scale = scale, alpha = alpha)
    if (cc$status != "ok") bad <- bad + 1L
    cov[b] <- (0 >= cc$ci[1]) && (0 <= cc$ci[2])
  }
  list(cov = mean(cov), bad = bad)
}

test_that("meanabs (scale = none) covers near nominal, never empty/unbounded", {
  skip_on_cran()
  ## scale = none is stable (~0.90-0.92 at 200 reps); loose bound guards gross
  ## miscalibration without flaking on Monte-Carlo noise. Precise calibration is
  ## the Phase 5 simulation study, not this unit test.
  r <- .conf_cover(20260608, "none")
  expect_equal(r$bad, 0L)
  expect_gt(r$cov, 0.80); expect_lte(r$cov, 1.0)
})

test_that("studentized-mean (scale = sd) runs, never empty/unbounded, ballpark coverage", {
  skip_on_cran()
  ## scale = sd is calibrated on average (~0.90) but has higher coverage variance
  ## than none, because its width is random via the treated pre-period sd. The
  ## strict assertion is bad == 0 (Family A can never empty); coverage is a wide
  ## gross-guard here, with the precise check deferred to the Phase 5 study.
  r <- .conf_cover(424242, "sd")
  expect_equal(r$bad, 0L)
  expect_gt(r$cov, 0.75); expect_lte(r$cov, 1.0)
})

## ---- fect(vartype = "conformal") end-to-end integration ---------------------

test_that("fect(vartype = 'conformal') returns a complete, printable object", {
  skip_on_cran()
  g <- .conf_dgp(1)
  Y <- g$mu + matrix(rnorm(g$N * g$T), g$T, g$N) + 1.5 * g$D   # true ATT = 1.5
  dat <- data.frame(id = rep(1:g$N, each = g$T), time = rep(1:g$T, g$N),
                    Y = as.vector(Y), D = as.vector(g$D))
  f <- suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
            method = "gsynth", force = 3, CV = FALSE, r = 2, se = TRUE,
            vartype = "conformal"))
  expect_equal(f$vartype, "conformal")
  expect_identical(f$conformal$status, "ok")
  ## est.avg: finite, ordered CI with the right columns
  expect_true(all(c("ATT.avg", "CI.lower", "CI.upper", "p.value") %in% colnames(f$est.avg)))
  expect_lt(f$est.avg[1, "CI.lower"], f$est.avg[1, "CI.upper"])
  ## est.att: event-time rownames, post-period CIs finite
  expect_equal(rownames(f$est.att), as.character(f$time))
  post <- f$est.att[as.numeric(rownames(f$est.att)) >= 0, , drop = FALSE]
  expect_true(all(is.finite(post[, c("CI.lower", "CI.upper")])))
  expect_silent(invisible(capture.output(print(f))))
})

test_that("conformal.scale threads through fect() and changes the interval", {
  skip_on_cran()
  g <- .conf_dgp(7)
  Y <- g$mu + matrix(rnorm(g$N * g$T), g$T, g$N) + 1.5 * g$D
  dat <- data.frame(id = rep(1:g$N, each = g$T), time = rep(1:g$T, g$N),
                    Y = as.vector(Y), D = as.vector(g$D))
  run <- function(scale) suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
            method = "gsynth", force = 3, CV = FALSE, r = 2, se = TRUE,
            vartype = "conformal", conformal.scale = scale))
  w_none <- diff(run("none")$est.avg[1, c("CI.lower", "CI.upper")])
  w_diff <- diff(run("diff")$est.avg[1, c("CI.lower", "CI.upper")])
  expect_true(is.finite(w_none) && is.finite(w_diff))
  expect_false(isTRUE(all.equal(w_none, w_diff)))   # the scale actually does something
  ## not-yet-implemented knobs error clearly
  expect_error(run("model-se"), "not yet implemented")
})

## ---- multi-treated weighting (cell / unit / precision) ----------------------

test_that("conformal weights run for multiple treated units; cell matches att.avg", {
  skip_on_cran()
  set.seed(3); N <- 30; T <- 20; T0 <- 15; r <- 2
  Fm <- matrix(rnorm(T * r), T, r); L <- matrix(rnorm(N * r), N, r)
  D <- matrix(0, T, N); D[(T0 + 1):T, 1:4] <- 1
  Y <- Fm %*% t(L) + matrix(rnorm(N * T), T, N) + 1.5 * D
  dat <- data.frame(id = rep(1:N, each = T), time = rep(1:T, N),
                    Y = as.vector(Y), D = as.vector(D))
  pf <- suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
            method = "gsynth", force = 3, CV = FALSE, r = 2, se = FALSE))
  for (w in c("cell", "unit", "precision")) {
    f <- suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
              method = "gsynth", force = 3, CV = FALSE, r = 2, se = TRUE,
              vartype = "conformal", conformal.weight = w))
    expect_identical(f$conformal$status, "ok")
    expect_true(all(is.finite(f$est.avg[1, c("CI.lower", "CI.upper")])))
    if (w == "cell") {
      expect_equal(unname(f$est.avg[1, "ATT.avg"]), pf$att.avg, tolerance = 1e-6)
    }
  }
})

test_that("precision weight down-weights a noisy treated unit", {
  skip_on_cran()
  set.seed(5); N <- 30; T <- 20; T0 <- 15; r <- 2
  Fm <- matrix(rnorm(T * r), T, r); L <- matrix(rnorm(N * r), N, r)
  D <- matrix(0, T, N); D[(T0 + 1):T, 1:4] <- 1
  E <- matrix(rnorm(N * T), T, N); E[, 1] <- E[, 1] * 4   # treated unit 1 noisy
  Y <- Fm %*% t(L) + E + 1.5 * D
  dat <- data.frame(id = rep(1:N, each = T), time = rep(1:T, N),
                    Y = as.vector(Y), D = as.vector(D))
  att <- function(w) suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
            method = "gsynth", force = 3, CV = FALSE, r = 2, se = TRUE,
            vartype = "conformal", conformal.weight = w))$est.avg[1, "ATT.avg"]
  expect_false(isTRUE(all.equal(att("unit"), att("precision"))))
})

## ---- bands: pointwise / simultaneous / pooled cutoff ------------------------

test_that("simultaneous band is wider than pointwise; est.att.sim always present", {
  skip_on_cran()
  g <- .conf_dgp(1)
  Y <- g$mu + matrix(rnorm(g$N * g$T), g$T, g$N) + 1.5 * g$D
  dat <- data.frame(id = rep(1:g$N, each = g$T), time = rep(1:g$T, g$N),
                    Y = as.vector(Y), D = as.vector(g$D))
  run <- function(...) suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
            method = "gsynth", force = 3, CV = FALSE, r = 2, se = TRUE,
            vartype = "conformal", ...))
  fp <- run(); fs <- run(conformal.band = "simultaneous")
  post_w <- function(m) {
    pm <- m[as.numeric(rownames(m)) >= 0, c("CI.lower", "CI.upper"), drop = FALSE]
    mean(pm[, 2] - pm[, 1])
  }
  expect_gte(post_w(fs$est.att), post_w(fp$est.att))      # uniform band is wider
  expect_false(is.null(fp$est.att.sim))                   # always computed
  expect_gte(post_w(fp$est.att.sim), post_w(fp$est.att))
  ## pooled cutoff runs and is finite
  fpool <- run(conformal.cutoff = "pooled")
  expect_identical(fpool$conformal$status, "ok")
  expect_true(all(is.finite(fpool$est.att[as.numeric(rownames(fpool$est.att)) >= 0,
                                          c("CI.lower", "CI.upper")])))
})

## ---- plotting ---------------------------------------------------------------

test_that("plot() works on conformal fits (block + staggered, gap + counterfactual)", {
  skip_on_cran()
  mkfit <- function(D, seed, ...) {
    set.seed(seed); N <- ncol(D); T <- nrow(D); r <- 2
    Fm <- matrix(rnorm(T * r), T, r); L <- matrix(rnorm(N * r), N, r)
    Y <- Fm %*% t(L) + matrix(rnorm(N * T), T, N) + 1.5 * D
    dat <- data.frame(id = rep(1:N, each = T), time = rep(1:T, N),
                      Y = as.vector(Y), D = as.vector(D))
    suppressWarnings(suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
          method = "gsynth", force = 3, CV = FALSE, r = r, se = TRUE,
          vartype = "conformal", parallel = FALSE, ...)))
  }
  Db <- matrix(0, 20, 25); Db[16:20, 1] <- 1
  Ds <- matrix(0, 20, 40); Ds[12:20, 1:4] <- 1; Ds[16:20, 5:8] <- 1
  fb  <- mkfit(Db, 1)
  fbs <- mkfit(Db, 1, conformal.band = "simultaneous")
  fs  <- mkfit(Ds, 2)
  for (f in list(fb, fbs, fs)) {
    expect_s3_class(plot(f, type = "gap"), "gg")
    expect_s3_class(plot(f, type = "counterfactual"), "gg")
  }
  ## the uniform band is stored and wider than the pointwise band in the post window
  post <- as.numeric(rownames(fb$est.att)) >= 0
  w_pt <- mean(fb$est.att[post, "CI.upper"] - fb$est.att[post, "CI.lower"])
  w_sim <- mean(fb$est.att.sim[post, "CI.upper"] - fb$est.att.sim[post, "CI.lower"])
  expect_gte(w_sim, w_pt)
  ## the inner (1 - 2*alpha) band (est.att90) is narrower than the main band
  w_in <- mean(fb$est.att90[post, "CI.upper"] - fb$est.att90[post, "CI.lower"])
  expect_lt(w_in, w_pt)
})

## ---- staggered adoption -----------------------------------------------------

test_that("conformal supports staggered adoption: runs, aligns event time, covers", {
  skip_on_cran()
  N <- 40; T <- 20; r <- 2
  D <- matrix(0, T, N); D[12:T, 1:4] <- 1; D[16:T, 5:8] <- 1   # cohorts at 12 and 16
  set.seed(20260609); reps <- 40; alpha <- 0.10
  cov <- logical(reps); bad <- 0L; aligned <- TRUE
  for (b in seq_len(reps)) {
    Fm <- matrix(rnorm(T * r), T, r); L <- matrix(rnorm(N * r), N, r)
    Y <- Fm %*% t(L) + matrix(rnorm(N * T), T, N)              # true effect 0
    dat <- data.frame(id = rep(1:N, each = T), time = rep(1:T, N),
                      Y = as.vector(Y), D = as.vector(D))
    f <- suppressWarnings(suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
            method = "gsynth", force = 3, CV = FALSE, r = r, se = TRUE,
            vartype = "conformal", conformal.scale = "sd", alpha = alpha, parallel = FALSE)))
    if (f$conformal$status != "ok") bad <- bad + 1L
    cov[b] <- (f$est.avg[1, "CI.lower"] <= 0) && (0 <= f$est.avg[1, "CI.upper"])
    if (!identical(rownames(f$est.att), as.character(f$time))) aligned <- FALSE
  }
  expect_true(aligned)              # est.att indexed by event time, matched to fit$time
  expect_equal(bad, 0L)             # never empty/unbounded
  expect_gt(mean(cov), 0.80)        # scalar coverage near nominal under staggering
})

test_that("simultaneous band gives materially better joint coverage than pointwise", {
  skip_on_cran()
  ## Joint (whole-post-path) coverage. Pointwise undercovers jointly (the
  ## multiple-comparison problem, ~0.44 at 150 reps); the simultaneous band is
  ## much better (~0.77). NOTE: it undercovers nominal 1 - alpha in this small-N
  ## regime, sitting near the jackknife+ floor 1 - 2*alpha plus the treated-vs-
  ## control LOO asymmetry; the precise characterization (and a possible
  ## symmetric-treated-gap fix) is a Phase 5 item. Here we assert the robust
  ## qualitative property only.
  g <- .conf_dgp(909); reps <- 40; alpha <- 0.10
  joint <- function(band) {
    post <- band[as.numeric(rownames(band)) >= 0, c("CI.lower", "CI.upper"), drop = FALSE]
    all(post[, 1] <= 0 & 0 <= post[, 2])
  }
  jp <- js <- logical(reps)
  for (b in seq_len(reps)) {
    Y <- g$mu + matrix(rnorm(g$N * g$T), g$T, g$N)   # true effect 0 everywhere
    dat <- data.frame(id = rep(1:g$N, each = g$T), time = rep(1:g$T, g$N),
                      Y = as.vector(Y), D = as.vector(g$D))
    f <- suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
              method = "gsynth", force = 3, CV = FALSE, r = 2, se = TRUE,
              vartype = "conformal", alpha = alpha, parallel = FALSE))
    jp[b] <- joint(f$est.att)        # pointwise
    js[b] <- joint(f$est.att.sim)    # simultaneous
  }
  expect_gt(mean(js) - mean(jp), 0.15)   # simultaneous materially better jointly
  expect_gt(mean(js), 0.65)              # and well above pointwise's joint coverage
})
