## ---------------------------------------------------------------
## Regression test: standard rolling-window CV in r.cv.rolling.
##
## After the rewrite, r.cv.rolling implements the standard time-series
## rolling-window design (random anchors per unit per fold, with
## training restricted to observations strictly before each anchor;
## an optional cv.buffer excludes a buffer immediately before the anchor;
## all observations after anchor + cv.nobs are dropped from training).
##
## Tests cover: new parameters (cv.buffer, k, seed), reproducibility under
## a fixed seed, fold-aggregated MSPE structure, identical CV behavior
## across method = "ife" and method = "gsynth".
## ---------------------------------------------------------------

.make_unbalanced_panel <- function(seed = 11) {
  set.seed(seed)
  N  <- 30L; T <- 25L
  ids   <- rep(seq_len(N), each = T)
  times <- rep(seq_len(T), N)
  F1 <- rnorm(T); L1 <- rnorm(N)
  err <- rnorm(N * T, sd = 0.3)
  mu  <- rep(L1, each = T) * rep(F1, N)
  Y   <- mu + err
  treated      <- 26L:30L
  treat_starts <- c(18L, 19L, 20L, 21L, 22L)
  D <- rep(0L, N * T)
  for (k in seq_along(treated)) {
    u <- treated[k]
    D[ids == u & times >= treat_starts[k]] <- 1L
  }
  Y <- Y + D * 1.5
  df <- data.frame(id = ids, time = times, Y = Y, D = D)
  drop_keys <- c(
    paste(27L, c(24L, 25L), sep = "_"),
    paste(28L, c(22L, 23L, 24L, 25L), sep = "_"),
    paste(1L,  c(23L, 24L, 25L), sep = "_"),
    paste(2L,  c(24L, 25L), sep = "_")
  )
  df_keys <- paste(df$id, df$time, sep = "_")
  df[!df_keys %in% drop_keys, ]
}

test_that("k and cv.buffer parameters are honored", {

  skip_on_cran()
  df <- .make_unbalanced_panel(seed = 11)

  set.seed(99)
  res_k5 <- suppressWarnings(suppressMessages(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r.max = 2, cv.nobs = 2, cv.buffer = 1,
      k = 5L, seed = 99,
      verbose = FALSE
    )
  ))
  expect_equal(res_k5$k, 5L)
  expect_equal(res_k5$cv.buffer, 1L)
  expect_equal(res_k5$cv.nobs, 2L)
  expect_equal(dim(res_k5$mspe.per.fold),
               c(3L, 5L))    # (r.max + 1) rows x k columns
  expect_true(all(res_k5$mspe$n_folds_used <= 5L))
})

test_that("seed gives reproducible per-fold MSPE", {

  skip_on_cran()
  df <- .make_unbalanced_panel(seed = 11)

  res_a <- suppressWarnings(suppressMessages(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r.max = 2, cv.nobs = 2, cv.buffer = 1,
      k = 3L, seed = 42L, verbose = FALSE
    )
  ))
  res_b <- suppressWarnings(suppressMessages(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r.max = 2, cv.nobs = 2, cv.buffer = 1,
      k = 3L, seed = 42L, verbose = FALSE
    )
  ))
  ## The seed is per-fold; identical inputs should give identical fold MSPE.
  expect_equal(res_a$mspe.per.fold, res_b$mspe.per.fold)
  expect_equal(res_a$r.cv, res_b$r.cv)
})

test_that("rolling-window CV identical for method = 'ife' and 'gsynth' inner-call paths", {

  skip_on_cran()
  ## The two methods route to fect with the same masking logic;
  ## with the same seed they should sample the same per-fold anchors.
  ## Numerical results may differ because of different inner estimators,
  ## but the anchor sets and units_seen should match.
  df <- .make_unbalanced_panel(seed = 11)

  res_ife <- suppressWarnings(suppressMessages(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r.max = 1, cv.nobs = 2, cv.buffer = 0,
      k = 2L, seed = 7L, verbose = FALSE
    )
  ))
  res_gsc <- suppressWarnings(suppressMessages(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "gsynth", r.max = 1, cv.nobs = 2, cv.buffer = 0,
      k = 2L, seed = 7L, verbose = FALSE
    )
  ))
  expect_equal(res_ife$n.units.masked, res_gsc$n.units.masked)
  expect_equal(res_ife$k, res_gsc$k)
  expect_equal(res_ife$cv.nobs, res_gsc$cv.nobs)
  expect_equal(res_ife$cv.buffer, res_gsc$cv.buffer)
})

test_that("invalid arguments are rejected", {

  skip_on_cran()
  df <- .make_unbalanced_panel(seed = 11)

  expect_error(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r.max = 1, cv.nobs = 0, verbose = FALSE
    ),
    regexp = "cv\\.nobs must be"
  )
  expect_error(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r.max = 1, cv.buffer = -1, verbose = FALSE
    ),
    regexp = "cv\\.buffer must be"
  )
  expect_error(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r.max = 1, k = 0, verbose = FALSE
    ),
    regexp = "k must be"
  )
  expect_error(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r.max = 1, cv.prop = 0, verbose = FALSE
    ),
    regexp = "cv\\.prop"
  )
  expect_error(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r.max = 1, cv.prop = 1.5, verbose = FALSE
    ),
    regexp = "cv\\.prop"
  )
})

test_that("cv.prop controls per-fold unit sampling", {

  skip_on_cran()
  ## With 30 units (25 controls + 5 treated, all 5 treated have at
  ## least min.T0 + cv.nobs pre-treatment cells in the fixture), the
  ## eligible pool covers controls AND treated pre-treatment.
  df <- .make_unbalanced_panel(seed = 11)

  ## Small cv.prop -> few units sampled per fold; capped at >= 1.
  res_small <- suppressWarnings(suppressMessages(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r.max = 1, cv.nobs = 2, cv.buffer = 1,
      k = 2L, cv.prop = 0.05, seed = 17L, verbose = FALSE
    )
  ))
  expect_equal(res_small$cv.prop, 0.05)
  ## At cv.prop = 0.05 with ~30 eligible units -> max(1, round(1.5)) = 2
  ## sampled per fold; with k = 2 folds, units_seen <= 4.
  expect_lte(res_small$n.units.masked, 4L)
  expect_gte(res_small$n.units.masked, 1L)

  ## Full cv.prop -> every eligible unit sampled every fold.
  res_full <- suppressWarnings(suppressMessages(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r.max = 1, cv.nobs = 2, cv.buffer = 1,
      k = 2L, cv.prop = 1.0, seed = 17L, verbose = FALSE
    )
  ))
  expect_equal(res_full$cv.prop, 1.0)
  ## n.units.masked at cv.prop = 1.0 should exceed n.units.masked at
  ## cv.prop = 0.05 by a wide margin.
  expect_gt(res_full$n.units.masked, res_small$n.units.masked)
})

test_that("treated pre-treatment cells are eligible for masking", {

  skip_on_cran()
  ## With cv.prop = 1.0, every eligible unit gets masked. We then
  ## check that at least one treated unit appears in the holdout
  ## (would be impossible under the old controls-only design).
  df <- .make_unbalanced_panel(seed = 11)
  treated_ids <- 26:30  # from .make_unbalanced_panel

  res <- suppressWarnings(suppressMessages(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r.max = 1, cv.nobs = 2, cv.buffer = 1,
      k = 2L, cv.prop = 1.0, seed = 21L, verbose = FALSE
    )
  ))
  ## n.units.masked should include at least one treated unit (their
  ## pre-treatment cells are eligible). Total eligible >= 25
  ## controls + the treated units with sufficient pre-period.
  expect_gt(res$n.units.masked, 25L)
})
