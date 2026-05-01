## ---------------------------------------------------------------
## Integration tests: cv.method = "rolling" wired into the main
## fect() CV dispatcher (additive change; block CV defaults
## untouched).
##
## Coverage:
##   1. fect(method = "ife", CV = TRUE, cv.method = "rolling", ...)
##   2. fect(method = "gsynth", CV = TRUE, cv.method = "rolling", ...)
##      (routes through fect_nevertreated)
##   3. fect(method = "cfe", CV = TRUE, cv.method = "rolling", ...)
##      with extra-FE column in index
##   4. fect_mspe(out, cv.method = "rolling", cv.buffer = 1, ...)
##   5. Regression: existing block CV defaults still work and produce
##      a finite r.cv.
## ---------------------------------------------------------------

.make_rolling_panel <- function(seed = 17) {
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
  ## Make panel slightly unbalanced (so masking some control units'
  ## tail observations doesn't drop entire time points).
  drop_keys <- c(
    paste(27L, c(24L, 25L), sep = "_"),
    paste(28L, c(22L, 23L, 24L, 25L), sep = "_"),
    paste(1L,  c(23L, 24L, 25L), sep = "_"),
    paste(2L,  c(24L, 25L), sep = "_")
  )
  df_keys <- paste(df$id, df$time, sep = "_")
  df[!df_keys %in% drop_keys, ]
}

test_that("fect dispatcher: cv.method = 'rolling' works for method = 'ife'", {

  skip_on_cran()

  df <- .make_rolling_panel(seed = 17)

  set.seed(101)
  out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", force = "two-way",
      CV = TRUE, cv.method = "rolling",
      cv.buffer = 1, cv.nobs = 2, k = 2, cv.prop = 0.2,
      r = c(0, 2), min.T0 = 5,
      se = FALSE
    )
  ))

  expect_true(is.list(out))
  expect_true(!is.null(out$r.cv))
  expect_true(is.finite(as.numeric(out$r.cv)))
})

test_that("fect dispatcher: cv.method = 'rolling' works for method = 'gsynth'", {

  skip_on_cran()

  df <- .make_rolling_panel(seed = 19)

  set.seed(202)
  out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "gsynth", force = "two-way",
      CV = TRUE, cv.method = "rolling",
      cv.buffer = 1, cv.nobs = 2, k = 2, cv.prop = 0.2,
      r = c(0, 2), min.T0 = 5,
      se = FALSE
    )
  ))

  expect_true(is.list(out))
  expect_true(!is.null(out$r.cv))
  expect_true(is.finite(as.numeric(out$r.cv)))
})

test_that("fect dispatcher: cv.method = 'rolling' works for method = 'cfe' with extra-FE", {

  skip_on_cran()

  if (!exists("sim_region", inherits = TRUE)) {
    suppressWarnings(try(utils::data("sim_region", package = "fect"),
                         silent = TRUE))
  }
  skip_if_not(exists("sim_region", inherits = TRUE),
              "sim_region dataset not available")

  set.seed(303)
  out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = sim_region,
      index = c("id", "time", "region_time"),
      method = "cfe", force = "two-way",
      CV = TRUE, cv.method = "rolling",
      cv.buffer = 1, cv.nobs = 2, k = 2, cv.prop = 0.2,
      r = c(0, 1), min.T0 = 5,
      se = FALSE
    )
  ))

  expect_true(is.list(out))
  expect_true(!is.null(out$r.cv))
  expect_true(is.finite(as.numeric(out$r.cv)))
})

test_that("fect_mspe: cv.method = 'rolling' produces a finite MSPE", {

  skip_on_cran()

  df <- .make_rolling_panel(seed = 23)

  set.seed(404)
  ## A single fect output is the input.
  out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", force = "two-way",
      CV = FALSE, r = 1, se = FALSE
    )
  ))

  res <- suppressWarnings(suppressMessages(
    fect::fect_mspe(
      out, cv.method = "rolling",
      cv.buffer = 1, cv.nobs = 2, k = 2, cv.prop = 0.2,
      min.T0 = 5, seed = 7
    )
  ))

  expect_true(is.list(res))
  expect_true(is.data.frame(res$summary))
  expect_true(any(is.finite(res$summary$MSPE)))
})

test_that("fect dispatcher: block-CV defaults still work (regression)", {

  skip_on_cran()

  df <- .make_rolling_panel(seed = 31)

  set.seed(505)
  ## Use the existing block CV path with all defaults.
  out_block <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", force = "two-way",
      CV = TRUE, cv.method = "all_units",
      cv.nobs = 2, k = 2, cv.prop = 0.2,
      r = c(0, 2), min.T0 = 5,
      se = FALSE
    )
  ))

  expect_true(is.list(out_block))
  expect_true(!is.null(out_block$r.cv))
  expect_true(is.finite(as.numeric(out_block$r.cv)))
})
