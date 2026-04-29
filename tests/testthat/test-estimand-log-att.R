## ---------------------------------------------------------------
## Tests for estimand(fit, "log.att", "event.time").
##
## logATT_t = mean_{(t,i): D=1} (log(Y) - log(Y0_hat))
##
## Validates against a manual computation.
## ---------------------------------------------------------------

## simdata's Y can be negative, which makes log undefined for many
## cells. Build a small synthetic positive-Y panel for testing.

.fit_positive_Y <- function(nboots = 30, keep.sims = TRUE, seed = 7) {
  set.seed(seed)
  N <- 30; TT <- 20
  df <- expand.grid(id = 1:N, time = 1:TT)
  treat_start <- sample(c(NA, 8:15), N, replace = TRUE)
  df$D <- ifelse(is.na(treat_start[df$id]) | df$time < treat_start[df$id],
                 0, 1)
  df$Y <- exp(0.5 + 0.05 * df$time + 0.3 * df$D + rnorm(nrow(df), sd = 0.2))
  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, nboots = nboots, parallel = FALSE,
      keep.sims = keep.sims
    )
  ))
}


## -- LA.1  Schema -----------------------------------------------------

test_that("LA.1: estimand(fit, 'log.att', 'event.time') has the documented schema", {

  skip_on_cran()

  fit <- .fit_positive_Y()
  est <- suppressWarnings(fect::estimand(fit, "log.att", "event.time"))

  expect_s3_class(est, "data.frame")
  expected_cols <- c("event.time", "estimate", "se",
                     "ci.lo", "ci.hi", "n_cells", "vartype")
  expect_setequal(names(est), expected_cols)

  expect_true(is.numeric(est$estimate))
})


## -- LA.2  Manual recipe match ----------------------------------------

test_that("LA.2: estimate matches manual log-ATT", {

  skip_on_cran()

  fit <- .fit_positive_Y()
  est <- suppressWarnings(fect::estimand(fit, "log.att", "event.time"))

  Don <- !is.na(fit$D.dat) & fit$D.dat == 1 & !is.na(fit$T.on)
  ets <- sort(unique(fit$T.on[Don]))

  manual <- vapply(ets, function(et) {
    m <- Don & fit$T.on == et
    eff_t <- fit$eff[m]
    Y_t   <- fit$Y.dat[m]
    Y0_t  <- Y_t - eff_t
    ok <- !is.na(Y_t) & !is.na(Y0_t) & Y_t > 0 & Y0_t > 0
    if (sum(ok) == 0) return(NA_real_)
    mean(log(Y_t[ok]) - log(Y0_t[ok]), na.rm = TRUE)
  }, numeric(1))

  est_sorted <- est[order(est$event.time), ]
  expect_equal(est_sorted$estimate, manual, tolerance = 1e-10)
})


## -- LA.3  keep.sims = FALSE error path -------------------------------

test_that("LA.3: log.att without keep.sims errors with the locked wording", {

  skip_on_cran()

  fit_no_keep <- .fit_positive_Y(keep.sims = FALSE)
  expect_error(
    fect::estimand(fit_no_keep, "log.att", "event.time"),
    "keep\\.sims = TRUE",
    fixed = FALSE
  )
})


## -- LA.4  Negative Y cells: warning + drop -------------------------

test_that("LA.4: negative Y cells trigger a one-time warning and are dropped", {

  skip_on_cran()

  ## simdata has negative Y values.
  data("simdata", package = "fect")
  set.seed(42)
  fit_neg <- suppressWarnings(suppressMessages(
    fect::fect(Y ~ D, data = simdata, index = c("id", "time"),
               method = "fe", force = "two-way",
               se = TRUE, nboots = 20, parallel = FALSE,
               keep.sims = TRUE)
  ))

  expect_warning(
    fect::estimand(fit_neg, "log.att", "event.time"),
    "dropped",
    fixed = FALSE
  )
})
