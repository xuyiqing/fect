test_that("fect_sens attaches sensitivity results when inputs present", {
  skip_on_cran()
  skip_if_not_installed("HonestDiDFEct")
  suppressWarnings(try(data("simdata", package = "fect"), silent = TRUE))
  expect_true(exists("simdata"))

  set.seed(3)
  out <- fect::fect(
    Y ~ D + X1 + X2,
    data = simdata,
    index = c("id", "time"),
    method = "ife",
    r = 1,
    CV = FALSE,
    se = TRUE,
    nboots = 20,
    parallel = FALSE
  )

  # ensure names and shapes needed by fect_sens exist
  expect_true(!is.null(out$est.att))
  expect_true(!is.null(out$att.vcov))

  # Ensure placebo.period exists to avoid length-0 error inside fect_sens
  if (is.null(out$placebo.period)) {
    out$placebo.period <- c(-3, -1)
  }
  # restrict periods so numPre+numPost matches available row count
  out$placebo.period <- c(-3, -1)
  out2 <- fect_sens(out, post.periods = 1:10, parallel = FALSE)
  expect_true(!is.null(out2$sensitivity.rm) || !is.null(out2$sensitivity.smooth))
})

## 2.4.8: the default post.periods, parallel = TRUE with cores, and NULL
## vectors no longer stop fect_sens() or drop the overall results.
test_that("fect_sens defaults, parallel with cores, and NULL vectors", {
  skip_on_cran()
  skip_if_not_installed("HonestDiDFEct")
  data("simdata", package = "fect", envir = environment())
  set.seed(3)
  out <- suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
    method = "fe", se = TRUE, nboots = 40, parallel = FALSE,
    placeboTest = TRUE, placebo.period = c(-2, 0)
  ))
  post <- 1:5

  ## default post.periods: event times 1..max, no error
  s0 <- fect_sens(out, Mbarvec = c(0, 0.5), Mvec = NULL,
                  periodMbarvec = NULL, periodMvec = NULL)
  expect_false(is.null(s0$sensitivity.rm$results))
  expect_null(s0$sensitivity.rm$periods)        # periodMbarvec = NULL keeps the overall result
  expect_null(s0$sensitivity.smooth)

  ## Mbarvec = NULL with periodMbarvec set: periods only
  s1 <- fect_sens(out, post.periods = post, Mbarvec = NULL, Mvec = NULL,
                  periodMbarvec = c(0, 0.5), periodMvec = NULL)
  expect_null(s1$sensitivity.rm$results)
  expect_true(is.data.frame(s1$sensitivity.rm$periods))
  expect_setequal(unique(s1$sensitivity.rm$periods$postPeriod), post)

  ## parallel = TRUE with cores given
  s2 <- fect_sens(out, post.periods = post, Mbarvec = c(0, 0.5), Mvec = NULL,
                  periodMbarvec = NULL, periodMvec = NULL,
                  parallel = TRUE, cores = 2)
  expect_false(is.null(s2$sensitivity.rm$results))
})
