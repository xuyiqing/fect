test_that("fect_sens stops on a fit without a placebo test", {
  skip_on_cran()
  skip_unless_honest_solvers()
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
  expect_true(!is.null(out$est.att))
  expect_true(!is.null(out$att.vcov))

  out$placebo.period <- c(-3, -1)
  expect_error(fect_sens(out, post.periods = 1:5, parallel = FALSE), "placebo periods")
})

test_that("fect_sens attaches sensitivity results when inputs present", {
  skip_on_cran()
  skip_unless_honest_solvers()
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
    parallel = FALSE,
    placeboTest = TRUE,
    placebo.period = c(-2, 0)
  )
  expect_true(!is.null(out$est.att))
  expect_true(!is.null(out$att.vcov))

  # grids are reduced from the defaults to keep the test short
  out2 <- fect_sens(out, post.periods = 1:4, Mbarvec = c(0, 0.5), Mvec = c(0, 0.1),
                    periodMbarvec = 0.5, periodMvec = 0.1, parallel = FALSE)
  expect_false(is.null(out2$sensitivity.rm$results))
  expect_false(is.null(out2$sensitivity.smooth$results))
  expect_equal(nrow(out2$sensitivity.rm$results), 2L)
  expect_equal(nrow(out2$sensitivity.smooth$results), 2L)
})

## 2.4.8: the default post.periods, parallel = TRUE with cores, and NULL
## vectors no longer stop fect_sens() or drop the overall results.
test_that("fect_sens defaults, parallel with cores, and NULL vectors", {
  skip_on_cran()
  skip_unless_honest_solvers()
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

  ## serial and parallel give the same bounds, in both restrictions
  args <- list(out, post.periods = post, Mbarvec = c(0, 0.5), Mvec = c(0, 0.1), periodMbarvec = NULL, periodMvec = NULL)
  ser <- do.call(fect_sens, c(args, list(parallel = FALSE)))
  par <- do.call(fect_sens, c(args, list(parallel = TRUE, cores = 2)))
  expect_close(par$sensitivity.rm$results$lb,  ser$sensitivity.rm$results$lb,  1e-10)
  expect_close(par$sensitivity.rm$results$ub,  ser$sensitivity.rm$results$ub,  1e-10)
  expect_close(par$sensitivity.smooth$results$lb, ser$sensitivity.smooth$results$lb, 1e-10)
  expect_close(par$sensitivity.smooth$results$ub, ser$sensitivity.smooth$results$ub, 1e-10)
  expect_identical(par$sensitivity.rm$results$Mbar, ser$sensitivity.rm$results$Mbar)
  expect_identical(par$sensitivity.smooth$results$M, ser$sensitivity.smooth$results$M)

  ## a second parallel call works
  par2 <- do.call(fect_sens, c(args, list(parallel = TRUE, cores = 2)))
  expect_close(par2$sensitivity.rm$results$lb, ser$sensitivity.rm$results$lb, 1e-10)
})

test_that("fect_sens names the missing solver packages", {
  expect_error(
    fect:::.fect_sens_check_packages(c("tibble", "notInstalledPkgA", "notInstalledPkgB")),
    "^fect_sens\\(\\) needs the package\\(s\\) 'notInstalledPkgA', 'notInstalledPkgB'\\. Install with:\\ninstall\\.packages\\(c\\(\"notInstalledPkgA\", \"notInstalledPkgB\"\\)\\)$")
  expect_silent(fect:::.fect_sens_check_packages(character(0)))
  expect_identical(fect:::.fect_sens_needs(FALSE, FALSE, 3L), character(0))
  expect_true(all(c("tibble", "purrr", "TruncatedNormal", "matrixStats", "Rglpk") %in% fect:::.fect_sens_needs(FALSE, TRUE, 1L)))
  expect_false("lpSolveAPI" %in% fect:::.fect_sens_needs(FALSE, TRUE, 1L))
  expect_true(all(c("lpSolveAPI", "Matrix", "pracma") %in% fect:::.fect_sens_needs(FALSE, TRUE, 3L)))
  expect_true("pracma" %in% fect:::.fect_sens_needs(TRUE, FALSE, 1L))
})
