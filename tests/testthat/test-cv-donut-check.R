## `cv.donut` is checked against `cv.nobs` (added in 2.4.7). Block
## cross-validation holds out cv.nobs consecutive periods and scores only the
## middle cv.nobs - 2 * cv.donut of them, so 2 * cv.donut >= cv.nobs left
## nothing to score. Before the check, fect() stopped deep inside the CV with
## the internal message "No residuals to score." (e.g. cv.method = "block",
## cv.nobs = 3, cv.donut = 2). Rolling folds do not use cv.donut.

skip_on_cran()

.donut_panel <- function() {
  withr::with_seed(11, {
    N <- 30; TT <- 12; Ntr <- 10
    a <- rnorm(N); xi <- rnorm(TT)
    D <- matrix(0L, TT, N); D[8:TT, seq_len(Ntr)] <- 1L
    Y <- outer(xi, rep(1, N)) + outer(rep(1, TT), a) + 2 * D +
      matrix(rnorm(N * TT), TT, N)
    data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
               Y = c(Y), D = c(D))
  })
}

.fit_donut <- function(...) {
  suppressMessages(suppressWarnings(
    fect(Y ~ D, data = .donut_panel(), index = c("id", "time"), CV = TRUE,
         r = c(0, 1), k = 2, parallel = FALSE, seed = 1, ...)
  ))
}

test_that("block CV stops early with a clear message when cv.donut leaves nothing to score", {
  expect_error(.fit_donut(method = "ife", cv.method = "block", cv.nobs = 3, cv.donut = 2),
               "leaves no held-out period to score")
  expect_error(.fit_donut(method = "ife", cv.method = "block", cv.nobs = 4, cv.donut = 2),
               "Use cv.donut <= 1")
  ## never-treated CV (gsynth) takes the same check
  expect_error(.fit_donut(method = "gsynth", cv.method = "block", cv.nobs = 3, cv.donut = 2),
               "leaves no held-out period to score")
})

test_that("cv.donut must be a non-negative whole number", {
  expect_error(.fit_donut(method = "ife", cv.method = "block", cv.donut = -1),
               "non-negative whole number")
  expect_error(.fit_donut(method = "ife", cv.method = "block", cv.donut = 0.5),
               "non-negative whole number")
})

test_that("valid cv.donut values, and rolling CV, are unaffected", {
  ## largest valid donut for cv.nobs = 5 is 2 (one scored period per run)
  fit <- .fit_donut(method = "ife", cv.method = "block", cv.nobs = 5, cv.donut = 2)
  expect_true(fit$r.cv %in% 0:1)
  ## rolling folds ignore cv.donut
  fit <- .fit_donut(method = "ife", cv.method = "rolling", cv.nobs = 3, cv.donut = 2)
  expect_true(fit$r.cv %in% 0:1)
})
