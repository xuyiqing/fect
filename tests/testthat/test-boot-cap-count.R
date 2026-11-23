## boot.diag$n.cap / n.cap.unknown: bootstrap fits that stopped at max.iteration.
## Written from the behaviour description only.

skip_on_cran()

cap_text <- "bootstrap fits stopped at max.iteration"

cap_data <- function() {
  data(simdata, package = "fect", envir = environment())
  simdata[simdata$id %in% unique(simdata$id)[1:40], ]
}

## fit and collect every warning raised
cap_fit <- function(max.iteration = 1000, parallel = FALSE, ...) {
  d <- cap_data()
  w <- character(0)
  out <- withCallingHandlers(
    suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"),
      method = "ife", r = 1, CV = FALSE, force = "two-way",
      se = TRUE, nboots = 20, seed = 123, parallel = parallel,
      cores = if (parallel) 2 else NULL,
      max.iteration = max.iteration, ...)),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  list(fit = out, warnings = w, n.cap.warn = sum(grepl(cap_text, w, fixed = TRUE)))
}

print_lines <- function(fit) capture.output(print(fit))

test_that("a fit forced to the cap counts, warns once and prints a line", {
  r <- cap_fit(max.iteration = 1)
  bd <- r$fit$boot.diag
  expect_false(is.null(bd))
  expect_gt(bd$n.cap, 0L)
  expect_lte(bd$n.cap + bd$n.cap.unknown, bd$n.used)
  expect_equal(r$n.cap.warn, 1L)
  expect_true(any(grepl("max.iteration", print_lines(r$fit), fixed = TRUE)))
})

test_that("parallel and serial fits give the same cap counts", {
  skip_on_os("windows")
  s <- cap_fit(max.iteration = 1, parallel = FALSE)
  p <- cap_fit(max.iteration = 1, parallel = TRUE)
  expect_equal(p$fit$boot.diag$n.cap, s$fit$boot.diag$n.cap)
  expect_equal(p$fit$boot.diag$n.cap.unknown, s$fit$boot.diag$n.cap.unknown)
})

test_that("a converging fit has n.cap = 0, no warning and no print line", {
  r <- cap_fit()
  expect_equal(r$fit$boot.diag$n.cap, 0L)
  expect_equal(r$n.cap.warn, 0L)
  expect_false(any(grepl("stopped at max.iteration", print_lines(r$fit), fixed = TRUE)))
})

test_that("n.cap and n.cap.unknown are length-one non-negative whole numbers", {
  check <- function(bd) {
    for (f in c("n.cap", "n.cap.unknown")) {
      x <- bd[[f]]
      expect_length(x, 1L)
      expect_true(is.numeric(x))
      expect_true(is.finite(x) && x >= 0 && x == round(x))
    }
  }
  check(cap_fit(max.iteration = 1)$fit$boot.diag)
  check(cap_fit()$fit$boot.diag)
  d <- cap_data()
  j <- suppressWarnings(suppressMessages(fect(Y ~ D, data = d,
    index = c("id", "time"), method = "ife", r = 1, CV = FALSE, se = TRUE,
    vartype = "jackknife", parallel = FALSE)))
  if (!is.null(j$boot.diag)) {
    check(j$boot.diag)
  } else {
    succeed("no boot.diag for jackknife")
  }
})

test_that("the count does not change results: two runs agree", {
  a <- cap_fit()$fit
  b <- cap_fit()$fit
  expect_true(all(is.finite(a$att.avg)))
  expect_true(all(is.finite(a$est.avg)))
  expect_true(all(is.finite(a$est.att[, "ATT"])))
  expect_identical(a$att.avg, b$att.avg)
  expect_identical(a$est.avg, b$est.avg)
  expect_identical(a$est.att, b$est.att)
})
