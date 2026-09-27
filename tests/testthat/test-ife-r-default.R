## `r` defaults to NULL (since 2.4.7). For method = "ife" / "gsynth" a missing r
## means: cross-validate the number of factors over 0:5, as method = "mc" does
## for a missing lambda. With CV = FALSE and no r, fall back to r = 0 (FEct)
## with a message. An explicit r = 0 is honoured silently. Before 2.4.7 the
## default was r = 0, so method = "ife" without r silently ran the two-way FE
## model: the manual's own LOO example did (reader report, 2026-09-21).

skip_on_cran()

.small_panel <- function() {
  data(sim_base, package = "fect")
  sim_base[sim_base$id %in% unique(sim_base$id)[1:40], ]
}

test_that("ife without r cross-validates the number of factors", {
  d <- .small_panel()
  out <- suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"),
                               method = "ife", se = FALSE, k = 3, parallel = FALSE))
  expect_false(is.null(out$r.cv))
  expect_true(out$r.cv %in% 0:5)
})

test_that("ife with CV = FALSE and no r falls back to r = 0 with a message", {
  d <- .small_panel()
  expect_message(
    out <- fect(Y ~ D, data = d, index = c("id", "time"), method = "ife",
                CV = FALSE, se = FALSE, parallel = FALSE),
    "No r is supplied and CV = FALSE"
  )
  expect_equal(out$r.cv, 0)
})

test_that("an explicit r, fe, and CV = TRUE stay quiet", {
  d <- .small_panel()
  expect_no_message(
    fect(Y ~ D, data = d, index = c("id", "time"), method = "ife", r = 0,
         se = FALSE, parallel = FALSE),
    message = "No r is supplied"
  )
  expect_no_message(
    fect(Y ~ D, data = d, index = c("id", "time"), method = "ife", r = 2,
         se = FALSE, parallel = FALSE),
    message = "No r is supplied"
  )
  expect_no_message(
    fect(Y ~ D, data = d, index = c("id", "time"), se = FALSE, parallel = FALSE),
    message = "No r is supplied"
  )
  expect_no_message(
    fect(Y ~ D, data = d, index = c("id", "time"), method = "ife",
         CV = TRUE, r = c(0, 2), k = 3, se = FALSE, parallel = FALSE),
    message = "No r is supplied"
  )
})
