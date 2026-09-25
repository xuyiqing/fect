## ---------------------------------------------------------------
## fix246 section B: CV, weights in CV, implied weights, plots, guards,
## never-treated se = TRUE, cumulative ATT (items B1-B13).
##
## Regression tests for the 2.4.6 correctness run. Each block fails on
## fect dev @ 412d7ae and passes after the fix. Fits are small and serial
## (parallel = FALSE); helpers use pkg::fn for non-base functions.
## ---------------------------------------------------------------

.fb_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## Run `expr`, returning its value (or the error) and every message text.
.fb_msgs <- function(expr) {
  m <- character(0)
  res <- tryCatch(
    withCallingHandlers(
      suppressWarnings(expr),
      message = function(cnd) {
        m <<- c(m, conditionMessage(cnd))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) e
  )
  list(res = res, msgs = m)
}

.fb_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "fect", envir = e)
  e[[name]]
}


## -- B1  criterion = "pc" selects the r with the lowest PC -----------------

test_that("B1: criterion = 'pc' picks argmin PC for gsynth and cfe + never-treated", {
  skip_on_cran()
  turnout <- .fb_data("turnout")
  for (m in c("gsynth", "cfe")) {
    args <- list(turnout ~ policy_edr + policy_mail_in + policy_motor,
                 data = turnout, index = c("abb", "year"), method = m,
                 force = "two-way", CV = TRUE, r = c(0, 5), criterion = "pc",
                 se = FALSE, parallel = FALSE, seed = 1)
    if (m == "cfe") args$time.component.from <- "nevertreated"
    fit <- .fb_quiet(do.call(fect::fect, args))
    cv <- fit$CV.out
    ## PC is lowest at r = 4 on turnout; the MSPE rules pick r = 2
    expect_equal(unname(cv[which.min(cv[, "PC"]), "r"]), 4, info = m)
    expect_equal(as.numeric(fit$r.cv), 4, info = m)
  }
})

test_that("B1: fect's IFE CV table under 'pc' has the PC columns, correctly labelled", {
  skip_on_cran()
  simdata <- .fb_data("simdata")
  fit <- .fb_quiet(fect::fect(Y ~ D + X1 + X2, data = simdata,
                              index = c("id", "time"), method = "ife",
                              CV = TRUE, r = c(0, 3), criterion = "pc",
                              se = FALSE, parallel = FALSE, seed = 1))
  cv <- fit$CV.out.ife
  expect_identical(colnames(cv), c("r", "sigma2", "IC", "PC", "MSPTATT", "MSE"))
  ## the two last columns hold the values computed for them (finite), and
  ## no MSPE-family column carries MSPTATT / MSE values
  expect_true(all(is.finite(cv[, "MSPTATT"])))
  expect_true(all(is.finite(cv[, "MSE"])))
  expect_true(all(cv[, "MSE"] > cv[, "MSPTATT"]))
  ## fect's own IFE path already selects argmin PC
  expect_equal(as.numeric(fit$r.cv), unname(cv[which.min(cv[, "PC"]), "r"]))
})
