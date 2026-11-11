## vartype = "conformal" on the code merged with fect 2.4.8 (#181-#184):
## alpha guard, failure counts and messages, non-converged fold fits, print(),
## summary(), boot.diag, plots, fect_sens(), and documentation facts.
## Spec: statsclaw-workspace/fect/runs/2026-10-28-fect-248-assembly/conformal/
## test-spec.md, sections 4.5 to 4.10. Every fit uses parallel = FALSE (mocks do
## not reach future workers) and CV = FALSE.

.c8_sim <- function(J = 20, Ntr = 1, TT = 20, T0 = 15, r = 1, seed = 1,
                    adopt = NULL, tau = 1) {
  set.seed(seed)
  N <- J + Ntr
  F <- matrix(rnorm(TT * r), TT, r)
  L <- matrix(rnorm(N * r), N, r)
  Y0 <- outer(rnorm(TT), rep(1, N)) + outer(rep(1, TT), rnorm(N)) +
    F %*% t(L) + matrix(rnorm(TT * N), TT, N)
  if (is.null(adopt)) adopt <- rep(T0 + 1, Ntr)
  D <- matrix(0, TT, N)
  for (i in seq_len(Ntr)) D[seq_len(TT) >= adopt[i], i] <- 1
  Y <- Y0 + tau * D
  data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
             Y = c(Y), D = c(D))
}

.c8_fit <- function(df, ...) {
  suppressMessages(fect(Y ~ D, data = df, index = c("id", "time"),
                        se = TRUE, vartype = "conformal", parallel = FALSE,
                        CV = FALSE, ...))
}

## fit with failures injected into .conformal_component_fit; collects warnings
.c8_run <- function(df, fail.calls = integer(0), nan.calls = integer(0),
                    nan.fun = NULL, ...) {
  orig <- fect:::.conformal_component_fit
  st <- new.env()
  st$n <- 0L
  mock <- function(...) {
    st$n <- st$n + 1L
    if (st$n %in% fail.calls) return(NULL)
    m <- orig(...)
    if (st$n %in% nan.calls) m <- nan.fun(m)
    m
  }
  testthat::local_mocked_bindings(.conformal_component_fit = mock,
                                  .package = "fect")
  ws <- character(0)
  fit <- withCallingHandlers(
    .c8_fit(df, ...),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, ws = ws, calls = st$n)
}

.c8_collect <- function(expr) {
  ws <- character(0)
  val <- withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = val, ws = ws)
}

## ---- 4.5 alpha guard ------------------------------------------------------

test_that("alpha must be a single number in (0, 1); 0.5 and 0.9 run", {
  skip_on_cran()
  df <- .c8_sim(J = 20, Ntr = 1, seed = 5)
  for (a in list(0, 1, -0.1, 1.5, NA_real_, c(0.05, 0.1))) {
    expect_error(.c8_fit(df, method = "gsynth", r = 1, alpha = a), "alpha",
                 info = paste(format(a), collapse = ","))
  }
  f5 <- .c8_fit(df, method = "gsynth", r = 1, alpha = 0.5)
  u <- f5$conformal$unit
  expect_identical(unique(u$N), 21)
  expect_identical(unique(u$q), 5)
  expect_identical(as.integer(unique(u$q)), fect:::.conformal_q(21, 0.5))
  expect_true(all(is.finite(u$CI.lower) & is.finite(u$CI.upper)))
  f9 <- .c8_fit(df, method = "gsynth", r = 1, alpha = 0.9)
  expect_identical(unique(f9$conformal$unit$q), 9)
  f05 <- .c8_fit(df, method = "gsynth", r = 1)
  expect_identical(unique(f05$conformal$unit$q), 1)
  ## one treated unit: no caveat line for aggregates
  expect_false(any(grepl("^Aggregates:", capture.output(print(f05)))))
  expect_identical(unique(f05$conformal$unit$N), 21)
})

## ---- 4.6 failure counts, warning and status ------------------------------

test_that("failure counts, the warning and the status agree", {
  skip_on_cran()
  df <- .c8_sim(J = 40, Ntr = 1, TT = 20, T0 = 15, r = 1, seed = 21, tau = 6)
  args <- list(method = "gsynth", r = 1, alpha = 0.1)
  ## one fit returns NULL
  r1 <- do.call(.c8_run, c(list(df, fail.calls = 7), args))
  expect_length(r1$ws, 1)
  expect_match(r1$ws, "^vartype = \"conformal\": 1 fold\\(s\\) failed \\(1 fit\\(s\\) returned no prediction\\)")
  expect_match(r1$ws, "widened")
  cf <- r1$fit$conformal
  expect_equal(cf$n.failed.folds, 1)
  expect_equal(cf$n.failed.fits, 1)
  expect_match(cf$status, "; 1 failed fold\\(s\\), affected intervals widened$")
  expect_true(all(cf$unit$N == 41) && all(cf$unit$q == 2))
  ## the fit returns a matrix with NaN in the control's query column
  r2 <- do.call(.c8_run, c(list(df, nan.calls = 7,
                                nan.fun = function(m) { m[, ncol(m)] <- NaN; m }),
                           args))
  expect_length(r2$ws, 1)
  expect_match(r2$ws, "^vartype = \"conformal\": 1 fold\\(s\\) failed \\(0 fit\\(s\\) returned no prediction\\)")
  cf2 <- r2$fit$conformal
  expect_equal(cf2$n.failed.folds, 1)
  expect_equal(cf2$n.failed.fits, 0)
  expect_true(all(cf2$unit$n.failed == 1))
  expect_match(cf2$status, "1 failed fold\\(s\\)")
  ## no failure
  r0 <- do.call(.c8_run, c(list(df), args))
  expect_length(r0$ws, 0)
  expect_identical(r0$fit$conformal$status, "ok")
  expect_equal(r0$fit$conformal$n.failed.folds, 0)
})

## ---- 4.7 non-converged fold fits are not failures -------------------------

test_that("fold fits stopped by max.iteration are used, not counted as failed", {
  skip_on_cran()
  df <- .c8_sim(J = 25, Ntr = 1, seed = 5)
  for (args in list(list(method = "ife", r = 1),
                    list(method = "mc", lambda = 0.01))) {
    res <- .c8_collect(do.call(.c8_fit, c(list(df, max.iteration = 2, tol = 1e-12),
                                          args)))
    cf <- res$value$conformal
    expect_false(any(grepl("vartype = \"conformal\"", res$ws, fixed = TRUE)),
                 info = args$method)
    expect_equal(cf$n.failed.fits, 0)
    expect_equal(cf$n.failed.folds, 0)
    expect_true(all(cf$unit$n.failed == 0))
    expect_identical(cf$status, "ok")
    expect_true(all(is.finite(cf$unit$CI.lower) & is.finite(cf$unit$CI.upper)))
  }
})

## ---- 4.8 print(), summary(), boot.diag and plots --------------------------

test_that("conformal fits carry no bootstrap slots and print cleanly", {
  skip_on_cran()
  df <- .c8_sim(J = 40, Ntr = 2, seed = 3)
  fits <- list(gsynth = .c8_fit(df, method = "gsynth", r = 1),
               ife = .c8_fit(df, method = "ife", r = 1))
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    expect_false(any(c("boot.diag", "n.boot.fail", "para.sim.check") %in% names(fit)),
                 info = nm)
    expect_null(fit$boot.diag)
    expect_identical(fit$vartype, "conformal")
    out <- NULL
    expect_no_warning(out <- capture.output(print(fit)))
    expect_no_warning(capture.output(summary(fit)))
    expect_true(any(grepl(paste0("^Inference: +conformal, Jackknife\\+ for full-path ",
                                 "estimators \\(95% level; 40 never-treated ",
                                 "controls; no S\\.E\\.\\)"), out)), info = nm)
    expect_false(any(grepl("Bootstrap draws", out)))
    expect_false(any(grepl("^Conformal:", out)))
    ## two treated units: the aggregated intervals carry a printed caveat
    ag <- grep("^Aggregates:", out, value = TRUE)
    expect_length(ag, 1)
    expect_match(ag, "averages over 2 treated units have no coverage guarantee",
                 fixed = TRUE)
    if (nm == "ife") {
      em <- grep("^EM convergence:", out, value = TRUE)
      expect_length(em, 1)
      expect_match(em, "stopped by tol")
    }
    s <- summary(fit)
    expect_s3_class(s, "summary.fect")
    expect_identical(capture.output(print(s)), capture.output(print(fit)))
    expect_s3_class(suppressMessages(plot(fit)), "ggplot")
    expect_s3_class(suppressMessages(plot(fit, type = "calendar")), "ggplot")
    expect_error(plot(fit, type = "equiv"), "conformal")
  }
  ## a failed fold at alpha = 0.05 (N = 41, q = 1) gives the whole real line
  r <- .c8_run(df, fail.calls = 3, method = "ife", r = 1)
  out <- capture.output(print(r$fit))
  cl <- grep("^Conformal:", out, value = TRUE)
  expect_length(cl, 1)
  expect_match(cl, "1 failed fold(s), affected intervals widened", fixed = TRUE)
  expect_identical(capture.output(print(summary(r$fit))), out)
  expect_true(any(grepl("-Inf", out, fixed = TRUE)))
  expect_true(any(grepl("\\bInf\\b", sub("-Inf", "", out))))
})

## ---- 4.9 fect_sens --------------------------------------------------------

test_that("fect_sens() stops for a conformal fit before any cluster starts", {
  skip_on_cran()
  df <- .c8_sim(J = 20, Ntr = 1, seed = 5)
  fit <- .c8_fit(df, method = "gsynth", r = 1)
  expect_error(suppressMessages(fect_sens(fit)), "conformal")
  expect_error(suppressMessages(fect_sens(fit, parallel = TRUE)), "conformal")
})

## ---- 4.10 documentation and version facts ---------------------------------

test_that("DESCRIPTION, NEWS and ?fect describe conformal as 2.4.9", {
  skip_on_cran()
  root <- normalizePath(testthat::test_path("..", ".."))
  skip_if_not(file.exists(file.path(root, "NEWS.md")) &&
                file.exists(file.path(root, "man", "fect.Rd")),
              "source tree not available")
  desc <- readLines(file.path(root, "DESCRIPTION"))
  ## the conformal features are 2.4.9; later versions keep the 2.4.9 section
  v <- sub("^Version: ", "", grep("^Version: ", desc, value = TRUE))
  expect_true(package_version(v) >= "2.4.9")
  news <- readLines(file.path(root, "NEWS.md"))
  expect_false(any(news == "# fect 2.4.7.9000"))
  h9 <- which(news == "# fect 2.4.9"); h8 <- which(news == "# fect 2.4.8")
  expect_length(h9, 1)
  expect_length(h8, 1)
  expect_lt(h9, h8)
  heads <- grep("^# fect ", news)
  s9 <- paste(news[h9:(min(heads[heads > h9]) - 1)], collapse = "\n")
  s8 <- paste(news[h8:(min(heads[heads > h8]) - 1)], collapse = "\n")
  expect_identical(lengths(regmatches(s9, gregexpr("\n\\* New `vartype = \"conformal\"`", s9))), 1L)
  for (ph in c("finite reference completion", "observation pattern",
               "not established by the unit-level theorem", "widened",
               "n.excluded", "n.failed", "CV = FALSE")) {
    expect_true(grepl(ph, s9, fixed = TRUE), info = ph)
  }
  expect_false(grepl("conformal", s8, fixed = TRUE))
  rd <- file.path(root, "man", "fect.Rd")
  expect_length(tools::checkRd(rd), 0)
  rdt <- paste(readLines(rd), collapse = "\n")
  for (ph in c("finite reference completion", "observation pattern",
               "not established by the unit-level theorem", "widened",
               "n.excluded", "n.failed", "n.failed.folds",
               "n.intervals.affected", "CV = FALSE", "no tolerance",
               "max.iteration", "whole real line")) {
    expect_true(grepl(ph, rdt, fixed = TRUE), info = ph)
  }
  expect_false(grepl("Three values are supported", rdt, fixed = TRUE))
  expect_false(grepl("conservativ", rdt, ignore.case = TRUE))
  df <- .c8_sim(J = 20, Ntr = 1, seed = 5)
  ext <- .c8_fit(df, method = "gsynth", r = 1)$conformal$extension
  expect_match(ext, "not established by the unit-level theorem", fixed = TRUE)
  expect_false(grepl("conservative", ext, fixed = TRUE))
})
