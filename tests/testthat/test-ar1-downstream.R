## serial = "ar1": downstream consumers (WORKS / STOP matrix), argument validation, and the
## byte-identical default path (serial = "none" fits against the 2.4.8 base build 6c9c16c).
##
## WORKS rows must return a value without error; STOP rows must stop with a message that contains
## the literal text  serial = "ar1"  (never a cryptic error). Byte-identical comparisons need a
## second library with the base build in FECT_BASE_LIB; they skip when it is not set.

skip_on_cran()

lit <- 'serial = "ar1"'
P2 <- ar1_P2()
fit_se <- function(d = P2, formula = y ~ d, ...) ar1_fect(d, formula, serial = "ar1", r = 1, se = TRUE, vartype = "parametric",
                                                           nboots = 40, keep.sims = TRUE, seed = 7, ...)
render <- function(p) { grDevices::pdf(NULL); on.exit(grDevices::dev.off()); print(p); invisible(p) }

test_that("WORKS: print, plots, post-hoc estimands, imputed outcomes", {
  f <- fit_se()
  expect_gt(length(capture.output(print(f))), 3)
  for (ty in c("gap", "status", "loading.overlap", "loadings")) {
    p <- NULL
    expect_no_error(p <- plot(f, type = ty), message = ty)
    expect_false(is.null(p), info = ty)
    expect_no_error(render(p))
  }
  expect_no_error(att.cumu(f)); expect_no_error(effect(f))
  expect_no_error(estimand(f, "att")); expect_no_error(estimand(f, "att.cumu")); expect_no_error(estimand(f, "aptt"))
  expect_no_error(imputed_outcomes(f)); expect_no_error(imputed_outcomes(f, replicates = TRUE))
  ## log.att on a positive outcome. The spec adds 10 to y; on this panel that leaves Y0_b <= 0 in >= 10%
  ## of the draws for the worst cell of the unfiltered fit too (fect's own guard, any serial), so the
  ## panel is shifted by 20 here.
  dp <- P2; dp$y <- dp$y + 20
  fp <- fit_se(dp)
  expect_no_error(estimand(fp, "log.att"))
})

test_that("WORKS: plot(type = \"gap\") shows the post-treatment points with the est.att estimates", {
  f <- fit_se()
  x <- plot(f, type = "gap", return.data = TRUE)
  est <- x$data$estimate
  expect_true(all(1:5 %in% est$Period))
  post <- est[est$Period > 0, ]
  expect_close(post$ATT, unname(f$est.att[as.character(post$Period), "ATT"]), 1e-10, "plot data ATT at post event times")
})

test_that("WORKS: printing a fit twice gives identical output", {
  f <- fit_se()
  a <- tempfile(); b <- tempfile()
  sink(a); print(f); sink(); sink(b); print(f); sink()
  expect_identical(readLines(a), readLines(b))
})

test_that("STOP: features that need the levels-scale pretreatment residuals refuse with the argument in the message", {
  f <- fit_se()
  for (ty in c("equiv", "exit", "factors", "calendar", "box", "counterfactual", "sens", "heterogeneous"))
    expect_error(plot(f, type = ty), lit, fixed = TRUE, info = paste("plot", ty))
  expect_error(fect:::diagtest(f), lit, fixed = TRUE)
  expect_error(fect_sens(f), lit, fixed = TRUE)
  expect_error(fect_iden(f), lit, fixed = TRUE)
  expect_error(fect_mspe(f), lit, fixed = TRUE)
  expect_error(estimand(f, "att", test = "placebo"), lit, fixed = TRUE)
  expect_error(estimand(f, "att", test = "carryover"), lit, fixed = TRUE)
  for (a in list(list(placeboTest = TRUE), list(carryoverTest = TRUE), list(loo = TRUE), list(dloo = TRUE),
                 list(permute = TRUE), list(cm = TRUE)))
    expect_error(do.call(ar1_fect, c(list(P2, serial = "ar1", r = 1), a)), lit, fixed = TRUE, info = names(a))
})

test_that("plot(type = \"cumul\") returns a plot or stops naming the argument", {
  f <- fit_se()
  out <- tryCatch({ render(plot(f, type = "cumul")); "plot" }, error = function(e) conditionMessage(e))
  expect_true(identical(out, "plot") || grepl(lit, out, fixed = TRUE), info = out)
})

## Section 7: argument validation ---------------------------------------------------------------------------------------------
test_that("refusals under serial = \"ar1\" name the argument", {
  d <- P2; d$w <- 1; d$g <- rep(1:5, length.out = nrow(d))
  chk <- function(...) expect_error(ar1_fect(d, serial = "ar1", r = 1, ...), lit, fixed = TRUE)
  chk(method = "ife"); chk(method = "fe"); chk(method = "mc"); chk(method = "cfe")
  chk(W = "w"); chk(group = "g"); chk(cfe = list(c("g"))); chk(group.fe = "g")
  chk(loading.bound = "simplex"); chk(balance.period = c(1, 15))
  for (cm in c("block", "all_units", "treated_units", "loo")) chk(cv.method = cm, CV = TRUE, r = c(0, 3))
  chk(cv.rule = "1se", CV = TRUE, r = c(0, 3)); chk(criterion = "pc", CV = TRUE, r = c(0, 3))
  chk(se = TRUE, vartype = "jackknife")
  expect_error(ar1_fect(d, serial = "bogus", r = 1))
})

test_that("accepted: rolling (default), pretreatment, normalize, every force, r single / range / NULL", {
  expect_no_error(f <- ar1_fect(P2, serial = "ar1", r = c(0, 3), cv.method = "rolling", quiet = FALSE) |> suppressMessages())
  m <- testthat::capture_messages(ar1_fect(P2, serial = "ar1", r = c(0, 3), cv.method = "rolling", quiet = FALSE))
  expect_gte(sum(grepl(lit, m, fixed = TRUE)), 1)                         # one message about the rank rule used
  expect_no_error(ar1_fect(P2, serial = "ar1", r = c(0, 3), cv.method = "pretreatment"))
  expect_no_error(ar1_fect(P2, serial = "ar1", r = 1, normalize = TRUE))
  for (fo in c("none", "unit", "time", "two-way")) {
    f <- ar1_fect(P2, serial = "ar1", r = 1, force = fo)
    expect_true(is.finite(f$att.avg), info = fo)
    expect_identical(f$force, c("none" = 0, "unit" = 1, "time" = 2, "two-way" = 3)[[fo]], info = fo)
  }
  expect_no_error(ar1_fect(P2, serial = "ar1", r = 2))
  expect_no_error(ar1_fect(P2, serial = "ar1", r = c(1, 3)))
  expect_no_error(ar1_fect(P2, serial = "ar1"))
})

test_that("the three signatures carry the same default for serial", {
  d <- function(f) eval(formals(f)$serial)
  expect_identical(d(fect::fect), c("none", "ar1"))
  expect_identical(d(getS3method("fect", "formula")), d(fect::fect))
  expect_identical(d(getS3method("fect", "default")), d(fect::fect))
})

## Section 8: serial = "none" fits are untouched ---------------------------------------------------------------------------
test_that("serial = \"none\" fits carry no ar1 state", {
  f <- ar1_fect(P2, CV = FALSE, r = 2)
  expect_true(is.null(f$ar1)); expect_false(identical(f$serial, "ar1"))
})

test_that("byte-identical default path: five fits equal the 6c9c16c base", {
  base_lib <- ar1_base_lib()
  skip_if(!nzchar(base_lib), "set FECT_BASE_LIB to a library with fect 6c9c16c")
  code <- '
    source(Sys.getenv("AR1_HELPER"))
    P2 <- ar1_P2()
    out <- list(
      gs_rolling = ar1_fect(P2, CV = TRUE, r = c(0, 3), cv.method = "rolling", seed = 1),
      gs_loo     = ar1_fect(P2, CV = TRUE, r = c(0, 3), cv.method = "loo", seed = 1),
      ife_tu     = ar1_fect(P2, method = "ife", time.component.from = "nevertreated", CV = TRUE, r = c(0, 3),
                            cv.method = "treated_units", seed = 1),
      fe_para    = ar1_fect(P2, method = "fe", CV = FALSE, se = TRUE, vartype = "parametric", nboots = 30, seed = 1,
                            time.component.from = "nevertreated"),
      mc         = ar1_fect(P2, method = "mc", CV = TRUE, seed = 1, lambda = c(1, 0.1)))
    out <- lapply(out, function(f) f[c("att.avg", "r.cv", "est.att", "est.avg", "CV.out", "eff")])
  '
  withr::local_envvar(AR1_HELPER = normalizePath(test_path("helper-ar1.R")))
  new <- local({ out <- NULL; eval(parse(text = code), envir = environment()); out })
  old <- ar1_run_in_lib(base_lib, code)$out
  expect_identical(names(new), names(old))
  for (nm in names(old)) {
    expect_equal(new[[nm]]$att.avg, old[[nm]]$att.avg, tolerance = 1e-12, info = nm)
    expect_identical(new[[nm]]$r.cv, old[[nm]]$r.cv, info = nm)
    expect_equal(new[[nm]]$est.att, old[[nm]]$est.att, tolerance = 1e-12, info = nm)
    expect_equal(new[[nm]]$est.avg, old[[nm]]$est.avg, tolerance = 1e-12, info = nm)
    expect_equal(new[[nm]]$CV.out, old[[nm]]$CV.out, tolerance = 1e-12, info = nm)
    expect_equal(new[[nm]]$eff, old[[nm]]$eff, tolerance = 1e-12, info = nm)
  }
})
