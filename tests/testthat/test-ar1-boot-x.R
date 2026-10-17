## Phase 3c: the parametric bootstrap under serial = "ar1" re-estimates the filter coefficient x in every draw
## (the rank is held at the fit's chosen rank; x = 0 is allowed). test-spec.md "Amendments 3" (T1-T9, R-a, R-d, R-e).
##
## Test hooks used here (set with options(); restored by withr): fect.ar1.boot = "fixed" is the bootstrap of the build
## before 3c (draws at the chosen pair, mapping with the chosen x); fect.ar1.x fixes x everywhere, draws included.
## T4 (the existing test files) is a run, not a test file. T5 uses the seam fect:::.fect_ar1_search_x. The cost study
## (T11), the draw-spread study (T10) and the spot check of three stored simulation cells (R-b) are scripts, not package
## tests (they need minutes of compute or the simulation folder).
##
## Departures from the literal text of test-spec Amendments 3, each reported in audit.md Round 4:
##  * T2 compares the draws of the effect with 0 (the effect in the parametric draws), not with fit$att.avg: the draws are
##    generated under no effect, so att.avg.boot centres on 0 in the old bootstrap too (SE 2.5e-16 there as well).
##  * T7/T8 "att.avg.boot is NA where x.boot is NA": fect_boot removes failed draws from att.avg.boot, att.boot and
##    eff.boot (as it always has); x.boot keeps NA in draw order. The tests assert the equivalent statement.
##  * T3(c) exact equality across parallel and serial mode is not required (different random-number streams, also on
##    serial = "none"); the test states what the spec states.

skip_on_cran()

ref_slots <- c("att.avg", "est.avg", "est.avg.unit", "est.att", "est.att90", "att.bound", "att.boot", "att.avg.boot",
               "att.avg.unit.boot", "att.count.boot", "att.vcov", "eff.boot", "D.boot", "I.boot", "colnames.boot", "eff", "att")
b3c <- list(se = TRUE, vartype = "parametric", nboots = 40, keep.sims = TRUE, seed = 7, tol = 1e-10)

## One fect(serial = "ar1") call, quiet; boot = "x" (default) or "fixed" (the build before 3c). `opts` are extra options().
## Returns the fit, the warnings it raised and .Random.seed right after the call.
run3c <- function(d, formula = y ~ d, boot = "x", opts = list(), ..., r = 2) {
  args <- c(list(d, formula, serial = "ar1", r = r, CV = FALSE), b3c)
  dots <- list(...); args[names(dots)] <- dots
  o <- c(list(fect.ar1.boot = if (identical(boot, "x")) NULL else boot), opts)      # NULL unsets: "x" is the default
  res <- withr::with_options(o, ar1_collect(do.call(ar1_fect, args)))
  list(fit = res$value, warnings = res$warnings, seed = get0(".Random.seed", envir = globalenv(), inherits = FALSE))
}
P2 <- ar1_P2(); P3 <- ar1_P3(); P4 <- ar1_P4(); P5 <- ar1_P5()
P5g <- P5; P5g$x[P5g$id == 32 & P5g$time %in% c(15, 16)] <- NA          # the E8 gap: unit 32 cannot be rebuilt at x != 0
panels <- list(P2 = list(d = P2, f = y ~ d), P3 = list(d = P3, f = y ~ d), P4 = list(d = P4, f = y ~ d), P5 = list(d = P5, f = y ~ d + x),
               P5g = list(d = P5g, f = y ~ d + x))
memo <- new.env()
fit_memo <- function(key, ...) { if (is.null(memo[[key]])) memo[[key]] <- run3c(...); memo[[key]] }
## Known warnings: the E8 gap warns at x != 0 (both bootstraps alike).
e8_msg <- "could not be rebuilt on the outcome scale"
only_known_warnings <- function(w) all(grepl(e8_msg, w, fixed = TRUE))

## ---- T1 identity with the old bootstrap when x is fixed -------------------------------------------------------------
tolT1 <- function(panel, xval) if (xval == 0) 1e-12 else if (panel %in% c("P2", "P3")) 1e-9 else 1e-6
check_t1 <- function(label, panel, xval, ...) {
  p <- panels[[panel]]; tol <- tolT1(panel, xval)
  a <- run3c(p$d, p$f, "x", list(fect.ar1.x = xval), ...); b <- run3c(p$d, p$f, "fixed", list(fect.ar1.x = xval), ...)
  worst <- 0
  for (s in ref_slots) {
    e <- ar1_slot_diff(a$fit[[s]], b$fit[[s]]); worst <- max(worst, e)
    expect_true(e <= tol, info = sprintf("%s: slot %s differs by %.3e (tolerance %.0e)", label, s, e, tol))
  }
  ar1_log(paste("T1", label), "max over 17 slots", worst, tol)
  expect_identical(a$seed, b$seed, info = paste(label, ": .Random.seed after the call"))
  expect_identical(a$fit$ar1$boot, "x"); expect_identical(b$fit$ar1$boot, "fixed")
  expect_length(a$fit$ar1$x.boot, 40)
  expect_true(all(a$fit$ar1$x.boot[!is.na(a$fit$ar1$x.boot)] == a$fit$ar1$x), info = paste(label, ": x.boot equals the held x exactly"))
  expect_equal(a$fit$ar1$x, xval)
  expect_true(only_known_warnings(a$warnings) && only_known_warnings(b$warnings), info = paste(label, ":", paste(c(a$warnings, b$warnings), collapse = " | ")))
  expect_identical(sort(unique(a$warnings)), sort(unique(b$warnings)), info = paste(label, ": the two bootstraps warn alike"))
  invisible(list(a = a, b = b))
}

test_that("T1: with x fixed the new bootstrap equals the old one (P2, P3, P4, P5 with the E8 gap; x = 0, 0.5, 1)", {
  for (pn in c("P2", "P3", "P4", "P5g")) for (xv in c(0, 0.5, 1)) check_t1(sprintf("%s x=%g", pn, xv), pn, xv)
})

test_that("T1 (additional cells): ci.method = 'basic', normalize = TRUE, para.error = 'ar' and 'auto'", {
  check_t1("P2 x=0.5 ci.method=basic", "P2", 0.5, ci.method = "basic")
  check_t1("P5 x=0.5 normalize=TRUE", "P5", 0.5, normalize = TRUE)
  check_t1("P4 x=0.5 para.error=ar", "P4", 0.5, para.error = "ar")
  check_t1("P4 x=0.5 para.error=auto", "P4", 0.5, para.error = "auto")
})

## ---- T2 noiseless panel ---------------------------------------------------------------------------------------------------
test_that("T2: noiseless panel -- every draw's x equals x-hat, the interval equals the old one", {
  d <- ar1_noiseless_3c()
  args <- list(r = 1, para.error = "empirical", tol = 1e-12)
  a <- do.call(run3c, c(list(d, y ~ d, "x"), args)); b <- do.call(run3c, c(list(d, y ~ d, "fixed"), args))
  f <- a$fit; g <- b$fit
  ## preconditions (recorded; the test stops here with a clear message if the generator is not noiseless)
  expect_lt(abs(f$ar1$x - 0.5), 1e-3, label = "precondition |x-hat - 0.5| < 1e-3")
  expect_lt(max(abs(f$res), na.rm = TRUE), 1e-6, label = "precondition max|res| < 1e-6")
  if (abs(f$ar1$x - 0.5) >= 1e-3 || max(abs(f$res), na.rm = TRUE) >= 1e-6) stop("T2 preconditions failed: generator not noiseless; report to the leader")
  xb <- f$ar1$x.boot; ok <- !is.na(xb)
  expect_equal(f$ar1$boot.failed, 0L); expect_true(all(ok))
  e_x <- max(abs(xb[ok] - f$ar1$x)); ar1_log("T2 max |x.boot - x-hat|", sprintf("%.3e", e_x), e_x, 3e-4)
  expect_lte(e_x, 3e-4)
  expect_lt(abs(f$att.avg - 1), 1e-3)
  ## the draws of the parametric bootstrap carry no effect: the estimate in every draw is the draw's own truth, 0
  ## (departure from the literal text, which compares with att.avg; the old bootstrap's draws are 0 as well)
  e_d <- max(abs(f$att.avg.boot)); ar1_log("T2 max |att.avg.boot - 0|", sprintf("%.3e", e_d), e_d, 5e-3)
  expect_lte(e_d, 5e-3)
  expect_lte(max(abs(g$att.avg.boot)), 5e-3)
  se_new <- f$est.avg[, "S.E."]; se_old <- g$est.avg[, "S.E."]
  ar1_log("T2 S.E. new", sprintf("%.3e", se_new), se_new, 5e-3); ar1_log("T2 S.E. fixed", sprintf("%.3e", se_old), se_old, 1e-6)
  expect_lt(se_new, 5e-3); expect_lt(se_old, 1e-6)
  e_i <- max(abs(f$est.avg[, c("CI.lower", "CI.upper")] - g$est.avg[, c("CI.lower", "CI.upper")]))
  ar1_log("T2 interval ends, new vs fixed", sprintf("%.3e", e_i), e_i, 1e-2)
  expect_lt(e_i, 1e-2)
  expect_true(only_known_warnings(c(a$warnings, b$warnings)))
})

## ---- T3 reproducibility ----------------------------------------------------------------------------------------------------
test_that("T3(a): same seed, serial mode: identical in every slot; another seed changes x.boot", {
  a <- run3c(P2); b <- run3c(P2)
  expect_identical(ar1_fit_identical(a$fit, b$fit), character(0))
  expect_identical(a$fit$ar1$x.boot, b$fit$ar1$x.boot)
  expect_identical(a$seed, b$seed)
  d <- run3c(P2, seed = 8)
  expect_false(isTRUE(all.equal(a$fit$ar1$x.boot, d$fit$ar1$x.boot, tolerance = 0)))
  expect_identical(a$fit$ar1$x, d$fit$ar1$x)                                   # the point estimate does not depend on the seed
})

skip_par <- function() {
  skip_if(ar1_is_dev_load(), "parallel = TRUE needs an installed fect (workers load the installed package)")
}
test_that("T3(b): same seed, parallel = TRUE with 2 workers: two calls identical in every slot", {
  skip_par()
  a <- run3c(P2, parallel = TRUE, cores = 2); b <- run3c(P2, parallel = TRUE, cores = 2)
  expect_identical(ar1_fit_identical(a$fit, b$fit), character(0))
  expect_identical(a$fit$ar1$x.boot, b$fit$ar1$x.boot)
  expect_gt(length(unique(round(a$fit$ar1$x.boot, 10))), 5)                    # the workers ran the new draws
  ## 2 against 3 workers: required for ar1 only if serial = "none" already agrees across 2 and 3 workers
  n2 <- do.call(ar1_fect, c(list(P2, serial = "none", r = 2, CV = FALSE, parallel = TRUE, cores = 2), b3c))
  n3 <- do.call(ar1_fect, c(list(P2, serial = "none", r = 2, CV = FALSE, parallel = TRUE, cores = 3), b3c))
  base_same <- length(ar1_fit_identical(n2, n3)) == 0
  ar1_note("T3(b): serial = 'none', 2 against 3 workers identical in every slot: ", base_same)
  c3 <- run3c(P2, parallel = TRUE, cores = 3)
  same23 <- length(ar1_fit_identical(a$fit, c3$fit)) == 0
  ar1_note("T3(b): serial = 'ar1' default, 2 against 3 workers identical in every slot: ", same23)
  if (base_same) expect_true(same23, info = "the base case agrees across 2 and 3 workers, so the ar1 default must as well")
})

test_that("T3(c): parallel against serial (200 draws): same point estimates; SE ratio, x.boot mean, failures in range", {
  skip_par()
  for (pn in c("P2", "P5")) {
    p <- panels[[pn]]
    ser <- run3c(p$d, p$f, nboots = 200)$fit
    par <- run3c(p$d, p$f, nboots = 200, parallel = TRUE, cores = 2)$fit
    expect_close(par$att.avg, ser$att.avg, 1e-12, paste(pn, "att.avg"))
    expect_close(par$ar1$x, ser$ar1$x, 1e-12, paste(pn, "x-hat"))
    expect_close(par$r.cv, ser$r.cv, 1e-12, paste(pn, "r.cv"))
    ratio <- par$est.avg[, "S.E."] / ser$est.avg[, "S.E."]
    dm <- abs(mean(par$ar1$x.boot, na.rm = TRUE) - mean(ser$ar1$x.boot, na.rm = TRUE)); sdm <- sd(ser$ar1$x.boot, na.rm = TRUE)
    ar1_note(sprintf("T3(c) %s: SE ratio par/ser %.4f; |mean x.boot diff| %.4f vs 0.5 sd %.4f; failed par %d ser %d", pn, ratio, dm, 0.5 * sdm,
                     par$ar1$boot.failed, ser$ar1$boot.failed))
    expect_gte(ratio, 0.7); expect_lte(ratio, 1.43)
    expect_lte(dm, 0.5 * sdm)
    expect_lte(abs(par$ar1$boot.failed - ser$ar1$boot.failed), 5)
    ## baseline: the same two numbers for serial = "none"
    nser <- do.call(ar1_fect, c(list(p$d, p$f, serial = "none", r = 2, CV = FALSE), modifyList(b3c, list(nboots = 200))))
    npar <- do.call(ar1_fect, c(list(p$d, p$f, serial = "none", r = 2, CV = FALSE, parallel = TRUE, cores = 2), modifyList(b3c, list(nboots = 200))))
    ar1_note(sprintf("T3(c) %s baseline serial = 'none': SE ratio par/ser %.4f", pn, npar$est.avg[, "S.E."] / nser$est.avg[, "S.E."]))
  }
})

## ---- T5 the seam ---------------------------------------------------------------------------------------------------------------
test_that("T5: fect:::.fect_ar1_search_x reproduces the main fit's x and profile (P2, P4, P5); boundary behaviour", {
  for (pn in c("P2", "P4", "P5")) {
    p <- panels[[pn]]
    fit <- ar1_fect(p$d, p$f, serial = "ar1", r = 2, CV = FALSE, se = FALSE, normalize = FALSE, tol = 1e-10)
    if (nrow(fit$ar1$dropped) > 0) { message("T5: units dropped on ", pn, "; skipped"); next }
    A <- ar1_levels_arrays(p$d, covar = pn == "P5")
    res <- fect:::.fect_ar1_search_x(A$Y, A$D, A$X, r = 2, force = 3L, tol = 1e-10)
    expect_true(is.finite(res$x) && res$x >= 0 && res$x <= 1)
    expect_close(res$x, fit$ar1$x, 1e-6, paste(pn, "x"))
    ## the 21 grid points seq(0, 1, by = 0.05) (the profile also lists the refinement evaluations)
    grid <- round(seq(0, 1, by = 0.05), 8)
    m1 <- match(grid, round(res$profile$x, 8)); m2 <- match(grid, round(fit$ar1$profile$x, 8))
    expect_false(anyNA(m1)); expect_false(anyNA(m2))
    expect_close(res$profile$ssr[m1], fit$ar1$profile$ssr[m2], 1e-8, paste(pn, "profile ssr at the 21 grid points"), rel = TRUE)
    ar1_note(sprintf("T5 %s: x seam %.8f main %.8f; profile rows seam %d main %d", pn, res$x, fit$ar1$x, nrow(res$profile), nrow(fit$ar1$profile)))
  }
  wn <- ar1_levels_arrays(ar1_panel(30, 5, 15, 5, rho = 0, seed = 5001))
  expect_lte(fect:::.fect_ar1_search_x(wn$Y, wn$D, NULL, r = 2, force = 3L, tol = 1e-10)$x, 0.15)
  rw <- ar1_levels_arrays(ar1_rw_panel(30, 5, 15, 5, 6001))
  expect_gte(fect:::.fect_ar1_search_x(rw$Y, rw$D, NULL, r = 2, force = 3L, tol = 1e-10)$x, 0.7)
})

## ---- T6 invariances of the new bootstrap ----------------------------------------------------------------------------------------
test_that("T6: scale, shift, unit constants (P2, P4) and normalize (P5) leave x.boot, the S.E. and the estimate alone", {
  for (pn in c("P2", "P4")) {
    p <- panels[[pn]]; base <- fit_memo(paste0("t6_", pn), p$d, p$f)$fit
    set.seed(5); cst <- rnorm(length(unique(p$d$id)), sd = 3)
    sc <- local({ d <- p$d; d$y <- d$y * 1000; run3c(d, p$f)$fit })
    sh <- local({ d <- p$d; d$y <- d$y + 100; run3c(d, p$f)$fit })
    un <- local({ d <- p$d; d$y <- d$y + cst[match(d$id, sort(unique(d$id)))]; run3c(d, p$f)$fit })
    expect_close(sc$att.avg, 1000 * base$att.avg, 1e-6, paste(pn, "(a) att.avg x1000"), rel = TRUE)
    expect_close(sc$est.avg[, "S.E."], 1000 * base$est.avg[, "S.E."], 1e-6, paste(pn, "(a) S.E. x1000"), rel = TRUE)
    expect_close(sc$est.avg[, c("CI.lower", "CI.upper")], 1000 * base$est.avg[, c("CI.lower", "CI.upper")], 1e-6, paste(pn, "(a) CI x1000"), rel = TRUE)
    expect_close(sc$ar1$x.boot, base$ar1$x.boot, 1e-6, paste(pn, "(a) x.boot"), rel = TRUE)
    expect_close(sc$ar1$x, base$ar1$x, 1e-6, paste(pn, "(a) x"), rel = TRUE)
    for (nm in c("sh", "un")) {
      z <- get(nm); lb <- if (nm == "sh") "(b) +100" else "(c) unit constants"
      expect_close(z$ar1$x.boot, base$ar1$x.boot, 1e-6, paste(pn, lb, "x.boot"), rel = TRUE)
      expect_close(z$est.avg[, "S.E."], base$est.avg[, "S.E."], 1e-6, paste(pn, lb, "S.E."), rel = TRUE)
      expect_close(z$att.avg, base$att.avg, 1e-6, paste(pn, lb, "att.avg"), rel = TRUE)
    }
  }
  nt <- run3c(P5, y ~ d + x, normalize = TRUE)$fit; nf <- run3c(P5, y ~ d + x, normalize = FALSE)$fit
  expect_close(nt$est.avg[, "S.E."], nf$est.avg[, "S.E."], 1e-5, "(d) normalize: S.E.", rel = TRUE)
  expect_close(nt$ar1$x.boot, nf$ar1$x.boot, 1e-5, "(d) normalize: x.boot", rel = TRUE)
})

## ---- T7 contract of the new fields -----------------------------------------------------------------------------------------------
test_that("T7: boot, x.boot, boot.failed and options on the four panels; se = FALSE; keep.sims; bad option; ci.method", {
  for (pn in c("P2", "P3", "P4", "P5")) {
    p <- panels[[pn]]; f <- fit_memo(paste0("t7_", pn), p$d, p$f)$fit
    expect_identical(f$ar1$boot, "x")
    xb <- f$ar1$x.boot
    expect_length(xb, 40); expect_true(all(xb[!is.na(xb)] >= 0 & xb[!is.na(xb)] <= 1))
    ## failed draws are removed from att.avg.boot (as fect_boot always does) and are NA in x.boot
    expect_equal(sum(is.na(xb)), f$ar1$boot.failed); expect_equal(length(f$att.avg.boot), 40 - f$ar1$boot.failed)
    expect_false(anyNA(f$att.avg.boot))
    expect_identical(f$ar1$boot.failed, as.integer(sum(is.na(xb))))
    expect_identical(f$ar1$options$boot, "x")
    expect_true(is.na(f$ar1$options$x) && is.na(f$ar1$options$rounds))
    expect_true(isTRUE(f$ar1$options$exclude.nonconverged)); expect_false(isTRUE(f$ar1$options$mask.zero))
    expect_gt(length(unique(xb)), 1)
    ## the old bootstrap
    g <- run3c(p$d, p$f, "fixed")$fit
    expect_identical(g$ar1$boot, "fixed"); expect_null(g$ar1$x.boot); expect_identical(g$ar1$options$boot, "fixed")
  }
  ## options recorded
  o <- run3c(P2, opts = list(fect.ar1.x = 0.5, fect.ar1.rounds = 3, fect.ar1.exclude.nonconverged = FALSE, fect.ar1.mask.zero = TRUE))$fit$ar1$options
  expect_equal(o$x, 0.5); expect_equal(o$rounds, 3); expect_false(o$exclude.nonconverged); expect_true(o$mask.zero); expect_identical(o$boot, "x")
  ## se = FALSE
  n <- run3c(P2, se = FALSE)$fit
  expect_null(n$ar1$x.boot); expect_identical(n$ar1$boot.failed, NA_integer_)
  ## keep.sims = FALSE
  k <- run3c(P2, keep.sims = FALSE)$fit
  expect_length(k$ar1$x.boot, 40); expect_null(k$eff.boot)
  expect_equal(k$ar1$x.boot, fit_memo("t7_P2", P2, y ~ d)$fit$ar1$x.boot)
  ## a bad value
  expect_error(run3c(P2, boot = "bad"), "fect.ar1.boot", fixed = TRUE)
  ## ci.method
  for (cm in c("basic", "normal")) {
    f <- run3c(P2, ci.method = cm)$fit
    expect_close(f$est.avg[, "S.E."], sd(f$att.avg.boot, na.rm = TRUE), 1e-12, paste("ci.method", cm, "S.E. equals sd(att.avg.boot)"))
    expect_true(all(is.finite(f$est.avg)))
  }
})

## T7b (added by the tester from the decision's "the search must allow x = 0"): x = 0 is reachable inside the draws.
test_that("T7b: on white-noise errors x-hat is 0 and most draws land on exactly 0; no draw fails", {
  d <- ar1_panel(30, 5, 15, 5, rho = 0, seed = 5001)
  f <- run3c(d, nboots = 60, keep.sims = FALSE)$fit
  xb <- f$ar1$x.boot
  expect_equal(f$ar1$x, 0); expect_equal(f$ar1$boot.failed, 0L)
  expect_true(all(xb >= 0 & xb <= 1)); expect_gt(sum(xb == 0), 0)
  ar1_note(sprintf("T7b: seed 5001, share of draws at exactly 0: %.2f; mean %.4f; max %.4f", mean(xb == 0), mean(xb), max(xb)))
})

## ---- T8 failure handling (mock) ---------------------------------------------------------------------------------------------------
test_that("T8: draws whose search throws are dropped, counted, warned about once; all fail -> error; 1 of 30 -> no warning", {
  orig <- fect:::.fect_ar1_search_x
  calls <- new.env(); calls$n <- 0L
  n_main <- 0L                                                                 # search calls made by the main fit (set below)
  nd <- function() calls$n - n_main                                           # index of the current draw-level call
  count_only <- function(...) { calls$n <- calls$n + 1L; orig(...) }
  every3 <- function(...) { calls$n <- calls$n + 1L; if (nd() > 0L && nd() %% 3L == 0L) stop("injected failure"); orig(...) }
  once <- function(...) { calls$n <- calls$n + 1L; if (nd() == 7L) stop("injected failure"); orig(...) }
  always <- function(...) { calls$n <- calls$n + 1L; if (nd() > 0L) stop("injected failure"); orig(...) }
  a30 <- list(nboots = 30)
  ## does the main fit use the function? (the spec counts only the draw calls)
  testthat::local_mocked_bindings(.fect_ar1_search_x = count_only, .package = "fect")
  invisible(run3c(P2, se = FALSE)); n_main <- calls$n
  calls$n <- 0L; clean <- run3c(P2, nboots = 30)
  n_draw <- calls$n - n_main
  expect_equal(n_draw, 30L)                                                    # one search call per draw
  ## every third draw call throws
  testthat::local_mocked_bindings(.fect_ar1_search_x = every3, .package = "fect")
  calls$n <- 0L
  r3 <- run3c(P2, nboots = 30); f <- r3$fit
  expect_equal(f$ar1$boot.failed, 10L)
  expect_equal(which(is.na(f$ar1$x.boot)), seq(3L, 30L, by = 3L))
  expect_false(anyNA(f$att.avg.boot)); expect_equal(length(f$att.avg.boot), 20L)    # fect_boot removes failed draws
  expect_equal(dim(f$eff.boot)[3], 20L)
  expect_close(f$est.avg[, "S.E."], sd(f$att.avg.boot), 1e-12, "S.E. equals sd of the surviving draws")
  ## the surviving draws are the clean run's draws 1, 2, 4, 5, ... (the hook draws no random numbers)
  keep <- !is.na(f$ar1$x.boot)
  expect_close(f$ar1$x.boot[keep], clean$fit$ar1$x.boot[keep], 1e-12, "x.boot of the surviving draws equals the clean run's")
  expect_close(as.vector(f$att.avg.boot), as.vector(clean$fit$att.avg.boot)[keep], 1e-12, "att.avg.boot of the surviving draws equals the clean run's")
  expect_length(r3$warnings, 1L)
  expect_match(r3$warnings, 'serial = "ar1"', fixed = TRUE); expect_match(r3$warnings, "bootstrap draws failed", fixed = TRUE)
  ## one of 30
  testthat::local_mocked_bindings(.fect_ar1_search_x = once, .package = "fect")
  calls$n <- 0L
  r1 <- run3c(P2, nboots = 30)
  expect_equal(r1$fit$ar1$boot.failed, 1L); expect_length(r1$warnings, 0L)
  ## all of them
  testthat::local_mocked_bindings(.fect_ar1_search_x = always, .package = "fect")
  w <- character(0)
  expect_error(withCallingHandlers(run3c(P2, nboots = 30), warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") }),
               'serial = "ar1": all bootstrap draws failed', fixed = TRUE)
  ar1_note("T8: warnings raised before the all-failed error: ", paste(w, collapse = " | "))
})

## ---- T9 the old bootstrap is reachable and unchanged ------------------------------------------------------------------------------
test_that("T9: fect.ar1.boot = 'fixed' equals the build before 3c (3600f95) with x searched; serial = 'none' unchanged", {
  lib3a <- Sys.getenv("FECT_3A_LIB", "")
  skip_if(!nzchar(lib3a) || !dir.exists(file.path(lib3a, "fect")), "set FECT_3A_LIB to a library with fect 3600f95")
  helper <- normalizePath(test_path("helper-ar1.R"))
  code <- sprintf('
    source(%s)
    P2 <- ar1_P2(); P3 <- ar1_P3(); P4 <- ar1_P4(); P5 <- ar1_P5()
    pan <- list(P2 = list(P2, y ~ d, 2), P3 = list(P3, y ~ d, 2), P4 = list(P4, y ~ d, 2), P5 = list(P5, y ~ d + x, 2),
                P2r = list(P2, y ~ d, c(0, 5)), P5r = list(P5, y ~ d + x, c(0, 5)))
    slots <- c("att.avg","est.avg","est.avg.unit","est.att","est.att90","att.bound","att.boot","att.avg.boot","att.avg.unit.boot",
               "att.count.boot","att.vcov","eff.boot","D.boot","I.boot","colnames.boot","eff","att")
    out <- lapply(pan, function(p) {
      f <- ar1_fect(p[[1]], p[[2]], serial = "ar1", r = p[[3]], se = TRUE, vartype = "parametric", nboots = 40, keep.sims = TRUE,
                    seed = 7, tol = 1e-10)
      c(f[slots], list(x = f$ar1$x, path = f$ar1$path, r.cv = f$r.cv))
    })
    nn <- ar1_fect(P2, y ~ d, serial = "none", r = 2, CV = FALSE, se = TRUE, vartype = "parametric", nboots = 40, keep.sims = TRUE, seed = 7, tol = 1e-10)
    out$none <- nn[slots]
    out', deparse(helper))
  old <- ar1_run_in_lib(lib3a, code)$out
  cmp <- function(f, o, label) {
    worst <- 0
    for (s in c(ref_slots, "x", "path", "r.cv")) {
      a <- if (s == "x") f$ar1$x else if (s == "path") f$ar1$path[, c("round", "x", "r", "att")] else if (s == "r.cv") f$r.cv else f[[s]]
      b <- o[[s]]; if (s == "path") b <- b[, c("round", "x", "r", "att")]
      e <- ar1_slot_diff(a, b); worst <- max(worst, e)
      expect_true(e <= 1e-12, info = sprintf("%s: slot %s differs by %.3e", label, s, e))
    }
    ar1_log(paste("T9", label), "max over 20 slots", worst, 1e-12)
  }
  for (nm in names(old)) {
    if (nm == "none") next
    p <- switch(sub("r$", "", nm), P2 = list(P2, y ~ d), P3 = list(P3, y ~ d), P4 = list(P4, y ~ d), P5 = list(P5, y ~ d + x))
    rr <- if (grepl("r$", nm)) c(0, 5) else 2
    x <- do.call(run3c, c(list(p[[1]], p[[2]], "x"), list(r = rr, CV = if (length(rr) > 1) TRUE else FALSE)))
    fx <- do.call(run3c, c(list(p[[1]], p[[2]], "fixed"), list(r = rr, CV = if (length(rr) > 1) TRUE else FALSE)))
    cmp(fx$fit, old[[nm]], paste(nm, "fixed against 3600f95"))
    expect_identical(x$seed, fx$seed, info = paste(nm, ": .Random.seed after the x and fixed calls"))
    expect_equal(x$fit$ar1$x, fx$fit$ar1$x); expect_identical(x$fit$ar1$boot, "x"); expect_false(x$fit$ar1$x.fixed)
  }
  nn <- ar1_fect(P2, y ~ d, serial = "none", r = 2, CV = FALSE, se = TRUE, vartype = "parametric", nboots = 40, keep.sims = TRUE, seed = 7, tol = 1e-10)
  worst <- 0
  for (s in ref_slots[ref_slots %in% names(old$none)]) {
    e <- ar1_slot_diff(nn[[s]], old$none[[s]]); worst <- max(worst, e)
    expect_true(e <= 1e-12, info = sprintf("serial = 'none': slot %s differs by %.3e", s, e))
  }
  ar1_log("T9 serial = none against 3600f95", paste(sum(ref_slots %in% names(old$none)), "slots compared"), worst, 1e-12)
})

## ---- R-a: the x = 0 covariate fix (57bb66f), the reproduction in the round-3 BLOCK ------------------------------------------------------
test_that("R-a: at x = 0 a covariate gap at (32, 16) costs only its own cell; att.avg equals gsynth on the panel without period 1", {
  d <- ar1_panel(30, 5, 15, 5, rho = .6, seed = 104, cov = TRUE)
  d <- d[!(d$id == 32 & d$time == 16), ]
  f <- withr::with_options(list(fect.ar1.x = 0), ar1_fect(d, y ~ d + x, serial = "ar1", r = 2))
  g <- ar1_fit_levels(d, 2, y ~ d + x)
  expect_close(f$att.avg, g$att.avg, 1e-8, "att.avg equals gsynth (r = 2, CV = FALSE, min.T0 = 4) on the data without period 1")
  expect_lt(abs(f$att.avg - 0.030904), 1e-6)                                    # the audit's reference value, quoted to 6 digits
  e <- f$eff[match(17, f$rawtime), match("32", as.character(f$id))]
  expect_true(is.finite(e))
  expect_close(e, g$eff[match(17, g$rawtime), match("32", as.character(g$id))], 1e-8, "unit 32 at t = 17 equals the levels fit")
})

## ---- R-d: the #178 gate (did before 2.5.1 stops) ---------------------------------------------------------------------------------------
test_that("R-d: did_wrapper(method = 'cs_never') stops with a message when did < 2.5.1 and passes the gate at 2.5.1", {
  msg <- "needs did >= 2.5.1"
  d <- ar1_P2()[, c("id", "time", "y", "d")]
  call <- function() fect::did_wrapper(d, Y = "y", D = "d", index = c("id", "time"), method = "cs_never", parallel = FALSE)
  ## the gate reads requireNamespace("did") and utils::packageVersion("did"); both are mocked so the test needs no did
  orig_rn <- base::requireNamespace
  testthat::local_mocked_bindings(requireNamespace = function(package, ..., quietly = FALSE) {
    if (identical(package, "did")) TRUE else orig_rn(package, ..., quietly = quietly)
  }, .package = "base")
  testthat::local_mocked_bindings(packageVersion = function(pkg, ...) package_version("2.5.0"), .package = "utils")
  expect_error(call(), msg, fixed = TRUE)
  testthat::local_mocked_bindings(packageVersion = function(pkg, ...) package_version("2.5.1"), .package = "utils")
  out <- tryCatch({ suppressWarnings(suppressMessages(call())); "ran" }, error = function(e) conditionMessage(e))
  expect_false(grepl(msg, out, fixed = TRUE), info = out)
  ar1_note("R-d: with did 2.5.1 the gate is passed; the call then ends with: ", substr(out, 1, 120))
})

## ---- R-e: plot(type = "cumul") under serial = "ar1" returns a ggplot ---------------------------------------------------------------------
test_that("R-e: plot(type = 'cumul') on an ar1 fit: a ggplot from the output of effect(cumu = TRUE); the fit itself stops naming the argument", {
  f <- ar1_fect(P2, serial = "ar1", r = 1, se = TRUE, vartype = "parametric", nboots = 40, keep.sims = TRUE, seed = 7)
  expect_error(plot(f, type = "cumul"), 'serial = "ar1"', fixed = TRUE)
  ce <- effect(f, cumu = TRUE)
  p <- plot(ce, type = "cumul")
  expect_s3_class(p, "ggplot")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off()); expect_no_error(print(p))
})
