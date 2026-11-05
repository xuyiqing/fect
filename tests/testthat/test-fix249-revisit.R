## fect 2.4.8 fixes from the gsynth-revisit reanalysis (2026-10-27).
##
## B1. Exploding parametric-bootstrap draws. In 02 Gilens "tort" (r = 2, 26
##     controls, 18 of them with a constant outcome) 5 of 1,000 simulated-error
##     draws (Step 1) came from refits whose pre-treatment factors were nearly
##     collinear with the unit intercept (condition number 1e11-1e13). Their
##     prediction errors were 5-19 million, and the 94 of 1,000 replicates that
##     reused one of them gave an SE of 158,668 for an ATT of -3.24. fect now
##     reports a robust spread and the count of extreme draws in boot.diag and
##     warns when SD > 5 * IQR/1.349. The reported SE does not change.
## B3. IFE-EM: print() and summary() report the rounds and the final relative
##     change (fit$convergence); the C++ EM loops return the final change.

skip_on_cran()

test_that(".fect_boot_diag flags a heavy-tailed set of draws and leaves a normal one alone", {
  set.seed(1)
  ok <- rnorm(1000)
  bd <- fect:::.fect_boot_diag(ok, nboots = 1000)
  expect_false(bd$breakdown)
  expect_equal(bd$sd, sd(ok))
  expect_equal(bd$iqr.sd, IQR(ok) / 1.349)
  expect_lt(bd$sd.ratio, 1.2)
  expect_identical(bd$n.extreme, 0L)
  expect_identical(bd$n.used, 1000L)

  bad <- ok
  bad[1:30] <- 2e5   # 3% of the draws explode, as in Gilens tort
  bd <- fect:::.fect_boot_diag(bad, nboots = 1000, n.failed = 2L)
  expect_true(bd$breakdown)
  expect_identical(bd$n.extreme, 30L)
  expect_identical(bd$n.failed, 2L)
  expect_gt(bd$sd.ratio, 1000)
  expect_lt(abs(bd$iqr.sd - 1), 0.15)   # the robust spread ignores the 30 draws
  expect_match(fect:::.fect_boot_diag_message(bd), "break down")

  ## non-finite draws are counted and left out of the spread
  bd <- fect:::.fect_boot_diag(c(ok, Inf, -Inf), nboots = 1002)
  expect_identical(bd$n.nonfinite, 2L)
  expect_identical(bd$n.used, 1000L)
  expect_false(bd$breakdown)
})

test_that(".fect_para_sim_check counts failed and extreme simulated-error draws", {
  set.seed(2)
  e <- array(rnorm(10 * 3 * 200), dim = c(10, 3, 200))
  e[, , 5] <- NA                 # a failed draw
  e[, , 7] <- e[, , 7] * 1e6     # an exploding draw
  sc <- fect:::.fect_para_sim_check(e)
  expect_identical(sc$n.sim, 200L)
  expect_identical(sc$n.sim.failed, 1L)
  expect_identical(sc$n.sim.extreme, 1L)
})

test_that("parametric bootstrap: boot.diag sits beside an unchanged S.E.", {
  data(simgsynth, package = "fect")
  out <- NULL
  expect_no_warning(out <- suppressMessages(fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    method = "gsynth", force = "two-way", r = 2, CV = FALSE, se = TRUE,
    vartype = "parametric", nboots = 100, parallel = FALSE, seed = 11)))
  bd <- out$boot.diag
  expect_type(bd, "list")
  b <- as.numeric(out$att.avg.boot)
  expect_equal(bd$sd, unname(out$est.avg[1, "S.E."]), tolerance = 1e-12)   # the S.E. is the SD of the draws
  expect_equal(bd$iqr.sd, IQR(b) / 1.349)
  expect_false(bd$breakdown)
  expect_identical(bd$nboots, 100L)
  expect_identical(bd$n.sim, 100L)
  expect_identical(bd$n.sim.extreme, 0L)
  ## working fields of fect_boot do not leak onto the fit
  expect_null(out$n.boot.fail)
  expect_null(out$para.sim.check)
  ## no breakdown, so print() does not show the check; summary() does
  expect_false(any(grepl("Bootstrap draws of the average effect", capture.output(print(out)))))
  expect_true(any(grepl("Bootstrap draws of the average effect", capture.output(print(summary(out))))))
})

test_that("case bootstrap gets boot.diag too; jackknife and se = FALSE do not", {
  data(simdata, package = "fect")
  d <- simdata[simdata$id %in% unique(simdata$id)[1:60], ]
  b <- suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"), method = "fe",
                             se = TRUE, nboots = 30, parallel = FALSE, seed = 3))
  expect_equal(b$boot.diag$sd, unname(b$est.avg[1, "S.E."]), tolerance = 1e-12)
  expect_null(b$boot.diag$n.sim)
  j <- suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"), method = "fe",
                             se = TRUE, vartype = "jackknife", parallel = FALSE))
  expect_null(j$boot.diag)
  n <- suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"), method = "fe", se = FALSE))
  expect_null(n$boot.diag)
})

## The real case. Its data are not shipped with fect; the test runs where the
## gsynth-revisit triage folder exists (FECT_TRIAGE_DIR, or Yiqing's Dropbox).
test_that("Gilens tort (r = 2): the breakdown warning fires and the S.E. stays 158,668", {
  dir <- Sys.getenv("FECT_TRIAGE_DIR",
    "~/Dropbox/Research_Hub/Packages/fect/triage/2026-10-14-gsynth-revisit-bugs")
  f <- file.path(path.expand(dir), "cases.rds")
  skip_if_not(file.exists(f), "gsynth-revisit triage cases not available")
  x <- readRDS(f)$explode_tort
  w <- character(0)
  out <- withCallingHandlers(
    suppressMessages(do.call(fect, c(list(data = x$data), x$args))),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("break down", w)))
  expect_equal(unname(out$est.avg[1, "S.E."]), 158667.8, tolerance = 1e-6)   # unchanged (bb435b9)
  expect_equal(out$att.avg, -3.244229, tolerance = 1e-6)
  bd <- out$boot.diag
  expect_true(bd$breakdown)
  expect_identical(bd$n.extreme, 94L)
  expect_gt(bd$n.sim.extreme, 0L)
  expect_lt(bd$iqr.sd, 2)
  expect_true(any(grepl("the bootstrap has broken down", capture.output(print(out)))))
})

test_that("IFE-EM fits report how they stopped (convergence, print, summary)", {
  data(simdata, package = "fect")
  out <- suppressMessages(fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                               method = "ife", r = 2, CV = FALSE, se = FALSE))
  cv <- out$convergence
  expect_type(cv, "list")
  expect_identical(cv$niter, as.integer(out$niter))
  expect_true(cv$converged)
  expect_lte(cv$rel.change, cv$tol)
  expect_gt(cv$rel.change, 0)
  expect_identical(cv$tol, 1e-5)
  expect_identical(cv$max.iteration, 5000L)
  expect_identical(out$rel.change, cv$rel.change)
  txt <- capture.output(print(out))
  expect_true(any(grepl(sprintf("EM convergence: stopped by tol after %d rounds", cv$niter), txt)))
  s <- summary(out)
  expect_s3_class(s, "summary.fect")
  expect_true(any(grepl("EM convergence: stopped by tol", capture.output(print(s)))))

  ## the final relative change is the one of the last round: one round short
  ## of convergence it is still above tol
  short <- NULL
  expect_warning(short <- suppressMessages(fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                                       method = "ife", r = 2, CV = FALSE, se = FALSE,
                                       max.iteration = cv$niter - 2)),
                 "did not converge")
  expect_false(short$convergence$converged)
  expect_gt(short$convergence$rel.change, cv$tol)
  expect_true(any(grepl("EM convergence: NOT converged", capture.output(print(short)))))
})

test_that("IFE-EM with cross-validation and r = 0 also report convergence; gsynth does not", {
  data(simdata, package = "fect")
  d <- simdata[simdata$id %in% unique(simdata$id)[1:80], ]
  cvfit <- suppressMessages(fect(Y ~ D + X1 + X2, data = d, index = c("id", "time"),
                                 method = "ife", r = c(0, 2), CV = TRUE, k = 3, se = FALSE,
                                 parallel = FALSE))
  expect_identical(cvfit$convergence$niter, as.integer(cvfit$niter))
  expect_true(is.finite(cvfit$convergence$rel.change))
  r0 <- suppressMessages(fect(Y ~ D + X1 + X2, data = d, index = c("id", "time"),
                              method = "ife", r = 0, CV = FALSE, se = FALSE))
  expect_true(is.finite(r0$convergence$rel.change))
  data(simgsynth, package = "fect")
  g <- suppressMessages(fect(Y ~ D, data = simgsynth, index = c("id", "time"),
                             method = "gsynth", r = 1, CV = FALSE, se = FALSE))
  expect_null(g$convergence)
  expect_false(any(grepl("EM convergence", capture.output(print(g)))))
})

## Review of PR #183 (2026-10-28): edge cases of the check and the convergence line.

test_that("the breakdown check needs 20 usable draws and a positive IQR", {
  set.seed(4)
  few <- c(rnorm(9), 1e4)                      # 10 draws, one wild
  bd <- fect:::.fect_boot_diag(few, nboots = 10)
  expect_gt(bd$sd.ratio, 5)
  expect_false(bd$breakdown)                   # too few draws to read
  expect_identical(bd$min.draws, 20L)
  many <- c(rnorm(29), 1e4)
  expect_true(fect:::.fect_boot_diag(many, nboots = 30)$breakdown)

  flat <- c(rep(1, 80), rnorm(20))             # IQR = 0
  bd <- fect:::.fect_boot_diag(flat, nboots = 100)
  expect_identical(bd$iqr.sd, 0)
  expect_true(is.na(bd$sd.ratio))
  expect_true(is.na(bd$n.extreme))             # not every draw off the median
  expect_false(bd$breakdown)
  txt <- capture.output(fect:::.print_boot_diag(bd))
  expect_true(any(grepl("not defined", txt)))
  expect_false(any(grepl("Inf", txt)))
})

test_that("the warning describes the first-step pool by error mode", {
  bd <- fect:::.fect_boot_diag(c(rnorm(99), 1e5), nboots = 100,
                               sim = list(n.sim = 100L, n.sim.failed = 0L, n.sim.extreme = 2L),
                               para.error = "empirical")
  m <- fect:::.fect_boot_diag_message(bd, r.cv = 2)
  expect_match(m, "reused by every replicate")
  bd$para.error <- "ar"
  m <- fect:::.fect_boot_diag_message(bd, r.cv = 2)
  expect_match(m, "error covariance")
  expect_false(grepl("reused by every replicate", m))
  expect_match(m, "stopped far from its solution")
})

test_that("summary() shows the first-step counts of a parametric bootstrap", {
  data(simgsynth, package = "fect")
  out <- suppressMessages(fect(Y ~ D, data = simgsynth, index = c("id", "time"),
                               method = "gsynth", force = "two-way", r = 1, CV = FALSE, se = TRUE,
                               vartype = "parametric", nboots = 50, parallel = FALSE, seed = 2))
  expect_identical(out$boot.diag$para.error, out$para.error)
  txt <- capture.output(print(summary(out)))
  expect_true(any(grepl("simulated-error draws \\(first step\\): 50;  failed: 0", txt)))
})

test_that("the convergence line reports the cap the EM loop applied", {
  data(simdata, package = "fect")
  fe <- suppressMessages(fect(Y ~ D, data = simdata, index = c("id", "time"),
                              method = "fe", se = FALSE))
  ## since 2.4.8 fe_ad_iter() (no factors, no covariates) obeys max.iteration
  ## too; before, it stopped at 500 rounds
  expect_identical(fe$convergence$iter.cap, 5000L)
  expect_null(fe$iter.cap)
  ife <- suppressMessages(fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                               method = "ife", r = 1, CV = FALSE, se = FALSE))
  expect_identical(ife$convergence$iter.cap, 5000L)
  capped <- list(niter = 5001L, rel.change = 2e-5, tol = 1e-5, max.iteration = 5000L,
                 iter.cap = 5000L, converged = FALSE)
  txt <- capture.output(fect:::.print_convergence(capped))
  expect_match(txt, "NOT converged; stopped at max.iteration = 5000")
  expect_false(any(grepl("routine's own cap", txt)))
})
