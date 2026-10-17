## Phases 2 and 3a: fect(serial = "ar1") -- GSC on an AR(1)-filtered panel with a joint (x, r).
##
## The estimator replaces each unit's outcome and covariates by V_it = Y_it - x Y_i,t-1 (previous
## period on the observed grid; a cell whose previous cell is missing is missing), fits gsynth to V,
## and turns the counterfactual back to the outcome scale starting at each treated unit's last
## observed pretreatment outcome. x = 0 is gsynth on the panel without its first period.
##
## References (fixtures/ar1-reference.rds from fixtures/ar1-reference-gen.R): the prototype's own
## filtered_att(), ssr_filtered() and rank_scores() (xuyiqing/gsynth-revisit wip/template-v2-2026-10-16,
## e1eebae, sim/ar1_options.R, sim/rank_rules.R) run unchanged on P1 =
## make_design(30, 5, 15, 5, 0.5, seed = 7); set.seed(1); ar1_paths(35, 20, 0.6). Pinned numbers below
## are from the same run (2026-10-16). Other references are analytic (noiseless oracle, invariances) or
## internal identities (x = 0 against gsynth on the panel without its first period, the recursion on the
## bootstrap draws). Test hooks: options(fect.ar1.x = ), options(fect.ar1.rounds = ).

skip_on_cran()

ref <- readRDS(test_path("fixtures", "ar1-reference.rds"))
P1 <- ar1_p1_long(ref)
ar1_fit <- function(d, formula = y ~ d, ...) ar1_fect(d, formula, serial = "ar1", ...)
fix_x <- function(x) withr::local_options(fect.ar1.x = x, .local_envir = parent.frame())
post_cells <- function(f) f$D.dat == 1 & f$I.dat == 1

## E1: x = 0 equals gsynth on the panel without its first period -----------------------------------
test_that("E1: x = 0 equals levels GSC on the panel without its first period", {
  cases <- list(P2 = list(d = ar1_P2(), fm = y ~ d, tol = 1e-10),
                P3 = list(d = ar1_P3(), fm = y ~ d, tol = 1e-10),
                P5 = list(d = ar1_P5(), fm = y ~ d + x, tol = 1e-8))
  for (nm in names(cases)) for (r in 0:2) {
    cs <- cases[[nm]]; fix_x(0)
    f <- ar1_fit(cs$d, cs$fm, r = r); g <- ar1_fit_levels(cs$d, r, cs$fm)
    lab <- sprintf("%s r=%d", nm, r)
    expect_identical(f$ar1$x, 0); expect_identical(as.integer(f$ar1$first.period), 1L)
    expect_close(f$att.avg, g$att.avg, cs$tol, paste(lab, "att.avg"))
    expect_close(f$att.avg.unit, g$att.avg.unit, cs$tol, paste(lab, "att.avg.unit"))
    expect_close(f$eff, g$eff, cs$tol, paste(lab, "eff"))
    expect_close(f$att, g$att, cs$tol, paste(lab, "att"))
    expect_identical(as.numeric(f$count), as.numeric(g$count))
  }
})

## P4 has cells whose previous cell is missing. Section 1 of the spec defines V as missing there,
## so x = 0 is gsynth on the hand-filtered panel (not on the panel merely without its first period).
test_that("E1 (P4): x = 0 equals gsynth on the hand-filtered panel", {
  P4 <- ar1_P4()
  for (r in 0:2) {
    fix_x(0)
    f <- ar1_fit(P4, r = r); g <- ar1_fect(ar1_hand_filter(P4, 0), CV = FALSE, r = r, min.T0 = 4)
    expect_close(f$att.avg, g$att.avg, 1e-8, sprintf("P4 hand-filter r=%d att.avg", r))
    expect_close(f$eff, g$eff, 1e-8, sprintf("P4 hand-filter r=%d eff", r))
    expect_close(f$att, g$att, 1e-8, sprintf("P4 hand-filter r=%d att", r))
  }
})

## Spec E1 for P4 as written compares with the panel without its first period only. The two panels
## differ on P4 (cells with a missing previous cell are kept there and dropped by the filter), so the
## literal reference cannot match the filter's definition in section 1; recorded for the planner.
test_that("E1 [spec literal, P4]: x = 0 equals the panel without its first period (1e-8)", {
  P4 <- ar1_P4()
  for (r in 0:2) {
    fix_x(0)
    f <- ar1_fit(P4, r = r); g <- ar1_fit_levels(P4, r)
    expect_close(f$att.avg, g$att.avg, 1e-8, sprintf("P4 literal r=%d att.avg", r))
  }
})

## E2: fixed x against the prototype's filtered_att on P1 ------------------------------------------
test_that("E2: fixed x equals the prototype's filtered_att on P1", {
  xs <- c(0, 0.3, 0.6, 1)
  pinned <- rbind(  # prototype filtered_att(Y, co, tr, 15, r, x), r = 0..3 by x = 0, 0.3, 0.6, 1 (12 decimals)
    c(0.115295052654, 0.005621577564, -0.204948065207, -0.867960311177),
    c(0.458182771194, 0.411075101080, 0.281077781387, 0.015016207046),
    c(-0.050194455724, -0.146103382435, -0.242886460870, -0.353363836938),
    c(0.249111928246, 0.066922076736, -0.198143227104, -0.362029704756))
  expect_equal(unname(ref$P1$filtered), pinned, tolerance = 1e-10)           # the stored reference is the pinned table
  expect_equal(unname(ref$P1$levels_att), c(0.264004371078, 0.598640713046, -0.018260323129, -0.149502370658),
               tolerance = 1e-10)
  got <- matrix(NA_real_, 4, 4)
  for (i in 1:4) for (j in 1:4) {
    fix_x(xs[j]); got[i, j] <- ar1_fit(P1, r = i - 1)$att.avg
  }
  expect_close(got, pinned, 1e-8, "P1 fixed-x att.avg (r = 0..3 by x)")
  ## levels GSC at r = 0..3 (the reference for the unfiltered fits)
  lev <- vapply(0:3, function(r) ar1_fect(P1, CV = FALSE, r = r)$att.avg, 0)
  expect_close(lev, ref$P1$levels_att, 1e-8, "P1 levels att.avg")
})

## E3: joint x ------------------------------------------------------------------------------------------
test_that("E3: the joint x equals the prototype's on P1 (fixed r), att.avg at that x, profile", {
  pinned_x <- c(0.5107110861, 0.5419791563, 0.4259964905, 0.3806170575)   # r = 0..3
  expect_close(ref$P1$x_joint, pinned_x, 1e-9, "stored x_joint vs pinned")
  for (r in 0:3) {
    f <- ar1_fit(P1, r = r)
    expect_close(f$ar1$x, ref$P1$x_joint[r + 1], 1e-6, sprintf("P1 r=%d joint x", r))
    expect_close(f$ar1$x, pinned_x[r + 1], 1e-6, sprintf("P1 r=%d joint x (pinned)", r))
    expect_close(f$att.avg, ref$P1$att_at_joint[r + 1], 1e-6, sprintf("P1 r=%d att.avg at joint x", r))
    pr <- f$ar1$profile
    expect_true(all(round(seq(0, 1, by = 0.05), 10) %in% round(pr$x, 10)), info = "the 21 grid points are in the profile")
    expect_close(pr$x[which.min(pr$ssr)], f$ar1$x, 1e-5, sprintf("P1 r=%d profile minimum at x", r))
  }
})

test_that("E3 (boundaries): white-noise errors give x at or next to 0; random-walk errors give x near 1", {
  bx <- ref$boundary
  for (s in bx$wn$seed) {
    x <- ar1_fit(ar1_panel(30, 5, 15, 5, rho = 0, seed = 5000 + s), r = 2)$ar1$x
    expect_close(x, bx$wn$x[s], 1e-6, sprintf("white noise seed %d: x vs the prototype's", s))
    expect_lte(x, max(0.05, bx$wn$x[s] + 1e-6))     # at most 0.05 unless the prototype's own x is higher
  }
  for (s in bx$rw30$seed) {
    x <- ar1_fit(ar1_rw_panel(30, 5, 15, 5, 6000 + s), r = 2)$ar1$x
    expect_close(x, bx$rw30$x[s], 1e-6, sprintf("random walk (30 controls) seed %d", s))
    expect_gte(x, 0.8)
  }
  for (s in bx$rw100$seed) {
    x <- ar1_fit(ar1_rw_panel(100, 20, 40, 5, 7000 + s), r = 2)$ar1$x
    expect_close(x, bx$rw100$x[s], 1e-6, sprintf("random walk (100 controls) seed %d", s))
    expect_gte(x, 0.9)
  }
})

## Spec E3 as written: the joint x on white-noise panels is at most 0.05 (ten seeds). The prototype's
## own x on the seeds used here (5001..5010) is 0.0649 at seed 10, so the threshold is not stable.
test_that("E3 [spec threshold]: white-noise panels give x <= 0.05", {
  x <- vapply(1:10, function(s) ar1_fit(ar1_panel(30, 5, 15, 5, rho = 0, seed = 5000 + s), r = 2)$ar1$x, 0)
  expect_lte(max(x), 0.05)
})

## E4: alternation ----------------------------------------------------------------------------------------
test_that("E4: the alternation reproduces the reference loop on P1", {
  f <- ar1_fit(P1, r = c(0, 5))
  pth <- f$ar1$path; rp <- ref$P1$path
  expect_identical(nrow(pth), nrow(rp))
  expect_identical(as.integer(pth$r), as.integer(rp$r))
  expect_identical(as.integer(pth$round), as.integer(rp$round))
  expect_true(is.na(pth$x[1]))
  expect_close(pth$x[-1], rp$x[-1], 1e-6, "path x")
  ## pinned: round 0 r = 3; round 1 x = 0.3806170575, r = 2; round 2 x = 0.4259964905, r = 2 (stop)
  expect_identical(as.integer(pth$r), c(3L, 2L, 2L))
  expect_close(pth$x[-1], c(0.3806170575, 0.4259964905), 1e-6, "path x (pinned)")
  expect_close(f$att.avg, -0.189588027968, 1e-6, "final att.avg (pinned)")
  expect_close(f$att.avg, rp$att[3], 1e-6, "final att.avg (reference loop)")
  expect_true(isTRUE(f$ar1$converged)); expect_identical(as.integer(f$r.cv), 2L)
  expect_identical(f$serial, "ar1")
  ## one round
  withr::local_options(fect.ar1.rounds = 1)
  f1 <- ar1_fit(P1, r = c(0, 5))
  expect_close(f1$ar1$x, 0.3806170575, 1e-6, "one-round x")
  expect_identical(as.integer(f1$ar1$r), 2L)
  expect_close(f1$att.avg, ref$P1$one_round$att, 1e-6, "one-round att.avg (filtered_att at x_1, r_1)")
  ## a single fixed r: no alternation
  withr::local_options(fect.ar1.rounds = NULL)
  f2 <- ar1_fit(P1, r = 2)
  expect_true(is.na(f2$ar1$converged))
})

## E5: noiseless oracle ------------------------------------------------------------------------------------
noiseless_variants <- function(tau) list(
  balanced  = list(d = ar1_panel(30, 5, 15, 5, rho = 0.6, seed = 201, noiseless = TRUE, tau = tau), fm = y ~ d),
  staggered = list(d = ar1_panel(30, 6, NA, 20, rho = 0.6, seed = 202, noiseless = TRUE, tau = tau,
                                 T0.unit = c(11, 11, 13, 13, 15, 15)), fm = y ~ d),
  missing   = list(d = ar1_drop_cells(ar1_panel(30, 5, 15, 5, rho = 0.6, seed = 203, noiseless = TRUE, tau = tau),
                                      seed = 204), fm = y ~ d),
  covariate = list(d = ar1_panel(30, 5, 15, 5, rho = 0.6, seed = 205, noiseless = TRUE, tau = tau, cov = TRUE),
                   fm = y ~ d + x))

test_that("E5: noiseless panels recover the effect exactly at x = rho", {
  fix_x(0.6)
  for (tau in c(0, 2)) {
    v <- noiseless_variants(tau)
    f <- ar1_fit(v$balanced$d, r = 2)
    cells <- post_cells(f)
    expect_close(f$att.avg, tau, 1e-8, sprintf("balanced tau=%g att.avg", tau))
    expect_close(f$eff[cells], rep(tau, sum(cells)), 1e-8, sprintf("balanced tau=%g every post eff", tau))
    for (nm in c("staggered", "missing", "covariate")) {
      f <- ar1_fit(v[[nm]]$d, v[[nm]]$fm, r = 2)
      cells <- post_cells(f) & is.finite(f$eff)
      expect_gt(sum(cells), 15)
      expect_close(f$eff[cells], rep(tau, sum(cells)), 1e-6, sprintf("%s tau=%g every rebuilt post eff", nm, tau))
      ## every observed treated post cell is rebuilt
      expect_identical(sum(cells), sum(post_cells(f)), info = sprintf("%s tau=%g: all observed post cells rebuilt", nm, tau))
    }
  }
})

test_that("E5: the oracle can fail (x = 0 and x = 0.3 are not exact)", {
  d <- ar1_panel(30, 5, 15, 5, rho = 0.6, seed = 201, noiseless = TRUE, tau = 0)
  for (x in c(0, 0.3)) { fix_x(x); expect_gt(abs(ar1_fit(d, r = 2)$att.avg), 1e-3) }
})

## E6: invariances ------------------------------------------------------------------------------------------
test_that("E6: unit shifts, unit drift at x = 1, and scale", {
  P2 <- ar1_P2(); fix_x(0.5)
  base <- ar1_fit(P2, r = 2)
  for (u in c(32, 7)) {                                   # one treated unit, one control
    d <- P2; d$y[d$id == u] <- d$y[d$id == u] + 5
    expect_close(ar1_fit(d, r = 2)$eff, base$eff, 1e-8, sprintf("unit %d shifted by 5", u))
  }
  fix_x(1)
  b1 <- ar1_fit(P2, r = 2)
  d <- P2; d$y[d$id == 32] <- d$y[d$id == 32] + 0.7 * d$time[d$id == 32]
  expect_close(ar1_fit(d, r = 2)$eff, b1$eff, 1e-8, "treated unit with drift 0.7 t at x = 1")
  ## scale
  fix_x(0.5)
  P5 <- ar1_P5(); s <- 1e8
  for (nm in c("P2", "P5")) {
    dd <- if (nm == "P2") P2 else P5; fm <- if (nm == "P2") y ~ d else y ~ d + x
    b <- ar1_fit(dd, fm, r = 2); d2 <- dd; d2$y <- s * d2$y
    g <- ar1_fit(d2, fm, r = 2)
    expect_close(g$eff, s * b$eff, 1e-6, paste(nm, "eff scales with y"), rel = TRUE)
    expect_close(g$att.avg, s * b$att.avg, 1e-6, paste(nm, "att.avg scales with y"), rel = TRUE)
    ## normalize
    nz <- ar1_fit(dd, fm, r = 2, normalize = TRUE)
    expect_close(nz$att.avg, b$att.avg, 1e-6, paste(nm, "normalize = TRUE equals FALSE"))
  }
})

## E7: staggered, unbalanced, covariate runs -----------------------------------------------------------------
test_that("E7: staggered, unbalanced and covariate panels run and rebuild", {
  for (nm in c("P3", "P4", "P5")) {
    d <- switch(nm, P3 = ar1_P3(), P4 = ar1_P4(), P5 = ar1_P5()); fm <- if (nm == "P5") y ~ d + x else y ~ d
    expect_no_warning(f <- ar1_fect(d, fm, serial = "ar1", r = c(0, 5)))
    rb <- f$ar1$rebuilt
    expect_identical(nrow(rb), length(f$tr))
    expect_named(rb, c("id", "anchor", "n.post", "ok"), ignore.order = TRUE)
    if (nm != "P4") expect_true(all(rb$ok), info = nm)
    expect_true(is.finite(f$att.avg))
    cells <- post_cells(f)
    expect_equal(sum(f$count[f$time > 0]), sum(cells & is.finite(f$eff)), info = nm)
  }
  ## the two treated units whose last pretreatment cell is missing are rebuilt from an earlier anchor
  f <- ar1_fit(ar1_P4(), r = c(0, 5))
  rb <- f$ar1$rebuilt; rb <- rb[rb$id %in% c("31", "32"), ]
  expect_true(all(rb$ok))
  expect_true(all(rb$anchor < 15), info = "anchor earlier than T0 = 15 for units 31, 32")
})

## E8: units that cannot be rebuilt -------------------------------------------------------------------------------
test_that("E8: a treated unit with a covariate gap is reported, warned about, and does not affect the others", {
  P5 <- ar1_P5(); fix_x(0.5)
  d <- P5; bad <- 32
  d$x[d$id == bad & d$time %in% c(15, 16)] <- NA
  expect_warning(f <- ar1_fit(d, y ~ d + x, r = 2))
  expect_true(as.character(bad) %in% f$ar1$not.rebuilt)
  g <- ar1_fit(P5[P5$id != bad, ], y ~ d + x, r = 2)
  expect_close(f$att.avg, g$att.avg, 1e-8, "att.avg without the unit that cannot be rebuilt")
})

## E9: dropped first period and units ------------------------------------------------------------------------------
test_that("E9: the grid omits the first period; units with too few filtered pretreatment cells are dropped", {
  P2 <- ar1_P2()
  f <- ar1_fit(P2, r = 1)
  expect_identical(as.integer(f$rawtime), 2:20)
  expect_identical(nrow(f$Y.dat), 19L)
  Ymat <- matrix(P2$y, 20, 35)
  expect_close(f$Y.dat, Ymat[-1, ][, as.integer(f$id)], 1e-12, "Y.dat is the original outcome on the grid")
  ## unit 31 observed in alternate periods before treatment: five original pretreatment cells, no two adjacent
  d <- P2[!(P2$id == 31 & P2$time <= 15 & !(P2$time %in% c(1, 3, 5, 7, 9))), ]
  f2 <- ar1_fit(d, r = 1)
  expect_true("31" %in% f2$ar1$dropped$id)
  expect_s3_class(f2$ar1$dropped, "data.frame"); expect_named(f2$ar1$dropped, c("id", "reason"))
  expect_true(is.finite(f2$att.avg))
})

## E10: pre-treatment rows ----------------------------------------------------------------------------------------------
test_that("E10: pre-treatment rows are one-step-ahead errors; event time 1 does not depend on the recursion", {
  P2 <- ar1_P2(); fix_x(0.5)
  f <- ar1_fit(P2, r = 2)
  expect_true(any(f$time <= 0)); expect_identical(f$ar1$pre.scale, "one-step-ahead")
  g <- ar1_fect(ar1_hand_filter(P2, 0.5), CV = FALSE, r = 2, min.T0 = 4)
  expect_close(f$att[f$time == 1], g$att[g$time == 1], 1e-8, "att at event time 1 equals gsynth on the filtered data")
  expect_close(f$att[f$time <= 0], g$att[g$time <= 0], 1e-8, "pre-treatment rows equal gsynth on the filtered data")
  ## noiseless, x = rho: the pre-treatment rows are 0
  fix_x(0.6)
  dn <- ar1_panel(30, 5, 15, 5, rho = 0.6, seed = 201, noiseless = TRUE)
  fn <- ar1_fit(dn, r = 2)
  expect_close(fn$att[fn$time <= 0], rep(0, sum(fn$time <= 0)), 1e-8, "noiseless pre-treatment rows at x = rho")
})

## E11: a gap in the time index -----------------------------------------------------------------------------------------------
test_that("E11: the lag crosses a gap in the time index", {
  d <- ar1_panel(30, 5, 6, 4, rho = 0.6, seed = 301)
  years <- c(2000:2002, 2004:2010)                         # no unit observed in 2003
  d2 <- d; d2$time <- years[d$time]
  fix_x(0.4)
  f <- ar1_fit(d2, r = 1); g <- ar1_fit(d, r = 1)
  expect_close(f$eff, g$eff, 1e-10, "eff, years 2000..2010 without 2003 against times 1..10")
  expect_close(f$att.avg, g$att.avg, 1e-10, "att.avg")
  expect_close(f$att, g$att, 1e-10, "att")
})

## E12/E13: determinism, slots, print ----------------------------------------------------------------------------------------
test_that("E12: determinism", {
  P2 <- ar1_P2()
  a <- ar1_fit(P2, r = c(0, 5), se = TRUE, vartype = "parametric", nboots = 20, seed = 11)
  b <- ar1_fit(P2, r = c(0, 5), se = TRUE, vartype = "parametric", nboots = 20, seed = 11)
  expect_identical(a$ar1, b$ar1); expect_identical(a$est.att, b$est.att); expect_identical(a$est.avg, b$est.avg)
})

test_that("E13: slots and print", {
  P2 <- ar1_P2()
  f <- ar1_fit(P2, r = c(0, 5))
  expect_identical(f$serial, "ar1")
  a <- f$ar1
  expect_true(is.list(a))
  expect_true(all(c("x", "x.fixed", "r", "path", "converged", "profile", "rebuilt", "not.rebuilt", "dropped",
                    "first.period", "min.T0", "pre.scale", "v.scale.slots") %in% names(a)))
  expect_true(is.numeric(a$x) && length(a$x) == 1); expect_true(is.logical(a$x.fixed))
  expect_true(is.data.frame(a$path)); expect_named(a$path, c("round", "x", "r", "att"), ignore.order = TRUE)
  expect_true(is.na(a$path$x[1]) && a$path$round[1] == 0)
  expect_true(is.data.frame(a$profile)); expect_named(a$profile, c("x", "ssr"), ignore.order = TRUE)
  expect_true(is.data.frame(a$rebuilt)); expect_true(is.character(a$not.rebuilt)); expect_true(is.character(a$v.scale.slots))
  expect_identical(a$r, f$r.cv)
  for (nm in c("att.avg", "att.avg.unit", "eff", "Y.dat", "D.dat", "I.dat", "time", "att", "count"))
    expect_false(is.null(f[[nm]]), info = nm)
  expect_equal(nrow(f$eff), 19)
  out <- capture.output(print(f)); expect_true(any(grepl("AR(1)", out, fixed = TRUE)))
  expect_true(any(grepl(sprintf("%.2f", f$ar1$x), out, fixed = TRUE)))
  fs <- ar1_fit(P2, r = c(0, 5), se = TRUE, vartype = "parametric", nboots = 20, seed = 1)
  out2 <- capture.output(print(fs)); expect_true(any(grepl("AR(1)", out2, fixed = TRUE)))
  expect_true(any(grepl(sprintf("%.2f", fs$ar1$x), out2, fixed = TRUE)))
  ## serial = "none" fits carry no ar1 state
  g <- ar1_fect(P2, CV = FALSE, r = 2)
  expect_true(is.null(g$ar1)); expect_false(identical(g$serial, "ar1"))
})

test_that("E14: performance", {
  skip_on_cran()
  d <- ar1_panel(100, 20, 30, 5, rho = 0.6, seed = 401)
  t1 <- system.time(ar1_fit(d, r = c(0, 5)))[["elapsed"]]
  t2 <- system.time(ar1_fit(ar1_P2(), r = c(0, 5)))[["elapsed"]]
  message(sprintf("E14 timings: 100 controls / 20 treated / T = 35: %.2f s; 30 / 5 / T = 20: %.2f s", t1, t2))
  expect_lt(t1, 10); expect_lt(t2, 3)
})

## Phase 3a: inference ---------------------------------------------------------------------------------------------------------
boot_args <- list(se = TRUE, vartype = "parametric", nboots = 40, keep.sims = TRUE, seed = 7)
ar1_boot <- function(d, formula = y ~ d, ..., serial = "ar1") do.call(ar1_fect, c(list(d, formula, serial = serial), boot_args, list(...)))

test_that("I1: draws run, shapes", {
  for (nm in c("P2", "P3", "P5")) {
    d <- switch(nm, P2 = ar1_P2(), P3 = ar1_P3(), P5 = ar1_P5()); fm <- if (nm == "P5") y ~ d + x else y ~ d
    f <- ar1_boot(d, fm, r = 1)
    Tn <- attr(d, "T") - 1; N <- length(unique(d$id))
    expect_equal(dim(f$eff.boot), c(Tn, N, 40))
    expect_length(f$att.avg.boot, 40); expect_equal(dim(f$att.boot)[2], 40)
    expect_true(all(c("S.E.", "CI.lower", "CI.upper", "p.value") %in% colnames(f$est.avg)))
    expect_true(all(c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value") %in% colnames(f$est.att)))
    expect_true(all(is.finite(f$est.avg)) && all(is.finite(f$est.att[, c("ATT", "S.E.")])))
    ## in each draw the treated units come first (columns 1..Ntr), in the order of f$tr
    ntr <- length(f$tr)
    expect_true(all(vapply(f$colnames.boot, function(cn) identical(as.character(cn[seq_len(ntr)]), as.character(f$id[f$tr])), TRUE)))
    cells <- post_cells(f)[, f$tr, drop = FALSE] & is.finite(f$eff[, f$tr, drop = FALSE])
    eb <- f$eff.boot[, seq_len(ntr), , drop = FALSE]
    ok <- vapply(seq_len(40), function(b) all(is.finite(eb[, , b][cells])), TRUE)
    expect_true(all(ok), info = paste(nm, ": treated post cells of eff.boot finite where eff is"))
  }
})

test_that("I2: x = 0 reproduces gsynth's parametric bootstrap on the panel without its first period", {
  P2 <- ar1_P2(); fix_x(0)
  f <- ar1_boot(P2, r = 1)
  d1 <- P2[P2$time > 1, ]
  g <- do.call(ar1_fect, c(list(d1, y ~ d, CV = FALSE, r = 1, min.T0 = 4), boot_args))
  expect_close(f$est.avg, g$est.avg, 1e-12, "est.avg")
  expect_close(f$est.att, g$est.att, 1e-12, "est.att")
  expect_close(f$att.boot, g$att.boot, 1e-12, "att.boot")
  expect_close(f$att.avg.boot, g$att.avg.boot, 1e-12, "att.avg.boot")
  expect_close(f$eff.boot, g$eff.boot, 1e-12, "eff.boot")
})

test_that("I3: the recursion on the draws", {
  P2 <- ar1_P2(); fix_x(0.5); x <- 0.5
  f <- ar1_boot(P2, r = 2)
  g <- do.call(ar1_fect, c(list(ar1_hand_filter(P2, x), CV = FALSE, r = 2, min.T0 = 4), boot_args))
  ntr <- length(f$tr); expect_identical(f$tr, g$tr)
  expect_identical(f$colnames.boot, g$colnames.boot)                       # same draws of units
  eo <- f$eff.boot[, seq_len(ntr), , drop = FALSE]; ev <- g$eff.boot[, seq_len(ntr), , drop = FALSE]
  tg <- f$rawtime; anchor <- 15; rows <- which(tg > anchor); pre <- which(tg <= anchor)
  worst <- 0; worst_pre <- 0
  for (b in 1:40) for (j in seq_len(ntr)) {
    e <- numeric(length(tg)); e[pre] <- ev[pre, j, b]                      # cells at or before the anchor: one-step-ahead
    prev <- 0
    for (k in rows) { e[k] <- x * prev + ev[k, j, b]; prev <- e[k] }
    worst <- max(worst, abs(eo[, j, b] - e), na.rm = TRUE)
    worst_pre <- max(worst_pre, abs(eo[pre, j, b] - ev[pre, j, b]), na.rm = TRUE)
  }
  expect_true(worst <= 1e-10, info = sprintf("max |outer draw - recursion| = %.3e", worst))
  expect_true(worst_pre <= 1e-10, info = sprintf("cells at or before the anchor differ from eV by %.3e", worst_pre))
  ## the draws of att.avg are the means of the mapped treated post cells
  m <- vapply(1:40, function(b) mean(eo[rows, , b]), 0)
  expect_close(f$att.avg.boot, m, 1e-10, "att.avg.boot is the mean of the mapped treated post cells")
})

test_that("I4: intervals", {
  P2 <- ar1_P2(); fix_x(0.5)
  f <- ar1_boot(P2, r = 2)
  expect_close(f$est.avg[, "S.E."], sd(f$att.avg.boot), 1e-12, "S.E. equals sd(att.avg.boot)")
  est <- f$est.avg[, 1]
  expect_close(f$est.avg[, "CI.lower"], est - qnorm(0.975) * f$est.avg[, "S.E."], 1e-10, "normal CI lower")
  expect_close(f$est.avg[, "CI.upper"], est + qnorm(0.975) * f$est.avg[, "S.E."], 1e-10, "normal CI upper")
  fb <- ar1_boot(P2, r = 2, ci.method = "basic")
  eb <- estimand(fb, "att", "overall", ci.method = "basic")
  expect_close(c(eb$ci.lo, eb$ci.hi), unname(fb$est.avg[, c("CI.lower", "CI.upper")]), 1e-10, "basic interval equals estimand(..., ci.method = 'basic')")
})

test_that("I5: post-hoc estimands agree with the fit", {
  P2 <- ar1_P2(); fix_x(0.5)
  f <- ar1_boot(P2, r = 2)
  ov <- estimand(f, "att", "overall")
  expect_close(c(ov$ci.lo, ov$ci.hi), unname(f$est.avg[, c("CI.lower", "CI.upper")]), 1e-8, "estimand overall CI")
  expect_close(ov$estimate, unname(f$est.avg[, 1]), 1e-8, "estimand overall estimate")
  ev <- estimand(f, "att", "event.time")
  post <- f$est.att[as.numeric(rownames(f$est.att)) > 0, , drop = FALSE]
  expect_close(ev$estimate[ev$event.time > 0], unname(post[, "ATT"]), 1e-8, "estimand event.time estimate")
  cu <- att.cumu(f); ef <- effect(f)
  expect_false(is.null(cu)); expect_false(is.null(ef))
  expect_close(cu[1, "catt"], unname(f$att[f$time == 1]), 1e-8, "cumulative effect in the first post period equals att there")
  io <- imputed_outcomes(f)
  expect_true(is.data.frame(io)); expect_true(all(is.finite(io$Y0_hat)))
  expect_close(io$Y_obs - io$eff, io$Y0_hat, 1e-10, "Y_obs - eff equals the imputed counterfactual")
  ij <- cbind(match(io$time, f$rawtime), match(as.character(io$id), as.character(f$id)))
  expect_close(io$Y0_hat, f$Y.ct[ij], 1e-10, "imputed counterfactual equals Y.ct at the treated post cells")
  expect_equal(nrow(io), sum(post_cells(f)))
  ir <- imputed_outcomes(f, replicates = TRUE)
  expect_equal(nrow(ir), 40 * nrow(io))
  expect_close(ir$Y_obs - ir$eff, ir$Y0_hat, 1e-10, "replicates: Y_obs - eff equals the imputed counterfactual")
})

test_that("I6: vartype handling", {
  P2 <- ar1_P2()
  msgs <- testthat::capture_messages(f <- do.call(ar1_fect, c(list(P2, y ~ d, serial = "ar1", r = 1, quiet = FALSE),
                                                               modifyList(boot_args, list(vartype = "bootstrap")))))
  expect_identical(f$vartype, "parametric")
  expect_equal(sum(grepl('vartype', msgs, fixed = TRUE)), 1)          # one message about the parametric bootstrap
  expect_error(do.call(ar1_fect, c(list(P2, y ~ d, serial = "ar1", r = 1), modifyList(boot_args, list(vartype = "jackknife")))),
               'serial = "ar1"', fixed = TRUE)
})

test_that("I7: missing data and staggered adoption", {
  f <- ar1_boot(ar1_P4(), r = 1, para.error = "auto")
  expect_true(all(is.finite(f$est.avg)))
  f3 <- ar1_boot(ar1_P3(), r = 1)
  expect_true(all(is.finite(f3$est.avg)))
})

test_that("I8: requirements", {
  d <- ar1_panel(1, 5, 15, 5, rho = 0.6, seed = 801)
  expect_error(ar1_boot(d, r = 0), "needs at least two never-treated units")      # the parametric bootstrap's own error
  P2 <- ar1_P2()
  f <- ar1_boot(P2, r = 1, keep.sims = FALSE)
  expect_true(is.null(f$eff.boot)); expect_false(is.null(f$att.boot))
})
