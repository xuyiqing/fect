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
## bootstrap draws). Test hooks: options(fect.ar1.x = ), options(fect.ar1.rounds = ),
## options(fect.ar1.exclude.nonconverged = ), options(fect.ar1.mask.zero = ).
##
## Round 3 (test-spec "Amendments 2"): at x exactly 0 the outcome is used unchanged (V = Y) on the panel
## without its first period, with no lag mask, so x = 0 equals gsynth on that panel on every panel (M1);
## the hand-filtered (masked) definition applies at x = 1e-6 (M2). The alternation stops at the first repeated
## rank (A1-A10, X1), and an inner EM that reaches max.iteration is reported and ineligible (X2, N1-N4).

## A: the alternation logic with a stub (no data, no fits) -------------------------------------------------
## step(r) maps the rank fed to the x-step to list(x, r); every call is recorded.
alt_stub <- function(tab, r0, K = 5L) {
  calls <- integer(0)
  step <- function(r) {
    calls <<- c(calls, as.integer(r))
    e <- tab[[as.character(r)]]
    if (is.null(e)) stop("step called with a rank outside the table: ", r)
    list(x = e[[1]], r = as.integer(e[[2]]))
  }
  res <- fect:::.fect_ar1_alternate(step, r0, K = K)
  res$calls <- calls
  res
}
alt_pair <- function(res, round) { row <- res$path[res$path$round == round, ]; c(x = row$x, r = row$r) }
tab_A3 <- list("2" = c(0.1, 3), "3" = c(0.2, 2))

test_that("A1: a repeat of the previous rank is convergence", {
  res <- alt_stub(list("2" = c(0.3, 2)), 2)
  expect_identical(as.integer(res$path$round), 0:1)
  expect_identical(res$final, 1L); expect_true(res$converged); expect_false(res$cycled)
  expect_identical(res$cycle, integer(0)); expect_false(res$capped); expect_identical(res$calls, 2L)
  expect_true(is.na(res$path$x[1])); expect_identical(as.integer(res$path$r[1]), 2L)
})

test_that("A2: convergence after a change of rank", {
  res <- alt_stub(list("3" = c(0.4, 2), "2" = c(0.5, 2)), 3)
  expect_identical(res$final, 2L); expect_equal(alt_pair(res, res$final), c(x = 0.5, r = 2))
  expect_true(res$converged); expect_false(res$cycled); expect_identical(res$calls, c(3L, 2L))
})

test_that("A3/A4: a repeat of an earlier rank is a cycle; the pair with the lowest rank is reported", {
  a3 <- alt_stub(tab_A3, 2)
  expect_equal(unname(a3$path$x[-1]), c(0.1, 0.2)); expect_identical(as.integer(a3$path$r[-1]), c(3L, 2L))
  expect_true(a3$cycled); expect_false(a3$converged); expect_identical(as.integer(a3$cycle), 1:2)
  expect_identical(a3$final, 2L); expect_equal(alt_pair(a3, a3$final), c(x = 0.2, r = 2))
  expect_identical(a3$calls, c(2L, 3L))
  a4 <- alt_stub(tab_A3, 3)
  expect_equal(unname(a4$path$x[-1]), c(0.2, 0.1)); expect_identical(as.integer(a4$path$r[-1]), c(2L, 3L))
  expect_true(a4$cycled); expect_false(a4$converged); expect_identical(as.integer(a4$cycle), 1:2)
  expect_identical(a4$final, 1L); expect_equal(alt_pair(a4, a4$final), c(x = 0.2, r = 2))     # the same pair as A3
})

test_that("A5/A6: the cycle starts where the repeated rank first appeared", {
  a5 <- alt_stub(list("0" = c(0.05, 2), "2" = c(0.1, 3), "3" = c(0.2, 2)), 0)
  expect_identical(as.integer(a5$path$r[-1]), c(2L, 3L, 2L)); expect_identical(nrow(a5$path), 4L)
  expect_identical(as.integer(a5$cycle), 2:3)            # the repeat is of round 1, so the cycle is rounds 2 and 3
  expect_identical(a5$final, 3L); expect_equal(alt_pair(a5, 3L), c(x = 0.2, r = 2))
  expect_true(a5$cycled); expect_false(a5$converged); expect_identical(a5$calls, c(0L, 2L, 3L))
  a6 <- alt_stub(list("1" = c(0.1, 4), "4" = c(0.2, 2), "2" = c(0.3, 1)), 1)
  expect_identical(as.integer(a6$cycle), 1:3); expect_identical(as.integer(a6$path$r[-1]), c(4L, 2L, 1L))
  expect_identical(a6$final, 3L); expect_equal(alt_pair(a6, 3L), c(x = 0.3, r = 1)); expect_true(a6$cycled)
})

test_that("A7/A8: the cap on rounds", {
  tab <- list("0" = c(0.1, 1), "1" = c(0.1, 2), "2" = c(0.1, 3), "3" = c(0.1, 4), "4" = c(0.1, 5), "5" = c(0.1, 5))
  a7 <- alt_stub(tab, 0, K = 5L)
  expect_true(a7$capped); expect_false(a7$converged); expect_false(a7$cycled)
  expect_identical(a7$final, 5L); expect_length(a7$calls, 5L)
  a7b <- alt_stub(tab, 0, K = 6L)
  expect_identical(a7b$final, 6L); expect_true(a7b$converged); expect_false(a7b$capped); expect_length(a7b$calls, 6L)
  a8 <- alt_stub(list("2" = c(0.1, 3)), 2, K = 1L)
  expect_true(a8$capped); expect_false(a8$converged); expect_identical(a8$final, 1L)
  a8b <- alt_stub(list("2" = c(0.1, 2)), 2, K = 1L)
  expect_true(a8b$converged); expect_false(a8b$capped)
})

test_that("A9: the reported pair does not depend on the cap (nor on its parity)", {
  for (K in 2:7) {
    res <- alt_stub(tab_A3, 2, K = K)
    expect_equal(alt_pair(res, res$final), c(x = 0.2, r = 2), info = paste("K =", K))
    expect_length(res$calls, 2L)
  }
})

test_that("A10: step is called with the previous round's rank, never after the decision, and traced once per call", {
  for (res in list(alt_stub(tab_A3, 2), alt_stub(list("0" = c(0.05, 2), "2" = c(0.1, 3), "3" = c(0.2, 2)), 0))) {
    expect_identical(res$calls, as.integer(res$path$r[seq_along(res$calls)]))     # c(r0, r_1, ...)
    expect_length(res$trace, length(res$calls))
    expect_identical(length(res$calls), nrow(res$path) - 1L)                      # no call after the decision
  }
})

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
    expect_false(f$ar1$lag.mask)                                           # M1: no lag mask at x = 0
    expect_close(f$att.avg, g$att.avg, cs$tol, paste(lab, "att.avg"))
    expect_close(f$att.avg.unit, g$att.avg.unit, cs$tol, paste(lab, "att.avg.unit"))
    expect_close(f$eff, g$eff, cs$tol, paste(lab, "eff"))
    expect_close(f$att, g$att, cs$tol, paste(lab, "att"))
    expect_identical(as.numeric(f$count), as.numeric(g$count))
  }
})

## M1 (amended E1 on P4, test-spec Amendments 2): at x exactly 0 the outcome is used unchanged (V = Y) on the panel
## without its first period and no cell is removed because its predecessor is missing, so x = 0 equals gsynth on the
## panel without its first period on every panel, P4 included (the literal reference of the original E1 is valid again).
test_that("M1: x = 0 equals gsynth on P4 without its first period (all cells, 1e-8), no lag mask", {
  P4 <- ar1_P4()
  for (r in 0:2) {
    fix_x(0)
    f <- ar1_fit(P4, r = r); g <- ar1_fit_levels(P4, r)
    lab <- sprintf("P4 r=%d", r)
    expect_false(f$ar1$lag.mask); expect_identical(f$ar1$x, 0)
    expect_identical(unname(is.finite(f$eff)), unname(is.finite(g$eff)))      # identical finite sets, incl. predecessor-missing post cells
    expect_close(f$att.avg, g$att.avg, 1e-8, paste(lab, "att.avg"))
    expect_close(f$att.avg.unit, g$att.avg.unit, 1e-8, paste(lab, "att.avg.unit"))
    expect_close(f$eff, g$eff, 1e-8, paste(lab, "eff (all cells)"))
    expect_close(f$att, g$att, 1e-8, paste(lab, "att"))
    expect_identical(as.numeric(f$count), as.numeric(g$count))
    ## every observed treated post cell is rebuilt, including those whose predecessor is missing (chain through gaps)
    expect_identical(sum(post_cells(f) & is.finite(f$eff)), sum(post_cells(f)))
    ar1_note(sprintf("M1 P4 r=%d: finite eff cells %d, finite treated post cells %d", r, sum(is.finite(f$eff)), sum(post_cells(f) & is.finite(f$eff))))
  }
})

## M1b: the contract says x = 0 equals gsynth on the panel without its first period on EVERY panel; here with a covariate
## and missing cells (controls, pretreatment cells, and the treated unit's own post cells), where V_x = x - 0 * lag(x) must
## not turn a missing predecessor into a missing covariate.
test_that("M1b: x = 0 equals gsynth on the panel without its first period with a covariate and missing cells", {
  P5m <- ar1_drop_cells(ar1_P5(), seed = 103, force.units = integer(0))
  P5m <- P5m[!(P5m$id == 32 & P5m$time == 16), ]                 # a treated post cell missing (outcome and covariate)
  for (r in 0:2) {
    fix_x(0)
    f <- ar1_fit(P5m, y ~ d + x, r = r); g <- ar1_fit_levels(P5m, r, y ~ d + x)
    expect_false(f$ar1$lag.mask)
    expect_identical(unname(is.finite(f$eff)), unname(is.finite(g$eff)), info = sprintf("r=%d finite cells", r))
    expect_close(f$att.avg, g$att.avg, 1e-8, sprintf("P5 with missing cells r=%d att.avg", r))
    expect_close(f$eff, g$eff, 1e-8, sprintf("P5 with missing cells r=%d eff", r))
  }
})

## M2 (the former E1 on P4): the masked definition applies at x = 1e-6. H is the hand-filtered panel (first period and cells
## with a missing predecessor removed). The chain passes through unobserved cells (E5, E7), so the build also rebuilds the
## observed treated post cells whose predecessor is missing; gsynth on H leaves them NA.
test_that("M2: x = 1e-6 keeps the lag mask; equals gsynth on the hand-filtered panel cell by cell (1e-4)", {
  P4 <- ar1_P4()
  obs <- with(P4, tapply(!is.na(y), list(time, id), any)); obs[is.na(obs)] <- FALSE
  for (r in 0:2) {
    fix_x(1e-6)
    f <- ar1_fit(P4, r = r); g <- ar1_fect(ar1_hand_filter(P4, 0), CV = FALSE, r = r, min.T0 = 4)
    expect_true(f$ar1$lag.mask)
    expect_identical(dim(f$eff), dim(g$eff)); expect_identical(as.character(f$id), as.character(g$id))
    fe <- is.finite(f$eff); ge <- is.finite(g$eff); both <- fe & ge
    expect_true(all(fe[ge]), info = sprintf("r=%d: the build keeps every cell of the reference", r))
    expect_close(f$eff[both], g$eff[both], 1e-4, sprintf("P4 hand-filter r=%d eff on shared cells", r))
    pred_miss <- matrix(FALSE, nrow(f$eff), ncol(f$eff))
    for (j in seq_len(ncol(f$eff))) for (i in seq_len(nrow(f$eff)))
      pred_miss[i, j] <- !obs[as.character(f$rawtime[i] - 1), as.character(f$id[j])]
    expected_extra <- post_cells(f) & pred_miss
    expect_identical(unname(fe & !ge), unname(expected_extra & fe), info = sprintf("r=%d: extra cells are predecessor-missing post cells", r))
    expect_gt(sum(fe & !ge), 0)
  }
})

ar1_P4c <- function() {
  ## controls and treated pretreatment cells missing at random; no treated post cell has a missing predecessor
  P4c <- ar1_drop_cells(ar1_P2(), seed = 103, force.units = integer(0))
  full <- ar1_P2()                                   # put back any dropped treated cell at t >= T0 (T0 = 15)
  add <- full[full$id > 30 & full$time >= 15 & !(paste(full$id, full$time) %in% paste(P4c$id, P4c$time)), ]
  P4c <- rbind(P4c, add); P4c[order(P4c$id, P4c$time), ]
}

test_that("M2 (P4c, no treated post gaps): x = 1e-6 equals gsynth on the hand-filtered panel (att.avg, eff, att, 1e-4)", {
  P4c <- ar1_P4c(); full <- ar1_P2()
  expect_lt(nrow(P4c), nrow(full))                   # controls and earlier cells are still missing
  for (r in 0:2) {
    fix_x(1e-6)
    f <- ar1_fit(P4c, r = r); g <- ar1_fect(ar1_hand_filter(P4c, 0), CV = FALSE, r = r, min.T0 = 4)
    expect_true(f$ar1$lag.mask)
    expect_identical(unname(is.finite(f$eff)), unname(is.finite(g$eff)))
    expect_close(f$att.avg, g$att.avg, 1e-4, sprintf("P4c hand-filter r=%d att.avg", r))
    expect_close(f$eff, g$eff, 1e-4, sprintf("P4c hand-filter r=%d eff", r))
    expect_close(f$att, g$att, 1e-4, sprintf("P4c hand-filter r=%d att", r))
  }
})

## M3 (informational): the two definitions differ, so M1 and M2 can fail (the tester measured 0.32 between the two references).
test_that("M3: on P4 at r = 1 the fits at x = 0 and x = 1e-6 differ by more than 1e-3", {
  P4 <- ar1_P4()
  fix_x(0); a0 <- ar1_fit(P4, r = 1)$att.avg
  fix_x(1e-6); a1 <- ar1_fit(P4, r = 1)$att.avg
  message(sprintf("M3: |att.avg(x = 0) - att.avg(x = 1e-6)| = %.4f on P4 at r = 1", abs(a0 - a1)))
  expect_gt(abs(a0 - a1), 1e-3)
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
  wn_x <- numeric(length(bx$wn$seed))
  for (s in bx$wn$seed) {
    x <- ar1_fit(ar1_panel(30, 5, 15, 5, rho = 0, seed = 5000 + s), r = 2)$ar1$x
    expect_close(x, bx$wn$x[s], 1e-6, sprintf("white noise seed %d: x vs the prototype's", s))
    expect_lte(x, 0.15)                              # amended bound (above the prototype's Monte Carlo maximum 0.111)
    wn_x[s] <- x
  }
  expect_gte(sum(wn_x == 0), 6)                      # amended: x = 0 for at least 6 of the 10 seeds
  rw_x <- numeric(length(bx$rw30$seed))
  for (s in bx$rw30$seed) {
    x <- ar1_fit(ar1_rw_panel(30, 5, 15, 5, 6000 + s), r = 2)$ar1$x
    expect_close(x, bx$rw30$x[s], 1e-6, sprintf("random walk (30 controls) seed %d", s))
    expect_gte(x, 0.7)                               # amended per-seed bound (Monte Carlo minimum 0.731)
    rw_x[s] <- x
  }
  expect_gte(median(rw_x), 0.8)                      # amended
  for (s in bx$rw100$seed) {
    x <- ar1_fit(ar1_rw_panel(100, 20, 40, 5, 7000 + s), r = 2)$ar1$x
    expect_close(x, bx$rw100$x[s], 1e-6, sprintf("random walk (100 controls) seed %d", s))
    expect_gte(x, 0.9)
  }
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

## Chain through gaps (round 2): the recursion starts at the last observed pretreatment outcome and passes
## through unobserved cells on the model's prediction. Noiseless P6 at x = rho: every rebuilt post cell is tau.
test_that("E5 (gaps): one missing cell does not break or restart the chain", {
  fix_x(0.6)
  d0 <- ar1_panel(30, 5, 15, 5, rho = 0.6, seed = 201, noiseless = TRUE, tau = 2)
  ## (i) Y(31, t = 15) missing: the last observed pretreatment outcome of unit 31 is at t = 14
  d <- d0[!(d0$id == 31 & d0$time == 15), ]
  f <- ar1_fit(d, r = 2)
  rb <- f$ar1$rebuilt; rb <- rb[rb$id == "31", ]
  expect_true(rb$ok); expect_lt(rb$anchor, 15)
  expect_identical(as.integer(rb$anchor), 14L)
  post <- f$eff[f$rawtime > 15, match("31", as.character(f$id))]
  expect_true(all(is.finite(post)), info = "unit 31: all five post cells rebuilt")
  expect_close(post, rep(2, length(post)), 1e-6, "unit 31 post effects equal tau")
  ## (ii) Y(32, t = 18) missing, a post-treatment cell: later cells are rebuilt through it
  d <- d0[!(d0$id == 32 & d0$time == 18), ]
  f <- ar1_fit(d, r = 2)
  rb <- f$ar1$rebuilt; rb <- rb[rb$id == "32", ]
  expect_true(rb$ok); expect_identical(as.integer(rb$anchor), 15L)
  obs <- f$rawtime > 15 & f$rawtime != 18
  post <- f$eff[obs, match("32", as.character(f$id))]
  expect_true(all(is.finite(post)), info = "unit 32: observed post cells rebuilt through the gap")
  expect_close(post, rep(2, length(post)), 1e-6, "unit 32 post effects equal tau")
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
                    "first.period", "min.T0", "pre.scale", "v.scale.slots",
                    "cycled", "cycle", "final.round", "em.converged", "niter", "excluded", "fallback.rounds",
                    "search.nonconverged", "lag.mask") %in% names(a)))
  expect_true(is.numeric(a$x) && length(a$x) == 1); expect_true(is.logical(a$x.fixed))
  expect_true(is.data.frame(a$path)); expect_named(a$path, c("round", "x", "r", "att", "em.converged"), ignore.order = TRUE)
  expect_true(is.na(a$path$x[1]) && a$path$round[1] == 0)
  expect_true(is.data.frame(a$profile)); expect_named(a$profile, c("x", "ssr", "em.converged"), ignore.order = TRUE)
  expect_true(is.logical(a$path$em.converged) && is.logical(a$profile$em.converged))
  expect_true(is.data.frame(a$rebuilt)); expect_true(is.character(a$not.rebuilt)); expect_true(is.character(a$v.scale.slots))
  expect_true(is.logical(a$cycled) && length(a$cycled) == 1 && !a$cycled)
  expect_true(is.data.frame(a$cycle)); expect_named(a$cycle, c("round", "x", "r", "chosen"), ignore.order = TRUE); expect_identical(nrow(a$cycle), 0L)
  expect_true(is.numeric(a$final.round) && length(a$final.round) == 1)
  expect_true(is.logical(a$em.converged) && length(a$em.converged) == 1); expect_true(is.numeric(a$niter))
  expect_true(is.data.frame(a$excluded)); expect_named(a$excluded, c("round", "x", "r"), ignore.order = TRUE)
  expect_true(is.numeric(a$fallback.rounds)); expect_true(is.numeric(a$search.nonconverged) && length(a$search.nonconverged) == 1)
  expect_true(is.logical(a$lag.mask) && length(a$lag.mask) == 1)
  expect_identical(a$lag.mask, a$x != 0)                                  # FALSE iff x = 0
  expect_true(is.logical(a$converged) && !a$cycled)
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

## M4: the x = 0 bootstrap identity on P4 and P4c ----------------------------------------------------------------------------
test_that("M4: x = 0 reproduces gsynth's parametric bootstrap on P4 and P4c without the first period (1e-12)", {
  fix_x(0)
  for (nm in c("P4", "P4c")) {
    d <- if (nm == "P4") ar1_P4() else ar1_P4c()
    f <- ar1_boot(d, r = 2)
    expect_false(f$ar1$lag.mask)
    g <- do.call(ar1_fect, c(list(d[d$time > 1, ], y ~ d, CV = FALSE, r = 2, min.T0 = 4), boot_args))
    expect_close(f$est.avg, g$est.avg, 1e-12, paste(nm, "est.avg"))
    expect_close(f$est.att, g$est.att, 1e-12, paste(nm, "est.att"))
    expect_close(f$att.boot, g$att.boot, 1e-12, paste(nm, "att.boot"))
    expect_close(f$att.avg.boot, g$att.avg.boot, 1e-12, paste(nm, "att.avg.boot"))
    expect_close(f$eff.boot, g$eff.boot, 1e-12, paste(nm, "eff.boot"))
  }
})

## M5: E8 at x = 0 -----------------------------------------------------------------------------------------------------------
test_that("M5: at x = 0 a covariate gap costs only its own cells; the unit is rebuilt", {
  P5 <- ar1_P5(); bad <- 32
  d <- P5; d$x[d$id == bad & d$time %in% c(15, 16)] <- NA
  fix_x(0)
  full <- ar1_fit(P5, y ~ d + x, r = 2)
  f <- ar1_collect(ar1_fit(d, y ~ d + x, r = 2))
  ar1_note("M5 warnings at x = 0: ", if (length(f$warnings)) paste(f$warnings, collapse = " | ") else "none")
  f <- f$value
  expect_false(as.character(bad) %in% f$ar1$not.rebuilt)
  rb <- f$ar1$rebuilt[f$ar1$rebuilt$id == as.character(bad), ]; rb0 <- full$ar1$rebuilt[full$ar1$rebuilt$id == as.character(bad), ]
  expect_true(rb$ok)
  expect_identical(as.integer(rb$n.post), as.integer(rb0$n.post) - 1L)       # the one treated post cell with a missing covariate (t = 16)
  g <- ar1_fit_levels(d, 2, y ~ d + x)
  expect_close(f$att.avg, g$att.avg, 1e-8, "att.avg equals gsynth on the data without the first period")
  ## x = 0.5: E8 as before
  fix_x(0.5)
  expect_warning(f5 <- ar1_fit(d, y ~ d + x, r = 2))
  expect_true(as.character(bad) %in% f5$ar1$not.rebuilt)
})

## M6/M7: the mask rule leaves the search and complete panels alone ------------------------------------------------------------
test_that("M6: the profile and the joint x do not depend on the mask rule (P4, fixed r)", {
  P4 <- ar1_P4()
  for (r in 1:2) {
    withr::local_options(fect.ar1.mask.zero = FALSE); a <- ar1_fit(P4, r = r)
    withr::local_options(fect.ar1.mask.zero = TRUE);  b <- ar1_fit(P4, r = r)
    expect_close(a$ar1$profile$x, b$ar1$profile$x, 1e-10, sprintf("r=%d profile x", r))
    expect_close(a$ar1$profile$ssr, b$ar1$profile$ssr, 1e-10, sprintf("r=%d profile ssr", r), rel = TRUE)
    expect_close(a$ar1$x, b$ar1$x, 1e-10, sprintf("r=%d x-hat", r))
    ar1_note(sprintf("M6 r=%d: x-hat %.6f, lag.mask %s (new) vs %s (earlier)", r, a$ar1$x, a$ar1$lag.mask, b$ar1$lag.mask))
  }
})

test_that("M7: complete panels give identical results with and without the mask-at-zero hook", {
  for (nm in c("P2", "P3", "P5")) {
    d <- switch(nm, P2 = ar1_P2(), P3 = ar1_P3(), P5 = ar1_P5()); fm <- if (nm == "P5") y ~ d + x else y ~ d
    withr::local_options(fect.ar1.mask.zero = FALSE); a <- ar1_fit(d, fm, r = c(0, 5))
    withr::local_options(fect.ar1.mask.zero = TRUE);  b <- ar1_fit(d, fm, r = c(0, 5))
    ar1_note(sprintf("M7 %s: x %.6f/%.6f r %d/%d cycled %s/%s excluded rows %d/%d", nm, a$ar1$x, b$ar1$x, a$r.cv, b$r.cv,
                     a$ar1$cycled, b$ar1$cycled, nrow(a$ar1$excluded), nrow(b$ar1$excluded)))
    if (!identical(a$ar1$cycled, b$ar1$cycled) || nrow(a$ar1$excluded) != nrow(b$ar1$excluded)) next    # reported, not compared
    expect_close(a$att.avg, b$att.avg, 1e-10, paste(nm, "att.avg")); expect_close(a$eff, b$eff, 1e-10, paste(nm, "eff"))
    expect_close(a$ar1$x, b$ar1$x, 1e-10, paste(nm, "x")); expect_identical(a$r.cv, b$r.cv)
  }
})

## D2 (Amendments 3): the "3c is not built" check is deleted. What stays is the part that is still true and
## independent of 3c: the point estimate (x, rank, ATT) does not depend on whether intervals are requested. The
## new bootstrap (T1-T11) is tested in test-ar1-boot-x.R.
test_that("D2: the point estimate and the chosen pair do not depend on se", {
  P2 <- ar1_P2()
  a <- ar1_fit(P2, r = c(0, 5)); b <- ar1_boot(P2, r = c(0, 5))
  expect_identical(a$r.cv, b$r.cv); expect_close(a$ar1$x, b$ar1$x, 1e-12, "x with and without se")
  expect_close(a$att.avg, b$att.avg, 1e-12, "att.avg with and without se")
})

## X: the real-data cycle and the runaway fit (H1 replication 14) ---------------------------------------------------------------
## Fixture fixtures/ar1-h1-rep14.rds: the long data frame (id, time, y, d) of cell 22 (H1) replication 14 of the simulation
## study, seed 20261017 + 100000 * 22 + 14 = 22461031, regenerated with the harness dgp.R (scratchpad ar1sim; provenance in
## the rds). 35 units, 20 periods, 633 of 700 cells, 5 treated, true effect 0. Pinned numbers are the diagnostic's
## (fect 2.4.8 build d4c3b7f, 2026-10-16).
h1 <- readRDS(test_path("fixtures", "ar1-h1-rep14.rds"))
H1 <- h1$data
h1_fit <- function(serial = "ar1", ...) {
  ar1_collect(ar1_fect(H1, serial = serial, r = c(0, 5), CV = TRUE, cv.method = "pretreatment", tol = 1e-5, max.iteration = 5000, ...))
}

test_that("X0: the H1 replication-14 fixture matches the regenerated panel; serial = 'none' gives -0.3283", {
  expect_identical(nrow(H1), 633L); expect_identical(length(unique(H1$id)), 35L); expect_identical(length(unique(H1$time)), 20L)
  expect_equal(sum(H1$d), 23)
  expect_close(sum(H1$y), h1$provenance$sum_y, 1e-9, "sum(y) equals the regenerated panel's")
  expect_close(sum(H1$y), 247.936024812346, 1e-9, "sum(y) pinned")
  n <- h1_fit("none")
  expect_equal(round(n$value$att.avg, 4), -0.3283)
})

test_that("X1: the cycle on real data, with the two hooks reproducing the earlier behaviour", {
  withr::local_options(fect.ar1.exclude.nonconverged = FALSE, fect.ar1.mask.zero = TRUE)
  res <- h1_fit(); f <- res$value; a <- f$ar1
  expect_identical(as.integer(a$path$round), 0:2); expect_identical(as.integer(a$path$r), c(2L, 3L, 2L))
  expect_identical(a$path$x[2], 0); expect_close(a$path$x[3], 0.0063, 1e-3, "round 2 x")
  expect_true(a$cycled); expect_false(a$converged)
  expect_identical(as.integer(a$cycle$round), 1:2); expect_identical(a$cycle$chosen, c(FALSE, TRUE))
  expect_identical(as.integer(f$r.cv), 2L); expect_close(a$x, 0.0063, 1e-3, "reported x")
  expect_close(f$att.avg, -0.4370, 1e-3, "att.avg on the masked panel")
  expect_equal(sum(grepl("cycle", res$warnings)), 1)
  expect_true(all(grepl('serial = "ar1"', res$warnings[grepl("cycle", res$warnings)], fixed = TRUE)))
  expect_true(a$em.converged)
  ## the cap on rounds does not change the reported pair (as long as the cap lets the cycle be seen)
  for (K in c(2, 4, 5, 6)) {
    withr::local_options(fect.ar1.rounds = K)
    g <- h1_fit()$value
    expect_close(g$att.avg, f$att.avg, 1e-10, sprintf("cap %d att.avg", K))
    expect_close(g$ar1$x, a$x, 1e-10, sprintf("cap %d x", K)); expect_identical(g$r.cv, f$r.cv)
    expect_true(g$ar1$cycled, info = paste("cap", K))
  }
  withr::local_options(fect.ar1.rounds = 1)
  g1 <- h1_fit()$value
  expect_false(g1$ar1$converged); expect_false(g1$ar1$cycled)
})

test_that("X2: the runaway rank-3 fit is passed over", {
  withr::local_options(fect.ar1.mask.zero = TRUE)
  res <- h1_fit(); f <- res$value; a <- f$ar1
  ex <- a$excluded
  expect_true(any(ex$round == 1 & ex$x == 0 & ex$r == 3), info = "round 1, x = 0, r = 3 is excluded (its EM reached max.iteration)")
  expect_true(a$em.converged)
  expect_false(any(grepl("did not converge", res$warnings)))
  expect_false(a$cycled)
  expect_false(a$x == 0 && f$r.cv == 3)
  expect_gt(f$att.avg, -3); expect_lt(f$att.avg, 3)
  ar1_note(sprintf("X2 (informational): x %.6f r %d rounds %d att.avg %.4f (planner's expectation: x 0, r 2, -0.4318); excluded: %s",
                   a$x, f$r.cv, max(a$path$round), f$att.avg, paste0("(", ex$round, ",", ex$x, ",", ex$r, ")", collapse = " ")))
})

test_that("X3: the new defaults on the same panel", {
  res <- h1_fit(); f <- res$value; a <- f$ar1
  expect_identical(a$lag.mask, a$x != 0)
  expect_true(a$em.converged)
  expect_gt(f$att.avg, -3); expect_lt(f$att.avg, 3)
  ar1_note(sprintf("X3 (informational): x %.6f r %d rounds %d att.avg %.4f cycled %s; warnings: %s", a$x, f$r.cv, max(a$path$round),
                   f$att.avg, a$cycled, if (length(res$warnings)) paste(res$warnings, collapse = " | ") else "none"))
  if (a$x == 0) {
    g <- ar1_fit_levels(H1, f$r.cv, tol = 1e-5, max.iteration = 5000, min.T0 = a$min.T0)
    expect_close(f$att.avg, g$att.avg, 1e-8, "x = 0: att.avg equals gsynth on the panel without period 1 at r.cv")
  }
})

## N: non-converged inner fits (P4, max.iteration = 3, tol = 1e-10, x fixed at 0.5) -------------------------------------------------
test_that("N1: a fixed rank whose EM stops at max.iteration is reported and warned about once", {
  P4 <- ar1_P4(); fix_x(0.5)
  res <- ar1_collect(ar1_fit(P4, r = 2, max.iteration = 3)); f <- res$value
  expect_false(f$ar1$em.converged); expect_gt(f$ar1$niter, 3); expect_true(is.finite(f$att.avg))
  w <- res$warnings
  expect_length(w, 1L)                                                           # N4: exactly one warning reaches the caller
  expect_true(grepl('serial = "ar1"', w, fixed = TRUE) && grepl("did not converge", w, fixed = TRUE))
  res2 <- ar1_collect(ar1_fit(P4, r = 2))                                        # default max.iteration (5000)
  expect_true(res2$value$ar1$em.converged); expect_false(any(grepl("did not converge", res2$warnings)))
})

test_that("N2/N4: a rank range with non-converged ranks: ineligible ranks, one warning at most", {
  P4 <- ar1_P4(); fix_x(0.5)
  res <- ar1_collect(ar1_fit(P4, r = c(0, 3), CV = TRUE, cv.method = "pretreatment", max.iteration = 3)); f <- res$value; a <- f$ar1
  expect_gt(a$search.nonconverged, 0)
  last <- max(a$path$round)
  if (isTRUE(a$em.converged)) {
    ar1_note("N2 branch (a): the final pair is converged; r.cv ", f$r.cv, "; excluded ranks round 1: ", paste(a$excluded$r[a$excluded$round == 1], collapse = ","))
    expect_false(f$r.cv %in% a$excluded$r[a$excluded$round == last])
    if (last == 1) {
      scored <- f$cv.treated$table$r[is.finite(f$cv.treated$table$loo) & f$cv.treated$table$r >= 1]
      expect_true(all(scored %in% a$excluded$r[a$excluded$round == 1]), info = "every scored rank >= 1 of round 1 is listed")
    }
    expect_length(res$warnings, 0L)
  } else {
    ar1_note("N2 branch (b): every rank was ineligible")
    expect_true(last %in% a$fallback.rounds)
    expect_length(res$warnings, 1L)
    expect_true(grepl('serial = "ar1"', res$warnings, fixed = TRUE) && grepl("did not converge", res$warnings, fixed = TRUE))
  }
})

test_that("N3: with the exclusion hook off the pick is the plain rule's on the same filtered data", {
  P4 <- ar1_P4(); fix_x(0.5)
  withr::local_options(fect.ar1.exclude.nonconverged = FALSE)
  res <- ar1_collect(ar1_fit(P4, r = c(0, 3), CV = TRUE, cv.method = "pretreatment", max.iteration = 3)); f <- res$value
  expect_identical(nrow(f$ar1$excluded), 0L)
  plain <- suppressWarnings(ar1_fect(ar1_hand_filter(P4, 0.5), CV = TRUE, r = c(0, 3), cv.method = "pretreatment",
                                     max.iteration = 3, min.T0 = f$ar1$min.T0))
  expect_identical(as.integer(f$r.cv), as.integer(plain$cv.treated$r_loo))
  expect_close(f$cv.treated$table$loo, plain$cv.treated$table$loo, 1e-8, "rule scores equal the plain rule on the hand-filtered V", rel = TRUE)
})

test_that("N5: the plain path (serial = 'none') is untouched by the eligibility rule", {
  P4 <- ar1_P4(); ref3 <- readRDS(test_path("fixtures", "ar1-reference.rds"))$rules$P4_maxit3
  expect_false(is.null(ref3))
  res <- ar1_collect(ar1_fect(P4, CV = TRUE, r = c(0, 5), cv.method = "pretreatment", max.iteration = 3)); f <- res$value; ct <- f$cv.treated
  expect_identical(ct$table$r, ref3$table$r)
  expect_close(ct$table$loo, ref3$table$loo, 1e-10, "N5 loo scores equal the prototype's treated_rank on fits from the same call", rel = TRUE)
  expect_close(ct$table$fwd, ref3$table$fwd, 1e-10, "N5 fwd scores", rel = TRUE)
  expect_identical(as.integer(ct$r_loo), as.integer(ref3$r_loo)); expect_identical(as.integer(f$r.cv), as.integer(ref3$r_loo))
  expect_identical(as.integer(ct$r_fwd), as.integer(ref3$r_fwd)); expect_identical(as.integer(ct$rbar), as.integer(ref3$rbar))
  expect_false(any(grepl('serial = "ar1"', res$warnings, fixed = TRUE)))
  expect_true(is.null(f$ar1))
})
