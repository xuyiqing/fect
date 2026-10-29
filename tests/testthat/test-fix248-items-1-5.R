## ---------------------------------------------------------------
## fect 2.4.8, fix items 1-5 (run 2026-10-28-fect-248-assembly).
##
## 1. carryover.rm removes the marked cells from the fit (before, they
##    stayed in as untreated controls, so nothing but carry.att changed);
##    method = "both" runs at defaults; plot(type = "hte") without a
##    covariate; carryoverTest = TRUE without carryover.period.
## 2. boot.diag gains sd.trim and n.dropped; the warning text says what
##    not to report and how to cross-check.
## 4. The breakdown warning also reads the first-step simulated-error
##    pool of the parametric bootstrap (cause "pool" / "both").
## 3, 5. The FE loop with neither factors nor covariates (fe_ad_iter)
##    obeys max.iteration (before, 501 rounds), and every EM path warns
##    at the cap.
##
## Helpers call non-base functions with their namespace (testthat
## parallel workers source files before stats/utils are attached).
## ---------------------------------------------------------------

suppressWarnings(utils::data("simdata", package = "fect"))
suppressWarnings(utils::data("simgsynth", package = "fect"))

## Cells that are not treated and lie 1..k periods after a treated period
## of the same unit (by time value; simdata has step 1).
mark_cells <- function(df, k = 2) {
  out <- logical(nrow(df))
  for (u in unique(df$id)) {
    pos <- which(df$id == u)
    tr <- df$time[pos][df$D[pos] == 1]
    cr <- setdiff(unique(unlist(lapply(seq_len(k), function(j) tr + j))), tr)
    out[pos] <- df$time[pos] %in% cr
  }
  out
}

## A T x N logical matrix (dimnames of obs.missing) TRUE on the flagged rows.
cells_to_mat <- function(df, flag, ref) {
  m <- matrix(FALSE, nrow(ref), ncol(ref), dimnames = dimnames(ref))
  m[cbind(as.character(df$time[flag]), as.character(df$id[flag]))] <- TRUE
  m
}

quiet_fit <- function(..., data = simdata, formula = Y ~ D + X1 + X2) {
  suppressWarnings(suppressMessages(fect::fect(
    formula, data = data, index = c("id", "time"), parallel = FALSE, ...
  )))
}

## all warning messages of an expression, muffled
collect_warnings <- function(expr) {
  ws <- character(0)
  val <- withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = val, warnings = ws)
}

## absolute-difference check (testthat 3e's expect_equal has no `scale`)
expect_close <- function(object, expected, tol, label = NULL) {
  d <- max(abs(as.numeric(object) - as.numeric(expected)))
  expect_true(length(object) == length(expected) && is.finite(d) && d <= tol,
              label = paste0(if (is.null(label)) "" else paste0(label, ": "),
                             "max abs difference ", format(d), " <= ", format(tol)))
}

first_units <- function(df, n) df[df$id %in% unique(df$id)[seq_len(n)], ]

mk <- mark_cells(simdata, 2)

## =============================================================
## 1. carryover.rm
## =============================================================

test_that("helper: k = 2 marks 237 cells of simdata", {
  expect_identical(sum(mk), 237L)
})

## -- R7 validation (fast: stops before any fit) --
test_that("R7: invalid carryover.rm values stop with a message naming it", {
  for (bad in list("a", c(1, 2), 2.5, -1, NA_real_, Inf)) {
    expect_error(
      suppressMessages(fect::fect(Y ~ D, data = simdata, index = c("id", "time"),
                                  method = "fe", se = FALSE, carryover.rm = bad)),
      "carryover.rm", info = paste(deparse(bad), collapse = "")
    )
  }
})

test_that("R1: carryover.rm changes the estimates (ife, fe, cfe, mc)", {
  skip_on_cran()
  ife0 <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE)
  ife2 <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE,
                    carryover.rm = 2)
  expect_close(ife2$att.avg, 3.9295, tol = 1e-3)
  expect_close(ife0$att.avg, 3.3886, tol = 1e-3)
  expect_equal(as.numeric(ife2$carry.time), c(1, 2))
  expect_close(as.numeric(ife2$carry.att), c(-0.084, 0.079), tol = 1e-3)
  expect_null(ife0$carry.att)

  fe0 <- quiet_fit(method = "fe", se = FALSE)
  fe2 <- quiet_fit(method = "fe", se = FALSE, carryover.rm = 2)
  expect_close(fe2$att.avg, 8.137, tol = 1e-3)
  expect_close(fe0$att.avg, 7.811, tol = 1e-3)

  cfe0 <- quiet_fit(method = "cfe", r = 1, CV = FALSE, se = FALSE)
  cfe2 <- quiet_fit(method = "cfe", r = 1, CV = FALSE, se = FALSE, carryover.rm = 2)
  expect_close(cfe2$att.avg, 3.938, tol = 1e-3)
  expect_close(cfe0$att.avg, 3.392, tol = 1e-3)

  mc0 <- quiet_fit(method = "mc", lambda = 0.05, CV = FALSE, se = FALSE)
  mc2 <- quiet_fit(method = "mc", lambda = 0.05, CV = FALSE, se = FALSE, carryover.rm = 2)
  expect_close(mc2$att.avg, 8.137, tol = 1e-3)
  expect_close(mc0$att.avg, 7.811, tol = 1e-3)
})

## -- R2 exclusion invariance: the marked cells take no part in the fit --
perturbed <- simdata
perturbed$Y[mk] <- perturbed$Y[mk] + 50

expect_invariant <- function(a, b, label) {
  expect_equal(a$att.avg, b$att.avg, tolerance = 1e-12, label = paste(label, "att.avg"))
  expect_equal(a$beta, b$beta, tolerance = 1e-12, label = paste(label, "beta"))
  expect_equal(a$att.avg.unit, b$att.avg.unit, tolerance = 1e-12, label = paste(label, "att.avg.unit"))
  expect_identical(a$r.cv, b$r.cv, label = paste(label, "r.cv"))
  expect_identical(a$lambda.cv, b$lambda.cv, label = paste(label, "lambda.cv"))
  expect_identical(a$method, b$method, label = paste(label, "method"))
  expect_close(as.numeric(b$carry.att - a$carry.att), rep(50, length(a$carry.att)), tol = 1e-8, label = paste(label, "carry.att shift"))
}

test_that("R2: adding 50 to Y on the marked cells moves nothing but carry.att (fixed r / lambda)", {
  skip_on_cran()
  specs <- list(
    ife = list(method = "ife", r = 1, CV = FALSE, force = "two-way"),
    fe = list(method = "fe"),
    cfe = list(method = "cfe", r = 1, CV = FALSE),
    mc = list(method = "mc", lambda = 0.05, CV = FALSE)
  )
  for (nm in names(specs)) {
    a <- do.call(quiet_fit, c(specs[[nm]], list(se = FALSE, carryover.rm = 2)))
    b <- do.call(quiet_fit, c(specs[[nm]], list(se = FALSE, carryover.rm = 2, data = perturbed)))
    expect_invariant(a, b, nm)
  }
})

test_that("R2: invariance holds through cross-validation and method = 'both'", {
  skip_on_cran()
  specs <- list(
    ife.cv = list(method = "ife", r = c(0, 3), CV = TRUE, k = 5, seed = 1),
    mc.cv = list(method = "mc", CV = TRUE, nlambda = 3, k = 3, seed = 2),
    both = list(method = "both", k = 5, seed = 1)
  )
  for (nm in names(specs)) {
    a <- do.call(quiet_fit, c(specs[[nm]], list(se = FALSE, carryover.rm = 2)))
    b <- do.call(quiet_fit, c(specs[[nm]], list(se = FALSE, carryover.rm = 2, data = perturbed)))
    expect_invariant(a, b, nm)
  }
})

test_that("R2: invariance holds for the bootstrap S.E.", {
  skip_on_cran()
  for (spec in list(list(method = "ife", r = 1, CV = FALSE, force = "two-way"),
                    list(method = "fe"))) {
    a <- do.call(quiet_fit, c(spec, list(se = TRUE, nboots = 15, seed = 3, carryover.rm = 2)))
    b <- do.call(quiet_fit, c(spec, list(se = TRUE, nboots = 15, seed = 3, carryover.rm = 2,
                                         data = perturbed)))
    expect_equal(a$est.avg[, "S.E."], b$est.avg[, "S.E."], tolerance = 1e-12, label = spec$method)
    expect_equal(a$att.avg, b$att.avg, tolerance = 1e-12)
  }
})

## -- R3 deleted-rows equivalence --
test_that("R3: carryover.rm = 2 equals fitting the panel without the marked rows", {
  skip_on_cran()
  specs <- list(
    ife = list(method = "ife", r = 1, CV = FALSE, force = "two-way"),
    fe = list(method = "fe"),
    cfe = list(method = "cfe", r = 1, CV = FALSE),
    ife.unit = list(method = "ife", r = 2, CV = FALSE, force = "unit")
  )
  for (nm in names(specs)) {
    a <- do.call(quiet_fit, c(specs[[nm]], list(se = FALSE, carryover.rm = 2,
                                                tol = 1e-10, max.iteration = 20000)))
    b <- do.call(quiet_fit, c(specs[[nm]], list(se = FALSE, data = simdata[!mk, ],
                                                tol = 1e-10, max.iteration = 20000)))
    expect_close(a$att.avg, b$att.avg, tol = 1e-6, label = paste(nm, "att.avg"))
    expect_close(as.numeric(a$beta), as.numeric(b$beta), tol = 1e-6,
                 label = paste(nm, "beta"))
    expect_close(as.numeric(a$att.avg.unit), as.numeric(b$att.avg.unit), tol = 1e-6, label = paste(nm, "att.avg.unit"))
  }
})

## -- R4 interactions --
test_that("R4: placebo test with carryover.rm codes 5 and 6", {
  skip_on_cran()
  f <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE,
                 carryover.rm = 2, placeboTest = TRUE, placebo.period = c(-2, 0))
  tab <- table(f$obs.missing)
  expect_true(all(c("5", "6") %in% names(tab)))
  expect_true(is.finite(f$att.placebo))
})

test_that("R4: carryover test counts from the first period after the removed window", {
  skip_on_cran()
  a <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE,
                 carryover.rm = 2, carryoverTest = TRUE, carryover.period = c(1, 3))
  b <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE,
                 carryoverTest = TRUE, carryover.period = c(1, 3))
  expect_true(is.finite(a$att.carryover))
  expect_false(isTRUE(all.equal(a$att.carryover, b$att.carryover)))
  expect_true(all(c("5", "6") %in% names(table(a$obs.missing))))
  ## carryover.period = 1 with a window of 2: the held-out cells are those
  ## exactly 3 periods after a treated period (and not within 2 of one)
  c1 <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE,
                  carryover.rm = 2, carryoverTest = TRUE, carryover.period = c(1, 1))
  expected <- cells_to_mat(simdata, mark_cells(simdata, 3) & !mk, c1$obs.missing)
  kept <- c1$obs.missing != 4L
  expect_identical((c1$obs.missing == 5L)[kept], expected[kept])
  expect_identical(c1$obs.missing == 6L,
                   cells_to_mat(simdata, mk, c1$obs.missing) & kept)
})

test_that("R4: jackknife, loo, balance.period and permute with carryover.rm", {
  skip_on_cran()
  j <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = TRUE,
                 vartype = "jackknife", carryover.rm = 2, data = first_units(simdata, 60))
  expect_true(is.finite(j$est.avg[1, "S.E."]))
  l <- NULL
  expect_no_error(l <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = first_units(simdata, 100), index = c("id", "time"),
    method = "ife", r = 1, CV = FALSE, force = "two-way", se = TRUE, nboots = 10,
    parallel = FALSE, seed = 1, loo = TRUE, carryover.rm = 2))))
  expect_s3_class(l, "fect")
  expect_error(suppressMessages(fect::fect(
    Y ~ D, data = simdata, index = c("id", "time"), method = "fe", se = FALSE,
    carryover.rm = 2, balance.period = c(-2, 2))))
  expect_error(suppressMessages(fect::fect(
    Y ~ D, data = simdata, index = c("id", "time"), method = "fe", se = TRUE,
    nboots = 10, parallel = FALSE, permute = TRUE, carryover.rm = 2)),
    "permute.*carryover.rm")
})

## -- R5 obs.missing and sample --
test_that("R5: code 6 marks the removed cells, and they are TRUE in sample", {
  skip_on_cran()
  f <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE,
                 carryover.rm = 2)
  removed <- f$obs.missing == 4L
  expect_identical(sum(f$obs.missing == 6L), 237L - sum(cells_to_mat(simdata, mk, f$obs.missing) & removed))
  expect_true(all(f$sample[f$obs.missing == 6L]))
  expect_identical(f$sample, matrix(f$obs.missing %in% c(1L, 2L, 5L, 6L),
                                    nrow(f$obs.missing), dimnames = dimnames(f$obs.missing)))
})

## -- R6 never-treated routes ignore carryover.rm with a warning --
test_that("R6: carryover.rm has no effect, with a warning, under nevertreated", {
  skip_on_cran()
  w <- collect_warnings(fect::fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                                   method = "ife", r = 0, CV = FALSE, se = FALSE,
                                   time.component.from = "nevertreated", carryover.rm = 2))
  expect_true(any(grepl("no effect", w$warnings) & grepl("nevertreated", w$warnings)))
  b <- quiet_fit(method = "ife", r = 0, CV = FALSE, se = FALSE,
                 time.component.from = "nevertreated")
  expect_identical(w$value$att.avg, b$att.avg)
  expect_identical(w$value$beta, b$beta)
  expect_null(w$value$carryover.rm)
  expect_null(w$value$carry.att)

  g <- collect_warnings(fect::fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
                                   method = "gsynth", r = 1, CV = FALSE, se = FALSE,
                                   carryover.rm = 2))
  expect_true(any(grepl("no effect", g$warnings) & grepl("nevertreated", g$warnings)))
  expect_null(g$value$carryover.rm)
})

## -- R7 values that run --
test_that("R7: carryover.rm = 0 equals none; 3L works; a k that marks nothing says so", {
  skip_on_cran()
  a <- quiet_fit(method = "fe", se = FALSE, carryover.rm = 0)
  b <- quiet_fit(method = "fe", se = FALSE)
  expect_identical(a$att.avg, b$att.avg)
  expect_null(a$carry.att)
  f3 <- quiet_fit(method = "fe", se = FALSE, carryover.rm = 3L)
  expect_equal(as.numeric(f3$carry.time), 1:3)
  ## units with no reversal: the window marks nothing
  norev <- vapply(split(simdata$D, simdata$id), function(d) all(diff(d) >= 0), logical(1))
  sub <- simdata[simdata$id %in% as.numeric(names(norev)[norev]), ]
  expect_true(any(sub$D == 1))
  expect_message(
    suppressWarnings(fect::fect(Y ~ D, data = sub, index = c("id", "time"),
                                method = "fe", se = FALSE, carryover.rm = 2)),
    "marked no cell"
  )
})

## -- R8 period position --
test_that("R8: the window counts periods by position, not by time value", {
  skip_on_cran()
  base <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE,
                    carryover.rm = 2)
  s5 <- simdata
  s5$time <- s5$time * 5
  f5 <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE,
                  carryover.rm = 2, data = s5)
  expect_close(f5$att.avg, base$att.avg, tol = 1e-8)
  expect_close(as.numeric(f5$carry.att), as.numeric(base$carry.att), tol = 1e-8)
  expect_identical(sum(f5$obs.missing == 6L), 237L)

  drop <- c(10, 500, 1200, 3000)
  u1 <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE,
                  carryover.rm = 2, data = simdata[-drop, ])
  u5 <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE,
                  carryover.rm = 2, data = s5[-drop, ])
  expect_close(u5$att.avg, u1$att.avg, tol = 1e-8)
  expect_close(as.numeric(u5$carry.att), as.numeric(u1$carry.att), tol = 1e-8)
})

## -- R9 step-1 balanced panel: exactly the helper's cells --
test_that("R9: on a step-1 balanced panel the window is exactly the helper's cells", {
  skip_on_cran()
  f <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE,
                 carryover.rm = 2)
  expect_identical(f$obs.missing == 6L, cells_to_mat(simdata, mk, f$obs.missing))
  ## a control cell just outside the window does move the fit
  mk3 <- mark_cells(simdata, 3) & !mk
  p3 <- simdata
  p3$Y[mk3] <- p3$Y[mk3] + 50
  g <- quiet_fit(method = "ife", r = 1, CV = FALSE, force = "two-way", se = FALSE,
                 carryover.rm = 2, data = p3)
  expect_gt(abs(g$att.avg - f$att.avg), 1e-3)
})

## =============================================================
## 2. method = "both"
## =============================================================

test_that("B1: method = 'both' runs at defaults", {
  skip_on_cran()
  d80 <- first_units(simdata, 80)
  f <- NULL
  expect_no_error(f <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = d80, index = c("id", "time"), method = "both"))))
  expect_s3_class(f, "fect")
  expect_identical(f$method, "ife")
  expect_identical(as.integer(f$r.cv), 2L)
  expect_close(f$att.avg, 2.33941, tol = 1e-5)
  g <- quiet_fit(method = "both", r = 0:2, k = 3, se = FALSE, data = d80)
  expect_s3_class(g, "fect")
  h <- quiet_fit(method = "both", se = TRUE, nboots = 10, data = d80)
  expect_true(is.finite(h$est.avg[1, "S.E."]))
})

## =============================================================
## 3. small bugs
## =============================================================

test_that("S1: plot(type = 'hte') without covariate", {
  skip_on_cran()
  f1 <- quiet_fit(method = "fe", se = FALSE, formula = Y ~ D + X1)
  expect_no_error(p <- plot(f1, type = "hte"))
  expect_no_error(plot(f1, type = "heterogeneous"))
  expect_no_error(plot(f1, type = "hte", covariate = "X1"))
  expect_error(plot(f1, type = "hte", covariate = "X9"), "valid covariate")
  f2 <- quiet_fit(method = "fe", se = FALSE)
  expect_error(plot(f2, type = "hte"), "Please provide a covariate")
  f0 <- quiet_fit(method = "fe", se = FALSE, formula = Y ~ D)
  e <- tryCatch(plot(f0, type = "hte"), error = function(e) conditionMessage(e))
  expect_match(e, "no covariates")
  expect_false(grepl("coercion to 'logical\\(1\\)'", e))
})

test_that("S2: carryoverTest without carryover.period stops with a clear message", {
  expect_error(
    suppressMessages(fect::fect(Y ~ D, data = simdata, index = c("id", "time"),
                                method = "fe", se = FALSE, carryoverTest = TRUE)),
    "carryover.period"
  )
  e <- tryCatch(suppressMessages(fect::fect(Y ~ D, data = simdata, index = c("id", "time"),
                                            method = "fe", se = FALSE, carryoverTest = TRUE)),
                error = function(e) conditionMessage(e))
  expect_false(grepl("II.origin", e))
  expect_error(
    suppressMessages(fect::fect(Y ~ D, data = simdata, index = c("id", "time"),
                                method = "fe", se = FALSE, placeboTest = TRUE)),
    "placeboTest = TRUE needs \"placebo.period\""
  )
})

## =============================================================
## 4. bootstrap diagnostics
## =============================================================

test_that("D1: sd.trim and n.dropped on the final draws", {
  set.seed(1)
  b <- c(stats::rnorm(95), 1e5, 2e5, 3e5, 4e5, 5e5)
  bd <- fect:::.fect_boot_diag(b, nboots = 100)
  expect_identical(bd$n.extreme, 5L)
  expect_identical(bd$n.dropped, 5L)
  expect_equal(bd$sd.trim, stats::sd(b[1:95]), tolerance = 1e-12)
  expect_lt(abs(bd$sd.trim - 0.89), 0.05)
  expect_equal(bd$sd, stats::sd(b))
  expect_gt(bd$sd, 7e4)
  expect_true(bd$breakdown)
  expect_identical(bd$cause, "draws")
  expect_null(bd$sim.breakdown)

  set.seed(1)
  ok <- stats::rnorm(1000)
  bd <- fect:::.fect_boot_diag(ok, nboots = 1000)
  expect_identical(bd$n.dropped, 0L)
  expect_identical(bd$sd.trim, bd$sd)
  expect_identical(bd$cause, "none")

  set.seed(1)
  flat <- c(rep(1, 80), stats::rnorm(20))
  bd <- fect:::.fect_boot_diag(flat, nboots = 100)
  expect_true(is.na(bd$sd.trim))
  expect_true(is.na(bd$n.dropped))
  expect_true(is.na(bd$n.extreme))
  expect_false(bd$breakdown)
  expect_identical(bd$cause, "none")
})

test_that("D3: the final-draw warning text", {
  set.seed(1)
  bd <- fect:::.fect_boot_diag(c(stats::rnorm(95), 1e5, 2e5, 3e5, 4e5, 5e5), nboots = 100)
  m <- fect:::.fect_boot_diag_message(bd)
  for (p in c("bootstrap draws of the average effect break down", "should not be reported",
              "by period", "jackknife", "vartype = \"bootstrap\"", "sd.trim")) {
    expect_true(grepl(p, m, fixed = TRUE), info = p)
  }
  expect_false(grepl("simulated-error draws of the first step", m, fixed = TRUE))
})

## first-step pool
pool_wild <- function(seed = 2, B = 200) {
  set.seed(seed)
  e <- array(stats::rnorm(10 * 3 * B), dim = c(10, 3, B))
  e[, , 7] <- e[, , 7] * 1e6
  e[, , 5] <- NA
  e
}

test_that("P1: pool fields of .fect_para_sim_check", {
  sc <- fect:::.fect_para_sim_check(pool_wild())
  expect_identical(sc$n.sim, 200L)
  expect_identical(sc$n.sim.failed, 1L)
  expect_identical(sc$n.sim.used, 199L)
  expect_identical(sc$n.sim.extreme, 1L)
  expect_true(is.finite(sc$sim.sd.ratio))
  expect_gt(sc$sim.sd.ratio, 5)
  m <- sc$sim.max.error[is.finite(sc$sim.max.error)]
  expect_equal(sc$sim.sd, stats::sd(m))
  expect_equal(sc$sim.iqr.sd, stats::IQR(m) / 1.349)
})

test_that("P2: the pool check needs 20 usable draws", {
  set.seed(3)
  e <- array(stats::rnorm(10 * 3 * 15), dim = c(10, 3, 15))
  e[, , 2] <- e[, , 2] * 1e6
  sc <- fect:::.fect_para_sim_check(e)
  expect_identical(sc$n.sim.used, 15L)
  set.seed(1)
  bd <- fect:::.fect_boot_diag(stats::rnorm(200), nboots = 200, sim = sc, para.error = "ar")
  expect_false(bd$sim.breakdown)
  expect_identical(bd$cause, "none")
  ## and with 185 failed draws out of 200
  e <- pool_wild()
  e[, , 21:200] <- NA
  sc <- fect:::.fect_para_sim_check(e)
  expect_lt(sc$n.sim.used, 20L)
  bd <- fect:::.fect_boot_diag(stats::rnorm(200), nboots = 200, sim = sc, para.error = "ar")
  expect_false(bd$sim.breakdown)
  expect_identical(bd$cause, "none")
})

test_that("P3: cause 'pool', 'both', and the ratio-only trigger", {
  sc <- fect:::.fect_para_sim_check(pool_wild())
  set.seed(1)
  ok <- stats::rnorm(200)
  bd <- fect:::.fect_boot_diag(ok, nboots = 200, sim = sc, para.error = "ar")
  expect_true(bd$sim.breakdown)
  expect_false(bd$breakdown)
  expect_identical(bd$cause, "pool")
  m <- fect:::.fect_boot_diag_message(bd, r.cv = 2, has.missing = TRUE)
  for (p in c("simulated-error draws of the first step", "break down", "do not show it",
              "error covariance", "should not be reported", "by period", "sd.trim")) {
    expect_true(grepl(p, m, fixed = TRUE), info = p)
  }
  expect_false(grepl("bootstrap draws of the average effect break down", m, fixed = TRUE))
  expect_false(grepl("para.error = \"empirical\"", m, fixed = TRUE))  # panel has missing cells

  bd$para.error <- "empirical"
  m <- fect:::.fect_boot_diag_message(bd, r.cv = 2, has.missing = FALSE)
  expect_true(grepl("reused by every replicate", m, fixed = TRUE))
  bd$para.error <- "wild"
  m <- fect:::.fect_boot_diag_message(bd, r.cv = 2, has.missing = FALSE)
  expect_true(grepl("para.error = \"empirical\"", m, fixed = TRUE))

  bad <- ok
  bad[1:5] <- 1e6
  bd <- fect:::.fect_boot_diag(bad, nboots = 200, sim = sc, para.error = "ar")
  expect_identical(bd$cause, "both")
  m <- fect:::.fect_boot_diag_message(bd, r.cv = 2)
  expect_true(grepl("bootstrap draws of the average effect break down", m, fixed = TRUE))
  expect_true(grepl("simulated-error draws of the first step", m, fixed = TRUE))
  expect_false(grepl("do not show it", m, fixed = TRUE))

  ## ratio only: no draw beyond 20 IQRs, SD > 5 * IQR/1.349
  vals <- c(seq(1, 2, length.out = 186), rep(11, 14))
  e <- array(0.5, dim = c(10, 3, 200))
  for (j in 1:200) e[, , j] <- e[, , j] * vals[j] / 0.5 * c(1, rep(0.5, 29))
  sc <- fect:::.fect_para_sim_check(e)
  expect_equal(sc$sim.max.error, vals)
  expect_identical(sc$n.sim.extreme, 0L)
  expect_gt(sc$sim.sd.ratio, 5)
  bd <- fect:::.fect_boot_diag(ok, nboots = 200, sim = sc, para.error = "ar")
  expect_true(bd$sim.breakdown)
  expect_identical(bd$cause, "pool")
  m <- fect:::.fect_boot_diag_message(bd, r.cv = 1)
  expect_true(grepl("times their IQR/1.349", m, fixed = TRUE))
})

gs_fit <- function(...) {
  fect::fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
             method = "gsynth", force = "two-way", r = 2, CV = FALSE, se = TRUE,
             vartype = "parametric", parallel = FALSE, seed = 2, ...)
}

test_that("D2/P3: fect() warns once on the pool and leaves the S.E. alone", {
  skip_on_cran()
  plain <- collect_warnings(gs_fit(nboots = 50))
  expect_identical(plain$value$boot.diag$cause, "none")
  expect_false(any(grepl("break down", plain$warnings)))
  expect_equal(plain$value$boot.diag$sd, unname(plain$value$est.avg[1, "S.E."]), tolerance = 1e-12)

  orig <- fect:::.fect_para_sim_check
  testthat::local_mocked_bindings(
    .fect_para_sim_check = function(error.tr) {
      sc <- orig(error.tr)
      sc$n.sim.extreme <- 2L
      sc
    },
    .package = "fect"
  )
  w <- collect_warnings(gs_fit(nboots = 50))
  expect_identical(w$value$boot.diag$cause, "pool")
  hits <- grepl("break down", w$warnings)
  expect_identical(sum(hits), 1L)
  expect_true(grepl("simulated-error draws of the first step", w$warnings[hits], fixed = TRUE))
  expect_true(grepl("reused by every replicate", w$warnings[hits], fixed = TRUE))
  expect_identical(w$value$est.avg, plain$value$est.avg)
  expect_equal(w$value$boot.diag$sd, unname(w$value$est.avg[1, "S.E."]), tolerance = 1e-12)

  ## both causes: one warning holding both sentences
  orig_bd <- fect:::.fect_boot_diag
  testthat::local_mocked_bindings(
    .fect_boot_diag = function(att.avg.boot, ...) {
      att.avg.boot[1:5] <- att.avg.boot[1:5] + 1e6
      orig_bd(att.avg.boot, ...)
    },
    .package = "fect"
  )
  w2 <- collect_warnings(gs_fit(nboots = 50))
  expect_identical(w2$value$boot.diag$cause, "both")
  hits <- grepl("break down", w2$warnings)
  expect_identical(sum(hits), 1L)
  expect_true(grepl("bootstrap draws of the average effect break down", w2$warnings[hits], fixed = TRUE))
  expect_true(grepl("simulated-error draws of the first step", w2$warnings[hits], fixed = TRUE))
})

test_that("P2: a parametric fit with nboots = 5 does not warn", {
  skip_on_cran()
  w <- collect_warnings(gs_fit(nboots = 5))
  expect_false(any(grepl("break down", w$warnings)))
  expect_identical(w$value$boot.diag$cause, "none")
})

test_that("P4: print() shows the block when either cause fired, summary() once", {
  skip_on_cran()
  out <- suppressWarnings(suppressMessages(gs_fit(nboots = 50)))
  header <- "Bootstrap draws of the average effect"
  n_hdr <- function(txt) sum(grepl(header, txt, fixed = TRUE))
  expect_identical(n_hdr(capture.output(print(out))), 0L)
  s <- capture.output(print(summary(out)))
  expect_identical(n_hdr(s), 1L)
  expect_true(any(grepl("largest error of each simulated-error draw: SD", s, fixed = TRUE)))

  flagged <- out
  flagged$boot.diag$n.sim.extreme <- 2L
  flagged$boot.diag$sim.breakdown <- TRUE
  flagged$boot.diag$cause <- "pool"
  p <- capture.output(print(flagged))
  expect_identical(n_hdr(p), 1L)
  expect_true(any(grepl("simulated-error draws (first step)", p, fixed = TRUE)))
  expect_true(any(grepl("first-step simulated errors have broken down", p, fixed = TRUE)))
  expect_identical(n_hdr(capture.output(print(summary(flagged)))), 1L)

  ## a fit saved before `cause` existed falls back to `breakdown`
  old <- out
  old$boot.diag$cause <- NULL
  expect_identical(n_hdr(capture.output(print(old))), 0L)
  old$boot.diag$breakdown <- TRUE
  expect_identical(n_hdr(capture.output(print(old))), 1L)
})

## =============================================================
## 5. FE loop cap and cap warnings
## =============================================================

test_that("C1: the no-factor no-covariate FE loop obeys max.iteration", {
  skip_on_cran()
  w <- collect_warnings(fect::fect(Y ~ D, data = simdata, index = c("id", "time"),
                                   method = "fe", se = FALSE, max.iteration = 3, tol = 1e-12))
  hits <- grepl("did not converge within max.iteration = 3", w$warnings, fixed = TRUE)
  expect_identical(sum(hits), 1L)
  expect_identical(as.integer(w$value$niter), 4L)
  expect_identical(w$value$convergence$iter.cap, 3L)
  expect_false(w$value$convergence$converged)
  f <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = simdata, index = c("id", "time"), method = "fe", se = FALSE,
    max.iteration = 2000L, tol = 1e-300)))
  expect_identical(f$convergence$iter.cap, 2000L)
})

make_chain <- function(N, w, seed = 1) {
  set.seed(seed)
  TT <- N + w
  rows <- do.call(rbind, lapply(1:N, function(i) data.frame(id = i, time = i:(i + w))))
  al <- stats::rnorm(N)
  xi <- stats::rnorm(TT)
  rows$Y <- al[rows$id] + xi[rows$time] + stats::rnorm(nrow(rows), 0, 0.3)
  rows$D <- 0L
  tr <- rows$id > N * 0.75 & rows$time >= rows$id + w - 1
  rows$D[tr] <- 1L
  rows$Y[tr] <- rows$Y[tr] + 2
  rows
}

test_that("C2: a binding fixture runs to max.iteration and warns", {
  skip_on_cran()
  cases <- list(list(N = 40, w = 3, att = 2.176713),
                list(N = 60, w = 3, att = 2.134410),
                list(N = 80, w = 2, att = 1.955957))
  for (cs in cases) {
    d <- make_chain(cs$N, cs$w)
    w <- collect_warnings(fect::fect(Y ~ D, data = d, index = c("id", "time"), method = "fe",
                                     se = FALSE, min.T0 = 1, tol = 1e-13))
    expect_identical(as.integer(w$value$niter), 5001L)
    expect_identical(sum(grepl("did not converge", w$warnings)), 1L)
    expect_close(w$value$att.avg, cs$att, tol = 1e-5)
    dflt <- suppressMessages(fect::fect(Y ~ D, data = d, index = c("id", "time"), method = "fe",
                                        se = FALSE, min.T0 = 1))
    expect_identical(as.integer(dflt$niter), 1L)
  }
})

test_that("C3: every EM path warns at the cap", {
  skip_on_cran()
  cap <- list(max.iteration = 3, tol = 1e-12, se = FALSE, parallel = FALSE)
  sd_ <- list(data = simdata, index = c("id", "time"))
  sg <- list(data = simgsynth, index = c("id", "time"))
  paths <- list(
    ife.r1.x = c(list(Y ~ D + X1 + X2, method = "ife", r = 1, CV = FALSE), sd_),
    ife.r0.x = c(list(Y ~ D + X1 + X2, method = "ife", r = 0, CV = FALSE), sd_),
    ife.r2 = c(list(Y ~ D, method = "ife", r = 2, CV = FALSE), sd_),
    fe = c(list(Y ~ D, method = "fe"), sd_),
    fe.x = c(list(Y ~ D + X1 + X2, method = "fe"), sd_),
    cfe.r1 = c(list(Y ~ D + X1 + X2, method = "cfe", r = 1, CV = FALSE), sd_),
    cfe.r0 = c(list(Y ~ D + X1 + X2, method = "cfe", r = 0, CV = FALSE), sd_),
    mc.x = c(list(Y ~ D + X1 + X2, method = "mc", lambda = 0.5, CV = FALSE), sd_),
    mc = c(list(Y ~ D, method = "mc", lambda = 0.5, CV = FALSE), sd_),
    mc.small = c(list(Y ~ D + X1 + X2, method = "mc", lambda = 0.001, CV = FALSE), sd_),
    gsynth = c(list(Y ~ D + X1 + X2, method = "gsynth", r = 1, CV = FALSE), sg),
    ife.nt = c(list(Y ~ D + X1 + X2, method = "ife", r = 1, CV = FALSE,
                    time.component.from = "nevertreated"), sg),
    cfe.nt = c(list(Y ~ D + X1 + X2, method = "cfe", r = 1, CV = FALSE,
                    time.component.from = "nevertreated"), sg),
    ife.cv = list(Y ~ D + X1 + X2, method = "ife", r = c(0, 2), CV = TRUE, k = 3, seed = 1,
                  data = first_units(simdata, 80), index = c("id", "time"))
  )
  for (nm in names(paths)) {
    w <- collect_warnings(do.call(fect::fect, c(paths[[nm]], cap)))
    expect_true(any(grepl("did not converge", w$warnings)), info = nm)
  }
})
