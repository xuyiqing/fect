## Phase 1: cv.method = "pretreatment" (the treated-unit rank rule, method = "gsynth").
##
## Contract (observable): fit$cv.treated is a list(table = data.frame(r, loo, fwd), r_loo,
## r_fwd, rbar, n_tr, n_loo, n_fwd, share_loo, share_fwd, note) with the meaning and values of the
## prototype's treated_rank() (xuyiqing/gsynth-revisit wip/template-v2-2026-10-16, e1eebae, R/rank_rule.R);
## fit$r.cv = r_loo and CV.out[, "MSPE"] holds the leave-one-out scores by rank (NA above rbar).
##
## Reference values: fixtures/ar1-reference.rds, written by fixtures/ar1-reference-gen.R, which runs
## the prototype's own treated_rank() on gsynth fits (CV = FALSE, r = 0..5, tol = 1e-10) of each
## test panel (helper-ar1.R). The prototype gate tests/check_rank_rule.R is run by hand (see the
## test-spec); it is not part of the package tests because it needs the prototype checkout.
##
## Byte-identical comparisons with the base build (R8) need a second library with fect 6c9c16c in
## FECT_BASE_LIB; they skip when it is not set.

skip_on_cran()

ref <- readRDS(test_path("fixtures", "ar1-reference.rds"))
P1 <- ar1_p1_long(ref)
pre_fit <- function(d, formula = y ~ d, ...) ar1_fect(d, formula, CV = TRUE, r = c(0, 5), cv.method = "pretreatment", ...)

expect_same_rule <- function(fit, rule, label) {
  ct <- fit$cv.treated
  expect_true(is.list(ct), info = label)
  expect_named(ct, c("table", "r_loo", "r_fwd", "rbar", "n_tr", "n_loo", "n_fwd", "share_loo", "share_fwd", "note"),
               ignore.order = TRUE)
  expect_identical(ct$table$r, rule$table$r)
  expect_close(ct$table$loo, rule$table$loo, 1e-10, paste(label, "loo scores"), rel = TRUE)
  expect_close(ct$table$fwd, rule$table$fwd, 1e-10, paste(label, "fwd scores"), rel = TRUE)
  expect_identical(as.integer(ct$r_loo), as.integer(rule$r_loo))
  expect_identical(as.integer(ct$r_fwd), as.integer(rule$r_fwd))
  expect_identical(as.integer(ct$rbar), as.integer(rule$rbar))
  expect_identical(as.integer(c(ct$n_tr, ct$n_loo, ct$n_fwd)), as.integer(c(rule$n_tr, rule$n_loo, rule$n_fwd)))
  expect_close(c(ct$share_loo, ct$share_fwd), c(rule$share_loo, rule$share_fwd), 1e-12, paste(label, "shares"))
  expect_identical(ct$note, rule$note)
  expect_identical(as.integer(fit$r.cv), as.integer(rule$r_loo))
  ## CV.out[, "MSPE"]: the loo scores by rank, NA above rbar
  cv <- fit$CV.out
  expect_close(cv[match(rule$table$r, cv[, "r"]), "MSPE"], rule$table$loo, 1e-10, paste(label, "CV.out MSPE"), rel = TRUE)
  above <- !(cv[, "r"] %in% rule$table$r)
  if (any(above)) expect_true(all(is.na(cv[above, "MSPE"])), info = label)
}

test_that("R2: scores and picks equal the prototype's treated_rank (balanced panels)", {
  f1 <- pre_fit(P1)
  expect_same_rule(f1, ref$rules$P1, "P1")
  ## pinned values (prototype treated_rank and simulation code rank_scores, which agree to 1e-15;
  ## test-spec amendment 4: 15 significant digits, relative tolerance 1e-8)
  pin_loo <- c(3.81734114141702, 2.21407309400062, 1.31595821309278, 1.26553950347434, 1.34476455369255, 1.44049830784233)
  pin_fwd <- c(2.84734647127381, 2.12808394805313, 1.57099346192948, 1.76512502816633, 1.89850059990273, 3.31582367958158)
  expect_close(f1$cv.treated$table$loo, pin_loo, 1e-8, "P1 pinned loo", rel = TRUE)
  expect_close(f1$cv.treated$table$fwd, pin_fwd, 1e-8, "P1 pinned fwd", rel = TRUE)
  expect_identical(f1$cv.treated$r_loo, 3L); expect_identical(f1$cv.treated$r_fwd, 2L)
  expect_same_rule(pre_fit(ar1_P2()), ref$rules$P2, "P2")
})

test_that("R3: scores equal the prototype's on staggered, unbalanced and covariate panels", {
  f3 <- pre_fit(ar1_P3()); expect_same_rule(f3, ref$rules$P3, "P3")
  f4 <- pre_fit(ar1_P4()); expect_same_rule(f4, ref$rules$P4, "P4")
  f5 <- pre_fit(ar1_P5(), y ~ d + x); expect_same_rule(f5, ref$rules$P5, "P5")
  for (f in list(f3, f4, f5)) {
    tab <- f$cv.treated$table
    expect_true(all(is.finite(tab$loo)))                       # finite for every rank up to rbar
    expect_identical(max(tab$r), f$cv.treated$rbar)            # and none above
    expect_true(all(is.na(f$CV.out[f$CV.out[, "r"] > f$cv.treated$rbar, "MSPE"])))
  }
})

test_that("R4: the forward window selects r.cv through the test hook", {
  withr::local_options(fect.cv.treated.window = "fwd")
  f <- pre_fit(P1)
  expect_identical(as.integer(f$r.cv), as.integer(ref$rules$P1$r_fwd))
  expect_identical(as.integer(f$r.cv), 2L)
})

test_that("R5: candidate ranks are fixed before scoring", {
  ## (a) one short treated unit (< 20% of the ATT weight) cannot support rank 5; the rest can.
  ##     8 units treated after t = 12 (8 post cells each) and one after t = 6 (14 post cells; 14/78 = 18%).
  d <- ar1_panel(30, 9, NA, 20, rho = 0.6, seed = 501, T0.unit = c(rep(12, 8), 6))
  f <- pre_fit(d)
  ct <- f$cv.treated
  expect_identical(ct$rbar, 5L)
  expect_identical(ct$n_tr, 9L); expect_identical(ct$n_loo, 8L)
  expect_gte(ct$share_loo, 0.8)
  expect_close(ct$share_loo, 64 / 78, 1e-12, "R5a share_loo = 64/78")        # pinned instance (test-spec R5)
  expect_true(all(is.finite(ct$table$loo))); expect_identical(ct$table$r, 0:5)
  expect_match(ct$note, "8 of 9")
  ## (b) units holding more than 20% of the weight are too short for rank 5: three units treated after
  ##     t = 6 (14 post cells each, 42 of 90 = 47%) support ranks up to 4 only.
  d2 <- ar1_panel(30, 9, NA, 20, rho = 0.6, seed = 502, T0.unit = c(rep(12, 6), rep(6, 3)))
  f2 <- pre_fit(d2)
  ct2 <- f2$cv.treated
  expect_identical(ct2$rbar, 4L)
  expect_identical(ct2$n_loo, 9L)
  expect_true(all(is.finite(ct2$table$loo))); expect_identical(max(ct2$table$r), 4L)
  cv <- f2$CV.out
  expect_true(all(is.na(cv[cv[, "r"] > 4, "MSPE"])))
  expect_true(all(is.finite(cv[cv[, "r"] <= 4, "MSPE"])))
  expect_lte(f2$r.cv, 4)
})

test_that("R6: panels too short to score anything stop with the argument in the message", {
  d <- ar1_panel(30, 5, 2, 8, rho = 0.6, seed = 503)
  expect_error(ar1_fect(d, CV = TRUE, r = c(2, 4), min.T0 = 2, cv.method = "pretreatment"),
               'cv.method = "pretreatment"', fixed = TRUE)
})

test_that("R7: refusals name cv.method = \"pretreatment\"; a fixed r runs silently", {
  P2 <- ar1_P2(); P2$w <- 1
  lit <- 'cv.method = "pretreatment"'
  chk <- function(...) expect_error(ar1_fect(P2, CV = TRUE, r = c(0, 5), cv.method = "pretreatment", ...), lit, fixed = TRUE)
  chk(method = "ife"); chk(method = "fe"); chk(method = "mc"); chk(method = "cfe")
  chk(W = "w"); chk(loading.bound = "simplex"); chk(cv.rule = "1se"); chk(cv.rule = "1pct"); chk(criterion = "pc")
  ## a fixed single r runs without cross-validation and without output
  expect_silent(f <- ar1_fect(P2, r = 2, cv.method = "pretreatment", quiet = FALSE))
  expect_identical(as.integer(f$r.cv), 2L)
})

test_that("R8b: \"pretreatment\" refuses elsewhere and never runs another cross-validation", {
  P2 <- ar1_P2(); lit <- 'cv.method = "pretreatment"'
  for (m in c("ife", "fe", "mc", "cfe", "both"))
    expect_error(ar1_fect(P2, method = m, CV = TRUE, r = c(0, 3), cv.method = "pretreatment"), lit, fixed = TRUE, info = m)
  expect_error(ar1_fect(P2, method = "ife", time.component.from = "nevertreated", CV = TRUE, r = c(0, 3),
                        cv.method = "pretreatment"), lit, fixed = TRUE)
  fit <- ar1_fect(P2, CV = TRUE, r = c(0, 3))
  expect_error(fect_mspe(fit, cv.method = "pretreatment"))
})

test_that("R8: cv.method = \"treated\", \"treated_units\" and \"loo\" are unchanged from the base build", {
  base_lib <- ar1_base_lib()
  skip_if(!nzchar(base_lib), "set FECT_BASE_LIB to a library with fect 6c9c16c")
  code <- '
    source(Sys.getenv("AR1_HELPER"))
    P2 <- ar1_P2()
    fits <- list()
    run <- function(nm, ...) {
      fits[[nm]] <<- tryCatch(ar1_fect(P2, CV = TRUE, r = c(0, 3), seed = 1, quiet = FALSE, ...), error = function(e) conditionMessage(e))
    }
    run("gs_treated", cv.method = "treated")
    run("gs_treated_units", cv.method = "treated_units")
    run("gs_loo", cv.method = "loo")
    run("gs_loo_se", cv.method = "loo", se = TRUE, vartype = "bootstrap", nboots = 20)
    run("ife_nt_treated", method = "ife", time.component.from = "nevertreated", cv.method = "treated")
    run("ife_nt_treated_units", method = "ife", time.component.from = "nevertreated", cv.method = "treated_units")
    out <- lapply(fits, function(f) if (is.character(f)) f else f[c("att.avg", "r.cv", "CV.out", "est.att", "est.avg")])
  '
  helper <- normalizePath(test_path("helper-ar1.R"))
  withr::local_envvar(AR1_HELPER = helper)
  run_here <- function() {
    msgs <- character(0)
    out <- NULL
    withCallingHandlers({ eval(parse(text = code), envir = environment()) },
                        message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
    list(out = out, msgs = msgs)
  }
  new <- run_here()
  old <- ar1_run_in_lib(base_lib, code)
  expect_identical(names(new$out), names(old$out))
  for (nm in names(old$out)) {
    if (is.character(old$out[[nm]])) {
      expect_identical(new$out[[nm]], old$out[[nm]], info = paste(nm, "error message"))
    } else {
      expect_equal(new$out[[nm]]$att.avg, old$out[[nm]]$att.avg, tolerance = 1e-12, info = nm)
      expect_identical(new$out[[nm]]$r.cv, old$out[[nm]]$r.cv, info = nm)
      expect_equal(new$out[[nm]]$CV.out, old$out[[nm]]$CV.out, tolerance = 1e-12, info = nm)
      expect_equal(new$out[[nm]]$est.att, old$out[[nm]]$est.att, tolerance = 1e-12, info = nm)
      expect_equal(new$out[[nm]]$est.avg, old$out[[nm]]$est.avg, tolerance = 1e-12, info = nm)
    }
  }
  expect_identical(new$msgs, old$msgs)                              # same messages, including the deprecation note
  expect_true(any(grepl("deprecated", old$msgs)))
})

test_that("R9: inference with the pretreatment rule", {
  P2 <- ar1_P2()
  f0 <- pre_fit(P2)
  f1 <- pre_fit(P2, se = TRUE, vartype = "parametric", nboots = 30, seed = 1)
  expect_false(is.null(f1$cv.treated))
  expect_identical(f1$r.cv, f0$r.cv)
  f2 <- pre_fit(P2, se = TRUE, vartype = "bootstrap", nboots = 30, seed = 1)
  expect_false(is.null(f2$cv.treated))
  expect_identical(f2$r.cv, f0$r.cv)
  expect_true(all(is.finite(f1$est.avg)) && all(is.finite(f2$est.avg)))
})

## R10 (amended): (a) the build's pick equals the prototype's leave-one-out pick on every one of the 40
## draws (the gate; picks stored in ar1-reference.rds$r10); (b) r.cv never above rbar; (c) r.cv >= 2 on at
## least 38 of 40 draws (the prototype itself: 40/40). The old requirement "share of r.cv == 2 >= 0.35" is
## withdrawn by the planner (the prototype's share is 7/40).
test_that("R10: on the motivating design the picks equal the prototype's", {
  picks <- vapply(1:40, function(s) {
    f <- pre_fit(ar1_mot_panel(s)); c(r = f$r.cv, rbar = f$cv.treated$rbar)
  }, c(r = 0, rbar = 0))
  expect_identical(as.integer(picks["r", ]), as.integer(ref$r10$proto_pick))   # (a) every draw
  expect_true(all(picks["r", ] <= picks["rbar", ]))                            # (b)
  expect_gte(sum(picks["r", ] >= 2), 38)                                       # (c)
})
