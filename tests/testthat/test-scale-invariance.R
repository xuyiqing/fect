## ---------------------------------------------------------------
## Scale invariance of a fect fit (a permanent homogeneity check; fect
## #166).
##
## Multiplying the outcome by c = 100 must multiply every slot by its power
## of c: estimates, standard errors, intervals, effects, counterfactuals,
## residuals, fixed effects, loadings and the equivalence threshold by c;
## variances (sigma2, sigma2.fect, PC, est$VNT, att.vcov) and squared-error
## CV scores (MSPE, WMSPE, GMSPE, WGMSPE, MAD of the squared residuals,
## MSPTATT, MSE) by c^2; p-values, counts, indices, the chosen r, the
## chosen lambda's index, lambda.norm and the factors not at all;
## information criteria (IC) shift by 2 log(c). A given mc lambda scales by
## c, and so do lambda.cv, lambda.seq and eigen.all.
##
## Each fit is checked with normalize = FALSE and with normalize = TRUE.
## With normalize = TRUE the fit runs on Y / sd(Y), which is the same
## whatever the scale, so a slot put back with the wrong power of sd(Y)
## breaks homogeneity here even though the normalize = FALSE fit is exact
## by construction. Before the #166 fix that caught sigma2.fect (so
## tost.threshold), est$VNT, the gsynth IC, data.long$Y, the CV tables and
## the mc lambda slots.
##
## Precision: two things in the solvers depend on the scale without being
## results. The balanced-panel solver with covariates (beta_iter in
## src/ife_sub.cpp, used by gsynth for its control panel) stops on the
## absolute change of beta; and the starting values come from
## fixest::feols(), whose demeaning converges on the covariates' scale,
## which normalize = TRUE changes (X is divided by sd(Y)). At the default
## tol = 1e-5 the two scales can stop one EM iteration apart. So every
## fit here uses tol = 1e-10 (the estimates then agree to about 1e-9),
## the iteration counts are left out, and the CV tables, whose fold fits
## stop at cv_tol = max(tol, 1e-3), are checked at that precision, which
## still tells a wrong power (a factor 100 or 10000) from the right one.
## fect_cv() also returns the fit made in its loop at cv_tol, and the mc
## grid comes from the starting fit's residual singular values (feols
## precision, about 1e-7), so the estimates of a CV = TRUE mc fit carry
## about 1e-5 across scales: that block is checked at 1e-4. Rows a CV loop
## never visited are NA and are skipped.
## ---------------------------------------------------------------

## every numeric slot of a fit, named by its path; matrices are kept whole
## so that their columns can be classified
.si_slots <- function(x, prefix = "") {
  out <- list()
  if (is.data.frame(x)) {
    for (nm in names(x)) {
      if (is.numeric(x[[nm]])) out[[paste0(prefix, "$", nm)]] <- x[[nm]]
    }
  } else if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- rep("", length(x))
    for (i in seq_along(x)) {
      nm <- if (is.na(nms[i]) || nms[i] == "") paste0("[[", i, "]]") else nms[i]
      out <- c(out, .si_slots(x[[i]], paste0(prefix, "$", nm)))
    }
  } else if (is.numeric(x)) {
    out[[prefix]] <- x
  }
  out
}

## the power of c a slot scales with, by the slot's last name; "ic" is an
## additive shift by 2 log(c), "abs0"/"abs1" compare absolute values
## (factors and loadings are identified up to sign)
.si_leaf_power <- c(
  ## degree 2: variances and squared-error scores
  sigma2 = 2, sigma2.fect = 2, PC = 2, VNT = 2, att.vcov = 2,
  ## additive
  IC = "ic",
  ## degree 1: everything on the outcome's scale
  Y = 1, Y.dat = 1, Y.ct = 1, Y.ct.full = 1, Y.tr.cnt = 1, Y.ct.cnt = 1,
  eff = 1, eff.tr = 1, eff.equiv = 1, res = 1, res.full = 1, residuals = 1,
  fit = 1, mu = 1, alpha = 1, xi = 1, alpha.tr = 1, alpha.co = 1,
  beta = 1, att = 1, att.avg = 1, att.avg.unit = 1, att.avg.W = 1,
  att.avg.balance = 1, rmse = 1, tost.threshold = 1, att.bound = 1,
  att.avg.boot = 1, att.avg.unit.boot = 1, att.boot = 1,
  att.boot.original = 1, beta.boot = 1, eff.boot = 1,
  lambda.cv = 1, lambda.seq = 1, eigen.all = 1,
  lambda = "abs1", lambda.tr = "abs1", lambda.co = "abs1",
  factor = "abs0",
  ## degree 0: indicators, indices, counts, p-values, the covariates
  D = 0, D.dat = 0, I = 0, I.dat = 0, II = 0, T.on = 0, T.off = 0,
  N = 0, T = 0, Nco = 0, Ntr = 0, p = 0, tr = 0, co = 0, id = 0,
  rawtime = 0, time = 0, count = 0, N.calendar = 0, calendar.enp = 0,
  ci.alpha = 0, proportion = 0, pre.periods = 0, unit.type = 0,
  obs.missing = 0, hasRevs = 0, r.cv = 0, validF = 0, validX = 0,
  niter = 0, force = 0, lambda.norm = 0, wgt.implied = 0, X = 0,
  att.count.boot = 0, N_bar = 0, df1 = 0, df2 = 0, f.threshold = 0,
  tost.equiv.p = 0, f.p = 0, placebo.p = 0, placebo.equiv.p = 0,
  carryover.p = 0, carryover.equiv.p = 0
)

## the power of c a named column scales with (CV tables, estimate tables,
## eff.pre, pre.sd, eff.calendar)
.si_col_power <- c(
  sigma2 = 2, PC = 2, MSPE = 2, WMSPE = 2, GMSPE = 2, WGMSPE = 2, MAD = 2,
  MSPTATT = 2, MSE = 2,
  IC = "ic",
  Moment = 1, GMoment = 1, RMSE = 1, Bias = 1,
  ATT = 1, ATT.avg = 1, ATT.avg.unit = 1, `S.E.` = 1, CI.lower = 1,
  CI.upper = 1, Coef = 1, `ATT-calendar` = 1, `ATT-calendar Fitted` = 1,
  eff = 1, eff.equiv = 1, sd = 1,
  r = 0, lambda.norm = 0, count = 0, period = 0, unit = 0, p.value = 0,
  n.Treated = 0
)

.si_power <- function(path, col = NULL) {
  if (!is.null(col) && col %in% names(.si_col_power)) return(.si_col_power[[col]])
  leaf <- sub(".*\\$", "", path)
  if (grepl("^\\$(estCV|rmCV)\\$", path)) return(0)   # fold indices
  ## the CFE coefficients on Z (gamma) and loadings on Q (kappa) are lists
  ## of matrices on the outcome's scale (Z and Q are not divided by sd(Y))
  if (grepl("\\$(gamma|kappa)(\\$?\\[\\[[0-9]+\\]\\])?$", path)) return(1)
  if (leaf %in% names(.si_leaf_power)) return(.si_leaf_power[[leaf]])
  NA
}

## NULL when v equals u scaled by c^k (relative to max(1, |value|)), else
## a message; an unknown slot passes when it is homogeneous of degree 0,
## 1 or 2
.si_compare <- function(u, v, k, c, tol) {
  u <- as.numeric(u)
  v <- as.numeric(v)
  if (!identical(is.finite(u), is.finite(v))) return("(NA or Inf pattern)")
  ok <- is.finite(u)
  if (!any(ok)) return(NULL)
  u <- u[ok]
  v <- v[ok]
  dev <- function(kk) {
    if (identical(kk, "ic")) return(max(abs((v - u) - 2 * log(c))))
    if (identical(kk, "abs0")) return(max(abs(abs(v) - abs(u)) / pmax(1, abs(u))))
    if (identical(kk, "abs1")) return(max(abs(abs(v) - c * abs(u)) / pmax(1, c * abs(u))))
    kk <- as.numeric(kk)
    max(abs(v - u * c^kk) / pmax(1, abs(u * c^kk)))
  }
  if (is.na(k)) {
    d <- vapply(c(0, 1, 2), dev, numeric(1))
    if (min(d) <= tol) return(NULL)
    return(sprintf("(unknown slot; not homogeneous of degree 0, 1 or 2: %.3g, %.3g, %.3g)",
                   d[1], d[2], d[3]))
  }
  d <- dev(k)
  if (d <= tol) return(NULL)
  sprintf("(expected degree %s, max rel diff %.3g)", as.character(k), d)
}

## the slots of f100 (outcome times c) that do not scale as they should
## against f1; character(0) when all do
.si_differing <- function(f1, f100, c = 100, tol = 1e-6, ignore = character()) {
  a <- .si_slots(f1)
  b <- .si_slots(f100)
  bad <- character()
  for (nm in setdiff(union(names(a), names(b)), ignore)) {
    if (!nm %in% names(a) || !nm %in% names(b)) {
      bad <- c(bad, paste(nm, "(only in one fit)"))
      next
    }
    u <- a[[nm]]
    v <- b[[nm]]
    if (length(u) != length(v)) {
      bad <- c(bad, paste(nm, "(length)"))
      next
    }
    if (is.matrix(u) && !is.null(colnames(u)) &&
        any(colnames(u) %in% names(.si_col_power))) {
      for (j in seq_len(ncol(u))) {
        cn <- colnames(u)[j]
        msg <- .si_compare(u[, j], v[, j], .si_power(nm, cn), c, tol)
        if (!is.null(msg)) bad <- c(bad, sprintf("%s[, \"%s\"] %s", nm, cn, msg))
      }
    } else {
      msg <- .si_compare(u, v, .si_power(nm), c, tol)
      if (!is.null(msg)) bad <- c(bad, paste(nm, msg))
    }
  }
  bad
}

.si_fit <- function(data, ...) {
  suppressWarnings(suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = data, index = c("id", "time"),
    force = "two-way", parallel = FALSE, tol = 1e-10, ...
  )))
}

.si_cv_tables <- c("$CV.out", "$CV.out.ife", "$CV.out.mc")
.si_cv_tol <- 1e-3

## one specification, fitted on sim_gsynth and on sim_gsynth with the
## outcome times 100, with normalize = FALSE and TRUE; `scaled` names the
## arguments that are on the outcome's scale (a given lambda)
.si_expect_homogeneous <- function(..., scaled = character(),
                                   ignore = c("$niter", "$est$niter"),
                                   tol = 1e-6, c = 100) {
  data("sim_gsynth", package = "fect", envir = environment())
  d100 <- sim_gsynth
  d100$Y <- d100$Y * c
  args <- list(...)
  args100 <- args
  for (nm in scaled) args100[[nm]] <- args[[nm]] * c
  fits <- list()
  for (nz in c(FALSE, TRUE)) {
    f1 <- do.call(.si_fit, c(list(data = sim_gsynth, normalize = nz), args))
    f100 <- do.call(.si_fit, c(list(data = d100, normalize = nz), args100))
    bad <- .si_differing(f1, f100, c = c, tol = tol, ignore = c(ignore, .si_cv_tables))
    expect_identical(bad, character(0), label = paste0("normalize = ", nz, ": ", deparse(bad)))
    ## the CV tables at the fold fits' precision (see the note above)
    cv <- intersect(sub("^\\$", "", .si_cv_tables), names(f1))
    if (length(cv) > 0) {
      bad <- .si_differing(f1[cv], f100[cv], c = c, tol = .si_cv_tol)
      expect_identical(bad, character(0), label = paste0("CV tables, normalize = ", nz, ": ", deparse(bad)))
    }
    ## the headline slots by name, so a failure names them
    expect_equal(f100$att.avg, c * f1$att.avg, tolerance = 1e-8,
                 label = paste0("att.avg, normalize = ", nz))
    expect_equal(f100$sigma2.fect, c^2 * f1$sigma2.fect, tolerance = 1e-8,
                 label = paste0("sigma2.fect, normalize = ", nz))
    expect_equal(f100$tost.threshold, c * f1$tost.threshold, tolerance = 1e-8,
                 label = paste0("tost.threshold, normalize = ", nz))
    expect_equal(range(f100$data.long$Y), c * range(f1$data.long$Y), tolerance = 1e-10,
                 label = paste0("data.long$Y, normalize = ", nz))
    fits[[paste0("normalize.", nz)]] <- list(f1, f100)
  }
  invisible(fits)
}

test_that("scale invariance: fe", {
  skip_on_cran()
  fits <- .si_expect_homogeneous(method = "fe", se = FALSE)
  f <- fits[["normalize.TRUE"]]
  expect_equal(f[[2]]$IC, f[[1]]$IC + 2 * log(100), tolerance = 1e-8)
})

test_that("scale invariance: ife and cfe with r = 2", {
  skip_on_cran()
  fits <- .si_expect_homogeneous(method = "ife", r = 2, CV = FALSE, se = FALSE)
  f <- fits[["normalize.TRUE"]]
  expect_equal(f[[2]]$est$VNT, 100^2 * f[[1]]$est$VNT, tolerance = 1e-8)
  expect_equal(f[[2]]$sigma2, 100^2 * f[[1]]$sigma2, tolerance = 1e-8)
  expect_equal(f[[2]]$IC, f[[1]]$IC + 2 * log(100), tolerance = 1e-8)
  fits <- .si_expect_homogeneous(method = "cfe", r = 2, CV = FALSE, se = FALSE)
  f <- fits[["normalize.TRUE"]]
  expect_equal(f[[2]]$est$VNT, 100^2 * f[[1]]$est$VNT, tolerance = 1e-8)
})

test_that("scale invariance: cfe with Z and with Q.type = \"linear\" (gamma, kappa)", {
  skip_on_cran()
  ## Z and Q are not divided by sd(Y), so gamma and kappa come out of the
  ## solver on the Y / sd(Y) scale under normalize = TRUE; the fit must put
  ## them back (found in the review of #166). Tolerance 1e-4: the CFE
  ## solver with Z or Q stops on an absolute change, so where it stops
  ## depends on the scale (alpha, xi, gamma and kappa differ by about 1e-5
  ## relative between Y and 100 Y); a wrong power would be a factor 100.
  fits <- .si_expect_homogeneous(method = "cfe", Z = "L1", CV = FALSE, se = FALSE, tol = 1e-4)
  f <- fits[["normalize.TRUE"]]
  expect_equal(f[[2]]$gamma[[1]], 100 * f[[1]]$gamma[[1]], tolerance = 1e-4)
  expect_equal(f[[2]]$est$gamma[[1]], 100 * f[[1]]$est$gamma[[1]], tolerance = 1e-4)
  fits <- .si_expect_homogeneous(method = "cfe", Q.type = "linear", CV = FALSE, se = FALSE, tol = 1e-4)
  f <- fits[["normalize.TRUE"]]
  expect_equal(f[[2]]$kappa[[1]], 100 * f[[1]]$kappa[[1]], tolerance = 1e-4)
  expect_equal(f[[2]]$est$kappa[[1]], 100 * f[[1]]$est$kappa[[1]], tolerance = 1e-4)
})

test_that("scale invariance: ife with CV = TRUE", {
  skip_on_cran()
  fits <- .si_expect_homogeneous(method = "ife", r = c(0, 3), CV = TRUE, seed = 1, se = FALSE)
  for (f in fits) {
    expect_equal(f[[2]]$r.cv, f[[1]]$r.cv)
    ## IC differences across r do not depend on the scale
    expect_equal(diff(f[[2]]$CV.out.ife[, "IC"]), diff(f[[1]]$CV.out.ife[, "IC"]), tolerance = .si_cv_tol)
    expect_equal(f[[2]]$CV.out.ife[, "MSPE"], 100^2 * f[[1]]$CV.out.ife[, "MSPE"], tolerance = .si_cv_tol)
    expect_equal(f[[2]]$CV.out.ife[, "Moment"], 100 * f[[1]]$CV.out.ife[, "Moment"], tolerance = .si_cv_tol)
  }
})

test_that("scale invariance: gsynth with CV = FALSE and CV = TRUE", {
  skip_on_cran()
  fits <- .si_expect_homogeneous(method = "gsynth", r = 2, CV = FALSE, se = FALSE)
  f <- fits[["normalize.TRUE"]]
  expect_equal(f[[2]]$IC, f[[1]]$IC + 2 * log(100), tolerance = 1e-8)
  expect_equal(f[[2]]$est$VNT, 100^2 * f[[1]]$est$VNT, tolerance = 1e-8)
  fits <- .si_expect_homogeneous(method = "gsynth", r = c(0, 3), CV = TRUE, seed = 1, se = FALSE)
  for (f in fits) {
    expect_equal(f[[2]]$r.cv, f[[1]]$r.cv)
    expect_equal(diff(f[[2]]$CV.out[, "IC"]), diff(f[[1]]$CV.out[, "IC"]), tolerance = .si_cv_tol)
    expect_equal(f[[2]]$CV.out[, "MSPE"], 100^2 * f[[1]]$CV.out[, "MSPE"], tolerance = .si_cv_tol)
    expect_equal(f[[2]]$CV.out[, "Bias"], 100 * f[[1]]$CV.out[, "Bias"], tolerance = .si_cv_tol)
  }
})

test_that("scale invariance: mc with a given lambda and with CV = TRUE", {
  skip_on_cran()
  fits <- .si_expect_homogeneous(method = "mc", lambda = 0.01, CV = FALSE, se = FALSE,
                                 scaled = "lambda")
  ## eigen.all (the starting fit's residual singular values) has feols
  ## precision, about 1e-7
  for (f in fits) {
    expect_equal(f[[2]]$lambda.cv, 100 * f[[1]]$lambda.cv)
    expect_equal(f[[2]]$lambda.norm, f[[1]]$lambda.norm, tolerance = 1e-6)
    expect_equal(f[[2]]$eigen.all[1], 100 * f[[1]]$eigen.all[1], tolerance = 1e-6)
  }
  ## the estimates of a CV = TRUE mc fit carry about 1e-5 across scales
  ## (see the note above)
  fits <- .si_expect_homogeneous(method = "mc", CV = TRUE, seed = 1, se = FALSE, tol = 1e-4)
  for (f in fits) {
    expect_equal(which(f[[2]]$lambda.seq == f[[2]]$lambda.cv),
                 which(f[[1]]$lambda.seq == f[[1]]$lambda.cv))
    expect_equal(f[[2]]$lambda.cv, 100 * f[[1]]$lambda.cv, tolerance = 1e-6)
    expect_equal(f[[2]]$lambda.seq, 100 * f[[1]]$lambda.seq, tolerance = 1e-6)
    expect_equal(f[[2]]$lambda.norm, f[[1]]$lambda.norm, tolerance = 1e-6)
    expect_equal(f[[2]]$CV.out.mc[, "lambda.norm"], f[[1]]$CV.out.mc[, "lambda.norm"], tolerance = 1e-6)
    visited <- which(is.finite(f[[1]]$CV.out.mc[, "MSPE"]))
    expect_gte(length(visited), 3L)
    expect_equal(f[[2]]$CV.out.mc[visited, "MSPE"], 100^2 * f[[1]]$CV.out.mc[visited, "MSPE"],
                 tolerance = .si_cv_tol)
  }
})

test_that("scale invariance: a parametric fit (standard errors, intervals, p-values, the equivalence test)", {
  skip_on_cran()
  fits <- .si_expect_homogeneous(method = "gsynth", r = 2, CV = FALSE, se = TRUE, nboots = 20,
                                 vartype = "parametric", seed = 1)
  for (f in fits) {
    expect_equal(f[[2]]$est.avg[, "S.E."], 100 * f[[1]]$est.avg[, "S.E."], tolerance = 1e-6)
    expect_equal(f[[2]]$est.avg[, "p.value"], f[[1]]$est.avg[, "p.value"], tolerance = 1e-6)
    expect_equal(f[[2]]$est.att[, "CI.lower"], 100 * f[[1]]$est.att[, "CI.lower"], tolerance = 1e-6)
    expect_equal(f[[2]]$test.out$tost.equiv.p, f[[1]]$test.out$tost.equiv.p, tolerance = 1e-6)
    expect_equal(f[[2]]$tost.threshold, 100 * f[[1]]$tost.threshold, tolerance = 1e-8)
  }
})
