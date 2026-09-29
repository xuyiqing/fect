## ---------------------------------------------------------
## Normalization regression test
## Verifies that sigma2 (est.fect) is consistent with and
## without normalize=TRUE (the sigma2 normalization bug fix)
## ---------------------------------------------------------

test_that("sigma2 is consistent with and without normalization (FE)", {

  skip_on_cran()

    set.seed(6001)

    N <- 30
    TT <- 15
    T0 <- 10
    Ntr <- 10
    tau <- 3.0

    alpha_i <- rnorm(N, 0, 1)
    xi_t <- rnorm(TT, 0, 0.5)

    Y_vec <- numeric(N * TT)
    D_vec <- integer(N * TT)
    id_vec <- integer(N * TT)
    time_vec <- integer(N * TT)

    idx <- 1
    for (i in 1:N) {
        for (t in 1:TT) {
            treated <- (i <= Ntr) && (t > T0)
            D_vec[idx] <- as.integer(treated)
            eps <- rnorm(1, 0, 1)
            Y_vec[idx] <- alpha_i[i] + xi_t[t] + tau * D_vec[idx] + eps
            id_vec[idx] <- i
            time_vec[idx] <- t
            idx <- idx + 1
        }
    }

    simdf <- data.frame(
        id = id_vec,
        time = time_vec,
        Y = Y_vec,
        D = D_vec
    )

    out_raw <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "fe",
        force = "two-way",
        normalize = FALSE,
        se = FALSE,
        parallel = FALSE
    ))

    out_norm <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "fe",
        force = "two-way",
        normalize = TRUE,
        se = FALSE,
        parallel = FALSE
    ))

    ## ATT should be essentially the same
    expect_equal(out_raw$att.avg, out_norm$att.avg, tolerance = 0.01,
        label = "ATT should match with/without normalization")

    ## sigma2.fect should be close (within 10% relative error)
    ## This is the regression test for the normalization bug fix
    if (!is.null(out_raw$sigma2.fect) && !is.null(out_norm$sigma2.fect)) {
        rel_diff <- abs(out_raw$sigma2.fect - out_norm$sigma2.fect) /
            max(abs(out_raw$sigma2.fect), 1e-10)
        expect_lt(rel_diff, 0.1,
            label = paste0("sigma2.fect relative difference = ",
                           round(rel_diff, 6),
                           " (raw=", round(out_raw$sigma2.fect, 4),
                           ", norm=", round(out_norm$sigma2.fect, 4), ")"))
    }
})

test_that("sigma2 is consistent with and without normalization (IFE)", {

  skip_on_cran()

    set.seed(6002)

    N <- 30
    TT <- 15
    T0 <- 10
    Ntr <- 10
    tau <- 2.0

    alpha_i <- rnorm(N, 0, 1)
    xi_t <- rnorm(TT, 0, 0.5)
    lambda_i <- rnorm(N, 0, 0.5)
    f_t <- rnorm(TT, 0, 0.5)

    Y_vec <- numeric(N * TT)
    D_vec <- integer(N * TT)
    id_vec <- integer(N * TT)
    time_vec <- integer(N * TT)

    idx <- 1
    for (i in 1:N) {
        for (t in 1:TT) {
            treated <- (i <= Ntr) && (t > T0)
            D_vec[idx] <- as.integer(treated)
            eps <- rnorm(1, 0, 1)
            Y_vec[idx] <- alpha_i[i] + xi_t[t] +
                lambda_i[i] * f_t[t] +
                tau * D_vec[idx] + eps
            id_vec[idx] <- i
            time_vec[idx] <- t
            idx <- idx + 1
        }
    }

    simdf <- data.frame(
        id = id_vec,
        time = time_vec,
        Y = Y_vec,
        D = D_vec
    )

    out_raw <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "ife",
        r = 1,
        CV = FALSE,
        force = "two-way",
        normalize = FALSE,
        se = FALSE,
        parallel = FALSE
    ))

    out_norm <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "ife",
        r = 1,
        CV = FALSE,
        force = "two-way",
        normalize = TRUE,
        se = FALSE,
        parallel = FALSE
    ))

    ## ATT should be essentially the same
    expect_equal(out_raw$att.avg, out_norm$att.avg, tolerance = 0.05,
        label = "IFE ATT should match with/without normalization")

    ## sigma2.fect should be close
    if (!is.null(out_raw$sigma2.fect) && !is.null(out_norm$sigma2.fect)) {
        rel_diff <- abs(out_raw$sigma2.fect - out_norm$sigma2.fect) /
            max(abs(out_raw$sigma2.fect), 1e-10)
        expect_lt(rel_diff, 0.15,
            label = paste0("IFE sigma2.fect relative difference = ",
                           round(rel_diff, 6)))
    }
})


## ---------------------------------------------------------------
## normalize = TRUE must not change any result (fect #166).
##
## The fit divides the outcome and the covariates by sd(Y) before fitting
## and puts every result back on the outcome's scale. Every numeric slot of
## a fit with normalize = TRUE must equal the normalize = FALSE fit's slot
## to relative 1e-6. Before the fix: sigma2.fect of gsynth and CV = TRUE
## fits was multiplied by sd(Y) once instead of sd(Y)^2, so the default
## equivalence threshold 0.36 * sqrt(sigma2.fect) was too small by
## sqrt(sd(Y)); data.long$Y stayed divided by sd(Y) (panelview(fit) reads
## it); a given mc lambda penalized the normalized outcome, so the
## estimates moved; est$VNT, the gsynth IC, the CV tables and the mc
## lambda.cv, lambda.seq and eigen.all were on the wrong scale.
## ---------------------------------------------------------------

## every numeric slot of a fit, named by its path, recursively
.n166_num_slots <- function(x, prefix = "") {
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
      out <- c(out, .n166_num_slots(x[[i]], paste0(prefix, "$", nm)))
    }
  } else if (is.numeric(x)) {
    out[[prefix]] <- x
  }
  out
}

## the numeric slots on which two fits differ (relative to max(1, |value|)),
## each with its largest difference; character(0) when they agree
.n166_differing <- function(f0, f1, tol = 1e-6, ignore = character()) {
  a <- .n166_num_slots(f0)
  b <- .n166_num_slots(f1)
  bad <- character()
  for (nm in setdiff(union(names(a), names(b)), ignore)) {
    if (!nm %in% names(a) || !nm %in% names(b)) {
      bad <- c(bad, paste(nm, "(only in one fit)"))
      next
    }
    u <- as.numeric(a[[nm]])
    v <- as.numeric(b[[nm]])
    if (length(u) != length(v)) {
      bad <- c(bad, paste(nm, "(length)"))
      next
    }
    if (!identical(is.finite(u), is.finite(v)) || !identical(is.na(u), is.na(v))) {
      bad <- c(bad, paste(nm, "(NA or Inf pattern)"))
      next
    }
    ok <- is.finite(u)
    if (!any(ok)) next
    d <- max(abs(u[ok] - v[ok]) / pmax(1, abs(u[ok])))
    if (d > tol) bad <- c(bad, sprintf("%s (max rel diff %.3g)", nm, d))
  }
  bad
}

.n166_fit <- function(normalize, data = sim_gsynth, ...) {
  suppressWarnings(suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = data, index = c("id", "time"),
    force = "two-way", parallel = FALSE, normalize = normalize, ...
  )))
}

## one pair of fits; every numeric slot must agree, and the three headline
## slots are checked by name so a failure names them
.n166_expect_same <- function(..., ignore = character()) {
  f0 <- .n166_fit(FALSE, ...)
  f1 <- .n166_fit(TRUE, ...)
  expect_identical(.n166_differing(f0, f1, ignore = ignore), character(0))
  expect_equal(f1$att.avg, f0$att.avg, tolerance = 1e-8)
  expect_equal(f1$sigma2.fect, f0$sigma2.fect, tolerance = 1e-8)
  expect_equal(f1$tost.threshold, f0$tost.threshold, tolerance = 1e-8)
  ## data.long holds the input outcome (panelview(fit) draws it)
  expect_equal(range(f1$data.long$Y), range(f0$data.long$Y), tolerance = 1e-10)
  invisible(list(f0, f1))
}

test_that("normalize = TRUE: fe", {
  skip_on_cran()
  data("sim_gsynth", package = "fect")
  fits <- .n166_expect_same(method = "fe", se = FALSE)
  expect_equal(range(fits[[2]]$data.long$Y), range(sim_gsynth$Y), tolerance = 1e-10)
})

test_that("normalize = TRUE: ife and cfe with r = 2", {
  skip_on_cran()
  data("sim_gsynth", package = "fect")
  fits <- .n166_expect_same(method = "ife", r = 2, CV = FALSE, se = FALSE)
  ## VNT holds singular values of E E' / (N T): it scales with sd(Y)^2
  expect_equal(fits[[2]]$est$VNT, fits[[1]]$est$VNT, tolerance = 1e-8)
  fits <- .n166_expect_same(method = "cfe", r = 2, CV = FALSE, se = FALSE)
  expect_equal(fits[[2]]$est$VNT, fits[[1]]$est$VNT, tolerance = 1e-8)
})

test_that("normalize = TRUE: cfe with Z and with Q.type = \"linear\" (gamma, kappa)", {
  skip_on_cran()
  data("sim_gsynth", package = "fect")
  ## Z and Q are not divided by sd(Y): gamma and kappa must be put back on
  ## the outcome's scale (found in the review of #166)
  fits <- .n166_expect_same(method = "cfe", Z = "L1", CV = FALSE, se = FALSE)
  expect_equal(fits[[2]]$gamma[[1]], fits[[1]]$gamma[[1]], tolerance = 1e-8)
  fits <- .n166_expect_same(method = "cfe", Q.type = "linear", CV = FALSE, se = FALSE)
  expect_equal(fits[[2]]$kappa[[1]], fits[[1]]$kappa[[1]], tolerance = 1e-8)
})

test_that("normalize = TRUE: ife with CV = TRUE", {
  skip_on_cran()
  data("sim_gsynth", package = "fect")
  fits <- .n166_expect_same(method = "ife", r = c(0, 3), CV = TRUE, seed = 1, se = FALSE)
  expect_equal(fits[[2]]$r.cv, fits[[1]]$r.cv)
  expect_equal(fits[[2]]$CV.out.ife, fits[[1]]$CV.out.ife, tolerance = 1e-8)
})

test_that("normalize = TRUE: gsynth with CV = FALSE and CV = TRUE", {
  skip_on_cran()
  data("sim_gsynth", package = "fect")
  fits <- .n166_expect_same(method = "gsynth", r = 2, CV = FALSE, se = FALSE)
  ## the IC was not put back at all (off by 2 log(sd(Y)))
  expect_equal(fits[[2]]$IC, fits[[1]]$IC, tolerance = 1e-8)
  expect_equal(fits[[2]]$est$IC, fits[[1]]$est$IC, tolerance = 1e-8)
  expect_equal(fits[[2]]$est$VNT, fits[[1]]$est$VNT, tolerance = 1e-8)
  fits <- .n166_expect_same(method = "gsynth", r = c(0, 3), CV = TRUE, seed = 1, se = FALSE)
  expect_equal(fits[[2]]$r.cv, fits[[1]]$r.cv)
  expect_equal(fits[[2]]$CV.out, fits[[1]]$CV.out, tolerance = 1e-8)
})

test_that("normalize = TRUE: mc with a given lambda means the same model", {
  skip_on_cran()
  data("sim_gsynth", package = "fect")
  fits <- .n166_expect_same(method = "mc", lambda = 0.01, CV = FALSE, se = FALSE)
  ## before the fix the ATT moved from 5.118 to 5.085 (the two-way fixed
  ## effects ATT: the penalty on the normalized outcome removed the whole
  ## low-rank part)
  expect_equal(fits[[2]]$att.avg, fits[[1]]$att.avg, tolerance = 1e-8)
  expect_equal(fits[[2]]$eff, fits[[1]]$eff, tolerance = 1e-8)
  expect_equal(fits[[2]]$lambda.cv, 0.01)
  expect_equal(fits[[2]]$lambda.norm, fits[[1]]$lambda.norm, tolerance = 1e-8)
  expect_equal(fits[[2]]$eigen.all, fits[[1]]$eigen.all, tolerance = 1e-8)
})

test_that("normalize = TRUE: mc with CV = TRUE reports lambda on the outcome's scale", {
  skip_on_cran()
  data("sim_gsynth", package = "fect")
  fits <- .n166_expect_same(method = "mc", CV = TRUE, seed = 1, se = FALSE)
  expect_equal(fits[[2]]$lambda.cv, fits[[1]]$lambda.cv, tolerance = 1e-8)
  expect_equal(fits[[2]]$lambda.seq, fits[[1]]$lambda.seq, tolerance = 1e-8)
  expect_equal(fits[[2]]$eigen.all, fits[[1]]$eigen.all, tolerance = 1e-8)
  expect_equal(fits[[2]]$CV.out.mc, fits[[1]]$CV.out.mc, tolerance = 1e-8)
  ## refitting without normalize at the reported lambda.cv gives the CV fit
  ## (fect_cv's fits stop at cv_tol = max(tol, 1e-3), so the refit uses that
  ## tolerance; before the fix lambda.cv was on the normalized scale and the
  ## refit gave 5.385 instead of 5.238)
  refit <- .n166_fit(FALSE, method = "mc", lambda = fits[[2]]$lambda.cv, CV = FALSE, se = FALSE,
                     tol = 1e-3)
  expect_equal(refit$att.avg, fits[[2]]$att.avg, tolerance = 1e-6)
  ## a user-supplied grid is on the outcome's scale too
  grid <- c(0.05, 0.01, 0.002)
  fits <- .n166_expect_same(method = "mc", lambda = grid, CV = TRUE, seed = 1, se = FALSE)
  expect_equal(fits[[2]]$lambda.seq, grid)
  expect_true(fits[[2]]$lambda.cv %in% grid)
})

test_that("normalize = TRUE: the equivalence test of a parametric fit", {
  skip_on_cran()
  data("sim_gsynth", package = "fect")
  ## tol = 1e-10: the parametric replicates refit the balanced control panel,
  ## whose solver starts from feols() values that depend on the covariates'
  ## scale (normalize divides X by sd(Y)) and stops on an absolute change of
  ## beta; at the default tol the replicates differ by about 1e-5 between the
  ## two legs, at 1e-10 by about 1e-10
  fits <- .n166_expect_same(method = "gsynth", r = 2, CV = FALSE, se = TRUE,
                            nboots = 20, vartype = "parametric", seed = 1, tol = 1e-10)
  ## the default threshold 0.36 * sqrt(sigma2.fect): 0.489 both ways (it was
  ## 0.202 with normalize = TRUE, and the TOST p-value moved with it)
  expect_equal(fits[[2]]$tost.threshold, fits[[1]]$tost.threshold, tolerance = 1e-8)
  expect_equal(fits[[2]]$test.out$tost.equiv.p, fits[[1]]$test.out$tost.equiv.p,
               tolerance = 1e-6)
  expect_equal(fits[[2]]$est.avg, fits[[1]]$est.avg, tolerance = 1e-6)
})

test_that("normalize = TRUE: the outcome in hundredths (the threshold gap was 24-fold)", {
  skip_on_cran()
  data("sim_gsynth", package = "fect")
  d100 <- sim_gsynth
  d100$Y <- d100$Y * 100
  fits <- .n166_expect_same(data = d100, method = "gsynth", r = 2, CV = FALSE, se = FALSE)
  base <- .n166_fit(FALSE, method = "gsynth", r = 2, CV = FALSE, se = FALSE)
  ## before the fix: 48.9 vs 2.02
  expect_equal(fits[[2]]$tost.threshold, 100 * base$tost.threshold, tolerance = 1e-8)
  expect_equal(range(fits[[2]]$data.long$Y), range(d100$Y), tolerance = 1e-10)
})
