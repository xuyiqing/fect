## fect 2.5.0: the gsynth-revisit steps against the paper's replication code
## (github.com/xuyiqing/gsynth-revisit, private). These tests run only where
## that repository is checked out (GSYNTH_REVISIT_REPO, default
## ~/GitHub/gsynth-revisit) and never on CRAN. The replication code is read,
## never changed: each file is sourced into a local environment. The run over
## the paper's 95 cases (statsclaw-workspace/fect/runs/, results/match.md) is
## the full comparison; these tests keep a small version of it in the suite.

.gr_repo <- function() {
  path.expand(Sys.getenv("GSYNTH_REVISIT_REPO", "~/GitHub/gsynth-revisit"))
}
.skip_no_ref <- function(files) {
  skip_on_cran()
  repo <- .gr_repo()
  skip_if_not(all(file.exists(file.path(repo, files))), "gsynth-revisit replication code not available")
}
.ref_env <- function(files, level = "0.90") {
  old <- Sys.getenv("OF_LEVEL", unset = NA)
  Sys.setenv(OF_LEVEL = level)
  on.exit(if (is.na(old)) Sys.unsetenv("OF_LEVEL") else Sys.setenv(OF_LEVEL = old))
  env <- new.env(parent = globalenv())
  for (f in files) sys.source(file.path(.gr_repo(), f), envir = env)
  env
}
.rel <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b) || any(is.na(a) != is.na(b))) return(Inf)
  k <- !is.na(a)
  if (!any(k)) return(0)
  max(abs(a[k] - b[k]) / pmax(1, abs(b[k])))
}

data(simgsynth, package = "fect")
data(turnout, package = "fect")

## a staggered panel with gaps and many controls: 3 adoption cohorts, 260
## never-treated units (above the 200 held-out donors and the dense solver)
.big_panel <- function(seed = 9) {
  set.seed(seed)
  N <- 280; T <- 16; tr <- 1:20
  F <- matrix(rnorm(T * 2), T, 2); L <- matrix(rnorm(N * 2), N, 2)
  Y <- F %*% t(L) + outer(rep(1, T), rnorm(N)) + matrix(rnorm(T * N), T, N)
  onset <- rep(c(8, 10, 12), length.out = length(tr))
  D <- matrix(0, T, N)
  for (k in seq_along(tr)) D[onset[k]:T, tr[k]] <- 1
  d <- data.frame(id = rep(1:N, each = T), time = rep(1:T, N), Y = c(Y), D = c(D))
  d[-sample(nrow(d), 150), ]
}

test_that("overlap_check() = overlap_free_panel() at level 0.90", {
  skip_if_not_installed("digest"); skip_if_not_installed("osqp"); skip_if_not_installed("Matrix")
  .skip_no_ref("R/overlap_free.R")
  ref <- .ref_env("R/overlap_free.R")
  expect_identical(ref$OF_LEVEL, 0.90)
  for (d in list(simgsynth, .big_panel())) {
    w <- fect:::.gsc_wide(Y ~ D, d, c("id", "time"))
    a <- ref$overlap_free_panel(w$Y, w$D, ids = w$ids)
    b <- overlap_check(Y ~ D, data = d, index = c("id", "time"), level = 0.90)
    expect_identical(a$units$unit, b$units$unit)
    expect_identical(w$time[a$units$onset], b$units$onset)
    expect_identical(a$units$status, b$units$status)
    expect_lte(.rel(b$units$error, a$units$err), 1e-8)
    expect_lte(.rel(b$units$yardstick, a$units$q95), 1e-8)
    expect_identical(a$units$pct, b$units$pct)
    expect_identical(a$units$flag, b$units$flag)
    expect_lte(max(mapply(.rel, b$yard, a$yard)), 1e-8)
    expect_identical(a$cohorts$flag, b$cohorts$flag)
    expect_identical(a$cohorts$pct, b$cohorts$pct)
    expect_lte(.rel(b$cohorts$error, a$cohorts$err), 1e-8)
  }
})

test_that("fect_rank() scores = treated_rank() and window_minimum() on the same fits", {
  .skip_no_ref("R/rank_rule.R")
  ref <- .ref_env("R/rank_rule.R")
  cases <- list(
    list(Y ~ D + X1 + X2, simgsynth, c("id", "time"), list(force = "two-way")),
    list(turnout ~ policy_edr + policy_mail_in + policy_motor, turnout, c("abb", "year"), list(force = "two-way")),
    list(Y ~ D, .big_panel(), c("id", "time"), list(force = "unit"))
  )
  for (cs in cases) {
    rk <- do.call(fect_rank, c(list(formula = cs[[1]], data = cs[[2]], index = cs[[3]], r = 0:5,
                                    keep.fits = TRUE), cs[[4]]))
    a <- ref$treated_rank(rk$fits)
    b <- fect:::.rank_scores(rk$fits)
    expect_identical(b$table$r, a$table$r)
    expect_lte(.rel(b$table$fwd, a$table$fwd), 1e-8)
    expect_lte(.rel(b$table$loo, a$table$loo), 1e-8)
    expect_identical(b$r.fwd, a$r_fwd)
    expect_identical(b$r.loo, a$r_loo)
    expect_identical(c(b$rbar.fwd, b$rbar.loo, b$n.fwd, b$n.loo), c(a$rbar_fwd, a$rbar_loo, a$n_fwd, a$n_loo))
    wa <- ref$window_minimum(rk$fits[["0"]])
    expect_identical(rk$window$meets, wa$meets)
    expect_identical(rk$window$all8, wa$all8)
    expect_equal(rk$window$share8, wa$share8)
  }
})

test_that("loading_check() = overlap_hull(), and the extrapolation share", {
  skip_if_not_installed("Rglpk")
  .skip_no_ref("R/rank_rule.R")
  ref <- .ref_env("R/rank_rule.R")
  for (r in 1:3) {
    fit <- fect(Y ~ D + X1, data = simgsynth, index = c("id", "time"), method = "gsynth", r = r, se = FALSE)
    a <- ref$overlap_hull(fit)
    b <- loading_check(fit)
    expect_identical(a$table$in_hull, b$table$in.hull)
    expect_lte(.rel(b$table$d2, a$table$d2), 1e-6)
    expect_identical(a$table$flagged, b$table$flagged)
    expect_equal(b$Lco, a$Lco, tolerance = 1e-12)
  }
})

test_that("IFE-EM bootstrap = ife_noise_group_x() + ife_pb_x(): groups, draws, noise and intervals", {
  .skip_no_ref(c("sim/ife_em.R", "sim/ife_boot.R", "sim/ife_boot_x.R"))
  ref <- new.env(parent = globalenv())
  for (f in c("sim/ife_em.R", "sim/ife_boot.R", "sim/ife_boot_x.R"))
    sys.source(file.path(.gr_repo(), f), envir = ref)
  ref$IFE_TOL <- 1e-15; ref$IFE_MAXIT <- 500000L         # the reference converged well past fect's tol
  set.seed(12)
  T <- 20; N <- 43; Ntr <- 3; r <- 2; B <- 15
  F <- matrix(rnorm(T * r), T, r); L <- matrix(rnorm(N * r), N, r)
  X <- matrix(rnorm(T * N), T, N) + 0.5 * F %*% t(L)
  Y <- outer(rnorm(T), rep(1, N)) + outer(rep(1, T), rnorm(N)) + F %*% t(L) + 1.5 * X + matrix(rnorm(T * N), T, N)
  D <- matrix(0, T, N); D[16:T, 1:Ntr] <- 1
  d <- data.frame(id = rep(1:N, each = T), time = rep(1:T, N), Y = c(Y), D = c(D), X = c(X))
  d <- d[-sample(nrow(d), 60), ]
  fit <- suppressWarnings(suppressMessages(
    fect(Y ~ D + X, data = d, index = c("id", "time"), method = "ife", r = r, CV = FALSE, se = TRUE,
         vartype = "parametric", nboots = B, seed = 5, parallel = FALSE, keep.sims = TRUE,
         tol = 1e-12, max.iteration = 20000)))
  I <- fit$I.dat == 1
  Yr <- fit$Y.dat; Yr[!I] <- NA
  Dm <- fit$D.dat
  Xa <- array(0, c(T, N, 1))
  Xa[cbind(match(d$time, fit$rawtime), match(d$id, fit$id), 1)] <- d$X
  co <- which(colSums(Dm) == 0)
  f0 <- ref$ife_em_x(Yr, I & Dm == 0, r, Xa, "two-way", tol = 1e-15, maxit = 500000)
  pool <- ref$ife_noise_group_x(Yr, Dm, co, ref$noise_rank(r, length(co)), seed = 5, obs = I, X = Xa)
  e <- ref$ife_pb_x(Yr, Dm, r, B, pool, f0$fit, obs = I, X = Xa)
  ## the groups and draws, replayed in the reference's order
  set.seed(5)
  grp <- sample(rep_len(1:30, length(co)))
  idx <- t(vapply(seq_len(B), function(b) {
    k <- sample.int(ncol(pool), N, replace = TRUE)
    invisible(ref$patch_noise(pool[, k, drop = FALSE], pool, I))
    k
  }, integer(N)))
  expect_identical(fit$ife.boot$groups, grp)
  expect_identical(unname(fit$ife.boot$draws), idx)
  expect_identical(is.na(fit$ife.boot$noise), is.na(pool))
  expect_lte(.rel(fit$ife.boot$noise, pool), 1e-5)
  expect_lte(.rel(fit$att.avg, mean((Yr - f0$fit)[I & Dm == 1])), 1e-5)
  expect_lte(.rel(fit$att.avg.boot, e), 1e-5)
  expect_lte(.rel(fit$att.avg - quantile(fit$att.avg.boot, c(0.975, 0.025)), ref$basic_ci(fit$att.avg, e)), 1e-5)
})

test_that("overlap_check() on one of the paper's cases (needs the case data)", {
  skip_if_not_installed("digest"); skip_if_not_installed("osqp")
  .skip_no_ref(c("R/overlap_free.R", "R/paths.R", "papers/02_gilens2021/00_spec.R"))
  skip_if_not(dir.exists(Sys.getenv("GSYNTH_REVISIT_DATA", "")), "GSYNTH_REVISIT_DATA not set")
  ref <- .ref_env(c("R/paths.R", "R/overlap_free.R"))
  sp <- new.env(parent = ref); sp$ID <- "02_gilens2021"
  sys.source(file.path(.gr_repo(), "papers/02_gilens2021/00_spec.R"), envir = sp)
  cell <- sp$cells[[1]]
  s <- sp$spec
  for (k in intersect(names(cell), c("data", "treat", "covariates", "force", "min.T0", "index"))) s[[k]] <- cell[[k]]
  vars <- c(cell$y, s$treat, s$covariates, s$index)
  d <- s$data[, vars]; d <- d[stats::complete.cases(d), ]
  fm <- stats::as.formula(paste(cell$y, "~", paste(c(s$treat, s$covariates), collapse = " + ")))
  w <- fect:::.gsc_wide(fm, d, s$index)
  a <- ref$overlap_free_panel(w$Y, w$D, ids = w$ids)
  b <- overlap_check(fm, data = d, index = s$index, level = 0.90)
  expect_identical(a$units$flag, b$units$flag)
  expect_identical(a$units$pct, b$units$pct)
  expect_lte(.rel(b$units$error, a$units$err), 1e-8)
})
