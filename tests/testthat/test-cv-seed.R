## `seed` controls the cross-validation fold draws (fixed in 2.4.7).
## Before the fix, fect.default() called set.seed(seed) only when se = TRUE or
## permute = TRUE. With se = FALSE the folds came from the session RNG, so two
## identical calls with the same seed could give different CV tables, and on
## some panels a different r (reported 2026-09-23). With se = FALSE,
## fect(..., seed = s) now draws the same folds as set.seed(s); fect(...).
## Without a seed the session RNG still drives the folds, and runs with
## se = TRUE or permute = TRUE are unchanged.

skip_on_cran()

## N = 30, T = 16, 10 treated units adopting at t = 11..13, two factors (the
## second one weak), one covariate. On the unfixed code, the session RNG state
## alone changed the chosen r on this panel: r.cv = 1 after set.seed(1) and
## 2 after set.seed(2), both with seed = 1.
.seed_panel <- function() {
  withr::with_seed(2, {
    N <- 30; TT <- 16; Ntr <- 10
    F <- matrix(rnorm(TT * 2), TT, 2)
    L <- matrix(rnorm(N * 2), N, 2); L[, 2] <- 0.5 * L[, 2]
    a <- rnorm(N); xi <- rnorm(TT); X <- matrix(rnorm(N * TT), TT, N)
    T0 <- rep(Inf, N); T0[seq_len(Ntr)] <- sample(11:13, Ntr, replace = TRUE)
    D <- sapply(seq_len(N), function(i) as.integer(seq_len(TT) >= T0[i]))
    Y <- outer(xi, rep(1, N)) + outer(rep(1, TT), a) + F %*% t(L) + 0.5 * X +
      2 * D + matrix(rnorm(N * TT), TT, N)
    data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
               Y = c(Y), D = c(D), X = c(X))
  })
}

## Set the session RNG to `session`, then fit with CV. Two calls with different
## `session` values start from different RNG states, so only `seed` can make
## their folds agree.
.cv_fit <- function(d, session, ...) {
  set.seed(session)
  suppressMessages(suppressWarnings(
    fect(Y ~ D + X, data = d, index = c("id", "time"), CV = TRUE, ...)
  ))
}

test_that("the same seed gives the same CV table and r, whatever the session RNG", {
  d <- .seed_panel()
  f1 <- .cv_fit(d, session = 1, method = "ife", r = c(0, 3), seed = 1)
  f2 <- .cv_fit(d, session = 2, method = "ife", r = c(0, 3), seed = 1)
  expect_identical(f1$CV.out.ife, f2$CV.out.ife)
  expect_identical(f1$r.cv, f2$r.cv)
})

test_that("seed = s draws the same folds as set.seed(s) before the call", {
  d <- .seed_panel()
  seeded <- .cv_fit(d, session = 99, method = "ife", r = c(0, 3), seed = 5)
  manual <- .cv_fit(d, session = 5, method = "ife", r = c(0, 3))
  expect_identical(seeded$CV.out.ife, manual$CV.out.ife)
  ## A different seed from the same session state gives different folds.
  other <- .cv_fit(d, session = 99, method = "ife", r = c(0, 3), seed = 6)
  expect_false(identical(seeded$CV.out.ife, other$CV.out.ife))
})

test_that("without a seed, the session RNG still drives the folds", {
  d <- .seed_panel()
  a <- .cv_fit(d, session = 3, method = "ife", r = c(0, 3))
  b <- .cv_fit(d, session = 3, method = "ife", r = c(0, 3))
  c <- .cv_fit(d, session = 4, method = "ife", r = c(0, 3))
  expect_identical(a$CV.out.ife, b$CV.out.ife)
  expect_false(identical(a$CV.out.ife, c$CV.out.ife))
})

test_that("seed also fixes block, gsynth and mc cross-validation", {
  d <- .seed_panel()
  ## gsynth runs its CV in fect_nevertreated() and stores the table in CV.out.
  cases <- list(
    block  = list(args = list(method = "ife", r = c(0, 3), cv.method = "block"),
                  table = "CV.out.ife", pick = "r.cv"),
    gsynth = list(args = list(method = "gsynth", r = c(0, 3)),
                  table = "CV.out", pick = "r.cv"),
    mc     = list(args = list(method = "mc"),
                  table = "CV.out.mc", pick = "lambda.cv")
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    f1 <- do.call(.cv_fit, c(list(d = d, session = 1, seed = 1, parallel = FALSE), cs$args))
    f2 <- do.call(.cv_fit, c(list(d = d, session = 2, seed = 1, parallel = FALSE), cs$args))
    expect_false(is.null(f1[[cs$table]]), info = nm)
    expect_identical(f1[[cs$table]], f2[[cs$table]], info = nm)
    expect_identical(f1[[cs$pick]], f2[[cs$pick]], info = nm)
  }
})

test_that("seed fixes the folds when CV runs in parallel", {
  d <- .seed_panel()
  f1 <- .cv_fit(d, session = 1, method = "ife", r = c(0, 3), seed = 1,
                parallel = "cv", cores = 2)
  f2 <- .cv_fit(d, session = 2, method = "ife", r = c(0, 3), seed = 1,
                parallel = "cv", cores = 2)
  expect_identical(f1$CV.out.ife, f2$CV.out.ife)
  expect_identical(f1$r.cv, f2$r.cv)
})

test_that("se = FALSE and a parallel bootstrap run pick the same folds for one seed", {
  ## The parallel bootstrap seeds with set.seed(seed) before its own CV, so the
  ## same seed gives the same CV table and r with or without standard errors.
  d <- .seed_panel()
  nose <- .cv_fit(d, session = 1, method = "ife", r = c(0, 3), seed = 7)
  boot <- .cv_fit(d, session = 2, method = "ife", r = c(0, 3), seed = 7,
                  se = TRUE, nboots = 10, parallel = TRUE, cores = 2)
  expect_identical(nose$CV.out.ife, boot$CV.out.ife)
  expect_identical(nose$r.cv, boot$r.cv)
})
