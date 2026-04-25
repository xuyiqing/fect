## Unit tests for .solve_bounded_loading and .cv_gamma_loading internals.
## See statsclaw-workspace/fect/runs/REQ-bounded-loadings/test-spec.md section S.

test_that("S.1: solver output lies on the simplex", {
    set.seed(1)
    T0 <- 20; r <- 2; Nco <- 15
    F_pre     <- matrix(rnorm(T0 * r), T0, r)
    Lambda_co <- matrix(rnorm(Nco * r), Nco, r)
    u_pre     <- F_pre %*% c(0.3, 0.5) + rnorm(T0, sd = 0.1)
    sol       <- fect:::.solve_bounded_loading(u_pre, F_pre, Lambda_co, gamma = 1)
    expect_equal(sum(sol$w), 1, tolerance = 1e-10)
    expect_true(all(sol$w >= -1e-12))
    expect_true(sol$converged)
    expect_equal(length(sol$lambda_hat), r)
})

test_that("S.2: lambda_hat = t(Lambda_co) %*% w", {
    set.seed(2)
    T0 <- 20; r <- 3; Nco <- 20
    F_pre     <- matrix(rnorm(T0 * r), T0, r)
    Lambda_co <- matrix(rnorm(Nco * r), Nco, r)
    u_pre     <- F_pre %*% c(0.2, 0.4, -0.1) + rnorm(T0, sd = 0.1)
    sol       <- fect:::.solve_bounded_loading(u_pre, F_pre, Lambda_co, gamma = 1)
    expect_equal(as.numeric(sol$lambda_hat),
                 as.numeric(crossprod(Lambda_co, sol$w)),
                 tolerance = 1e-8)
})

test_that("S.3: small gamma -> weights shrink toward uniform; large gamma -> near-unregularized", {
    set.seed(3)
    T0 <- 30; r <- 2; Nco <- 12
    F_pre     <- matrix(rnorm(T0 * r), T0, r)
    Lambda_co <- matrix(rnorm(Nco * r), Nco, r)
    u_pre     <- as.numeric(F_pre %*% c(0.4, 0.2))   # no noise; OLS target inside hull
    sol_small <- fect:::.solve_bounded_loading(u_pre, F_pre, Lambda_co, gamma = 1e-4)
    sol_large <- fect:::.solve_bounded_loading(u_pre, F_pre, Lambda_co, gamma = 1e6)
    ## Small gamma: entropy dominates, weights pulled very close to uniform
    expect_lt(max(abs(sol_small$w - 1 / Nco)), 0.05)
    ## Large gamma: entropy negligible, projection fits the noise-free target well
    resid_large <- u_pre - as.numeric(F_pre %*% sol_large$lambda_hat)
    expect_lt(sum(resid_large^2), 1e-4)
})

test_that("S.4: exterior OLS target -> projection on hull, nontrivial residual", {
    set.seed(4)
    T0 <- 20; r <- 2; Nco <- 10
    Lambda_co <- matrix(rnorm(Nco * r, sd = 0.5), Nco, r)
    F_pre     <- matrix(rnorm(T0 * r), T0, r)
    lambda_tr_true <- c(3, 3)              # far outside hull
    u_pre <- as.numeric(F_pre %*% lambda_tr_true) + rnorm(T0, sd = 0.1)
    sol   <- fect:::.solve_bounded_loading(u_pre, F_pre, Lambda_co, gamma = 1e6)
    ## Projected loading should lie within the bounding box of Lambda_co
    expect_lte(max(sol$lambda_hat), max(Lambda_co) + 1e-6)
    expect_gte(min(sol$lambda_hat), min(Lambda_co) - 1e-6)
    ## Residual should be substantially nonzero (cannot fit (3,3) with hull)
    resid <- u_pre - as.numeric(F_pre %*% sol$lambda_hat)
    expect_gt(sum(resid^2), 0.1)
})

test_that("S.5: Nco < 2 raises informative error", {
    expect_error(
        fect:::.solve_bounded_loading(
            u_pre = rnorm(10), F_pre = matrix(rnorm(10), 10, 1),
            Lambda_co = matrix(0.5, 1, 1), gamma = 1
        ),
        regexp = "at least 2 control units"
    )
})

test_that("S.6: all-zero Lambda_co returns uniform w, zero lambda_hat", {
    sol <- fect:::.solve_bounded_loading(
        u_pre     = rnorm(10),
        F_pre     = matrix(rnorm(20), 10, 2),
        Lambda_co = matrix(0, 15, 2),
        gamma     = 1
    )
    expect_equal(sol$w, rep(1 / 15, 15), tolerance = 1e-10)
    expect_equal(as.numeric(sol$lambda_hat), c(0, 0), tolerance = 1e-12)
    expect_true(sol$converged)
    expect_equal(sol$method, "degenerate-return-uniform")
})

test_that("S.7: mirror-descent fallback can run on ill-posed inputs", {
    set.seed(7)
    T0 <- 5; r <- 4; Nco <- 30       # r > T0 -- ill-posed regression
    F_pre     <- matrix(rnorm(T0 * r), T0, r)
    Lambda_co <- matrix(rnorm(Nco * r), Nco, r)
    u_pre     <- rnorm(T0)
    ## fallback = TRUE (default) ensures we don't silently fail
    sol <- fect:::.solve_bounded_loading(u_pre, F_pre, Lambda_co,
                                         gamma = 1, fallback = TRUE)
    expect_true(sol$converged)
    expect_true(sol$method %in% c("lbfgs", "mirror-descent"))
    expect_equal(sum(sol$w), 1, tolerance = 1e-6)
    expect_true(all(sol$w >= -1e-8))
})

test_that("S.8: .cv_gamma_loading picks from grid, returns valid selection", {
    set.seed(8)
    T0 <- 30; r <- 2; Nco <- 20; Ntr <- 5
    F_pre     <- matrix(rnorm(T0 * r), T0, r)
    Lambda_co <- matrix(rnorm(Nco * r), Nco, r)
    ## Treated loadings sampled from same dist; interior case
    Lambda_tr <- matrix(rnorm(Ntr * r), Ntr, r)
    U_tr_pre  <- F_pre %*% t(Lambda_tr) + matrix(rnorm(T0 * Ntr, sd = 0.3), T0, Ntr)
    grid      <- c(0.1, 1.0, 10.0)
    res <- fect:::.cv_gamma_loading(
        U_tr_pre   = U_tr_pre,
        F_hat_pre  = F_pre,
        Lambda_co  = Lambda_co,
        gamma_grid = grid,
        cv_k       = 3L
    )
    expect_true(res$gamma_cv %in% grid)
    expect_length(res$cv_mspe, length(grid))
    expect_true(all(is.finite(res$cv_mspe)))
})
