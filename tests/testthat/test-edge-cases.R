## ---------------------------------------------------------
## Edge case and boundary condition tests
## ---------------------------------------------------------

test_that("fect handles single treated unit", {

  skip_on_cran()

    set.seed(7001)

    N <- 20
    TT <- 15
    T0 <- 10
    Ntr <- 1    # only one treated unit
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
            eps <- rnorm(1, 0, 0.5)
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

    ## Should not error with a single treated unit
    out <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "fe",
        force = "two-way",
        se = FALSE,
        parallel = FALSE
    ))

    expect_true(is.numeric(out$att.avg))
    expect_true(!is.na(out$att.avg))
    expect_true(is.matrix(out$eff))
})

test_that("fect handles missing data (scattered NAs in Y)", {

  skip_on_cran()

    set.seed(7002)

    N <- 25
    TT <- 15
    T0 <- 10
    Ntr <- 8
    tau <- 2.0

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

    ## Introduce ~5% missing values in Y (but NOT in treated post-period)
    miss_candidates <- which(D_vec == 0)
    n_miss <- min(floor(0.05 * length(Y_vec)), length(miss_candidates))
    miss_idx <- sample(miss_candidates, n_miss)
    Y_vec[miss_idx] <- NA

    simdf <- data.frame(
        id = id_vec,
        time = time_vec,
        Y = Y_vec,
        D = D_vec
    )

    ## Should handle missing data without error
    out <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "fe",
        force = "two-way",
        se = FALSE,
        parallel = FALSE,
        na.rm = TRUE
    ))

    expect_true(is.numeric(out$att.avg))
    expect_true(!is.na(out$att.avg))
    ## ATT should still be roughly correct
    expect_lt(abs(out$att.avg - tau), 1.5,
        label = paste0("Missing data FE: att.avg = ",
                       round(out$att.avg, 4), ", tau = ", tau))
})

test_that("fect output structure is complete (FE)", {

  skip_on_cran()

    set.seed(7003)

    N <- 20
    TT <- 10
    T0 <- 6
    Ntr <- 8

    Y_vec <- numeric(N * TT)
    D_vec <- integer(N * TT)
    id_vec <- integer(N * TT)
    time_vec <- integer(N * TT)

    idx <- 1
    for (i in 1:N) {
        for (t in 1:TT) {
            treated <- (i <= Ntr) && (t > T0)
            D_vec[idx] <- as.integer(treated)
            Y_vec[idx] <- rnorm(1) + 2.0 * D_vec[idx]
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

    out <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "fe",
        force = "two-way",
        se = FALSE,
        parallel = FALSE
    ))

    ## Check required output slots exist
    expect_true(!is.null(out$att.avg))
    expect_true(!is.null(out$att.avg.unit))
    expect_true(!is.null(out$eff))
    expect_true(!is.null(out$Y.ct))
    expect_true(!is.null(out$time))
    expect_true(!is.null(out$att))
    expect_true(!is.null(out$count))
    expect_true(!is.null(out$sigma2))
})

test_that("fect with many treated units still works", {

  skip_on_cran()

    set.seed(7004)

    N <- 25
    TT <- 12
    T0 <- 8
    Ntr <- 20  # 80% treated
    tau <- 1.5

    alpha_i <- rnorm(N, 0, 0.5)
    xi_t <- rnorm(TT, 0, 0.3)

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

    ## Should work even with many treated / few control
    out <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "fe",
        force = "two-way",
        se = FALSE,
        parallel = FALSE
    ))

    expect_true(is.numeric(out$att.avg))
    expect_true(!is.na(out$att.avg))
})
