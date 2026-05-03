## CFE coverage at K=200 to resolve whether 0.917 was real bias or MC noise.
## If true coverage ~ 0.95, K=200 will give SE-of-coverage ~ 0.015 ≈ 1.5pp.

suppressPackageStartupMessages({ library(fect) })
make_factor_data <- function(N = 200, T = 35, Ntr = 100, tau = 3.0, r = 2,
                             T_treat_start = 18, sigma = 1.0, seed = 1) {
    set.seed(seed)
    F_mat <- matrix(rnorm(T * r), T, r)
    L_mat <- matrix(rnorm(N * r), N, r)
    Y_factor <- F_mat %*% t(L_mat)
    Y_noise <- matrix(rnorm(T * N, sd = sigma), T, N)
    D_mat <- matrix(0, T, N)
    D_mat[T_treat_start:T, 1:Ntr] <- 1
    Y_mat <- Y_factor + Y_noise + tau * D_mat
    data.frame(id = rep(1:N, each = T), time = rep(1:T, N),
               Y = as.vector(Y_mat), D = as.vector(D_mat))
}

K <- 200
cat(sprintf("=== CFE coverage K=%d at new defaults (tol=1e-5) ===\n", K))
cover <- numeric(K); att <- numeric(K); se <- numeric(K)
t0 <- Sys.time()
for (k in 1:K) {
    s <- 5000 + k
    df <- make_factor_data(seed = s)
    fit <- fect(Y ~ D, data = df, index = c("id","time"),
                method = "cfe", r = 2, force = "two-way",
                se = TRUE, vartype = "bootstrap", nboots = 100,
                parallel = FALSE, CV = FALSE, seed = s)
    att[k] <- fit$est.avg[1, "ATT.avg"]
    se[k] <- fit$est.avg[1, "S.E."]
    cover[k] <- (fit$est.avg[1, "CI.lower"] < 3) && (fit$est.avg[1, "CI.upper"] > 3)
    if (k %% 50 == 0) {
        cat(sprintf("  k=%d: cumulative coverage = %.4f   bias = %+.4f   elapsed=%.0fs\n",
                    k, mean(cover[1:k]), mean(att[1:k] - 3.0),
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
}
cat(sprintf("\nFinal: coverage = %.4f   bias = %+.4f   emp.SE = %.4f   boot.SE = %.4f\n",
            mean(cover), mean(att - 3.0), sd(att), mean(se)))
cat(sprintf("Coverage 95%% CI (Wilson): [%.4f, %.4f]\n",
            (function(p, n) {
                z <- 1.96; den <- 1 + z^2/n
                ctr <- (p + z^2/(2*n)) / den
                hw <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / den
                c(ctr - hw, ctr + hw)
            })(mean(cover), K)))
saveRDS(list(cover = cover, att = att, se = se),
        "/tmp/fect-tol-char/cfe_high_K_K200.rds")
