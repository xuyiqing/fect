## ============================================================================
## Extended tol coverage: CFE + AR(1) noise + multi-seed stability
##
## After Phase 1c showed inference is preserved at tol=1e-3 on a clean IFE
## DGP, test the harder cases:
##   1. CFE r=2 (was 40% gap in Phase 1a)
##   2. AR(1) noise (real-world serial correlation)
##   3. Multi-seed att.avg stability (does the deterministic output drift
##      with initialization perturbations?)
## ============================================================================

suppressPackageStartupMessages({ library(fect); library(dplyr) })
OUT_DIR <- "/tmp/fect-tol-char"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

make_factor_data <- function(N = 200, T = 35, Ntr = 100, tau = 3.0, r = 2,
                             T_treat_start = 18, sigma = 1.0, ar1_rho = 0,
                             seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    F_mat <- matrix(rnorm(T * r), T, r)
    L_mat <- matrix(rnorm(N * r), N, r)
    Y_factor <- F_mat %*% t(L_mat)
    if (ar1_rho == 0) {
        Y_noise <- matrix(rnorm(T * N, sd = sigma), T, N)
    } else {
        Y_noise <- matrix(0, T, N)
        for (i in 1:N) {
            eps <- rnorm(T, sd = sigma * sqrt(1 - ar1_rho^2))
            y <- numeric(T); y[1] <- rnorm(1, sd = sigma)
            for (t in 2:T) y[t] <- ar1_rho * y[t-1] + eps[t]
            Y_noise[, i] <- y
        }
    }
    D_mat <- matrix(0, T, N)
    treat_units <- 1:Ntr
    D_mat[T_treat_start:T, treat_units] <- 1
    Y_mat <- Y_factor + Y_noise + tau * D_mat
    df <- data.frame(
        id = rep(1:N, each = T),
        time = rep(1:T, N),
        Y = as.vector(Y_mat),
        D = as.vector(D_mat)
    )
    df
}

run_one_rep <- function(rep_seed, tol, max_iter, method = "ife", ar1_rho = 0,
                        nboots = 80) {
    df <- make_factor_data(N = 200, T = 35, Ntr = 100, tau = 3.0, r = 2,
                           ar1_rho = ar1_rho, seed = rep_seed)
    fit <- tryCatch(
        fect(Y ~ D, data = df, index = c("id","time"),
             method = method, r = 2, force = "two-way",
             se = TRUE, vartype = "bootstrap", nboots = nboots,
             parallel = FALSE, keep.sims = FALSE, CV = FALSE,
             tol = tol, max.iteration = max_iter, seed = rep_seed),
        error = function(e) e)
    if (inherits(fit, "error")) return(list(att = NA, se = NA, ci_l = NA, ci_u = NA, cover = NA, niter = NA))
    att_post <- fit$est.avg
    if (is.null(att_post) || nrow(att_post) == 0) return(list(att = NA, se = NA, ci_l = NA, ci_u = NA, cover = NA, niter = fit$niter))
    point <- att_post[1, "ATT.avg"]
    se <- att_post[1, "S.E."]
    ci_l <- att_post[1, "CI.lower"]
    ci_u <- att_post[1, "CI.upper"]
    cover <- (3.0 >= ci_l) && (3.0 <= ci_u)
    list(att = point, se = se, ci_l = ci_l, ci_u = ci_u, cover = cover, niter = as.integer(fit$niter))
}

run_coverage <- function(tol, max_iter, method = "ife", ar1_rho = 0, K = 60, label) {
    cat(sprintf("\n=== %s (method=%s, tol=%.0e, max_iter=%d, ar1=%.1f, K=%d) ===\n",
                label, method, tol, max_iter, ar1_rho, K))
    t0 <- Sys.time()
    res <- vector("list", K)
    for (k in seq_len(K)) {
        res[[k]] <- run_one_rep(rep_seed = 5000 + k, tol = tol, max_iter = max_iter,
                                method = method, ar1_rho = ar1_rho)
        if (k %% 20 == 0) {
            cov_so_far <- mean(sapply(res[seq_len(k)], "[[", "cover"), na.rm = TRUE)
            cat(sprintf("  k=%d: cumulative coverage = %.3f\n", k, cov_so_far))
        }
    }
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cover_vec <- sapply(res, "[[", "cover")
    att_vec <- sapply(res, "[[", "att")
    se_vec <- sapply(res, "[[", "se")
    niter_vec <- sapply(res, "[[", "niter")
    n_valid <- sum(!is.na(cover_vec))
    coverage <- mean(cover_vec, na.rm = TRUE)
    se_emp <- sd(att_vec, na.rm = TRUE)
    se_avg <- mean(se_vec, na.rm = TRUE)
    bias <- mean(att_vec - 3.0, na.rm = TRUE)
    cat(sprintf("\n  Coverage: %.3f   bias: %+.4f   emp.SE: %.4f   avg.boot.SE: %.4f\n",
                coverage, bias, se_emp, se_avg))
    cat(sprintf("  median niter: %s   n_valid: %d/%d   elapsed: %.1fs\n",
                as.character(median(niter_vec, na.rm = TRUE)),
                n_valid, K, elapsed))
    list(label = label, method = method, tol = tol, max_iter = max_iter, ar1 = ar1_rho,
         K = K, coverage = coverage, bias = bias, se_emp = se_emp, se_avg = se_avg,
         elapsed = elapsed,
         niter_med = median(niter_vec, na.rm = TRUE))
}

cat("=== Extended coverage validation (v2.4.3) ===\n")
cat("Started:", as.character(Sys.time()), "\n")

results <- list(
    ## CFE on IID DGP (Phase 1a worst case at tol=1e-3)
    run_coverage(tol = 1e-3, max_iter = 1000, method = "cfe", ar1_rho = 0, K = 60,
                 label = "CFE iid old default"),
    run_coverage(tol = 1e-5, max_iter = 5000, method = "cfe", ar1_rho = 0, K = 60,
                 label = "CFE iid NEW default"),

    ## IFE on AR(1) noise (real-world serial correlation)
    run_coverage(tol = 1e-3, max_iter = 1000, method = "ife", ar1_rho = 0.5, K = 60,
                 label = "IFE AR(0.5) old default"),
    run_coverage(tol = 1e-5, max_iter = 5000, method = "ife", ar1_rho = 0.5, K = 60,
                 label = "IFE AR(0.5) NEW default")
)

cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
for (r in results) {
    cat(sprintf("  %-32s  cover=%.3f   bias=%+.4f   emp.SE=%.4f   boot.SE=%.4f   niter_med=%s\n",
                r$label, r$coverage, r$bias, r$se_emp, r$se_avg,
                as.character(r$niter_med)))
}

ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
saveRDS(results, file.path(OUT_DIR, sprintf("coverage_extended_%s.rds", ts)))
cat("\nSaved to:", OUT_DIR, "\n")
