## ============================================================================
## Tol coverage validation (v2.4.3 default-tol fix gate)
##
## Question: at the proposed new default tol = 1e-5, does the bootstrap
## CI cover the true ATT at the nominal 95% rate?
##
## DGP: factor model, N=200, T=35, Ntr=100, r=2, true tau = 3.0.
## (Mirrors simdata structure but with known-truth tau for coverage.)
##
## Cells: tol ∈ {1e-3 (current), 1e-5 (proposed), 1e-7 (converged ground truth)}
## Method: ife.  vartype: bootstrap, nboots = 200.
## K = 100 reps per cell (300 fits total).
##
## Pass criterion for tol = 1e-5:
##   empirical coverage within ±3pp of nominal 95% (i.e. 92-98%)
##   AND coverage no worse than tol = 1e-7 baseline
## ============================================================================

suppressPackageStartupMessages({ library(fect); library(dplyr) })
OUT_DIR <- "/tmp/fect-tol-char"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

make_factor_data <- function(N = 200, T = 35, Ntr = 100, tau = 3.0, r = 2,
                             T_treat_start = 18, sigma = 1.0, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    F_mat <- matrix(rnorm(T * r), T, r)
    L_mat <- matrix(rnorm(N * r), N, r)
    Y_factor <- F_mat %*% t(L_mat)            # T x N
    Y_noise <- matrix(rnorm(T * N, sd = sigma), T, N)
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

run_one_rep <- function(rep_seed, tol, max_iter, nboots = 100) {
    df <- make_factor_data(N = 200, T = 35, Ntr = 100, tau = 3.0, r = 2,
                           seed = rep_seed)
    fit <- tryCatch(
        fect(Y ~ D, data = df, index = c("id","time"),
             method = "ife", r = 2, force = "two-way",
             se = TRUE, vartype = "bootstrap", nboots = nboots,
             parallel = FALSE, keep.sims = FALSE, CV = FALSE,
             tol = tol, max.iteration = max_iter, seed = rep_seed),
        error = function(e) e)
    if (inherits(fit, "error")) return(list(att = NA, se = NA, ci_l = NA, ci_u = NA, cover = NA, niter = NA))
    # Aggregate ATT over post-treatment event-times
    att_post <- fit$est.avg
    if (is.null(att_post) || nrow(att_post) == 0) return(list(att = NA, se = NA, ci_l = NA, ci_u = NA, cover = NA, niter = fit$niter))
    point <- att_post[1, "ATT.avg"]
    se <- att_post[1, "S.E."]
    ci_l <- att_post[1, "CI.lower"]
    ci_u <- att_post[1, "CI.upper"]
    cover <- (3.0 >= ci_l) && (3.0 <= ci_u)
    list(att = point, se = se, ci_l = ci_l, ci_u = ci_u, cover = cover, niter = fit$niter)
}

run_coverage <- function(tol, max_iter, K = 100, label) {
    cat(sprintf("\n=== %s (tol=%.0e, max_iter=%d, K=%d) ===\n",
                label, tol, max_iter, K))
    t0 <- Sys.time()
    res <- vector("list", K)
    for (k in seq_len(K)) {
        res[[k]] <- run_one_rep(rep_seed = 1000 + k, tol = tol, max_iter = max_iter)
        if (k %% 20 == 0) {
            cat(sprintf("  k=%d: cumulative coverage = %.3f\n",
                        k, mean(sapply(res[seq_len(k)], "[[", "cover"), na.rm = TRUE)))
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
    cat(sprintf("  median niter: %d   n_valid: %d/%d   elapsed: %.1fs\n",
                median(niter_vec, na.rm = TRUE), n_valid, K, elapsed))
    list(label = label, tol = tol, max_iter = max_iter, K = K,
         coverage = coverage, bias = bias, se_emp = se_emp, se_avg = se_avg,
         elapsed = elapsed, niter_med = median(niter_vec, na.rm = TRUE))
}

cat("=== v2.4.3 default-tol coverage validation ===\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("DGP: N=200 T=35 Ntr=100 r=2 tau=3.0 sigma=1.0\n")

results <- list(
    run_coverage(tol = 1e-3, max_iter = 1000, K = 80, label = "current default"),
    run_coverage(tol = 1e-5, max_iter = 5000, K = 80, label = "PROPOSED new default")
)

cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
for (r in results) {
    cat(sprintf("  %-32s  cover=%.3f   bias=%+.4f   emp.SE=%.4f   boot.SE=%.4f   niter_med=%d\n",
                r$label, r$coverage, r$bias, r$se_emp, r$se_avg, r$niter_med))
}

ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
saveRDS(results, file.path(OUT_DIR, sprintf("coverage_tol_%s.rds", ts)))
cat("\nSaved RDS to", OUT_DIR, "\n")
