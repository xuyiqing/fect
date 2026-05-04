## Targeted coverage test: gsynth-style factor-model regime.
##
## Purpose: diagnose whether fect's parametric coverage approaches 0.95 when
## the DGP matches gsynth-note's `sim_coverage.R` *structure* (factor model
## r=2, T0/T = ~72%, Ntr/Nco = 1:3, fixed treated units, gsynth-style
## estimator with `time.component.from = "nevertreated"`).
##
## Quarter-scale of gsynth-note (N=40 vs 160) so K=500 finishes in ~2 hours
## at cores=10. ATT held constant (gsynth-note has heterogeneous ATT 1:10),
## no covariates (gsynth-note has X1, X2), IID Gaussian residuals.
##
## Single configuration (per user request):
##   - vartype = "parametric", para.error = "auto" (resolves to "empirical"
##     on this fully-observed panel)
##   - ci.method = "normal" (the att default; matches what fit$est.att uses)
##   - K = 500 outer reps, nboots = 200 inner

suppressPackageStartupMessages({
  library(devtools)
})
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
pkg_root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(pkg_root)
cat("Loading fect from:", pkg_root, "\n")
suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))

## ------------------------------------------------------------ DGP --------

## Factor-model DGP, gsynth-note structure (quarter-scale).
##   Ntr = 10 treated (units 1:10, fixed across replications)
##   Nco = 30 controls (units 11:40)
##   T   = 25 periods (T0 = 18 pre-treatment, 7 post-treatment)
##   r   = 2 latent factors
##   ATT = 3.0 constant
make_factor_data <- function(seed) {
    set.seed(seed)
    Ntr <- 10; Nco <- 30; N <- Ntr + Nco
    TT <- 25; T0 <- 18
    r <- 2

    F_mat <- matrix(rnorm(TT * r), TT, r)
    L_mat <- matrix(rnorm(N * r), N, r)
    alpha_i <- rnorm(N, 0, 1)
    xi_t    <- rnorm(TT, 0, 1)
    eps     <- matrix(rnorm(N * TT, 0, 1), TT, N)

    ## Treatment block: units 1:Ntr treated from period T0+1 onward
    D <- matrix(0L, TT, N); D[(T0 + 1):TT, 1:Ntr] <- 1L

    Y <- outer(xi_t, rep(1, N)) + outer(rep(1, TT), alpha_i) +
         F_mat %*% t(L_mat) + 3 * D + eps
    data.frame(id = rep(1:N, each = TT), time = rep(1:TT, N),
               Y  = c(Y), D = c(D))
}

## ------------------------------------------------------------ FIT --------

fect_gsynth_style <- function(df, nboots = 200, cores = 10) {
    suppressMessages(suppressWarnings(
        fect(Y ~ D, data = df, index = c("id", "time"),
             method = "ife", r = 2, force = "two-way",
             se = TRUE, vartype = "parametric", para.error = "auto",
             nboots = nboots, parallel = TRUE, cores = cores,
             time.component.from = "nevertreated",
             keep.sims = TRUE, CV = FALSE)
    ))
}

## ------------------------------------------------------------ COVERAGE ---

run_coverage <- function(K = 500, nboots = 200, cores = 10) {
    cat(sprintf("\n=== gsynth-style coverage test ===\n"))
    cat(sprintf("DGP: N=40 (Ntr=10, Nco=30), TT=25, T0=18, r=2 latent factors\n"))
    cat(sprintf("Estimator: method='ife' r=2 force='two-way' time.component.from='nevertreated'\n"))
    cat(sprintf("Inference: vartype='parametric' para.error='auto' ci.method='normal'\n"))
    cat(sprintf("Sims: K=%d outer reps, nboots=%d inner, cores=%d\n", K, nboots, cores))

    cover <- numeric(K)
    width <- numeric(K)
    estim <- numeric(K)

    t0 <- Sys.time()
    for (r in seq_len(K)) {
        df <- make_factor_data(seed = 1000 + r)
        fit <- fect_gsynth_style(df, nboots = nboots, cores = cores)
        res <- suppressWarnings(estimand(fit, "att", "overall",
                                         window = c(1, 7),
                                         ci.method = "normal"))
        estim[r] <- res$estimate
        width[r] <- res$ci.hi - res$ci.lo
        cover[r] <- (res$ci.lo <= 3.0) && (3.0 <= res$ci.hi)

        if (r %% 10 == 0 || r == K) {
            elapsed <- as.numeric(Sys.time() - t0, units = "mins")
            eta <- elapsed / r * (K - r)
            cat(sprintf("  rep %4d / %d  elapsed %.1f min  eta %.1f min  running coverage %.3f\n",
                        r, K, elapsed, eta, mean(cover[1:r])))
        }
    }

    coverage <- mean(cover)
    se <- sqrt(coverage * (1 - coverage) / K)
    cat(sprintf("\n=== RESULT ===\n"))
    cat(sprintf("Coverage: %.3f  (MC SE: %.3f, 95%% CI: [%.3f, %.3f])\n",
                coverage, se, coverage - 1.96 * se, coverage + 1.96 * se))
    cat(sprintf("Mean CI width: %.4f\n", mean(width)))
    cat(sprintf("Mean estimate: %.4f  (true ATT = 3.0; bias %.4f)\n",
                mean(estim), mean(estim) - 3.0))

    out <- data.frame(rep = seq_len(K), estimate = estim,
                      width = width, cover = cover)
    csv_path <- sprintf("/tmp/fect-coverage-study/gsynth_style_K%d_nb%d_%s.csv",
                        K, nboots, format(Sys.time(), "%Y%m%d-%H%M%S"))
    dir.create(dirname(csv_path), showWarnings = FALSE, recursive = TRUE)
    write.csv(out, csv_path, row.names = FALSE)
    cat(sprintf("Per-rep raw saved to: %s\n", csv_path))

    invisible(list(coverage = coverage, se = se,
                   width = mean(width), bias = mean(estim) - 3.0))
}

## ------------------------------------------------------------ RUN --------

run_coverage(K = 500, nboots = 200, cores = 10)
