## ============================================================================
## Tol-convergence characterization: where does fect's EM actually converge?
## ============================================================================
## v2.4.3 investigation (2026-05-02): default tol = 1e-3 produces
## under-converged IFE/CFE estimates. Question: what's the right new default?
##
## Approach: tol-sweep across methods x DGPs x panel sizes. For each cell,
## measure att.avg, niter, wall time. Derive: at what tol does att.avg
## stabilize to within X% of the truly-converged value?
##
## Pass criterion for a candidate default tol:
##   max relative att.avg gap (vs tol=1e-7) < 1% across all tested cells
##   AND wall time penalty vs current default <= 5x
##
## RUN: ~10-15 min wall time. Output:
##   /tmp/fect-tol-char/tol_sweep_<ts>.csv
## ============================================================================

suppressPackageStartupMessages({ library(fect); library(dplyr) })
OUT_DIR <- "/tmp/fect-tol-char"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

run_tol_cell <- function(method, formula, data, index, r = 0, lambda = NULL,
                         force = "two-way", time.component.from = "notyettreated",
                         tol_seq = c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7),
                         max_iter_seq = c(1000, 1000, 5000, 10000, 50000),
                         label = "") {
    cells <- list()
    for (i in seq_along(tol_seq)) {
        t <- tol_seq[i]
        mi <- max_iter_seq[i]
        t0 <- Sys.time()
        args <- list(formula = formula, data = data, index = index,
                     method = method, force = force,
                     time.component.from = time.component.from,
                     se = FALSE, CV = FALSE,
                     tol = t, max.iteration = mi)
        if (!is.null(r)) args$r <- r
        if (!is.null(lambda)) args$lambda <- lambda
        fit <- tryCatch(do.call(fect, args), error = function(e) e)
        elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        if (inherits(fit, "error")) {
            cells[[i]] <- data.frame(
                label = label, method = method, tol = t, max_iter = mi,
                niter = NA, att.avg = NA, elapsed = elapsed,
                error = conditionMessage(fit), stringsAsFactors = FALSE)
        } else {
            cells[[i]] <- data.frame(
                label = label, method = method, tol = t, max_iter = mi,
                niter = fit$niter %||% NA,
                att.avg = fit$att.avg, elapsed = elapsed,
                error = NA_character_, stringsAsFactors = FALSE)
        }
    }
    do.call(rbind, cells)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

cat("=== Tol-convergence characterization ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

results <- list()

cat("[1/8] simdata + IFE r=2 (canonical IFE benchmark)\n")
data(simdata)
results[[1]] <- run_tol_cell(
    method = "ife", formula = Y ~ D + X1 + X2, data = simdata,
    index = c("id","time"), r = 2,
    label = "simdata-ife-r2"
)
print(results[[1]][, c("tol","niter","att.avg","elapsed")])

cat("\n[2/8] simdata + CFE r=2 (canonical CFE benchmark)\n")
results[[2]] <- run_tol_cell(
    method = "cfe", formula = Y ~ D + X1 + X2, data = simdata,
    index = c("id","time"), r = 2,
    label = "simdata-cfe-r2"
)
print(results[[2]][, c("tol","niter","att.avg","elapsed")])

cat("\n[3/8] simdata + IFE r=1 (smaller rank)\n")
results[[3]] <- run_tol_cell(
    method = "ife", formula = Y ~ D + X1 + X2, data = simdata,
    index = c("id","time"), r = 1,
    label = "simdata-ife-r1"
)
print(results[[3]][, c("tol","niter","att.avg","elapsed")])

cat("\n[4/8] sim_gsynth + GSC r=2 (gsynth canonical)\n")
data(sim_gsynth)
results[[4]] <- run_tol_cell(
    method = "gsynth", formula = Y ~ D + X1 + X2, data = sim_gsynth,
    index = c("id","time"), r = 2,
    label = "sim_gsynth-gsynth-r2"
)
print(results[[4]][, c("tol","niter","att.avg","elapsed")])

cat("\n[5/8] simdata + MC lambda=0.05 (typical regularization)\n")
results[[5]] <- run_tol_cell(
    method = "mc", formula = Y ~ D + X1 + X2, data = simdata,
    index = c("id","time"), r = NULL, lambda = 0.05,
    label = "simdata-mc-lam05"
)
print(results[[5]][, c("tol","niter","att.avg","elapsed")])

cat("\n[6/8] simdata + MC lambda=0.01 (light regularization)\n")
results[[6]] <- run_tol_cell(
    method = "mc", formula = Y ~ D + X1 + X2, data = simdata,
    index = c("id","time"), r = NULL, lambda = 0.01,
    label = "simdata-mc-lam01"
)
print(results[[6]][, c("tol","niter","att.avg","elapsed")])

cat("\n[7/8] turnout + IFE r=2 (real-world panel, Liu Wang Xu 2024)\n")
data(turnout)
results[[7]] <- run_tol_cell(
    method = "ife", formula = turnout ~ policy_edr,
    data = turnout, index = c("abb","year"), r = 2,
    label = "turnout-ife-r2"
)
print(results[[7]][, c("tol","niter","att.avg","elapsed")])

cat("\n[8/8] hh2019 + IFE r=2 (real-world, smaller)\n")
data(hh2019)
results[[8]] <- run_tol_cell(
    method = "ife", formula = hr ~ indirect,
    data = hh2019, index = c("bfs","year"), r = 2,
    label = "hh2019-ife-r2"
)
print(results[[8]][, c("tol","niter","att.avg","elapsed")])

all_results <- do.call(rbind, results)
ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_path <- file.path(OUT_DIR, sprintf("tol_sweep_%s.csv", ts))
write.csv(all_results, out_path, row.names = FALSE)

cat("\n=== Summary: gap from converged (tol=1e-7) per cell ===\n")
summary_df <- all_results %>%
    group_by(label) %>%
    mutate(att_converged = att.avg[tol == min(tol)],
           gap_pct = 100 * (att.avg - att_converged) / att_converged,
           speedup_vs_1e7 = elapsed[tol == min(tol)] / elapsed) %>%
    select(label, tol, niter, att.avg, gap_pct, elapsed, speedup_vs_1e7)
print(as.data.frame(summary_df), digits = 4)

cat("\nSaved:", out_path, "\n")
cat("\n=== Verdict per candidate default ===\n")
for (cand in c(1e-3, 1e-4, 1e-5, 1e-6)) {
    sub <- subset(summary_df, tol == cand)
    if (nrow(sub) == 0) next
    max_gap <- max(abs(sub$gap_pct), na.rm = TRUE)
    median_speedup <- median(sub$speedup_vs_1e7, na.rm = TRUE)
    cat(sprintf("  tol = %.0e:  max |gap| = %.2f%%   median speedup vs 1e-7 = %.2fx\n",
                cand, max_gap, median_speedup))
}
