## ============================================================================
## MC partial warm-start validation (v2.4.3)
##
## Per Yiqing's request 2026-05-02: "If MC fails, you can also try partial
## warm-start with MC."
##
## MC differs from IFE in two relevant ways:
##   1. Optimization is convex (nuclear-norm regularized) -- no basin
##      identification issue at all. Warm-start should be variance-neutral
##      in theory regardless of tol.
##   2. EM converges differently due to soft-thresholding step.
##
## Test: cold vs partial-warm bootstrap on simdata + MC.
## ============================================================================

suppressPackageStartupMessages({ library(devtools) })
setwd("/Users/xyq/GitHub/fect-warmstart")
devtools::load_all(quiet = TRUE)

cat("\n=== MC partial warm-start (n=100, tol=1e-5) ===\n")

data(simdata)

common <- list(
    formula = Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
    method = "mc", lambda = 0.05, force = "two-way",
    se = TRUE, vartype = "bootstrap", nboots = 100,
    parallel = FALSE, keep.sims = TRUE, CV = FALSE,
    tol = 1e-5, max.iteration = 5000
)

set.seed(42)
t0 <- Sys.time()
fit_cold <- tryCatch(do.call(fect, common), error = function(e) e)
ct <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (inherits(fit_cold, "error")) {
    cat("cold FAILED:", conditionMessage(fit_cold), "\n")
    quit(save = "no", status = 1)
}
cat(sprintf("MC cold time: %.2fs   att.avg: %.4f\n", ct, fit_cold$att.avg))

set.seed(42)
t0 <- Sys.time()
fit_warm <- tryCatch(do.call(fect, c(common, list(warm.start = "linear"))),
                     error = function(e) e)
wt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (inherits(fit_warm, "error")) {
    cat("warm FAILED:", conditionMessage(fit_warm), "\n")
    quit(save = "no", status = 1)
}
cat(sprintf("MC warm time: %.2fs   att.avg: %.4f\n", wt, fit_warm$att.avg))

cat(sprintf("\nMC speedup:        %.2fx %s\n",
            ct/wt, if (ct/wt >= 2) "PASS" else "FAIL"))
cat(sprintf("MC att.avg diff:   %.2e %s\n",
            abs(fit_cold$att.avg - fit_warm$att.avg),
            if (abs(fit_cold$att.avg - fit_warm$att.avg) < 1e-8) "PASS" else "FAIL"))
cat(sprintf("MC eff.boot diff:  %.4e\n",
            max(abs(fit_cold$eff.boot - fit_warm$eff.boot), na.rm=TRUE)))
se_diff <- max(abs(fit_cold$est.att[, "S.E."] - fit_warm$est.att[, "S.E."]) /
               pmax(fit_cold$est.att[, "S.E."], 1e-10))
cat(sprintf("MC rel S.E. diff:  %.2f%% %s\n",
            100 * se_diff, if (se_diff < 0.05) "PASS" else "FAIL"))
