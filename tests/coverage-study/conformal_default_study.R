## -----------------------------------------------------------------------------
## conformal_default_study.R
##
## Picks the default for vartype = "conformal" on evidence. Sweeps the `scale`
## knob (and, separately, the `weight` knob on a staggered DGP) across five data-
## generating processes, recording, for a true effect of zero:
##   - scalar att.avg coverage (target 1 - alpha) and median CI width,
##   - status != "ok" count (must be 0; Family A never empties),
##   - pointwise and simultaneous JOINT (whole-post-path) coverage.
##
## Decision rule: the option that holds scalar coverage across ALL five DGPs at
## the smallest median width is the default. Run:
##   Rscript tests/coverage-study/conformal_default_study.R [reps] [alpha]
## Writes a markdown table to tests/coverage-study/results/.
## -----------------------------------------------------------------------------

suppressMessages(devtools::load_all("/Users/xyq/GitHub/fect", quiet = TRUE))

args  <- commandArgs(trailingOnly = TRUE)
reps  <- if (length(args) >= 1) as.integer(args[1]) else 150L
alpha <- if (length(args) >= 2) as.numeric(args[2]) else 0.10
N <- 25L; T <- 20L; T0 <- 15L; r <- 2L
post.terms <- 1:(T - T0)

## --- DGPs: return a TT x N outcome matrix with true effect 0 (controls only) --
## Shared latent structure regenerated per rep via the seed stream.
gen <- function(dgp, Fm, L) {
  mu <- Fm %*% t(L)
  E  <- matrix(rnorm(N * T), T, N)
  if (dgp == "iid") {
    NULL
  } else if (dgp == "ar1") {
    rho <- 0.5
    for (t in 2:T) E[t, ] <- rho * E[t - 1, ] + sqrt(1 - rho^2) * E[t, ]
  } else if (dgp == "hetero") {
    s <- rep(c(1, 3), length.out = N)             # half the units 3x noisier
    E <- E * matrix(s, T, N, byrow = TRUE)
  } else if (dgp == "nonstat") {
    ramp <- 1 + 1.2 * (seq_len(T) - 1) / (T - 1)  # variance grows over time
    g    <- runif(N, 0.5, 1.5)                     # unit-specific drift rate
    E <- E * outer(ramp, g)
  }
  mu + E
}

dgps   <- c("iid", "ar1", "hetero", "nonstat")
scales <- c("none", "sd", "rmspe", "mad", "diff")

joint <- function(m) {
  p <- m[as.numeric(rownames(m)) >= 0, c("CI.lower", "CI.upper"), drop = FALSE]
  all(p[, 1] <= 0 & 0 <= p[, 2])
}

fit_one <- function(Y, D, scale) {
  dat <- data.frame(id = rep(1:N, each = T), time = rep(1:T, N),
                    Y = as.vector(Y), D = as.vector(D))
  suppressWarnings(suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
        method = "gsynth", force = 3, CV = FALSE, r = r, se = TRUE,
        vartype = "conformal", conformal.scale = scale, alpha = alpha,
        parallel = FALSE)))
}

cat(sprintf("conformal default study: reps=%d alpha=%.2f (target cov %.2f)\n\n",
            reps, alpha, 1 - alpha))

## ===== Part 1: scale sweep, single treated unit (block) ======================
D1 <- matrix(0, T, N); D1[(T0 + 1):T, 1] <- 1
rows <- list()
t0 <- Sys.time()
for (dgp in dgps) {
  for (sc in scales) {
    set.seed(20260609)                       # same stream across scales -> paired
    cov <- wid <- jp <- js <- numeric(reps); bad <- 0L
    for (b in seq_len(reps)) {
      Fm <- matrix(rnorm(T * r), T, r); L <- matrix(rnorm(N * r), N, r)
      Y  <- gen(dgp, Fm, L)
      f  <- fit_one(Y, D1, sc)
      e  <- f$est.avg[1, ]
      if (f$conformal$status != "ok") bad <- bad + 1L
      cov[b] <- (e["CI.lower"] <= 0) && (0 <= e["CI.upper"])
      wid[b] <- e["CI.upper"] - e["CI.lower"]
      jp[b]  <- joint(f$est.att); js[b] <- joint(f$est.att.sim)
    }
    rows[[length(rows) + 1]] <- data.frame(
      dgp = dgp, scale = sc, cov = mean(cov), med.width = median(wid),
      bad = bad, joint.point = mean(jp), joint.sim = mean(js))
    cat(sprintf("  [%-8s %-6s] cov=%.3f width=%.2f bad=%d joint(pt/sim)=%.2f/%.2f\n",
                dgp, sc, mean(cov), median(wid), bad, mean(jp), mean(js)))
  }
}
tab <- do.call(rbind, rows)
cat(sprintf("\nPart 1 done in %.0fs\n", as.numeric(Sys.time() - t0, units = "secs")))

## worst-case scalar coverage per scale across DGPs, and mean width
summ <- do.call(rbind, lapply(scales, function(sc) {
  s <- tab[tab$scale == sc, ]
  data.frame(scale = sc, worst.cov = min(s$cov), mean.cov = mean(s$cov),
             mean.width = mean(s$med.width), max.bad = max(s$bad))
}))
summ <- summ[order(-summ$worst.cov, summ$mean.width), ]
cat("\n=== scale ranking (worst-case scalar coverage, then width) ===\n")
print(summ, row.names = FALSE, digits = 3)
best <- summ$scale[1]
cat(sprintf("\nProvisional default scale (most robust at smallest width): %s\n", best))

## ===== Part 2: weight comparison on a staggered DGP ==========================
set.seed(424242)
cat("\n=== Part 2: weight sweep, staggered + heteroskedastic treated ===\n")
Ds <- matrix(0, T, N); Ds[(T0 - 2 + 1):T, 1] <- 1; Ds[(T0 + 1):T, 2] <- 1; Ds[(T0 + 3):T, 3] <- 1
wrows <- list()
for (w in c("cell", "unit", "precision")) {
  cov <- wid <- numeric(reps); bad <- 0L
  set.seed(13)
  for (b in seq_len(reps)) {
    Fm <- matrix(rnorm(T * r), T, r); L <- matrix(rnorm(N * r), N, r)
    E <- matrix(rnorm(N * T), T, N); E[, 1] <- E[, 1] * 3      # treated 1 noisy
    Y <- Fm %*% t(L) + E
    dat <- data.frame(id = rep(1:N, each = T), time = rep(1:T, N),
                      Y = as.vector(Y), D = as.vector(Ds))
    f <- tryCatch(suppressWarnings(suppressMessages(fect(Y ~ D, data = dat,
              index = c("id", "time"), method = "gsynth", force = 3, CV = FALSE,
              r = r, se = TRUE, vartype = "conformal", conformal.weight = w,
              alpha = alpha, parallel = FALSE))), error = function(e) NULL)
    if (is.null(f)) { bad <- bad + 1L; cov[b] <- NA; wid[b] <- NA; next }
    e <- f$est.avg[1, ]
    cov[b] <- (e["CI.lower"] <= 0) && (0 <= e["CI.upper"]); wid[b] <- e["CI.upper"] - e["CI.lower"]
  }
  wrows[[length(wrows) + 1]] <- data.frame(weight = w, cov = mean(cov, na.rm = TRUE),
                                           med.width = median(wid, na.rm = TRUE), bad = bad)
  cat(sprintf("  [weight=%-9s] cov=%.3f width=%.2f bad=%d\n", w,
              mean(cov, na.rm = TRUE), median(wid, na.rm = TRUE), bad))
}
wtab <- do.call(rbind, wrows)

## ===== write results =========================================================
out <- "/Users/xyq/GitHub/fect/tests/coverage-study/results/conformal_default_study.md"
con <- file(out, "w")
writeLines(c(
  "# Conformal default study",
  sprintf("reps=%d, alpha=%.2f, target coverage=%.2f, N=%d T=%d T0=%d r=%d, true effect 0.",
          reps, alpha, 1 - alpha, N, T, T0, r),
  "", "## Part 1: scale sweep (single treated, block)", "",
  "| dgp | scale | cov | med.width | bad | joint.point | joint.sim |",
  "|---|---|---|---|---|---|---|"), con)
for (i in seq_len(nrow(tab))) writeLines(sprintf("| %s | %s | %.3f | %.2f | %d | %.2f | %.2f |",
  tab$dgp[i], tab$scale[i], tab$cov[i], tab$med.width[i], tab$bad[i],
  tab$joint.point[i], tab$joint.sim[i]), con)
writeLines(c("", "## Scale ranking (worst-case scalar coverage, then width)", "",
  "| scale | worst.cov | mean.cov | mean.width | max.bad |",
  "|---|---|---|---|---|"), con)
for (i in seq_len(nrow(summ))) writeLines(sprintf("| %s | %.3f | %.3f | %.2f | %d |",
  summ$scale[i], summ$worst.cov[i], summ$mean.cov[i], summ$mean.width[i], summ$max.bad[i]), con)
writeLines(c("", sprintf("**Provisional default scale: %s** (most robust at smallest width).", best),
  "", "## Part 2: weight sweep (staggered + heteroskedastic treated)", "",
  "| weight | cov | med.width | bad |", "|---|---|---|---|"), con)
for (i in seq_len(nrow(wtab))) writeLines(sprintf("| %s | %.3f | %.2f | %d |",
  wtab$weight[i], wtab$cov[i], wtab$med.width[i], wtab$bad[i]), con)
close(con)
cat(sprintf("\nResults written to %s\n", out))
