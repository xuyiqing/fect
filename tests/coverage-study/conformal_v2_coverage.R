## Coverage of vartype = "conformal" (pooled information-matched Jackknife+,
## Liu write-up v2) for method = "gsynth".
##
## One treated unit, factor DGP with r = 2 factors, T = 30 (T0 = 20, 10
## post-treatment dates), J = 30 or 60 never-treated controls, iid N(0, 1)
## errors or stationary AR(1) errors (rho = 0.5, unit variance). The treated
## unit's loadings, fixed effects and errors come from the controls'
## distribution, so the units are exchangeable. r is fixed at the truth
## (CV = FALSE). Reported: date-specific coverage of the 90% interval
## (gamma = alpha = 0.1) averaged over the 10 post-treatment dates, its range
## over dates, the finite-sample bound 1 - (2q - 1)/N, the mean width, and
## the coverage of the interval for the average effect over the 10 dates
## (an extension beyond v2). A second table covers five treated units (block
## and staggered adoption; the event-time and average-ATT intervals, also an
## extension beyond v2).
##
## Run from the package root (not part of testthat):
##   FECT_LIB=<library with this branch installed> \
##     Rscript tests/coverage-study/conformal_v2_coverage.R
## Uses 4 cores (forked) and writes results/conformal_v2_coverage.md next to
## this script.

lib <- Sys.getenv("FECT_LIB", "")
if (nzchar(lib)) {
  library(fect, lib.loc = lib)
} else {
  library(fect)
}
stopifnot("conformal_jplus" %in% ls(asNamespace("fect")))

NREP <- as.integer(Sys.getenv("NREP", "200"))
CORES <- 4L
TT <- 30L
T0 <- 20L
R <- 2L
GAMMA <- 0.1

sim_one <- function(J, err, seed) {
  set.seed(seed)
  N <- J + 1L
  F <- matrix(rnorm(TT * R), TT, R)
  L <- matrix(rnorm(N * R), N, R)
  alpha <- rnorm(N)
  xi <- rnorm(TT)
  if (err == "iid") {
    e <- matrix(rnorm(TT * N), TT, N)
  } else {
    rho <- 0.5
    e <- matrix(NA_real_, TT, N)
    e[1, ] <- rnorm(N)
    for (t in 2:TT) e[t, ] <- rho * e[t - 1, ] + sqrt(1 - rho^2) * rnorm(N)
  }
  Y0 <- outer(xi, rep(1, N)) + outer(rep(1, TT), alpha) + F %*% t(L) + e
  tau <- seq(1, 3, length.out = TT - T0)
  D <- matrix(0, TT, N)
  D[(T0 + 1):TT, 1] <- 1
  Y <- Y0
  Y[(T0 + 1):TT, 1] <- Y0[(T0 + 1):TT, 1] + tau
  df <- data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
                   Y = c(Y), D = c(D))
  fit <- suppressMessages(fect(Y ~ D, data = df, index = c("id", "time"),
                               method = "gsynth", r = R, CV = FALSE,
                               se = TRUE, vartype = "conformal",
                               alpha = GAMMA, parallel = FALSE))
  u <- fit$conformal$unit
  cover <- (u$CI.lower <= tau) & (tau <= u$CI.upper)
  avg <- fit$est.avg[1, ]
  data.frame(
    date = seq_along(tau), cover = cover, width = u$CI.upper - u$CI.lower,
    N = u$N, q = u$q,
    cover.avg = (avg["CI.lower"] <= mean(tau)) & (mean(tau) <= avg["CI.upper"]),
    width.avg = avg["CI.upper"] - avg["CI.lower"],
    row.names = NULL
  )
}

cells <- expand.grid(J = c(30L, 60L), err = c("iid", "ar1"),
                     stringsAsFactors = FALSE)
t.start <- Sys.time()
res <- list()
for (k in seq_len(nrow(cells))) {
  J <- cells$J[k]
  err <- cells$err[k]
  seeds <- 1e5 * k + seq_len(NREP)
  out <- parallel::mclapply(seq_along(seeds), function(b) {
    tryCatch(cbind(rep = b, sim_one(J, err, seeds[b])),
             error = function(e) NULL)
  }, mc.cores = CORES)
  n.fail <- sum(vapply(out, is.null, logical(1)))
  d <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  n.ok <- NREP - n.fail
  by.date <- tapply(d$cover, d$date, mean)
  N <- d$N[1]
  q <- d$q[1]
  avg.rows <- d[d$date == 1, ]
  res[[k]] <- data.frame(
    J = J, errors = err, reps = n.ok, failed = n.fail, N = N, q = q,
    bound = 1 - (2 * q - 1) / N,
    coverage = mean(d$cover),
    min.date = min(by.date), max.date = max(by.date),
    mc.se = stats::sd(tapply(d$cover, d$rep, mean)) / sqrt(n.ok),
    width = mean(d$width),
    coverage.avg = mean(avg.rows$cover.avg),
    width.avg = mean(avg.rows$width.avg)
  )
  cat(sprintf("J = %d, %s: coverage %.3f (%d reps, %.1f min so far)\n",
              J, err, mean(d$cover), n.ok,
              as.numeric(difftime(Sys.time(), t.start, units = "mins"))))
}
tab <- do.call(rbind, res)

## ---- extension beyond v2: several treated units -------------------------
## Five treated units, J = 40, iid errors, r = 2 fixed. "block": all adopt at
## T0 + 1 = 21; "staggered": three adopt at 17 and two at 21 (T = 30). Each
## unit's effect at event time l is tau_l, so the ATT at l is tau_l.
## Reported: coverage of the event-time intervals (est.att, averaged over the
## post-treatment event times with an interval) and of est.avg, their widths,
## and the mean width of the per-unit, per-date intervals.
sim_multi <- function(design, seed, J = 40L, Ntr = 5L) {
  set.seed(seed)
  N <- J + Ntr
  F <- matrix(rnorm(TT * R), TT, R)
  L <- matrix(rnorm(N * R), N, R)
  Y0 <- outer(rnorm(TT), rep(1, N)) + outer(rep(1, TT), rnorm(N)) +
    F %*% t(L) + matrix(rnorm(TT * N), TT, N)
  adopt <- if (design == "block") rep(T0 + 1L, Ntr) else c(17L, 17L, 17L, 21L, 21L)
  D <- matrix(0, TT, N)
  for (i in seq_len(Ntr)) D[seq_len(TT) >= adopt[i], i] <- 1
  ev <- matrix(0, TT, N)
  for (i in seq_len(Ntr)) ev[, i] <- seq_len(TT) - adopt[i] + 1
  tau.fun <- function(l) 1 + 0.2 * l
  Y <- Y0 + D * tau.fun(ev)
  df <- data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
                   Y = c(Y), D = c(D))
  fit <- suppressMessages(fect(Y ~ D, data = df, index = c("id", "time"),
                               method = "gsynth", r = R, CV = FALSE,
                               se = TRUE, vartype = "conformal",
                               alpha = GAMMA, parallel = FALSE))
  post <- fit$time >= 1
  l <- fit$time[post]
  est <- fit$est.att[post, , drop = FALSE]
  cov.att <- (est[, "CI.lower"] <= tau.fun(l)) & (tau.fun(l) <= est[, "CI.upper"])
  cells <- D[, seq_len(Ntr)] == 1
  true.avg <- mean(tau.fun(ev[, seq_len(Ntr)][cells]))
  avg <- fit$est.avg[1, ]
  data.frame(
    cover.att = mean(cov.att), width.att = mean(est[, "CI.upper"] - est[, "CI.lower"]),
    cover.avg = (avg["CI.lower"] <= true.avg) & (true.avg <= avg["CI.upper"]),
    width.avg = avg["CI.upper"] - avg["CI.lower"],
    width.unit = mean(fit$conformal$unit$CI.upper - fit$conformal$unit$CI.lower),
    row.names = NULL
  )
}
res.multi <- list()
for (design in c("block", "staggered")) {
  seeds <- 9e5 + 1e4 * (design == "staggered") + seq_len(NREP)
  out <- parallel::mclapply(seeds, function(s) {
    tryCatch(sim_multi(design, s), error = function(e) NULL)
  }, mc.cores = CORES)
  ok <- !vapply(out, is.null, logical(1))
  d <- do.call(rbind, out[ok])
  res.multi[[design]] <- data.frame(
    design = design, reps = sum(ok), failed = sum(!ok),
    coverage.att = mean(d$cover.att), width.att = mean(d$width.att),
    coverage.avg = mean(d$cover.avg), width.avg = mean(d$width.avg),
    width.unit = mean(d$width.unit)
  )
  cat(sprintf("multi %s: coverage.att %.3f, coverage.avg %.3f\n", design,
              mean(d$cover.att), mean(d$cover.avg)))
}
tab.multi <- do.call(rbind, res.multi)
elapsed <- as.numeric(difftime(Sys.time(), t.start, units = "mins"))

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
lines <- c(
  "# Coverage of vartype = \"conformal\" (gsynth)",
  "",
  sprintf("Generated by `conformal_v2_coverage.R` on %s with fect %s (%s reps per cell, %d cores, %.1f minutes).",
          format(Sys.Date()), as.character(utils::packageVersion("fect")),
          NREP, CORES, elapsed),
  "",
  paste0("DGP: one treated unit, two factors, T = 30 (T0 = 20), J never-treated controls, ",
         "iid N(0, 1) errors or AR(1) errors (rho = 0.5, unit variance); r fixed at 2. ",
         "Nominal level 90% (gamma = 0.1). `coverage` is the date-specific coverage averaged ",
         "over the 10 post-treatment dates; `min.date` and `max.date` are its range over dates; ",
         "`bound` is the finite-sample bound 1 - (2q - 1)/N of v2; `mc.se` is the Monte Carlo ",
         "standard error of `coverage` (from the per-rep coverage over the 10 dates). `coverage.avg` and `width.avg` are for the ",
         "interval of the average effect over the 10 dates (extension beyond v2)."),
  "",
  "| J | errors | reps | failed | N | q | bound | coverage | min.date | max.date | mc.se | width | coverage.avg | width.avg |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"
)
for (k in seq_len(nrow(tab))) {
  x <- tab[k, ]
  lines <- c(lines, sprintf("| %d | %s | %d | %d | %d | %d | %s | %s | %s | %s | %s | %s | %s | %s |",
                            x$J, x$errors, x$reps, x$failed, x$N, x$q,
                            fmt(x$bound), fmt(x$coverage), fmt(x$min.date),
                            fmt(x$max.date), fmt(x$mc.se), fmt(x$width, 2),
                            fmt(x$coverage.avg), fmt(x$width.avg, 2)))
}
lines <- c(lines, "",
  "## Extension beyond v2: five treated units",
  "",
  paste0("J = 40, iid errors, r = 2 fixed, 90% level. `block`: all five adopt at date 21; ",
         "`staggered`: three at 17, two at 21. `coverage.att` is the coverage of the event-time ",
         "intervals (`est.att`) for the ATT at each event time, averaged over event times; ",
         "`coverage.avg` is for `est.avg`. `width.unit` is the mean width of the per-unit, ",
         "per-date intervals, for comparison: the averaged intervals stay about as wide as one ",
         "unit's (a single control's error does not average down). In these designs they cover ",
         "at or above the nominal level. That is a simulation result, not a guarantee: coverage ",
         "for aggregated targets is not established by the unit-level theorem, and a constructed ",
         "example (nine controls, 100 treated units, iid Bernoulli outcomes) has coverage 0.42 ",
         "at nominal 0.90."),
  "",
  "| design | reps | failed | coverage.att | width.att | coverage.avg | width.avg | width.unit |",
  "|---|---|---|---|---|---|---|---|")
for (k in seq_len(nrow(tab.multi))) {
  x <- tab.multi[k, ]
  lines <- c(lines, sprintf("| %s | %d | %d | %s | %s | %s | %s | %s |",
                            x$design, x$reps, x$failed, fmt(x$coverage.att),
                            fmt(x$width.att, 2), fmt(x$coverage.avg),
                            fmt(x$width.avg, 2), fmt(x$width.unit, 2)))
}
out.dir <- file.path("tests", "coverage-study", "results")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
writeLines(lines, file.path(out.dir, "conformal_v2_coverage.md"))
cat(lines, sep = "\n")
