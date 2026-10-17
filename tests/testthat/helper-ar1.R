## Helpers for the serial = "ar1" and cv.method = "pretreatment" tests
## (test-cv-treated.R, test-ar1-joint.R, test-ar1-downstream.R).
##
## ar1_panel() is the test generator of the validation spec: 
##   Y_it(0) = alpha_i + xi_t + lambda_i' F_t + beta X_it + e_it,
## alpha_i, xi_t ~ N(0, 1); F is T x 2 iid N(0, 1); lambda_i ~ N(0, I_2);
## e_it stationary AR(1) with coefficient rho and variance 1
## (e_i1 ~ N(0, 1), e_it = rho e_i,t-1 + sqrt(1 - rho^2) z_it); X_it stationary
## AR(1) 0.7, variance 1 (when present), beta = 1. Units 1..Nco are never
## treated, the next Ntr are treated; D = 1 for t > T0 (common adoption) or for
## t > T0[i] (staggered). A constant effect tau is added to treated post cells.
## Draw order from the seed: alpha, xi, F, lambda, then the errors, then X.
## noiseless = TRUE sets e_it = rho^(t - 1) e_i1 (no innovations after period 1).

ar1_panel <- function(Nco, Ntr, T0, q, rho, seed, cov = FALSE, tau = 0,
                      noiseless = FALSE, T0.unit = NULL) {
  set.seed(seed)
  N <- Nco + Ntr; T <- if (is.null(T0.unit)) T0 + q else q
  if (!is.null(T0.unit)) T0 <- NA
  alpha <- rnorm(N); xi <- rnorm(T)
  Fm <- matrix(rnorm(T * 2), T, 2); Lam <- matrix(rnorm(N * 2), N, 2)
  ar1_paths <- function(n, coef, innov = TRUE) {
    E <- matrix(0, T, n); E[1, ] <- rnorm(n)
    for (t in 2:T) E[t, ] <- if (innov) coef * E[t - 1, ] + sqrt(1 - coef^2) * rnorm(n) else coef * E[t - 1, ]
    E
  }
  E <- ar1_paths(N, rho, innov = !noiseless)
  Y <- outer(xi, rep(1, N)) + outer(rep(1, T), alpha) + Fm %*% t(Lam) + E
  X <- NULL
  if (cov) { X <- ar1_paths(N, 0.7); Y <- Y + X }
  t0 <- if (is.null(T0.unit)) rep(T0, Ntr) else T0.unit
  D <- matrix(0L, T, N)
  for (j in seq_len(Ntr)) D[seq_len(T) > t0[j], Nco + j] <- 1L
  Y <- Y + tau * D
  d <- data.frame(id = rep(seq_len(N), each = T), time = rep(seq_len(T), N),
                  y = as.vector(Y), d = as.vector(D))
  if (cov) d$x <- as.vector(X)
  attr(d, "Nco") <- Nco; attr(d, "Ntr") <- Ntr; attr(d, "T0") <- t0; attr(d, "T") <- T
  d
}

## P2 balanced, P3 staggered, P5 covariate, P6 noiseless variants.
ar1_P2 <- function(...) ar1_panel(30, 5, 15, 5, rho = 0.6, seed = 101, ...)
ar1_P3 <- function(...) ar1_panel(30, 6, NA, 20, rho = 0.6, seed = 102,
                                  T0.unit = c(11, 11, 13, 13, 15, 15), ...)
ar1_P5 <- function(...) ar1_panel(30, 5, 15, 5, rho = 0.6, seed = 104, cov = TRUE, ...)

## Drop cells at random from t >= 2 (probability p each), redrawing until every
## unit keeps >= 8 observed pretreatment cells and every period keeps an
## observed control; then force the last pretreatment cell (t = T0) of `force.units`
## (treated units) missing. Rows are removed (not set to NA).
ar1_drop_cells <- function(d, seed, p = 0.08, force.units = c(31, 32), min.pre = 8) {
  Nco <- attr(d, "Nco"); T0 <- attr(d, "T0"); Tn <- attr(d, "T")
  set.seed(seed)
  repeat {
    drop <- d$time >= 2 & runif(nrow(d)) < p
    k <- d[!drop, ]
    npre <- tapply(k$time <= max(T0) & k$d == 0, k$id, sum)
    ctrl_ok <- all(tabulate(k$time[k$id <= Nco], Tn) > 0)
    if (all(npre >= min.pre) && ctrl_ok) break
  }
  k <- k[!(k$id %in% force.units & k$time == max(T0)), ]
  attr(k, "Nco") <- Nco; attr(k, "Ntr") <- attr(d, "Ntr"); attr(k, "T0") <- T0; attr(k, "T") <- Tn
  k
}
ar1_P4 <- function() ar1_drop_cells(ar1_P2(), seed = 103)

## Reference quasi-differencing on the observed time grid: V = v - x * v_{previous grid period}.
ar1_hand_filter <- function(d, x, vars = c("y", "x"), id = "id", time = "time") {
  vars <- intersect(vars, names(d))
  tg <- sort(unique(d[[time]])); d$.pos <- match(d[[time]], tg)
  lag <- d[, c(id, ".pos", vars)]; lag$.pos <- lag$.pos + 1L
  names(lag)[-(1:2)] <- paste0(vars, ".lag")
  m <- merge(d, lag, by = c(id, ".pos"))
  for (v in vars) m[[v]] <- m[[v]] - x * m[[paste0(v, ".lag")]]
  m <- m[order(m[[id]], m[[time]]), setdiff(names(d), ".pos")]
  rownames(m) <- NULL
  m
}


## P1 (prototype reference panel) in long format from the stored matrix.
ar1_p1_long <- function(ref) {
  Y <- ref$P1$Y
  data.frame(id = rep(1:35, each = 20), time = rep(1:20, 35), y = as.vector(Y),
             d = as.integer(rep(1:35, each = 20) > 30 & rep(1:20, 35) > 15))
}

## Random-walk-error panel (e_i1 ~ N(0, 1), innovation SD 0.3122), two factors, no effect.
ar1_rw_panel <- function(Nco, Ntr, T0, q, seed) {
  set.seed(seed); T <- T0 + q; N <- Nco + Ntr
  al <- rnorm(N); xi <- rnorm(T); Fm <- matrix(rnorm(T * 2), T, 2); Lam <- matrix(rnorm(N * 2), N, 2)
  E <- matrix(0, T, N); E[1, ] <- rnorm(N); for (t in 2:T) E[t, ] <- E[t - 1, ] + 0.3122 * rnorm(N)
  Y <- outer(xi, rep(1, N)) + outer(rep(1, T), al) + Fm %*% t(Lam) + E
  D <- matrix(0L, T, N); D[(T0 + 1):T, Nco + seq_len(Ntr)] <- 1L
  data.frame(id = rep(seq_len(N), each = T), time = rep(seq_len(T), N), y = as.vector(Y), d = as.vector(D))
}

## The motivating design: 30 controls, 5 treated, T0 = 15, q = 5, unit and time effects, two
## factors (a standardized linear trend and white noise), loadings N(0, 1), AR(1) 0.7 errors, no effect.
ar1_mot_panel <- function(seed, rho = 0.7) {
  set.seed(seed); N <- 35; T <- 20; T0 <- 15
  f1 <- (1:T - (T + 1) / 2); f1 <- f1 / sd(f1); Fm <- cbind(f1, rnorm(T))
  Lam <- matrix(rnorm(N * 2), N, 2); al <- rnorm(N); xi <- rnorm(T)
  E <- matrix(0, T, N); E[1, ] <- rnorm(N); for (t in 2:T) E[t, ] <- rho * E[t - 1, ] + sqrt(1 - rho^2) * rnorm(N)
  Y <- outer(xi, rep(1, N)) + outer(rep(1, T), al) + Fm %*% t(Lam) + E
  D <- matrix(0L, T, N); D[(T0 + 1):T, 31:35] <- 1L
  data.frame(id = rep(seq_len(N), each = T), time = rep(seq_len(T), N), y = as.vector(Y), d = as.vector(D))
}

## Quiet fect() call.
ar1_fect <- function(d, formula = y ~ d, ..., quiet = TRUE) {
  args <- list(formula = formula, data = d, index = c("id", "time"), method = "gsynth", force = "two-way",
               se = FALSE, parallel = FALSE, tol = 1e-10)
  dots <- list(...); args[names(dots)] <- dots
  f <- function() do.call(fect::fect, args)
  if (quiet) suppressMessages(f()) else f()
}

## Reference for the serial = "none" fit on the panel without its first period.
ar1_fit_levels <- function(d, r, formula = y ~ d, min.T0 = 4, ...) {
  d1 <- d[d$time > min(d$time), ]
  args <- list(d1, formula = formula, CV = FALSE, r = r, ...)
  if (!is.null(min.T0)) args$min.T0 <- min.T0
  do.call(ar1_fect, args)
}

## Absolute (or relative) closeness, with the observed worst difference reported
## and, when AR1_TEST_LOG is set to a file, appended to it (audit evidence).
expect_close <- function(actual, expected, tol, label = "", rel = FALSE) {
  a <- as.numeric(actual); e <- as.numeric(expected)
  ok_na <- identical(is.na(a), is.na(e))
  keep <- !is.na(a) & !is.na(e)
  err <- if (rel) max(abs(a[keep] - e[keep]) / pmax(abs(e[keep]), 1e-300), 0) else max(abs(a[keep] - e[keep]), 0)
  lg <- Sys.getenv("AR1_TEST_LOG", "")
  if (nzchar(lg)) {
    cat(sprintf("%s\t%s\t%s\t%.3e\t%.1e\t%s\n", testthat::get_reporter()$.context %||% "", label,
                if (length(a) == 1) sprintf("%.12g", a) else sprintf("n=%d", length(a)),
                err, tol, if (rel) "rel" else "abs"), file = lg, append = TRUE)
  }
  testthat::expect_true(ok_na, info = paste(label, ": NA pattern differs"))
  testthat::expect_true(err <= tol, info = sprintf("%s: max %s difference %.3e exceeds %.1e", label,
                                                    if (rel) "relative" else "absolute", err, tol))
  invisible(err)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## Run code in a clean subprocess against another fect library (the base build),
## returning the object it saves. `expr` is a character string of R code that
## must assign its result to `out`. Used by the byte-identical-default tests.
ar1_base_lib <- function() {
  p <- Sys.getenv("FECT_BASE_LIB", "")
  if (nzchar(p) && dir.exists(file.path(p, "fect"))) p else ""
}
ar1_run_in_lib <- function(lib, code) {
  out <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R")
  writeLines(c("suppressPackageStartupMessages(library(fect))",
               "msgs <- character(0)",
               "out <- withCallingHandlers({", code, "}, message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart('muffleMessage') })",
               sprintf("saveRDS(list(out = out, msgs = msgs), %s)", deparse(out))), script)
  rs <- file.path(R.home("bin"), "Rscript")
  st <- system2(rs, c("--vanilla", shQuote(script)), env = c(paste0("R_LIBS=", lib)), stdout = TRUE, stderr = TRUE)
  if (!file.exists(out)) stop("subprocess failed:\n", paste(utils::tail(st, 20), collapse = "\n"))
  readRDS(out)
}

## Evaluate `expr`, collecting (and muffling) its warnings; returns list(value, warnings = character()).
ar1_collect <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  list(value = v, warnings = w)
}
## Write a line of informational evidence to the AR1_TEST_LOG file when set (and always to the message stream).
ar1_note <- function(...) {
  txt <- paste0(...)
  message(txt)
  lg <- Sys.getenv("AR1_TEST_LOG", "")
  if (nzchar(lg)) cat("NOTE\t", txt, "\n", sep = "", file = lg, append = TRUE)
}

## ---- Phase 3c (re-estimating x in every bootstrap draw): helpers for test-ar1-boot-x.R ----

## T2 generator (test-spec Amendments 3, seed 11): N_co = 20, N_tr = 3, T = 14, treatment on the last 4 periods of the
## treated units; x0 = 0.5; xi_t, alpha_i, f_t, lambda_i, Y_i1 ~ N(0, 1) (drawn in that order);
## m_it = 0.3 + xi_t + alpha_i + lambda_i f_t; Y0_it = x0 Y0_i,t-1 + m_it (t >= 2); observed = Y0 + tau D (the effect is added
## to the observed outcome only, not carried in the recursion). The filtered model at x0 is exactly rank 1 plus two-way
## effects: there is no noise.
ar1_noiseless_3c <- function(seed = 11, Nco = 20, Ntr = 3, T = 14, q = 4, x0 = 0.5, tau = 1) {
  set.seed(seed); N <- Nco + Ntr
  xi <- rnorm(T); alpha <- rnorm(N); f <- rnorm(T); lam <- rnorm(N); Y1 <- rnorm(N)
  M <- 0.3 + outer(xi, rep(1, N)) + outer(rep(1, T), alpha) + outer(f, lam)
  Y0 <- matrix(0, T, N); Y0[1, ] <- Y1
  for (t in 2:T) Y0[t, ] <- x0 * Y0[t - 1, ] + M[t, ]
  D <- matrix(0L, T, N); D[(T - q + 1):T, Nco + seq_len(Ntr)] <- 1L
  Y <- Y0 + tau * D
  data.frame(id = rep(seq_len(N), each = T), time = rep(seq_len(T), N), y = as.vector(Y), d = as.vector(D))
}

## Largest difference between two slots, relative to max(1, max |b|); Inf when the structure or the NA pattern differs.
ar1_slot_diff <- function(a, b) {
  if (is.null(a) || is.null(b)) return(if (is.null(a) && is.null(b)) 0 else Inf)
  if (is.list(a)) {
    if (!is.list(b) || length(a) != length(b)) return(Inf)
    return(max(0, unlist(Map(ar1_slot_diff, a, b))))
  }
  if (is.numeric(a) || is.logical(a)) {
    if (!(is.numeric(b) || is.logical(b)) || !identical(dim(a), dim(b)) || length(a) != length(b)) return(Inf)
    if (!identical(is.na(a), is.na(b))) return(Inf)
    k <- !is.na(a); if (!any(k)) return(0)
    sc <- max(1, max(abs(as.numeric(b[k]))))
    return(max(abs(as.numeric(a[k]) - as.numeric(b[k]))) / sc)
  }
  if (identical(a, b)) 0 else Inf
}

## Every slot but `call` equal with tolerance 0.
ar1_fit_identical <- function(a, b) {
  nm <- setdiff(union(names(a), names(b)), "call")
  bad <- nm[!vapply(nm, function(s) isTRUE(all.equal(a[[s]], b[[s]], tolerance = 0)), NA)]
  bad
}

## Log an evidence line (label, value, error, tolerance) when AR1_TEST_LOG is set.
ar1_log <- function(label, value, err, tol) {
  lg <- Sys.getenv("AR1_TEST_LOG", "")
  if (nzchar(lg)) cat(sprintf("3c\t%s\t%s\t%.3e\t%.1e\tabs\n", label, value, err, tol), file = lg, append = TRUE)
}

## Levels matrices on the full time grid for the low-level search seam fect:::.fect_ar1_search_x (NA = unobserved).
ar1_levels_arrays <- function(d, covar = FALSE) {
  ids <- sort(unique(d$id)); tg <- sort(unique(d$time))
  Y <- D <- matrix(NA_real_, length(tg), length(ids))
  i <- match(d$time, tg); j <- match(d$id, ids)
  Y[cbind(i, j)] <- d$y; D[cbind(i, j)] <- d$d
  X <- NULL
  if (covar) { X <- array(NA_real_, c(length(tg), length(ids), 1)); X[cbind(i, j, 1L)] <- d$x }
  list(Y = Y, D = D, X = X)
}

## TRUE when the loaded fect is a source tree (load_all), not an installed library: parallel workers load the installed one.
ar1_is_dev_load <- function() {
  p <- tryCatch(getNamespaceInfo("fect", "path"), error = function(e) "")
  nzchar(p) && length(list.files(file.path(p, "R"), pattern = "\\.R$")) > 0
}
