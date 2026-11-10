## vartype = "conformal", the procedures added on 2026-11-05 after the paper
## (Liu and Xu 2026): pre-treatment diagnostics (conformal.diag), the block
## Jackknife+ for averages (conformal.block), the relabeling test of the sharp
## null of no effect (conformal.relabel), the message when cross-validation
## chose the tuning, the gap-plot caption and the print() lines. Each procedure
## is checked against a reference built from public fect() calls on a panel
## rearranged by hand, so the check does not reuse the code it tests.

.cf_sim <- function(J = 20, Ntr = 1, TT = 20, T0 = 15, r = 1, seed = 1,
                    adopt = NULL, tau = 1) {
  set.seed(seed)
  N <- J + Ntr
  F <- matrix(rnorm(TT * r), TT, r)
  L <- matrix(rnorm(N * r), N, r)
  Y0 <- outer(rnorm(TT), rep(1, N)) + outer(rep(1, TT), rnorm(N)) +
    F %*% t(L) + matrix(rnorm(TT * N), TT, N)
  if (is.null(adopt)) adopt <- rep(T0 + 1, Ntr)
  D <- matrix(0, TT, N)
  for (i in seq_len(Ntr)) D[seq_len(TT) >= adopt[i], i] <- 1
  Y <- Y0 + tau * D
  data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
             Y = c(Y), D = c(D))
}

.cf_fit <- function(df, ...) {
  suppressMessages(fect(Y ~ D, data = df, index = c("id", "time"),
                        se = TRUE, vartype = "conformal", parallel = FALSE,
                        CV = FALSE, ...))
}

## point fit through the public API, no inference
.cf_point <- function(df, ...) {
  suppressMessages(fect(Y ~ D, data = df, index = c("id", "time"),
                        se = FALSE, parallel = FALSE, CV = FALSE, ...))
}

## Jackknife+ endpoints by hand (the paper's rule)
.cf_jplus <- function(M, R, y, gamma) {
  y <- unname(y)
  N <- length(M) + 1
  qs <- 0:floor(N / 2)
  q <- max(qs[qs == 0 | (2 * qs - 1) / N <= gamma])
  if (q == 0) return(c(lo = -Inf, hi = Inf, N = N, q = 0))
  a <- sort(M - R)
  b <- sort(M + R)
  out <- c(y - b[N - q], y - a[q], N, q)
  names(out) <- c("lo", "hi", "N", "q")
  out
}

## move unit i's adoption to date s (column D of the long data)
.cf_move <- function(df, i, s) {
  sel <- df$id == i
  df$D[sel] <- as.numeric(df$time[sel] >= s)
  df
}

## ---- the CV message ---------------------------------------------------------

test_that("a conformal fit with cross-validated tuning prints a message", {
  df <- .cf_sim(J = 15, Ntr = 1, seed = 2)
  expect_message(
    fect(Y ~ D, data = df, index = c("id", "time"), method = "gsynth",
         r = c(0, 2), CV = TRUE, se = TRUE, vartype = "conformal",
         parallel = FALSE),
    "chosen by cross-validation on the never-treated units' outcomes only.*refit with r = [0-9]+ and CV = FALSE"
  )
  ## ife's cross-validation also scores the treated units' pre-period cells
  expect_message(
    fect(Y ~ D, data = df, index = c("id", "time"), method = "ife",
         r = c(0, 2), CV = TRUE, se = TRUE, vartype = "conformal",
         parallel = FALSE),
    "chosen by cross-validation on a fit that includes the treated units' pre-treatment outcomes"
  )
  ## fixed tuning: no such message
  msgs <- character(0)
  withCallingHandlers(
    fect(Y ~ D, data = df, index = c("id", "time"), method = "gsynth",
         r = 1, CV = FALSE, se = TRUE, vartype = "conformal",
         parallel = FALSE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("chosen by cross-validation", msgs)))
  ## a single r turns CV off by default, and the message stays away
  msgs <- character(0)
  withCallingHandlers(
    fect(Y ~ D, data = df, index = c("id", "time"), method = "ife", r = 1,
         se = TRUE, vartype = "conformal", parallel = FALSE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("chosen by cross-validation", msgs)))
})

## ---- argument checks --------------------------------------------------------

test_that("the new arguments are checked", {
  df <- .cf_sim(J = 10, Ntr = 1, seed = 3)
  expect_error(.cf_fit(df, method = "gsynth", r = 1, conformal.diag = -1),
               "conformal.diag must be a whole number")
  expect_error(.cf_fit(df, method = "gsynth", r = 1, conformal.relabel = 2.5),
               "conformal.relabel must be a whole number")
  expect_error(.cf_fit(df, method = "gsynth", r = 1, conformal.block = NA),
               "conformal.block must be TRUE or FALSE")
  expect_warning(
    suppressMessages(fect(Y ~ D, data = df, index = c("id", "time"),
                          method = "gsynth", r = 1, CV = FALSE, se = FALSE,
                          conformal.diag = 3)),
    "apply only with se = TRUE and vartype = \"conformal\""
  )
  ## too few complete controls for one block of K
  df2 <- .cf_sim(J = 3, Ntr = 4, seed = 4)
  expect_error(.cf_fit(df2, method = "gsynth", r = 0, conformal.block = TRUE),
               "needs at least as many never-treated controls")
})

## ---- pre-treatment diagnostics ----------------------------------------------

test_that("diagnostics equal the Jackknife+ with the adoption moved (gsynth)", {
  df <- .cf_sim(J = 20, Ntr = 3, T0 = 14, seed = 6, r = 1)
  m <- 3
  fit <- .cf_fit(df, method = "gsynth", r = 1, conformal.diag = m)
  du <- fit$conformal$diag$unit
  expect_equal(nrow(du), 3 * m)
  expect_setequal(du$event.time, (1 - m):0)
  for (e in (1 - m):0) {
    s <- 15 + e - 1
    ref.df <- df
    for (i in 1:3) ref.df <- .cf_move(ref.df, i, s)
    ref <- .cf_fit(ref.df, method = "gsynth", r = 1)
    ru <- ref$conformal$unit[ref$conformal$unit$time == s, ]
    dd <- du[du$event.time == e, ]
    dd <- dd[order(dd$unit), ]
    ru <- ru[order(ru$unit), ]
    expect_equal(dd$time, rep(s, 3))
    for (v in c("CI.lower", "CI.upper", "p.value", "N", "q", "n.failed")) {
      expect_equal(dd[[v]], ru[[v]], tolerance = 1e-8, info = paste(e, v))
    }
    ## the point estimate is the prediction error from a fit with all controls
    expect_equal(dd$eff, ru$eff, tolerance = 1e-6, info = paste(e, "eff"))
  }
})

test_that("diagnostics equal the Jackknife+ with the adoption moved (ife)", {
  df <- .cf_sim(J = 12, Ntr = 1, T0 = 14, TT = 18, seed = 7, r = 1)
  fit <- .cf_fit(df, method = "ife", r = 1, conformal.diag = 2)
  du <- fit$conformal$diag$unit
  for (e in -1:0) {
    s <- 15 + e - 1
    ref <- .cf_fit(.cf_move(df, 1, s), method = "ife", r = 1)
    ru <- ref$conformal$unit[ref$conformal$unit$time == s, ]
    dd <- du[du$event.time == e, ]
    for (v in c("CI.lower", "CI.upper", "p.value", "N", "q")) {
      expect_equal(dd[[v]], ru[[v]], tolerance = 1e-8, info = paste(e, v))
    }
  }
})

test_that("diagnostics fill the pre-treatment rows of est.att, and only those", {
  df <- .cf_sim(J = 20, Ntr = 1, T0 = 14, seed = 8)
  fit <- .cf_fit(df, method = "gsynth", r = 1, alpha = 0.1,
                 conformal.diag = 4)
  ea <- fit$est.att
  et <- as.numeric(rownames(ea))
  pre.diag <- et %in% (-3:0)
  expect_true(all(is.finite(ea[pre.diag, "CI.lower"])))
  expect_true(all(is.na(ea[et < -3, "CI.lower"])))
  ## one treated unit: the event-time row is that unit's diagnostic
  du <- fit$conformal$diag$unit
  for (e in -3:0) {
    row <- ea[et == e, ]
    u <- du[du$event.time == e, ]
    expect_equal(unname(row["ATT"]), u$eff)
    expect_equal(unname(row[c("CI.lower", "CI.upper", "p.value")]),
                 c(u$CI.lower, u$CI.upper, u$p.value))
  }
  ## the 90% rows use 2 * alpha
  da <- fit$conformal$diag$att
  expect_equal(unname(fit$est.att90[match(da$event.time, et), "CI.lower"]),
               da$CI.lower.90)
  ## without diagnostics the pre-treatment rows have no interval
  fit0 <- .cf_fit(df, method = "gsynth", r = 1, alpha = 0.1)
  expect_true(all(is.na(fit0$est.att[et <= 0, "CI.lower"])))
  expect_null(fit0$conformal$diag)
  ## the post-treatment rows do not change
  expect_equal(fit$est.att[et >= 1, ], fit0$est.att[et >= 1, ])
})

test_that("a diagnostic date without earlier observed history is skipped", {
  df <- .cf_sim(J = 15, Ntr = 1, T0 = 3, TT = 10, seed = 9)
  fit <- .cf_fit(df, method = "gsynth", r = 0, conformal.diag = 4,
                 min.T0 = 3)
  ## dates 0, -1 have history; date -2 is time 1 with none before it
  expect_setequal(fit$conformal$diag$unit$event.time, c(-1, 0))
})

## ---- block Jackknife+ ---------------------------------------------------------

## reference: for block g, mark its controls as treated with their matched
## adopters' dates and fit gsynth through fect(); the never-treated rest are
## the controls
.cf_block_ref <- function(df, fit, adopt, gamma, r) {
  part <- fit$conformal$block$partition
  G <- max(part$block)
  tr <- unique(part$matched.to)
  et <- as.numeric(rownames(fit$est.att))
  Ms <- Rs <- NULL
  for (g in seq_len(G)) {
    pg <- part[part$block == g, ]
    d2 <- df
    for (k in seq_len(nrow(pg))) {
      d2 <- .cf_move(d2, pg$unit[k], adopt[match(pg$matched.to[k], tr)])
    }
    pf <- .cf_point(d2, method = "gsynth", r = r)
    Yct <- pf$Y.ct
    Yd <- pf$Y.dat
    ids <- as.numeric(pf$id)
    M.g <- R.g <- c()
    for (l in et[et >= 1]) {
      pr <- er <- c()
      for (k in seq_along(tr)) {
        t <- adopt[k] + l - 1
        if (t > nrow(Yd)) next
        pr <- c(pr, Yct[t, ids == tr[k]])
        cj <- pg$unit[pg$matched.to == tr[k]]
        er <- c(er, Yd[t, ids == cj] - Yct[t, ids == cj])
      }
      M.g <- c(M.g, mean(pr))
      R.g <- c(R.g, abs(mean(er)))
    }
    Ms <- cbind(Ms, M.g)
    Rs <- cbind(Rs, R.g)
  }
  list(M = Ms, R = Rs)
}

test_that("block Jackknife+ matches a reference built from fect() calls", {
  df <- .cf_sim(J = 22, Ntr = 3, T0 = 14, seed = 10, r = 1)
  gamma <- 0.2
  fit <- .cf_fit(df, method = "gsynth", r = 1, alpha = gamma,
                 conformal.block = TRUE, seed = 11)
  b <- fit$conformal$block
  expect_equal(b$G, 7)
  expect_equal(b$K, 3)
  expect_equal(b$N, 8)
  expect_equal(b$n.training, 1)
  expect_equal(nrow(b$partition), 21)
  expect_false(any(duplicated(b$partition$unit)))
  expect_true(all(b$partition$unit > 3))
  expect_identical(fit$conformal$aggregate, "block")
  ref <- .cf_block_ref(df, fit, adopt = rep(15, 3), gamma = gamma, r = 1)
  et <- as.numeric(rownames(fit$est.att))
  post <- which(et >= 1)
  pf <- .cf_point(df, method = "gsynth", r = 1)
  ybar <- rowMeans(pf$Y.dat[15:20, match(1:3, pf$id), drop = FALSE])
  for (ll in seq_along(post)) {
    h <- .cf_jplus(ref$M[ll, ], ref$R[ll, ], ybar[ll], gamma)
    expect_equal(unname(fit$est.att[post[ll], c("CI.lower", "CI.upper")]),
                 unname(h[c("lo", "hi")]), tolerance = 1e-6,
                 info = paste("event time", et[post[ll]]))
    expect_equal(fit$conformal$att$N[post[ll]], 8)
  }
})

test_that("block Jackknife+ under staggered adoption matches the reference", {
  adopt <- c(13, 13, 16)
  df <- .cf_sim(J = 19, Ntr = 3, TT = 20, seed = 12, r = 1, adopt = adopt)
  gamma <- 0.25
  fit <- .cf_fit(df, method = "gsynth", r = 1, alpha = gamma,
                 conformal.block = TRUE, seed = 13)
  ref <- .cf_block_ref(df, fit, adopt = adopt, gamma = gamma, r = 1)
  et <- as.numeric(rownames(fit$est.att))
  post <- which(et >= 1)
  pf <- .cf_point(df, method = "gsynth", r = 1)
  for (ll in seq_along(post)) {
    l <- et[post[ll]]
    ys <- c()
    for (k in 1:3) {
      t <- adopt[k] + l - 1
      if (t <= 20) ys <- c(ys, pf$Y.dat[t, match(k, pf$id)])
    }
    h <- .cf_jplus(ref$M[ll, ], ref$R[ll, ], mean(ys), gamma)
    expect_equal(unname(fit$est.att[post[ll], c("CI.lower", "CI.upper")]),
                 unname(h[c("lo", "hi")]), tolerance = 1e-6,
                 info = paste("event time", l))
  }
  ## event times reached by one adopter only still use blocks of size K
  expect_equal(unname(fit$est.att[et == 6, "count"]), 2)
})

test_that("with one treated unit, gsynth's block Jackknife+ is the pooled one", {
  df <- .cf_sim(J = 18, Ntr = 1, seed = 14)
  f0 <- .cf_fit(df, method = "gsynth", r = 1, alpha = 0.1)
  f1 <- .cf_fit(df, method = "gsynth", r = 1, alpha = 0.1,
                conformal.block = TRUE, seed = 1)
  expect_equal(f1$conformal$block$G, 18)
  expect_equal(f1$est.att, f0$est.att)
  expect_equal(f1$est.avg, f0$est.avg)
  expect_equal(f1$conformal$unit, f0$conformal$unit)
})

test_that("the partition follows the seed, and q = 0 gives the whole line", {
  df <- .cf_sim(J = 20, Ntr = 2, seed = 15)
  fa <- .cf_fit(df, method = "gsynth", r = 1, alpha = 0.2,
                conformal.block = TRUE, seed = 3)
  fb <- .cf_fit(df, method = "gsynth", r = 1, alpha = 0.2,
                conformal.block = TRUE, seed = 3)
  expect_identical(fa$conformal$block$partition, fb$conformal$block$partition)
  expect_equal(fa$est.att, fb$est.att)
  ## G + 1 = 11 < 1/0.05: no finite interval
  fz <- .cf_fit(df, method = "gsynth", r = 1, alpha = 0.05,
                conformal.block = TRUE, seed = 3)
  et <- as.numeric(rownames(fz$est.att))
  expect_true(all(is.infinite(fz$est.att[et >= 1, "CI.lower"])))
  expect_equal(fz$conformal$block$q, 0L)
})

test_that("controls with missing outcomes stay training units in the block fits", {
  df <- .cf_sim(J = 14, Ntr = 2, seed = 16)
  df <- df[!(df$id == 10 & df$time == 3), ]
  fit <- .cf_fit(df, method = "gsynth", r = 1, alpha = 0.25,
                 conformal.block = TRUE, seed = 2)
  expect_false(10 %in% fit$conformal$block$partition$unit)
  expect_equal(fit$conformal$block$G, 6)
})

test_that("block Jackknife+ and diagnostics run with ife and mc", {
  df <- .cf_sim(J = 12, Ntr = 2, TT = 16, T0 = 12, seed = 17)
  for (m in c("ife", "mc")) {
    a <- list(df, method = m, alpha = 0.3, conformal.block = TRUE,
              conformal.diag = 2, seed = 4)
    if (m == "ife") a$r <- 1 else a$lambda <- 0.05
    fit <- do.call(.cf_fit, a)
    et <- as.numeric(rownames(fit$est.att))
    expect_true(all(is.finite(fit$est.att[et >= -1, "CI.lower"])), info = m)
    expect_equal(fit$conformal$block$G, 6, info = m)
    expect_equal(nrow(fit$conformal$diag$unit), 4, info = m)
  }
})

## ---- relabeling test --------------------------------------------------------

test_that("relabeling statistics match fect() fits on the relabeled panel", {
  df <- .cf_sim(J = 15, Ntr = 2, seed = 18, tau = 2)
  fit <- .cf_fit(df, method = "gsynth", r = 1, conformal.relabel = 19,
                 seed = 5)
  rl <- fit$conformal$relabel
  expect_equal(rl$B, 19L)
  expect_equal(dim(rl$units), c(2L, 19L))
  ## observed statistic: the point fit's average effect
  pf <- .cf_point(df, method = "gsynth", r = 1)
  expect_equal(rl$ATT.avg, pf$att.avg, tolerance = 1e-6)
  ## three draws by hand: give unit units[k, b] the adoption of treated unit k
  for (b in 1:3) {
    d2 <- df
    d2$D <- 0
    for (k in 1:2) d2 <- .cf_move(d2, rl$units[k, b], 16)
    rf <- .cf_point(d2, method = "gsynth", r = 1)
    expect_equal(rl$draws.avg[b], rf$att.avg, tolerance = 1e-6,
                 info = paste("draw", b))
    expect_equal(unname(rl$paths[, b]),
                 unname(rf$att[rf$time >= 1]),
                 tolerance = 1e-6, info = paste("path", b))
  }
  ## p-value
  expect_equal(rl$p.value,
               (1 + sum(abs(rl$draws.avg) >= abs(rl$ATT.avg))) / 20)
  for (l in seq_len(nrow(rl$att))) {
    expect_equal(rl$att$p.value[l],
                 (1 + sum(abs(rl$paths[l, ]) >= abs(rl$att$ATT[l]))) / 20)
  }
  ## same seed, same draws
  fit2 <- .cf_fit(df, method = "gsynth", r = 1, conformal.relabel = 19,
                  seed = 5)
  expect_identical(fit2$conformal$relabel$units, rl$units)
})

test_that("the relabeling p-value is valid under exchangeable units", {
  skip_on_cran()
  ## no effect, iid units, one treated unit: P(p <= 0.1) <= 0.1
  pv <- vapply(1:60, function(s) {
    df <- .cf_sim(J = 14, Ntr = 1, TT = 12, T0 = 8, seed = 100 + s, tau = 0,
                  r = 1)
    .cf_fit(df, method = "gsynth", r = 1, conformal.relabel = 19,
            seed = s)$conformal$relabel$p.value
  }, numeric(1))
  expect_true(all(pv >= 1 / 20 & pv <= 1))
  expect_lte(mean(pv <= 0.1), 0.2)
})

test_that("a failed relabeling fit counts as extreme and warns", {
  df <- .cf_sim(J = 12, Ntr = 1, seed = 19, tau = 3)
  orig <- fect:::.conformal_component_fit
  st <- new.env()
  st$n <- 0L
  ## the relabeling fits are the last 1 + B calls; fail the last two
  mock <- function(...) {
    st$n <- st$n + 1L
    if (isTRUE(st$fail) && st$n > st$cut) return(NULL)
    orig(...)
  }
  ok <- .cf_fit(df, method = "gsynth", r = 1, conformal.relabel = 9, seed = 1)
  testthat::local_mocked_bindings(.conformal_component_fit = mock,
                                  .package = "fect")
  st$fail <- FALSE
  .cf_fit(df, method = "gsynth", r = 1, conformal.relabel = 9, seed = 1)
  total <- st$n
  st$n <- 0L
  st$fail <- TRUE
  st$cut <- total - 2L
  expect_warning(
    fit <- .cf_fit(df, method = "gsynth", r = 1, conformal.relabel = 9,
                   seed = 1),
    "2 of 9 relabeling fits failed"
  )
  rl <- fit$conformal$relabel
  expect_equal(rl$n.failed, 2L)
  expect_equal(rl$p.value,
               (1 + 2 + sum(abs(rl$draws.avg) >= abs(rl$ATT.avg),
                            na.rm = TRUE)) / 10)
  expect_gte(rl$p.value, ok$conformal$relabel$p.value)
})

## ---- print() and the gap plot --------------------------------------------------

test_that("print() and the gap plot describe the new procedures", {
  df <- .cf_sim(J = 20, Ntr = 2, seed = 20)
  fit <- .cf_fit(df, method = "gsynth", r = 1, alpha = 0.2,
                 conformal.block = TRUE, conformal.diag = 2,
                 conformal.relabel = 9, seed = 6)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Aggregates:    block Jackknife+, 10 control blocks of 2 (N = 11",
               fixed = TRUE)
  expect_match(out, "Diagnostics:   last 2 pre-treatment dates; [0-9]+ of 4")
  expect_match(out, "Relabeling:    sharp null of no effect, 9 draws, p = ")
  expect_false(grepl("no coverage guarantee", out))
  p <- plot(fit)
  expect_match(p$labels$caption, "Block Jackknife+ intervals: 10 control blocks of 2",
               fixed = TRUE)
  expect_match(p$labels$caption, "diagnostic intervals at the last 2 dates",
               fixed = TRUE)
  ## pooled aggregates over two treated units: the caveat
  f0 <- .cf_fit(df, method = "gsynth", r = 1, alpha = 0.2)
  expect_match(plot(f0)$labels$caption, "no coverage guarantee", fixed = TRUE)
  ## one treated unit, no options: no caption
  f1 <- .cf_fit(.cf_sim(J = 15, Ntr = 1, seed = 21), method = "gsynth", r = 1)
  expect_null(plot(f1)$labels$caption)
})

test_that("the new tables carry unit and time labels", {
  df <- .cf_sim(J = 12, Ntr = 2, seed = 22)
  df$id <- sprintf("u%02d", df$id)
  df$time <- df$time + 2000
  fit <- .cf_fit(df, method = "gsynth", r = 1, alpha = 0.3,
                 conformal.diag = 2, conformal.block = TRUE,
                 conformal.relabel = 5, seed = 1)
  du <- fit$conformal$diag$unit
  expect_setequal(du$unit, c("u01", "u02"))
  expect_setequal(as.numeric(du$time), c(2014, 2015))
  bp <- fit$conformal$block$partition
  expect_true(all(bp$unit %in% sprintf("u%02d", 3:14)))
  expect_setequal(bp$matched.to, c("u01", "u02"))
  expect_true(all(fit$conformal$relabel$units %in% unique(df$id)))
  expect_equal(dim(fit$conformal$relabel$units), c(2L, 5L))
})
