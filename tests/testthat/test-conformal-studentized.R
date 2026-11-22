## The studentized score for full conformal inference and the Jackknife+ for
## target-date estimators (conformal_sc(score = "studentized"), and
## fect(method = "fe", vartype = "conformal", conformal.score = "studentized")).
## Each unit's absolute prediction error is divided by its pre-treatment root
## mean squared residual from the fit that holds it out (Liu and Xu 2026).
## The checks here are by hand and by direct refits, not against the code.

stud_panel <- function(seed = 11, J = 9, T0 = 8, T1 = 3, effect = 0) {
  set.seed(seed)
  sdv <- c(0.3, 1, 0.5, 2, 0.8, 1.5, 0.4, 1.2, 0.7, 1.8)[seq_len(J + 1L)]
  TT <- T0 + T1
  fe.u <- rnorm(J + 1L)
  fe.t <- rnorm(TT)
  Y <- outer(fe.t, fe.u, "+") + sapply(sdv, function(s) rnorm(TT, sd = s))
  Y[(T0 + 1L):TT, 1L] <- Y[(T0 + 1L):TT, 1L] + effect
  d <- data.frame(id = rep(seq_len(J + 1L), each = TT),
                  time = rep(seq_len(TT), J + 1L), Y = as.vector(Y))
  d$D <- as.integer(d$id == 1L & d$time > T0)
  list(data = d, Y = Y, T0 = T0, TT = TT, J = J)
}

## difference in differences rule, written out
stud_did <- function(q, bp, bt) mean(q) + mean(bt - colMeans(bp))

## root mean squared in-sample residual of unit y against training columns bp
stud_sigma <- function(y, bp) {
  fit <- vapply(seq_len(nrow(bp)), function(s) stud_did(y, bp, bp[s, ]),
                numeric(1))
  sqrt(mean((y - fit)^2))
}

## studentized full conformal p-value of candidate c0 at target date t, by
## refitting every control with the treated outcome completed
stud_p_direct <- function(P, t, c0) {
  Ypre <- P$Y[seq_len(P$T0), ]
  ytg <- P$Y[t, ]
  N <- ncol(Ypre)
  sg <- vapply(seq_len(N), function(j) stud_sigma(Ypre[, j], Ypre[, -j, drop = FALSE]),
               numeric(1))
  yc <- ytg
  yc[1L] <- ytg[1L] - c0
  s <- numeric(N)
  s[1L] <- abs(ytg[1L] - stud_did(Ypre[, 1L], Ypre[, -1L, drop = FALSE],
                                  ytg[-1L]) - c0) / sg[1L]
  for (j in 2:N) {
    s[j] <- abs(ytg[j] - stud_did(Ypre[, j], Ypre[, -j, drop = FALSE],
                                  yc[-j])) / sg[j]
  }
  mean(s >= s[1L])
}

test_that("studentized full conformal sets match direct refits (DID rule)", {
  P <- stud_panel()
  out <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"),
                      rule = "did", alpha = 0.2, score = "studentized")
  expect_equal(out$score, "studentized")
  expect_equal(nrow(out$unit), 3L)
  for (k in seq_len(3L)) {
    t <- P$T0 + k
    lo <- out$unit$CI.lower[k]
    hi <- out$unit$CI.upper[k]
    expect_equal(out$unit$p.value[k], stud_p_direct(P, t, 0), tolerance = 1e-12)
    ## the set is {c: p(c) > alpha}: inside just in, outside just out
    h <- 1e-6 * max(abs(c(lo, hi)), 1)
    if (is.finite(lo)) {
      expect_gt(stud_p_direct(P, t, lo + h), 0.2)
      expect_lte(stud_p_direct(P, t, lo - h), 0.2)
    }
    if (is.finite(hi)) {
      expect_gt(stud_p_direct(P, t, hi - h), 0.2)
      expect_lte(stud_p_direct(P, t, hi + h), 0.2)
    }
    mid <- if (is.finite(lo) && is.finite(hi)) (lo + hi) / 2 else out$unit$eff[k]
    if (!out$unit$disconnected[k]) expect_gt(stud_p_direct(P, t, mid), 0.2)
  }
})

test_that("out$sigma holds the hand-computed root mean squared residuals", {
  P <- stud_panel()
  out <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"),
                      rule = "did", score = "studentized")
  expect_named(out$sigma, c("treated", "time", "unit", "sigma"))
  expect_equal(nrow(out$sigma), 3L * (P$J + 1L))
  Ypre <- P$Y[seq_len(P$T0), ]
  s1 <- out$sigma[out$sigma$time == P$T0 + 1L, ]
  expect_equal(s1$unit, 1:(P$J + 1L))
  for (j in seq_len(P$J + 1L)) {
    expect_equal(s1$sigma[j], stud_sigma(Ypre[, j], Ypre[, -j, drop = FALSE]),
                 tolerance = 1e-12)
  }
  ## the weights use history dates only: the same at every date
  s2 <- out$sigma[out$sigma$time == P$T0 + 2L, ]
  expect_equal(s1$sigma, s2$sigma)
  ## absolute: no sigma
  oa <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"),
                     rule = "did")
  expect_null(oa$sigma)
  expect_equal(oa$score, "absolute")
})

test_that("the default score is the absolute score, unchanged", {
  P <- stud_panel()
  a <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"), rule = "did")
  b <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"), rule = "did",
                    score = "absolute")
  expect_identical(a$unit, b$unit)
  expect_identical(a$avg, b$avg)
  ## heterogeneous residuals make the two scores differ
  s <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"), rule = "did",
                    score = "studentized")
  expect_false(isTRUE(all.equal(a$unit$CI.lower, s$unit$CI.lower)))
})

test_that("the studentized Jackknife+ radius is E_j * s0_j / sj_j", {
  P <- stud_panel()
  out <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"),
                      rule = "did", procedure = "jackknife+", alpha = 0.2,
                      score = "studentized")
  t <- P$T0 + 1L
  Ypre <- P$Y[seq_len(P$T0), ]
  ytg <- P$Y[t, ]
  J <- P$J
  M0 <- Rp <- numeric(J)
  for (k in seq_len(J)) {
    j <- k + 1L
    bg <- setdiff(2:(J + 1L), j)
    bp <- Ypre[, bg, drop = FALSE]
    M0[k] <- stud_did(Ypre[, 1L], bp, ytg[bg])
    E <- abs(ytg[j] - stud_did(Ypre[, j], bp, ytg[bg]))
    Rp[k] <- E * stud_sigma(Ypre[, 1L], bp) / stud_sigma(Ypre[, j], bp)
  }
  ## the interval of the paper: a_(q) and b_(N - q) of the fold endpoints
  N <- J + 1L
  qq <- 0
  while ((2 * (qq + 1) - 1) / N <= 0.2) qq <- qq + 1
  a <- sort(M0 - Rp)
  b <- sort(M0 + Rp)
  expect_equal(out$unit$CI.lower[1], ytg[1L] - b[N - qq], tolerance = 1e-12)
  expect_equal(out$unit$CI.upper[1], ytg[1L] - a[qq], tolerance = 1e-12)
  expect_equal(out$unit$eff[1], ytg[1L] - stud_did(Ypre[, 1L], Ypre[, -1L],
                                                    ytg[-1L]), tolerance = 1e-12)
  ## `$sigma` would partially match `$sigma.treated`: use exact indexing
  expect_null(out[["sigma"]])
})

test_that("a user rule that is DID gives the DID result with the same weights", {
  P <- stud_panel()
  a <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"), rule = "did",
                    score = "studentized")
  b <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"),
                    rule = stud_did, score = "studentized")
  expect_equal(b$unit$CI.lower, a$unit$CI.lower, tolerance = 1e-9)
  expect_equal(b$unit$CI.upper, a$unit$CI.upper, tolerance = 1e-9)
  expect_equal(b$sigma$sigma, a$sigma$sigma, tolerance = 1e-12)
  aj <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"), rule = "did",
                     procedure = "jackknife+", score = "studentized")
  bj <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"),
                     rule = stud_did, procedure = "jackknife+",
                     score = "studentized")
  expect_equal(bj$unit$CI.lower, aj$unit$CI.lower, tolerance = 1e-9)
  expect_equal(bj$unit$CI.upper, aj$unit$CI.upper, tolerance = 1e-9)
})

test_that("a rescaled outcome rescales the studentized sets", {
  P <- stud_panel()
  d10 <- P$data
  d10$Y <- 10 * d10$Y
  a <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"), rule = "did",
                    score = "studentized")
  b <- conformal_sc(Y ~ D, data = d10, index = c("id", "time"), rule = "did",
                    score = "studentized")
  expect_equal(b$unit$CI.lower, 10 * a$unit$CI.lower, tolerance = 1e-9)
  expect_equal(b$unit$CI.upper, 10 * a$unit$CI.upper, tolerance = 1e-9)
  expect_equal(b$sigma$sigma, 10 * a$sigma$sigma, tolerance = 1e-9)
})

test_that("the average over dates has no studentized set", {
  P <- stud_panel()
  out <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"),
                      rule = "did", score = "studentized")
  expect_equal(nrow(out$avg), 1L)
  expect_true(all(is.na(out$avg[, c("eff", "CI.lower", "CI.upper", "p.value",
                                    "N", "level")])))
  expect_match(out$avg$note, "not defined for the studentized score")
  txt <- capture.output(print(out))
  expect_true(any(grepl("^Score: studentized", txt)))
  expect_true(any(grepl("not computed with the studentized score", txt)))
  oa <- capture.output(print(conformal_sc(Y ~ D, data = P$data,
                                          index = c("id", "time"),
                                          rule = "did")))
  expect_false(any(grepl("Score:", oa)))
})

test_that("unsupported combinations stop with a clear message", {
  P <- stud_panel()
  ix <- c("id", "time")
  expect_error(conformal_sc(Y ~ D, data = P$data, index = ix, rule = "did",
                            score = "t"), "absolute")
  ## an exact pre-treatment fit has no finite weight: one pre-treatment date
  d1 <- P$data[P$data$time >= P$T0, ]
  expect_error(conformal_sc(Y ~ D, data = d1, index = ix, rule = "did",
                            score = "studentized"),
               "root mean squared residual")
  expect_error(conformal_sc(Y ~ D, data = d1, index = ix, rule = "did",
                            score = "studentized", procedure = "jackknife+"),
               "score = \"absolute\"")
  ## the absolute score still runs there
  expect_s3_class(conformal_sc(Y ~ D, data = d1, index = ix, rule = "did"),
                  "conformal_sc")
})

test_that("a failed history fit gives a failed set or an infinite interval", {
  P <- stud_panel()
  bad <- function(q, bp, bt) {
    if (identical(unname(bt), unname(bp[1L, ]))) NA_real_ else stud_did(q, bp, bt)
  }
  f <- suppressMessages(conformal_sc(Y ~ D, data = P$data,
                                     index = c("id", "time"), rule = bad,
                                     score = "studentized"))
  expect_true(all(f$unit$method == "failed"))
  expect_true(all(is.na(f$unit$CI.lower)))
  j <- suppressMessages(conformal_sc(Y ~ D, data = P$data,
                                     index = c("id", "time"), rule = bad,
                                     procedure = "jackknife+",
                                     score = "studentized"))
  expect_true(all(is.infinite(j$unit$CI.lower) & j$unit$CI.lower < 0))
  expect_true(all(is.infinite(j$unit$CI.upper) & j$unit$CI.upper > 0))
})

test_that("fect(method = \"fe\") accepts conformal.score and refuses the rest", {
  skip_on_cran()
  P <- stud_panel(J = 9)
  d <- P$data
  fit <- function(...) {
    suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"),
                          se = TRUE, vartype = "conformal", CV = FALSE,
                          alpha = 0.2, parallel = FALSE, ...))
  }
  f <- suppressWarnings(fit(method = "fe", conformal.score = "studentized"))
  g <- conformal_sc(Y ~ D, data = d, index = c("id", "time"), rule = "did",
                    alpha = 0.2, score = "studentized")
  expect_equal(f$conformal$score, "studentized")
  expect_equal(f$conformal$unit$CI.lower, g$unit$CI.lower, tolerance = 1e-8)
  expect_equal(f$conformal$unit$CI.upper, g$unit$CI.upper, tolerance = 1e-8)
  expect_true(all(is.na(f$est.avg[, c("CI.lower", "CI.upper", "p.value")])))
  expect_true(is.finite(f$est.avg[1, "ATT.avg"]))
  expect_true(any(grepl("Score:         studentized",
                        capture.output(print(f)))))
  fa <- suppressWarnings(fit(method = "fe"))
  expect_equal(fa$conformal$score, "absolute")
  expect_false(any(grepl("Score:", capture.output(print(fa)))))
  expect_error(fit(method = "fe", conformal.score = "z"), "studentized")
  expect_error(fit(method = "ife", r = 1, conformal.score = "studentized"),
               "method = \"fe\"")
  expect_error(fit(method = "gsynth", r = 1,
                   conformal.score = "studentized"), "conformal.score")
  expect_warning(
    suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"),
                          method = "fe", se = FALSE, CV = FALSE,
                          conformal.score = "studentized")),
    "conformal.score"
  )
})

## ---- pre-treatment diagnostics with the studentized score ---------------------

## studentized full conformal p-value of candidate c0 at target date t when the
## history is the rows `hist` (a placebo date: t is a pre-treatment date)
stud_p_hist <- function(P, hist, t, c0) {
  Ypre <- P$Y[hist, , drop = FALSE]
  ytg <- P$Y[t, ]
  N <- ncol(Ypre)
  sg <- vapply(seq_len(N), function(j) stud_sigma(Ypre[, j], Ypre[, -j, drop = FALSE]),
               numeric(1))
  yc <- ytg
  yc[1L] <- ytg[1L] - c0
  s <- numeric(N)
  s[1L] <- abs(ytg[1L] - stud_did(Ypre[, 1L], Ypre[, -1L, drop = FALSE],
                                  ytg[-1L]) - c0) / sg[1L]
  for (j in 2:N) {
    s[j] <- abs(ytg[j] - stud_did(Ypre[, j], Ypre[, -j, drop = FALSE],
                                  yc[-j])) / sg[j]
  }
  mean(s >= s[1L])
}

test_that("studentized placebo dates match direct refits on the shortened history", {
  P <- stud_panel()
  out <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"),
                      rule = "did", alpha = 0.2, score = "studentized", diag = 3)
  expect_equal(nrow(out$diag), 3L)
  expect_equal(out$diag$time, (P$T0 - 2L):P$T0)
  expect_equal(out$diag$event.time, -2:0)
  expect_equal(out$diag$n.history, (P$T0 - 3L):(P$T0 - 1L))
  for (k in seq_len(3L)) {
    s <- out$diag$time[k]
    hist <- seq_len(s - 1L)
    expect_equal(out$diag$p.value[k], stud_p_hist(P, hist, s, 0),
                 tolerance = 1e-12)
    ## the reported pseudo-effect is the held-out prediction error
    expect_equal(out$diag$eff[k],
                 P$Y[s, 1L] - stud_did(P$Y[hist, 1L], P$Y[hist, -1L], P$Y[s, -1L]),
                 tolerance = 1e-12)
    lo <- out$diag$CI.lower[k]
    hi <- out$diag$CI.upper[k]
    h <- 1e-6 * max(abs(c(lo, hi)), 1)
    expect_gt(stud_p_hist(P, hist, s, lo + h), 0.2)
    expect_lte(stud_p_hist(P, hist, s, lo - h), 0.2)
    expect_gt(stud_p_hist(P, hist, s, hi - h), 0.2)
    expect_lte(stud_p_hist(P, hist, s, hi + h), 0.2)
    ## each unit's sigma-hat is over the history before that date only
    sd <- out$sigma[out$sigma$time == s, ]
    expect_equal(nrow(sd), P$J + 1L)
    Ypre <- P$Y[hist, ]
    for (j in seq_len(P$J + 1L)) {
      expect_equal(sd$sigma[j], stud_sigma(Ypre[, j], Ypre[, -j, drop = FALSE]),
                   tolerance = 1e-12)
    }
  }
  ## the post-treatment rows are what they are without diagnostics
  o0 <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"),
                     rule = "did", alpha = 0.2, score = "studentized")
  expect_identical(out$unit, o0$unit)
  expect_true(any(grepl("Pre-treatment placebo dates",
                        capture.output(print(out)))))
})

test_that("studentized Jackknife+ placebo dates use the history before the date", {
  P <- stud_panel()
  out <- conformal_sc(Y ~ D, data = P$data, index = c("id", "time"),
                      rule = "did", procedure = "jackknife+", alpha = 0.2,
                      score = "studentized", diag = 2)
  J <- P$J
  N <- J + 1L
  qq <- 0
  while ((2 * (qq + 1) - 1) / N <= 0.2) qq <- qq + 1
  for (k in 1:2) {
    s <- out$diag$time[k]
    Ypre <- P$Y[seq_len(s - 1L), ]
    ytg <- P$Y[s, ]
    M0 <- Rp <- numeric(J)
    for (m in seq_len(J)) {
      j <- m + 1L
      bg <- setdiff(2:N, j)
      bp <- Ypre[, bg, drop = FALSE]
      M0[m] <- stud_did(Ypre[, 1L], bp, ytg[bg])
      E <- abs(ytg[j] - stud_did(Ypre[, j], bp, ytg[bg]))
      Rp[m] <- E * stud_sigma(Ypre[, 1L], bp) / stud_sigma(Ypre[, j], bp)
    }
    expect_equal(out$diag$CI.lower[k], ytg[1L] - sort(M0 + Rp)[N - qq],
                 tolerance = 1e-12)
    expect_equal(out$diag$CI.upper[k], ytg[1L] - sort(M0 - Rp)[qq],
                 tolerance = 1e-12)
  }
})

test_that("an exact fit on a shortened history is an error naming the date", {
  P <- stud_panel(T0 = 3L)
  ix <- c("id", "time")
  ## the first placebo date (time 2) has a single history date: DID fits it exactly
  e <- tryCatch(conformal_sc(Y ~ D, data = P$data, index = ix, rule = "did",
                             score = "studentized", diag = 2),
                error = function(e) conditionMessage(e))
  expect_match(e, "at diagnostic date 2", fixed = TRUE)
  expect_match(e, "use a smaller diag", fixed = TRUE)
  expect_match(e, "root mean squared residual", fixed = TRUE)
  ej <- tryCatch(conformal_sc(Y ~ D, data = P$data, index = ix, rule = "did",
                              score = "studentized", diag = 2,
                              procedure = "jackknife+"),
                 error = function(e) conditionMessage(e))
  expect_match(ej, "at diagnostic date 2", fixed = TRUE)
  ## a smaller diag runs, and so does the absolute score
  expect_s3_class(conformal_sc(Y ~ D, data = P$data, index = ix, rule = "did",
                               score = "studentized", diag = 1),
                  "conformal_sc")
  expect_s3_class(conformal_sc(Y ~ D, data = P$data, index = ix, rule = "did",
                               diag = 2), "conformal_sc")
})

test_that("sigma.treated holds the treated unit's divisor for both procedures", {
  P <- stud_panel()
  ix <- c("id", "time")
  fc <- conformal_sc(Y ~ D, data = P$data, index = ix, rule = "did",
                     score = "studentized", diag = 2)
  jp <- conformal_sc(Y ~ D, data = P$data, index = ix, rule = "did",
                     procedure = "jackknife+", score = "studentized", diag = 2)
  for (o in list(fc, jp)) {
    st <- o$sigma.treated
    expect_named(st, c("treated", "time", "sigma"))
    expect_equal(nrow(st), nrow(o$unit) + nrow(o$diag))
    expect_equal(st$time, c(o$unit$time, o$diag$time))
    expect_identical(class(st$time), class(o$unit$time))
    expect_identical(class(st$treated), class(o$unit$unit))
  }
  expect_equal(jp$sigma.treated$sigma, fc$sigma.treated$sigma, tolerance = 1e-12)
  ## full conformal: the treated unit's row of $sigma
  s1 <- fc$sigma[fc$sigma$unit == 1L, ]
  expect_equal(s1$sigma, fc$sigma.treated$sigma, tolerance = 1e-12)
  expect_identical(class(fc$sigma$time), class(fc$unit$time))
  ## by hand at a post date and at a placebo date
  expect_equal(fc$sigma.treated$sigma[1L],
               stud_sigma(P$Y[seq_len(P$T0), 1L], P$Y[seq_len(P$T0), -1L]),
               tolerance = 1e-12)
  s <- fc$diag$time[1L]
  expect_equal(fc$sigma.treated$sigma[nrow(fc$unit) + 1L],
               stud_sigma(P$Y[seq_len(s - 1L), 1L], P$Y[seq_len(s - 1L), -1L]),
               tolerance = 1e-12)
  ## absolute score: nothing stored
  oa <- conformal_sc(Y ~ D, data = P$data, index = ix, rule = "did", diag = 2)
  expect_null(oa[["sigma.treated"]])
  expect_null(oa[["sigma"]])
  ## the time column keeps its type for other kinds of time variable
  d <- P$data
  d$time <- as.Date("2000-01-01") + d$time
  o <- conformal_sc(Y ~ D, data = d, index = ix, rule = "did",
                    score = "studentized", diag = 1)
  expect_s3_class(o$sigma.treated$time, "Date")
  expect_s3_class(o$sigma$time, "Date")
  expect_s3_class(o$unit$time, "Date")
})

test_that("fect(method = \"fe\") runs studentized diagnostics and stores sigma.treated", {
  skip_on_cran()
  P <- stud_panel(J = 9)
  d <- P$data
  f <- suppressWarnings(suppressMessages(
    fect(Y ~ D, data = d, index = c("id", "time"), se = TRUE,
         vartype = "conformal", CV = FALSE, alpha = 0.2, parallel = FALSE,
         method = "fe", conformal.score = "studentized", conformal.diag = 2)))
  g <- conformal_sc(Y ~ D, data = d, index = c("id", "time"), rule = "did",
                    alpha = 0.2, score = "studentized", diag = 2)
  expect_equal(f$conformal$diag$unit$CI.lower, g$diag$CI.lower, tolerance = 1e-8)
  expect_equal(f$conformal$diag$unit$CI.upper, g$diag$CI.upper, tolerance = 1e-8)
  expect_equal(f$conformal$diag$unit$p.value, g$diag$p.value, tolerance = 1e-10)
  expect_equal(f$conformal$sigma.treated$sigma, g$sigma.treated$sigma,
               tolerance = 1e-10)
  expect_identical(class(f$conformal$sigma.treated$time),
                   class(f$conformal$unit$time))
  expect_equal(f$conformal$sigma.treated$time,
               c(f$conformal$unit$time, f$conformal$diag$unit$time))
  fa <- suppressWarnings(suppressMessages(
    fect(Y ~ D, data = d, index = c("id", "time"), se = TRUE,
         vartype = "conformal", CV = FALSE, alpha = 0.2, parallel = FALSE,
         method = "fe", conformal.diag = 2)))
  expect_null(fa$conformal[["sigma.treated"]])
})


test_that("fect(method = \"fe\"): the exact-fit error names the user's time label", {
  P <- stud_panel(T0 = 3L)
  d <- P$data
  d$time <- d$time + 1999L
  e <- tryCatch(suppressWarnings(suppressMessages(
    fect(Y ~ D, data = d, index = c("id", "time"), se = TRUE,
         vartype = "conformal", CV = FALSE, parallel = FALSE, method = "fe",
         force = "two-way", conformal.score = "studentized",
         conformal.diag = 2))), error = function(e) conditionMessage(e))
  expect_match(e, "at diagnostic date 2001", fixed = TRUE)
  expect_match(e, "use a smaller diag", fixed = TRUE)
})

## ---- the event-study plot on the studentized scale ------------------------------

stud_fits <- function(diag = 0L, rule = "did") {
  P <- stud_panel()
  ix <- c("id", "time")
  list(P = P,
       fc = conformal_sc(Y ~ D, data = P$data, index = ix, rule = rule,
                         alpha = 0.2, score = "studentized", diag = diag),
       jp = conformal_sc(Y ~ D, data = P$data, index = ix, rule = rule,
                         alpha = 0.2, score = "studentized", diag = diag,
                         procedure = "jackknife+"),
       fa = conformal_sc(Y ~ D, data = P$data, index = ix, rule = rule,
                         alpha = 0.2, diag = diag),
       ja = conformal_sc(Y ~ D, data = P$data, index = ix, rule = rule,
                         alpha = 0.2, diag = diag, procedure = "jackknife+"))
}

test_that("plot(y.axis = \"outcome\") is the plot of dev", {
  F <- stud_fits(diag = 2L)
  for (o in list(F$fa, F$fc)) {
    p <- plot(o)
    expect_s3_class(p, "ggplot")
    expect_identical(p$data, rbind(o$unit, o$diag))
    expect_identical(p$labels$y, "Effect")
    expect_identical(p$labels$x, "Periods since treatment")
    expect_null(p$labels$caption)
    expect_null(p$labels$title)
    expect_length(p$layers, 3L)
    expect_identical(plot(o, y.axis = "outcome")$data, p$data)
  }
  p2 <- plot(F$fa, xlab = "X", ylab = "Y", main = "M")
  expect_identical(p2$labels$x, "X")
  expect_identical(p2$labels$y, "Y")
  expect_identical(p2$labels$title, "M")
  ## `...` is ignored in the default plot, as in dev
  expect_no_warning(plot(F$fa, foo = 1))
})

test_that("plot(y.axis = \"t\") divides by the treated unit's divisor", {
  for (pr in c("fc", "jp")) {
    F <- stud_fits(diag = 0L)
    o <- F[[pr]]
    p <- plot(o, y.axis = "t")
    d <- p$data
    expect_named(d, c("unit", "time", "event.time", "eff", "CI.lower",
                      "CI.upper", "sigma"))
    expect_equal(d$sigma, o$sigma.treated$sigma)
    expect_equal(d$eff * d$sigma, o$unit$eff, tolerance = 1e-12)
    expect_equal(d$CI.lower * d$sigma, o$unit$CI.lower, tolerance = 1e-12)
    expect_equal(d$CI.upper * d$sigma, o$unit$CI.upper, tolerance = 1e-12)
    ## the question "does the set contain 0" is unchanged
    expect_identical(d$CI.lower <= 0 & d$CI.upper >= 0,
                     o$unit$CI.lower <= 0 & o$unit$CI.upper >= 0)
    expect_identical(p$labels$y, "Effect / pre-treatment residual size")
    expect_lt(nchar(p$labels$y), 45L)
    cp <- gsub("\n", " ", p$labels$caption, fixed = TRUE)
    expect_match(cp, "divided by the treated unit's pre-treatment root mean squared residual.",
                 fixed = TRUE)
    expect_false(grepl("t-test|t-stat", cp))
    expect_lte(max(nchar(strsplit(p$labels$caption, "\n", fixed = TRUE)[[1]])), 90L)
    expect_equal(grepl("he divisor uses all controls", cp), pr == "jp")
    expect_false(grepl("+", p$labels$y, fixed = TRUE))
    expect_length(p$layers, 3L)
  }
  expect_identical(plot(F$fc, y.axis = "t", ylab = "Z", main = "M")$labels$y, "Z")
})

test_that("plot(y.axis = \"t\") draws the placebo dates with their own divisor", {
  F <- stud_fits(diag = 2L)
  p <- plot(F$fc, y.axis = "t")
  d <- p$data
  expect_equal(nrow(d), nrow(F$fc$unit) + 2L)
  expect_equal(d$time, c(F$fc$unit$time, F$fc$diag$time))
  expect_equal(d$eff * d$sigma, c(F$fc$unit$eff, F$fc$diag$eff), tolerance = 1e-12)
  expect_equal(d$CI.upper * d$sigma, c(F$fc$unit$CI.upper, F$fc$diag$CI.upper),
               tolerance = 1e-12)
  expect_true(all(d$event.time[d$time %in% F$fc$diag$time] <= 0))
  expect_equal(length(unique(d$sigma)), 3L)
  expect_match(gsub("\n", " ", p$labels$caption, fixed = TRUE),
               "Dates before the line are placebo dates", fixed = TRUE)
})

test_that("plot(y.axis = \"t\") refuses an absolute-score fit and an old object", {
  F <- stud_fits()
  expect_error(plot(F$fa, y.axis = "t"),
               "needs a fit with the studentized score", fixed = TRUE)
  expect_error(plot(F$fa, y.axis = "t"), "score = \"studentized\"", fixed = TRUE)
  old <- F$fc
  old$sigma.treated <- NULL
  old$sigma <- NULL
  expect_error(plot(old, y.axis = "t"), "refit with the current version")
  expect_error(plot(F$fc, y.axis = "z"), "outcome")
  expect_warning(plot(F$fc, y.axis = "t", foo = 1),
                 'foo: not used with y.axis = "t"', fixed = TRUE)
})

test_that("plot(fc, jp) draws both procedures in either order", {
  F <- stud_fits(diag = 2L)
  for (ax in c("outcome", "t")) {
    p1 <- plot(F$fc, F$jp, y.axis = ax)
    p2 <- plot(F$jp, F$fc, y.axis = ax)
    expect_identical(p1$data, p2$data)
    d <- p1$data
    expect_identical(levels(d$procedure), c("Full conformal", "Jackknife+"))
    expect_equal(sum(d$procedure == "Full conformal"), nrow(F$fc$unit) + 2L)
    expect_equal(sum(d$procedure == "Jackknife+"), nrow(F$jp$unit) + 2L)
    expect_identical(p1$theme$legend.position, "bottom")
    fcrows <- d[d$procedure == "Full conformal", ]
    jprows <- d[d$procedure == "Jackknife+", ]
    if (ax == "t") {
      expect_named(d, c("unit", "time", "event.time", "procedure", "eff",
                        "CI.lower", "CI.upper", "sigma"))
      expect_equal(fcrows$sigma, jprows$sigma)
      expect_equal(fcrows$eff * fcrows$sigma, c(F$fc$unit$eff, F$fc$diag$eff),
                   tolerance = 1e-12)
      expect_equal(jprows$CI.upper * jprows$sigma,
                   c(F$jp$unit$CI.upper, F$jp$diag$CI.upper), tolerance = 1e-12)
      expect_match(gsub("\n", " ", p1$labels$caption, fixed = TRUE),
                   "the divisor uses all controls", fixed = TRUE)
    } else {
      expect_named(d, c("unit", "time", "event.time", "procedure", "eff",
                        "CI.lower", "CI.upper"))
      expect_equal(fcrows$CI.lower, c(F$fc$unit$CI.lower, F$fc$diag$CI.lower))
      expect_match(gsub("\n", " ", p1$labels$caption, fixed = TRUE),
                   "Dates before the line", fixed = TRUE)
      expect_identical(p1$labels$y, "Effect")
    }
  }
  ## absolute-score fits pair in outcome units; no caption without diagnostics
  pa <- plot(F$fa, F$ja)
  expect_null(plot(stud_fits()$fa, stud_fits()$ja)$labels$caption)
  expect_identical(levels(pa$data$procedure), c("Full conformal", "Jackknife+"))
  ## the Jackknife+ rows come first so the full conformal set is on top
  expect_identical(as.character(pa$data$procedure[1L]), "Jackknife+")
})

test_that("plot(x, y) names what differs between the two fits", {
  F <- stud_fits(diag = 2L)
  ix <- c("id", "time")
  expect_error(plot(F$fc, F$fc),
               "one full conformal fit and one Jackknife+ fit", fixed = TRUE)
  expect_error(plot(F$fc, 1), "conformal_sc")
  expect_error(plot(F$fc, F$ja), "differ in score")
  jp.alpha <- conformal_sc(Y ~ D, data = F$P$data, index = ix, rule = "did",
                           alpha = 0.1, score = "studentized", diag = 2,
                           procedure = "jackknife+")
  expect_error(plot(F$fc, jp.alpha), "differ in alpha")
  jp.diag <- conformal_sc(Y ~ D, data = F$P$data, index = ix, rule = "did",
                          alpha = 0.2, score = "studentized", diag = 1,
                          procedure = "jackknife+")
  expect_error(plot(F$fc, jp.diag), "differ in diagnostic dates")
  jp.rule <- conformal_sc(Y ~ D, data = F$P$data, index = ix, rule = stud_did,
                          alpha = 0.2, score = "studentized", diag = 2,
                          procedure = "jackknife+")
  expect_error(plot(F$fc, jp.rule), "differ in rule")
  ## the same data under another rule: the message names the rule alone
  e.rule <- tryCatch(plot(F$fc, jp.rule), error = function(e) conditionMessage(e))
  expect_false(grepl("data", sub("plot them together.*", "", e.rule)))
  jp.sc <- conformal_sc(Y ~ D, data = F$P$data, index = ix, rule = "sc",
                        alpha = 0.2, score = "studentized", diag = 2,
                        procedure = "jackknife+")
  e.sc <- tryCatch(plot(F$fc, jp.sc), error = function(e) conditionMessage(e))
  expect_match(e.sc, "differ in rule;", fixed = TRUE)
  d2 <- F$P$data
  ## a control's outcome at one post-treatment date moves the DID estimate there
  d2$Y <- d2$Y + 1 * (d2$id == 3L & d2$time == F$P$T0 + 1L)
  jp.data <- conformal_sc(Y ~ D, data = d2, index = ix, rule = "did",
                          alpha = 0.2, score = "studentized", diag = 2,
                          procedure = "jackknife+")
  expect_error(plot(F$fc, jp.data), "plot them together only when they differ in procedure alone")
  expect_warning(plot(F$fc, F$jp, foo = 1), "foo")
})

test_that("plot.fect handles y.axis = \"t\" for a studentized fe fit and one fit only", {
  skip_on_cran()
  P <- stud_panel(J = 9)
  d <- P$data
  fit <- function(...) suppressWarnings(suppressMessages(
    fect(Y ~ D, data = d, index = c("id", "time"), se = TRUE,
         vartype = "conformal", CV = FALSE, alpha = 0.2, parallel = FALSE, ...)))
  f <- fit(method = "fe", conformal.score = "studentized", conformal.diag = 2)
  g <- conformal_sc(Y ~ D, data = d, index = c("id", "time"), rule = "did",
                    alpha = 0.2, score = "studentized", diag = 2)
  p <- plot(f, y.axis = "t")
  expect_s3_class(p, "ggplot")
  expect_equal(p$data, plot(g, y.axis = "t")$data, tolerance = 1e-8)
  expect_s3_class(plot(f, y.axis = "t", return.data = TRUE), "data.frame")
  expect_identical(plot(f, y.axis = "t", main = "M", xlab = "X", ylab = "Y",
                        xlim = c(-2, 3), ylim = c(-5, 5))$labels$title, "M")
  expect_warning(plot(f, y.axis = "t", show.count = FALSE),
                 'show.count: not used with y.axis = "t"', fixed = TRUE)
  expect_error(plot(f, y.axis = "t", type = "status"), "type = \"gap\"")
  expect_error(plot(f, g), "conformal_sc objects", fixed = TRUE)
  expect_error(plot(f, f), "conformal_sc objects", fixed = TRUE)
  fa <- fit(method = "fe")
  expect_error(plot(fa, y.axis = "t"), "conformal.score = \"studentized\"",
               fixed = TRUE)
  fi <- fit(method = "ife", r = 1)
  expect_error(plot(fi, y.axis = "t"), "method = \"fe\"", fixed = TRUE)
  expect_error(plot(f, y.axis = "z"), "outcome")
})
