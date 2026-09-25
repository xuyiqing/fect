## ---------------------------------------------------------------
## fect 2.4.6: fixes ported onto PR #151 (5afa708) from the home-mac run
## 2026-09-24-fix246-correctness: B10-cfe, B11, B12, B13 and M1. Every
## test_that block below fails on 5afa708 and passes after the port.
## Self-contained: the helpers are defined here (prefix .fix246p_).
## ---------------------------------------------------------------

.fix246p_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## One of fect's datasets, loaded into a local environment.
.fix246p_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "fect", envir = e)
  e[[name]]
}

## simgsynth: 50 units x 30 periods; units 101-105 treated from period 21,
## 45 never-treated units. Y ~ D + X1 + X2, two-way FE, no CV, seed 1.
.fix246p_nt_fit <- function(...) {
  d <- .fix246p_data("simgsynth")
  .fix246p_quiet(fect::fect(Y ~ D + X1 + X2, data = d,
                            index = c("id", "time"), force = "two-way",
                            CV = FALSE, parallel = FALSE, seed = 1,
                            nboots = 10, ...))
}


## -- B10-cfe  cfe + never-treated: replicates, loo SEs, dloo ----------------

test_that("B10-cfe: bootstrap and jackknife replicates of cfe + never-treated use the never-treated model", {
  skip_on_cran()
  c0 <- .fix246p_nt_fit(method = "cfe", time.component.from = "nevertreated",
                        r = 2, se = FALSE)
  for (vt in c("bootstrap", "jackknife")) {
    cc <- .fix246p_nt_fit(method = "cfe", time.component.from = "nevertreated",
                          r = 2, se = TRUE, vartype = vt)
    g <- .fix246p_nt_fit(method = "gsynth", r = 2, se = TRUE, vartype = vt)
    ## the point estimate is the se = FALSE fit's, as before
    expect_identical(cc$att.avg, c0$att.avg, info = vt)
    ## cfe without extra fixed effects is gsynth's model; the two solvers
    ## agree to about 5e-5. On 5afa708 the replicates were not-yet-treated
    ## cfe fits (S.E. off by up to 0.19).
    expect_lt(max(abs(cc$est.att[, "S.E."] - g$est.att[, "S.E."]),
                  na.rm = TRUE), 1e-3)
    expect_lt(max(abs(cc$att.avg.boot - g$att.avg.boot)), 1e-3)
  }
})

test_that("B10-cfe: leave-one-period-out SEs of cfe + never-treated use the never-treated model", {
  skip_on_cran()
  lo <- function(...) .fix246p_nt_fit(r = 2, se = TRUE, loo = TRUE, ...)
  g  <- lo(method = "gsynth")
  a  <- lo(method = "ife", time.component.from = "nevertreated")
  cc <- lo(method = "cfe", time.component.from = "nevertreated")
  ## ife + never-treated is relabelled gsynth (PR #151, B8e): identical
  expect_identical(a$pre.est.att, g$pre.est.att)
  ## ATT column: loo refits keep time.component.from (B8c, on 5afa708);
  ## S.E. column: the refits' replicates (0.16 off on 5afa708)
  cols <- c("ATT", "S.E.")
  expect_lt(max(abs(cc$pre.est.att[, cols] - g$pre.est.att[, cols]),
                na.rm = TRUE), 1e-3)
})

test_that("B10-cfe: dloo with never-treated fixed effects stops before the bootstrap", {
  skip_on_cran()
  ## test-dloo.R's panel: 8 periods; cohorts adopt at 4 and 6; 25 never treated
  dat <- withr::with_seed(42, {
    rows <- list()
    id <- 0
    for (g in c(4, 6, Inf)) for (u in seq_len(25)) {
      id <- id + 1
      y <- stats::rnorm(1) + 0.1 * (1:8) + stats::rnorm(8, 0, 0.4)
      d <- as.integer(!is.infinite(g) & (1:8) >= g)
      rows[[length(rows) + 1]] <- data.frame(id = id, time = 1:8,
                                             Y = y + 2 * d, D = d)
    }
    do.call(rbind, rows)
  })
  ## 5afa708: ran the bootstrap, then stopped with a message about `method`
  expect_error(
    .fix246p_quiet(fect::fect(Y ~ D, data = dat, index = c("id", "time"),
                              method = "fe", force = "two-way", dloo = TRUE,
                              se = TRUE, vartype = "bootstrap", nboots = 10,
                              time.component.from = "nevertreated",
                              parallel = FALSE, seed = 1)),
    "time.component.from")
})
