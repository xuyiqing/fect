## v2.3.1 regression: when W is supplied, fit$est.att / fit$est.avg /
## plot(fit) / print(fit) all agree (W-weighted), and the redundant *.W
## slots are stripped from the user-facing fit object.

test_that("weighted fit: plot trajectory == fit$est.att (W-weighted)", {
  skip_on_cran()

  suppressWarnings(try(data("turnout", package = "fect"), silent = TRUE))
  expect_true(exists("turnout"))

  set.seed(123456)
  ub <- turnout[
    -c(which(turnout$abb == "WY")[1:15],
       sample(seq_len(nrow(turnout)), 50, replace = FALSE)),
  ]
  ub$weight <- runif(nrow(ub), 1, 10)

  fit <- suppressMessages(suppressWarnings(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor,
    data = ub, index = c("abb", "year"),
    method = "mc", CV = TRUE, k = 5, se = TRUE,
    nboots = 30, min.T0 = 5, seed = 42,
    W = "weight", na.rm = TRUE, parallel = FALSE
  )))

  ## the plot data layer must match est.att exactly
  p <- plot(fit, show.points = FALSE)
  plot_dat <- p$data
  expect_true(all(c("Period", "ATT", "CI.lower", "CI.upper") %in% colnames(plot_dat)))

  table <- fit$est.att
  table_periods <- as.numeric(rownames(table))

  ## intersect on shared periods (plot may restrict event-time window)
  shared <- intersect(plot_dat$Period, table_periods)
  expect_gt(length(shared), 0L)

  plot_atts  <- plot_dat$ATT[match(shared, plot_dat$Period)]
  table_atts <- table[as.character(shared), "ATT"]
  expect_equal(unname(plot_atts), unname(table_atts), tolerance = 1e-12)
})

test_that("weighted fit: *.W parallel slots are stripped", {
  skip_on_cran()

  suppressWarnings(try(data("turnout", package = "fect"), silent = TRUE))
  set.seed(123456)
  ub <- turnout[
    -c(which(turnout$abb == "WY")[1:15],
       sample(seq_len(nrow(turnout)), 50, replace = FALSE)),
  ]
  ub$weight <- runif(nrow(ub), 1, 10)

  fit <- suppressMessages(suppressWarnings(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor,
    data = ub, index = c("abb", "year"),
    method = "mc", CV = TRUE, k = 5, se = TRUE,
    nboots = 30, min.T0 = 5, seed = 42,
    W = "weight", na.rm = TRUE, parallel = FALSE
  )))

  ## Inferential *.W slots removed
  expect_null(fit$est.att.W)
  expect_null(fit$est.avg.W)
  expect_null(fit$att.W.bound)
  expect_null(fit$att.W.boot)
  expect_null(fit$att.W.vcov)
  expect_null(fit$est.placebo.W)
  expect_null(fit$est.carryover.W)

  ## Per-method *.W slots removed
  expect_null(fit$att.on.W)
  expect_null(fit$att.avg.W)
  expect_null(fit$time.on.W)
  expect_null(fit$count.on.W)
  expect_null(fit$att.on.sum.W)
  expect_null(fit$W.on.sum)
  expect_null(fit$att.placebo.W)
  expect_null(fit$att.carryover.W)
})

test_that("weighted fit: print() aggregate equals est.avg (W-weighted)", {
  skip_on_cran()

  suppressWarnings(try(data("turnout", package = "fect"), silent = TRUE))
  set.seed(123456)
  ub <- turnout[
    -c(which(turnout$abb == "WY")[1:15],
       sample(seq_len(nrow(turnout)), 50, replace = FALSE)),
  ]
  ub$weight <- runif(nrow(ub), 1, 10)

  fit <- suppressMessages(suppressWarnings(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor,
    data = ub, index = c("abb", "year"),
    method = "mc", CV = TRUE, k = 5, se = TRUE,
    nboots = 30, min.T0 = 5, seed = 42,
    W = "weight", na.rm = TRUE, parallel = FALSE
  )))

  out_lines <- capture.output(print(fit))
  ## label change confirms est.avg is now W-weighted
  expect_true(any(grepl("Tr obs sample-weighted \\(W\\)", out_lines)))
  expect_false(any(grepl("Tr obs equally weighted", out_lines, fixed = TRUE)))
})

test_that("W.est alone: aggregation unweighted, label unweighted", {
  skip_on_cran()

  suppressWarnings(try(data("turnout", package = "fect"), silent = TRUE))
  set.seed(123456)
  ub <- turnout[
    -c(which(turnout$abb == "WY")[1:15],
       sample(seq_len(nrow(turnout)), 50, replace = FALSE)),
  ]
  ub$weight <- runif(nrow(ub), 1, 10)

  fit <- suppressMessages(suppressWarnings(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor,
    data = ub, index = c("abb", "year"),
    method = "mc", CV = TRUE, k = 5, se = TRUE,
    nboots = 30, min.T0 = 5, seed = 42,
    W.est = "weight",
    na.rm = TRUE, parallel = FALSE
  )))

  ## Print label is unweighted
  out_lines <- capture.output(print(fit))
  expect_true(any(grepl("Tr obs equally weighted", out_lines, fixed = TRUE)))
  expect_false(any(grepl("sample-weighted", out_lines, fixed = TRUE)))

  ## per-role flags stashed on the fit
  expect_true(isTRUE(fit$W.in.fit))
  expect_false(isTRUE(fit$W.in.agg))

  ## Slots stripped
  expect_null(fit$est.att.W)
  expect_null(fit$att.on.W)
})

test_that("W.agg alone: fit matches W=NULL fit; aggregation is W-weighted", {
  skip_on_cran()

  suppressWarnings(try(data("turnout", package = "fect"), silent = TRUE))
  set.seed(123456)
  ub <- turnout[
    -c(which(turnout$abb == "WY")[1:15],
       sample(seq_len(nrow(turnout)), 50, replace = FALSE)),
  ]
  ub$weight <- runif(nrow(ub), 1, 10)

  fit_agg <- suppressMessages(suppressWarnings(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor,
    data = ub, index = c("abb", "year"),
    method = "mc", CV = TRUE, k = 5, se = TRUE,
    nboots = 30, min.T0 = 5, seed = 42,
    W.agg = "weight",
    na.rm = TRUE, parallel = FALSE
  )))

  fit_null <- suppressMessages(suppressWarnings(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor,
    data = ub, index = c("abb", "year"),
    method = "mc", CV = TRUE, k = 5, se = TRUE,
    nboots = 30, min.T0 = 5, seed = 42,
    na.rm = TRUE, parallel = FALSE
  )))

  ## Same outcome model (both fit unweighted): residuals/eff matrix identical
  expect_equal(fit_agg$eff, fit_null$eff, tolerance = 1e-10)

  ## But aggregation differs (W-weighted vs unweighted)
  expect_false(isTRUE(all.equal(
    unname(fit_agg$est.avg[1, "ATT.avg"]),
    unname(fit_null$est.avg[1, "ATT.avg"]),
    tolerance = 1e-6
  )))

  ## Print label is sample-weighted
  out_lines <- capture.output(print(fit_agg))
  expect_true(any(grepl("Tr obs sample-weighted \\(W\\)", out_lines)))

  expect_false(isTRUE(fit_agg$W.in.fit))
  expect_true(isTRUE(fit_agg$W.in.agg))
})

test_that("Distinct W.est and W.agg columns: errors in v2.3.1", {
  skip_on_cran()

  suppressWarnings(try(data("turnout", package = "fect"), silent = TRUE))
  set.seed(123456)
  ub <- turnout[
    -c(which(turnout$abb == "WY")[1:15],
       sample(seq_len(nrow(turnout)), 50, replace = FALSE)),
  ]
  ub$weight     <- runif(nrow(ub), 1, 10)
  ub$weight_alt <- runif(nrow(ub), 1, 10)

  expect_error(
    suppressMessages(suppressWarnings(fect::fect(
      turnout ~ policy_edr + policy_mail_in + policy_motor,
      data = ub, index = c("abb", "year"),
      method = "mc", se = FALSE, CV = FALSE,
      W.est = "weight", W.agg = "weight_alt",
      na.rm = TRUE, parallel = FALSE
    ))),
    "v2\\.4\\.0"
  )
})


test_that("unweighted fit (W = NULL): unchanged label, no *.W slots either", {
  skip_on_cran()

  suppressWarnings(try(data("turnout", package = "fect"), silent = TRUE))
  set.seed(123456)
  ub <- turnout[
    -c(which(turnout$abb == "WY")[1:15],
       sample(seq_len(nrow(turnout)), 50, replace = FALSE)),
  ]

  fit <- suppressMessages(suppressWarnings(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor,
    data = ub, index = c("abb", "year"),
    method = "mc", CV = TRUE, k = 5, se = TRUE,
    nboots = 30, min.T0 = 5, seed = 42,
    na.rm = TRUE, parallel = FALSE
  )))

  out_lines <- capture.output(print(fit))
  expect_true(any(grepl("Tr obs equally weighted", out_lines, fixed = TRUE)))
  expect_false(any(grepl("sample-weighted", out_lines, fixed = TRUE)))

  expect_null(fit$est.att.W)
  expect_null(fit$att.on.W)
})
