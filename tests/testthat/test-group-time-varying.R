## Cohort (group) effects with a time-varying group indicator.
##
## A time-invariant `group` (one index per unit) reports the cohort ATT and
## cohort-specific dynamic effects. A time-varying `group` (index changes
## within a unit over time, e.g. regime type by country-year) reports only
## the cohort ATT, defined cell-wise as the mean effect over treated
## observations whose group index equals each level (fect 1.1.x behaviour).

manual_group_att <- function(fit) {
    eff <- fit$eff; D <- fit$D.dat; G <- fit$G
    sel <- D == 1 & !is.na(eff)
    as.numeric(tapply(eff[sel], G[sel], mean))
}

test_that("time-varying group: cohort ATT only, cell-level definition", {
    suppressWarnings(try(data("simdata", package = "fect"), silent = TRUE))
    d <- simdata
    ## group switches over time within units
    d$grp <- ifelse(d$time <= 20, "early", "late")
    fit <- suppressMessages(fect(Y ~ D + X1 + X2, data = d,
        index = c("id", "time"), method = "fe", CV = FALSE, se = FALSE,
        parallel = FALSE, group = "grp"))
    expect_true(isTRUE(fit$group.time.varying))
    expect_equal(names(fit$g.level), c("early", "late"))
    expect_equal(names(fit$group.att), c("early", "late"))
    expect_equal(unname(fit$group.att), manual_group_att(fit))
    expect_null(fit$group.output)
    expect_error(plot(fit, show.group = "early"), "varies over time")
    expect_output(print(fit), "Cohort \\(group\\) ATT")
})

test_that("time-varying group: jackknife inference returns est.group.att only", {
    skip_on_cran()
    suppressWarnings(try(data("simdata", package = "fect"), silent = TRUE))
    d <- simdata[simdata$id %in% unique(simdata$id)[1:60], ]
    d$grp <- ifelse(d$time <= 20, 1, 2)
    fit <- suppressMessages(fect(Y ~ D + X1 + X2, data = d,
        index = c("id", "time"), method = "fe", CV = FALSE, se = TRUE,
        vartype = "jackknife", parallel = FALSE, group = "grp"))
    expect_true(is.matrix(fit$est.group.att))
    expect_equal(rownames(fit$est.group.att), c("1", "2"))
    expect_equal(unname(fit$est.group.att[, "ATT"]), unname(fit$group.att))
    expect_true(all(is.finite(fit$est.group.att[, "S.E."])))
    expect_null(fit$est.group.output)
})

test_that("time-invariant group still reports cohort dynamics", {
    suppressWarnings(try(data("simdata", package = "fect"), silent = TRUE))
    d <- simdata
    d$grp <- ifelse(d$id <= 150, "A", "B")
    fit <- suppressMessages(fect(Y ~ D + X1 + X2, data = d,
        index = c("id", "time"), method = "fe", CV = FALSE, se = FALSE,
        parallel = FALSE, group = "grp"))
    expect_false(isTRUE(fit$group.time.varying))
    expect_equal(names(fit$group.output), c("A", "B"))
    expect_true(length(fit$group.output[["A"]]$att.on) > 0)
    expect_equal(unname(fit$group.att), manual_group_att(fit))
})

test_that("time-varying group with missing cells and binary outcome", {
    set.seed(410)
    N <- 20L; TT <- 22L
    d <- expand.grid(time = seq_len(TT), id = seq_len(N))
    d$D <- as.integer(d$id > 15 & d$time >= 15 + (d$id %% 3))
    d$X <- rnorm(nrow(d), sd = 0.3)
    u <- rnorm(N, sd = 0.25); tt <- rnorm(TT, sd = 0.2)
    d$Y <- rbinom(nrow(d), 1, pnorm(d$X + u[d$id] + tt[d$time] + 0.3 * d$D))
    d$grp <- ifelse(d$time <= 17, "g1", "g2")
    ## drop a few cells so that the panel is unbalanced (missing cells -> 0)
    d <- d[-c(3, 40, 200, 430), ]
    fit <- suppressMessages(fect(Y ~ D + X, data = d, index = c("id", "time"),
        binary = TRUE, force = "two-way", QR = TRUE, parallel = FALSE,
        tol = 1e-3, CV = FALSE, r = 0, se = FALSE, group = "grp"))
    expect_true(isTRUE(fit$group.time.varying))
    expect_equal(names(fit$g.level), c("g1", "g2"))
    expect_equal(unname(fit$group.att), manual_group_att(fit))
    expect_null(fit$group.output)
})
