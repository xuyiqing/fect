## Presentation changes must preserve the plotted estimates and intervals.
publication_fixture <- function() {
  set.seed(314)
  d <- expand.grid(time = 1:12, id = 1:30)
  d$D <- as.integer(d$id <= 8 & d$time >= 8)
  d$Y <- 0.2 * d$id + 0.1 * d$time + 2 * d$D + rnorm(nrow(d))
  fit <- suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"),
    method = "fe", force = "two-way", se = FALSE, parallel = FALSE))
  list(fit = fit, moderator = matrix(rep(rep(0:1, 15), each = 12), 12, 30))
}

test_that("publication HTE preserves estimates and separates counts from CIs", {
  f <- publication_fixture()
  before <- plot(f$fit, type = "heterogeneous", covariate = f$moderator,
    return.data = TRUE)
  after <- plot(f$fit, type = "heterogeneous", covariate = f$moderator,
    preset = "publication", return.data = TRUE)
  expect_equal(after$data$discrete, before$data$discrete)
  expect_lt(max(after$data$count_bars$ymax), min(after$data$discrete$lower))
  expect_true(is.null(after$p$labels$title))
  expect_identical(after$p$labels$y, "Estimate")
  expect_equal(after$p$theme$axis.title$size, 18)
  built <- ggplot2::ggplot_build(after$p)
  expect_equal(built$layout$panel_params[[1]]$x.range, c(0.5, 2.5))
  hlines <- vapply(after$p$layers, function(l) inherits(l$geom, "GeomHline"), logical(1))
  expect_equal(sum(hlines), 1)
  expect_equal(built$data[[which(hlines)]]$yintercept, 0)
})

test_that("mean annotation is optional and never changes heterogeneous estimates", {
  f <- publication_fixture()
  a <- plot(f$fit, type = "heterogeneous", covariate = f$moderator,
    preset = "publication", show.count = FALSE, return.data = TRUE)
  b <- plot(f$fit, type = "heterogeneous", covariate = f$moderator,
    preset = "publication", show.count = FALSE, average.label = TRUE,
    show.average = TRUE, return.data = TRUE)
  expect_equal(a$data$discrete, b$data$discrete)
  texts <- ggplot2::ggplot_build(b$p)$data
  labels <- unlist(lapply(texts, function(d) if ("label" %in% names(d)) d$label))
  expect_true(any(grepl("Mean estimate =", labels, fixed = TRUE)))
  expect_equal(sum(vapply(b$p$layers, function(l) inherits(l$geom, "GeomHline"), logical(1))), 2)
})

test_that("publication event studies retain supplied estimates and supported ticks", {
  d <- data.frame(time = c(-1, 0, 1, 2), ATT = c(0.02, 0, -0.1, -0.3),
    CI.lower = c(-0.03, -0.02, -0.2, -0.5),
    CI.upper = c(0.07, 0.02, 0, -0.1))
  p <- esplot(d, Period = "time", show.count = FALSE, preset = "publication")
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$layout$panel_params[[1]]$x$breaks, d$time)
  expect_true(is.null(p$labels$title))
  expect_identical(p$labels$y, "Estimate")
  ci_layers <- vapply(p$layers, function(l) {
    inherits(l$geom, "GeomPointrange") || inherits(l$geom, "GeomLinerange")
  }, logical(1))
  actual <- do.call(rbind, lapply(built$data[ci_layers], function(x) x[, c("x", "y", "ymin", "ymax")]))
  actual <- actual[order(actual$x), ]
  expect_equal(unname(actual$y), d$ATT)
  expect_equal(unname(actual$ymin), d$CI.lower)
  expect_equal(unname(actual$ymax), d$CI.upper)
  expect_equal(p$theme$axis.title$size, 18)
})

test_that("publication options preserve explicit appearance overrides", {
  d <- data.frame(time = -1:1, ATT = c(0, 0, 1),
    CI.lower = c(-0.1, -0.1, 0.7), CI.upper = c(0.1, 0.1, 1.3))
  p <- esplot(d, Period = "time", show.count = FALSE, preset = "publication",
    main = "Requested title", ylab = "Custom label", cex.lab = 16)
  expect_equal(p$labels$title, "Requested title")
  expect_equal(p$labels$y, "Custom label")
  expect_equal(p$theme$axis.title$size, 16)
})
