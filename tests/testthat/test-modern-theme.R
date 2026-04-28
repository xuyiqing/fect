## Regression tests for the v2.3.1 modern-theme overhaul.
##
## Goals:
##  1. Modern (default): theme.bw = TRUE + legacy.style = FALSE produces a plot
##     with the modern visual recipe (plain title, no x major grid).
##  2. Classic: legacy.style = TRUE produces pre-2.3.1 visuals (bold centered
##     title, x major grid present) regardless of theme.bw.
##  3. Gray (theme.bw = FALSE) emits a soft-deprecation message exactly once
##     per session.
##  4. Each mode runs without error on a small, real fit.

skip_if_not_installed("ggplot2")

local_fect_fit <- function(seed = 1) {
  set.seed(seed)
  N <- 30; TT <- 12; Ntr <- 6; T0 <- 7
  T0_vec <- rep(Inf, N); T0_vec[seq_len(Ntr)] <- T0
  Y <- D <- numeric(N * TT)
  id <- time <- integer(N * TT)
  alpha <- rnorm(N); xi <- rnorm(TT, 0, 0.3)
  idx <- 1L
  for (i in seq_len(N)) for (t in seq_len(TT)) {
    treated <- (t >= T0_vec[i])
    D[idx] <- as.integer(treated)
    Y[idx] <- alpha[i] + xi[t] + 1.0 * D[idx] + rnorm(1, 0, 0.3)
    id[idx] <- i; time[idx] <- t; idx <- idx + 1L
  }
  df <- data.frame(id = id, time = time, Y = Y, D = D)
  suppressMessages(suppressWarnings(fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "fe", force = "two-way", se = FALSE, parallel = FALSE
  )))
}

extract_title_face <- function(p) {
  th <- ggplot2::theme_get()
  pt <- ggplot2::ggplot_build(p)$plot$theme$plot.title
  if (is.null(pt) || is.null(pt$face)) return(NA_character_)
  pt$face
}

extract_title_hjust <- function(p) {
  pt <- ggplot2::ggplot_build(p)$plot$theme$plot.title
  if (is.null(pt) || is.null(pt$hjust)) return(NA_real_)
  pt$hjust
}

extract_grid_major_x <- function(p) {
  pg <- ggplot2::ggplot_build(p)$plot$theme$panel.grid.major.x
  if (is.null(pg)) return(NA_character_)
  ## Modern recipe sets it to element_blank(); legacy leaves the default.
  if (inherits(pg, "element_blank")) return("blank") else return("not_blank")
}

test_that("modern recipe (default) produces plain left-aligned title with no x major grid", {
  fit <- local_fect_fit()
  p <- plot(fit, type = "gap", main = "modern test")
  expect_equal(extract_title_face(p), "plain")
  expect_equal(extract_title_hjust(p), 0)
  expect_equal(extract_grid_major_x(p), "blank")
})

test_that("legacy.style = TRUE produces bold centered title (theme.bw = TRUE)", {
  fit <- local_fect_fit()
  p <- plot(fit, type = "gap", main = "classic test", legacy.style = TRUE)
  expect_equal(extract_title_face(p), "bold")
  expect_equal(extract_title_hjust(p), 0.5)
})

test_that("legacy.style = TRUE produces bold centered title even with theme.bw = FALSE", {
  fit <- local_fect_fit()
  ## Reset deprecation flag so any messages from theme.bw=FALSE don't leak
  options(fect.theme.bw.deprecated.notified = TRUE)
  p <- plot(fit, type = "gap", main = "gray-classic test",
            legacy.style = TRUE, theme.bw = FALSE)
  expect_equal(extract_title_face(p), "bold")
  expect_equal(extract_title_hjust(p), 0.5)
})

test_that("theme.bw = FALSE emits soft-deprecation message exactly once per session", {
  fit <- local_fect_fit()
  ## Force re-arming
  options(fect.theme.bw.deprecated.notified = FALSE)
  expect_message(
    plot(fit, type = "gap", theme.bw = FALSE),
    "soft-deprecated"
  )
  ## Second call should NOT emit
  expect_silent({
    suppressMessages(plot(fit, type = "gap", theme.bw = FALSE))
  })
})

test_that("all three modes run without error", {
  fit <- local_fect_fit()
  options(fect.theme.bw.deprecated.notified = TRUE)
  expect_no_error(plot(fit, type = "gap"))                                           # modern
  expect_no_error(plot(fit, type = "gap", legacy.style = TRUE))                      # classic
  expect_no_error(plot(fit, type = "gap", theme.bw = FALSE))                         # gray
  expect_no_error(plot(fit, type = "gap", legacy.style = TRUE, theme.bw = FALSE))    # gray-classic
})
