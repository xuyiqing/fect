# Code adapted from the 'HonestDiD' package, version 0.2.8, by Ashesh Rambachan and Jonathan Roth,
# through 'HonestDiDFEct' 1.0.0 by Ziyi Liu (github.com/lzy318/HonestDiDFEct, commit 29c831d).
# Original file: tests/testthat/helper-path.R
# Copyright (c) 2023 HonestDiD authors. MIT licence; the full notice is in inst/COPYRIGHTS.
# Copied into fect 2.5.1 without change to any numerical line; the edits are listed in the fect repository history.
#
# A fect-style path: placebo periods -2, -1, 0, then post periods 1 to 4, no reference period.
path_s <- c(-2:0, 1:4)
path_b <- stats::setNames(c(0.020, -0.012, 0.015, 0.071, 0.088, 0.064, 0.102), path_s)
path_V <- 0.0004 * 0.5^abs(outer(1:7, 1:7, "-")) + 0.0001
path_l <- rep(0.25, 4)

# The same steps by hand, for HonestDiD's engine.
by_hand <- function(b, V, p, q, l) {
  k <- p + q
  J <- diag(k)[-p, , drop = FALSE]
  J[, p] <- J[, p] - 1
  d <- as.numeric(J %*% unname(b))
  S <- J %*% unname(V) %*% t(J)
  S <- (S + t(S)) / 2
  post <- p:(k - 1L)
  shift <- sum(l * d[post])
  d[post] <- d[post] - shift * l / sum(l^2)
  list(betahat = d, sigma = S, numPrePeriods = p - 1L, numPostPeriods = q, shift = shift)
}

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
