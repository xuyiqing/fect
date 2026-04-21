# Regression test for the "Unsupported bootstrap method: fe" crash when CV
# selects r = 0 under method = "gsynth".
#
# Chain of the original bug (fect <= 2.2.0, diagnosed 2026-04-21):
#   1. fect(method = "gsynth", CV = TRUE, r = c(0, ...), se = TRUE) dispatches
#      to fect_nevertreated.
#   2. When CV picks r.cv == 0, fect_nevertreated used to relabel the outgoing
#      `method` slot to "fe".
#   3. fect_boot reads out$method and dispatches per-iteration bootstrap. It
#      has branches for gsynth / ife / mc / cfe but no "fe" branch, so it
#      would stop() with "Unsupported bootstrap method: fe".
#
# The fix preserves method = "gsynth" (or "cfe") regardless of r.cv. The
# fixture below is a pure two-way-FE DGP that deterministically makes CV
# prefer r = 0 over r >= 1, so this test exercises the r.cv == 0 code path.

test_that("gsynth + CV + se survives when r.cv == 0 (regression for dispatcher crash)", {
  set.seed(20260421)
  N <- 20
  TT <- 15
  id <- rep(seq_len(N), each = TT)
  time <- rep(seq_len(TT), times = N)
  alpha <- rnorm(N)
  xi <- rnorm(TT)
  D <- as.integer(id <= 5 & time >= 10)
  Y <- alpha[id] + xi[time] + 2 * D + rnorm(N * TT, sd = 0.5)
  dat <- data.frame(id = id, time = time, Y = Y, D = D)

  # Before the fix this call raised:
  #   "Unsupported bootstrap method: fe"
  out <- fect::fect(
    Y ~ D,
    data = dat,
    index = c("id", "time"),
    method = "gsynth",
    CV = TRUE,
    r = c(0, 3),
    se = TRUE,
    nboots = 5,
    parallel = FALSE
  )

  expect_s3_class(out, "fect")
  expect_equal(out$r.cv, 0)
  # Method must remain "gsynth" after CV; relabel to "fe" was the bug.
  expect_equal(out$method, "gsynth")
  expect_true(is.numeric(out$att.avg))
  expect_true(!is.null(out$est.avg))
})
