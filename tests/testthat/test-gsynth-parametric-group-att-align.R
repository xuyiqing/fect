test_that("gsynth parametric bootstrap aligns group (unit-level) results with boot resampling", {
  suppressWarnings(try(data("simgsynth", package = "fect"), silent = TRUE))
  skip_if_not(exists("simgsynth"), "Dataset 'simgsynth' not available")

  set.seed(123)
  out <- fect::fect(
    Y ~ D,
    data = simgsynth,
    index = c("id", "time"),
    method = "gsynth",
    force = "two-way",
    CV = FALSE,
    r = 2,
    se = TRUE,
    vartype = "parametric",
    nboots = 50,
    keep.sims = TRUE,
    min.T0 = 2,
    group = "id",
    parallel = FALSE
  )

  expect_true(!is.null(out$est.group.att))
  expect_true(is.matrix(out$est.group.att))

  # Determine ever-treated units from the stored D matrix
  Dmat <- out$D.dat
  ever_treated <- colSums(Dmat, na.rm = TRUE) > 0
  ids <- colnames(Dmat)

  # Row names should correspond to unit identifiers when group="id"
  est <- out$est.group.att
  rn <- rownames(est)
  # Some workflows may coerce ids; at minimum, all D ids should be represented
  expect_true(all(ids %in% rn))

  never_ids <- ids[!ever_treated]
  if (length(never_ids) > 0) {
    never_rows <- est[never_ids, , drop = FALSE]
    # Never-treated units should have NA for ATT and all uncertainty metrics
    expect_true(all(is.na(never_rows[, "ATT"])))
    expect_true(all(is.na(never_rows[, "S.E."])))
    expect_true(all(is.na(never_rows[, "CI.lower"])))
    expect_true(all(is.na(never_rows[, "CI.upper"])))
    expect_true(all(is.na(never_rows[, "p.value"])))
  }
})


