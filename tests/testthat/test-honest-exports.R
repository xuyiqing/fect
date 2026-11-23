# The honest-sensitivity engine is internal: fect exports what dev exports and nothing else.

test_that("fect exports exactly the dev names and no honest-engine name", {
  # The 18 names exported by dev (commit 45ec576), from `export()` lines of its NAMESPACE.
  dev_exports <- c("att.cumu", "conformal_sc", "did_wrapper", "effect", "esplot", "estimand", "fect",
                   "fect_iden", "fect_mspe", "fect_rank", "fect_sens", "get.cohort", "imputed_outcomes",
                   "interFE", "loading_check", "overlap_check", "plot.fect", "r.cv.rolling")
  expect_setequal(getNamespaceExports("fect"), dev_exports)
  expect_false(any(c("createSensitivityResults", "createSensitivityResults_relativeMagnitudes",
                     "constructOriginalCS", "createSensitivityPlot", "createSensitivityPlot_relativeMagnitudes",
                     "createEventStudyPlot", "basisVector") %in% getNamespaceExports("fect")))
})

test_that("fect_sens() has the arguments and defaults of dev", {
  fm <- formals(fect::fect_sens)
  expect_identical(names(fm), c("fect.out", "post.periods", "l_vec", "Mbarvec", "Mvec", "periodMbarvec",
                                "periodMvec", "parallel", "cores"))
  expect_identical(deparse(fm$post.periods), "NA")
  expect_identical(deparse(fm$l_vec), "NA")
  expect_identical(deparse(fm$Mbarvec), "seq(0, 1, by = 0.1)")
  expect_identical(deparse(fm$Mvec), "seq(0, 0.25, 0.05)")
  expect_identical(deparse(fm$periodMbarvec), "c(0, 0.5)")
  expect_identical(deparse(fm$periodMvec), "c(0, 0.1)")
  expect_identical(deparse(fm$parallel), "FALSE")
  expect_identical(deparse(fm$cores), "NULL")
})

test_that("no Suggested solver package is attached or imported at load time", {
  loadNamespace("fect")
  expect_false("package:CVXR" %in% search())
  solvers <- c("CVXR", "ECOSolveR", "Rglpk", "lpSolveAPI", "TruncatedNormal", "matrixStats", "Matrix",
               "pracma", "purrr", "tibble")
  expect_false(any(solvers %in% setdiff(names(getNamespaceImports("fect")), "base")))
})
