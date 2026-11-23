# Stored bounds of the honest-sensitivity engine (tests/testthat/fixtures/honest-sens-ref.rds).

fixture_path <- testthat::test_path("fixtures", "honest-sens-ref.rds")

test_that("the stored fixture exists and has the expected format", {
  expect_true(file.exists(fixture_path))
  fx <- readRDS(fixture_path)
  expect_identical(fx$version, 1L)
  expect_identical(fx$tolerance, 1e-8)
  expect_true(is.list(fx$cases) && length(fx$cases) > 0L)
  expect_false(is.null(names(fx$cases)))
  expect_false(any(names(fx$cases) == ""))
})

test_that("the engine reproduces every stored case", {
  skip_unless_honest_solvers()
  fx <- readRDS(fixture_path)
  for (p in unique(unlist(lapply(fx$cases, `[[`, "needs")))) testthat::skip_if_not_installed(p)
  for (i in seq_along(fx$cases)) {
    cs <- fx$cases[[i]]
    id <- names(fx$cases)[i]
    res <- suppressWarnings(suppressMessages(do.call(get(cs$fn, envir = asNamespace("fect")), cs$args)))
    res <- as.data.frame(res)
    expect_identical(names(res), names(cs$expected), info = id)
    expect_identical(nrow(res), nrow(cs$expected), info = id)
    expect_close(res$lb, cs$expected$lb, fx$tolerance, info = id)
    expect_close(res$ub, cs$expected$ub, fx$tolerance, info = id)
    for (v in setdiff(names(cs$expected), c("lb", "ub"))) {
      if (is.numeric(cs$expected[[v]])) {
        expect_close(res[[v]], cs$expected[[v]], 1e-12, info = paste(id, v))
      } else {
        expect_identical(res[[v]], cs$expected[[v]], info = paste(id, v))
      }
    }
  }
})

test_that("no stored result is all-NA", {
  fx <- readRDS(fixture_path)
  for (i in seq_along(fx$cases)) {
    cs <- fx$cases[[i]]
    expect_true(any(is.finite(c(cs$expected$lb, cs$expected$ub))), info = names(fx$cases)[i])
  }
})
