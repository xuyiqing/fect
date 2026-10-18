## fect >= 2.4.8 needs fixest >= 0.14.0 (DESCRIPTION). fixest 0.13.2's
## cpp_to_index() allocated `first_obs` and `table` without PROTECT, so a
## garbage collection during the call freed them and the returned list held
## dangling pointers. R then crashed at random in did_wrapper()'s bootstrap
## (one feols() per replicate). gctorture() forces a collection at every
## allocation, which makes the bug deterministic: with 0.13.2, `first_obs`
## comes back as the same object as `table`.
test_that("fixest's indexer keeps its results under gctorture", {
  skip_on_cran()
  skip_if_not_installed("fixest")
  to_index <- tryCatch(utils::getFromNamespace("cpp_to_index", "fixest"),
                       error = function(e) NULL)
  skip_if(is.null(to_index), "fixest no longer has cpp_to_index()")

  x <- rep(c("b", "a", "c", "d"), times = c(5, 3, 2, 4))
  for (i in 1:5) {
    gctorture(TRUE)
    res <- to_index(x)
    gctorture(FALSE)
    expect_identical(res$index, rep(1:4, times = c(5, 3, 2, 4)))
    expect_identical(res$first_obs, c(1L, 6L, 9L, 11L))
    expect_identical(res$table, c(5L, 3L, 2L, 4L))
  }
})
