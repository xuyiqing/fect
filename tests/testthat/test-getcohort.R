test_that("get.cohort adds cohort columns and handles options", {
  suppressWarnings(try(data("simdata", package = "fect"), silent = TRUE))
  expect_true(exists("simdata"))
  df <- simdata

  out <- get.cohort(
    data = df,
    D = "D",
    index = c("id", "time"),
    start0 = TRUE,
    drop.always.treat = TRUE
  )

  expect_true(all(c("Cohort", "Time_to_Treatment") %in% names(out)))
  expect_true(all(out$Cohort %in% c("Control", unique(out$Cohort))))
})
