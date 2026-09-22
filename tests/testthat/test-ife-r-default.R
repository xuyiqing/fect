## Regression test: `method = "ife"` (or "gsynth") with the signature default
## r = 0 and no cross-validation estimates no factors --- it is the FEct model.
## fect must say so (2026-09-21: the manual's LOO example ran exactly this call
## and presented a two-way FE pre-trend as an IFE result).

test_that("ife with r = 0 and no CV announces that it is FEct", {
  data(sim_base)
  d <- sim_base[sim_base$id %in% unique(sim_base$id)[1:40], ]
  expect_message(
    fect(Y ~ D, data = d, index = c("id", "time"), method = "ife", se = FALSE),
    "no factors are estimated"
  )
  ## CV = FALSE with a vector r uses its smallest element
  expect_message(
    fect(Y ~ D, data = d, index = c("id", "time"), method = "ife",
         r = c(0, 5), CV = FALSE, se = FALSE),
    "no factors are estimated"
  )
})

test_that("the r = 0 message stays quiet for fe, ife with r > 0, and ife with CV", {
  data(sim_base)
  d <- sim_base[sim_base$id %in% unique(sim_base$id)[1:40], ]
  expect_no_message(
    fect(Y ~ D, data = d, index = c("id", "time"), se = FALSE),
    message = "no factors are estimated"
  )
  expect_no_message(
    fect(Y ~ D, data = d, index = c("id", "time"), method = "ife", r = 2, se = FALSE),
    message = "no factors are estimated"
  )
  expect_no_message(
    fect(Y ~ D, data = d, index = c("id", "time"), method = "ife",
         CV = TRUE, r = c(0, 2), k = 3, se = FALSE),
    message = "no factors are estimated"
  )
})
