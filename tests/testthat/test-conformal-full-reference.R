## Full conformal inference and the Jackknife+ against the reference results
## of Liu and Xu (2026): Proposition 99 and German reunification, DID and
## ridge-simplex synthetic control, 5 pre-treatment placebo dates and the
## post-treatment dates. fixtures/scconf-adh.rds holds the two panels
## (Abadie, Diamond and Hainmueller 2010, 2015; West German GDP in thousands)
## and the paper's reference table for these designs.

scconf <- readRDS(test_path("fixtures", "scconf-adh.rds"))
scconf_tr <- list(prop99 = c("California", 1989), germany = c("West Germany", 1990))

scconf_data <- function(nm) {
  d <- scconf[[nm]]
  d$D <- as.integer(d$id == scconf_tr[[nm]][1] & d$time >= as.integer(scconf_tr[[nm]][2]))
  d
}

test_that("conformal_sc reproduces the paper's reference results", {
  skip_on_cran()
  skip_if_not_installed("quadprog")
  ref <- scconf$ref
  for (nm in names(scconf_tr)) {
    d <- scconf_data(nm)
    for (est in c("did_twfe", "synthetic_control")) {
      rule <- if (est == "did_twfe") "did" else "sc"
      for (pr in c("full", "jackknife+")) {
        design <- if (pr == "full") "pointwise_full_conformal" else
          "pointwise_jackknife_plus"
        f <- conformal_sc(Y ~ D, data = d, index = c("id", "time"),
                          rule = rule, procedure = pr, alpha = 0.1, diag = 5)
        mine <- rbind(f$diag, f$unit)
        r <- ref[ref$scenario == nm & ref$estimator == est &
                   ref$design == design, ]
        r <- r[match(as.integer(mine$time), r$target_year), ]
        expect_false(anyNA(r$target_year))
        expect_equal(mine$n.history, r$n_history)
        expect_equal(mine$eff, r$estimate, tolerance = 1e-8)
        expect_equal(mine$CI.lower, r$ci_lower, tolerance = 1e-8)
        expect_equal(mine$CI.upper, r$ci_upper, tolerance = 1e-8)
        pl <- r$period_type == "placebo"
        expect_equal(mine$p.value[pl], r$p_value[pl], tolerance = 1e-12)
      }
    }
  }
})

test_that("the average over post-treatment dates matches the paper (exact)", {
  skip_on_cran()
  skip_if_not_installed("quadprog")
  ## Liu and Xu (2026), Sec. 6.2 and the 2026-11-10 computation:
  ## Prop 99 [-35.0540, -3.9733] (N = 39, rank 3, p = 0.077);
  ## West Germany [-4.7079, 2.1129] (N = 17, p = 0.412); at gamma = 0.05
  ## California [-58.8, 19.8] and West Germany the whole line.
  fp <- conformal_sc(Y ~ D, data = scconf_data("prop99"),
                     index = c("id", "time"), rule = "sc", alpha = 0.1)
  expect_equal(round(c(fp$avg$CI.lower, fp$avg$CI.upper), 4),
               c(-35.0540, -3.9733))
  expect_equal(fp$avg$p.value, 3 / 39)
  fg <- conformal_sc(Y ~ D, data = scconf_data("germany"),
                     index = c("id", "time"), rule = "sc", alpha = 0.1)
  expect_equal(round(c(fg$avg$CI.lower, fg$avg$CI.upper), 4),
               c(-4.7079, 2.1129))
  expect_equal(fg$avg$p.value, 7 / 17)
  fp5 <- conformal_sc(Y ~ D, data = scconf_data("prop99"),
                      index = c("id", "time"), rule = "sc", alpha = 0.05)
  expect_equal(round(c(fp5$avg$CI.lower, fp5$avg$CI.upper), 1),
               c(-58.8, 19.8))
  fg5 <- conformal_sc(Y ~ D, data = scconf_data("germany"),
                      index = c("id", "time"), rule = "sc", alpha = 0.05)
  expect_equal(c(fg5$avg$CI.lower, fg5$avg$CI.upper), c(-Inf, Inf))
})

test_that("fect(method = \"fe\") reproduces the paper's DID full conformal sets", {
  skip_on_cran()
  ref <- scconf$ref
  for (nm in names(scconf_tr)) {
    d <- scconf_data(nm)
    f <- suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"),
                               method = "fe", se = TRUE,
                               vartype = "conformal", alpha = 0.1,
                               conformal.diag = 5, parallel = FALSE))
    u <- f$conformal$unit
    r <- ref[ref$scenario == nm & ref$estimator == "did_twfe" &
               ref$design == "pointwise_full_conformal", ]
    r <- r[match(as.integer(u$time), r$target_year), ]
    expect_equal(u$eff, u$eff.fit, tolerance = 1e-8)
    expect_equal(u$eff, r$estimate, tolerance = 1e-8)
    expect_equal(u$CI.lower, r$ci_lower, tolerance = 1e-8)
    expect_equal(u$CI.upper, r$ci_upper, tolerance = 1e-8)
  }
})
