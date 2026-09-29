## ---------------------------------------------------------------
## fect 2.4.7: plot(type = "gap", id = ...) draws the chosen treated units'
## own effects (fect #162), with a band when the fit kept parametric
## bootstrap draws. Before, `id` was ignored and the gap plot showed the
## average effect over all treated units. Each block fails on f39791c and
## passes after the change. Self-contained: helpers prefixed .gi_.
## ---------------------------------------------------------------

## One of fect's datasets, loaded into a local environment.
.gi_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "fect", envir = e)
  e[[name]]
}

## Run a plot call; return the plot and the messages it emitted.
.gi_plot <- function(expr) {
  msgs <- character(0)
  p <- withCallingHandlers(expr, message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  list(p = p, messages = msgs)
}

.gi_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], "")

## The estimates the plot draws (default unconnected style): x, y, ymin and
## ymax of its point-range layers, sorted by x.
.gi_points <- function(p) {
  idx <- which(.gi_geoms(p) == "GeomPointrange")
  d <- do.call(rbind, lapply(idx, function(i)
    ggplot2::layer_data(p, i)[, c("x", "y", "ymin", "ymax")]))
  d <- d[order(d$x), ]
  rownames(d) <- NULL
  d
}

## All drawn estimates, highlighted periods included: x and y of the
## point-range and point layers, sorted by x.
.gi_xy <- function(p) {
  idx <- which(.gi_geoms(p) %in% c("GeomPointrange", "GeomPoint"))
  d <- do.call(rbind, lapply(idx, function(i)
    ggplot2::layer_data(p, i)[, c("x", "y")]))
  d <- d[order(d$x), ]
  rownames(d) <- NULL
  d
}

## The labels of the text layers (count label, test statistics).
.gi_text <- function(p) {
  idx <- which(.gi_geoms(p) == "GeomText")
  unlist(lapply(idx, function(i) ggplot2::layer_data(p, i)$label),
         use.names = FALSE)
}

## By hand from the fit: at each time relative to the treatment onset
## (fit$T.on), the mean of fit$eff over the chosen units' cells, and the
## number of cells.
.gi_by_hand <- function(fit, ids) {
  j <- match(as.character(ids), as.character(fit$id))
  rel <- fit$T.on[, j, drop = FALSE]
  eff <- fit$eff[, j, drop = FALSE]
  ok <- !is.na(rel) & !is.na(eff)
  data.frame(x = sort(unique(rel[ok])),
             y = as.numeric(tapply(eff[ok], rel[ok], mean)),
             n = as.numeric(table(rel[ok])))
}


## -- G1  one treated unit ------------------------------------------------

test_that("G1: the gap plot with one treated id draws that unit's own effects", {
  skip_on_cran()
  simgsynth <- .gi_data("simgsynth")
  fit <- suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE, parallel = FALSE))
  j <- which(fit$id == 101)
  r <- .gi_plot(plot(fit, type = "gap", id = 101))
  pts <- .gi_points(r$p)
  ## unit 101 is observed in all 30 periods: one point per relative time,
  ## its own effect fit$eff[, j]
  expect_equal(pts$x, as.numeric(fit$T.on[, j]))
  expect_equal(pts$y, as.numeric(fit$eff[, j]), tolerance = 1e-12)
  expect_equal(pts$y[pts$x == 1], 0.3381493, tolerance = 1e-6)
  expect_true(all(is.na(pts$ymin)) && all(is.na(pts$ymax)))
  expect_identical(r$p$labels$title, "id = 101")
  ## no count bars for one unit
  expect_false("GeomRect" %in% .gi_geoms(r$p))
  ## a fit without SEs: the same message as the average gap plot
  expect_identical(r$messages, "Uncertainty estimates not available.\n\n")
  ## a character id gives the same plot; a user title wins
  expect_equal(.gi_points(.gi_plot(plot(fit, type = "gap", id = "101"))$p), pts)
  pm <- suppressMessages(plot(fit, type = "gap", id = 101, main = "Unit 101"))
  expect_identical(pm$labels$title, "Unit 101")
  ## the default type is the gap plot
  expect_equal(.gi_points(suppressMessages(plot(fit, id = 101))), pts)
  ## an explicit plot.ci on a fit without SEs does not stop for an id
  ## (the average gap plot stops there: "No uncertainty estimates")
  expect_equal(.gi_points(suppressMessages(
    plot(fit, type = "gap", id = 101, plot.ci = "0.95"))), pts)
  ## id = NULL: the average effect over all treated units, as before
  r0 <- .gi_plot(plot(fit, type = "gap"))
  expect_equal(.gi_points(r0$p)$y, as.numeric(fit$att), tolerance = 1e-12)
  expect_identical(r0$p$labels$title, "Estimated Dynamic Treatment Effects")
  expect_identical(r0$messages, "Uncertainty estimates not available.\n\n")
  expect_true("GeomRect" %in% .gi_geoms(r0$p))
})


## -- G2  several treated units, staggered ------------------------------

test_that("G2: several ids draw the average of their effects at each relative time, with their count bars", {
  skip_on_cran()
  simdata <- .gi_data("simdata")
  fit <- suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
    method = "fe", force = "two-way", se = FALSE, parallel = FALSE))
  ## each switches on once and stays on, in periods 19, 13 and 9
  ids <- c(101, 130, 213)
  ref <- .gi_by_hand(fit, ids)
  r <- .gi_plot(plot(fit, type = "gap", id = ids))
  pts <- .gi_points(r$p)
  ## with 3 units, proportion = 0.3 keeps every relative time that has a cell
  expect_equal(pts$x, ref$x)
  expect_equal(pts$y, ref$y, tolerance = 1e-12)
  expect_equal(range(pts$x), c(-17, 27))
  expect_setequal(ref$n, 1:3)
  expect_true(all(is.na(pts$ymin)))
  expect_identical(r$p$labels$title, "Average over 3 units")
  expect_length(r$messages, 1)
  ## count bars: one per relative time, height proportional to the number
  ## of chosen cells there
  rect <- ggplot2::layer_data(r$p, which(.gi_geoms(r$p) == "GeomRect"))
  rect <- rect[order(rect$xmin), ]
  h <- rect$ymax - rect$ymin
  expect_equal((rect$xmin + rect$xmax) / 2, ref$x)
  expect_equal(h / max(h), ref$n / max(ref$n), tolerance = 1e-12)
  ## show.count = FALSE drops them
  pn <- suppressMessages(plot(fit, type = "gap", id = ids, show.count = FALSE))
  expect_false("GeomRect" %in% .gi_geoms(pn))
  expect_equal(.gi_points(pn)$y, ref$y, tolerance = 1e-12)
})


## -- G3  a bootstrap fit: estimates only ---------------------------------

test_that("G3: with a bootstrap fit, the gap plot for an id draws no interval, statistics or bounds, and says why", {
  skip_on_cran()
  simgsynth <- .gi_data("simgsynth")
  fit <- suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = TRUE, nboots = 20,
    parallel = FALSE, seed = 1))
  j <- which(fit$id == 102)
  ## the average gap plot of this fit has intervals and an F test
  p0 <- suppressMessages(plot(fit, type = "gap", stats = "F.p"))
  expect_true(all(is.finite(.gi_points(p0)$ymin)))
  expect_true(any(grepl("F test p-value", .gi_text(p0), fixed = TRUE)))
  for (args in list(list(), list(plot.ci = "0.95"), list(plot.ci = "0.9"),
                    list(stats = "F.p"), list(bound = "equiv"))) {
    r <- .gi_plot(do.call(plot, c(list(fit, type = "gap", id = 102), args)))
    pts <- .gi_points(r$p)
    expect_equal(pts$y, as.numeric(fit$eff[, j]), tolerance = 1e-12)
    expect_true(all(is.na(pts$ymin)) && all(is.na(pts$ymax)))
    expect_false("GeomLine" %in% .gi_geoms(r$p)) # no bound lines
    expect_true(all(.gi_text(r$p) == ""))        # no test statistics
    expect_length(r$messages, 1)
    expect_match(r$messages, "a band needs parametric bootstrap draws",
                 fixed = TRUE)
    expect_match(r$messages, "keep each unit's own outcomes fixed",
                 fixed = TRUE)
  }
  ## return.test: the fit's tests are about the average, so the unit plot,
  ## which shows none, returns none
  expect_false(is.null(suppressMessages(
    plot(fit, type = "gap", return.test = TRUE))$test.out))
  expect_null(suppressMessages(
    plot(fit, type = "gap", id = 102, return.test = TRUE))$test.out)
  ## connected style: the ribbons carry no interval
  pc <- suppressMessages(plot(fit, type = "gap", id = 102, connected = TRUE))
  rib <- which(.gi_geoms(pc) == "GeomRibbon")
  expect_true(length(rib) > 0)
  for (i in rib) {
    expect_true(all(is.na(ggplot2::layer_data(pc, i)$ymin)))
  }
})


## -- G4  placebo fits ------------------------------------------------------

test_that("G4: placebo fits: the gap plot for an id marks the placebo periods, shows no p-value and needs no SEs", {
  skip_on_cran()
  simgsynth <- .gi_data("simgsynth")
  args <- list(Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
               method = "ife", r = 2, CV = FALSE, parallel = FALSE,
               placeboTest = TRUE, placebo.period = c(-2, 0))
  fit_nose <- suppressMessages(do.call(fect::fect, c(args, list(se = FALSE))))
  ## the average gap plot needs SEs for the placebo test
  expect_error(suppressMessages(plot(fit_nose, type = "gap")),
               "No uncertainty estimates", fixed = TRUE)
  j <- which(fit_nose$id == 101)
  r <- .gi_plot(plot(fit_nose, type = "gap", id = 101))
  xy <- .gi_xy(r$p)
  expect_equal(xy$x, as.numeric(fit_nose$T.on[, j]))
  expect_equal(xy$y, as.numeric(fit_nose$eff[, j]), tolerance = 1e-12)
  ## the placebo periods -2, -1, 0 are drawn as triangles
  tri <- which(.gi_geoms(r$p) == "GeomPoint")
  tri_x <- unlist(lapply(tri, function(i) ggplot2::layer_data(r$p, i)$x),
                  use.names = FALSE)
  tri_shape <- unlist(lapply(tri, function(i) ggplot2::layer_data(r$p, i)$shape),
                      use.names = FALSE)
  expect_equal(sort(tri_x), c(-2, -1, 0))
  expect_true(all(tri_shape == 17))
  expect_length(r$messages, 1)
  ## with SEs: the average plot shows the placebo p-values, the unit plot none
  fit_se <- suppressMessages(do.call(fect::fect,
                                     c(args, list(se = TRUE, nboots = 20, seed = 1))))
  p0 <- suppressMessages(plot(fit_se, type = "gap"))
  expect_true(any(grepl("Placebo test p-value", .gi_text(p0), fixed = TRUE)))
  p1 <- suppressMessages(plot(fit_se, type = "gap", id = 101))
  expect_true(all(.gi_text(p1) == ""))
  expect_equal(.gi_xy(p1)$y, as.numeric(fit_se$eff[, j]), tolerance = 1e-12)
  ## a one-period placebo window without SEs: the unit's effects, no triangles
  fit_one <- suppressMessages(do.call(fect::fect,
    c(args[names(args) != "placebo.period"],
      list(placebo.period = 0, se = FALSE))))
  r1 <- .gi_plot(plot(fit_one, type = "gap", id = 101))
  expect_equal(.gi_xy(r1$p)$y, as.numeric(fit_one$eff[, j]), tolerance = 1e-12)
  expect_false("GeomPoint" %in% .gi_geoms(r1$p))
})


## -- G5  ids that cannot be drawn --------------------------------------

test_that("G5: control ids and unknown ids stop with a message that names them", {
  skip_on_cran()
  simgsynth <- .gi_data("simgsynth")
  fit <- suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE, parallel = FALSE))
  expect_error(plot(fit, type = "gap", id = 106),
               "Unit(s) in \"id\" never treated (control units): 106.", fixed = TRUE)
  expect_error(plot(fit, type = "gap", id = c(101, 106, 107)),
               "never treated (control units): 106, 107.", fixed = TRUE)
  expect_error(plot(fit, id = 106), "never treated (control units): 106.",
               fixed = TRUE)
  expect_error(plot(fit, type = "gap", id = 999),
               "Unit(s) in \"id\" not in the data: 999.", fixed = TRUE)
  expect_error(plot(fit, type = "gap", id = c(101, 999)),
               "Unit(s) in \"id\" not in the data: 999.", fixed = TRUE)
  expect_error(plot(fit, type = "gap", id = "Wisconsin"),
               "not in the data: Wisconsin.", fixed = TRUE)
  ## the counterfactual plot keeps its own check
  expect_error(suppressMessages(plot(fit, type = "counterfactual", id = 106)),
               "is not treated", fixed = TRUE)
})


## -- G6  treatment reversals ---------------------------------------------

test_that("G6: with reversals, a unit's cells at the same relative time are averaged, as in the average gap plot", {
  skip_on_cran()
  simdata <- .gi_data("simdata")
  fit <- suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
    method = "fe", force = "two-way", se = FALSE, parallel = FALSE))
  ## unit 103 switches on in period 15, off in 16 and on again in 18, so
  ## relative times -1, 0 and 1 each occur twice
  j <- which(fit$id == 103)
  ref <- .gi_by_hand(fit, 103)
  expect_equal(ref$n[ref$x %in% c(-1, 0, 1)], c(2, 2, 2))
  r <- .gi_plot(plot(fit, type = "gap", id = 103))
  pts <- .gi_points(r$p)
  expect_equal(pts$x, ref$x)
  expect_equal(pts$y, ref$y, tolerance = 1e-12)
  expect_equal(pts$y[pts$x == 1], mean(fit$eff[c(15, 18), j]), tolerance = 1e-12)
  expect_false("GeomRect" %in% .gi_geoms(r$p))
  ## all treated units: the same estimates and count bars as the average plot
  pall <- suppressMessages(plot(fit, type = "gap", id = fit$id[fit$tr]))
  p0 <- suppressMessages(plot(fit, type = "gap"))
  expect_equal(.gi_points(pall)[, c("x", "y")], .gi_points(p0)[, c("x", "y")],
               tolerance = 1e-12)
  expect_equal(.gi_points(p0)$y, as.numeric(fit$att[fit$time %in% .gi_points(p0)$x]),
               tolerance = 1e-12)
  cols <- c("xmin", "xmax", "ymin", "ymax")
  expect_equal(ggplot2::layer_data(pall, which(.gi_geoms(pall) == "GeomRect"))[, cols],
               ggplot2::layer_data(p0, which(.gi_geoms(p0) == "GeomRect"))[, cols],
               tolerance = 1e-12)
  expect_identical(pall$labels$title, paste0("Average over ", length(fit$tr), " units"))
})


## -- G7  options: start0, xlim, return.data, loo ----------------------------

test_that("G7: start0, xlim and return.data apply to the gap plot for an id; loo stops", {
  skip_on_cran()
  simgsynth <- .gi_data("simgsynth")
  fit <- suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    method = "fe", force = "two-way", se = FALSE, parallel = FALSE))
  j <- which(fit$id == 104)
  rel <- as.numeric(fit$T.on[, j])
  eff <- as.numeric(fit$eff[, j])
  ## start0: the first treated period is 0, the onset line at -0.5
  p <- suppressMessages(plot(fit, type = "gap", id = 104, start0 = TRUE))
  expect_equal(.gi_points(p)$x, rel - 1)
  expect_equal(.gi_points(p)$y, eff, tolerance = 1e-12)
  expect_equal(ggplot2::layer_data(p, which(.gi_geoms(p) == "GeomVline"))$xintercept, -0.5)
  ## xlim keeps the relative times in the window
  p <- suppressMessages(plot(fit, type = "gap", id = 104, xlim = c(-3, 5)))
  expect_equal(.gi_points(p)$x, -3:5)
  expect_equal(.gi_points(p)$y, eff[rel %in% -3:5], tolerance = 1e-12)
  ## return.data returns the plotted estimates
  out <- suppressMessages(plot(fit, type = "gap", id = 104, return.data = TRUE))
  expect_equal(out$data$estimate$Period, rel)
  expect_equal(out$data$estimate$ATT, eff, tolerance = 1e-12)
  ## loo = TRUE: the leave-one-out estimates exist for the average only
  fit_loo <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    method = "fe", force = "two-way", se = TRUE, nboots = 20, loo = TRUE,
    parallel = FALSE, seed = 1)))
  expect_error(plot(fit_loo, type = "gap", id = 104, loo = TRUE),
               "\"loo\" and \"dloo\" can't be used with \"id\" in the gap plot",
               fixed = TRUE)
  expect_equal(.gi_points(suppressMessages(plot(fit_loo, type = "gap", id = 104)))$y,
               as.numeric(fit_loo$eff[, j]), tolerance = 1e-12)
})


## -- G8-G12  the band from parametric bootstrap draws -----------------------

## By hand from the fit: at each relative time, the chosen cells' draws
## (fit$eff.boot keeps the treated units first, in fit$tr order) averaged
## within each replication, and the band with total tail probability `a`:
## the estimate -/+ qnorm(1 - a / 2) times the draws' sd ("normal"), or the
## basic interval of the draws shifted to the estimate ("basic").
.gi_band_by_hand <- function(fit, ids, a = 0.05, rule = "normal") {
  ref <- .gi_by_hand(fit, ids)
  j <- match(as.character(ids), as.character(fit$id))
  k <- match(j, fit$tr)
  rel <- fit$T.on[, j, drop = FALSE]
  eff <- fit$eff[, j, drop = FALSE]
  draws <- t(vapply(ref$x, function(s) {
    cells <- which(!is.na(rel) & !is.na(eff) & rel == s, arr.ind = TRUE)
    apply(fit$eff.boot, 3, function(b) mean(b[cbind(cells[, 1], k[cells[, 2]])]))
  }, numeric(dim(fit$eff.boot)[3])))
  if (rule == "normal") {
    h <- qnorm(1 - a / 2) * apply(draws, 1, sd)
    lo <- ref$y - h
    hi <- ref$y + h
  } else {
    q <- apply(draws - rowMeans(draws) + ref$y, 1, quantile,
               probs = c(a / 2, 1 - a / 2))
    lo <- 2 * ref$y - q[2, ]
    hi <- 2 * ref$y - q[1, ]
  }
  data.frame(x = ref$x, y = ref$y, ymin = lo, ymax = hi)
}

.gi_param_fit <- function(data, ...) {
  suppressWarnings(suppressMessages(fect::fect(
    Y ~ D + X1 + X2, data = data, index = c("id", "time"),
    method = "gsynth", r = 2, CV = FALSE, se = TRUE, vartype = "parametric",
    nboots = 50, keep.sims = TRUE, parallel = FALSE, seed = 1, ...)))
}

test_that("G8: with parametric draws, the gap plot for an id draws a band from the chosen units' draws", {
  skip_on_cran()
  fit <- .gi_param_fit(.gi_data("simgsynth"))
  expect_identical(fit$ci.method, "normal")
  expect_identical(fit$ci.alpha, 0.05)
  ## one unit: its effects and a band from its own draws, with no message
  r <- .gi_plot(plot(fit, type = "gap", id = 102))
  expect_length(r$messages, 0)
  pts <- .gi_points(r$p)
  ref <- .gi_band_by_hand(fit, 102)
  expect_equal(pts, ref, tolerance = 1e-10)
  expect_true(all(is.finite(pts$ymin)) && all(is.finite(pts$ymax)))
  ## one unit's band is wider than the band of the average
  p0 <- suppressMessages(plot(fit, type = "gap"))
  post <- pts$x >= 1
  expect_gt(mean((pts$ymax - pts$ymin)[post]),
            1.5 * mean((.gi_points(p0)$ymax - .gi_points(p0)$ymin)[post]))
  ## still no test statistics or bounds
  for (args in list(list(stats = "F.p"), list(bound = "equiv"))) {
    p <- suppressMessages(do.call(plot, c(list(fit, type = "gap", id = 102), args)))
    expect_true(all(.gi_text(p) == ""))
    expect_false("GeomLine" %in% .gi_geoms(p))
  }
  ## plot.ci = "0.9": the 90% band; "none": no band and no message
  p90 <- suppressMessages(plot(fit, type = "gap", id = 102, plot.ci = "0.9"))
  expect_equal(.gi_points(p90), .gi_band_by_hand(fit, 102, a = 0.1), tolerance = 1e-10)
  rn <- .gi_plot(plot(fit, type = "gap", id = 102, plot.ci = "none"))
  expect_true(all(is.na(.gi_points(rn$p)$ymin)))
  expect_length(rn$messages, 0)
  ## several ids: the chosen cells' draws averaged within each replication
  ids <- c(101, 102, 103)
  p3 <- suppressMessages(plot(fit, type = "gap", id = ids))
  expect_equal(.gi_points(p3), .gi_band_by_hand(fit, ids), tolerance = 1e-10)
  ## every treated unit: the default gap plot, band included (unweighted fit)
  pall <- suppressMessages(plot(fit, type = "gap", id = fit$id[fit$tr]))
  expect_equal(.gi_points(pall), .gi_points(p0), tolerance = 1e-10)
  ## return.data carries the band
  out <- suppressMessages(plot(fit, type = "gap", id = 102, return.data = TRUE))
  expect_equal(out$data$estimate$CI.lower, ref$ymin, tolerance = 1e-10)
  expect_equal(out$data$estimate$CI.upper, ref$ymax, tolerance = 1e-10)
  ## the connected style draws the band as ribbons
  pc <- suppressMessages(plot(fit, type = "gap", id = 102, connected = TRUE))
  rib <- which(.gi_geoms(pc) == "GeomRibbon")
  expect_true(any(vapply(rib, function(i)
    any(is.finite(ggplot2::layer_data(pc, i)$ymin)), TRUE)))
  ## the normal rule at another level: the band uses the fit's ci.alpha
  a10 <- fit
  a10$ci.alpha <- 0.1
  expect_equal(.gi_points(suppressMessages(plot(a10, type = "gap", id = 102))),
               .gi_band_by_hand(fit, 102, a = 0.1), tolerance = 1e-10)
  expect_equal(.gi_points(suppressMessages(plot(a10, type = "gap", id = 102,
                                                plot.ci = "0.9"))),
               .gi_band_by_hand(fit, 102, a = 0.2), tolerance = 1e-10)
  ## a fit object from before fect 2.4.7, without ci.method and ci.alpha:
  ## the normal 95% band
  old <- fit
  old$ci.method <- NULL
  old$ci.alpha <- NULL
  expect_equal(.gi_points(suppressMessages(plot(old, type = "gap", id = 102))),
               pts, tolerance = 1e-12)
  ## an older fit that does not record its variance type: no band, and a
  ## message that says so
  nov <- fit
  nov$vartype <- NULL
  r <- .gi_plot(plot(nov, type = "gap", id = 102))
  expect_true(all(is.na(.gi_points(r$p)$ymin)))
  expect_length(r$messages, 1)
  expect_match(r$messages, "does not record its variance type", fixed = TRUE)
})

test_that("G9: the band follows the fit's ci.method and alpha", {
  skip_on_cran()
  fit <- .gi_param_fit(.gi_data("simgsynth"), ci.method = "basic", alpha = 0.1)
  expect_identical(fit$ci.method, "basic")
  expect_identical(fit$ci.alpha, 0.1)
  ## plot.ci = "0.95" draws the fit's own level, 1 - alpha, as the default
  ## gap plot does; "0.9" the one-sided bound, 1 - 2 * alpha
  pts <- .gi_points(suppressMessages(plot(fit, type = "gap", id = 102)))
  expect_equal(pts, .gi_band_by_hand(fit, 102, a = 0.1, rule = "basic"),
               tolerance = 1e-10)
  p90 <- .gi_points(suppressMessages(plot(fit, type = "gap", id = 102, plot.ci = "0.9")))
  expect_equal(p90, .gi_band_by_hand(fit, 102, a = 0.2, rule = "basic"),
               tolerance = 1e-10)
  pall <- suppressMessages(plot(fit, type = "gap", id = fit$id[fit$tr]))
  expect_equal(.gi_points(pall), .gi_points(suppressMessages(plot(fit, type = "gap"))),
               tolerance = 1e-10)
})

test_that("G10: the band reads the chosen units' draws when the treated units come after the controls", {
  skip_on_cran()
  simgsynth <- .gi_data("simgsynth")
  ## treated units 101-105 renamed 901-905, so they sort after the controls:
  ## fit$tr is 46:50, while fit$eff.boot keeps the treated units first
  tr <- unique(simgsynth$id[simgsynth$D == 1])
  simgsynth$id[simgsynth$id %in% tr] <- simgsynth$id[simgsynth$id %in% tr] + 800
  fit <- .gi_param_fit(simgsynth)
  expect_equal(fit$tr, 46:50)
  pts <- .gi_points(suppressMessages(plot(fit, type = "gap", id = 902)))
  expect_equal(pts, .gi_band_by_hand(fit, 902), tolerance = 1e-10)
  ## the draws at the unit's position in the panel belong to a control unit
  j <- which(fit$id == 902)
  wrong <- fit$eff[, j] - qnorm(0.975) * apply(fit$eff.boot[, j, ], 1, sd)
  expect_gt(max(abs(pts$ymin - wrong)), 0.1)
  pall <- suppressMessages(plot(fit, type = "gap", id = fit$id[fit$tr]))
  expect_equal(.gi_points(pall), .gi_points(suppressMessages(plot(fit, type = "gap"))),
               tolerance = 1e-10)
})

test_that("G11: without kept parametric draws, the gap plot for an id draws estimates and says how to get a band", {
  skip_on_cran()
  simgsynth <- .gi_data("simgsynth")
  args <- list(Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
               method = "gsynth", r = 2, CV = FALSE, se = TRUE,
               parallel = FALSE, seed = 1)
  ## parametric draws not kept
  fit <- suppressMessages(do.call(fect::fect, c(args, list(vartype = "parametric",
                                                           nboots = 20))))
  j <- which(fit$id == 102)
  r <- .gi_plot(plot(fit, type = "gap", id = 102))
  pts <- .gi_points(r$p)
  expect_equal(pts$y, as.numeric(fit$eff[, j]), tolerance = 1e-12)
  expect_true(all(is.na(pts$ymin)))
  expect_length(r$messages, 1)
  expect_match(r$messages, "did not keep its bootstrap draws", fixed = TRUE)
  expect_match(r$messages, "keep.sims = TRUE", fixed = TRUE)
  ## the jackknife, draws kept
  fit_j <- suppressMessages(do.call(fect::fect, c(args, list(vartype = "jackknife",
                                                             keep.sims = TRUE))))
  r <- .gi_plot(plot(fit_j, type = "gap", id = 102))
  expect_true(all(is.na(.gi_points(r$p)$ymin)))
  expect_length(r$messages, 1)
  expect_match(r$messages, "a band needs parametric bootstrap draws", fixed = TRUE)
  ## with plot.ci = "none" the user asked for no interval: no message
  expect_length(.gi_plot(plot(fit, type = "gap", id = 102,
                              plot.ci = "none"))$messages, 0)
  expect_length(.gi_plot(plot(fit_j, type = "gap", id = 102,
                              plot.ci = "none"))$messages, 0)
  ## the average gap plot of these fits keeps its interval
  expect_true(all(is.finite(.gi_points(suppressMessages(
    plot(fit_j, type = "gap")))$ymin)))
})

test_that("G12: a placebo fit with parametric draws: the band covers the placebo periods, with no p-value", {
  skip_on_cran()
  fit <- .gi_param_fit(.gi_data("simgsynth"), placeboTest = TRUE,
                       placebo.period = c(-2, 0))
  r <- .gi_plot(plot(fit, type = "gap", id = 101))
  expect_length(r$messages, 0)
  expect_true(all(.gi_text(r$p) == "")) # the placebo p-value is the average's
  ref <- .gi_band_by_hand(fit, 101)
  pl <- ref$x %in% c(-2, -1, 0)
  ## the placebo periods: triangles with line ranges; the others: point ranges
  lr <- which(.gi_geoms(r$p) == "GeomLinerange")
  band_pl <- do.call(rbind, lapply(lr, function(i)
    ggplot2::layer_data(r$p, i)[, c("x", "ymin", "ymax")]))
  band_pl <- band_pl[order(band_pl$x), ]
  expect_equal(band_pl$x, c(-2, -1, 0))
  expect_equal(band_pl$ymin, ref$ymin[pl], tolerance = 1e-10)
  expect_equal(band_pl$ymax, ref$ymax[pl], tolerance = 1e-10)
  pts <- .gi_points(r$p)
  expect_equal(pts$ymin, ref$ymin[!pl], tolerance = 1e-10)
})
