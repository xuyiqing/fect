## Tests for the group.fe argument (added v2.4.5; closes #139).
## Design memo: statsclaw-workspace/fect/runs/2026-05-21-higher-level-fe.md

## ----------------------------------------------------------------------
## Helper: small nested panel (county within state, state-level treatment).
## Returns a balanced TT*N panel with non-trivial state structure.
## ----------------------------------------------------------------------
.make_nested_panel <- function(N = 40, TT = 10, n_states = 4, seed = 42) {
    set.seed(seed)
    df <- expand.grid(id = 1:N, time = 1:TT)
    df$state <- paste0("S", (df$id - 1) %% n_states + 1)
    df$D <- as.integer(df$state %in% c("S1", "S2") & df$time >= 6)
    df$Y <- 1 + 0.5 * df$D + rnorm(nrow(df), sd = 0.5)
    df
}

## ----------------------------------------------------------------------
## Point-estimate equivalence
## ----------------------------------------------------------------------

test_that("group.fe = 'state' is byte-equivalent to legacy index[3]", {
    df <- .make_nested_panel()
    fit_new <- fect(Y ~ D, data = df, index = c("id", "time"),
                    group.fe = "state", force = "time", se = FALSE)
    fit_old <- fect(Y ~ D, data = df, index = c("id", "time", "state"),
                    method = "cfe", r = 0, CV = FALSE,
                    force = "time", se = FALSE)
    expect_equal(fit_new$att.avg, fit_old$att.avg)
})

test_that("auto-route from method='fe' preserves the result", {
    df <- .make_nested_panel()
    fit_default <- fect(Y ~ D, data = df, index = c("id", "time"),
                        group.fe = "state", force = "time", se = FALSE)
    fit_explicit <- fect(Y ~ D, data = df, index = c("id", "time"),
                         group.fe = "state", method = "cfe", r = 0,
                         CV = FALSE, force = "time", se = FALSE)
    expect_equal(fit_default$att.avg, fit_explicit$att.avg)
})

## ----------------------------------------------------------------------
## Hard errors on unsupported methods (D3)
## ----------------------------------------------------------------------

test_that("method='ife' + group.fe hard-errors with guidance", {
    df <- .make_nested_panel()
    expect_error(
        fect(Y ~ D, data = df, index = c("id", "time"),
             group.fe = "state", method = "ife", se = FALSE),
        "not supported"
    )
})

test_that("method='mc' + group.fe hard-errors", {
    df <- .make_nested_panel()
    expect_error(
        fect(Y ~ D, data = df, index = c("id", "time"),
             group.fe = "state", method = "mc", se = FALSE),
        "not supported"
    )
})

test_that("method='both' + group.fe hard-errors", {
    df <- .make_nested_panel()
    expect_error(
        fect(Y ~ D, data = df, index = c("id", "time"),
             group.fe = "state", method = "both", se = FALSE),
        "not supported"
    )
})

test_that("method='gsynth' + group.fe hard-errors", {
    df <- .make_nested_panel()
    expect_error(
        fect(Y ~ D, data = df, index = c("id", "time"),
             group.fe = "state", method = "gsynth", se = FALSE),
        "not supported"
    )
})

## ----------------------------------------------------------------------
## Edge cases (D5)
## ----------------------------------------------------------------------

test_that("group.fe column missing from data hard-errors", {
    df <- .make_nested_panel()
    expect_error(
        fect(Y ~ D, data = df, index = c("id", "time"),
             group.fe = "nonexistent_col", se = FALSE),
        "not found"
    )
})

test_that("group.fe AND legacy index[3:] both used hard-errors", {
    df <- .make_nested_panel()
    expect_error(
        fect(Y ~ D, data = df, index = c("id", "time", "state"),
             group.fe = "state", se = FALSE),
        "OR extra index slots"
    )
})

test_that("non-character group.fe hard-errors", {
    df <- .make_nested_panel()
    expect_error(
        fect(Y ~ D, data = df, index = c("id", "time"),
             group.fe = 1L, se = FALSE),
        "character vector"
    )
})

test_that("group.fe overlapping index[1:2] is warned and dropped", {
    df <- .make_nested_panel()
    expect_warning(
        fect(Y ~ D, data = df, index = c("id", "time"),
             group.fe = c("id", "state"), force = "time", se = FALSE),
        "duplicate index"
    )
})

## ----------------------------------------------------------------------
## Nesting check (D5 edge 5, D3f) -- applies to BOTH group.fe and legacy
## ----------------------------------------------------------------------

test_that("non-nested group.fe hard-errors with offending units listed", {
    df <- .make_nested_panel()
    df$state[1:5] <- "BAD"  # county 1 now appears in two states across time
    expect_error(
        fect(Y ~ D, data = df, index = c("id", "time"),
             group.fe = "state", force = "time", se = FALSE),
        "not constant within"
    )
})

test_that("legacy index[3:] does NOT enforce nesting (supports cell-level interactions)", {
    ## Legacy index = c(unit, time, extra) syntax has historically supported
    ## both nested groupings AND cell-level interactions like region_time.
    ## The nesting check only fires for group.fe; legacy form keeps full scope.
    df <- .make_nested_panel()
    df$region_time <- as.numeric(df$state == "S1") + df$time / 100  # varies within id
    expect_no_error(
        fect(Y ~ D, data = df, index = c("id", "time", "region_time"),
             method = "cfe", force = "two-way", se = FALSE)
    )
})

## ----------------------------------------------------------------------
## cl auto-default + FALSE sentinel (D6)
## ----------------------------------------------------------------------

test_that("single-column group.fe auto-defaults cl to group.fe[1]", {
    df <- .make_nested_panel()
    fit <- fect(Y ~ D, data = df, index = c("id", "time"),
                group.fe = "state", force = "time", se = FALSE)
    expect_equal(fit$cl.label, "state")
})

test_that("cl = FALSE is rejected with a guiding error", {
    df <- .make_nested_panel()
    expect_error(
        fect(Y ~ D, data = df, index = c("id", "time"),
             group.fe = "state", force = "time", cl = FALSE, se = FALSE),
        "cl = FALSE is not supported"
    )
})

test_that("cl = index[1] explicitly clusters at the unit level", {
    df <- .make_nested_panel()
    fit <- fect(Y ~ D, data = df, index = c("id", "time"),
                group.fe = "state", force = "time", cl = "id", se = FALSE)
    expect_equal(fit$cl.label, "id")
})

test_that("cl = 'other_col' overrides the auto-default", {
    df <- .make_nested_panel()
    df$region <- substr(df$state, 1, 1)
    fit <- fect(Y ~ D, data = df, index = c("id", "time"),
                group.fe = "state", force = "time", cl = "region", se = FALSE)
    expect_equal(fit$cl.label, "region")
})

test_that("multi-column group.fe requires explicit cl", {
    df <- .make_nested_panel()
    df$region <- substr(df$state, 1, 1)
    expect_error(
        fect(Y ~ D, data = df, index = c("id", "time"),
             group.fe = c("state", "region"), force = "time", se = FALSE),
        "Multi-column group.fe requires explicit cl"
    )
})

## ----------------------------------------------------------------------
## Fit slots for print (D8)
## ----------------------------------------------------------------------

test_that("fit$group.fe and fit$cl.label are populated for downstream print", {
    df <- .make_nested_panel()
    fit <- fect(Y ~ D, data = df, index = c("id", "time"),
                group.fe = "state", force = "time", se = FALSE)
    expect_equal(fit$group.fe, "state")
    expect_equal(fit$cl.label, "state")
})

test_that("print(fit) surfaces Estimator + Fixed effects + Cluster SE", {
    df <- .make_nested_panel()
    fit <- fect(Y ~ D, data = df, index = c("id", "time"),
                group.fe = "state", force = "time", se = FALSE)
    out <- capture.output(print(fit))
    expect_true(any(grepl("^Estimator:", out)))
    expect_true(any(grepl("^Fixed effects:.*state", out)))
    expect_true(any(grepl("^Cluster SE:.*state", out)))
})

test_that("print(fit) shows the user's chosen cl in the Cluster SE line", {
    df <- .make_nested_panel()
    fit <- fect(Y ~ D, data = df, index = c("id", "time"),
                group.fe = "state", force = "time", cl = "id", se = FALSE)
    out <- capture.output(print(fit))
    expect_true(any(grepl("^Cluster SE:.*id", out)))
})
