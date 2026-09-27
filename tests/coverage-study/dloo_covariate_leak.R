## dloo / loo with covariates: does a late control cell move the pre-treatment
## placebos?
##
## Adds +50 to one never-treated unit's outcome in the last period (after every
## cohort has adopted, so outside every dloo comparison sample) and reports how
## much each pre-treatment placebo moves, for
##   * dloo with covariates   (Y ~ D + X1 + X2): should not move (fixed in 2.4.7)
##   * dloo without covariates (Y ~ D):         does not move
##   * loo with covariates:                      moves, by design (the imputation
##     fit estimates beta from all untreated observations; see ch. 4)
##
## Run from the package root against the installed fect, or after
## devtools::load_all(). Before the fix (dev f9346ae) the dloo placebo at event
## time -4 moved from -0.064 to -0.357 on this panel.

suppressMessages(library(fect))

make_panel_x <- function(TT = 8, cohorts = c(4, 6, Inf), n_per = 25, sd = 0.4,
                         seed = 42) {
    set.seed(seed)
    rows <- list(); id <- 0
    for (g in cohorts) for (u in seq_len(n_per)) {
        id <- id + 1
        tr <- !is.infinite(g)
        x1 <- rnorm(TT) + (if (tr) 0.9 / g else 0.05) * (1:TT) * 3
        x2 <- rnorm(TT, sd = 0.5) + (if (tr) 0.2 else 0) * (1:TT)
        y  <- rnorm(1) + 0.1 * (1:TT) + 0.8 * x1 - 0.5 * x2 + rnorm(TT, 0, sd)
        d  <- as.integer(tr & (1:TT) >= g)
        if (tr) y <- y + d * 2.0
        rows[[length(rows) + 1]] <- data.frame(id = id, time = 1:TT, Y = y,
                                               D = d, X1 = x1, X2 = x2)
    }
    do.call(rbind, rows)
}

dat  <- make_panel_x()
nt   <- min(dat$id[ave(dat$D, dat$id, FUN = max) == 0])
last <- dat$id == nt & dat$time == max(dat$time)
dat2 <- dat
dat2$Y[last] <- dat2$Y[last] + 50

fit <- function(d, f, ...) {
    suppressMessages(fect(f, data = d, index = c("id", "time"), method = "fe",
                          force = "two-way", parallel = FALSE, nboots = 50,
                          seed = 1, ...))
}
cases <- list(
    "dloo, Y ~ D + X1 + X2" = list(f = Y ~ D + X1 + X2, dloo = TRUE),
    "dloo, Y ~ D"           = list(f = Y ~ D, dloo = TRUE),
    "loo,  Y ~ D + X1 + X2" = list(f = Y ~ D + X1 + X2, loo = TRUE, se = TRUE)
)
cat("fect", as.character(packageVersion("fect")), "\n")
for (nm in names(cases)) {
    a  <- cases[[nm]]
    f1 <- do.call(fit, c(list(d = dat), a))
    f2 <- do.call(fit, c(list(d = dat2), a))
    et <- rownames(f1$pre.est.att)
    d1 <- f1$pre.est.att[, "ATT"]; d2 <- f2$pre.est.att[, "ATT"]
    cat(sprintf("\n%s\n", nm))
    print(round(data.frame(event.time = as.numeric(et), before = d1, after = d2,
                           moved = d2 - d1, row.names = NULL), 4))
}
