## ======================================================================
## serial = "ar1": gsynth on an AR(1)-filtered panel
## ======================================================================
##
## The estimator (Xu, gsynth-revisit memo "Two ways to add AR(1) errors to
## GSC", 2026-10-16; prototype sim/ar1_options.R, `filtered_att`):
##
##   1. quasi-difference every unit in calendar time with one coefficient x,
##      V_it = Y_it - x Y_i,t-1, covariates too (the first period is lost);
##   2. fit the gsynth model to V: the filtered model keeps the factor
##      structure and the loadings (V_it = (1-x) alpha_i + xi_t - x xi_t-1 +
##      lambda_i'(F_t - x F_t-1) + (X_it - x X_i,t-1)'b + w_it);
##   3. undo the filter for each treated unit from its last observed
##      pretreatment outcome s_i: Yhat_{i,s} = Y_{i,s}, Yhat_it = x Yhat_i,t-1 +
##      Vhat_it, so the effect e_it = Y_it - Yhat_it = x e_i,t-1 + (V_it -
##      Vhat_it), e_{i,s} = 0, is on the outcome's scale.
##
## x is chosen jointly with the rank: it minimizes the controls' residual sum
## of squares of the rank-r fit to V(x) over [0, 1] (the conditional Gaussian
## quasi-likelihood of factors plus AR(1) errors; the filter has unit
## Jacobian). The rank is chosen by the treated-unit rule (R/cv-treated.R) on
## the filtered panel, and x and r are alternated from the levels rule's rank
## until the rank repeats. x = 0 is gsynth in levels without the first
## period. Inference (se = TRUE) is the parametric bootstrap of the filtered
## panel at the chosen (x, r), each draw mapped through step 3.
##
## Outer fit slots live on the outcome's scale on the filtered grid (T - 1
## periods); the slots that cannot be put there stay on the filtered scale and
## are listed in fit$ar1$v.scale.slots.

## the value of serial this file implements (one constant, so a rename is one
## line)
.FECT_SERIAL_AR1 <- "ar1"

## x is searched on this grid, then refined between the neighbours of the grid
## minimum
.FECT_AR1_GRID <- seq(0, 1, by = 0.05)

## cap on the rank/x alternation
.FECT_AR1_ROUNDS <- 5L

## state read by the treated-unit rule inside fect_nevertreated(): set by
## .fect_ar1()'s rule_on() only, around one call, and reset on exit
.fect_ar1_state <- new.env(parent = emptyenv())
.fect_ar1_state$exclude <- FALSE

## filtered-scale slots kept from the inner fit
.FECT_AR1_V_SLOTS <- c(
    "res", "res.full", "sigma2", "sigma2.fect", "IC", "PC", "factor",
    "lambda.co", "lambda.tr", "alpha.co", "alpha.tr", "xi", "mu", "est",
    "niter", "eff.pre", "eff.pre.equiv", "pre.sd", "wgt.implied", "CV.out",
    "CV.out.se", "cv.treated", "rmse"
)

## ----------------------------------------------------------------------
## Entry point (called from fect.default())
## ----------------------------------------------------------------------
.fect_ar1 <- function(data, Y, D, X, index, W, W.est, W.agg, group, na.rm,
                      force, time.component.from, em, r, lambda, CV,
                      cv.method, criterion, binary, method, se, vartype,
                      para.error, cl, ci.method, quantile.CI, nboots, alpha,
                      parallel, cores, tol, max.iteration, seed, min.T0,
                      max.missing, group.fe, cfe, Z, gamma, Q, kappa,
                      Z.param, Q.param, Q.type, balance.period, fill.missing,
                      placeboTest, carryoverTest, carryover.rm, loo, dloo,
                      permute, normalize, keep.sims, cm, loading.bound,
                      gamma.loading, gamma.loading.grid, cv.rule, call) {
    ## the five internal options, read once (never documented)
    opt.x <- getOption("fect.ar1.x", NULL)
    exclude.nc <- !isFALSE(getOption("fect.ar1.exclude.nonconverged", TRUE))
    mask.zero <- isTRUE(getOption("fect.ar1.mask.zero", FALSE))
    opt.rounds <- getOption("fect.ar1.rounds", NULL)
    K <- if (is.null(opt.rounds)) .FECT_AR1_ROUNDS else as.integer(opt.rounds)
    if (!is.null(opt.x) && (!is.numeric(opt.x) || length(opt.x) != 1L ||
                            !is.finite(opt.x) || opt.x < 0 || opt.x > 1)) {
        stop("option fect.ar1.x must be a number in [0, 1].", call. = FALSE)
    }
    if (!is.numeric(K) || length(K) != 1L || is.na(K) || K < 1L) {
        stop("option fect.ar1.rounds must be a positive whole number.",
             call. = FALSE)
    }

    spec <- .fect_ar1_check_args(
        data = data, Y = Y, D = D, X = X, index = index, W = W, W.est = W.est,
        W.agg = W.agg, group = group, force = force, r = r, cv.method = cv.method,
        criterion = criterion, binary = binary, method = method, se = se,
        vartype = vartype, group.fe = group.fe, cfe = cfe, Z = Z, gamma = gamma,
        Q = Q, kappa = kappa, Z.param = Z.param, Q.param = Q.param,
        Q.type = Q.type, balance.period = balance.period,
        fill.missing = fill.missing, placeboTest = placeboTest,
        carryoverTest = carryoverTest, loo = loo, dloo = dloo,
        permute = permute, cm = cm, loading.bound = loading.bound,
        gamma.loading = gamma.loading,
        gamma.loading.grid = gamma.loading.grid, cv.rule = cv.rule,
        min.T0 = min.T0, tol = tol, nboots = nboots, ci.method = ci.method,
        quantile.CI = quantile.CI, para.error = para.error
    )
    se <- spec$se
    vartype <- spec$vartype
    force.code <- spec$force.code
    has_a <- force.code %in% c(1, 3)
    has_xi <- force.code %in% c(2, 3)
    min.T0 <- spec$min.T0
    ci.method <- spec$ci.method

    id <- index[1]
    time <- index[2]
    prep <- .fect_ar1_prep(data, Y, D, X, id, time)
    vars <- c(Y, X)
    ## messages are muffled, and so are fect's own "did not converge" warnings
    ## (counted: they would repeat dozens of times in the search); every other
    ## warning passes
    st.nc <- new.env(parent = emptyenv())
    st.nc$n <- 0L
    quiet <- function(expr, count = TRUE) {
        suppressMessages(withCallingHandlers(expr, warning = function(w) {
            if (grepl("did not converge", conditionMessage(w), fixed = TRUE)) {
                if (count) st.nc$n <- st.nc$n + 1L
                invokeRestart("muffleWarning")
            }
        }))
    }
    ## an inner EM fit that stopped at max.iteration did not converge (the
    ## condition of fect's own warning)
    em_ok <- function(f) {
        !(length(f$niter) == 1L && is.finite(f$niter) &&
              f$niter >= max.iteration)
    }
    fail <- function(e) {
        stop("serial = \"ar1\": ", conditionMessage(e), call. = FALSE)
    }

    ## one gsynth fit through the front end on a (filtered) panel
    run <- function(fd, ..., min.T0) {
        fect.default(
            data = fd, Y = Y, D = D, X = X, index = index, method = "gsynth",
            force = force, tol = tol, max.iteration = max.iteration,
            normalize = normalize, max.missing = max.missing, na.rm = na.rm,
            serial = "none", min.T0 = min.T0, ...
        )
    }
    fd.levels <- prep$data
    ## the profile is always on the lag-masked panel (the same cells at every
    ## x); the rule, the final fit and everything after use V = Y without the
    ## mask at x exactly 0
    filter <- function(x) .fect_ar1_filter(prep, x, vars)
    filter.fit <- function(x) {
        .fect_ar1_filter(prep, x, vars, mask = !(x == 0 && !mask.zero))
    }
    min.rule <- max(1, min.T0 - 1)  # the filter costs each unit one period
    min.fit <- function(rr) max(min.T0 - 1, rr + 1)

    ## ---- SSR(x; r), the controls' residual sum of squares at rank r ------
    ## the control cells that enter are identical at every x (asserted)
    st <- new.env(parent = emptyenv())
    st$sig <- NULL
    ssr_at <- function(x, rr) {
        f <- tryCatch(
            quiet(run(filter(x), r = rr, CV = FALSE, se = FALSE,
                      parallel = FALSE, min.T0 = min.fit(rr))),
            error = function(e) NULL
        )
        st$em <- NA
        if (is.null(f)) return(Inf)
        st$em <- em_ok(f)
        m <- f$II[, f$co, drop = FALSE] == 1
        s <- list(ids = as.character(f$id[f$co]), m = m)
        if (is.null(st$sig)) {
            st$sig <- s
        } else if (!identical(st$sig$ids, s$ids) || !identical(st$sig$m, s$m)) {
            stop("serial = \"ar1\": the control cells that enter the profile ",
                 "differ across x; the profile is not comparable.",
                 call. = FALSE)
        }
        v <- sum(f$res[m]^2)
        if (!is.finite(v)) Inf else v
    }
    ## x(r): grid, one optimize() between the grid minimum's neighbours, exact
    ## boundary. Returns x and the profile of every point evaluated.
    x_of_r <- function(rr) {
        st$sig <- NULL
        pts <- new.env(parent = emptyenv())
        pts$x <- numeric(0)
        pts$ssr <- numeric(0)
        pts$em <- logical(0)
        f <- function(x) {
            v <- ssr_at(x, rr)
            pts$x <- c(pts$x, x)
            pts$ssr <- c(pts$ssr, v)
            pts$em <- c(pts$em, st$em)
            v
        }
        g <- .FECT_AR1_GRID
        s <- vapply(g, f, 0)
        if (!any(is.finite(s))) {
            stop("serial = \"ar1\": the filtered panel cannot be fitted at ",
                 "rank ", rr, " at any x.", call. = FALSE)
        }
        k <- which.min(s)
        bracket <- c(g[max(1L, k - 1L)], g[min(length(g), k + 1L)])
        o <- stats::optimize(f, interval = bracket, tol = 1e-5)
        x.ref <- o$minimum
        xx <- if (is.finite(o$objective) && o$objective < s[k]) x.ref else g[k]
        o.ord <- order(pts$x)
        list(x = xx, profile = data.frame(x = pts$x[o.ord], ssr = pts$ssr[o.ord],
                                          em.converged = pts$em[o.ord]))
    }

    ## ---- the treated-unit rule on a panel ----------------------------------
    r.range <- spec$r.range
    rule_on <- function(fd, min.T0.rule, exclude = FALSE) {
        .fect_ar1_state$exclude <- exclude
        on.exit(.fect_ar1_state$exclude <- FALSE, add = TRUE)
        f <- tryCatch(
            quiet(run(fd, r = r.range, CV = TRUE, cv.method = .FECT_CV_PRETREATMENT,
                      se = FALSE, parallel = FALSE, min.T0 = min.T0.rule)),
            error = fail
        )
        list(r = as.integer(f$r.cv[1]), rule = f$cv.treated)
    }

    ## ---- choose (x, r) --------------------------------------------------------
    rule.levels <- NULL
    rule.last <- NULL
    profile <- data.frame(x = numeric(0), ssr = numeric(0),
                          em.converged = logical(0))
    profile.rank <- NA_integer_
    converged <- NA
    x.fixed <- !is.null(opt.x)
    cycle.df <- data.frame(round = integer(0), x = numeric(0), r = integer(0),
                           chosen = logical(0))
    cycled <- FALSE
    trace <- list()
    if (spec$fixed.r) {
        r.k <- spec$r
        if (x.fixed) {
            x.k <- as.numeric(opt.x)
        } else {
            xr <- x_of_r(r.k)
            x.k <- xr$x
            profile <- xr$profile
            profile.rank <- as.integer(r.k)
        }
        path <- data.frame(round = 1L, x = x.k, r = as.integer(r.k),
                           att = NA_real_, em.converged = NA)
        rounds <- 1L
        final.round <- 1L
        message("serial = \"ar1\": x = ", format(round(x.k, 4)), ", rank ", r.k)
    } else {
        r0 <- rule_on(fd.levels, min.T0)
        rule.levels <- r0$rule
        ## one round: the x-step at the rank fed in, then the rule on V(x)
        step <- function(r.prev) {
            if (x.fixed) {
                x.s <- as.numeric(opt.x)
                prof <- profile
                prof.r <- NA_integer_
            } else {
                xr <- x_of_r(r.prev)
                x.s <- xr$x
                prof <- xr$profile
                prof.r <- as.integer(r.prev)
            }
            rl <- rule_on(filter.fit(x.s), min.rule, exclude = exclude.nc)
            if (!is.null(rl$rule$eligible) && any(!rl$rule$eligible)) {
                bad <- rl$rule$table$r[!rl$rule$eligible]
                message("serial = \"ar1\": rank", if (length(bad) > 1L) "s " else " ",
                        paste(bad, collapse = ", "), " not eligible at x = ",
                        formatC(x.s, format = "f", digits = 2),
                        " (EM did not converge)")
            }
            message("serial = \"ar1\": x = ", format(round(x.s, 4)), ", rank ",
                    rl$r)
            list(x = x.s, r = rl$r, profile = prof, profile.rank = prof.r,
                 rule = rl$rule)
        }
        if (x.fixed) {
            ## the map does not depend on the rank: one round, no alternation
            s1 <- step(r0$r)
            alt <- list(
                path = data.frame(round = 0:1, x = c(NA, s1$x),
                                  r = c(r0$r, s1$r)),
                final = 1L, converged = NA, cycled = FALSE,
                cycle = integer(0), capped = FALSE, trace = list(s1)
            )
        } else {
            alt <- .fect_ar1_alternate(step, r0$r, K)
        }
        trace <- alt$trace
        fin <- trace[[alt$final]]
        x.k <- fin$x
        r.k <- fin$r
        profile <- fin$profile
        profile.rank <- fin$profile.rank
        rule.last <- fin$rule
        converged <- alt$converged
        cycled <- alt$cycled
        path <- data.frame(alt$path, att = NA_real_, em.converged = NA)
        rounds <- nrow(path) - 1L
        final.round <- as.integer(alt$final)
        if (cycled) {
            cy <- alt$path[match(alt$cycle, alt$path$round), ]
            cycle.df <- data.frame(round = as.integer(cy$round), x = cy$x,
                                   r = as.integer(cy$r),
                                   chosen = cy$round == final.round)
            fmt <- function(x, r) sprintf("(%.3f, %d)", x, r)
            prs <- fmt(cycle.df$x, cycle.df$r)
            lst <- if (length(prs) == 2L) {
                paste(prs, collapse = " and ")
            } else {
                paste0(paste(prs[-length(prs)], collapse = ", "), " and ",
                       prs[length(prs)])
            }
            ch <- cycle.df[cycle.df$chosen, ]
            warning("serial = \"ar1\": x and the rank cycle between (x, r) = ",
                    lst, "; using the pair with the lowest rank, ",
                    fmt(ch$x, ch$r), ".", call. = FALSE)
        } else if (isTRUE(alt$capped) && is.null(opt.rounds)) {
            warning("serial = \"ar1\": the rank and x did not settle in ",
                    K, " rounds; using the last pair (x = ",
                    format(round(x.k, 4)), ", rank ", r.k, ").",
                    call. = FALSE)
        }
    }
    if (x.k <= 0 || x.k >= 1) {
        message("serial = \"ar1\": x is at the boundary (x = ", x.k, ")",
                if (x.k <= 0) ": gsynth in levels without the first period"
                else ": the unit effect is free (a unit-specific drift)")
    }

    ## ---- final fit at (x, r) -----------------------------------------------
    min.final <- min.fit(r.k)
    inner <- tryCatch(
        quiet(
            if (se) {
                run(filter.fit(x.k), r = r.k, CV = FALSE, se = TRUE,
                    vartype = "parametric", nboots = nboots,
                    para.error = para.error, alpha = alpha, ci.method = ci.method,
                    parallel = parallel, cores = cores, seed = seed,
                    keep.sims = TRUE, min.T0 = min.final)
            } else {
                run(filter.fit(x.k), r = r.k, CV = FALSE, se = FALSE,
                    parallel = FALSE, min.T0 = min.final)
            },
            count = FALSE
        ),
        error = fail
    )
    em.final <- em_ok(inner)
    if (!em.final) {
        warning("serial = \"ar1\": the EM fit at x = ",
                formatC(x.k, format = "f", digits = 2), ", rank r = ", r.k,
                " did not converge within max.iteration = ", max.iteration,
                "; the estimate can change with max.iteration. Raise ",
                "max.iteration or choose r.", call. = FALSE)
    }

    ctx <- list(prep = prep, id = id, time = time, Y = Y, D = D, X = X,
                has_a = has_a, has_xi = has_xi, force.code = force.code,
                mask.zero = mask.zero)
    rb <- .fect_ar1_rebuild(inner, ctx, x.k)
    if (length(rb$not.rebuilt) > 0L) {
        warning("serial = \"ar1\": ", length(rb$not.rebuilt),
                " treated unit(s) could not be rebuilt on the outcome scale ",
                "and are left out: ",
                paste(utils::head(rb$not.rebuilt, 10), collapse = ", "),
                if (length(rb$not.rebuilt) > 10) ", ..." else "",
                call. = FALSE)
    }

    ## the outcome-scale ATT and the EM status at each pair of the path (cheap:
    ## one fit each)
    final.row <- which(path$round == final.round)
    for (i in seq_len(nrow(path))) {
        if (is.na(path$x[i])) next
        if (i == final.row) {
            path$att[i] <- rb$att.avg
            path$em.converged[i] <- em.final
            next
        }
        fi <- tryCatch(
            quiet(run(filter.fit(path$x[i]), r = path$r[i], CV = FALSE,
                      se = FALSE, parallel = FALSE,
                      min.T0 = min.fit(path$r[i]))),
            error = function(e) NULL
        )
        if (is.null(fi)) next
        path$em.converged[i] <- em_ok(fi)
        path$att[i] <- tryCatch(
            .fect_ar1_rebuild(fi, ctx, path$x[i])$att.avg,
            error = function(e) NA_real_)
    }

    ## ranks passed over, rounds without an eligible rank
    excluded <- do.call(rbind, c(list(data.frame(
        round = integer(0), x = numeric(0), r = integer(0))), lapply(
        seq_along(trace), function(k) {
            ru <- trace[[k]]$rule
            if (is.null(ru$eligible) || all(ru$eligible)) return(NULL)
            data.frame(round = as.integer(k), x = trace[[k]]$x,
                       r = as.integer(ru$table$r[!ru$eligible]))
        })))
    fallback.rounds <- as.integer(which(vapply(trace, function(t) {
        isTRUE(t$rule$fallback)
    }, TRUE)))

    out <- .fect_ar1_outer(inner, rb, ctx, x.k)
    ## the final round's rule output, as on a levels fit (also in fit$ar1$rule)
    if (!is.null(rule.last)) out$cv.treated <- rule.last
    out <- .fect_ar1_inference(out, inner, rb, x.k, se = se, alpha = alpha,
                               ci.method = ci.method, keep.sims = keep.sims)
    out$ar1 <- .fect_ar1_info(
        x = x.k, x.fixed = x.fixed, inner = inner, rb = rb, path = path,
        converged = converged, rounds = rounds, profile = profile,
        profile.rank = profile.rank, rule = rule.last, rule.levels = rule.levels,
        prep = prep, data = data, id = id, min.T0 = min.final,
        extra = list(
            cycled = cycled, cycle = cycle.df, final.round = final.round,
            em.converged = em.final, niter = as.integer(inner$niter[1]),
            excluded = excluded, fallback.rounds = fallback.rounds,
            search.nonconverged = st.nc$n, lag.mask = rb$lag.mask
        )
    )
    out$serial <- .FECT_SERIAL_AR1
    out$call <- call
    class(out) <- "fect"
    out
}

## ----------------------------------------------------------------------
## Argument checks (section 3.1 of the build spec)
## ----------------------------------------------------------------------
.fect_ar1_check_args <- function(data, Y, D, X, index, W, W.est, W.agg, group,
                                 force, r, cv.method, criterion, binary, method,
                                 se, vartype, group.fe, cfe, Z, gamma, Q, kappa,
                                 Z.param, Q.param, Q.type, balance.period,
                                 fill.missing, placeboTest, carryoverTest, loo,
                                 dloo, permute, cm, loading.bound, gamma.loading,
                                 gamma.loading.grid, cv.rule, min.T0, tol,
                                 nboots, ci.method, quantile.CI, para.error) {
    no <- function(what) {
        stop("serial = \"ar1\" does not support ", what, ".", call. = FALSE)
    }
    if (!identical(method, "gsynth")) {
        stop("serial = \"ar1\" requires method = \"gsynth\".", call. = FALSE)
    }
    if (!is.data.frame(data)) data <- as.data.frame(data)
    if (is.null(Y) || is.null(D) || is.null(index)) {
        stop("serial = \"ar1\" needs the outcome, the treatment and the index.",
             call. = FALSE)
    }
    if (length(index) != 2L || !all(index %in% colnames(data))) {
        stop("serial = \"ar1\": \"index\" option misspecified. Try, for ",
             "example, index = c(\"unit.id\", \"time\").", call. = FALSE)
    }
    miss.cols <- setdiff(c(Y, D, X), colnames(data))
    if (length(miss.cols) > 0L) {
        stop("serial = \"ar1\": variable \"", miss.cols[1],
             "\" is not in the data set.", call. = FALSE)
    }
    if (!is.null(W) || !is.null(W.est) || !is.null(W.agg)) no("weights (W, W.est, W.agg)")
    if (!is.null(group)) no("group")
    if (!is.null(cfe)) no("cfe")
    if (!is.null(Z) || !is.null(gamma) || !is.null(Q) || !is.null(kappa) ||
        !is.null(Z.param) || !is.null(Q.param) || !is.null(Q.type)) {
        no("Z, gamma, Q, kappa (complex fixed effects)")
    }
    if (!is.null(group.fe)) no("group.fe")
    if (isTRUE(binary)) no("binary outcomes")
    if (!identical(loading.bound, "none") || !is.null(gamma.loading) ||
        !is.null(gamma.loading.grid)) {
        no("loading.bound")
    }
    if (isTRUE(cm)) no("cm = TRUE")
    if (!is.null(balance.period)) no("balance.period")
    if (isTRUE(fill.missing)) no("fill.missing = TRUE")
    if (isTRUE(placeboTest)) no("placeboTest")
    if (isTRUE(carryoverTest)) no("carryoverTest")
    if (isTRUE(loo)) no("loo = TRUE")
    if (isTRUE(dloo)) no("dloo = TRUE")
    if (isTRUE(permute)) no("permute = TRUE")

    force.code <- switch(as.character(force), none = 0, unit = 1, time = 2,
                         `two-way` = 3, `0` = 0, `1` = 1, `2` = 2, `3` = 3,
                         NA)
    if (is.na(force.code)) {
        stop("\"force\" option misspecified; choose from c(\"none\", ",
             "\"unit\", \"time\", \"two-way\").", call. = FALSE)
    }

    se <- isTRUE(as.logical(se))
    if (se) {
        if (!vartype %in% c("bootstrap", "jackknife", "parametric")) {
            stop("\"vartype\" must be one of \"bootstrap\", \"jackknife\", ",
                 "or \"parametric\".", call. = FALSE)
        }
        if (identical(vartype, "jackknife")) {
            stop("serial = \"ar1\" supports vartype = \"parametric\" only ",
                 "(the parametric bootstrap of the filtered panel).",
                 call. = FALSE)
        }
        if (identical(vartype, "bootstrap")) {
            message("serial = \"ar1\": using vartype = \"parametric\", the ",
                    "bootstrap of the filtered panel.")
            vartype <- "parametric"
        }
        if (!is.numeric(nboots) || nboots <= 0) {
            stop("\"nboots\" option misspecified. Try, for example, ",
                 "nboots = 1000.", call. = FALSE)
        }
    }
    if (!is.null(quantile.CI)) {
        ci.method <- if (isTRUE(as.logical(quantile.CI))) "basic" else "normal"
    }
    if (!ci.method %in% c("normal", "basic")) {
        stop("`ci.method` must be one of \"normal\" or \"basic\"; got \"",
             ci.method, "\".", call. = FALSE)
    }

    ## the rank: a single number is fixed; a range or NULL is chosen by the
    ## rule (CV is ignored)
    fixed.r <- !is.null(r) && length(r) == 1L
    if (!is.null(r) && (!is.numeric(r) || anyNA(r) || any(r < 0))) {
        stop("\"r\" option misspecified. The number of factors must be ",
             "non-negative.", call. = FALSE)
    }
    r.range <- if (fixed.r) NULL else if (is.null(r)) c(0, 5) else range(r)
    if (!fixed.r && length(unique(r.range)) == 1L) {
        fixed.r <- TRUE
        r <- r.range[1]
    }
    if (!identical(cv.rule, "min")) {
        stop("serial = \"ar1\" ranks with one pooled loss per rank (no folds); ",
             "use cv.rule = \"min\".", call. = FALSE)
    }
    if (!identical(criterion, "mspe")) {
        stop("serial = \"ar1\" ranks with one pooled loss per rank; use ",
             "criterion = \"mspe\".", call. = FALSE)
    }
    if (!is.character(cv.method) || length(cv.method) != 1L) {
        stop("serial = \"ar1\": \"cv.method\" must be a single string.",
             call. = FALSE)
    }
    if (cv.method %in% c("block", "all_units", "treated_units", "loo") ||
        (!cv.method %in% c("rolling", .FECT_CV_PRETREATMENT) &&
         !is.na(pmatch(cv.method, c("block", "all_units", "treated_units",
                                    "loo"))))) {
        stop("serial = \"ar1\" ranks with cv.method = \"",
             .FECT_CV_PRETREATMENT, "\" only (not \"", cv.method, "\").",
             call. = FALSE)
    }
    if (identical(cv.method, "rolling") && !fixed.r) {
        message("serial = \"ar1\": ranking by cv.method = \"",
                .FECT_CV_PRETREATMENT, "\" (the treated-unit rule).")
    }
    if (is.null(min.T0)) {
        min.T0 <- 5
    } else if (!is.numeric(min.T0) || min.T0 <= 0) {
        stop("\"min.T0\" option should be larger than 0.", call. = FALSE)
    }
    if (!is.numeric(tol) || tol <= 0) {
        stop("\"tol\" option misspecified. Try using the default option.",
             call. = FALSE)
    }
    for (v in X) {
        xc <- data[[v]]
        if (!is.numeric(xc) && !is.logical(xc)) {
            stop("serial = \"ar1\": covariate \"", v, "\" must be numeric.",
                 call. = FALSE)
        }
    }
    list(se = se, vartype = vartype, force.code = force.code, fixed.r = fixed.r,
         r = if (fixed.r) as.integer(r[1]) else NULL, r.range = r.range,
         min.T0 = min.T0, ci.method = ci.method)
}

## ----------------------------------------------------------------------
## The filter (calendar time)
## ----------------------------------------------------------------------
## Everything that does not depend on x: the sorted panel, each row's position
## on the time grid, and the row of its lag. A period no unit has is not a
## period: the lag of a row is the same unit's row at the previous grid
## position.
.fect_ar1_prep <- function(data, Y, D, X, id, time) {
    data <- data[, unique(c(id, time, Y, D, X)), drop = FALSE]
    tcol <- data[[time]]
    if (is.factor(tcol)) {
        f <- droplevels(tcol)
        lev <- levels(f)
        num <- suppressWarnings(as.numeric(lev))
        data[[time]] <- if (length(lev) > 0 && !anyNA(num) && all(diff(num) > 0)) {
            num[as.integer(f)]
        } else {
            as.integer(f)
        }
    } else if (is.character(tcol)) {
        num <- suppressWarnings(as.numeric(tcol))
        if (any(!is.na(tcol) & is.na(num))) {
            stop("serial = \"ar1\": the time index \"", time, "\" is ",
                 "character and some values are not numbers.", call. = FALSE)
        }
        data[[time]] <- num
    }
    if (anyNA(data[[id]]) || anyNA(data[[time]])) {
        stop("serial = \"ar1\": \"index\" should not have missing values.",
             call. = FALSE)
    }
    ## infinite outcomes are missing, as fect treats them
    yv <- data[[Y]]
    if (any(is.infinite(yv))) {
        yv[is.infinite(yv)] <- NA
        data[[Y]] <- yv
        warning("Detect infinite values in outcome, automatically replace ",
                "them with NA.\n", call. = FALSE)
    }
    data <- data[order(data[[id]], data[[time]]), , drop = FALSE]
    rownames(data) <- NULL
    tg <- sort(unique(data[[time]]))
    pos <- match(data[[time]], tg)
    key <- paste(as.character(data[[id]]), pos, sep = "\r")
    if (anyDuplicated(key) > 0L) {
        stop("Observations are not uniquely defined by unit and time ",
             "indicators.", call. = FALSE)
    }
    lag.row <- match(paste(as.character(data[[id]]), pos - 1L, sep = "\r"), key)
    list(data = data, tg = tg, pos = pos, lag.row = lag.row,
         keep = which(!is.na(lag.row)), T = length(tg))
}

## V = v - x v_lag for the outcome and the covariates; a row exists iff its lag
## row exists in the data, and a value is missing if either is. D stays at its
## own period.
.fect_ar1_filter <- function(prep, x, vars, mask = TRUE) {
    if (!mask) {
        ## x exactly 0 (the lag is not used): V = Y on the panel without its
        ## first period, no cell removed for a missing lag
        k <- which(prep$pos >= 2L)
        out <- prep$data[k, , drop = FALSE]
        for (v in vars) out[[v]] <- as.numeric(prep$data[[v]][k])
        rownames(out) <- NULL
        return(out)
    }
    k <- prep$keep
    lr <- prep$lag.row[k]
    out <- prep$data[k, , drop = FALSE]
    for (v in vars) {
        out[[v]] <- as.numeric(prep$data[[v]][k]) - x * as.numeric(prep$data[[v]][lr])
    }
    rownames(out) <- NULL
    out
}

## ----------------------------------------------------------------------
## The alternation as a pure function of one map on ranks
## ----------------------------------------------------------------------
## step(r) runs the x-step at rank r and the rule at V(x): it returns list(x, r,
## ...) (extras are kept untouched in `trace`) and is called once per round, in
## order, never after the decision. Stops at a fixed point (converged), at the
## first repeated rank (cycle: the pair with the lowest rank in the cycle is
## the result) or at the cap K.
.fect_ar1_alternate <- function(step, r0, K = 5L) {
    hist <- as.integer(r0)
    xs <- numeric(0)
    trace <- list()
    converged <- FALSE
    cycled <- FALSE
    cycle <- integer(0)
    capped <- FALSE
    final <- NA_integer_
    done <- FALSE
    for (k in seq_len(K)) {
        s <- step(hist[k])
        trace[[k]] <- s
        xs <- c(xs, s$x)
        hist <- c(hist, as.integer(s$r))
        r.k <- hist[k + 1L]
        if (r.k == hist[k]) {
            converged <- TRUE
            final <- k
            done <- TRUE
            break
        }
        j <- if (k >= 2L) which(hist[seq_len(k - 1L)] == r.k)[1L] - 1L else NA_integer_
        if (!is.na(j)) {
            cycle <- (j + 1L):k
            cycled <- TRUE
            final <- cycle[which.min(hist[cycle + 1L])]
            done <- TRUE
            break
        }
    }
    if (!done) {
        capped <- TRUE
        final <- as.integer(K)
    }
    n <- length(trace)
    list(
        path = data.frame(round = 0:n, x = c(NA_real_, xs), r = hist[seq_len(n + 1L)]),
        final = as.integer(final), converged = converged, cycled = cycled,
        cycle = as.integer(cycle), capped = capped, trace = trace
    )
}

## ----------------------------------------------------------------------
## Rebuilding the counterfactual (the inverse recursion)
## ----------------------------------------------------------------------
## `inner` is the fit to the filtered panel at x. For each treated unit: the
## model's prediction Vhat at every cell where the covariates are available,
## the anchor s_i, the chain Yhat_it = x Yhat_i,t-1 + Vhat_it, the effects
## e_it = Y_it - Yhat_it at observed treated cells, and the one-step-ahead
## errors V - Vhat everywhere else.
.fect_ar1_rebuild <- function(inner, ctx, x) {
    prep <- ctx$prep
    ids <- as.character(inner$id)
    N1 <- length(ids)
    rows <- match(as.character(inner$rawtime), as.character(prep$tg))
    if (anyNA(rows) || min(rows) < 2L) {
        stop("serial = \"ar1\": the filtered fit's periods do not match the ",
             "panel's time grid.", call. = FALSE)
    }
    T1 <- length(rows)
    ui <- match(as.character(prep$data[[ctx$id]]), ids)
    sel <- !is.na(ui)
    wide <- function(v) {
        m <- matrix(NA_real_, prep$T, N1)
        m[cbind(prep$pos[sel], ui[sel])] <- as.numeric(v[sel])
        m
    }
    Yfull <- wide(prep$data[[ctx$Y]])
    Yo <- Yfull[rows, , drop = FALSE]
    p.all <- length(ctx$X)
    Xfull <- array(NA_real_, dim = c(prep$T, N1, p.all))
    for (k in seq_len(p.all)) Xfull[, , k] <- wide(prep$data[[ctx$X[k]]])
    ## observed on the outcome's scale: Y and every covariate present
    obs <- is.finite(Yo)
    for (k in seq_len(p.all)) obs <- obs & is.finite(Xfull[rows, , k])
    ## filtered covariates, on the inner grid
    Xf <- array(NA_real_, dim = c(T1, N1, p.all))
    for (k in seq_len(p.all)) {
        ## x exactly 0 without the lag mask: V = X, the lag is not used
        ## (0 * NA would be NA)
        Xf[, , k] <- if (x == 0 && !isTRUE(ctx$mask.zero)) {
            Xfull[rows, , k]
        } else {
            Xfull[rows, , k] - x * Xfull[rows - 1L, , k]
        }
    }

    tr <- inner$tr
    ## the treatment indicator as the data have it: a cell whose filtered row is
    ## missing (its lag, or the cell itself, is absent) still has its own D
    ## when the original row exists; a cell absent from the data carries the
    ## previous D, as fect's own relative-period rule does
    Dw <- wide(prep$data[[ctx$D]])[rows, , drop = FALSE]
    ID.w <- !is.na(Dw)
    Dw[!ID.w] <- 0
    Dm <- Dw
    Vobs <- inner$I.dat == 1
    F <- if (!is.null(inner$factor)) as.matrix(inner$factor) else matrix(0, T1, 0)
    beta <- inner$beta
    ## covariates the fit dropped (collinear) have NA coefficients and no
    ## part in the model
    bvec <- if (p.all == 0L) numeric(0) else as.numeric(beta)
    keepcov <- which(!is.na(bvec))
    mu <- if (is.null(inner$mu)) 0 else as.numeric(inner$mu)
    xi <- if (ctx$has_xi) as.numeric(inner$xi) else rep(0, T1)
    alpha <- if (ctx$has_a) as.numeric(inner$alpha.tr) else rep(0, length(tr))
    Lam <- if (ncol(F) > 0) as.matrix(inner$lambda.tr) else matrix(0, length(tr), 0)

    ## Vhat at every cell of the treated units (NA where a used covariate's
    ## filtered value is missing)
    Vhat <- matrix(NA_real_, T1, length(tr))
    for (w in seq_along(tr)) {
        j <- tr[w]
        v <- mu + xi + alpha[w]
        if (ncol(F) > 0) v <- v + drop(F %*% Lam[w, ])
        for (k in keepcov) v <- v + Xf[, j, k] * bvec[k]
        Vhat[, w] <- v
    }
    ## consistency: at every observed cell of a treated unit the prediction
    ## equals the fit's own counterfactual (this catches wrong offsets under
    ## covariates, force and normalize)
    scale <- max(1, stats::sd(as.numeric(Yo[obs]), na.rm = TRUE), na.rm = TRUE)
    Yct.tr <- inner$Y.ct[, tr, drop = FALSE]
    cmp <- Vobs[, tr, drop = FALSE] & is.finite(Yct.tr) & is.finite(Vhat)
    if (any(cmp)) {
        dev <- abs(Vhat[cmp] - Yct.tr[cmp])
        if (max(dev) > 1e-8 * scale) {
            wbad <- which(colSums(cmp & abs(Vhat - Yct.tr) > 1e-8 * scale,
                                  na.rm = TRUE) > 0)[1]
            stop("serial = \"ar1\": the rebuilt model prediction differs from ",
                 "the fit's counterfactual (max deviation ",
                 format(max(dev), digits = 3), ") for treated unit ",
                 ids[tr[wbad]], ".", call. = FALSE)
        }
    }

    ## relative periods from the data's D (fect's get_term; the inner fit's own
    ## T.on has the wrong clock after a cell dropped by the filter)
    Ton <- inner$T.on
    for (j in which(colSums(Dm != inner$D.dat) > 0)) {
        if (any(ID.w[, j])) Ton[, j] <- get_term(Dm[, j], as.numeric(ID.w[, j]), type = "on")
    }

    ## effects
    eff <- inner$eff
    anchor <- rep(NA_integer_, length(tr))
    n.post <- integer(length(tr))
    Yhat <- matrix(NA_real_, T1, length(tr))
    ## x exactly 0 on the unmasked panel: no recursion (0 * NA is NA); the
    ## effects are the inner fit's own at the observed cells with a finite
    ## prediction
    plain <- x == 0 && !isTRUE(ctx$mask.zero)
    for (w in seq_along(tr)) {
        j <- tr[w]
        pre.ok <- which(Dm[, j] == 0 & obs[, j])
        if (length(pre.ok) == 0L) next
        s <- max(pre.ok)
        anchor[w] <- s
        if (plain) {
            Yhat[, w] <- Vhat[, w]
            fin <- Vobs[, j] & obs[, j] & is.finite(Vhat[, w]) &
                is.finite(inner$eff[, j])
            n.post[w] <- sum(Dm[, j] == 1 & fin)
            col <- rep(NA_real_, T1)
            col[fin] <- inner$eff[fin, j]
            eff[, j] <- col
            next
        }
        yh <- rep(NA_real_, T1)
        yh[s] <- Yo[s, j]
        if (s < T1) {
            for (t in (s + 1L):T1) yh[t] <- x * yh[t - 1L] + Vhat[t, w]
        }
        Yhat[, w] <- yh
        post <- which(seq_len(T1) > s & Dm[, j] == 1 & obs[, j] & is.finite(yh))
        n.post[w] <- length(post)
        ## the unit's column: one-step-ahead errors up to the anchor, the
        ## rolled-forward effects after it
        col <- rep(NA_real_, T1)
        pre.cells <- which(seq_len(T1) <= s & is.finite(inner$eff[, j]) &
                               Vobs[, j] & obs[, j])
        col[pre.cells] <- inner$eff[pre.cells, j]
        col[post] <- Yo[post, j] - yh[post]
        eff[, j] <- col
    }
    ok <- !is.na(anchor) & n.post > 0L
    ## units that cannot be rebuilt contribute nothing
    for (w in which(!ok)) eff[, tr[w]] <- NA_real_
    eff[!obs] <- NA_real_
    D.mat <- Dm

    ## ATT as fect: observed cells with finite eff
    complete <- !is.na(eff)
    att.avg <- sum(eff[complete] * D.mat[complete]) / sum(D.mat[complete])
    tr.pos <- which(colSums(D.mat) > 0)
    att.unit <- vapply(tr.pos, function(j) {
        okc <- !is.na(eff[, j])
        d <- sum(D.mat[okc, j])
        if (d > 0) sum(eff[okc, j] * D.mat[okc, j]) / d else NA_real_
    }, 0)
    att.avg.unit <- mean(att.unit, na.rm = TRUE)

    ## event-time ATT over all cells with a finite effect
    eff.v <- c(eff)
    t.on <- c(Ton)
    use <- !is.na(eff.v) & !is.na(t.on)
    time.on <- sort(unique(t.on[use]))
    att.on <- as.numeric(tapply(eff.v[use], t.on[use], mean))
    count.on <- as.numeric(table(t.on[use]))

    ## counterfactuals: Y - eff where eff is defined; the chain value at the
    ## treated cells where Y is unobserved
    Yo.fill <- Yo
    Yo.fill[!obs] <- NA_real_
    Y.ct <- Yo.fill - eff
    for (w in seq_along(tr)) {
        j <- tr[w]
        cells <- which(!is.na(anchor[w]) & ok[w] &
                           (plain | seq_len(T1) > anchor[w]) &
                           !obs[, j] & is.finite(Yhat[, w]))
        if (length(cells) > 0L) Y.ct[cells, j] <- Yhat[cells, w]
    }

    list(
        eff = eff, Y.ct = Y.ct, Yo = Yo, obs = obs, Vhat = Vhat, Yhat = Yhat,
        D = Dm, T.on = Ton, plain = plain, lag.mask = !plain,
        anchor = anchor, n.post = n.post, ok = ok, rows = rows,
        att.avg = att.avg, att.avg.unit = att.avg.unit,
        time = time.on, att = att.on, count = count.on,
        not.rebuilt = ids[tr[!ok]]
    )
}

## ----------------------------------------------------------------------
## The outer fit: the inner fit with outcome-scale slots
## ----------------------------------------------------------------------
.fect_ar1_outer <- function(inner, rb, ctx, x) {
    out <- inner
    tr <- inner$tr
    co <- inner$co
    T1 <- nrow(rb$eff)
    Y.dat <- rb$Yo
    Y.dat[!rb$obs] <- 0
    out$Y.dat <- Y.dat
    obsm <- matrix(as.numeric(rb$obs), T1, ncol(rb$obs))
    out$I.dat <- obsm
    out$I <- obsm
    out$D.dat <- rb$D
    out$T.on <- rb$T.on
    out$eff <- rb$eff
    out$eff.tr <- rb$eff[, tr, drop = FALSE]
    out$Y.ct <- rb$Y.ct
    out$Y.tr.cnt <- rb$Y.ct[, tr, drop = FALSE]
    out$Y.ct.cnt <- rb$Y.ct[, co, drop = FALSE]
    out$Y.ct.full <- rb$Y.ct
    out$att.avg <- rb$att.avg
    out$att.avg.unit <- rb$att.avg.unit
    out$time <- rb$time
    out$att <- rb$att
    out$count <- rb$count
    ## the calendar-time slots have no outcome-scale counterpart here
    out["eff.calendar"] <- list(NULL)
    out["N.calendar"] <- list(NULL)
    out["eff.calendar.fit"] <- list(NULL)
    out["est.eff.calendar"] <- list(NULL)
    out["est.eff.calendar.fit"] <- list(NULL)
    ## no pre-trend tests (the pre-treatment cells are one-step-ahead errors)
    for (nm in c("test.out", "loo.test.out", "pre.est.att", "pre.att.bound",
                 "est.placebo", "est.carryover", "permute", "est.cm")) {
        out[nm] <- list(NULL)
    }
    out
}

## ----------------------------------------------------------------------
## Inference: the parametric bootstrap of the filtered panel, draws mapped
## through the recursion
## ----------------------------------------------------------------------
## The inner fit's parametric draws are V-scale treatment effects of pseudo
## treated units simulated with no effect (centred at zero). For a treated unit
## with anchor s, a draw's outcome-scale effect follows the same recursion as
## the point estimate: e_s = 0 and e_t = x e_{t-1} + effV_t, effV_t = 0 where V
## is unobserved (the draw has no error there); cells up to the anchor keep the
## draw's V-scale value (a one-step-ahead error). The recursion is linear, so
## the parametric H0 -> H1 location shift commutes with it. The draw statistics
## are then recomputed on the cells of the point fit, and the intervals are
## fect's own (normal, or basic with the parametric shift).
.fect_ar1_inference <- function(out, inner, rb, x, se, alpha, ci.method,
                                keep.sims) {
    ## the inner fit's V-scale inference slots never survive
    for (nm in c("est.avg", "est.avg.unit", "est.att", "est.att90", "att.bound",
                 "att.avg.boot", "att.avg.unit.boot", "att.boot",
                 "att.boot.original", "att.count.boot", "att.vcov", "eff.boot",
                 "D.boot", "I.boot", "colnames.boot")) {
        if (nm %in% names(out)) out[nm] <- list(NULL)
    }
    if (!isTRUE(se)) return(out)

    eb <- inner$eff.boot
    if (is.null(eb) || length(dim(eb)) != 3L) {
        stop("serial = \"ar1\": the parametric draws are not available.",
             call. = FALSE)
    }
    T1 <- dim(eb)[1]
    B <- dim(eb)[3]
    tr <- inner$tr
    Ntr <- length(tr)
    cn <- inner$colnames.boot
    same.tr <- vapply(cn, function(v) {
        identical(as.integer(v[seq_len(Ntr)]), as.integer(tr))
    }, TRUE)
    if (!all(same.tr)) {
        stop("serial = \"ar1\": the treated columns of the parametric draws ",
             "are not the fit's treated units.", call. = FALSE)
    }

    ## map each draw through the recursion
    Vobs <- inner$I.dat[, tr, drop = FALSE] == 1
    mapped <- eb
    for (w in seq_len(Ntr)) {
        j <- tr[w]
        if (!rb$ok[w]) {
            mapped[, w, ] <- NA_real_
            next
        }
        if (isTRUE(rb$plain)) {
            ## x = 0 without the lag mask: the draws are already on the
            ## outcome scale (V = Y)
            next
        }
        s <- rb$anchor[w]
        EV <- matrix(eb[, w, ], T1, B)
        EV0 <- EV
        EV0[!Vobs[, w], ] <- 0
        E <- matrix(NA_real_, T1, B)
        E[s, ] <- 0
        if (s < T1) {
            for (t in (s + 1L):T1) E[t, ] <- x * E[t - 1L, ] + EV0[t, ]
        }
        res <- matrix(NA_real_, T1, B)
        fin <- is.finite(rb$eff[, j])
        post <- which(seq_len(T1) > s & fin)
        pre <- which(seq_len(T1) <= s & fin)
        res[post, ] <- E[post, ]
        res[pre, ] <- EV[pre, ]
        mapped[, w, ] <- res
    }

    ## the draws' statistics, on the point fit's cells
    D.tr <- rb$D[, tr, drop = FALSE]
    d.v <- c(D.tr)
    t.on <- match(c(rb$T.on[, tr, drop = FALSE]), out$time)
    L <- length(out$time)
    att.avg.boot <- matrix(NA_real_, 1L, B)
    att.avg.unit.boot <- matrix(NA_real_, 1L, B)
    att.boot <- matrix(NA_real_, L, B)
    att.count.boot <- matrix(NA_real_, L, B)
    for (b in seq_len(B)) {
        eb.b <- matrix(mapped[, seq_len(Ntr), b], T1, Ntr)
        ev <- c(eb.b)
        fin <- !is.na(ev)
        if (any(fin)) {
            att.avg.boot[1L, b] <- sum(ev[fin] * d.v[fin]) / sum(d.v[fin])
        }
        unit <- vapply(seq_len(Ntr), function(w) {
            okc <- !is.na(eb.b[, w])
            dd <- sum(D.tr[okc, w])
            if (dd > 0) sum(eb.b[okc, w] * D.tr[okc, w]) / dd else NA_real_
        }, 0)
        if (any(!is.na(unit))) att.avg.unit.boot[1L, b] <- mean(unit, na.rm = TRUE)
        use <- fin & !is.na(t.on)
        if (any(use)) {
            cnt <- tabulate(t.on[use], nbins = L)
            sm <- numeric(L)
            agg <- rowsum(ev[use], t.on[use])
            sm[as.integer(rownames(agg))] <- agg[, 1]
            att.boot[, b] <- ifelse(cnt > 0, sm / pmax(cnt, 1), NA_real_)
            att.count.boot[, b] <- ifelse(cnt > 0, cnt, NA_real_)
        }
    }

    ## intervals, with fect's own rules (R/boot.R)
    basic <- identical(ci.method, "basic")
    pval.basic <- function(theta, boots) {
        one <- function(th, bb) {
            .pvalue_basic_dual(th, bb - mean(bb, na.rm = TRUE) + th)
        }
        if (is.matrix(boots)) {
            return(vapply(seq_len(nrow(boots)), function(k) {
                one(theta[k], boots[k, ])
            }, numeric(1)))
        }
        one(theta, boots)
    }
    z.half <- stats::qnorm(1 - alpha / 2)
    one_ci <- function(est, boots) {
        se1 <- stats::sd(boots, na.rm = TRUE)
        if (!basic) {
            list(se = se1, ci = c(est - se1 * z.half, est + se1 * z.half),
                 p = (1 - stats::pnorm(abs(est / se1))) * 2)
        } else {
            list(se = se1, ci = .basic_ci_shifted_one(est, boots, alpha, TRUE),
                 p = pval.basic(est, boots))
        }
    }
    a <- one_ci(rb$att.avg, att.avg.boot)
    est.avg <- t(as.matrix(c(rb$att.avg, a$se, a$ci, a$p)))
    colnames(est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")
    u <- one_ci(rb$att.avg.unit, att.avg.unit.boot)
    est.avg.unit <- t(as.matrix(c(rb$att.avg.unit, u$se, u$ci, u$p)))
    colnames(est.avg.unit) <- c("ATT.avg.unit", "S.E.", "CI.lower", "CI.upper",
                                "p.value")

    att <- out$att
    se.att <- apply(att.boot, 1, function(v) stats::sd(v, na.rm = TRUE))
    if (!basic) {
        CI.att <- cbind(att - se.att * z.half, att + se.att * z.half)
        p.att <- (1 - stats::pnorm(abs(att / se.att))) * 2
        att.bound <- cbind(att - se.att * stats::qnorm(1 - alpha),
                           att + se.att * stats::qnorm(1 - alpha))
    } else {
        CI.att <- .basic_ci_shifted(att, att.boot, alpha, TRUE)
        p.att <- pval.basic(att, att.boot)
        att.bound <- .basic_ci_shifted(att, att.boot, 2 * alpha, TRUE)
    }
    colnames(att.bound) <- c("CI.lower", "CI.upper")
    rownames(att.bound) <- out$time
    est.att <- cbind(att, se.att, CI.att, p.att, out$count)
    est.att90 <- cbind(att, se.att, att.bound, p.att, out$count)
    cols <- c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value", "count")
    colnames(est.att) <- cols
    colnames(est.att90) <- cols
    rownames(est.att) <- out$time
    rownames(est.att90) <- out$time
    vcov.att <- tryCatch(
        stats::cov(t(att.boot), use = "pairwise.complete.obs"),
        error = function(e) NA
    )

    out$est.avg <- est.avg
    out$att.avg.boot <- att.avg.boot
    out$est.avg.unit <- est.avg.unit
    out$att.avg.unit.boot <- att.avg.unit.boot
    out$est.att <- est.att
    out$est.att90 <- est.att90
    out$att.bound <- att.bound
    out$att.boot <- att.boot
    out$att.vcov <- vcov.att
    out$att.count.boot <- att.count.boot
    if (isTRUE(keep.sims)) {
        I.boot <- inner$I.boot
        for (w in seq_len(Ntr)) {
            I.boot[, w, ] <- as.numeric(rb$obs[, tr[w]])
        }
        out$eff.boot <- mapped
        out$D.boot <- inner$D.boot
        out$I.boot <- I.boot
        out$colnames.boot <- inner$colnames.boot
    }
    out
}

## ----------------------------------------------------------------------
## fit$ar1
## ----------------------------------------------------------------------
.fect_ar1_info <- function(x, x.fixed, inner, rb, path, converged, rounds,
                           profile, profile.rank, rule, rule.levels, prep, data,
                           id, min.T0, extra = list()) {
    ids <- as.character(inner$id)
    tr <- inner$tr
    anchor.time <- inner$rawtime[rb$anchor]
    rebuilt <- data.frame(
        id = ids[tr],
        anchor = anchor.time,
        n.post = rb$n.post,
        ok = rb$ok,
        stringsAsFactors = FALSE
    )
    all.ids <- unique(as.character(data[[id]]))
    gone <- setdiff(all.ids, ids)
    kept <- if (isFALSE(extra$lag.mask)) which(prep$pos >= 2L) else prep$keep
    in.filtered <- unique(as.character(prep$data[[id]][kept]))
    reason <- as.character(ifelse(
        gone %in% in.filtered,
        "removed by fect (min.T0 or max.missing on the filtered panel)",
        "no filtered period (fewer than two observed periods)"))
    dropped <- data.frame(id = gone, reason = reason, stringsAsFactors = FALSE)
    c(list(
        x = x, x.fixed = x.fixed, r = as.integer(inner$r.cv[1]),
        path = path, converged = converged, rounds = as.integer(rounds),
        profile = profile, profile.rank = profile.rank,
        rule = rule, rule.levels = rule.levels,
        rebuilt = rebuilt, not.rebuilt = rb$not.rebuilt, dropped = dropped,
        first.period = prep$tg[1], min.T0 = min.T0,
        pre.scale = "one-step-ahead", v.scale.slots = .FECT_AR1_V_SLOTS,
        x.bounds = c(0, 1)
    ), extra)
}

## ----------------------------------------------------------------------
## Features that cannot run on a serial = "ar1" fit
## ----------------------------------------------------------------------
## Returns silently for any other fit. Called at the top of each function that
## has no outcome-scale reading of a serial = "ar1" fit.
.fect_stop_ar1 <- function(fit, feature) {
    if (!is.list(fit) || is.null(fit[["serial"]]) ||
        !identical(fit[["serial"]], .FECT_SERIAL_AR1)) {
        return(invisible(NULL))
    }
    stop("serial = \"ar1\" fits do not support ", feature, ": the ",
         "pre-treatment cells are one-step-ahead prediction errors on the ",
         "outcome scale, not the residuals that feature reads.", call. = FALSE)
}
