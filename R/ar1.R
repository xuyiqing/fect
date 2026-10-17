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
    ## the three internal options, read once (never documented)
    opt.x <- getOption("fect.ar1.x", NULL)
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
    quiet <- function(expr) suppressWarnings(suppressMessages(expr))
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
    filter <- function(x) .fect_ar1_filter(prep, x, vars)
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
        if (is.null(f)) return(Inf)
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
        f <- function(x) {
            v <- ssr_at(x, rr)
            pts$x <- c(pts$x, x)
            pts$ssr <- c(pts$ssr, v)
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
        list(x = xx, profile = data.frame(x = pts$x[o.ord], ssr = pts$ssr[o.ord]))
    }

    ## ---- the treated-unit rule on a panel ----------------------------------
    r.range <- spec$r.range
    rule_on <- function(fd, min.T0.rule) {
        f <- tryCatch(
            quiet(run(fd, r = r.range, CV = TRUE, cv.method = .FECT_CV_PRETREATMENT,
                      se = FALSE, parallel = FALSE, min.T0 = min.T0.rule)),
            error = fail
        )
        list(r = as.integer(f$r.cv[1]), rule = f$cv.treated)
    }

    ## ---- choose (x, r) --------------------------------------------------------
    path <- data.frame(round = integer(0), x = numeric(0), r = integer(0),
                       att = numeric(0))
    rule.levels <- NULL
    rule.last <- NULL
    profile <- data.frame(x = numeric(0), ssr = numeric(0))
    profile.rank <- NA_integer_
    converged <- NA
    x.fixed <- !is.null(opt.x)
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
        path <- rbind(path, data.frame(round = 1L, x = x.k, r = as.integer(r.k),
                                       att = NA_real_))
        rounds <- 1L
        message("serial = \"ar1\": x = ", format(round(x.k, 4)), ", rank ", r.k)
    } else {
        r0 <- rule_on(fd.levels, min.T0)
        rule.levels <- r0$rule
        path <- rbind(path, data.frame(round = 0L, x = NA_real_, r = r0$r,
                                       att = NA_real_))
        r.prev <- r0$r
        rounds <- 0L
        converged <- FALSE
        for (k in seq_len(K)) {
            if (x.fixed) {
                x.k <- as.numeric(opt.x)
            } else {
                xr <- x_of_r(r.prev)
                x.k <- xr$x
                profile <- xr$profile
                profile.rank <- as.integer(r.prev)
            }
            rl <- rule_on(filter(x.k), min.rule)
            r.k <- rl$r
            rule.last <- rl$rule
            rounds <- k
            path <- rbind(path, data.frame(round = k, x = x.k, r = r.k,
                                           att = NA_real_))
            message("serial = \"ar1\": round ", k, ": x = ",
                    format(round(x.k, 4)), ", rank ", r.k)
            if (x.fixed) {
                converged <- NA
                break
            }
            if (r.k == r.prev) {
                converged <- TRUE
                break
            }
            r.prev <- r.k
        }
        if (isFALSE(converged)) {
            if (is.null(opt.rounds)) {
                warning("serial = \"ar1\": the rank and x did not settle in ",
                        K, " rounds; using the last pair (x = ",
                        format(round(x.k, 4)), ", rank ", r.k, ").",
                        call. = FALSE)
            }
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
        suppressMessages(
            if (se) {
                run(filter(x.k), r = r.k, CV = FALSE, se = TRUE,
                    vartype = "parametric", nboots = nboots,
                    para.error = para.error, alpha = alpha, ci.method = ci.method,
                    parallel = parallel, cores = cores, seed = seed,
                    keep.sims = TRUE, min.T0 = min.final)
            } else {
                run(filter(x.k), r = r.k, CV = FALSE, se = FALSE,
                    parallel = FALSE, min.T0 = min.final)
            }
        ),
        error = fail
    )

    ctx <- list(prep = prep, id = id, time = time, Y = Y, D = D, X = X,
                has_a = has_a, has_xi = has_xi, force.code = force.code)
    rb <- .fect_ar1_rebuild(inner, ctx, x.k)
    if (length(rb$not.rebuilt) > 0L) {
        warning("serial = \"ar1\": ", length(rb$not.rebuilt),
                " treated unit(s) could not be rebuilt on the outcome scale ",
                "and are left out: ",
                paste(utils::head(rb$not.rebuilt, 10), collapse = ", "),
                if (length(rb$not.rebuilt) > 10) ", ..." else "",
                call. = FALSE)
    }

    ## the outcome-scale ATT at each pair of the path (cheap: one fit each)
    for (i in seq_len(nrow(path))) {
        if (is.na(path$x[i])) next
        if (i == nrow(path)) {
            path$att[i] <- rb$att.avg
            next
        }
        path$att[i] <- tryCatch({
            fi <- quiet(run(filter(path$x[i]), r = path$r[i], CV = FALSE,
                            se = FALSE, parallel = FALSE,
                            min.T0 = min.fit(path$r[i])))
            .fect_ar1_rebuild(fi, ctx, path$x[i])$att.avg
        }, error = function(e) NA_real_)
    }

    out <- .fect_ar1_outer(inner, rb, ctx, x.k)
    out <- .fect_ar1_inference(out, inner, rb, x.k, se = se, alpha = alpha,
                               ci.method = ci.method, keep.sims = keep.sims)
    out$ar1 <- .fect_ar1_info(
        x = x.k, x.fixed = x.fixed, inner = inner, rb = rb, path = path,
        converged = converged, rounds = rounds, profile = profile,
        profile.rank = profile.rank, rule = rule.last, rule.levels = rule.levels,
        prep = prep, data = data, id = id, min.T0 = min.final
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
.fect_ar1_filter <- function(prep, x, vars) {
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
        Xf[, , k] <- Xfull[rows, , k] - x * Xfull[rows - 1L, , k]
    }

    tr <- inner$tr
    Dm <- inner$D.dat
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

    ## effects
    eff <- inner$eff
    Dtr <- Dm[, tr, drop = FALSE]
    anchor <- rep(NA_integer_, length(tr))
    n.post <- integer(length(tr))
    Yhat <- matrix(NA_real_, T1, length(tr))
    for (w in seq_along(tr)) {
        j <- tr[w]
        pre.ok <- which(Dm[, j] == 0 & obs[, j])
        if (length(pre.ok) == 0L) next
        s <- max(pre.ok)
        anchor[w] <- s
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
    Ton <- inner$T.on
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
        cells <- which(!is.na(anchor[w]) & ok[w] & seq_len(T1) > anchor[w] &
                           !obs[, j] & is.finite(Yhat[, w]))
        if (length(cells) > 0L) Y.ct[cells, j] <- Yhat[cells, w]
    }

    list(
        eff = eff, Y.ct = Y.ct, Yo = Yo, obs = obs, Vhat = Vhat, Yhat = Yhat,
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

## placeholder until the inference block is added
.fect_ar1_inference <- function(out, inner, rb, x, se, alpha, ci.method,
                                keep.sims) {
    if (isTRUE(se)) {
        stop("serial = \"ar1\" inference is not available in this build yet.",
             call. = FALSE)
    }
    out
}

## ----------------------------------------------------------------------
## fit$ar1
## ----------------------------------------------------------------------
.fect_ar1_info <- function(x, x.fixed, inner, rb, path, converged, rounds,
                           profile, profile.rank, rule, rule.levels, prep, data,
                           id, min.T0) {
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
    in.filtered <- unique(as.character(prep$data[[id]][prep$keep]))
    reason <- as.character(ifelse(
        gone %in% in.filtered,
        "removed by fect (min.T0 or max.missing on the filtered panel)",
        "no filtered period (fewer than two observed periods)"))
    dropped <- data.frame(id = gone, reason = reason, stringsAsFactors = FALSE)
    list(
        x = x, x.fixed = x.fixed, r = as.integer(inner$r.cv[1]),
        path = path, converged = converged, rounds = as.integer(rounds),
        profile = profile, profile.rank = profile.rank,
        rule = rule, rule.levels = rule.levels,
        rebuilt = rebuilt, not.rebuilt = rb$not.rebuilt, dropped = dropped,
        first.period = prep$tg[1], min.T0 = min.T0,
        pre.scale = "one-step-ahead", v.scale.slots = .FECT_AR1_V_SLOTS,
        x.bounds = c(0, 1)
    )
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
