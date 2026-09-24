###############################################################################
## Double leave-one-out (dloo) pre-trend placebos computed directly from
## the default in-sample (fect-default) imputation fit.
##
## These are INTERNAL helpers used by fect() when it is called with
## `dloo = TRUE` (and, optionally, `dloo.adjust = TRUE`). They are not
## exported. The user-facing surface is the pair of fect() flags, which mirror
## the existing `loo` flag: like `loo`, `dloo` fills the object's
## `pre.est.att` / `pre.att.bound` / `pre.att.boot` slots with corrected
## pre-treatment placebo estimates. Unlike `loo`, it does so WITHOUT re-fitting
## the imputation model --- everything is recovered in closed form from the
## in-sample fit (its `eff` matrix and the cohort sizes).
##
## Background (Li and Strezhnev, 2026, "Benchmarking parallel trends violations in
## regression imputation difference-in-differences").
##   fect's default pre-trend placebo is the in-sample TWFE residual
##   ATT^is_g(t) = Ybar[g,t] - Yhat0[g,t] at a pre-treatment control cell
##   (cohort g, period t < g). This carries an
##   *attenuation* bias (from structurally-zero / own-cohort comparisons) and a
##   *contamination* bias (under staggered adoption, from using early-adopting cohorts).
##
##   The double leave-one-out (dloo) estimator removes
##   both: for cohort g and pre-period t it drops the tested period t and
##   keeps only strictly-later-adopting (and never-treated) units as controls,
##
##     ATT^dloo_g(t) = 1 / ((g-2) N^(g)) *
##         sum_{g' > g} sum_{t' != t, t' < g} N_g'
##             ( Ybar[g,t] - Ybar[g,t'] - Ybar[g',t] + Ybar[g',t'] ),      (D)
##
##   where N^(g) = # units in cohorts adopting strictly after g, and g-1 is the
##   number of pre-treatment periods for cohort g (balanced panel).
##
##   Li and Strezhnev (2026) shows (D) is a fixed
##   LINEAR function of the per-cohort in-sample placebos that fect already
##   computes: because the unit/time fixed effects cancel in a DiD taken over
##   control cells, the 2x2 DiD of group-time means equals the 2x2 DiD of the
##   in-sample placebos, so no per-cell re-estimation is required. Equivalently:
##
##     ATT^dloo_g(t) = (g-1)/(g-2) * ( ATT^is_g(t)
##                       - (1/N^(g)) sum_{g' > g} N_g' ATT^is_g'(t) ).
##
##   We implement the algebraically identical form actually used below: an
##   average of the tested period t against the retained baseline periods,
##   written directly in terms of the in-sample placebos ATT^is_g(.). Each
##   bracketed quantity is the placebo at one period, decontaminated by the
##   N-weighted mean of the later-adopting controls (cohorts g' > g):
##
##     ATT^dloo_g(t) =
##         [ ATT^is_g(t)  - (1/N^(g)) sum_{g' > g} N_g' ATT^is_g'(t)  ]
##       - mean_{t' in B} [ ATT^is_g(t') - (1/N^(g)) sum_{g' > g} N_g' ATT^is_g'(t') ],
##
##   where B is the set of retained pre-treatment baseline periods of g (t' < g):
##   for dloo, B excludes the tested period t; for dloo.adjust, B keeps it (see
##   below).
##
##   Two reasons we average this way rather than applying the (g-1)/(g-2) closed
##   form directly: (i) no rescaling constant is ever formed, so it degrades
##   gracefully when a cohort has few baselines (a short B) instead of dividing
##   by g-2; and (ii) switching between dloo and dloo.adjust is just a change to
##   which periods are in B. This reproduces (D) exactly: because the in-sample
##   pre-period placebos sum to zero over t' < g, subtracting the leave-one-out
##   baseline mean rescales the tested-period term by (g-1)/(g-2), recovering the
##   closed form above.
##
## The dloo.adjust correction (Liu's "pre-treatment average" baseline).
## (Liu, "Cohort-Anchored Robust Inference for Event-Study with Staggered Adoption")
##   The *difference* of two dloo estimates over-states the period-to-period
##   change in the parallel-trends violation by the factor (g-1)/(g-2) (no
##   common baseline). Liu (2025) / Roth (2024) instead average over ALL
##   pre-treatment baselines including t' = t; the t' = t term is a mechanically
##   zero DiD, so this is exactly
##
##     ATT^avg_g(t) = ((g-2)/(g-1)) * ATT^dloo_g(t),
##
##   i.e. dloo with the tested period KEPT in the baseline set B. That is what
##   `dloo.adjust = TRUE` selects here --- a minimal change to which in-sample
##   estimates enter the aggregation (no rescaling constant is ever formed).
##   See the paper for the equality ATT^avg_g(t) - ATT^avg_g(t*) = the 2x2 DiD.
##
## Exactness. The overlay matches the literal "re-fit the imputation model on
## the restricted control pool per (cohort, period)" dloo to machine precision
## for a balanced control panel and additive two-way FE (method = "fe"). With
## latent factors / matrix completion the FE do not cleanly cancel; fect()
## stops before calling in that case.
##
## Covariates. In the DiD above the unit and time effects cancel, but X beta
## does not: X differs across the four cells. fit$eff carries the full-sample
## beta, which fect estimates from EVERY untreated cell, including control
## units' cells in periods after cohort g has adopted. Taken over fit$eff, the
## placebo would therefore move when a late control cell moves (found
## 2026-07-14: +50 on one never-treated cell in the last period shifted the
## event-time -1 placebo by 0.09). So with covariates we do what the re-fit
## does: for each (cohort g, tested period t) we estimate beta_(g,t) by
## two-way FE OLS on the same pool the re-fit uses (cohort g at its
## pre-periods other than t; the later-adopting / never-treated units at all
## periods before g), and take the DiD of Y - X beta_(g,t). The DiD of the
## adjusted outcome equals the re-fit's placebo exactly (the FE part is the
## no-covariate identity applied to Y - X beta_(g,t)). beta_(g,t) has a closed
## form: the pool is balanced except for one block (cohort g at period t);
## filling that block with its own fitted values makes the panel balanced
## without changing the fit (.dloo_fill_block), after which the usual two-way
## within transformation applies. dloo.adjust uses the same beta_(g,t), so the
## (g-2)/(g-1) relation to dloo holds with covariates too. Without covariates
## the code below is unchanged (it reads fit$eff).
###############################################################################


## ---------------------------------------------------------------------------
## Validate the (balanced, staggered) panel a fect fit must satisfy for the
## closed-form dloo overlay to be exact. Called both early in fect() (fail
## fast, before the expensive fit) and defensively here.
##
## Requirements, per Li & Strezhnev:
##   * staggered adoption / no reversal   -> checked via hasRevs upstream
##   * balanced panel with NO missing PRE-treatment cells (post-treatment
##     missingness is irrelevant to a pre-trend placebo, so it is allowed).
##
## D, I are TT x N matrices (treatment, inclusion). Returns TRUE invisibly if
## the pre-treatment panel is complete; otherwise a character string naming the
## violation (so the caller can stop() with a specific message).
## ---------------------------------------------------------------------------
.dloo_check_prebalance <- function(D, I) {
    TT <- nrow(D)
    N  <- ncol(D)
    ## cumulative treatment -> first-treated row per unit (Inf if never)
    for (i in seq_len(N)) {
        d <- D[, i]
        d[is.na(d)] <- 0
        wt <- which(d == 1)
        first <- if (length(wt) == 0) TT + 1L else min(wt)
        if (first <= 1L) next               # treated from period 1: no pre cells
        pre_rows <- seq_len(first - 1L)
        obs <- I[pre_rows, i]
        ## a pre-treatment cell is "missing" when it is not included (I != 1)
        if (any(is.na(obs)) || any(obs != 1)) {
            return("unbalanced pre-treatment panel (a unit has a missing pre-treatment period)")
        }
    }
    invisible(TRUE)
}


## ---------------------------------------------------------------------------
## Pre-compute the (weight-independent) panel structure from a fitted, in-sample
## fect object: cohort labels, per-unit cohort index, the control-cell mask, the
## set of estimable (cohort, period) cells, and per-cohort event times.
## Mirrors the slots the object already carries (eff / D.dat / I.dat / rawtime).
##
## With covariates, pass the outcome and covariates the fit used: `Y` (TT x N)
## and `X` (TT x N x p), aligned with fit$eff's columns, and `scale` (fect's
## normalization constant, 1 if none) to return them to the units of fit$eff.
## The per-cell coefficients beta_(g,t) are then computed here (see header).
## ---------------------------------------------------------------------------
.dloo_prepare <- function(fit, controls = "not-yet-treated",
                          Y = NULL, X = NULL, scale = 1) {

    eff     <- fit$eff
    D       <- fit$D.dat
    I       <- fit$I.dat
    TT      <- nrow(eff)
    N       <- ncol(eff)

    ## Work in ROW-INDEX time, not calendar time. fect's time axis is unique and
    ## sorted (default.R builds tname via unique(sort(.))), so ranks are an
    ## order-preserving relabelling: every use below is a comparison (which(t < g)),
    ## a lookup (match(g, .)) or a difference (tcol - gcol + 1), all invariant to
    ## it. Two reasons to normalise rather than carry fit$rawtime:
    ##   * fect accepts Date and character time indices. vapply(., numeric(1))
    ##     below strips a Date's class, after which match(<numeric>, <Date>)
    ##     compares "12782" to "2000-12-31" and returns NA for every cell --- the
    ##     point estimates come back silently all-NA while the bootstrap path
    ##     (which already ranks) still fills the S.E. column.
    ##   * it makes this path byte-identical to .dloo_boot_draw's, which passes
    ##     seq_len(TT). The two must never disagree.
    ## The original calendar labels are kept for cell_cohort_label only.
    rawtime_orig <- fit$rawtime
    rawtime      <- seq_len(TT)

    ## Control (= not-treated, observed) cell mask.
    ctrl_mask <- (!is.na(I) & I == 1) & (!is.na(D) & D == 0)

    ## First-treated row per unit; Inf for never-treated.
    first_treat <- vapply(seq_len(N), function(i) {
        d <- D[, i]
        w <- which(!is.na(d) & d == 1)
        if (length(w) == 0) Inf else rawtime[min(w)]
    }, numeric(1))

    cohorts    <- sort(unique(first_treat))            # incl. Inf if present
    n_cohorts  <- length(cohorts)
    cohort_idx <- match(first_treat, cohorts)          # per-unit index
    never_idx  <- which(!is.finite(cohorts))

    ## Enumerate (cohort, period) cells for treated cohorts with at least two
    ## pre-periods (the minimum for a dloo baseline). Cohorts with no
    ## later-adopting control pool are kept but evaluate to NA.
    ## NB: this >= 2 requirement is enforced here, so it applies to dloo.adjust
    ## too --- a single-pre-period cohort is excluded outright, NOT entered as
    ## the mechanically-zero placebo the dloo.adjust note below might suggest.
    ## fect() fails fast upstream if NO cohort clears this.
    treated_pos <- which(is.finite(cohorts))
    cohort_v <- integer(0); tcol_v <- integer(0)
    gcol_v <- integer(0); et_v <- numeric(0)
    for (gi in treated_pos) {
        g    <- cohorts[gi]
        gcol <- match(g, rawtime)
        pre  <- which(rawtime < g)
        if (length(pre) < 2) next
        for (tcol in pre) {
            cohort_v <- c(cohort_v, gi)
            tcol_v   <- c(tcol_v, tcol)
            gcol_v   <- c(gcol_v, gcol)
            et_v     <- c(et_v, tcol - gcol + 1)        # fect event time (T.on)
        }
    }

    ## Optional subgroup breakdown. fit$G is the TT x N group-code matrix
    ## (constant within unit, aligned to eff's columns); fit$group maps codes to
    ## raw labels. We take one code per unit and index it into the sorted set of
    ## codes so the overlay can partition the (linear) per-unit dloo placebo by
    ## subgroup --- shared controls, exactly mirroring how fect's own group-wise
    ## ATT partitions the treated units of a single imputation fit.
    unit_group <- NULL; group_labels <- NULL
    if (!is.null(fit$G) && !is.null(fit$group)) {
        G  <- fit$G
        ug <- vapply(seq_len(N), function(i) {
            u <- unique(G[, i][!is.na(G[, i])])
            if (length(u) == 0) NA_real_ else u[1]
        }, numeric(1))
        map    <- fit$group                            # cols: rawgroup, newgroup
        codes  <- sort(unique(map$newgroup))
        group_labels <- as.character(map$rawgroup[match(codes, map$newgroup)])
        unit_group   <- match(ug, codes)               # per-unit index; NA if none
    }

    ## Covariate adjustment (see header): outcome and covariates in the units
    ## of fit$eff, plus one coefficient vector per (cohort, tested period) cell.
    Ycov <- Xcov <- cell_beta <- NULL
    if (!is.null(X) && length(dim(X)) == 3 && dim(X)[3] > 0) {
        if (is.null(Y) || !identical(dim(as.matrix(Y)), c(TT, N)) ||
            !identical(dim(X)[1:2], c(TT, N))) {
            stop("Internal error: the outcome and covariates passed to the ",
                 "dloo overlay do not match the fitted panel.", call. = FALSE)
        }
        Ycov <- as.matrix(Y) * scale
        Xcov <- X * scale
        cell_beta <- .dloo_cell_betas(Ycov, Xcov, cohorts, cohort_idx,
                                      never_idx, controls, rawtime,
                                      cohort_v, tcol_v)
    }

    list(eff = eff, ctrl_mask = ctrl_mask, rawtime = rawtime,
         TT = TT, N = N, cohorts = cohorts, n_cohorts = n_cohorts,
         cohort_idx = cohort_idx, never_idx = never_idx,
         controls = controls,
         cell_cohort = cohort_v, cell_tcol = tcol_v,
         cell_event = et_v,
         cell_cohort_label = ifelse(is.finite(cohorts[cohort_v]),
                                    rawtime_orig[cohorts[cohort_v]], NA),
         unit_group = unit_group, group_labels = group_labels,
         Ycov = Ycov, Xcov = Xcov, cell_beta = cell_beta)
}


## ---------------------------------------------------------------------------
## Fill the one missing block of a balanced pool. `Z` is T_p x n (rows = the
## cohort's pre-periods, cols = pool units), `gmask` marks the tested cohort's
## columns and `tpos` the tested period's row. The tested cells are excluded
## from the fit; this replaces them by their two-way fitted values from the
## fit on the remaining cells. Those cells then have zero residual, so the
## two-way projection of the filled (balanced) panel equals the projection on
## the pool without them. The fill value solves
##   z_j = zbar_j. + zbar_.t - zbar_..   (means over the filled panel),
## which is linear in the unknowns; summing over the tested units first gives
## their total u in closed form, then each z_j.
## ---------------------------------------------------------------------------
.dloo_fill_block <- function(Z, gmask, tpos) {
    Tp <- nrow(Z)
    n  <- ncol(Z)
    nG <- sum(gmask)
    nL <- n - nG
    S_j <- colSums(Z[-tpos, gmask, drop = FALSE])   # tested units, other periods
    A   <- sum(Z[tpos, !gmask])                      # pool units at the tested period
    C   <- sum(Z[, !gmask]) + sum(S_j)               # total over the pool
    u   <- (sum(S_j) / Tp + nG * A / n - nG * C / (n * Tp)) /
        ((1 - 1 / Tp) * nL / n)
    Z[tpos, gmask] <- (S_j / Tp + (A + u) / n - (C + u) / (n * Tp)) /
        (1 - 1 / Tp)
    Z
}


## Two-way within transformation of a balanced T x n matrix.
.dloo_twoway_resid <- function(Z) {
    Z - outer(rowMeans(Z), rep(1, ncol(Z))) -
        outer(rep(1, nrow(Z)), colMeans(Z)) + mean(Z)
}


## ---------------------------------------------------------------------------
## beta_(g,t) for every enumerated cell: two-way FE OLS of Y on X over the
## pool of cohort g (its units at pre-periods other than t, plus the
## later-adopting / never-treated units at all pre-periods of g). Returns an
## ncell x p matrix; NA rows for cells without a control pool (those cells
## are NA in the overlay anyway). A covariate with no within variation in the
## pool (e.g. time-invariant) gets coefficient 0; its DiD is 0 too.
## ---------------------------------------------------------------------------
.dloo_cell_betas <- function(Ycov, Xcov, cohorts, cohort_idx, never_idx,
                             controls, rawtime, cell_cohort, cell_tcol) {
    p     <- dim(Xcov)[3]
    ncell <- length(cell_cohort)
    out   <- matrix(NA_real_, ncell, p)
    for (gi in unique(cell_cohort)) {
        g <- cohorts[gi]
        later_idx <- if (controls == "never-treated") never_idx else
                     which(cohorts > g)
        cols_G <- which(cohort_idx == gi)
        cols_L <- which(cohort_idx %in% later_idx)
        if (length(cols_L) == 0 || length(cols_G) == 0) next
        pre   <- which(rawtime < g)
        cols  <- c(cols_G, cols_L)
        gmask <- c(rep(TRUE, length(cols_G)), rep(FALSE, length(cols_L)))
        Yp    <- Ycov[pre, cols, drop = FALSE]
        Xp    <- lapply(seq_len(p), function(m) Xcov[pre, cols, m, drop = TRUE])
        Xp    <- lapply(Xp, function(x) matrix(x, length(pre), length(cols)))
        for (k in which(cell_cohort == gi)) {
            tpos <- match(cell_tcol[k], pre)
            ry <- c(.dloo_twoway_resid(.dloo_fill_block(Yp, gmask, tpos)))
            RX <- vapply(Xp, function(x)
                c(.dloo_twoway_resid(.dloo_fill_block(x, gmask, tpos))),
                numeric(length(ry)))
            RX <- matrix(RX, ncol = p)
            ## drop covariates without within variation in this pool
            raw_norm <- vapply(Xp, function(x) sqrt(sum(x^2)), numeric(1))
            ok <- sqrt(colSums(RX^2)) > 1e-10 * pmax(raw_norm, 1)
            b  <- rep(0, p)
            if (any(ok)) {
                bk <- qr.coef(qr(RX[, ok, drop = FALSE]), ry)
                bk[is.na(bk)] <- 0
                b[ok] <- bk
            }
            out[k, ] <- b
        }
    }
    out
}


## ---------------------------------------------------------------------------
## Evaluate the overlay for a vector of per-unit weights `w`. Returns the
## per-cell placebos and the size-weighted event-study series. `correct = TRUE`
## keeps the tested period in the baseline set B (Liu's pre-treatment-average /
## dloo.adjust); `correct = FALSE` drops it (the conventional dloo). Used both
## for the point estimate (w == 1) and for each bootstrap draw.
## ---------------------------------------------------------------------------
.dloo_eval <- function(prep, w, correct) {

    eff       <- prep$eff
    ctrl      <- prep$ctrl_mask
    cohorts   <- prep$cohorts
    n_cohorts <- prep$n_cohorts
    cohort_idx <- prep$cohort_idx
    rawtime   <- prep$rawtime
    TT        <- prep$TT

    ## With covariates the placebo is built from Y - X beta_(g,t) (header);
    ## without, from fit$eff. beta_(g,t) is unweighted, and every caller passes
    ## unit weights of 1.
    has_cov <- !is.null(prep$cell_beta)
    if (has_cov && any(w != 1)) {
        stop("Internal error: the covariate-adjusted dloo overlay assumes ",
             "unit weights of 1.", call. = FALSE)
    }
    mats <- if (has_cov) {
        c(list(prep$Ycov),
          lapply(seq_len(dim(prep$Xcov)[3]),
                 function(m) matrix(prep$Xcov[, , m], TT, prep$N)))
    } else {
        list(eff)
    }

    ## Weighted per-cohort x period control means M[ci, s] and weight totals
    ## Wt[ci, s] (for pooling later-adopter controls), for each matrix in
    ## `mats` (the weights, hence Wt, are the same for all of them).
    cohort_means <- function(mat) {
        M  <- matrix(NA_real_, n_cohorts, TT)
        Wt <- matrix(0, n_cohorts, TT)
        for (ci in seq_len(n_cohorts)) {
            cols <- which(cohort_idx == ci)
            if (length(cols) == 0) next
            e_sub <- mat[, cols, drop = FALSE]
            m_sub <- ctrl[, cols, drop = FALSE]
            wv    <- w[cols]
            wmat  <- matrix(wv, nrow = TT, ncol = length(cols), byrow = TRUE)
            wmat[!m_sub] <- 0
            denom <- rowSums(wmat)
            num   <- rowSums(wmat * ifelse(m_sub, e_sub, 0))
            ok    <- denom > 0
            M[ci, ok] <- num[ok] / denom[ok]
            Wt[ci, ]  <- denom
        }
        list(M = M, Wt = Wt)
    }
    cm <- lapply(mats, cohort_means)
    M  <- cm[[1]]$M
    Wt <- cm[[1]]$Wt

    ## Per-cohort size under weights w (for the size-weighted aggregate).
    wsize <- vapply(seq_len(n_cohorts),
                    function(ci) sum(w[cohort_idx == ci]), numeric(1))

    ncell <- length(prep$cell_cohort)
    val   <- rep(NA_real_, ncell)

    ## Later-adopter control mean ebar depends on the target cohort gi only
    ## through the "later" set; cache once per distinct gi (one vector per
    ## matrix in `mats`).
    uniq_gi <- unique(prep$cell_cohort)
    ebar_cache <- vector("list", n_cohorts)
    for (gi in uniq_gi) {
        g <- cohorts[gi]
        later_idx <- if (prep$controls == "never-treated") prep$never_idx else
                     which(cohorts > g)
        Wl <- Wt[later_idx, , drop = FALSE]
        ebar_cache[[gi]] <- lapply(cm, function(cmi) {
            Ml  <- cmi$M[later_idx, , drop = FALSE]
            num <- colSums(Wl * ifelse(is.na(Ml), 0, Ml))
            den <- colSums(Wl)
            ifelse(den > 0, num / den, NA_real_)
        })
    }

    ## Later-adopter mean for cell k: of fit$eff, or of Y - X beta_(g,t).
    cell_ebar <- function(k, gi) {
        eb <- ebar_cache[[gi]][[1]]
        if (has_cov) {
            b <- prep$cell_beta[k, ]
            for (m in seq_along(b)) eb <- eb - b[m] * ebar_cache[[gi]][[m + 1]]
        }
        eb
    }

    for (k in seq_len(ncell)) {
        gi   <- prep$cell_cohort[k]
        tcol <- prep$cell_tcol[k]
        g    <- cohorts[gi]
        pre  <- which(rawtime < g)
        if (has_cov && anyNA(prep$cell_beta[k, ])) next
        eb   <- cell_ebar(k, gi)
        Mg   <- M[gi, ]
        if (has_cov) {
            b <- prep$cell_beta[k, ]
            for (m in seq_along(b)) Mg <- Mg - b[m] * cm[[m + 1]]$M[gi, ]
        }
        Dvec <- Mg - eb                            # decontaminated placebo, all periods
        base <- if (correct) pre else setdiff(pre, tcol)
        if (length(base) == 0) next
        d_t  <- Dvec[tcol]
        d_b  <- Dvec[base]
        if (is.na(d_t) || any(is.na(d_b))) next
        val[k] <- d_t - mean(d_b)
    }

    ## Size-weighted aggregate over cohorts at each pre-treatment event time.
    ets <- sort(unique(prep$cell_event))
    att <- vapply(ets, function(et) {
        sel <- prep$cell_event == et & !is.na(val)
        if (!any(sel)) return(NA_real_)
        wk <- wsize[prep$cell_cohort[sel]]
        sum(wk * val[sel]) / sum(wk)
    }, numeric(1))
    n.units <- vapply(ets, function(et) {
        sel <- prep$cell_event == et & !is.na(val)
        if (!any(sel)) return(0)
        sum(wsize[prep$cell_cohort[sel]])
    }, numeric(1))

    ## Optional subgroup breakdown. The dloo cell value for cohort gi at period
    ## tcol is the cohort-average of a per-unit placebo
    ##   val_j = (eff[tcol,j] - eb[tcol]) - mean_{t' in B}(eff[t',j] - eb[t']),
    ## with eb the *shared* later-adopter control mean. Averaging val_j over a
    ## subgroup's units (weights w) gives that subgroup's dloo series; summing
    ## the weighted numerators over all subgroups reproduces the pooled series
    ## exactly (sum_{j in gi} w_j val_j = wsize[gi] * val[cell]).
    group_att <- NULL; group_count <- NULL
    if (!is.null(prep$unit_group)) {
        ng   <- length(prep$group_labels)
        nets <- length(ets)
        gnum <- matrix(0, ng, nets); gden <- matrix(0, ng, nets)
        gcnt <- matrix(0, ng, nets)
        for (k in seq_len(ncell)) {
            gi   <- prep$cell_cohort[k]
            tcol <- prep$cell_tcol[k]
            g    <- cohorts[gi]
            pre  <- which(rawtime < g)
            if (has_cov && anyNA(prep$cell_beta[k, ])) next
            eb   <- cell_ebar(k, gi)
            base <- if (correct) pre else setdiff(pre, tcol)
            if (length(base) == 0) next
            if (is.na(eb[tcol]) || any(is.na(eb[base]))) next
            eti  <- match(prep$cell_event[k], ets)
            for (j in which(cohort_idx == gi)) {
                grp <- prep$unit_group[j]
                if (is.na(grp)) next
                yj <- eff[, j]
                if (has_cov) {
                    b  <- prep$cell_beta[k, ]
                    yj <- prep$Ycov[, j]
                    for (m in seq_along(b)) yj <- yj - b[m] * prep$Xcov[, j, m]
                }
                dj_t <- yj[tcol] - eb[tcol]
                dj_b <- yj[base] - eb[base]
                if (is.na(dj_t) || any(is.na(dj_b))) next
                vj <- dj_t - mean(dj_b)
                gnum[grp, eti] <- gnum[grp, eti] + w[j] * vj
                gden[grp, eti] <- gden[grp, eti] + w[j]
                gcnt[grp, eti] <- gcnt[grp, eti] + 1
            }
        }
        group_att   <- ifelse(gden > 0, gnum / gden, NA_real_)
        group_count <- gcnt
        rownames(group_att) <- rownames(group_count) <- prep$group_labels
        colnames(group_att) <- colnames(group_count) <- ets
    }

    list(event.time = ets, att = att, n.units = n.units,
         group_att = group_att, group_count = group_count,
         group_labels = prep$group_labels)
}


## ---------------------------------------------------------------------------
## Per-replicate dloo draw, evaluated INSIDE fect's own bootstrap loop
## (`fect_boot`). Given one case-resampled replicate's panel --- eff / D / I on
## the resampled columns, plus the resampled group codes `G` --- it returns the
## pooled pre-trend vector aligned to `pre.term` and, when grouped, the per-group
## matrix (rows = sorted group codes, cols = pre.term). Because fect_boot calls
## this per replicate and keeps only the returned (small) vectors, inference
## rides the native resampling with no separate resampler and without retaining
## the full TT x N x nboots panels (no keep.sims needed). With covariates the
## replicate's Y and X are passed too, so beta_(g,t) is re-estimated on the
## replicate's own pools.
## ---------------------------------------------------------------------------
.dloo_boot_draw <- function(eff, D, I, rawtime, group.map, G,
                            controls, correct, pre.term,
                            Y = NULL, X = NULL, scale = 1) {
    pseudo <- list(eff = eff, D.dat = D, I.dat = I, rawtime = rawtime,
                   group = group.map, G = G)
    prep <- .dloo_prepare(pseudo, controls, Y = Y, X = X, scale = scale)
    ev   <- .dloo_eval(prep, rep(1, prep$N), correct)
    pos  <- match(pre.term, ev$event.time)
    list(att       = ev$att[pos],
         group_att = if (!is.null(ev$group_att))
                         ev$group_att[, pos, drop = FALSE] else NULL)
}


## ---------------------------------------------------------------------------
## Top-level: build the fect-format pre-treatment placebo slots (indexed by the
## object's own pre-treatment event times `fit$time[fit$time <= 0]`) from the
## dloo overlay. The point estimate is the deterministic overlay on the fitted
## object; inference is taken from fect's OWN bootstrap: `fect_boot` has already
## applied `.dloo_boot_draw` to each case-resampled replicate and handed back the
## draws (`boot.pre`, and per group `boot.pre.group`). We only turn those draws
## into SEs/CIs here, so the slots are drop-in compatible with the existing `loo`
## machinery (plot, print, diagtest).
##
##   pre.est.att  : matrix, rows = pre.term, cols
##                  ATT, S.E., CI.lower, CI.upper, p.value, count.on
##   pre.att.bound: matrix, rows = pre.term, cols CI.lower, CI.upper (one-sided,
##                  used by the equivalence test)
##   pre.att.boot : matrix, rows = pre.term, cols = surviving replicates
##
## `dloo.adjust = TRUE` selects Liu's pre-treatment-average baseline. SEs are
## normal-approximation (matching fect's default ci.method = "normal").
##
##   boot.pre       : matrix, rows = pre.term, cols = surviving replicates
##   boot.pre.group : array [n_group x length(pre.term) x replicates] or NULL,
##                    group rows aligned to `pt$group_labels`.
##   Y, X, scale    : the fitted outcome and covariates (fits with covariates
##                    only), passed on to .dloo_prepare().
## ---------------------------------------------------------------------------
.dloo_fill <- function(fit,
                       dloo.adjust    = FALSE,
                       alpha          = 0.05,
                       controls       = "not-yet-treated",
                       quantile.CI    = FALSE,
                       vartype        = "bootstrap",
                       boot.pre       = NULL,
                       boot.pre.group = NULL,
                       Y              = NULL,
                       X              = NULL,
                       scale          = 1) {

    prep <- .dloo_prepare(fit, controls, Y = Y, X = X, scale = scale)

    ## Point estimate (unit weights all 1), aligned to the object's pre.term.
    pt <- .dloo_eval(prep, rep(1, prep$N), dloo.adjust)

    pre.term <- fit$time[fit$time <= 0]
    pre.term <- sort(pre.term)
    pos <- match(pre.term, pt$event.time)
    att   <- pt$att[pos]
    count <- pt$n.units[pos]

    ## Per-subgroup point estimates, aligned to pre.term. Populated only when the
    ## fit carried a `group` breakdown.
    has_group  <- !is.null(pt$group_att)
    grp_labels <- pt$group_labels
    if (has_group) {
        att.group   <- pt$group_att[, pos, drop = FALSE]
        count.group <- pt$group_count[, pos, drop = FALSE]
    }

    if (is.null(boot.pre)) {
        stop("dloo inference requires the native bootstrap draws ",
             "(from fect_boot). None were supplied.", call. = FALSE)
    }
    ## boot.pre rows already correspond to pre.term (fect_boot aligned them).
    boot <- boot.pre

    ## Assemble the fect-format placebo triple (pre.est.att / pre.att.bound /
    ## pre.att.boot) from a point-estimate series, its bootstrap draws and cell
    ## counts. CIs / bounds / p-values follow EXACTLY what fect_boot uses for its
    ## own est.att (so dloo matches `loo`): jackknife-pseudovalue SEs under
    ## vartype = "jackknife"; otherwise sd-of-draws with either normal-approx
    ## (quantile.CI = FALSE, the default) or basic bootstrap CIs
    ## (quantile.CI = TRUE, via .basic_ci_shifted / an empirical two-sided
    ## p-value). Shared by the pooled series and each subgroup.
    .boot_pvalue <- function(vec) {          # mirrors fect_boot's get.pvalue
        ok <- !is.na(vec) & !is.nan(vec)
        n  <- sum(ok)
        if (n == 0) return(NA_real_)
        a <- sum(vec[ok] >= 0) / n * 2
        b <- sum(vec[ok] <= 0) / n * 2
        min(min(a, b), 1)
    }
    .assemble <- function(att, boot, count, terms = pre.term) {
        if (identical(vartype, "jackknife")) {
            j  <- jackknifed(att, boot, alpha, quantile.CI = quantile.CI)
            se <- j$se; ci.lower <- j$CI.l; ci.upper <- j$CI.u; p.value <- j$P
            bd.lower <- att + stats::qnorm(alpha) * se
            bd.upper <- att + stats::qnorm(1 - alpha) * se
        } else {
            is_param <- identical(vartype, "parametric")
            se <- apply(boot, 1, stats::sd, na.rm = TRUE)
            if (!isTRUE(quantile.CI)) {
                ci.lower <- att - se * stats::qnorm(1 - alpha / 2)
                ci.upper <- att + se * stats::qnorm(1 - alpha / 2)
                p.value  <- 2 * (1 - stats::pnorm(abs(att / se)))
                bd.lower <- att - se * stats::qnorm(1 - alpha)   # one-sided
                bd.upper <- att + se * stats::qnorm(1 - alpha)
            } else {
                cm <- .basic_ci_shifted(att, boot, alpha, is_param)
                ci.lower <- cm[, 1]; ci.upper <- cm[, 2]
                p.value  <- apply(boot, 1, .boot_pvalue)
                bm <- .basic_ci_shifted(att, boot, 2 * alpha, is_param)
                bd.lower <- bm[, 1]; bd.upper <- bm[, 2]
            }
        }
        pre.est.att <- cbind(att, se, ci.lower, ci.upper, p.value, count)
        colnames(pre.est.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                   "p.value", "count.on")
        rownames(pre.est.att) <- terms
        pre.att.bound <- cbind(bd.lower, bd.upper)
        colnames(pre.att.bound) <- c("CI.lower", "CI.upper")
        rownames(pre.att.bound) <- terms
        rownames(boot) <- terms
        list(pre.est.att = pre.est.att, pre.att.bound = pre.att.bound,
             pre.att.boot = boot)
    }

    pooled <- .assemble(att, boot, count)

    ## Per-subgroup slots, keyed by raw group label, mirroring the `loo`
    ## contract (list of pre.est.att / pre.att.bound / pre.att.boot per group)
    ## so they view through plot(fit, loo = TRUE, show.group = ...).
    ##
    ## Each subgroup is built on ITS OWN event-time axis, not the pooled one --
    ## exactly as `loo` does. A subgroup only spans event times its own cohorts
    ## reach: a subgroup made of cohort 6 has no cell at event time -6, so
    ## carrying the pooled axis would advertise rows that cannot exist, and
    ## plot()'s `which(t1 == t0[1])` (against the subgroup's own att.on axis)
    ## would return integer(0) and error with "argument of length 0". Note this
    ## only bites when the group correlates with adoption timing -- region,
    ## sector, state -- which is the usual reason to group at all.
    ##
    ## `keep` is always contiguous and ends at 0: each cohort contributes the
    ## run {1-g+1, ..., 0}, so a union over a subgroup's cohorts is the run of
    ## its latest-adopting cohort.
    pre.est.group.output <- NULL
    if (has_group && !is.null(boot.pre.group)) {
        pre.est.group.output <- stats::setNames(
            lapply(seq_along(grp_labels), function(gi) {
                keep <- which(count.group[gi, ] > 0 & !is.na(att.group[gi, ]))
                if (length(keep) == 0) return(NULL)
                .assemble(att.group[gi, keep],
                          matrix(boot.pre.group[gi, keep, ],
                                 nrow = length(keep)),
                          count.group[gi, keep],
                          terms = pre.term[keep])
            }),
            grp_labels)
    }

    list(pre.est.att   = pooled$pre.est.att,
         pre.att.bound = pooled$pre.att.bound,
         pre.att.boot  = pooled$pre.att.boot,
         pre.est.group.output = pre.est.group.output,
         pre.term      = pre.term)
}
