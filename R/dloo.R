###############################################################################
## Double leave-one-out (dloo) pre-trend placebos computed directly from
## the default in-sample (fect-default) imputation fit.
##
## These are INTERNAL helpers used by fect() when it is called with
## `dloo = TRUE` (and, optionally, `dloo_adjust = TRUE`). They are not
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
##   for dloo, B excludes the tested period t; for dloo_adjust, B keeps it (see
##   below).
##
##   Two reasons we average this way rather than applying the (g-1)/(g-2) closed
##   form directly: (i) no rescaling constant is ever formed, so it degrades
##   gracefully when a cohort has few baselines (a short B) instead of dividing
##   by g-2; and (ii) switching between dloo and dloo_adjust is just a change to
##   which periods are in B. This reproduces (D) exactly: because the in-sample
##   pre-period placebos sum to zero over t' < g, subtracting the leave-one-out
##   baseline mean rescales the tested-period term by (g-1)/(g-2), recovering the
##   closed form above.
##
## The dloo_adjust correction (Liu's "pre-treatment average" baseline).
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
##   `dloo_adjust = TRUE` selects here --- a minimal change to which in-sample
##   estimates enter the aggregation (no rescaling constant is ever formed).
##   See the paper for the equality ATT^avg_g(t) - ATT^avg_g(t*) = the 2x2 DiD.
##
## Exactness. The overlay matches the literal "re-fit the imputation model on
## the restricted control pool per (cohort, period)" dloo to machine precision
## for a balanced control panel and additive two-way FE (method = "fe"). With
## latent factors / matrix completion the FE do not cleanly cancel; fect()
## warns before calling in that case.
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
## ---------------------------------------------------------------------------
.dloo_prepare <- function(fit, controls = "not-yet-treated") {

    eff     <- fit$eff
    D       <- fit$D.dat
    I       <- fit$I.dat
    rawtime <- fit$rawtime
    TT      <- nrow(eff)
    N       <- ncol(eff)

    ## Control (= not-treated, observed) cell mask.
    ctrl_mask <- (!is.na(I) & I == 1) & (!is.na(D) & D == 0)

    ## First-treated calendar time per unit; Inf for never-treated.
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
    ## NB: this >= 2 requirement is enforced here, so it applies to dloo_adjust
    ## too --- a single-pre-period cohort is excluded outright, NOT entered as
    ## the mechanically-zero placebo the dloo_adjust note below might suggest.
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

    list(eff = eff, ctrl_mask = ctrl_mask, rawtime = rawtime,
         TT = TT, N = N, cohorts = cohorts, n_cohorts = n_cohorts,
         cohort_idx = cohort_idx, never_idx = never_idx,
         controls = controls,
         cell_cohort = cohort_v, cell_tcol = tcol_v,
         cell_event = et_v,
         cell_cohort_label = cohorts[cohort_v],
         unit_group = unit_group, group_labels = group_labels)
}


## ---------------------------------------------------------------------------
## Evaluate the overlay for a vector of per-unit weights `w`. Returns the
## per-cell placebos and the size-weighted event-study series. `correct = TRUE`
## keeps the tested period in the baseline set B (Liu's pre-treatment-average /
## dloo_adjust); `correct = FALSE` drops it (the conventional dloo). Used both
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

    ## Weighted per-cohort x period control means M[ci, s] and weight totals
    ## Wt[ci, s] (for pooling later-adopter controls).
    M  <- matrix(NA_real_, n_cohorts, TT)
    Wt <- matrix(0, n_cohorts, TT)
    for (ci in seq_len(n_cohorts)) {
        cols <- which(cohort_idx == ci)
        if (length(cols) == 0) next
        e_sub <- eff[, cols, drop = FALSE]
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

    ## Per-cohort size under weights w (for the size-weighted aggregate).
    wsize <- vapply(seq_len(n_cohorts),
                    function(ci) sum(w[cohort_idx == ci]), numeric(1))

    ncell <- length(prep$cell_cohort)
    val   <- rep(NA_real_, ncell)

    ## Later-adopter control mean ebar depends on the target cohort gi only
    ## through the "later" set; cache once per distinct gi.
    uniq_gi <- unique(prep$cell_cohort)
    ebar_cache <- vector("list", n_cohorts)
    for (gi in uniq_gi) {
        g <- cohorts[gi]
        later_idx <- if (prep$controls == "never-treated") prep$never_idx else
                     which(cohorts > g)
        Wl <- Wt[later_idx, , drop = FALSE]
        Ml <- M[later_idx, , drop = FALSE]
        num <- colSums(Wl * ifelse(is.na(Ml), 0, Ml))
        den <- colSums(Wl)
        ebar_cache[[gi]] <- ifelse(den > 0, num / den, NA_real_)
    }

    for (k in seq_len(ncell)) {
        gi   <- prep$cell_cohort[k]
        tcol <- prep$cell_tcol[k]
        g    <- cohorts[gi]
        pre  <- which(rawtime < g)
        eb   <- ebar_cache[[gi]]
        Dvec <- M[gi, ] - eb                       # decontaminated placebo, all periods
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
            eb   <- ebar_cache[[gi]]
            base <- if (correct) pre else setdiff(pre, tcol)
            if (length(base) == 0) next
            if (is.na(eb[tcol]) || any(is.na(eb[base]))) next
            eti  <- match(prep$cell_event[k], ets)
            for (j in which(cohort_idx == gi)) {
                grp <- prep$unit_group[j]
                if (is.na(grp)) next
                dj_t <- eff[tcol, j] - eb[tcol]
                dj_b <- eff[base, j] - eb[base]
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
## the full TT x N x nboots panels (no keep.sims needed).
## ---------------------------------------------------------------------------
.dloo_boot_draw <- function(eff, D, I, rawtime, group.map, G,
                            controls, correct, pre.term) {
    pseudo <- list(eff = eff, D.dat = D, I.dat = I, rawtime = rawtime,
                   group = group.map, G = G)
    prep <- .dloo_prepare(pseudo, controls)
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
## `dloo_adjust = TRUE` selects Liu's pre-treatment-average baseline. SEs are
## normal-approximation (matching fect's default ci.method = "normal").
##
##   boot.pre       : matrix, rows = pre.term, cols = surviving replicates
##   boot.pre.group : array [n_group x length(pre.term) x replicates] or NULL,
##                    group rows aligned to `pt$group_labels`.
## ---------------------------------------------------------------------------
.dloo_fill <- function(fit,
                       dloo_adjust    = FALSE,
                       alpha          = 0.05,
                       controls       = "not-yet-treated",
                       quantile.CI    = FALSE,
                       vartype        = "bootstrap",
                       boot.pre       = NULL,
                       boot.pre.group = NULL) {

    prep <- .dloo_prepare(fit, controls)

    ## Point estimate (unit weights all 1), aligned to the object's pre.term.
    pt <- .dloo_eval(prep, rep(1, prep$N), dloo_adjust)

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
    .assemble <- function(att, boot, count) {
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
        rownames(pre.est.att) <- pre.term
        pre.att.bound <- cbind(bd.lower, bd.upper)
        colnames(pre.att.bound) <- c("CI.lower", "CI.upper")
        rownames(pre.att.bound) <- pre.term
        rownames(boot) <- pre.term
        list(pre.est.att = pre.est.att, pre.att.bound = pre.att.bound,
             pre.att.boot = boot)
    }

    pooled <- .assemble(att, boot, count)

    ## Per-subgroup slots, keyed by raw group label, mirroring the `loo`
    ## contract (list of pre.est.att / pre.att.bound / pre.att.boot per group)
    ## so they view through plot(fit, loo = TRUE, show.group = ...).
    pre.est.group.output <- NULL
    if (has_group && !is.null(boot.pre.group)) {
        pre.est.group.output <- stats::setNames(
            lapply(seq_along(grp_labels), function(gi)
                .assemble(att.group[gi, ],
                          matrix(boot.pre.group[gi, , ],
                                 nrow = length(pre.term)),
                          count.group[gi, ])),
            grp_labels)
    }

    list(pre.est.att   = pooled$pre.est.att,
         pre.att.bound = pooled$pre.att.bound,
         pre.att.boot  = pooled$pre.att.boot,
         pre.est.group.output = pre.est.group.output,
         pre.term      = pre.term)
}
