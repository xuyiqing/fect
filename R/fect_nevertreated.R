fect_nevertreated <- function(Y, # Outcome variable, (T*N) matrix
                        X, # Explanatory variables:  (T*N*p) array
                        D, #  Indicator for treated unit (tr==1)
                        W,
                        I,
                        II,
                        cm = FALSE,
                        II.cm = NULL,
                        T.on,
                        T.off = NULL,
                        T.on.carry = NULL,
                        T.on.balance = NULL,
                        balance.period = NULL,
                        CV = TRUE,
                        criterion = "mspe",
                        cv.method = "rolling",
                        k = 20,
                        cv.prop = 0.1,
                        cv.nobs = 3,
                        cv.donut = 1,
                        cv.buffer = 1,
                        min.T0 = 5,
                        r = 0, # r.end when CV==TRUE
                        r.end = 3,
                        binary = FALSE,
                        QR = FALSE,
                        force,
                        hasRevs = 1,
                        tol, # tolerance level
                        max.iteration = 1000,
                        boot = FALSE, # bootstrapped sample
                        placeboTest = 0,
                        placebo.period = NULL,
                        carryoverTest = 0,
                        carryover.period = NULL,
                        norm.para = NULL,
                        calendar.enp.seq = NULL,
                        time.on.seq = NULL,
                        time.off.seq = NULL,
                        time.on.balance.seq = NULL,
                        time.on.seq.W = NULL,
                        time.off.seq.W = NULL,
                        group.level = NULL,
                        group = NULL,
                        time.on.seq.group = NULL,
                        time.off.seq.group = NULL,
                        ## CFE-specific parameters (new)
                        method = "ife",           ## "ife" (existing) or "cfe" (new)
                        X.extra.FE = NULL,        ## TT x N x n_extra array
                        X.Z = NULL,               ## TT x N x n_z array
                        X.Q = NULL,               ## TT x N x n_q array
                        X.gamma = NULL,           ## TT x N x n_gamma array
                        X.kappa = NULL,           ## TT x N x n_kappa array
                        Zgamma.id = NULL,         ## list mapping gamma groups to Z columns
                        kappaQ.id = NULL,         ## list mapping kappa groups to Q columns
                        parallel = TRUE,
                        cores = NULL,
                        do_parallel_cv = NULL,  ## pre-computed flag from default.R; NULL means derive from parallel
                        loading.bound = "none",
                        gamma.loading = NULL,
                        gamma.loading.grid = NULL,
                        cv.rule = "1se",
                        W.in.fit = TRUE
                        ) {
    ## -------------------------------##
    ## Parsing data
    ## -------------------------------##
    cv.rule <- .fect_validate_cv_rule(cv.rule)
    carryover.pos <- placebo.pos <- na.pos <- NULL
    res.sd1 <- res.sd2 <- NULL
    ## unit id and time
    TT <- dim(Y)[1]
    N <- dim(Y)[2]
    id <- 1:N
    time <- 1:TT
    if (is.null(X) == FALSE) {
        p <- dim(X)[3]
    } else {
        p <- 0
        X <- array(0, dim = c(1, 1, 0))
    }

    ## replicate data
    YY <- Y
    YY[which(II == 0)] <- 0 ## reset to 0

    ## once treated, always treated
    ## careful unbalanced case: calculate treated units
    ## treat reversals as always treated
    D <- apply(D, 2, function(vec) {
        cumsum(vec)
    })
    D <- ifelse(D > 0, 1, 0)
    D.sum <- colSums(D)
    tr <- which(D.sum >= 1)
    Ntr <- length(tr)
    co <- which(D.sum == 0)
    Nco <- length(co)
    r <- min(r, TT, Nco)

    I.tr <- as.matrix(I[, tr]) ## maybe only 1 treated unit
    I.co <- I[, co]
    II.tr <- as.matrix(II[, tr])
    II.co <- II[, co]
    Y.tr <- as.matrix(Y[, tr])
    Y.co <- as.matrix(Y[, co])
    YY.tr <- as.matrix(YY[, tr])
    YY.co <- as.matrix(YY[, co])
    Ntr <- length(tr)
    Nco <- length(co)
    if (p == 0) {
        X.tr <- array(0, dim = c(TT, Ntr, 0))
        X.co <- array(0, dim = c(TT, Nco, 0))
    } else {
        X.tr <- array(NA, dim = c(TT, Ntr, p))
        X.co <- array(NA, dim = c(TT, Nco, p))
        for (j in 1:p) {
            X.tr[, , j] <- X[, tr, j]
            X.co[, , j] <- X[, co, j]
        }
    }

    ## ---- CFE array subsetting (co/tr split) ----
    if (method == "cfe") {
        ## helper to split a TT x N x k array into co/tr slices
        .split_array <- function(arr, co_idx, tr_idx, TT, Nco, Ntr) {
            if (!is.null(arr) && length(dim(arr)) == 3 && dim(arr)[3] > 0) {
                list(co = arr[, co_idx, , drop = FALSE],
                     tr = arr[, tr_idx, , drop = FALSE])
            } else {
                list(co = array(0, dim = c(TT, Nco, 0)),
                     tr = array(0, dim = c(TT, Ntr, 0)))
            }
        }
        efe.split  <- .split_array(X.extra.FE, co, tr, TT, Nco, Ntr)
        X.extra.FE.co <- efe.split$co;  X.extra.FE.tr <- efe.split$tr
        xz.split   <- .split_array(X.Z, co, tr, TT, Nco, Ntr)
        X.Z.co <- xz.split$co;  X.Z.tr <- xz.split$tr
        xq.split   <- .split_array(X.Q, co, tr, TT, Nco, Ntr)
        X.Q.co <- xq.split$co;  X.Q.tr <- xq.split$tr
        xg.split   <- .split_array(X.gamma, co, tr, TT, Nco, Ntr)
        X.gamma.co <- xg.split$co;  X.gamma.tr <- xg.split$tr
        xk.split   <- .split_array(X.kappa, co, tr, TT, Nco, Ntr)
        X.kappa.co <- xk.split$co;  X.kappa.tr <- xk.split$tr

        ## ---- Extra FE Type A/B classification ----
        n_extra <- dim(X.extra.FE.co)[3]
        fe_type <- character(0)
        typeA_idx <- integer(0)
        typeB_idx <- integer(0)

        if (n_extra > 0) {
            fe_type <- character(n_extra)
            for (k in 1:n_extra) {
                co_levels <- unique(X.extra.FE[1, co, k])
                tr_levels <- unique(X.extra.FE[1, tr, k])
                overlap <- intersect(co_levels, tr_levels)
                if (length(overlap) == 0) {
                    fe_type[k] <- "A"
                } else {
                    fe_type[k] <- "B"
                    missing <- setdiff(tr_levels, co_levels)
                    if (length(missing) > 0) {
                        stop(paste0("Extra fixed effect dimension ", k,
                             " has levels in treated units not found in controls: ",
                             paste(missing, collapse = ", "),
                             ". Cannot estimate from controls."))
                    }
                }
            }
            typeA_idx <- which(fe_type == "A")
            typeB_idx <- which(fe_type == "B")
        }

        ## Build co-only Type-B FE array for complex_fe_ub
        if (length(typeB_idx) > 0) {
            X.extra.FE.co.B <- X.extra.FE.co[, , typeB_idx, drop = FALSE]
        } else {
            X.extra.FE.co.B <- array(0, dim = c(TT, Nco, 0))
        }
    }

    if (is.null(W) || !W.in.fit) {
        W.use <- as.matrix(0)
    } else {
        W.use <- as.matrix(W[, co, drop = FALSE])
        W.use[which(II.co == 0)] <- 0
    }

    ## ---- cv.method validation ---- ##
    cv.method <- .fect_normalize_cv_method(
        cv.method,
        allowed = c("rolling", "block", "all_units", "treated_units", "loo")
    )

    ## ---- W for treated units (scoring) ---- ##
    if (!is.null(W)) {
        W.tr <- as.matrix(W[, tr, drop = FALSE])
    } else {
        W.tr <- NULL
    }

    if (!0 %in% I.tr) {
        ## a (TT*Ntr) matrix, time dimension: before treatment
        pre <- as.matrix(D[, tr] == 0 & II[, tr] == 1)
    } else {
        pre <- as.matrix(D[, tr] == 0 & I[, tr] == 1 & II[, tr] == 1)
    }

    T0 <- apply(pre, 2, sum)
    T0.min <- min(T0)
    pre.v <- as.vector(pre) ## vectorized "pre-treatment" indicator
    id.tr.pre.v <- rep(id, each = TT)[which(pre.v == 1)] ## vectorized pre-treatment grouping variable for the treated
    time.pre <- split(rep(time, Ntr)[which(pre.v == 1)], id.tr.pre.v) ## a list of pre-treatment periods
    sameT0 <- length(unique(T0)) == 1

    ## ---- count.T.cv construction ---- ##
    count.T.cv <- NULL
    if (!is.null(T.on)) {
        t.on.full <- c(T.on)
        count.T.cv <- table(t.on.full)
        count.T.cv <- count.T.cv[which(as.numeric(names(count.T.cv)) <= 0)]
        if (length(count.T.cv) > 0) {
            count.T.cv <- count.T.cv / mean(count.T.cv)
            nm <- names(count.T.cv)
            count.T.cv <- c(count.T.cv, median(count.T.cv))
            names(count.T.cv) <- c(nm, "Control")
        }
    }

    ## ====================================================================
    ## Solver dispatch wrapper
    ## ====================================================================
    ## Dispatches to inter_fe (balanced) / inter_fe_ub (unbalanced) for IFE,
    ## or complex_fe_ub for CFE. Eliminates balanced/unbalanced branching
    ## at each call site.
    .estimate_co <- function(Y, Y0, X, I_obs, W_in, beta0_in, r, force,
                             tol, max.iteration, use_cfe = FALSE,
                             center = TRUE) {
        ## Center Y to improve convergence conditioning:
        ## removes grand mean from fit so tol applies to variation, not level.
        mu_init <- 0
        if (center && 0 %in% I_obs) {
            ## Only center for unbalanced panels (where EM iterates)
            mu_init <- sum(Y * I_obs) / sum(I_obs)
            Y  <- Y - mu_init * I_obs   ## observed positions centered, zeros stay
            Y0 <- Y0 - mu_init          ## initial fit centered
        }

        if (use_cfe) {
            out <- complex_fe_ub(Y, Y0, X,
                X.extra.FE.co.B, X.Z.co, X.Q.co, X.gamma.co, X.kappa.co,
                Zgamma.id, kappaQ.id,
                I_obs, W_in, beta0_in, r, force = force, tol, max.iteration)
        } else if (!0 %in% I_obs) {
            out <- inter_fe(Y, X, r, force = force, beta0_in = beta0_in, tol, max.iteration)
        } else {
            out <- inter_fe_ub(Y, Y0, X, I_obs, W_in, beta0_in, r,
                        force = force, tol, max.iteration)
        }

        ## Undo centering
        if (mu_init != 0) {
            out$mu  <- out$mu + mu_init
            out$fit <- out$fit + mu_init
        }
        out
    }

  if (method == "ife") {
    ## ====================================================================
    ## IFE PATH
    ## ====================================================================

    beta0 <- matrix(0, p, 1)
    # initial fit using Y.co
    if (0 %in% I.co) {
        data.ini <- matrix(NA, Nco * TT, (p + 3))
        data.ini[, 1] <- c(Y.co)
        data.ini[, 2] <- rep(1:Nco, each = TT)
        data.ini[, 3] <- rep(1:TT, Nco)
        if (p > 0) {
            for (i in 1:p) {
                data.ini[, (3 + i)] <- c(X.co[, , i])
            }
        }
        ## observed Y0 indicator:
        initialOut <- Y0.co <- beta0 <- FE0 <- xi0 <- factor0 <- NULL
        oci <- which(c(II.co) == 1)
        if (is.null(W) || !W.in.fit) {
            initialOut <- initialFit(data = data.ini, force = force, oci = oci)
        } else {
            initialOut <- initialFit(data = data.ini, force = force, w = c(W.use), oci = oci)
        }
        Y0.co <- initialOut$Y0
        beta0 <- initialOut$beta0
        if (p > 0 && sum(is.na(beta0)) > 0) {
            beta0[which(is.na(beta0))] <- 0
        }
    }


    validX <- 1 ## no multi-colinearity

    ## cross-validation only for gsynth
    if (CV == TRUE) {
        ## Two-tier tolerance: CV uses looser tol for r-selection speed,
        ## final estimation (after CV) uses the user's tol for precision.
        cv_tol <- max(tol, 1e-3)

        ## starting r
        if ((r > (T0.min - 1) & force %in% c(0, 2)) | (r > (T0.min - 2) & force %in% c(1, 3))) {
            message("r is too big compared with T0; reset to 0.")
            r <- 0
        }

        ## store all MSPE
        if (force %in% c(0, 2)) {
            r.max <- max(min((T0.min - 1), r.end), 0)
        } else {
            r.max <- max(min((T0.min - 2), r.end), 0)
        }

        if (r.max == 0) {
            r.cv <- 0
            message("Cross validation cannot be performed since available pre-treatment records of treated units are too few. So set r.cv = 0.")
            est.co.best <- .estimate_co(YY.co, Y0.co, X.co, I.co, W.use, beta0, 0, force, cv_tol, max.iteration)
        } else {
            r.old <- r ## save the minimal number of factors

            message("Cross-validating ...", "\r")
            score_names <- c("MSPE", "WMSPE", "GMSPE", "WGMSPE",
                             "MAD", "Moment", "GMoment", "RMSE", "Bias")
            CV.out <- matrix(NA, (r.max - r.old + 1), 4 + length(score_names))
            colnames(CV.out) <- c("r", "sigma2", "IC", "PC", score_names)
            CV.out[, "r"] <- c(r.old:r.max)
            CV.out[, score_names] <- 1e10
            CV.out[, "PC"] <- 1e10
            r.pc <- est.co.pc.best <- NULL

            ## Per-fold SE matrix parallel to CV.out (added v2.3.0). Populated
            ## below in both parallel and serial CV branches; consumed at the
            ## end of the IFE CV block to apply `cv.rule` (default "1se").
            CV.out.se <- matrix(NA_real_, nrow(CV.out), ncol(CV.out))
            colnames(CV.out.se) <- colnames(CV.out)
            CV.out.se[, "r"] <- CV.out[, "r"]

            crit_col <- switch(criterion,
                mspe = "MSPE", wmspe = "WMSPE", gmspe = "GMSPE", wgmspe = "WGMSPE",
                mad = "MAD", moment = "Moment", gmoment = "GMoment", "MSPE")

            ## ---- cv.sample pre-computation (IFE) ---- ##
            if (cv.method != "loo" && r.max > 0) {
                if (cv.method %in% c("all_units", "rolling")) {
                    rm.count.co <- floor(sum(II.co) * cv.prop)
                    if (rm.count.co == 0 && cv.method == "all_units") {
                        message("cv.prop too small for control panel; falling back to LOO.")
                        cv.method <- "loo"
                    } else {
                        D.co.fake <- matrix(0, TT, Nco)
                        oci.co <- which(c(II.co) == 1)
                        ## Ensure data.ini exists (balanced IFE case skips its creation)
                        if (!exists("data.ini", inherits = FALSE) || is.null(data.ini)) {
                            data.ini <- matrix(NA, Nco * TT, (p + 3))
                            data.ini[, 1] <- c(Y.co)
                            data.ini[, 2] <- rep(1:Nco, each = TT)
                            data.ini[, 3] <- rep(1:TT, Nco)
                            if (p > 0) {
                                for (i.tmp in 1:p) {
                                    data.ini[, (3 + i.tmp)] <- c(X.co[, , i.tmp])
                                }
                            }
                        }
                        rmCV <- list()
                        ociCV <- list()
                        estCV <- list()
                        Y0CV.co <- array(NA, dim = c(TT, Nco, k))
                        if (p > 0) {
                            beta0CV.co <- array(NA, dim = c(p, 1, k))
                        } else {
                            beta0CV.co <- array(0, dim = c(1, 0, k))
                        }
                        flag.cv <- 0

                        ## ---- rolling-window pre-computation ---- ##
                        ## When rolling, masks come from .build_cv_mask_rolling
                        ## (per-fold sampling of cv.prop of eligible control
                        ## units; per-unit random anchor + cv.nobs scored
                        ## holdout + cv.buffer past-side buffer + drop-future
                        ## as the rolling-window step). Skip the con1/con2
                        ## block-CV feasibility checks --- rolling preserves
                        ## per-time donor coverage by construction via
                        ## per-fold unit sampling.
                        rolling_folds <- NULL
                        if (cv.method == "rolling") {
                            rolling_folds <- .build_cv_mask_rolling(
                                II = II.co, D = D.co.fake, k = k,
                                cv.nobs = cv.nobs, cv.buffer = cv.buffer,
                                cv.prop = cv.prop, min.T0 = min.T0, seed = NULL
                            )
                        }

                        for (i.cv in 1:k) {
                            if (cv.method == "rolling") {
                                cv.n <- 0
                                cv.id  <- rolling_folds[[i.cv]]$cv.id
                                est.id <- rolling_folds[[i.cv]]$est.id
                            } else {
                                cv.n <- 0
                                repeat {
                                    cv.n <- cv.n + 1
                                    get.cv <- cv.sample(II.co, D.co.fake,
                                        count = rm.count.co,
                                        cv.count = cv.nobs,
                                        cv.treat = FALSE,
                                        cv.donut = cv.donut)
                                    cv.id <- get.cv$cv.id
                                    II.co.cv <- II.co
                                    II.co.cv[cv.id] <- 0
                                    II.co.cv.valid <- II.co
                                    II.co.cv.valid[cv.id] <- -1
                                    con1 <- sum(apply(II.co.cv, 1, sum) >= 1) == TT
                                    con2 <- sum(apply(II.co.cv, 2, sum) >= min.T0) == Nco
                                    if (con1 && con2) break
                                    if (cv.n >= 200) {
                                        flag.cv <- 1
                                        keep.1 <- which(apply(II.co.cv, 1, sum) < 1)
                                        keep.2 <- which(apply(II.co.cv, 2, sum) < min.T0)
                                        II.co.cv[keep.1, ] <- II.co[keep.1, ]
                                        II.co.cv[, keep.2] <- II.co[, keep.2]
                                        II.co.cv.valid[keep.1, ] <- II.co[keep.1, ]
                                        II.co.cv.valid[, keep.2] <- II.co[, keep.2]
                                        cv.id <- which(II.co.cv.valid != II.co)
                                        break
                                    }
                                }
                            }
                            rmCV[[i.cv]] <- cv.id
                            ociCV[[i.cv]] <- setdiff(oci.co, cv.id)
                            if (cv.method == "rolling") {
                                estCV[[i.cv]] <- est.id
                            } else if (cv.n < 200) {
                                estCV[[i.cv]] <- get.cv$est.id
                            } else {
                                cv.diff <- setdiff(get.cv$cv.id, cv.id)
                                estCV[[i.cv]] <- setdiff(get.cv$est.id, cv.diff)
                            }
                            if (is.null(W) || !W.in.fit) {
                                initialOutCv <- initialFit(data = data.ini, force = force, oci = ociCV[[i.cv]])
                            } else {
                                initialOutCv <- initialFit(data = data.ini, force = force, w = c(W.use), oci = ociCV[[i.cv]])
                            }
                            Y0CV.co[, , i.cv] <- initialOutCv$Y0
                            if (p > 0) {
                                beta0cv <- initialOutCv$beta0
                                if (sum(is.na(beta0cv)) > 0) {
                                    beta0cv[which(is.na(beta0cv))] <- 0
                                }
                                beta0CV.co[, , i.cv] <- beta0cv
                            }
                        }
                        if (flag.cv == 1) {
                            message("Some control units have too few observations. Removed automatically in CV.\n")
                        }
                    }
                } else if (cv.method == "treated_units") {
                    rm.count.tr <- floor(sum(pre) * cv.prop)
                    if (rm.count.tr == 0) {
                        message("cv.prop too small for treated pre-treatment panel; falling back to LOO.")
                        cv.method <- "loo"
                    } else {
                        D.tr.fake <- matrix(0, TT, Ntr)
                        rmCV.tr <- list()
                        estCV.tr <- list()
                        flag.cv <- 0
                        for (i.cv in 1:k) {
                            cv.n <- 0
                            repeat {
                                cv.n <- cv.n + 1
                                get.cv <- cv.sample(pre, D.tr.fake,
                                    count = rm.count.tr,
                                    cv.count = cv.nobs,
                                    cv.treat = FALSE,
                                    cv.donut = cv.donut)
                                cv.id <- get.cv$cv.id
                                pre.cv <- pre
                                pre.cv[cv.id] <- 0
                                con1 <- TRUE
                                pre.rows <- which(rowSums(pre) > 0)
                                if (length(pre.rows) > 0) {
                                    con1 <- all(rowSums(pre.cv[pre.rows, , drop = FALSE]) >= 1)
                                }
                                con2 <- all(colSums(pre.cv) >= min(min.T0, 2))
                                if (con1 && con2) break
                                if (cv.n >= 200) {
                                    flag.cv <- 1
                                    pre.cv.valid <- pre
                                    pre.cv.valid[cv.id] <- -1
                                    keep.1 <- pre.rows[rowSums(pre.cv[pre.rows, , drop = FALSE]) < 1]
                                    keep.2 <- which(colSums(pre.cv) < min(min.T0, 2))
                                    if (length(keep.1) > 0) {
                                        pre.cv[keep.1, ] <- pre[keep.1, ]
                                        pre.cv.valid[keep.1, ] <- pre[keep.1, ]
                                    }
                                    if (length(keep.2) > 0) {
                                        pre.cv[, keep.2] <- pre[, keep.2]
                                        pre.cv.valid[, keep.2] <- pre[, keep.2]
                                    }
                                    cv.id <- which(pre.cv.valid != pre)
                                    break
                                }
                            }
                            rmCV.tr[[i.cv]] <- cv.id
                            if (cv.n < 200) {
                                estCV.tr[[i.cv]] <- get.cv$est.id
                            } else {
                                cv.diff <- setdiff(get.cv$cv.id, cv.id)
                                estCV.tr[[i.cv]] <- setdiff(get.cv$est.id, cv.diff)
                            }
                        }
                        if (flag.cv == 1) {
                            message("Some treated units have too few pre-treatment observations. Removed automatically in CV.\n")
                        }
                    }
                }
            }

            ## ---- Parallel backend setup (IFE CV) — Phase 3 flat dispatch ---- ##
            ## do_parallel_cv: pre-computed flag from default.R/cv.R (NULL means derive from parallel).
            ## parallel=TRUE (default from fect): auto-enable for all_units
            ##   when control panel is large enough (Nco*TT > .CV_PARALLEL_THRESH$ife)
            ## parallel=FALSE: force sequential (user override)
            ## parallel="cv": force parallel regardless of threshold
            if (is.null(do_parallel_cv)) {
                ## backwards-compat: derive from parallel if not pre-computed
                do_parallel_cv <- isTRUE(parallel) || "cv" %in% as.character(parallel)
            }
            if (!do_parallel_cv) {
                cv_ife_parallel <- FALSE
            } else {
                ## Explicit "cv" override: bypass threshold; auto mode: threshold gates
                use_explicit_cv_ife <- "cv" %in% as.character(parallel) && !isTRUE(parallel)
                ## Centralized threshold gate — replaces bespoke Nco * TT > 20000 check
                ife_threshold_met <- (Nco * TT) > .CV_PARALLEL_THRESH$ife
                cv_ife_parallel <- (cv.method %in% c("all_units", "rolling")) &&
                                   (ife_threshold_met || use_explicit_cv_ife) &&
                                   (k > 1)
            }
            if (cv_ife_parallel) {
                if (is.null(cores)) {
                    cores <- max(1L, min(parallelly::availableCores(omit = 2L), 8L))
                }
                old.future.plan.ife <- future::plan()
                on.exit(future::plan(old.future.plan.ife), add = TRUE, after = FALSE)
                future::plan(future::cluster, workers = .fect_make_future_cluster(cores))
                ## doFuture::registerDoFuture() removed — not needed for future_lapply dispatch
                avail <- parallelly::availableCores()
                msg_line <- sprintf("Parallel CV: using %d of %d available cores.", cores, avail)
                pad <- strrep(" ", max(0, 56 - nchar(msg_line)))
                message("\n",
                    " +----------------------------------------------------------+\n",
                    " | ", msg_line, pad, " |\n",
                    " |                                                          |\n",
                    " | To change: set cores = <n> in fect().                    |\n",
                    " | Default: min(available - 2, 8).                          |\n",
                    " +----------------------------------------------------------+\n")
            }

            ## ---- Task list construction (IFE flat r×k dispatch) ---- ##
            if (cv_ife_parallel) {
                r_seq_ife <- CV.out[, "r"]
                tasks_ife <- vector("list", length(r_seq_ife) * k)
                idx <- 1L
                for (ri in seq_along(r_seq_ife)) {
                    for (ii in 1:k) {
                        tasks_ife[[idx]] <- list(r = r_seq_ife[ri], ii = ii, ri = ri)
                        idx <- idx + 1L
                    }
                }
                ## Capture helper in closure for worker serialization
                .score_fn_ife_all <- .fect_cv_score_one_ife_nt_all
            }

            if (cv_ife_parallel) {
            ## ---- PARALLEL BRANCH: flat r×k future_lapply dispatch (IFE all_units) ---- ##

                ## Step 1: dispatch all (r, fold) scoring tasks
                fold_scores_ife <- future.apply::future_lapply(
                    tasks_ife,
                    FUN = function(task) {
                        .score_fn_ife_all(
                            ii            = task$ii,
                            YY.co         = YY.co,
                            Y0CV.co       = Y0CV.co,
                            X.co          = X.co,
                            II.co         = II.co,
                            W.use         = W.use,
                            W             = W,
                            beta0CV.co    = beta0CV.co,
                            rmCV          = rmCV,
                            estCV         = estCV,
                            r             = task$r,
                            force         = force,
                            cv_tol        = cv_tol,
                            max.iteration = max.iteration
                        )
                    },
                    future.seed     = TRUE,
                    future.packages = "fect"
                )

                ## Step 2: sequential master walk — apply 1% rule in rank order
                n_r_ife <- length(r_seq_ife)
                for (i in seq_len(n_r_ife)) {
                    r <- unname(r_seq_ife[i])
                    ## Full-data fit (sequential in master) — needed for sigma2/IC/PC
                    est.co <- .estimate_co(YY.co, Y0.co, X.co, I.co, W.use, beta0, r, force, cv_tol, max.iteration)

                    if (p > 0) {
                        na.pos <- is.nan(est.co$beta)
                        beta <- est.co$beta
                        beta[is.nan(est.co$beta)] <- 0
                    }
                    if (is.null(norm.para)) {
                        sigma2 <- est.co$sigma2; IC <- est.co$IC; PC <- est.co$PC
                    } else {
                        sigma2 <- est.co$sigma2 * (norm.para[1]^2)
                        IC <- est.co$IC - log(est.co$sigma2) + log(sigma2)
                        PC <- est.co$PC * (norm.para[1]^2)
                    }

                    ## Aggregate fold scores for this rank
                    task_idx <- which(vapply(tasks_ife, function(t) t$ri == i, logical(1)))
                    all_resid    <- unlist(lapply(fold_scores_ife[task_idx], `[[`, "resid"))
                    all_time_idx <- unlist(lapply(fold_scores_ife[task_idx], `[[`, "time_idx"))
                    all_obs_w    <- if (!is.null(W)) unlist(lapply(fold_scores_ife[task_idx], `[[`, "obs_w")) else c()

                    if (length(all_resid) == 0) {
                        scores <- c(MSPE = Inf, WMSPE = Inf, GMSPE = Inf, WGMSPE = Inf,
                                    MAD = Inf, Moment = Inf, GMoment = Inf, RMSE = Inf, Bias = Inf)
                        se_v <- setNames(rep(NA_real_, length(scores)), names(scores))
                    } else {
                        agg <- .fect_cv_aggregate_folds(
                            fold_list  = fold_scores_ife[task_idx],
                            count.T.cv = count.T.cv,
                            use_weight = as.integer(!is.null(W)),
                            norm.para  = NULL
                        )
                        scores <- agg$pooled
                        se_v   <- agg$se
                    }

                    ## 1% rule — identical logic to serial path
                    if ((min(CV.out[, crit_col]) - scores[crit_col]) > 0.01 * min(CV.out[, crit_col])) {
                        est.co.best <- est.co
                        r.cv <- r
                    } else {
                        if (r == r.cv + 1) message("*")
                    }
                    if (PC < min(CV.out[, "PC"])) {
                        r.pc <- r
                        est.co.pc.best <- est.co
                    }
                    CV.out[i, 2:4] <- c(sigma2, IC, PC)
                    CV.out[i, score_names] <- scores[score_names]
                    ## Persist per-fold SEs for end-of-loop cv.rule application
                    for (cn in score_names) {
                        if (cn %in% names(se_v)) CV.out.se[i, cn] <- se_v[cn]
                    }
                    message("r = ", r, "; sigma2 = ",
                        sprintf("%.5f", sigma2), "; IC = ",
                        sprintf("%.5f", IC), "; PC = ",
                        sprintf("%.5f", PC), "; MSPE = ",
                        sprintf("%.5f", scores["MSPE"]), sep = "")
                } ## end per-r master walk (IFE parallel)

            } else {
            ## ---- SERIAL BRANCH (existing r-loop, all cv.method values) ---- ##

            for (i in 1:dim(CV.out)[1]) {
                r <- unname(CV.out[i, "r"])
                est.co <- .estimate_co(YY.co, Y0.co, X.co, I.co, W.use, beta0, r, force, cv_tol, max.iteration)

                if (p > 0) {
                    na.pos <- is.nan(est.co$beta)
                    beta <- est.co$beta
                    beta[is.nan(est.co$beta)] <- 0
                }

                if (is.null(norm.para)) {
                    sigma2 <- est.co$sigma2
                    IC <- est.co$IC
                    PC <- est.co$PC
                } else {
                    sigma2 <- est.co$sigma2 * (norm.para[1]^2)
                    IC <- est.co$IC - log(est.co$sigma2) + log(sigma2)
                    PC <- est.co$PC * (norm.para[1]^2)
                }

              if (cv.method == "loo") {
                ## ---- LOO CV (existing code) ---- ##
                if (r != 0) {
                    F.hat <- as.matrix(est.co$factor)
                    if (force %in% c(1, 3)) {
                        F.hat <- cbind(F.hat, rep(1, TT))
                    }
                }
                U.tr <- Y.tr
                if (p > 0) {
                    for (j in 1:p) {
                        U.tr <- U.tr - X.tr[, , j] * beta[j]
                    }
                }
                if (force != 0) {
                    U.tr <- U.tr - matrix(est.co$mu, TT, Ntr) ## grand mean
                }
                if (force %in% c(2, 3)) {
                    U.tr <- U.tr - matrix(est.co$xi, TT, Ntr, byrow = FALSE)
                }

                if (0 %in% I.tr) {
                    U.tr[which(I.tr == 0)] <- 0
                }

                U.sav <- U.tr

                resid_all <- c()
                for (lv in unique(unlist(time.pre))) {
                    U.tr <- U.sav
                    if (max(T0) == T0.min & (!0 %in% I.tr)) {
                        U.lv <- as.matrix(U.tr[setdiff(c(1:T0.min), lv), ]) ## setdiff : x
                    } else {
                        U.tr.pre.v <- as.vector(U.tr)[which(pre.v == 1)] ## pre-treatment residual in a vector
                        U.tr.pre <- split(U.tr.pre.v, id.tr.pre.v) ##  a list of pretreatment residuals
                        if (!0 %in% I.tr) {
                            U.lv <- lapply(U.tr.pre, function(vec) {
                                return(vec[-lv])
                            }) ## a list
                        } else {
                            ## U.tr.pre.sav <- U.tr.pre
                            for (i.tr in 1:Ntr) {
                                U.tmp <- U.tr.pre[[i.tr]]
                                U.tr.pre[[i.tr]] <- U.tmp[!time.pre[[i.tr]] == lv]
                            }
                            U.lv <- U.tr.pre
                        }
                    }

                    if (r == 0) {
                        if (force %in% c(1, 3)) { ## take out unit fixed effect
                            if (max(T0) == T0.min & (!0 %in% I.tr)) {
                                alpha.tr.lv <- colMeans(U.lv)
                                U.tr <- U.tr - matrix(alpha.tr.lv, TT, Ntr, byrow = TRUE)
                            } else {
                                alpha.tr.lv <- sapply(U.lv, mean)
                                U.tr <- U.tr - matrix(alpha.tr.lv, TT, Ntr, byrow = TRUE)
                            }
                        }
                        e <- U.tr[which(time == lv), ] ## that period
                    } else {
                        F.lv <- as.matrix(F.hat[which(time != lv), ])
                        if (max(T0) == T0.min & (!0 %in% I.tr)) {
                            F.lv.pre <- F.hat[setdiff(c(1:T0.min), lv), ]
                            lambda.lv <- try(
                                solve(t(F.lv.pre) %*% F.lv.pre) %*% t(F.lv.pre) %*% U.lv,
                                silent = TRUE
                            )
                            if ("try-error" %in% class(lambda.lv)) {
                                break
                            }
                        } else {
                            if (!0 %in% I.tr) {
                                lambda.lv <- try(as.matrix(sapply(U.lv, function(vec) {
                                    F.lv.pre <- as.matrix(F.lv[1:length(vec), ])
                                    l.lv.tr <- solve(t(F.lv.pre) %*% F.lv.pre) %*% t(F.lv.pre) %*% vec
                                    return(l.lv.tr)
                                })), silent = TRUE)
                                if ("try-error" %in% class(lambda.lv)) {
                                    break
                                } else {
                                    if ((r == 1) & (force %in% c(0, 2))) {
                                        lambda.lv <- t(lambda.lv)
                                    }
                                }
                            } else {
                                if (force %in% c(1, 3)) {
                                    lambda.lv <- matrix(NA, (r + 1), Ntr)
                                } else {
                                    lambda.lv <- matrix(NA, r, Ntr)
                                }
                                test <- try(
                                    for (i.tr in 1:Ntr) {
                                        F.lv.pre <- as.matrix(F.hat[setdiff(time.pre[[i.tr]], lv), ])
                                        lambda.lv[, i.tr] <- solve(t(F.lv.pre) %*% F.lv.pre) %*% t(F.lv.pre) %*% as.matrix(U.lv[[i.tr]])
                                    },
                                    silent = TRUE
                                )
                                if ("try-error" %in% class(test)) {
                                    break
                                }
                            }
                        }
                        lambda.lv <- t(lambda.lv) ## N*r
                        e <- U.tr[which(time == lv), ] - c(F.hat[which(time == lv), ] %*% t(lambda.lv))
                    }
                    if (sameT0 == FALSE | 0 %in% I.tr) {
                        e <- e[which(pre[which(time == lv), ] == TRUE)]
                    }
                    ## accumulate residuals
                    resid_all <- c(resid_all, e)
                } ## end of leave-one-out

                if (length(resid_all) == 0) {
                    scores <- c(MSPE = Inf, WMSPE = Inf, GMSPE = Inf, WGMSPE = Inf,
                                MAD = Inf, Moment = Inf, GMoment = Inf, RMSE = Inf, Bias = Inf)
                } else {
                    ## Build time indices for LOO residuals
                    time_idx_loo <- NULL
                    obs_w_loo <- NULL
                    if (!is.null(count.T.cv)) {
                        time_idx_loo <- c()
                        for (lv in unique(unlist(time.pre))) {
                            if (sameT0 == FALSE | 0 %in% I.tr) {
                                n_resid_lv <- sum(pre[which(time == lv), ] == TRUE)
                            } else {
                                n_resid_lv <- Ntr
                            }
                            if (n_resid_lv > 0) {
                                t.on.tr <- T.on[, tr, drop = FALSE]
                                t.on.lv <- unique(t.on.tr[lv, ])
                                t.on.lv <- t.on.lv[!is.na(t.on.lv)]
                                if (length(t.on.lv) > 0) {
                                    time_idx_loo <- c(time_idx_loo, rep(as.character(t.on.lv[1]), n_resid_lv))
                                } else {
                                    time_idx_loo <- c(time_idx_loo, rep("Control", n_resid_lv))
                                }
                            }
                        }
                    }
                    if (!is.null(W.tr)) {
                        obs_w_loo <- c()
                        for (lv in unique(unlist(time.pre))) {
                            if (sameT0 == FALSE | 0 %in% I.tr) {
                                w_lv <- W.tr[lv, which(pre[which(time == lv), ] == TRUE)]
                            } else {
                                w_lv <- W.tr[lv, ]
                            }
                            obs_w_loo <- c(obs_w_loo, w_lv)
                        }
                    }
                    scores <- .score_residuals(
                        resid_all,
                        obs_weights = obs_w_loo,
                        time_index = time_idx_loo,
                        count_weights = count.T.cv,
                        norm.para = norm.para
                    )
                }

              } else if (cv.method %in% c("all_units", "rolling")) {
                ## ---- cv.sample "all_units" / "rolling" IFE CV (serial path — lapply only) ---- ##
                ## Rolling reuses the all_units scoring helper; only the rmCV/estCV
                ## fold construction differs (built via .build_cv_mask_rolling above).
                fold_results <- lapply(1:k, function(ii) {
                    .fect_cv_score_one_ife_nt_all(
                        ii       = ii,
                        YY.co    = YY.co,
                        Y0CV.co  = Y0CV.co,
                        X.co     = X.co,
                        II.co    = II.co,
                        W.use    = W.use,
                        W        = W,
                        beta0CV.co = beta0CV.co,
                        rmCV     = rmCV,
                        estCV    = estCV,
                        r        = r,
                        force    = force,
                        cv_tol   = cv_tol,
                        max.iteration = max.iteration
                    )
                })
                if (length(fold_results) == 0L) {
                    scores <- c(MSPE = Inf, WMSPE = Inf, GMSPE = Inf, WGMSPE = Inf,
                                MAD = Inf, Moment = Inf, GMoment = Inf, RMSE = Inf, Bias = Inf)
                    se_v <- setNames(rep(NA_real_, length(scores)), names(scores))
                } else {
                    agg <- .fect_cv_aggregate_folds(
                        fold_list  = fold_results,
                        count.T.cv = count.T.cv,
                        use_weight = as.integer(!is.null(W)),
                        norm.para  = NULL
                    )
                    scores <- agg$pooled
                    se_v   <- agg$se
                }

              } else {
                ## ---- cv.sample "treated_units" IFE CV (serial path — lapply only) ---- ##
                if (r != 0) {
                    F.hat <- as.matrix(est.co$factor)
                    if (force %in% c(1, 3)) {
                        F.hat <- cbind(F.hat, rep(1, TT))
                    }
                }
                U.tr <- Y.tr
                if (p > 0) {
                    for (j in 1:p) {
                        U.tr <- U.tr - X.tr[, , j] * beta[j]
                    }
                }
                if (force != 0) {
                    U.tr <- U.tr - matrix(est.co$mu, TT, Ntr)
                }
                if (force %in% c(2, 3)) {
                    U.tr <- U.tr - matrix(est.co$xi, TT, Ntr, byrow = FALSE)
                }
                if (0 %in% I.tr) {
                    U.tr[which(I.tr == 0)] <- 0
                }

                fold_results <- lapply(1:k, function(ii) {
                    .fect_cv_score_one_ife_nt_tr(
                        ii       = ii,
                        U.tr     = U.tr,
                        F.hat    = F.hat,
                        pre      = pre,
                        r        = r,
                        force    = force,
                        rmCV.tr  = rmCV.tr,
                        estCV.tr = estCV.tr,
                        W.tr     = W.tr,
                        T.on     = T.on,
                        tr       = tr,
                        TT       = TT,
                        Ntr      = Ntr
                    )
                })
                if (length(fold_results) == 0L) {
                    scores <- c(MSPE = Inf, WMSPE = Inf, GMSPE = Inf, WGMSPE = Inf,
                                MAD = Inf, Moment = Inf, GMoment = Inf, RMSE = Inf, Bias = Inf)
                    se_v <- setNames(rep(NA_real_, length(scores)), names(scores))
                } else {
                    agg <- .fect_cv_aggregate_folds(
                        fold_list  = fold_results,
                        count.T.cv = count.T.cv,
                        use_weight = as.integer(!is.null(W.tr)),
                        norm.para  = norm.para
                    )
                    scores <- agg$pooled
                    se_v   <- agg$se
                }

              } ## end cv.method branching

                if ((min(CV.out[, crit_col]) - scores[crit_col]) > 0.01 * min(CV.out[, crit_col])) {
                    ## at least 1% improvement for selected criterion
                    est.co.best <- est.co ## interFE result with the best r
                    r.cv <- r
                } else {
                    if (r == r.cv + 1) message("*")
                }

                if (PC < min(CV.out[, "PC"])) {
                    r.pc <- r
                    est.co.pc.best <- est.co
                }
                CV.out[i, 2:4] <- c(sigma2, IC, PC)
                CV.out[i, score_names] <- scores[score_names]
                ## Persist per-fold SEs for end-of-loop cv.rule application.
                ## Some serial paths (e.g., LOO branch) may leave se_v undefined;
                ## existsCheck protects those cases.
                if (exists("se_v", inherits = FALSE) && !is.null(se_v)) {
                    for (cn in score_names) {
                        if (cn %in% names(se_v)) CV.out.se[i, cn] <- se_v[cn]
                    }
                }
                message("r = ", r, "; sigma2 = ",
                    sprintf("%.5f", sigma2), "; IC = ",
                    sprintf("%.5f", IC), "; PC = ",
                    sprintf("%.5f", PC), "; MSPE = ",
                    sprintf("%.5f", scores["MSPE"]),
                    sep = ""
                )
            } ## end of while: search for r_star over

            } ## end SERIAL BRANCH (IFE)

            MSPE.best <- min(CV.out[, "MSPE"])

            ## --- Apply cv.rule (added v2.3.0) -----------------------------
            ## Override r.cv based on the user-selected rule. The default "1se"
            ## picks the smallest r within one fold-SE of the minimum-CV-error r;
            ## "min" picks the argmin; "1pct" preserves the legacy 1% rule from
            ## the in-loop assignments above.
            if (criterion %in% c("mspe","wmspe","gmspe","wgmspe","mad","moment","gmoment")) {
                means <- CV.out[, crit_col]
                ses   <- CV.out.se[, crit_col]
                ## Treat sentinels (1e10, Inf, NA) as missing for selection.
                means[!is.finite(means) | means >= 1e9] <- NA_real_
                i_pick <- .fect_apply_cv_rule(means, ses, rule = cv.rule)
                if (!is.na(i_pick) && i_pick >= 1L && i_pick <= nrow(CV.out)) {
                    new_r_cv <- unname(CV.out[i_pick, "r"])
                    if (!is.null(new_r_cv) && is.finite(new_r_cv)) {
                        if (new_r_cv != as.integer(unname(r.cv))) {
                            message(sprintf(
                                "  [cv.rule = %s] r.cv adjusted from %d to %d (1-SE band)",
                                cv.rule,
                                as.integer(unname(r.cv)),
                                as.integer(new_r_cv)
                            ))
                            est.co.best <- .estimate_co(
                                YY.co, Y0.co, X.co, I.co, W.use, beta0,
                                as.integer(new_r_cv), force, cv_tol, max.iteration
                            )
                        }
                        ## Preserve the in-loop names convention: serial path
                        ## leaves r.cv unnamed; parallel path sets names "r".
                        had_name <- !is.null(names(r.cv))
                        r.cv <- new_r_cv
                        if (had_name) names(r.cv) <- "r"
                    }
                }
            }
            ## --------------------------------------------------------------

            if (r > (T0.min - 1)) {
                message(" (r hits maximum)")
            }
            message("\n r* = ", r.cv, sep = "")
            message("\n")
        }
    } else {
        r.cv <- r
        r.min <- r.max <- r
    }

    est.co.fect <- NULL

    est.co.best <- .estimate_co(YY.co, Y0.co, X.co, II.co, W.use, beta0, r.cv, force, tol, max.iteration)

    if (boot == FALSE) {
        if (r.cv == 0) {
            est.co.fect <- est.co.best
        } else {
            est.co.fect <- .estimate_co(YY.co, Y0.co, X.co, II.co, W.use, beta0, 0, force, tol, max.iteration)
        }
    }
    validX <- est.co.best$validX
    validF <- ifelse(r.cv > 0, 1, 0)

    # get the counterfactual of X
    ## ## take out the effect of X
    U.tr.r0 <- U.tr <- Y.tr
    if (p > 0) {
        beta <- est.co.best$beta

        if (est.co.best$validX == 0) {
            beta <- matrix(0, p, 1)
        } else {
            beta <- est.co.best$beta
            beta[is.nan(est.co.best$beta)] <- 0
        }
        for (j in 1:p) {
            U.tr <- U.tr - X.tr[, , j] * beta[j]
        }

        if (boot == FALSE) {
            beta.r0 <- est.co.fect$beta
            if (est.co.fect$validX == 0) {
                beta.r0 <- matrix(0, p, 1)
            } else {
                beta.r0 <- est.co.fect$beta
                beta.r0[is.nan(est.co.fect$beta)] <- 0
            }
            for (j in 1:p) {
                U.tr.r0 <- U.tr.r0 - X.tr[, , j] * beta.r0[j]
            }
        }
    } else {
        beta <- NA
        beta.r0 <- NA
    }

    mu <- est.co.best$mu
    U.tr <- U.tr - matrix(mu, TT, Ntr) ## grand mean
    Y.fe.bar <- rep(mu, TT)

    if (boot == FALSE) {
        mu.r0 <- est.co.fect$mu
        U.tr.r0 <- U.tr.r0 - matrix(mu.r0, TT, Ntr)
        Y.fe.bar.r0 <- rep(mu.r0, TT)
    }

    if (force %in% c(2, 3)) {
        xi <- est.co.best$xi ## a (TT*1) matrix
        U.tr <- U.tr - matrix(c(xi), TT, Ntr, byrow = FALSE) ## will be adjusted at last
        Y.fe.bar <- Y.fe.bar + xi
        if (boot == FALSE) {
            xi.r0 <- est.co.fect$xi ## a (TT*1) matrix
            U.tr.r0 <- U.tr.r0 - matrix(c(xi.r0), TT, Ntr, byrow = FALSE)
            Y.fe.bar.r0 <- Y.fe.bar.r0 + xi.r0
        }
    }

    if (max(T0) == T0.min & (!0 %in% I.tr)) {
        U.tr.pre <- as.matrix(U.tr[1:T0.min, ])
        if (boot == FALSE) {
            U.tr.pre.r0 <- as.matrix(U.tr.r0[1:T0.min, ])
        }
    } else {
        ## not necessary to reset utr for ub data for pre.v doesn't include them
        U.tr.pre.v <- as.vector(U.tr)[which(pre.v == 1)] # pre-treatment residual in a vector
        U.tr.pre <- split(U.tr.pre.v, id.tr.pre.v) ##  a list of pretreatment residuals

        if (boot == FALSE) {
            U.tr.pre.v.r0 <- as.vector(U.tr.r0)[which(pre.v == 1)] # pre-treatment residual in a vector
            U.tr.pre.r0 <- split(U.tr.pre.v.r0, id.tr.pre.v) ##  a list of pretreatment residuals
        }
    }

    ## the error structure

    # for r=0
    if (force %in% c(1, 3)) { ## take out unit fixed effect
        if ((max(T0) == T0.min) & (!0 %in% I.tr)) {
            if (boot == FALSE) {
                alpha.tr.r0 <- as.matrix(colMeans(U.tr.pre.r0))
                U.tr.r0 <- U.tr.r0 - matrix(alpha.tr.r0, TT, Ntr, byrow = TRUE)
            }
        } else {
            if (boot == FALSE) {
                alpha.tr.r0 <- as.matrix(sapply(U.tr.pre.r0, mean))
                U.tr.r0 <- U.tr.r0 - matrix(alpha.tr.r0, TT, Ntr, byrow = TRUE)
            }
        }
    }
    if (boot == FALSE) {
        eff.r0 <- U.tr.r0
    }



    if (r.cv == 0) {
        if (force %in% c(1, 3)) { ## take out unit fixed effect
            if ((max(T0) == T0.min) & (!0 %in% I.tr)) {
                alpha.tr <- as.matrix(colMeans(U.tr.pre))
                U.tr <- U.tr - matrix(alpha.tr, TT, Ntr, byrow = TRUE)
            } else {
                alpha.tr <- as.matrix(sapply(U.tr.pre, mean))
                U.tr <- U.tr - matrix(alpha.tr, TT, Ntr, byrow = TRUE)
            }
        }
        eff <- U.tr
        lambda.tr <- NULL
        lambda.co <- NULL
    } else { ## Factors
        F.hat <- as.matrix(est.co.best$factor)
        if (force %in% c(1, 3)) {
            F.hat <- cbind(F.hat, rep(1, TT))
        }

        ## Bounded-loading dispatch: resolves use_bounded, gamma_use, and
        ## optional W_tr / loading.proj.resid diagnostics. See statsclaw-
        ## workspace/fect/runs/REQ-bounded-loadings/spec.md for the design.
        use_bounded <- identical(loading.bound, "simplex")
        W_tr <- NULL
        loading.proj.resid <- NULL

        if (use_bounded) {
            Lambda.co.raw <- as.matrix(est.co.best$lambda)   # Nco x r.cv (no intercept col)
            F.hat.raw     <- as.matrix(est.co.best$factor)   # TT  x r.cv (no intercept col)
            Nco_b         <- nrow(Lambda.co.raw)
            gamma_grid_use <- if (!is.null(gamma.loading.grid)) gamma.loading.grid
                              else .default_gamma_grid()
            gamma_use <- gamma.loading
            if (is.null(gamma_use)) {
                balanced_pre <- max(T0) == T0.min & !0 %in% I.tr
                if (balanced_pre) {
                    cv_F_pre <- F.hat.raw[1:T0.min, , drop = FALSE]
                    cv_U_pre <- U.tr.pre
                } else if (!0 %in% I.tr) {
                    ## Different T0 per unit; U.tr.pre is a list
                    cv_F_pre <- lapply(U.tr.pre, function(vec)
                        F.hat.raw[seq_along(vec), , drop = FALSE]
                    )
                    cv_U_pre <- lapply(U.tr.pre, as.numeric)
                } else {
                    ## Some missing observations in treated pre-period
                    cv_F_pre <- lapply(seq_along(U.tr.pre), function(i.tr)
                        F.hat.raw[time.pre[[i.tr]], , drop = FALSE]
                    )
                    cv_U_pre <- lapply(U.tr.pre, as.numeric)
                }
                gamma_use <- .cv_gamma_loading(
                    U_tr_pre   = cv_U_pre,
                    F_hat_pre  = cv_F_pre,
                    Lambda_co  = Lambda.co.raw,
                    gamma_grid = gamma_grid_use,
                    cv_k       = 5L
                )$gamma_cv
            }
        }

        ## Lambda_tr (Ntr*r) or (Ntr*(r+1))
        if (max(T0) == T0.min & (!0 %in% I.tr)) {
            F.hat.pre <- F.hat[1:T0.min, ]
            if (!use_bounded) {
                lambda.tr <- try(solve(t(F.hat.pre) %*% F.hat.pre) %*% t(F.hat.pre) %*% U.tr.pre,
                    silent = TRUE
                )
                if ("try-error" %in% class(lambda.tr)) {
                    return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1)))
                }
            } else {
                ## Bounded: solve simplex QP per treated unit on the r-col F block only
                F.pre.r <- F.hat.raw[1:T0.min, , drop = FALSE]
                lambda.tr.r <- matrix(NA_real_, nrow = r.cv, ncol = Ntr)
                W_tr <- matrix(NA_real_, nrow = Ntr, ncol = Nco_b)
                loading.proj.resid <- numeric(Ntr)
                for (i.tr in seq_len(Ntr)) {
                    sol <- .solve_bounded_loading(
                        u_pre     = U.tr.pre[, i.tr],
                        F_pre     = F.pre.r,
                        Lambda_co = Lambda.co.raw,
                        gamma     = gamma_use
                    )
                    lambda.tr.r[, i.tr] <- sol$lambda_hat
                    W_tr[i.tr, ] <- sol$w
                    loading.proj.resid[i.tr] <- sqrt(sum(
                        (U.tr.pre[, i.tr] - F.pre.r %*% sol$lambda_hat)^2
                    ))
                }
                if (force %in% c(1, 3)) {
                    ## alpha.tr is residual-mean (NOT bounded; per locked decision)
                    resid_mat <- U.tr.pre - F.pre.r %*% lambda.tr.r   # T0.min x Ntr
                    alpha_row <- matrix(colMeans(resid_mat), nrow = 1L)
                    lambda.tr <- rbind(lambda.tr.r, alpha_row)
                } else {
                    lambda.tr <- lambda.tr.r
                }
            }
        } else {
            if (!0 %in% I.tr) {
                if (!use_bounded) {
                    lambda.tr <- try(as.matrix(sapply(U.tr.pre, function(vec) {
                        F.hat.pre <- as.matrix(F.hat[1:length(vec), ])
                        l.tr <- solve(t(F.hat.pre) %*% F.hat.pre) %*% t(F.hat.pre) %*% vec
                        return(l.tr) ## a vector of each individual lambdas
                    })), silent = TRUE)
                    if ("try-error" %in% class(lambda.tr)) {
                        return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1)))
                        ## stop("Error occurs. Please set a smaller value of factor number.")
                    }
                    if ((r.cv == 1) & (force %in% c(0, 2))) {
                        lambda.tr <- t(lambda.tr)
                    }
                } else {
                    lambda.tr.r <- matrix(NA_real_, nrow = r.cv, ncol = Ntr)
                    W_tr <- matrix(NA_real_, nrow = Ntr, ncol = Nco_b)
                    loading.proj.resid <- numeric(Ntr)
                    alpha_vec <- numeric(Ntr)
                    for (i.tr in seq_len(Ntr)) {
                        vec <- U.tr.pre[[i.tr]]
                        F.pre.i <- F.hat.raw[seq_along(vec), , drop = FALSE]
                        sol <- .solve_bounded_loading(
                            u_pre     = vec,
                            F_pre     = F.pre.i,
                            Lambda_co = Lambda.co.raw,
                            gamma     = gamma_use
                        )
                        lambda.tr.r[, i.tr] <- sol$lambda_hat
                        W_tr[i.tr, ] <- sol$w
                        loading.proj.resid[i.tr] <- sqrt(sum(
                            (vec - F.pre.i %*% sol$lambda_hat)^2
                        ))
                        alpha_vec[i.tr] <- mean(vec - F.pre.i %*% sol$lambda_hat)
                    }
                    if (force %in% c(1, 3)) {
                        lambda.tr <- rbind(lambda.tr.r, matrix(alpha_vec, nrow = 1L))
                    } else {
                        lambda.tr <- lambda.tr.r
                    }
                }
            } else {
                if (!use_bounded) {
                    if (force %in% c(1, 3)) {
                        lambda.tr <- matrix(NA, (r.cv + 1), Ntr)
                    } else {
                        lambda.tr <- matrix(NA, r.cv, Ntr)
                    }
                    test <- try(
                        for (i.tr in 1:Ntr) {
                            F.hat.pre <- as.matrix(F.hat[time.pre[[i.tr]], ])
                            lambda.tr[, i.tr] <- solve(t(F.hat.pre) %*% F.hat.pre) %*% t(F.hat.pre) %*% as.matrix(U.tr.pre[[i.tr]])
                        },
                        silent = TRUE
                    )
                    if ("try-error" %in% class(test)) {
                        return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1), eff = matrix(NA, TT, Ntr)))
                        ## stop("Error occurs. Please set a smaller value of factor number.")
                    }
                } else {
                    lambda.tr.r <- matrix(NA_real_, nrow = r.cv, ncol = Ntr)
                    W_tr <- matrix(NA_real_, nrow = Ntr, ncol = Nco_b)
                    loading.proj.resid <- numeric(Ntr)
                    alpha_vec <- numeric(Ntr)
                    for (i.tr in seq_len(Ntr)) {
                        vec    <- as.numeric(U.tr.pre[[i.tr]])
                        F.pre.i <- F.hat.raw[time.pre[[i.tr]], , drop = FALSE]
                        sol <- .solve_bounded_loading(
                            u_pre     = vec,
                            F_pre     = F.pre.i,
                            Lambda_co = Lambda.co.raw,
                            gamma     = gamma_use
                        )
                        lambda.tr.r[, i.tr] <- sol$lambda_hat
                        W_tr[i.tr, ] <- sol$w
                        loading.proj.resid[i.tr] <- sqrt(sum(
                            (vec - F.pre.i %*% sol$lambda_hat)^2
                        ))
                        alpha_vec[i.tr] <- mean(vec - F.pre.i %*% sol$lambda_hat)
                    }
                    if (force %in% c(1, 3)) {
                        lambda.tr <- rbind(lambda.tr.r, matrix(alpha_vec, nrow = 1L))
                    } else {
                        lambda.tr <- lambda.tr.r
                    }
                }
            }
        }

        lambda.tr <- t(lambda.tr)
        eff <- U.tr - F.hat %*% t(lambda.tr)
        if (force %in% c(1, 3)) {
            alpha.tr <- as.matrix(lambda.tr[, (r.cv + 1), drop = FALSE])
            lambda.tr <- lambda.tr[, 1:r.cv, drop = FALSE]
        }
        if (boot == 0) {
            if (use_bounded && !is.null(W_tr)) {
                wgt.implied <- W_tr
            } else {
                inv.tr <- try(
                    ginv(t(as.matrix(lambda.tr))),
                    silent = TRUE
                )

                if (!"try-error" %in% class(inv.tr)) {
                    wgt.implied <- t(inv.tr %*% t(as.matrix(est.co.best$lambda)))
                }
            }
        }
    } ## end of r!=0 case


    if (0 %in% I.tr) {
        eff[which(I.tr == 0)] <- 0 ## adjust
        if (boot == FALSE) {
            eff.r0[which(I.tr == 0)] <- 0
        }
    } ## missing data will be adjusted to NA finally

  } else if (method == "cfe") {
    ## ====================================================================
    ## CFE PATH (new code)
    ## ====================================================================

    ## ---- Validate sufficient control units ----
    if (Nco < 2) {
        stop("Too few never-treated (control) units for CFE estimation. ",
             "At least 2 control units are required, but only ", Nco, " found.")
    }

    ## ---- Initial fit for co-only data ----
    beta0 <- matrix(0, p, 1)
    data.ini <- matrix(NA, Nco * TT, (p + 3))
    data.ini[, 1] <- c(Y.co)
    data.ini[, 2] <- rep(1:Nco, each = TT)
    data.ini[, 3] <- rep(1:TT, Nco)
    if (p > 0) {
        for (i in 1:p) {
            data.ini[, (3 + i)] <- c(X.co[, , i])
        }
    }
    initialOut <- Y0.co <- NULL
    oci <- which(c(II.co) == 1)
    if (is.null(W)) {
        initialOut <- initialFit(data = data.ini, force = force, oci = oci)
    } else {
        initialOut <- initialFit(data = data.ini, force = force, w = c(W.use), oci = oci)
    }
    Y0.co <- initialOut$Y0
    beta0 <- initialOut$beta0
    if (p > 0 && sum(is.na(beta0)) > 0) {
        beta0[which(is.na(beta0))] <- 0
    }

    validX <- 1

    ## ---- CV loop for CFE ----
    if (CV == TRUE) {
        ## Two-tier tolerance for CFE CV
        cv_tol <- max(tol, 1e-3)

        ## starting r
        if ((r > (T0.min - 1) & force %in% c(0, 2)) | (r > (T0.min - 2) & force %in% c(1, 3))) {
            message("r is too big compared with T0; reset to 0.")
            r <- 0
        }
        if (force %in% c(0, 2)) {
            r.max <- max(min((T0.min - 1), r.end), 0)
        } else {
            r.max <- max(min((T0.min - 2), r.end), 0)
        }

        if (r.max == 0) {
            r.cv <- 0
            message("Cross validation cannot be performed since available pre-treatment records of treated units are too few. So set r.cv = 0.")
            est.co.best <- complex_fe_ub(YY.co, Y0.co, X.co,
                X.extra.FE.co.B, X.Z.co, X.Q.co, X.gamma.co, X.kappa.co,
                Zgamma.id, kappaQ.id,
                II.co, W.use, beta0, 0, force = force, cv_tol, max.iteration)
        } else {
            r.old <- r
            message("Cross-validating ...", "\r")
            score_names <- c("MSPE", "WMSPE", "GMSPE", "WGMSPE",
                             "MAD", "Moment", "GMoment", "RMSE", "Bias")
            CV.out <- matrix(NA, (r.max - r.old + 1), 4 + length(score_names))
            colnames(CV.out) <- c("r", "sigma2", "IC", "PC", score_names)
            CV.out[, "r"] <- c(r.old:r.max)
            CV.out[, score_names] <- 1e10
            CV.out[, "PC"] <- 1e10
            r.pc <- est.co.pc.best <- NULL

            crit_col <- switch(criterion,
                mspe = "MSPE", wmspe = "WMSPE", gmspe = "GMSPE", wgmspe = "WGMSPE",
                mad = "MAD", moment = "Moment", gmoment = "GMoment", "MSPE")

            ## ---- cv.sample pre-computation (CFE) ---- ##
            if (cv.method != "loo" && r.max > 0) {
                if (cv.method %in% c("all_units", "rolling")) {
                    rm.count.co <- floor(sum(II.co) * cv.prop)
                    if (rm.count.co == 0 && cv.method == "all_units") {
                        message("cv.prop too small for control panel; falling back to LOO.")
                        cv.method <- "loo"
                    } else {
                        D.co.fake <- matrix(0, TT, Nco)
                        oci.co <- which(c(II.co) == 1)
                        rmCV <- list()
                        ociCV <- list()
                        estCV <- list()
                        Y0CV.co <- array(NA, dim = c(TT, Nco, k))
                        if (p > 0) {
                            beta0CV.co <- array(NA, dim = c(p, 1, k))
                        } else {
                            beta0CV.co <- array(0, dim = c(1, 0, k))
                        }
                        flag.cv <- 0

                        ## ---- rolling-window pre-computation (CFE) ---- ##
                        rolling_folds <- NULL
                        if (cv.method == "rolling") {
                            rolling_folds <- .build_cv_mask_rolling(
                                II = II.co, D = D.co.fake, k = k,
                                cv.nobs = cv.nobs, cv.buffer = cv.buffer,
                                cv.prop = cv.prop, min.T0 = min.T0, seed = NULL
                            )
                        }

                        for (i.cv in 1:k) {
                            if (cv.method == "rolling") {
                                cv.n <- 0
                                cv.id  <- rolling_folds[[i.cv]]$cv.id
                                est.id <- rolling_folds[[i.cv]]$est.id
                            } else {
                                cv.n <- 0
                                repeat {
                                    cv.n <- cv.n + 1
                                    get.cv <- cv.sample(II.co, D.co.fake,
                                        count = rm.count.co,
                                        cv.count = cv.nobs,
                                        cv.treat = FALSE,
                                        cv.donut = cv.donut)
                                    cv.id <- get.cv$cv.id
                                    II.co.cv <- II.co
                                    II.co.cv[cv.id] <- 0
                                    II.co.cv.valid <- II.co
                                    II.co.cv.valid[cv.id] <- -1
                                    con1 <- sum(apply(II.co.cv, 1, sum) >= 1) == TT
                                    con2 <- sum(apply(II.co.cv, 2, sum) >= min.T0) == Nco
                                    if (con1 && con2) break
                                    if (cv.n >= 200) {
                                        flag.cv <- 1
                                        keep.1 <- which(apply(II.co.cv, 1, sum) < 1)
                                        keep.2 <- which(apply(II.co.cv, 2, sum) < min.T0)
                                        II.co.cv[keep.1, ] <- II.co[keep.1, ]
                                        II.co.cv[, keep.2] <- II.co[, keep.2]
                                        II.co.cv.valid[keep.1, ] <- II.co[keep.1, ]
                                        II.co.cv.valid[, keep.2] <- II.co[, keep.2]
                                        cv.id <- which(II.co.cv.valid != II.co)
                                        break
                                    }
                                }
                            }
                            rmCV[[i.cv]] <- cv.id
                            ociCV[[i.cv]] <- setdiff(oci.co, cv.id)
                            if (cv.method == "rolling") {
                                estCV[[i.cv]] <- est.id
                            } else if (cv.n < 200) {
                                estCV[[i.cv]] <- get.cv$est.id
                            } else {
                                cv.diff <- setdiff(get.cv$cv.id, cv.id)
                                estCV[[i.cv]] <- setdiff(get.cv$est.id, cv.diff)
                            }
                            if (is.null(W) || !W.in.fit) {
                                initialOutCv <- initialFit(data = data.ini, force = force, oci = ociCV[[i.cv]])
                            } else {
                                initialOutCv <- initialFit(data = data.ini, force = force, w = c(W.use), oci = ociCV[[i.cv]])
                            }
                            Y0CV.co[, , i.cv] <- initialOutCv$Y0
                            if (p > 0) {
                                beta0cv <- initialOutCv$beta0
                                if (sum(is.na(beta0cv)) > 0) {
                                    beta0cv[which(is.na(beta0cv))] <- 0
                                }
                                beta0CV.co[, , i.cv] <- beta0cv
                            }
                        }
                        if (flag.cv == 1) {
                            message("Some control units have too few observations. Removed automatically in CV.\n")
                        }
                    }
                } else if (cv.method == "treated_units") {
                    rm.count.tr <- floor(sum(pre) * cv.prop)
                    if (rm.count.tr == 0) {
                        message("cv.prop too small for treated pre-treatment panel; falling back to LOO.")
                        cv.method <- "loo"
                    } else {
                        D.tr.fake <- matrix(0, TT, Ntr)
                        rmCV.tr <- list()
                        estCV.tr <- list()
                        flag.cv <- 0
                        for (i.cv in 1:k) {
                            cv.n <- 0
                            repeat {
                                cv.n <- cv.n + 1
                                get.cv <- cv.sample(pre, D.tr.fake,
                                    count = rm.count.tr,
                                    cv.count = cv.nobs,
                                    cv.treat = FALSE,
                                    cv.donut = cv.donut)
                                cv.id <- get.cv$cv.id
                                pre.cv <- pre
                                pre.cv[cv.id] <- 0
                                con1 <- TRUE
                                pre.rows <- which(rowSums(pre) > 0)
                                if (length(pre.rows) > 0) {
                                    con1 <- all(rowSums(pre.cv[pre.rows, , drop = FALSE]) >= 1)
                                }
                                con2 <- all(colSums(pre.cv) >= min(min.T0, 2))
                                if (con1 && con2) break
                                if (cv.n >= 200) {
                                    flag.cv <- 1
                                    pre.cv.valid <- pre
                                    pre.cv.valid[cv.id] <- -1
                                    keep.1 <- pre.rows[rowSums(pre.cv[pre.rows, , drop = FALSE]) < 1]
                                    keep.2 <- which(colSums(pre.cv) < min(min.T0, 2))
                                    if (length(keep.1) > 0) {
                                        pre.cv[keep.1, ] <- pre[keep.1, ]
                                        pre.cv.valid[keep.1, ] <- pre[keep.1, ]
                                    }
                                    if (length(keep.2) > 0) {
                                        pre.cv[, keep.2] <- pre[, keep.2]
                                        pre.cv.valid[, keep.2] <- pre[, keep.2]
                                    }
                                    cv.id <- which(pre.cv.valid != pre)
                                    break
                                }
                            }
                            rmCV.tr[[i.cv]] <- cv.id
                            if (cv.n < 200) {
                                estCV.tr[[i.cv]] <- get.cv$est.id
                            } else {
                                cv.diff <- setdiff(get.cv$cv.id, cv.id)
                                estCV.tr[[i.cv]] <- setdiff(get.cv$est.id, cv.diff)
                            }
                        }
                        if (flag.cv == 1) {
                            message("Some treated units have too few pre-treatment observations. Removed automatically in CV.\n")
                        }
                    }
                }
            }

            ## ---- Parallel backend setup (CFE CV) — Phase 3 flat dispatch ---- ##
            ## Reuse do_parallel_cv derived at IFE setup; if not yet set (CFE-only call), derive now.
            if (is.null(do_parallel_cv)) {
                do_parallel_cv <- isTRUE(parallel) || "cv" %in% as.character(parallel)
            }
            if (!do_parallel_cv) {
                cv_cfe_parallel <- FALSE
            } else {
                ## Re-derive explicitly for CFE block (do NOT reuse use_explicit_cv_ife)
                use_explicit_cv_cfe <- "cv" %in% as.character(parallel) && !isTRUE(parallel)
                ## Centralized threshold gate: CFE threshold is 60000L (higher overhead per fit)
                cfe_threshold_met <- (Nco * TT) > .CV_PARALLEL_THRESH$cfe
                cv_cfe_parallel <- (cv.method %in% c("all_units", "rolling")) &&
                                   (cfe_threshold_met || use_explicit_cv_cfe) &&
                                   (k > 1)
            }
            if (cv_cfe_parallel) {
                if (is.null(cores)) {
                    cores <- max(1L, min(parallelly::availableCores(omit = 2L), 8L))
                }
                old.future.plan.cfe <- future::plan()
                on.exit(future::plan(old.future.plan.cfe), add = TRUE, after = FALSE)
                future::plan(future::cluster, workers = .fect_make_future_cluster(cores))
                ## doFuture::registerDoFuture() removed — not needed for future_lapply dispatch
                avail <- parallelly::availableCores()
                msg_line <- sprintf("Parallel CV (CFE): using %d of %d available cores.", cores, avail)
                pad <- strrep(" ", max(0, 56 - nchar(msg_line)))
                message("\n",
                    " +----------------------------------------------------------+\n",
                    " | ", msg_line, pad, " |\n",
                    " |                                                          |\n",
                    " | To change: set cores = <n> in fect().                    |\n",
                    " | Default: min(available - 2, 8).                          |\n",
                    " +----------------------------------------------------------+\n")
            }

            ## ---- Task list construction (CFE flat r×k dispatch) ---- ##
            if (cv_cfe_parallel) {
                r_seq_cfe <- CV.out[, "r"]
                tasks_cfe <- vector("list", length(r_seq_cfe) * k)
                idx <- 1L
                for (ri in seq_along(r_seq_cfe)) {
                    for (ii in 1:k) {
                        tasks_cfe[[idx]] <- list(r = r_seq_cfe[ri], ii = ii, ri = ri)
                        idx <- idx + 1L
                    }
                }
                ## Capture helper in closure for worker serialization
                .score_fn_cfe_all <- .fect_cv_score_one_cfe_nt_all
            }

            if (cv_cfe_parallel) {
            ## ---- PARALLEL BRANCH: flat r×k future_lapply dispatch (CFE all_units) ---- ##

                ## Step 1: dispatch all (r, fold) scoring tasks
                fold_scores_cfe <- future.apply::future_lapply(
                    tasks_cfe,
                    FUN = function(task) {
                        .score_fn_cfe_all(
                            ii              = task$ii,
                            YY.co           = YY.co,
                            Y0CV.co         = Y0CV.co,
                            X.co            = X.co,
                            II.co           = II.co,
                            W.use           = W.use,
                            W               = W,
                            beta0CV.co      = beta0CV.co,
                            X.extra.FE.co.B = X.extra.FE.co.B,
                            X.Z.co          = X.Z.co,
                            X.Q.co          = X.Q.co,
                            X.gamma.co      = X.gamma.co,
                            X.kappa.co      = X.kappa.co,
                            Zgamma.id       = Zgamma.id,
                            kappaQ.id       = kappaQ.id,
                            rmCV            = rmCV,
                            estCV           = estCV,
                            r               = task$r,
                            force           = force,
                            cv_tol          = cv_tol,
                            max.iteration   = max.iteration
                        )
                    },
                    future.seed     = TRUE,
                    future.packages = "fect"
                )

                ## Step 2: sequential master walk — apply 1% rule in rank order
                n_r_cfe <- length(r_seq_cfe)
                for (i in seq_len(n_r_cfe)) {
                    r <- unname(r_seq_cfe[i])
                    ## Full-data fit (sequential in master) — needed for sigma2/IC/PC
                    est.co <- complex_fe_ub(YY.co, Y0.co, X.co,
                        X.extra.FE.co.B, X.Z.co, X.Q.co, X.gamma.co, X.kappa.co,
                        Zgamma.id, kappaQ.id,
                        II.co, W.use, beta0, r, force = force, cv_tol, max.iteration)

                    if (p > 0) {
                        na.pos <- is.nan(est.co$beta)
                        beta <- est.co$beta
                        beta[is.nan(est.co$beta)] <- 0
                    }
                    if (is.null(norm.para)) {
                        sigma2 <- est.co$sigma2; IC <- est.co$IC; PC <- est.co$PC
                    } else {
                        sigma2 <- est.co$sigma2 * (norm.para[1]^2)
                        IC <- est.co$IC - log(est.co$sigma2) + log(sigma2)
                        PC <- est.co$PC * (norm.para[1]^2)
                    }

                    ## Aggregate fold scores for this rank
                    task_idx <- which(vapply(tasks_cfe, function(t) t$ri == i, logical(1)))
                    all_resid    <- unlist(lapply(fold_scores_cfe[task_idx], `[[`, "resid"))
                    all_time_idx <- unlist(lapply(fold_scores_cfe[task_idx], `[[`, "time_idx"))
                    all_obs_w    <- if (!is.null(W)) unlist(lapply(fold_scores_cfe[task_idx], `[[`, "obs_w")) else c()

                    if (length(all_resid) == 0) {
                        scores <- c(MSPE = Inf, WMSPE = Inf, GMSPE = Inf, WGMSPE = Inf,
                                    MAD = Inf, Moment = Inf, GMoment = Inf, RMSE = Inf, Bias = Inf)
                    } else {
                        scores <- .score_residuals(
                            all_resid,
                            obs_weights   = if (!is.null(W)) all_obs_w else NULL,
                            time_index    = all_time_idx,
                            count_weights = count.T.cv,
                            norm.para     = NULL
                        )
                    }

                    ## 1% rule — identical logic to serial path
                    if ((min(CV.out[, crit_col]) - scores[crit_col]) > 0.01 * min(CV.out[, crit_col])) {
                        est.co.best <- est.co
                        r.cv <- r
                    } else {
                        if (r == r.cv + 1) message("*")
                    }
                    if (PC < min(CV.out[, "PC"])) {
                        r.pc <- r
                        est.co.pc.best <- est.co
                    }
                    CV.out[i, 2:4] <- c(sigma2, IC, PC)
                    CV.out[i, score_names] <- scores[score_names]
                    message("r = ", r, "; sigma2 = ",
                        sprintf("%.5f", sigma2), "; IC = ",
                        sprintf("%.5f", IC), "; PC = ",
                        sprintf("%.5f", PC), "; MSPE = ",
                        sprintf("%.5f", scores["MSPE"]), sep = "")
                } ## end per-r master walk (CFE parallel)

            } else {
            ## ---- SERIAL BRANCH (existing r-loop, all cv.method values) ---- ##

            for (i in 1:dim(CV.out)[1]) {
                r <- unname(CV.out[i, "r"])
                est.co <- complex_fe_ub(YY.co, Y0.co, X.co,
                    X.extra.FE.co.B, X.Z.co, X.Q.co, X.gamma.co, X.kappa.co,
                    Zgamma.id, kappaQ.id,
                    II.co, W.use, beta0, r, force = force, cv_tol, max.iteration)

                if (p > 0) {
                    na.pos <- is.nan(est.co$beta)
                    beta <- est.co$beta
                    beta[is.nan(est.co$beta)] <- 0
                }

                if (is.null(norm.para)) {
                    sigma2 <- est.co$sigma2
                    IC <- est.co$IC
                    PC <- est.co$PC
                } else {
                    sigma2 <- est.co$sigma2 * (norm.para[1]^2)
                    IC <- est.co$IC - log(est.co$sigma2) + log(sigma2)
                    PC <- est.co$PC * (norm.para[1]^2)
                }

              if (cv.method == "loo") {
                ## ---- LOO CV (existing code) ---- ##
                ## Build U.tr: subtract Layer 1 from Y.tr
                if (r != 0) {
                    F.hat <- as.matrix(est.co$factor)
                    if (force %in% c(1, 3)) {
                        F.hat <- cbind(F.hat, rep(1, TT))
                    }
                }
                U.tr <- Y.tr
                if (p > 0) {
                    for (j in 1:p) {
                        U.tr <- U.tr - X.tr[, , j] * beta[j]
                    }
                }
                if (force != 0) {
                    U.tr <- U.tr - matrix(est.co$mu, TT, Ntr)
                }
                if (force %in% c(2, 3)) {
                    U.tr <- U.tr - matrix(est.co$xi, TT, Ntr, byrow = FALSE)
                }

                ## Subtract gamma for treated (Layer 1)
                if (!is.null(est.co$gamma) && length(est.co$gamma) > 0) {
                    for (k_g in seq_along(est.co$gamma)) {
                        gamma.fit.tr.k <- .reconstruct_gamma_fit_tr(
                            est.co$gamma[[k_g]], X.Z.tr, X.gamma.tr[, , k_g, drop = FALSE],
                            Zgamma.id[[k_g]], TT, Ntr)
                        U.tr <- U.tr - gamma.fit.tr.k
                    }
                }

                ## Subtract Type-B extra FE for treated (Layer 1)
                if (length(typeB_idx) > 0) {
                    typeB.fit.tr <- .extract_and_apply_typeB_fe(
                        est.co, X.co, X.extra.FE.co.B, X.extra.FE.tr,
                        typeB_idx, X.Z.co, X.gamma.co, X.Q.co, X.kappa.co,
                        Zgamma.id, kappaQ.id,
                        TT, Nco, Ntr, p, r, force)
                    U.tr <- U.tr - typeB.fit.tr
                }

                if (0 %in% I.tr) {
                    U.tr[which(I.tr == 0)] <- 0
                }

                U.sav <- U.tr

                ## Leave-one-out CV (same structure as IFE)
                resid_all <- c()
                for (lv in unique(unlist(time.pre))) {
                    U.tr <- U.sav
                    if (max(T0) == T0.min & (!0 %in% I.tr)) {
                        U.lv <- as.matrix(U.tr[setdiff(c(1:T0.min), lv), ])
                    } else {
                        U.tr.pre.v <- as.vector(U.tr)[which(pre.v == 1)]
                        U.tr.pre <- split(U.tr.pre.v, id.tr.pre.v)
                        if (!0 %in% I.tr) {
                            U.lv <- lapply(U.tr.pre, function(vec) {
                                return(vec[-lv])
                            })
                        } else {
                            for (i.tr in 1:Ntr) {
                                U.tmp <- U.tr.pre[[i.tr]]
                                U.tr.pre[[i.tr]] <- U.tmp[!time.pre[[i.tr]] == lv]
                            }
                            U.lv <- U.tr.pre
                        }
                    }

                    if (r == 0) {
                        if (force %in% c(1, 3)) {
                            if (max(T0) == T0.min & (!0 %in% I.tr)) {
                                alpha.tr.lv <- colMeans(U.lv)
                                U.tr <- U.tr - matrix(alpha.tr.lv, TT, Ntr, byrow = TRUE)
                            } else {
                                alpha.tr.lv <- sapply(U.lv, mean)
                                U.tr <- U.tr - matrix(alpha.tr.lv, TT, Ntr, byrow = TRUE)
                            }
                        }
                        e <- U.tr[which(time == lv), ]
                    } else {
                        F.lv <- as.matrix(F.hat[which(time != lv), ])
                        if (max(T0) == T0.min & (!0 %in% I.tr)) {
                            F.lv.pre <- F.hat[setdiff(c(1:T0.min), lv), ]
                            lambda.lv <- try(
                                solve(t(F.lv.pre) %*% F.lv.pre) %*% t(F.lv.pre) %*% U.lv,
                                silent = TRUE
                            )
                            if ("try-error" %in% class(lambda.lv)) {
                                break
                            }
                        } else {
                            if (!0 %in% I.tr) {
                                lambda.lv <- try(as.matrix(sapply(U.lv, function(vec) {
                                    F.lv.pre <- as.matrix(F.lv[1:length(vec), ])
                                    l.lv.tr <- solve(t(F.lv.pre) %*% F.lv.pre) %*% t(F.lv.pre) %*% vec
                                    return(l.lv.tr)
                                })), silent = TRUE)
                                if ("try-error" %in% class(lambda.lv)) {
                                    break
                                } else {
                                    if ((r == 1) & (force %in% c(0, 2))) {
                                        lambda.lv <- t(lambda.lv)
                                    }
                                }
                            } else {
                                if (force %in% c(1, 3)) {
                                    lambda.lv <- matrix(NA, (r + 1), Ntr)
                                } else {
                                    lambda.lv <- matrix(NA, r, Ntr)
                                }
                                test <- try(
                                    for (i.tr in 1:Ntr) {
                                        F.lv.pre <- as.matrix(F.hat[setdiff(time.pre[[i.tr]], lv), ])
                                        lambda.lv[, i.tr] <- solve(t(F.lv.pre) %*% F.lv.pre) %*%
                                            t(F.lv.pre) %*% as.matrix(U.lv[[i.tr]])
                                    },
                                    silent = TRUE
                                )
                                if ("try-error" %in% class(test)) {
                                    break
                                }
                            }
                        }
                        lambda.lv <- t(lambda.lv)
                        e <- U.tr[which(time == lv), ] - c(F.hat[which(time == lv), ] %*% t(lambda.lv))
                    }
                    if (sameT0 == FALSE | 0 %in% I.tr) {
                        e <- e[which(pre[which(time == lv), ] == TRUE)]
                    }
                    ## accumulate residuals
                    resid_all <- c(resid_all, e)
                } ## end of leave-one-out

                if (length(resid_all) == 0) {
                    scores <- c(MSPE = Inf, WMSPE = Inf, GMSPE = Inf, WGMSPE = Inf,
                                MAD = Inf, Moment = Inf, GMoment = Inf, RMSE = Inf, Bias = Inf)
                } else {
                    ## Build time indices for LOO residuals (CFE)
                    time_idx_loo <- NULL
                    obs_w_loo <- NULL
                    if (!is.null(count.T.cv)) {
                        time_idx_loo <- c()
                        for (lv in unique(unlist(time.pre))) {
                            if (sameT0 == FALSE | 0 %in% I.tr) {
                                n_resid_lv <- sum(pre[which(time == lv), ] == TRUE)
                            } else {
                                n_resid_lv <- Ntr
                            }
                            if (n_resid_lv > 0) {
                                t.on.tr <- T.on[, tr, drop = FALSE]
                                t.on.lv <- unique(t.on.tr[lv, ])
                                t.on.lv <- t.on.lv[!is.na(t.on.lv)]
                                if (length(t.on.lv) > 0) {
                                    time_idx_loo <- c(time_idx_loo, rep(as.character(t.on.lv[1]), n_resid_lv))
                                } else {
                                    time_idx_loo <- c(time_idx_loo, rep("Control", n_resid_lv))
                                }
                            }
                        }
                    }
                    if (!is.null(W.tr)) {
                        obs_w_loo <- c()
                        for (lv in unique(unlist(time.pre))) {
                            if (sameT0 == FALSE | 0 %in% I.tr) {
                                w_lv <- W.tr[lv, which(pre[which(time == lv), ] == TRUE)]
                            } else {
                                w_lv <- W.tr[lv, ]
                            }
                            obs_w_loo <- c(obs_w_loo, w_lv)
                        }
                    }
                    scores <- .score_residuals(
                        resid_all,
                        obs_weights = obs_w_loo,
                        time_index = time_idx_loo,
                        count_weights = count.T.cv,
                        norm.para = norm.para
                    )
                }

              } else if (cv.method %in% c("all_units", "rolling")) {
                ## ---- cv.sample "all_units" / "rolling" CFE CV (serial path — lapply only) ---- ##
                ## Rolling reuses the all_units CFE scoring helper; only rmCV/estCV
                ## fold construction differs (built via .build_cv_mask_rolling above).
                fold_results <- lapply(1:k, function(ii) {
                    .fect_cv_score_one_cfe_nt_all(
                        ii              = ii,
                        YY.co           = YY.co,
                        Y0CV.co         = Y0CV.co,
                        X.co            = X.co,
                        II.co           = II.co,
                        W.use           = W.use,
                        W               = W,
                        beta0CV.co      = beta0CV.co,
                        X.extra.FE.co.B = X.extra.FE.co.B,
                        X.Z.co          = X.Z.co,
                        X.Q.co          = X.Q.co,
                        X.gamma.co      = X.gamma.co,
                        X.kappa.co      = X.kappa.co,
                        Zgamma.id       = Zgamma.id,
                        kappaQ.id       = kappaQ.id,
                        rmCV            = rmCV,
                        estCV           = estCV,
                        r               = r,
                        force           = force,
                        cv_tol          = cv_tol,
                        max.iteration   = max.iteration
                    )
                })
                all_resid    <- unlist(lapply(fold_results, `[[`, "resid"))
                all_time_idx <- unlist(lapply(fold_results, `[[`, "time_idx"))
                all_obs_w    <- if (!is.null(W)) unlist(lapply(fold_results, `[[`, "obs_w")) else c()
                if (length(all_resid) == 0) {
                    scores <- c(MSPE = Inf, WMSPE = Inf, GMSPE = Inf, WGMSPE = Inf,
                                MAD = Inf, Moment = Inf, GMoment = Inf, RMSE = Inf, Bias = Inf)
                } else {
                    scores <- .score_residuals(
                        all_resid,
                        obs_weights = if (!is.null(W)) all_obs_w else NULL,
                        time_index = all_time_idx,
                        count_weights = count.T.cv,
                        norm.para = NULL
                    )
                }

              } else {
                ## ---- cv.sample "treated_units" CFE CV (serial path — lapply only) ---- ##
                ## Build U.tr: subtract Layer 1 from Y.tr
                if (r != 0) {
                    F.hat <- as.matrix(est.co$factor)
                    if (force %in% c(1, 3)) {
                        F.hat <- cbind(F.hat, rep(1, TT))
                    }
                }
                U.tr <- Y.tr
                if (p > 0) {
                    for (j in 1:p) {
                        U.tr <- U.tr - X.tr[, , j] * beta[j]
                    }
                }
                if (force != 0) {
                    U.tr <- U.tr - matrix(est.co$mu, TT, Ntr)
                }
                if (force %in% c(2, 3)) {
                    U.tr <- U.tr - matrix(est.co$xi, TT, Ntr, byrow = FALSE)
                }

                ## Subtract gamma for treated (Layer 1)
                if (!is.null(est.co$gamma) && length(est.co$gamma) > 0) {
                    for (k_g in seq_along(est.co$gamma)) {
                        gamma.fit.tr.k <- .reconstruct_gamma_fit_tr(
                            est.co$gamma[[k_g]], X.Z.tr, X.gamma.tr[, , k_g, drop = FALSE],
                            Zgamma.id[[k_g]], TT, Ntr)
                        U.tr <- U.tr - gamma.fit.tr.k
                    }
                }

                ## Subtract Type-B extra FE for treated (Layer 1)
                if (length(typeB_idx) > 0) {
                    typeB.fit.tr <- .extract_and_apply_typeB_fe(
                        est.co, X.co, X.extra.FE.co.B, X.extra.FE.tr,
                        typeB_idx, X.Z.co, X.gamma.co, X.Q.co, X.kappa.co,
                        Zgamma.id, kappaQ.id,
                        TT, Nco, Ntr, p, r, force)
                    U.tr <- U.tr - typeB.fit.tr
                }

                if (0 %in% I.tr) {
                    U.tr[which(I.tr == 0)] <- 0
                }

                fold_results <- lapply(1:k, function(ii) {
                    .fect_cv_score_one_cfe_nt_tr(
                        ii       = ii,
                        U.tr     = U.tr,
                        F.hat    = F.hat,
                        pre      = pre,
                        r        = r,
                        force    = force,
                        rmCV.tr  = rmCV.tr,
                        estCV.tr = estCV.tr,
                        W.tr     = W.tr,
                        T.on     = T.on,
                        tr       = tr,
                        TT       = TT,
                        Ntr      = Ntr
                    )
                })
                all_resid    <- unlist(lapply(fold_results, `[[`, "resid"))
                all_time_idx <- unlist(lapply(fold_results, `[[`, "time_idx"))
                all_obs_w    <- if (!is.null(W.tr)) unlist(lapply(fold_results, `[[`, "obs_w")) else c()
                if (length(all_resid) == 0) {
                    scores <- c(MSPE = Inf, WMSPE = Inf, GMSPE = Inf, WGMSPE = Inf,
                                MAD = Inf, Moment = Inf, GMoment = Inf, RMSE = Inf, Bias = Inf)
                } else {
                    scores <- .score_residuals(
                        all_resid,
                        obs_weights = if (!is.null(W.tr)) all_obs_w else NULL,
                        time_index = if (length(all_time_idx) > 0) all_time_idx else NULL,
                        count_weights = count.T.cv,
                        norm.para = norm.para
                    )
                }

              } ## end cv.method branching

                if ((min(CV.out[, crit_col]) - scores[crit_col]) > 0.01 * min(CV.out[, crit_col])) {
                    est.co.best <- est.co
                    r.cv <- r
                } else {
                    if (r == r.cv + 1) message("*")
                }

                if (PC < min(CV.out[, "PC"])) {
                    r.pc <- r
                    est.co.pc.best <- est.co
                }
                CV.out[i, 2:4] <- c(sigma2, IC, PC)
                CV.out[i, score_names] <- scores[score_names]
                message("r = ", r, "; sigma2 = ",
                    sprintf("%.5f", sigma2), "; IC = ",
                    sprintf("%.5f", IC), "; PC = ",
                    sprintf("%.5f", PC), "; MSPE = ",
                    sprintf("%.5f", scores["MSPE"]),
                    sep = ""
                )
            } ## end of CV loop

            } ## end SERIAL BRANCH (CFE)

            MSPE.best <- min(CV.out[, "MSPE"])
            if (r > (T0.min - 1)) {
                message(" (r hits maximum)")
            }
            message("\n r* = ", r.cv, sep = "")
            message("\n")
        }
    } else {
        r.cv <- r
        r.min <- r.max <- r
    }

    ## ---- Final fit with selected r.cv ----
    est.co.best <- complex_fe_ub(YY.co, Y0.co, X.co,
        X.extra.FE.co.B, X.Z.co, X.Q.co, X.gamma.co, X.kappa.co,
        Zgamma.id, kappaQ.id,
        II.co, W.use, beta0, r.cv, force = force, tol, max.iteration)

    ## Convergence check
    if (!is.null(est.co.best$niter) && est.co.best$niter >= max.iteration) {
        warning(paste0("CFE optimization did not converge within ", max.iteration,
            " iterations. Results may be unreliable."))
    }

    est.co.fect <- NULL
    if (boot == FALSE) {
        if (r.cv == 0) {
            est.co.fect <- est.co.best
        } else {
            est.co.fect <- complex_fe_ub(YY.co, Y0.co, X.co,
                X.extra.FE.co.B, X.Z.co, X.Q.co, X.gamma.co, X.kappa.co,
                Zgamma.id, kappaQ.id,
                II.co, W.use, beta0, 0, force = force, tol, max.iteration)
        }
    }

    validX <- est.co.best$validX
    validF <- ifelse(r.cv > 0, 1, 0)

    ## ---- Three-Layer Projection ----
    ## Layer 1: subtract shared parameters from Y.tr
    U.tr.r0 <- U.tr <- Y.tr
    if (p > 0) {
        beta <- est.co.best$beta
        if (est.co.best$validX == 0) {
            beta <- matrix(0, p, 1)
        } else {
            beta <- est.co.best$beta
            beta[is.nan(est.co.best$beta)] <- 0
        }
        for (j in 1:p) {
            U.tr <- U.tr - X.tr[, , j] * beta[j]
        }
        if (boot == FALSE) {
            beta.r0 <- est.co.fect$beta
            if (est.co.fect$validX == 0) {
                beta.r0 <- matrix(0, p, 1)
            } else {
                beta.r0 <- est.co.fect$beta
                beta.r0[is.nan(est.co.fect$beta)] <- 0
            }
            for (j in 1:p) {
                U.tr.r0 <- U.tr.r0 - X.tr[, , j] * beta.r0[j]
            }
        }
    } else {
        beta <- NA
        beta.r0 <- NA
    }

    mu <- est.co.best$mu
    U.tr <- U.tr - matrix(mu, TT, Ntr)
    Y.fe.bar <- rep(mu, TT)

    if (boot == FALSE) {
        mu.r0 <- est.co.fect$mu
        U.tr.r0 <- U.tr.r0 - matrix(mu.r0, TT, Ntr)
        Y.fe.bar.r0 <- rep(mu.r0, TT)
    }

    if (force %in% c(2, 3)) {
        xi <- est.co.best$xi
        U.tr <- U.tr - matrix(c(xi), TT, Ntr, byrow = FALSE)
        Y.fe.bar <- Y.fe.bar + xi
        if (boot == FALSE) {
            xi.r0 <- est.co.fect$xi
            U.tr.r0 <- U.tr.r0 - matrix(c(xi.r0), TT, Ntr, byrow = FALSE)
            Y.fe.bar.r0 <- Y.fe.bar.r0 + xi.r0
        }
    }

    ## Subtract gamma for treated (Layer 1)
    if (!is.null(est.co.best$gamma) && length(est.co.best$gamma) > 0) {
        for (k_g in seq_along(est.co.best$gamma)) {
            gamma.fit.tr.k <- .reconstruct_gamma_fit_tr(
                est.co.best$gamma[[k_g]], X.Z.tr, X.gamma.tr[, , k_g, drop = FALSE],
                Zgamma.id[[k_g]], TT, Ntr)
            U.tr <- U.tr - gamma.fit.tr.k
        }
    }

    ## Subtract Type-B extra FE for treated (Layer 1)
    if (length(typeB_idx) > 0) {
        typeB.fit.tr <- .extract_and_apply_typeB_fe(
            est.co.best, X.co, X.extra.FE.co.B, X.extra.FE.tr,
            typeB_idx, X.Z.co, X.gamma.co, X.Q.co, X.kappa.co,
            Zgamma.id, kappaQ.id,
            TT, Nco, Ntr, p, r.cv, force)
        U.tr <- U.tr - typeB.fit.tr
    }

    ## Layer 2: estimate alpha, kappa, Type-A FE, lambda from pre-treatment
    ## Save Layer 1 residuals (before any Layer 2 subtractions)
    U.tr.L1 <- U.tr

    has_kappa <- !is.null(est.co.best$kappa) && length(est.co.best$kappa) > 0

    ## --- Helper: compute kappa_fit (TT x Ntr) from current residuals ---
    .estimate_kappa_fit <- function(U.cur) {
        kappa_fit <- matrix(0, TT, Ntr)
        if (!has_kappa) return(kappa_fit)
        for (k_k in seq_along(est.co.best$kappa)) {
            q_cols <- kappaQ.id[[k_k]]
            kappa_groups <- X.kappa.tr[1, , k_k]
            unique_kgroups <- sort(unique(kappa_groups))
            for (g_idx in seq_along(unique_kgroups)) {
                g <- unique_kgroups[g_idx]
                units_in_group <- which(kappa_groups == g)
                Q.pre.list <- list()
                U.pre.list <- list()
                for (i_idx in seq_along(units_in_group)) {
                    ii <- units_in_group[i_idx]
                    pre_t <- which(pre[, ii])
                    if (length(pre_t) > 0) {
                        Q.mat <- matrix(0, length(pre_t), length(q_cols))
                        for (jj in seq_along(q_cols)) {
                            Q.mat[, jj] <- X.Q.tr[pre_t, ii, q_cols[jj]]
                        }
                        Q.pre.list[[length(Q.pre.list) + 1]] <- Q.mat
                        U.pre.list[[length(U.pre.list) + 1]] <- U.cur[pre_t, ii]
                    }
                }
                if (length(U.pre.list) > 0) {
                    Q.pre.all <- do.call(rbind, Q.pre.list)
                    U.pre.all <- unlist(U.pre.list)
                    if (length(U.pre.all) > length(q_cols)) {
                        kappa.hat <- try(solve(t(Q.pre.all) %*% Q.pre.all) %*%
                            t(Q.pre.all) %*% U.pre.all, silent = TRUE)
                        if (!"try-error" %in% class(kappa.hat)) {
                            for (ii in units_in_group) {
                                Q.full <- matrix(0, TT, length(q_cols))
                                for (jj in seq_along(q_cols)) {
                                    Q.full[, jj] <- X.Q.tr[, ii, q_cols[jj]]
                                }
                                kappa_fit[, ii] <- kappa_fit[, ii] + Q.full %*% kappa.hat
                            }
                        }
                    }
                }
            }
        }
        return(kappa_fit)
    }

    ## --- Helper: compute alpha (Ntr x 1) from current residuals ---
    .estimate_alpha <- function(U.cur) {
        if (!force %in% c(1, 3)) return(matrix(0, Ntr, 1))
        if (max(T0) == T0.min & (!0 %in% I.tr)) {
            return(as.matrix(colMeans(U.cur[1:T0.min, ])))
        } else {
            U.pre.v <- as.vector(U.cur)[which(pre.v == 1)]
            U.pre.l <- split(U.pre.v, id.tr.pre.v)
            return(as.matrix(sapply(U.pre.l, mean)))
        }
    }

    ## --- Helper: estimate Type-A extra FE fit (TT x Ntr) from residuals ---
    .estimate_typeA_fit <- function(U.cur) {
        fe_fit <- matrix(0, TT, Ntr)
        if (length(typeA_idx) == 0) return(fe_fit)
        for (k_a in typeA_idx) {
            labels.tr <- X.extra.FE.tr[1, , k_a]
            unique_levels <- sort(unique(labels.tr))
            for (g in unique_levels) {
                units_in_level <- which(labels.tr == g)
                pre_vals <- c()
                for (ii in units_in_level) {
                    pre_t <- which(pre[, ii])
                    pre_vals <- c(pre_vals, U.cur[pre_t, ii])
                }
                if (length(pre_vals) > 0) {
                    fe_mean <- mean(pre_vals)
                    for (ii in units_in_level) {
                        fe_fit[, ii] <- fe_fit[, ii] + fe_mean
                    }
                }
            }
        }
        return(fe_fit)
    }

    ## --- Helper: estimate lambda fit (TT x Ntr) from residuals ---
    ## When force %in% c(1,3), F.hat.aug includes intercept column;
    ## alpha is embedded as the last column of lambda.
    ## Returns list(fit = TT x Ntr, lambda = Ntr x ncol(F.hat.aug))
    .estimate_lambda_fit <- function(U.cur, F.hat.aug) {
        ncol_f <- ncol(F.hat.aug)
        if (max(T0) == T0.min & (!0 %in% I.tr)) {
            F.pre <- F.hat.aug[1:T0.min, , drop = FALSE]
            U.pre <- as.matrix(U.cur[1:T0.min, ])
            lam <- try(solve(t(F.pre) %*% F.pre) %*% t(F.pre) %*% U.pre,
                       silent = TRUE)
            if ("try-error" %in% class(lam)) {
                return(list(fit = matrix(0, TT, Ntr),
                            lambda = matrix(0, Ntr, ncol_f), ok = FALSE))
            }
        } else if (!0 %in% I.tr) {
            lam <- try(as.matrix(sapply(seq_len(Ntr), function(j) {
                pre_t <- which(pre[, j])
                F.pre <- as.matrix(F.hat.aug[pre_t, , drop = FALSE])
                solve(t(F.pre) %*% F.pre) %*% t(F.pre) %*% U.cur[pre_t, j]
            })), silent = TRUE)
            if ("try-error" %in% class(lam)) {
                return(list(fit = matrix(0, TT, Ntr),
                            lambda = matrix(0, Ntr, ncol_f), ok = FALSE))
            }
            if (ncol_f == 1) lam <- t(lam)
        } else {
            lam <- matrix(NA, ncol_f, Ntr)
            test <- try(
                for (i.tr in 1:Ntr) {
                    F.pre <- as.matrix(F.hat.aug[time.pre[[i.tr]], , drop = FALSE])
                    lam[, i.tr] <- solve(t(F.pre) %*% F.pre) %*%
                        t(F.pre) %*% as.matrix(U.cur[time.pre[[i.tr]], i.tr])
                }, silent = TRUE)
            if ("try-error" %in% class(test)) {
                return(list(fit = matrix(0, TT, Ntr),
                            lambda = matrix(0, Ntr, ncol_f), ok = FALSE))
            }
        }
        lam_mat <- t(lam)  ## Ntr x ncol_f
        fit <- F.hat.aug %*% t(lam_mat)  ## TT x Ntr
        return(list(fit = fit, lambda = lam_mat, ok = TRUE))
    }

    ## ============================================================
    ## Block coordinate descent: jointly estimate all unit-specific
    ## parameters (alpha, kappa, Type-A FE, lambda) from
    ## pre-treatment residuals. Mirrors C++ cfe_iter logic.
    ## ============================================================
    has_alpha  <- force %in% c(1, 3)
    has_typeA  <- length(typeA_idx) > 0
    has_factor <- r.cv > 0

    ## Build augmented factor matrix (factors + intercept for alpha)
    F.hat.aug <- NULL
    if (has_factor) {
        F.hat.aug <- as.matrix(est.co.best$factor)
        if (has_alpha) F.hat.aug <- cbind(F.hat.aug, rep(1, TT))
    }

    ## Initialize all component fits to zero
    ## NOTE: When has_factor && has_alpha, alpha is embedded in lambda_fit
    ## (via the augmented intercept column in F.hat.aug). In that case,
    ## alpha_mat is NOT used in residual computation to avoid double subtraction.
    ## alpha_mat is only used when has_alpha && !has_factor.
    alpha_mat  <- matrix(0, TT, Ntr)   ## alpha broadcast to TT x Ntr
    kappa_fit  <- matrix(0, TT, Ntr)
    typeA_fit  <- matrix(0, TT, Ntr)
    lambda_fit <- matrix(0, TT, Ntr)   ## F %*% t(lambda), includes alpha when augmented
    alpha.tr   <- matrix(0, Ntr, 1)
    lambda.tr  <- NULL
    lambda.co  <- NULL

    ## When has_factor && has_alpha, alpha lives inside lambda_fit.
    ## alpha_mat is separate only when !has_factor.
    alpha_in_lambda <- has_factor && has_alpha

    n_components <- has_kappa + has_typeA + (has_factor || has_alpha)

    if (n_components >= 2) {
        ## Multiple components: iterate to convergence
        max_iter_bcd <- 100
        tol_bcd <- 1e-8

        for (iter_bcd in 1:max_iter_bcd) {
            old_kappa  <- kappa_fit
            old_typeA  <- typeA_fit
            old_alpha  <- alpha_mat
            old_lambda <- lambda_fit

            ## Step 1: estimate kappa from residual
            if (has_kappa) {
                if (alpha_in_lambda) {
                    resid <- U.tr.L1 - typeA_fit - lambda_fit
                } else {
                    resid <- U.tr.L1 - alpha_mat - typeA_fit - lambda_fit
                }
                kappa_fit <- .estimate_kappa_fit(resid)
            }

            ## Step 2: estimate Type-A FE from residual
            if (has_typeA) {
                if (alpha_in_lambda) {
                    resid <- U.tr.L1 - kappa_fit - lambda_fit
                } else {
                    resid <- U.tr.L1 - alpha_mat - kappa_fit - lambda_fit
                }
                typeA_fit <- .estimate_typeA_fit(resid)
            }

            ## Step 3: estimate alpha (+lambda if r.cv > 0)
            if (has_factor) {
                ## alpha embedded in augmented F.hat → lambda_fit includes alpha
                resid <- U.tr.L1 - kappa_fit - typeA_fit
                result <- .estimate_lambda_fit(resid, F.hat.aug)
                if (!result$ok) {
                    return(list(att = rep(NA, TT), att.avg = NA,
                                beta = matrix(NA, p, 1)))
                }
                lambda_fit <- result$fit
                lambda.tr <- result$lambda
                if (has_alpha) {
                    alpha.tr <- as.matrix(lambda.tr[, ncol(F.hat.aug), drop = FALSE])
                    ## Do NOT update alpha_mat — alpha is inside lambda_fit
                }
            } else if (has_alpha) {
                resid <- U.tr.L1 - kappa_fit - typeA_fit
                alpha.tr <- .estimate_alpha(resid)
                alpha_mat <- matrix(alpha.tr, TT, Ntr, byrow = TRUE)
            }

            ## Check convergence of all components
            delta <- max(
                max(abs(kappa_fit  - old_kappa)),
                max(abs(typeA_fit  - old_typeA)),
                if (alpha_in_lambda) 0 else max(abs(alpha_mat - old_alpha)),
                max(abs(lambda_fit - old_lambda))
            )
            if (delta < tol_bcd) break
        }
    } else {
        ## Single component (or none): one-pass estimation
        if (has_kappa) {
            kappa_fit <- .estimate_kappa_fit(U.tr.L1)
        }
        if (has_factor) {
            resid <- U.tr.L1 - kappa_fit - typeA_fit
            result <- .estimate_lambda_fit(resid, F.hat.aug)
            if (!result$ok) {
                return(list(att = rep(NA, TT), att.avg = NA,
                            beta = matrix(NA, p, 1)))
            }
            lambda_fit <- result$fit
            lambda.tr <- result$lambda
            if (has_alpha) {
                alpha.tr <- as.matrix(lambda.tr[, ncol(F.hat.aug), drop = FALSE])
            }
        } else if (has_alpha) {
            resid <- U.tr.L1 - kappa_fit - typeA_fit
            alpha.tr <- .estimate_alpha(resid)
            alpha_mat <- matrix(alpha.tr, TT, Ntr, byrow = TRUE)
        }
        if (has_typeA) {
            if (alpha_in_lambda) {
                resid <- U.tr.L1 - kappa_fit - lambda_fit
            } else {
                resid <- U.tr.L1 - alpha_mat - kappa_fit - lambda_fit
            }
            typeA_fit <- .estimate_typeA_fit(resid)
        }
    }

    ## Final residual = treatment effect
    if (alpha_in_lambda) {
        ## alpha is inside lambda_fit — don't subtract alpha_mat
        U.tr <- U.tr.L1 - kappa_fit - typeA_fit - lambda_fit
    } else {
        U.tr <- U.tr.L1 - alpha_mat - kappa_fit - typeA_fit - lambda_fit
    }
    eff <- U.tr

    ## Extract clean lambda.tr (remove alpha column if embedded)
    if (has_factor && has_alpha && !is.null(lambda.tr)) {
        alpha.tr <- as.matrix(lambda.tr[, ncol(F.hat.aug), drop = FALSE])
        lambda.tr <- lambda.tr[, 1:r.cv, drop = FALSE]
    } else if (has_factor && !is.null(lambda.tr)) {
        ## lambda.tr already clean
    } else {
        lambda.tr <- NULL
    }
    lambda.co <- if (has_factor) est.co.best$lambda else NULL

    ## Implied weights
    if (has_factor && boot == 0 && !is.null(lambda.tr)) {
        inv.tr <- try(ginv(t(as.matrix(lambda.tr))), silent = TRUE)
        if (!"try-error" %in% class(inv.tr)) {
            wgt.implied <- t(inv.tr %*% t(as.matrix(est.co.best$lambda)))
        }
    }

    ## r=0 path (for equivalence test baseline — uses FE-only model)
    if (force %in% c(1, 3)) {
        if (boot == FALSE) {
            if ((max(T0) == T0.min) & (!0 %in% I.tr)) {
                alpha.tr.r0 <- as.matrix(colMeans(as.matrix(U.tr.r0[1:T0.min, ])))
            } else {
                U.tr.pre.v.r0 <- as.vector(U.tr.r0)[which(pre.v == 1)]
                U.tr.pre.r0 <- split(U.tr.pre.v.r0, id.tr.pre.v)
                alpha.tr.r0 <- as.matrix(sapply(U.tr.pre.r0, mean))
            }
            U.tr.r0 <- U.tr.r0 - matrix(alpha.tr.r0, TT, Ntr, byrow = TRUE)
        }
    }
    if (boot == FALSE) {
        eff.r0 <- U.tr.r0
    }

    if (0 %in% I.tr) {
        eff[which(I.tr == 0)] <- 0
        if (boot == FALSE) {
            eff.r0[which(I.tr == 0)] <- 0
        }
    }

  } ## end of method bifurcation

    ## -------------------------------##
    ## Summarize
    ## -------------------------------##

    ## counterfactuals for treated units
    Y.ct.tr <- as.matrix(Y.tr - eff)
    Y.ct.co <- Y.co - est.co.best$residuals
    # print(Y.ct.co)
    Y.ct <- Y
    Y.ct[, tr] <- Y.ct.tr
    Y.ct[, co] <- Y.ct.co

    if (boot == FALSE) {
        Y.ct.tr.r0 <- as.matrix(Y.tr - eff.r0)
        Y.ct.co.r0 <- Y.co - est.co.fect$residuals
        Y.ct.r0 <- Y
        Y.ct.r0[, tr] <- Y.ct.tr.r0
        Y.ct.r0[, co] <- Y.ct.co.r0
    }

    ## we first adjustment for normalization
    if (!is.null(norm.para)) {
        Y <- Y * norm.para[1]
        ## variance of the error term
        sigma2 <- est.co.best$sigma2 <- est.co.best$sigma2 * (norm.para[1]^2)
        IC <- est.co.best$IC <- est.co.best$IC - log(est.co.best$sigma2) + log(sigma2)
        PC <- est.co.best$PC <- est.co.best$PC * (norm.para[1]^2)

        ## output of estimates
        mu <- est.co.best$mu <- est.co.best$mu * norm.para[1]
        if (r.cv > 0) {
            est.co.best$lambda <- est.co.best$lambda * norm.para[1]
            lambda.tr <- lambda.tr * norm.para[1]
            est.co.best$VNT <- est.co.best$VNT * norm.para[1]
        }
        if (force %in% c(1, 3)) {
            est.co.best$alpha <- est.co.best$alpha * norm.para[1]
            alpha.tr <- alpha.tr * norm.para[1]
        }
        if (force %in% c(2, 3)) {
            xi <- est.co.best$xi <- est.co.best$xi * norm.para[1]
        }
        est.co.best$residuals <- est.co.best$residuals * norm.para[1]
        est.co.best$fit <- est.co.best$fit * norm.para[1]
        if (boot == FALSE) {
            est.co.fect$fit <- est.co.fect$fit * norm.para[1]
        }
        est.co.fect$sigma2 <- est.co.fect$sigma2 * norm.para[1]

        # estimated counterfactual
        Y.tr <- Y.tr * norm.para[1]
        Y.co <- Y.co * norm.para[1]

        Y.ct <- Y.ct * norm.para[1]
        Y.ct.tr <- Y.ct.tr * norm.para[1]
        Y.ct.co <- Y.ct.co * norm.para[1]
        eff <- eff * norm.para[1]

        if (boot == FALSE) {
            Y.ct.r0 <- Y.ct.r0 * norm.para[1]
            Y.ct.tr.r0 <- Y.ct.tr.r0 * norm.para[1]
            Y.ct.co.r0 <- Y.ct.co.r0 * norm.para[1]
            eff.r0 <- eff.r0 * norm.para[1]
        }
    }

    ## 0. relevant parameters
    IC <- est.co.best$IC
    sigma2 <- est.co.best$sigma2
    PC <- est.co.best$PC
    loglikelihood <- NULL

    if (p > 0) {
        na.pos <- is.nan(est.co.best$beta)
        beta <- est.co.best$beta
        if (sum(na.pos) > 0) {
            beta[na.pos] <- NA
        }
    } else {
        beta <- NA
    }

    ## 1. estimated att and counterfactuals
    if (boot == FALSE) {
        Y.ct.equiv <- Y.ct.r0
    } else {
        Y.ct.equiv <- NULL
    }

    eff <- Y - Y.ct

    missing.index <- which(is.na(eff))
    if (length(missing.index) > 0) {
        I[missing.index] <- 0
        II[missing.index] <- 0
    }
    if (0 %in% I) {
        eff[which(I == 0)] <- NA
    }
    complete.index <- which(!is.na(eff))
    att.avg <- sum(eff[complete.index] * D[complete.index]) / (sum(D[complete.index]))
    marginal <- NULL

    att.avg.balance <- NA
    if (!is.null(balance.period)) {
        complete.index2 <- which(!is.na(T.on.balance))
        att.avg.balance <- sum(eff[complete.index2] * D[complete.index2]) / (sum(D[complete.index2]))
    }

    # weighted effect
    att.avg.W <- NA
    if (!is.null(W)) {
        att.avg.W <- sum(eff[complete.index] * D[complete.index] * W[complete.index]) / (sum(D[complete.index] * W[complete.index]))
    }

    ## att.avg.unit
    tr.pos <- which(apply(D, 2, sum) > 0)
    att.unit <- sapply(1:length(tr.pos), function(vec) {
        return(sum(eff[, tr.pos[vec]] * D[, tr.pos[vec]]) / sum(D[, tr.pos[vec]]))
    })
    att.avg.unit <- mean(att.unit, na.rm = TRUE)

    equiv.att.avg <- eff.equiv <- NULL
    if (boot == FALSE) {
        eff.equiv <- Y - Y.ct.equiv
        if (0 %in% I) {
            eff.equiv[which(I == 0)] <- NA
        }
        complete.index <- which(!is.na(eff.equiv))
        equiv.att.avg <- sum(eff.equiv[complete.index] * D[complete.index]) / (sum(D[complete.index]))
    }

    ## 2. rmse for treated units' observations under control
    if (binary == 0) {
        tr <- which(apply(D, 2, sum) > 0)
        tr.co <- which((as.matrix(1 - D[, tr]) * as.matrix(II[, tr])) == 1)
        eff.tr <- as.matrix(eff[, tr])
        v.eff.tr <- eff.tr[tr.co]
        rmse <- sqrt(mean(v.eff.tr^2, na.rm = TRUE))
    }

    ## 3. unbalanced output
    Y.ct.full <- Y.ct
    ## Stage 2: populate Y.ct.full at control positions on the GSC path.
    ## On the IFE-EM (notyettreated) path EM imputation already fills in
    ## masked control positions, so user-space rolling CV can read off
    ## fit$Y.ct.full there. On GSC the residual recipe `Y.co - residuals`
    ## leaves NA at masked positions because Y.co is NA. Overwriting with
    ## the model-implied F * t(lambda_co) closes that gap. ATT, gap, and
    ## est.avg are computed from treated units only and are unaffected.
    if (method == "ife" &&
        !is.null(est.co.best$factor) &&
        !is.null(est.co.best$lambda)) {
        F.hat.ctf     <- as.matrix(est.co.best$factor)
        Lambda.co.ctf <- as.matrix(est.co.best$lambda)
        if (ncol(F.hat.ctf) > 0 &&
            ncol(F.hat.ctf) == ncol(Lambda.co.ctf) &&
            nrow(F.hat.ctf) == TT && nrow(Lambda.co.ctf) == Nco) {
            Y.ct.full[, co] <- F.hat.ctf %*% t(Lambda.co.ctf)
        }
    }
    res.full <- Y - Y.ct
    if (0 %in% I) {
        eff[which(I == 0)] <- NA
        Y.ct[which(I == 0)] <- NA
    }
    if (binary == FALSE) {
        res.full[which(II == 0)] <- NA
    }


    ## 4. dynamic effects
    t.on <- c(T.on)
    eff.v <- c(eff) ## a vector
    eff.equiv.v <- NULL
    if (binary == FALSE && boot == FALSE) {
        eff.equiv.v <- c(eff.equiv)
    }

    rm.pos1 <- which(is.na(eff.v))
    rm.pos2 <- which(is.na(t.on))

    eff.v.use1 <- eff.v
    t.on.use <- t.on
    n.on.use <- rep(1:N, each = TT)

    if (NA %in% eff.v | NA %in% t.on) {
        eff.v.use1 <- eff.v[-c(rm.pos1, rm.pos2)]
        t.on.use <- t.on[-c(rm.pos1, rm.pos2)]
        n.on.use <- n.on.use[-c(rm.pos1, rm.pos2)]
        if (binary == FALSE && boot == FALSE) {
            eff.equiv.v <- eff.equiv.v[-c(rm.pos1, rm.pos2)]
        }
    }

    pre.pos <- which(t.on.use <= 0)
    eff.pre <- cbind(eff.v.use1[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
    colnames(eff.pre) <- c("eff", "period", "unit")

    pre.sd <- eff.pre.equiv <- NULL
    if (binary == FALSE && boot == FALSE) {
        eff.pre.equiv <- cbind(eff.equiv.v[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
        colnames(eff.pre.equiv) <- c("eff.equiv", "period", "unit")

        pre.sd <- tapply(eff.pre.equiv[, 1], eff.pre.equiv[, 2], sd)
        pre.sd <- cbind(pre.sd, sort(unique(eff.pre.equiv[, 2])), table(eff.pre.equiv[, 2]))
        colnames(pre.sd) <- c("sd", "period", "count")
    }

    time.on <- sort(unique(t.on.use))
    att.on <- as.numeric(tapply(eff.v.use1, t.on.use, mean)) ## NA already removed
    count.on <- as.numeric(table(t.on.use))

    if (!is.null(time.on.seq)) {
        count.on.med <- att.on.med <- rep(NA, length(time.on.seq))
        att.on.med[which(time.on.seq %in% time.on)] <- att.on
        count.on.med[which(time.on.seq %in% time.on)] <- count.on
        att.on <- att.on.med
        count.on <- count.on.med
        time.on <- time.on.seq
    }

    if (!is.null(W)) {
        W.v <- c(W)
        rm.pos.W <- which(is.na(W))
        if (NA %in% eff.v | NA %in% t.on | NA %in% W.v) {
            eff.v.use.W <- eff.v[-c(rm.pos1, rm.pos2, rm.pos.W)]
            W.v.use <- W.v[-c(rm.pos1, rm.pos2, rm.pos.W)]
            t.on.use.W <- t.on[-c(rm.pos1, rm.pos2, rm.pos.W)]
            n.on.use.W <- n.on.use[-c(rm.pos1, rm.pos2, rm.pos.W)]
        } else {
            eff.v.use.W <- eff.v.use1
            t.on.use.W <- t.on.use
            n.on.use.W <- n.on.use
            W.v.use <- W.v
        }
        time.on.W <- sort(unique(t.on.use.W))
        att.on.sum.W <- as.numeric(tapply(eff.v.use.W * W.v.use, t.on.use.W, sum)) ## NA already removed
        W.on.sum <- as.numeric(tapply(W.v.use, t.on.use.W, sum))
        att.on.W <- att.on.sum.W / W.on.sum
        count.on.W <- as.numeric(table(t.on.use.W))

        if (!is.null(time.on.seq.W)) {
            att.on.sum.med.W <- W.on.sum.med <- count.on.med.W <- att.on.med.W <- rep(NA, length(time.on.seq.W))
            att.on.sum.med.W[which(time.on.seq.W %in% time.on.W)] <- att.on.sum.W
            att.on.med.W[which(time.on.seq.W %in% time.on.W)] <- att.on.W
            count.on.med.W[which(time.on.seq.W %in% time.on.W)] <- count.on.W
            W.on.sum.med[which(time.on.seq.W %in% time.on.W)] <- W.on.sum
            att.on.sum.W <- att.on.sum.med.W
            att.on.W <- att.on.med.W
            count.on.W <- count.on.med.W
            time.on.W <- time.on.seq.W
            W.on.sum <- W.on.sum.med
        }
    } else {
        att.on.sum.med.W <- att.on.sum.W <- count.on.med.W <- att.on.med.W <- W.on.sum.med <- att.on.W <- count.on.W <- time.on.W <- W.on.sum <- NULL
    }

    ## 4.2 balance effect
    balance.att <- NULL
    if (!is.null(balance.period)) {
        t.on.balance <- c(T.on.balance)
        rm.pos4 <- which(is.na(t.on.balance))
        t.on.balance.use <- t.on.balance

        if (NA %in% eff.v | NA %in% t.on.balance) {
            eff.v.use3 <- eff.v[-c(rm.pos1, rm.pos4)]
            t.on.balance.use <- t.on.balance[-c(rm.pos1, rm.pos4)]
        }

        balance.time <- sort(unique(t.on.balance.use))
        balance.att <- as.numeric(tapply(eff.v.use3, t.on.balance.use, mean)) ## NA already removed
        balance.count <- as.numeric(table(t.on.balance.use))

        if (!is.null(time.on.balance.seq)) {
            balance.att.med <- rep(NA, length(time.on.balance.seq))
            balance.count.med <- rep(0, length(time.on.balance.seq))
            balance.att.med[which(time.on.balance.seq %in% balance.time)] <- balance.att
            if (length(balance.count) > 0) {
                balance.count.med[which(time.on.balance.seq %in% balance.time)] <- balance.count
            }
            balance.count <- balance.count.med
            balance.att <- balance.att.med
            balance.time <- time.on.balance.seq
        }

        # placebo for balanced samples
        if (!is.null(placebo.period) && placeboTest == 1) {
            if (length(placebo.period) == 1) {
                balance.placebo.pos <- which(balance.time == placebo.period)
                balance.att.placebo <- balance.att[balance.placebo.pos]
            } else {
                balance.placebo.pos <- which(balance.time >= placebo.period[1] & balance.time <= placebo.period[2])
                balance.att.placebo <- sum(balance.att[balance.placebo.pos] * balance.count[balance.placebo.pos]) / sum(balance.count[balance.placebo.pos])
            }
        }
    }

    ## 5. placebo effect, if placeboTest == 1
    if (!is.null(placebo.period) && placeboTest == 1) {
        if (length(placebo.period) == 1) {
            placebo.pos <- which(time.on == placebo.period)
            att.placebo <- att.on[placebo.pos]
        } else {
            placebo.pos <- which(time.on >= placebo.period[1] & time.on <= placebo.period[2])
            att.placebo <- sum(att.on[placebo.pos] * count.on[placebo.pos]) / sum(count.on[placebo.pos])
        }

        if (!is.null(W)) {
            if (length(placebo.period) == 1) {
                placebo.pos.W <- which(time.on.W == placebo.period)
                att.placebo.W <- att.on.W[placebo.pos.W]
            } else {
                placebo.pos.W <- which(time.on.W >= placebo.period[1] & time.on.W <= placebo.period[2])
                att.placebo.W <- sum(att.on.sum.W[placebo.pos.W]) / sum(W.on.sum[placebo.pos.W])
            }
        }
    }

    ## 6. switch-off effects
    eff.off.equiv <- off.sd <- eff.off <- NULL
    if (hasRevs == 1) {
        t.off <- c(T.off)
        rm.pos3 <- which(is.na(t.off))
        eff.v.use2 <- eff.v
        t.off.use <- t.off
        if (NA %in% eff.v | NA %in% t.off) {
            eff.v.use2 <- eff.v[-c(rm.pos1, rm.pos3)]
            t.off.use <- t.off[-c(rm.pos1, rm.pos3)]
        }

        off.pos <- which(t.off.use > 0)
        eff.off <- cbind(eff.v.use2[off.pos], t.off.use[off.pos], n.on.use[off.pos])
        colnames(eff.off) <- c("eff", "period", "unit")

        if (binary == FALSE && boot == FALSE) {
            eff.off.equiv <- cbind(eff.equiv.v[off.pos], t.off.use[off.pos], n.on.use[off.pos])
            colnames(eff.off.equiv) <- c("off.equiv", "period", "unit")

            off.sd <- tapply(eff.off.equiv[, 1], eff.off.equiv[, 2], sd)
            off.sd <- cbind(off.sd, sort(unique(eff.off.equiv[, 2])), table(eff.off.equiv[, 2]))
            colnames(off.sd) <- c("sd", "period", "count")
        }

        time.off <- sort(unique(t.off.use))

        att.off <- as.numeric(tapply(eff.v.use2, t.off.use, mean)) ## NA already removed
        count.off <- as.numeric(table(t.off.use))

        if (!is.null(time.off.seq)) {
            count.off.med <- att.off.med <- rep(NA, length(time.off.seq))
            att.off.med[which(time.off.seq %in% time.off)] <- att.off
            count.off.med[which(time.off.seq %in% time.off)] <- count.off
            att.off <- att.off.med
            count.off <- count.off.med
            time.off <- time.off.seq
        }

        if (!is.null(W)) {
            if (NA %in% eff.v | NA %in% t.off | NA %in% W.v) {
                eff.v.use2.W <- eff.v[-c(rm.pos1, rm.pos3, rm.pos.W)]
                W.v.use2 <- W.v[-c(rm.pos1, rm.pos3, rm.pos.W)]
                t.off.use.W <- t.off[-c(rm.pos1, rm.pos3, rm.pos.W)]
            } else {
                eff.v.use2.W <- eff.v.use2
                t.off.use.W <- t.off.use
                W.v.use2 <- W.v
            }

            time.off.W <- sort(unique(t.off.use.W))
            att.off.sum.W <- as.numeric(tapply(eff.v.use2.W * W.v.use2, t.off.use.W, sum))
            W.off.sum <- as.numeric(tapply(W.v.use2, t.off.use.W, sum))
            att.off.W <- att.off.sum.W / W.off.sum ## NA already removed
            count.off.W <- as.numeric(table(t.off.use.W))

            if (!is.null(time.off.seq.W)) {
                att.off.sum.med.W <- W.off.sum.med <- count.off.med.W <- att.off.med.W <- rep(NA, length(time.off.seq.W))
                att.off.sum.med.W[which(time.off.seq.W %in% time.off.W)] <- att.off.sum.W
                att.off.med.W[which(time.off.seq.W %in% time.off.W)] <- att.off.W
                count.off.med.W[which(time.off.seq.W %in% time.off.W)] <- count.off.W
                W.off.sum.med[which(time.off.seq.W %in% time.off.W)] <- W.off.sum
                att.off.sum.W <- att.off.sum.med.W
                att.off.W <- att.off.med.W
                count.off.W <- count.off.med.W
                time.off.W <- time.off.seq.W
                W.off.sum <- W.off.sum.med
            }
        } else {
            W.off.sum.med <- W.off.sum <- att.off.sum.W <- att.off.sum.med.W <- count.off.med.W <- att.off.med.W <- count.off.med.W <- att.off.W <- count.off.W <- time.off.W <- NULL
        }
    }

    ## 7. carryover effects
    if (!is.null(carryover.period) && carryoverTest == 1 && hasRevs == 1) {
        ## construct att.carryover
        ## eff is derived from eff.v
        ## period and Num.Units are derived from T.off
        if (length(carryover.period) == 1) {
            carryover.pos <- which(time.off == carryover.period)
            att.carryover <- att.off[carryover.pos]
        } else {
            carryover.pos <- which(time.off >= carryover.period[1] & time.off <= carryover.period[2])
            att.carryover <- sum(att.off[carryover.pos] * count.off[carryover.pos]) / sum(count.off[carryover.pos])
        }
    }

    ## 9. loess HTE by time
    D.missing <- D
    D.missing[which(D == 0)] <- NA
    eff.calendar <- apply(eff * D.missing, 1, mean, na.rm = TRUE)
    N.calendar <- apply(!is.na(eff * D.missing), 1, sum)
    T.calendar <- c(1:TT)
    if (sum(!is.na(eff.calendar)) > 1) {
        # loess fit
        if (!is.null(calendar.enp.seq)) {
            if (length(calendar.enp.seq) == 1 & is.na(calendar.enp.seq)) {
                calendar.enp.seq <- NULL
            }
        }
        if (is.null(calendar.enp.seq)) {
            loess.fit <- suppressWarnings(try(loess(eff.calendar ~ T.calendar, weights = N.calendar), silent = TRUE))
        } else {
            loess.fit <- suppressWarnings(try(loess(eff.calendar ~ T.calendar, weights = N.calendar, enp.target = calendar.enp.seq), silent = TRUE))
        }
        if ("try-error" %in% class(loess.fit)) {
            eff.calendar.fit <- eff.calendar
            calendar.enp <- NULL
        } else {
            eff.calendar.fit <- eff.calendar
            eff.calendar.fit[which(!is.na(eff.calendar))] <- loess.fit$fit
            calendar.enp <- loess.fit$enp
        }
    } else {
        eff.calendar.fit <- eff.calendar
        calendar.enp <- NULL
    }


    ## 8. cohort effects
    if (!is.null(group)) {
        cohort <- cbind(c(group), c(D), c(eff.v))
        rm.pos <- unique(c(rm.pos1, which(cohort[, 2] == 0)))
        cohort <- cohort[-rm.pos, ]

        g.level <- sort(unique(cohort[, 1]))
        raw.group.att <- as.numeric(tapply(cohort[, 3], cohort[, 1], mean))

        group.att <- rep(NA, length(group.level))
        group.att[which(group.level %in% g.level)] <- raw.group.att

        # by-group dynamic effects
        group.level.name <- names(group.level)

        group.output <- list()
        for (i in c(1:length(group.level))) {
            sub.group <- group.level[i]
            sub.group.name <- group.level.name[i]

            ## by-group dynamic effects
            t.on.sub <- c(T.on[which(group == sub.group)])
            eff.v.sub <- c(eff[which(group == sub.group)]) ## a vector
            rm.pos1.sub <- which(is.na(eff.v.sub))
            rm.pos2.sub <- which(is.na(t.on.sub))
            eff.v.use1.sub <- eff.v.sub
            t.on.use.sub <- t.on.sub
            if (NA %in% eff.v.sub | NA %in% t.on.sub) {
                eff.v.use1.sub <- eff.v.sub[-c(rm.pos1.sub, rm.pos2.sub)]
                t.on.use.sub <- t.on.sub[-c(rm.pos1.sub, rm.pos2.sub)]
            }
            if (length(t.on.use.sub) > 0) {
                time.on.sub <- sort(unique(t.on.use.sub))
                att.on.sub <- as.numeric(tapply(
                    eff.v.use1.sub,
                    t.on.use.sub,
                    mean
                )) ## NA already removed
                count.on.sub <- as.numeric(table(t.on.use.sub))
            } else {
                time.on.sub <- att.on.sub <- count.on.sub <- NULL
            }

            if (!is.null(time.on.seq.group)) {
                count.on.med.sub <- att.on.med.sub <- rep(NA, length(time.on.seq.group[[sub.group.name]]))
                time.on.seq.sub <- time.on.seq.group[[sub.group.name]]
                att.on.med.sub[which(time.on.seq.sub %in% time.on.sub)] <- att.on.sub
                count.on.med.sub[which(time.on.seq.sub %in% time.on.sub)] <- count.on.sub
                att.on.sub <- att.on.med.sub
                count.on.sub <- count.on.med.sub
                time.on.sub <- time.on.seq.sub
            }
            if (length(att.on.sub) == 0) {
                att.on.sub <- NULL
            }
            if (length(time.on.sub) == 0) {
                time.on.sub <- NULL
            }
            if (length(count.on.sub) == 0) {
                count.on.sub <- NULL
            }
            suboutput <- list(
                att.on = att.on.sub,
                time.on = time.on.sub,
                count.on = count.on.sub
            )

            ## placebo effect, if placeboTest == 1
            if (!is.null(placebo.period) && placeboTest == 1) {
                if (length(placebo.period) == 1) {
                    placebo.pos.sub <- which(time.on.sub == placebo.period)
                    if (length(placebo.pos.sub) > 0) {
                        att.placebo.sub <- att.on.sub[placebo.pos.sub]
                    } else {
                        att.placebo.sub <- NULL
                    }
                } else {
                    placebo.pos.sub <- which(time.on.sub >= placebo.period[1] & time.on.sub <= placebo.period[2])
                    if (length(placebo.pos.sub) > 0) {
                        att.placebo.sub <- sum(att.on.sub[placebo.pos.sub] * count.on.sub[placebo.pos.sub]) / sum(count.on.sub[placebo.pos.sub])
                    } else {
                        att.placebo.sub <- NULL
                    }
                }
                if (length(att.placebo.sub) == 0) {
                    att.placebo.sub <- NULL
                }
                suboutput <- c(suboutput, list(att.placebo = att.placebo.sub))
            }

            ## T.off
            if (hasRevs == 1) {
                t.off.sub <- c(T.off[which(group == sub.group)])
                rm.pos3.sub <- which(is.na(t.off.sub))
                eff.v.use2.sub <- eff.v.sub
                t.off.use.sub <- t.off.sub
                if (NA %in% eff.v.sub | NA %in% t.off.sub) {
                    eff.v.use2.sub <- eff.v.sub[-c(rm.pos1.sub, rm.pos3.sub)]
                    t.off.use.sub <- t.off.sub[-c(rm.pos1.sub, rm.pos3.sub)]
                }
                if (length(t.off.use.sub) > 0) {
                    time.off.sub <- sort(unique(t.off.use.sub))
                    att.off.sub <- as.numeric(tapply(eff.v.use2.sub, t.off.use.sub, mean)) ## NA already removed
                    count.off.sub <- as.numeric(table(t.off.use.sub))
                } else {
                    time.off.sub <- att.off.sub <- count.off.sub <- NULL
                }

                if (!is.null(time.off.seq.group)) {
                    count.off.med.sub <- att.off.med.sub <- rep(NA, length(time.off.seq.group[[sub.group.name]]))
                    time.off.seq.sub <- time.off.seq.group[[sub.group.name]]
                    att.off.med.sub[which(time.off.seq.sub %in% time.off.sub)] <- att.off.sub
                    count.off.med.sub[which(time.off.seq.sub %in% time.off.sub)] <- count.off.sub
                    att.off.sub <- att.off.med.sub
                    count.off.sub <- count.off.med.sub
                    time.off.sub <- time.off.seq.sub
                }
                if (length(att.off.sub) == 0) {
                    att.off.sub <- NULL
                }
                if (length(time.off.sub) == 0) {
                    time.off.sub <- NULL
                }
                if (length(count.off.sub) == 0) {
                    count.off.sub <- NULL
                }
                suboutput <- c(suboutput, list(
                    att.off = att.off.sub,
                    count.off = count.off.sub,
                    time.off = time.off.sub
                ))

                if (!is.null(carryover.period) && carryoverTest == 1) {
                    if (length(carryover.period) == 1) {
                        carryover.pos.sub <- which(time.off.sub == carryover.period)
                        if (length(carryover.pos.sub) > 0) {
                            att.carryover.sub <- att.off.sub[carryover.pos.sub]
                        } else {
                            att.carryover.sub <- NULL
                        }
                    } else {
                        carryover.pos.sub <- which(time.off.sub >= carryover.period[1] & time.off.sub <= carryover.period[2])
                        if (length(carryover.pos.sub) > 0) {
                            att.carryover.sub <- sum(att.off.sub[carryover.pos.sub] * count.off.sub[carryover.pos.sub]) / sum(count.off.sub[carryover.pos.sub])
                        } else {
                            att.carryover.sub <- NULL
                        }
                    }
                    if (length(att.carryover.sub) == 0) {
                        att.carryover.sub <- NULL
                    }
                    suboutput <- c(suboutput, list(att.carryover = att.carryover.sub))
                }
            }
            group.output[[sub.group.name]] <- suboutput
        }
    }

    ## Preserve method regardless of r.cv: the bootstrap dispatcher in fect_boot
    ## has branches for gsynth/ife/mc/cfe but not for "fe", so relabelling to
    ## "fe" when r.cv == 0 caused stop("Unsupported bootstrap method: fe") whenever
    ## CV selected the two-way FE model. The gsynth/cfe branches already handle
    ## r = 0 correctly on the estimation side (inside fect_nevertreated itself).
    if (method == "cfe") {
        method <- "cfe"
    } else {
        method <- "gsynth"
    }

    ## -------------------------------##
    ##            Storage            ##
    ## -------------------------------##
    out <- list(
        ## main results
        method = method,
        Y.ct = Y.ct,
        Y.tr.cnt = Y.ct.tr,
        Y.ct.cnt = Y.ct.co,
        Y.ct.full = Y.ct.full,
        Y = Y,
        D = D,
        tr = tr,
        co = co,
        eff = eff,
        eff.tr = eff[, tr],
        I = I,
        II = II,
        att.avg = att.avg,
        att.avg.unit = att.avg.unit,
        ## supporting
        force = force,
        T = TT,
        N = N,
        Ntr = Ntr,
        Nco = Nco,
        p = p,
        r.cv = r.cv,
        IC = IC,
        beta = beta,
        est = est.co.best,
        mu = est.co.best$mu,
        niter = est.co.best$niter,
        validX = validX,
        validF = validF,
        time = time.on,
        att = att.on,
        count = count.on,
        eff.calendar = eff.calendar,
        N.calendar = N.calendar,
        eff.calendar.fit = eff.calendar.fit,
        eff.pre = eff.pre,
        eff.pre.equiv = eff.pre.equiv,
        pre.sd = pre.sd,
        gamma = est.co.best$gamma,
        kappa = est.co.best$kappa
    )


    out <- c(out, list(
        PC = PC,
        sigma2 = sigma2,
        sigma2.fect = est.co.fect$sigma2,
        res = est.co.best$residuals,
        res.full = res.full,
        rmse = rmse
    ))


    if (hasRevs == 1) {
        out <- c(out, list(
            time.off = time.off,
            att.off = att.off,
            count.off = count.off,
            eff.off = eff.off,
            eff.off.equiv = eff.off.equiv,
            off.sd = off.sd
        ))
    }
    ## Include CV.out in the return list when cross-validation was performed
    if (exists("CV.out", inherits = FALSE)) {
        out <- c(out, list(CV.out = CV.out))
    }

    if (r.cv > 0) {
        lambda.co <- as.matrix(est.co.best$lambda)
        rownames(lambda.co) <- co
        rownames(lambda.tr) <- tr
        out <- c(out, list(
            factor = as.matrix(est.co.best$factor),
            lambda.co = as.matrix(lambda.co),
            lambda.tr = as.matrix(lambda.tr)
        ))
        if (boot == 0) {
            if (exists("use_bounded", inherits = FALSE) && isTRUE(use_bounded)) {
                out <- c(out, list(wgt.implied = wgt.implied))
            } else if (exists("inv.tr", inherits = FALSE) &&
                       !inherits(inv.tr, "try-error")) {
                out <- c(out, list(wgt.implied = wgt.implied))
            }
        }
        if (exists("use_bounded", inherits = FALSE) && isTRUE(use_bounded)) {
            out <- c(out, list(
                loading.bound       = "simplex",
                gamma.loading       = gamma_use,
                loading.proj.resid  = loading.proj.resid
            ))
        } else {
            out <- c(out, list(loading.bound = loading.bound))
        }
    }

    if (force == 1) {
        out <- c(out, list(
            alpha.tr = alpha.tr,
            alpha.co = est.co.best$alpha
        ))
    } else if (force == 2) {
        out <- c(out, list(xi = est.co.best$xi))
    } else if (force == 3) {
        out <- c(out, list(
            alpha.tr = alpha.tr,
            alpha.co = est.co.best$alpha,
            xi = est.co.best$xi
        ))
    }

    if (!is.null(placebo.period) && placeboTest == 1) {
        out <- c(out, list(att.placebo = att.placebo))
    }

    if (!is.null(W)) {
        out <- c(out, list(
            W = W,
            att.avg.W = att.avg.W,
            att.on.sum.W = att.on.sum.W,
            att.on.W = att.on.W,
            count.on.W = count.on.W,
            time.on.W = time.on.W,
            W.on.sum = W.on.sum
        ))
        if (hasRevs == 1) {
            out <- c(out, list(
                att.off.sum.W = att.off.sum.W,
                att.off.W = att.off.W,
                count.off.W = count.off.W,
                time.off.W = time.off.W,
                W.off.sum = W.off.sum
            ))
        }
        if (!is.null(placebo.period) && placeboTest == 1) {
            out <- c(out, list(att.placebo.W = att.placebo.W))
        }
    }

    if (!is.null(balance.period)) {
        out <- c(out, list(balance.att = balance.att, balance.time = balance.time, balance.count = balance.count, balance.avg.att = att.avg.balance))
        if (!is.null(placebo.period) && placeboTest == 1) {
            out <- c(out, list(balance.att.placebo = balance.att.placebo))
        }
    }

    if (!is.null(carryover.period) && carryoverTest == 1) {
        out <- c(out, list(att.carryover = att.carryover))
    }

    if (!is.null(group)) {
        out <- c(out, list(
            group.att = group.att,
            group.output = group.output
        ))
    }
    return(out)
}


## ====================================================================
## Helper functions for CFE nevertreated path
## ====================================================================

## Reconstruct gamma fitted values for treated units using co-estimated gamma
## gamma_coef: coefficient matrix from C++ Gamma(), n_groups x n_z_per_gamma
## X.Z.tr: TT x Ntr x n_z array
## gamma_groups_mat: TT x Ntr x 1 array (or TT x Ntr matrix via drop)
## z_cols: integer vector of which Z columns belong to this gamma
## Returns: TT x Ntr matrix
.reconstruct_gamma_fit_tr <- function(gamma_coef, X.Z.tr, gamma_groups_arr,
                                       z_cols, TT, Ntr) {
    ## gamma_groups_arr: TT x Ntr x 1 array (from X.gamma.tr[,,k,drop=FALSE])
    ## extract the TT x Ntr matrix
    if (length(dim(gamma_groups_arr)) == 3) {
        gamma_groups_mat <- gamma_groups_arr[, , 1, drop = FALSE]
        dim(gamma_groups_mat) <- dim(gamma_groups_mat)[1:2]
    } else {
        gamma_groups_mat <- gamma_groups_arr
    }
    if (is.null(dim(gamma_groups_mat))) {
        gamma_groups_mat <- matrix(gamma_groups_mat, nrow = TT, ncol = Ntr)
    }
    gamma_groups <- gamma_groups_mat[, 1]  ## all units share same time grouping
    unique_groups <- sort(unique(gamma_groups))

    fit <- matrix(0, TT, Ntr)

    ## Extract Z values for treated (time-invariant, take t=1)
    Z.tr.sub <- matrix(0, Ntr, length(z_cols))
    for (j in seq_along(z_cols)) {
        Z.tr.sub[, j] <- X.Z.tr[1, , z_cols[j]]
    }

    for (g_idx in seq_along(unique_groups)) {
        t_in_group <- which(gamma_groups == unique_groups[g_idx])
        if (g_idx <= nrow(gamma_coef)) {
            coef_g <- gamma_coef[g_idx, , drop = FALSE]  ## 1 x n_z
            contribution <- Z.tr.sub %*% t(coef_g)  ## Ntr x 1
            for (tt in t_in_group) {
                fit[tt, ] <- c(contribution)
            }
        }
    }
    return(fit)
}

## Reconstruct kappa fitted values
## kappa_coef: coefficient matrix from C++ Kappa()
## X.Q: TT x N x n_q array
## kappa_groups_mat: from X.kappa[,,k] -- unit-level grouping
## q_cols: integer vector
## Returns: TT x N matrix
.reconstruct_kappa_fit <- function(kappa_coef, X.Q, kappa_groups_mat, q_cols, TT, N) {
    ## kappa_coef: q x N matrix from C++ Kappa().t() — each column is a unit's
    ## kappa coefficients (q basis dimensions); units in the same group share
    ## identical columns.
    ## kappa_groups_mat may be a TT x N x 1 array (from X.kappa[,,k,drop=FALSE])
    ## extract the TT x N matrix
    if (length(dim(kappa_groups_mat)) == 3) {
        kappa_groups_mat <- kappa_groups_mat[, , 1, drop = FALSE]
        dim(kappa_groups_mat) <- dim(kappa_groups_mat)[1:2]
    }
    if (is.null(dim(kappa_groups_mat))) {
        kappa_groups_mat <- matrix(kappa_groups_mat, nrow = TT, ncol = N)
    }
    kappa_groups <- kappa_groups_mat[1, ]  ## N-length, unit grouping
    unique_groups <- sort(unique(kappa_groups))

    ## Ensure kappa_coef is a matrix (handles q=1 case)
    kappa_coef <- as.matrix(kappa_coef)

    fit <- matrix(0, TT, N)
    for (g_idx in seq_along(unique_groups)) {
        g <- unique_groups[g_idx]
        units_in_group <- which(kappa_groups == g)
        ## Extract kappa coefficients for this group from any unit in the group
        ## kappa_coef is q x N; column = unit's coefficients
        rep_unit <- units_in_group[1]
        if (rep_unit <= ncol(kappa_coef)) {
            kappa_g <- kappa_coef[, rep_unit, drop = FALSE]  ## q x 1
            ## Q basis values (time-varying)
            Q.basis <- matrix(0, TT, length(q_cols))
            for (j in seq_along(q_cols)) {
                Q.basis[, j] <- X.Q[, units_in_group[1], q_cols[j]]
            }
            contribution <- Q.basis %*% kappa_g  ## TT x 1
            for (ii in units_in_group) {
                fit[, ii] <- c(contribution)
            }
        }
    }
    return(fit)
}

## Extract Type-B FE from co estimation and apply to treated
## Returns: TT x Ntr matrix of Type-B FE fitted values for treated
.extract_and_apply_typeB_fe <- function(est.co, X.co, X.extra.FE.co.B,
                                         X.extra.FE.tr, typeB_idx,
                                         X.Z.co, X.gamma.co, X.Q.co, X.kappa.co,
                                         Zgamma.id, kappaQ.id,
                                         TT, Nco, Ntr, p, r, force) {
    ## Reconstruct non-FE fit for co
    beta <- est.co$beta
    if (!is.null(beta) && any(is.nan(beta))) beta[is.nan(beta)] <- 0

    fit_no_fe.co <- matrix(est.co$mu, TT, Nco)
    if (force %in% c(1, 3))
        fit_no_fe.co <- fit_no_fe.co + matrix(c(est.co$alpha), TT, Nco, byrow = TRUE)
    if (force %in% c(2, 3))
        fit_no_fe.co <- fit_no_fe.co + matrix(c(est.co$xi), TT, Nco, byrow = FALSE)
    if (p > 0 && !is.null(beta)) {
        for (j in 1:p) fit_no_fe.co <- fit_no_fe.co + X.co[, , j] * beta[j]
    }
    ## Add gamma for co
    if (!is.null(est.co$gamma) && length(est.co$gamma) > 0) {
        for (k_g in seq_along(est.co$gamma)) {
            gamma.fit.co.k <- .reconstruct_gamma_fit_tr(
                est.co$gamma[[k_g]], X.Z.co,
                X.gamma.co[, , k_g, drop = FALSE],
                Zgamma.id[[k_g]], TT, Nco)
            fit_no_fe.co <- fit_no_fe.co + gamma.fit.co.k
        }
    }
    ## Add kappa for co
    if (!is.null(est.co$kappa) && length(est.co$kappa) > 0) {
        for (k_k in seq_along(est.co$kappa)) {
            kappa.fit.co.k <- .reconstruct_kappa_fit(
                est.co$kappa[[k_k]], X.Q.co,
                X.kappa.co[, , k_k, drop = FALSE],
                kappaQ.id[[k_k]], TT, Nco)
            fit_no_fe.co <- fit_no_fe.co + kappa.fit.co.k
        }
    }
    ## Add factors for co
    if (r > 0 && !is.null(est.co$factor)) {
        fit_no_fe.co <- fit_no_fe.co + est.co$factor %*% t(est.co$lambda)
    }

    fe_fit.co <- est.co$fit - fit_no_fe.co  ## TT x Nco: extra FE contribution

    ## For each Type-B FE dimension, compute group means and apply to treated
    typeB.fit.tr <- matrix(0, TT, Ntr)
    for (k_b_idx in seq_along(typeB_idx)) {
        k <- typeB_idx[k_b_idx]
        labels.co <- X.extra.FE.co.B[1, , k_b_idx]
        labels.tr <- X.extra.FE.tr[1, , k]
        levels.all <- sort(unique(labels.co))

        for (g in levels.all) {
            co_mask <- which(labels.co == g)
            tr_mask <- which(labels.tr == g)
            if (length(tr_mask) > 0 && length(co_mask) > 0) {
                fe_val <- mean(fe_fit.co[, co_mask])
                typeB.fit.tr[, tr_mask] <- typeB.fit.tr[, tr_mask] + fe_val
            }
        }
    }
    return(typeB.fit.tr)
}
