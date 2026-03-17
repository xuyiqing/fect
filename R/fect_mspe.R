fect_mspe <- function(
    out.fect,
    hide_mask   = NULL,
    hide_n      = 20,
    seed        = NULL,
    n_rep       = 1,
    ## ----- masking strategy ----- ##
    pre.trend   = NULL,           # DEPRECATED — use mask.method instead
    mask.method = NULL,           # "random" | "pre.trend" | "cv.sample" | "user"
    pre.trend.n = 2,              # for mask.method = "pre.trend" or legacy pre.trend = TRUE
    cv.treat    = TRUE,           # for mask.method = "cv.sample"
    cv.nobs     = 3,              # for mask.method = "cv.sample"
    cv.donut    = 1,              # for mask.method = "cv.sample"
    cv.prop     = 0.1,            # for mask.method = "cv.sample"
    min.T0      = 5,              # for mask.method = "cv.sample"
    k           = 5,              # for mask.method = "cv.sample" (number of folds)
    ## ----- scoring ----- ##
    criterion   = "mspe",         # "mspe","wmspe","gmspe","wgmspe","mad","moment","gmoment"
    W           = NULL,           # TT x N observation weight matrix, or NULL
    norm.para   = NULL,           # normalization vector, or NULL
    proportion  = 0,              # proportion cutoff for count.T.cv (same as fect_cv)
    ## ----- ground truth ----- ##
    actual      = NULL,           # TT x N ground-truth matrix (default: observed Y)
    control.only = TRUE           # mask only D==0 cells? FALSE = also mask treated
) {
    if (!is.null(seed)) set.seed(seed)
    hide_n_given <- !missing(hide_n)
    caller_env <- parent.frame()

    ## ---- validate criterion ---- ##
    valid_criteria <- c("mspe", "wmspe", "gmspe", "wgmspe", "mad", "moment", "gmoment")
    criterion <- match.arg(criterion, valid_criteria)

    ## ---- backward compat for pre.trend ---- ##
    if (!is.null(pre.trend) && is.null(mask.method)) {
        mask.method <- if (isTRUE(pre.trend)) "pre.trend" else "random"
        message("Note: 'pre.trend' is deprecated; use mask.method='",
                mask.method, "' instead.")
    } else if (!is.null(pre.trend) && !is.null(mask.method)) {
        message("Note: both 'pre.trend' and 'mask.method' given; using mask.method='",
                mask.method, "'.")
    }
    if (is.null(mask.method)) {
        mask.method <- if (!is.null(hide_mask)) "user" else "random"
    }
    mask.method <- match.arg(mask.method,
                             c("random", "pre.trend", "cv.sample", "user"))

    ## ---- helper functions (unchanged) ---- ##
    .build_rerun_args <- function(out_obj, formula_obj, data_obj, index_obj, caller_env) {
        formula_env <- environment(out_obj$call$formula)
        if (is.null(formula_env)) formula_env <- caller_env
        call_args <- as.list(out_obj$call)[-1]
        rerun_args <- list(
            formula = formula_obj,
            data = data_obj,
            index = index_obj
        )
        arg_names <- setdiff(names(call_args), c("", "formula", "data", "index"))
        for (nm in arg_names) {
            rerun_args[[nm]] <- eval(call_args[[nm]], envir = formula_env, enclos = caller_env)
        }
        rerun_args
    }
    .is_fect_output <- function(obj) {
        is.list(obj) && !is.null(obj$Y.ct.full)
    }

    .get_last_matrix <- function(out_obj, nm, TT, N) {
        idx <- which(names(out_obj) == nm)
        if (length(idx) == 0) return(NULL)
        for (kk in rev(idx)) {
            obj <- out_obj[[kk]]
            if (is.matrix(obj) && all(dim(obj) == c(TT, N))) return(obj)
        }
        NULL
    }

    .as_mask <- function(mask, TT, N) {
        if (is.null(mask)) return(NULL)
        if (is.matrix(mask)) {
            if (all(dim(mask) == c(TT, N))) return(matrix(as.logical(mask), nrow = TT, ncol = N))
            if (all(dim(mask) == c(N, TT))) return(matrix(as.logical(t(mask)), nrow = TT, ncol = N))
            stop("`hide_mask` must be TT x N or N x TT.")
        }
        if (is.vector(mask) && length(mask) == TT * N) {
            return(matrix(as.logical(mask), nrow = TT, ncol = N))
        }
        stop("`hide_mask` must be NULL, a matrix, or a vector of length TT*N.")
    }

    .pre_trend_candidates <- function(d_mat, y_mat, n_last) {
        cand <- integer(0)
        if (!isTRUE(n_last > 0)) return(cand)
        for (j in seq_len(ncol(d_mat))) {
            treated_idx <- which(d_mat[, j] == 1)
            if (length(treated_idx) == 0) next
            first_treat <- min(treated_idx)
            pre_idx <- which(d_mat[, j] == 0 & !is.na(y_mat[, j]) & seq_len(nrow(d_mat)) < first_treat)
            if (length(pre_idx) == 0) next
            keep_idx <- tail(pre_idx, min(n_last, length(pre_idx)))
            cand <- c(cand, keep_idx + (j - 1) * nrow(d_mat))
        }
        unique(cand)
    }

    ## ---- validate out.fect ---- ##
    out_list <- if (.is_fect_output(out.fect)) list(out.fect) else out.fect
    if (!is.list(out_list) || length(out_list) == 0 ||
        !all(vapply(out_list, .is_fect_output, logical(1)))) {
        stop("`out.fect` must be a fect output or list of fect outputs.")
    }
    if (is.null(names(out_list)) || any(names(out_list) == "")) {
        names(out_list) <- paste0("out.fect_", seq_along(out_list))
    }

    ref <- out_list[[1]]
    TT <- nrow(ref$Y.ct.full)
    N <- ncol(ref$Y.ct.full)
    y_ref <- .get_last_matrix(ref, "Y", TT, N)
    if (is.null(y_ref) && !is.null(ref$Y.dat)) y_ref <- ref$Y.dat
    d_ref <- .get_last_matrix(ref, "D", TT, N)
    if (is.null(d_ref) && !is.null(ref$D.dat)) d_ref <- ref$D.dat
    if (is.null(y_ref) || is.null(d_ref)) stop("`out.fect` must contain matrix Y and D.")

    ## validate pre.trend.n
    if (!is.numeric(pre.trend.n) || length(pre.trend.n) != 1 || is.na(pre.trend.n) ||
        pre.trend.n < 1) {
        stop("`pre.trend.n` must be a positive number.")
    }
    pre.trend.n <- as.integer(pre.trend.n)

    ## validate mask.method = "user" requires hide_mask
    if (mask.method == "user" && is.null(hide_mask)) {
        stop("mask.method='user' requires a non-NULL hide_mask.")
    }

    ## validate W dimensions
    if (!is.null(W)) {
        if (!is.matrix(W) || nrow(W) != TT || ncol(W) != N) {
            stop("`W` must be a TT x N matrix.")
        }
    }

    ## validate actual dimensions
    if (!is.null(actual)) {
        if (!is.matrix(actual) || nrow(actual) != TT || ncol(actual) != N) {
            stop("`actual` must be a TT x N matrix.")
        }
    }

    ## ---- count.T.cv construction ---- ##
    ## Build period-level weights when criterion needs them
    count.T.cv <- NULL
    needs_period_weights <- criterion %in% c("wmspe", "wgmspe", "moment", "gmoment") ||
        mask.method == "cv.sample"
    if (needs_period_weights && !is.null(ref$T.on)) {
        T_on <- ref$T.on
        count.T.cv <- table(T_on)
        count.T.cv <- count.T.cv[which(as.numeric(names(count.T.cv)) <= 0)]
        cv.prop.cut <- max(count.T.cv) * proportion
        cv.drop.index <- which(count.T.cv <= cv.prop.cut)
        count.T.cv <- count.T.cv / mean(count.T.cv)
        nm <- names(count.T.cv)
        count.T.cv <- c(count.T.cv, median(count.T.cv))
        names(count.T.cv) <- c(nm, "Control")
        count.T.cv[cv.drop.index] <- 0
    }

    ## ---- cv.sample masking setup ---- ##
    rmCV <- estCV <- NULL
    if (mask.method == "cv.sample") {
        I_mat <- if (!is.null(ref$I.dat)) ref$I.dat else {
            tmp <- .get_last_matrix(ref, "I", TT, N)
            if (is.null(tmp)) matrix(1L, TT, N)
            else tmp
        }
        D_mat <- d_ref
        II_mat <- I_mat
        II_mat[D_mat == 1] <- 0
        if (!is.null(ref$factors.from) && ref$factors.from == "nevertreated") {
            D_cum <- apply(D_mat, 2, function(vec) cummax(vec))
            ever_treated <- which(colSums(D_cum) > 0)
            II_mat[, ever_treated] <- 0
        }
        rm_count <- floor(sum(II_mat) * cv.prop)

        rmCV <- list()
        estCV <- list()
        for (ii in 1:k) {
            cv.n <- 0
            repeat {
                cv.n <- cv.n + 1
                get.cv <- cv.sample(II_mat, D_mat,
                    count = rm_count,
                    cv.count = cv.nobs,
                    cv.treat = cv.treat,
                    cv.donut = cv.donut
                )
                cv.id <- get.cv$cv.id
                II.cv <- II_mat
                II.cv[cv.id] <- 0

                con1 <- sum(apply(II.cv, 1, sum) >= 1) == TT
                con2 <- sum(apply(II.cv, 2, sum) >= min.T0) == N

                if (con1 && con2) break
                if (cv.n >= 200) {
                    keep.1 <- which(apply(II.cv, 1, sum) < 1)
                    keep.2 <- which(apply(II.cv, 2, sum) < min.T0)
                    II.cv[keep.1, ] <- II_mat[keep.1, ]
                    II.cv[, keep.2] <- II_mat[, keep.2]
                    II.cv.valid <- II_mat
                    II.cv.valid[cv.id] <- -1
                    II.cv.valid[keep.1, ] <- II_mat[keep.1, ]
                    II.cv.valid[, keep.2] <- II_mat[, keep.2]
                    cv.id <- which(II.cv.valid != II_mat)
                    break
                }
            }
            rmCV[[ii]] <- cv.id
            if (cv.n < 200) {
                estCV[[ii]] <- get.cv$est.id
            } else {
                cv.id.old <- get.cv$cv.id
                cv.diff <- setdiff(cv.id.old, cv.id)
                estCV[[ii]] <- setdiff(get.cv$est.id, cv.diff)
            }
        }
    }

    ## ---- records data.frame ---- ##
    mask_base <- .as_mask(hide_mask, TT, N)
    records <- data.frame(
        Rep = integer(0),
        Model = character(0),
        Hidden_N = integer(0),
        MSPE = numeric(0), WMSPE = numeric(0), GMSPE = numeric(0),
        WGMSPE = numeric(0), MAD = numeric(0),
        Moment = numeric(0), GMoment = numeric(0),
        RMSE = numeric(0), Bias = numeric(0),
        stringsAsFactors = FALSE
    )
    fits_last <- NULL
    mask_last <- NULL
    scores_last <- NULL

    ## ---- main repetition loop ---- ##
    for (rep_id in seq_len(n_rep)) {

        ## ---- build mask for this rep ---- ##
        if (mask.method == "cv.sample") {
            ## cv.sample masking is per-fold inside the model loop below
            mask_mat <- NULL
        } else if (mask.method == "user") {
            mask_mat <- mask_base
            if (!is.null(mask_mat) && isTRUE(control.only)) {
                ## restrict to D==0 cells
                mask_mat[d_ref == 1] <- FALSE
            }
            if (!is.null(mask_mat) && hide_n_given) {
                cand <- which(mask_mat)
                pick <- sample(cand, min(hide_n, length(cand)))
                mask_mat <- matrix(FALSE, nrow = TT, ncol = N)
                mask_mat[pick] <- TRUE
            }
        } else if (mask.method == "pre.trend") {
            cand <- .pre_trend_candidates(d_ref, y_ref, pre.trend.n)
            if (length(cand) == 0) {
                ## fall back to random among D==0
                if (isTRUE(control.only)) {
                    cand <- which(d_ref == 0 & !is.na(y_ref))
                } else {
                    cand <- which(!is.na(y_ref))
                }
            }
            pick <- sample(cand, min(hide_n, length(cand)))
            mask_mat <- matrix(FALSE, nrow = TT, ncol = N)
            mask_mat[pick] <- TRUE
        } else {
            ## mask.method == "random"
            if (isTRUE(control.only)) {
                cand <- which(d_ref == 0 & !is.na(y_ref))
            } else {
                cand <- which(!is.na(y_ref))
            }
            pick <- sample(cand, min(hide_n, length(cand)))
            mask_mat <- matrix(FALSE, nrow = TT, ncol = N)
            mask_mat[pick] <- TRUE
        }

        fits <- vector("list", length(out_list))
        names(fits) <- names(out_list)

        for (i in seq_along(out_list)) {
            out_i <- out_list[[i]]
            y_true_i <- .get_last_matrix(out_i, "Y", TT, N)
            if (is.null(y_true_i) && !is.null(out_i$Y.dat)) y_true_i <- out_i$Y.dat
            d_mat_i <- .get_last_matrix(out_i, "D", TT, N)
            if (is.null(d_mat_i) && !is.null(out_i$D.dat)) d_mat_i <- out_i$D.dat
            if (is.null(y_true_i) || is.null(d_mat_i)) {
                stop("Each out.fect must provide Y and D matrices.")
            }

            ## determine actual values for scoring
            actual_mat <- if (!is.null(actual)) actual else y_true_i

            formula_env_i <- environment(out_i$call$formula)
            if (is.null(formula_env_i)) formula_env_i <- caller_env
            data_i <- eval(out_i$call$data, envir = formula_env_i, enclos = caller_env)
            idx_i <- eval(out_i$call$index, envir = formula_env_i, enclos = caller_env)
            formula_obj_i <- tryCatch(
                eval(out_i$call$formula, envir = formula_env_i, enclos = caller_env),
                error = function(e) NULL
            )
            if (!inherits(formula_obj_i, "formula")) {
                formula_obj_i <- stats::reformulate(c(out_i$D, out_i$X), response = out_i$Y)
            }
            f_vars <- all.vars(formula_obj_i)
            y_col_i <- f_vars[1]

            rr_i <- match(data_i[[idx_i[2]]], out_i$rawtime)
            cc_i <- match(data_i[[idx_i[1]]], out_i$id)

            if (mask.method == "cv.sample") {
                ## ---- cv.sample: accumulate residuals across k folds ---- ##
                all_resid <- c()
                all_obs_w <- c()
                all_time_idx <- c()

                for (ii in 1:k) {
                    ## hide fold ii positions in long-form data
                    long_rm <- rep(FALSE, nrow(data_i))
                    valid_rm <- !is.na(rr_i) & !is.na(cc_i)
                    mat_idx_rm <- cbind(rr_i[valid_rm], cc_i[valid_rm])
                    in_rm <- rep(FALSE, nrow(data_i))
                    in_rm[valid_rm] <- (mat_idx_rm[, 1] + (mat_idx_rm[, 2] - 1) * TT) %in% rmCV[[ii]]
                    hide_rows_ii <- which(in_rm)

                    if (length(hide_rows_ii) == 0) next

                    data_hidden <- data_i
                    data_hidden[hide_rows_ii, y_col_i] <- NA_real_

                    rerun_args <- .build_rerun_args(
                        out_obj = out_i,
                        formula_obj = formula_obj_i,
                        data_obj = data_hidden,
                        index_obj = idx_i,
                        caller_env = caller_env
                    )
                    out_new <- do.call(fect, rerun_args)

                    ## collect residuals at estCV positions
                    ## est_pos are flat column-major indices into the
                    ## *original* TT x N matrix (II_mat).  out_new from
                    ## the rerun may have different dimensions (dropped
                    ## rows/cols when all Y are NA), so we map through
                    ## rawtime/id labels instead of using flat indices
                    ## directly on out_new$Y.ct.full.
                    est_pos <- estCV[[ii]]
                    est_row <- ((est_pos - 1) %% TT) + 1
                    est_col <- ((est_pos - 1) %/% TT) + 1

                    ## map original (row, col) to rawtime/id labels,
                    ## then match into the rerun's coordinate system
                    est_time_labels <- out_i$rawtime[est_row]
                    est_id_labels   <- out_i$id[est_col]
                    rr_new <- match(est_time_labels, out_new$rawtime)
                    cc_new <- match(est_id_labels,   out_new$id)
                    valid_map <- !is.na(rr_new) & !is.na(cc_new)

                    pred_vals <- rep(NA_real_, length(est_pos))
                    pred_vals[valid_map] <- out_new$Y.ct.full[cbind(rr_new[valid_map], cc_new[valid_map])]
                    ## actual_mat uses the original coordinate system (always TT x N)
                    actual_vals <- actual_mat[cbind(est_row, est_col)]
                    valid_score <- valid_map & !is.na(pred_vals) & !is.na(actual_vals)

                    resid_ii <- pred_vals[valid_score] - actual_vals[valid_score]
                    all_resid <- c(all_resid, resid_ii)

                    if (!is.null(W)) {
                        all_obs_w <- c(all_obs_w, W[cbind(est_row[valid_score], est_col[valid_score])])
                    }

                    ## time index for moment/gmoment
                    if (!is.null(out_i$T.on)) {
                        tidx <- as.character(out_i$T.on[cbind(est_row[valid_score], est_col[valid_score])])
                        tidx[which(is.na(tidx))] <- "Control"
                        all_time_idx <- c(all_time_idx, tidx)
                    }
                }

                if (length(all_resid) == 0) {
                    stop("No valid residuals collected from cv.sample folds.")
                }
                fits[[i]] <- NULL  ## cv.sample doesn't retain a single refit

                scores <- .score_residuals(
                    resid = all_resid,
                    obs_weights = if (length(all_obs_w) > 0) all_obs_w else NULL,
                    time_index = if (length(all_time_idx) > 0) all_time_idx else NULL,
                    count_weights = count.T.cv,
                    norm.para = norm.para
                )

            } else {
                ## ---- non-cv.sample: single refit with masked data ---- ##
                long_mask_i <- !is.na(rr_i) & !is.na(cc_i) &
                    as.logical(mask_mat[cbind(rr_i, cc_i)])
                hide_rows_i <- which(long_mask_i)
                if (length(hide_rows_i) == 0) {
                    stop("No valid hide positions mapped to source data rows.")
                }

                data_hidden <- data_i
                data_hidden[hide_rows_i, y_col_i] <- NA_real_

                rerun_args <- .build_rerun_args(
                    out_obj = out_i,
                    formula_obj = formula_obj_i,
                    data_obj = data_hidden,
                    index_obj = idx_i,
                    caller_env = caller_env
                )

                out_new <- do.call(fect, rerun_args)
                fits[[i]] <- out_new

                rr_new <- match(data_i[[idx_i[2]]][hide_rows_i], out_new$rawtime)
                cc_new <- match(data_i[[idx_i[1]]][hide_rows_i], out_new$id)
                valid_pred <- !is.na(rr_new) & !is.na(cc_new)
                if (!any(valid_pred)) {
                    stop("No hidden positions remain comparable after rerun; try smaller `hide_n` or `pre.trend.n`.")
                }
                pred <- out_new$Y.ct.full[cbind(rr_new[valid_pred], cc_new[valid_pred])]
                actual_vals <- actual_mat[cbind(
                    rr_i[hide_rows_i][valid_pred],
                    cc_i[hide_rows_i][valid_pred]
                )]
                resid_vec <- pred - actual_vals

                ## observation weights at evaluation positions
                obs_w_vec <- NULL
                if (!is.null(W)) {
                    obs_w_vec <- W[cbind(
                        rr_i[hide_rows_i][valid_pred],
                        cc_i[hide_rows_i][valid_pred]
                    )]
                }

                ## time index for moment/gmoment
                time_idx_vec <- NULL
                if (!is.null(out_i$T.on)) {
                    time_idx_vec <- as.character(out_i$T.on[cbind(
                        rr_i[hide_rows_i][valid_pred],
                        cc_i[hide_rows_i][valid_pred]
                    )])
                    time_idx_vec[which(is.na(time_idx_vec))] <- "Control"
                }

                scores <- .score_residuals(
                    resid = resid_vec,
                    obs_weights = obs_w_vec,
                    time_index = time_idx_vec,
                    count_weights = count.T.cv,
                    norm.para = norm.para
                )
            }

            scores_last <- scores

            records <- rbind(
                records,
                data.frame(
                    Rep = rep_id,
                    Model = names(out_list)[i],
                    Hidden_N = if (mask.method == "cv.sample") length(all_resid) else sum(valid_pred),
                    MSPE = scores["MSPE"],
                    WMSPE = scores["WMSPE"],
                    GMSPE = scores["GMSPE"],
                    WGMSPE = scores["WGMSPE"],
                    MAD = scores["MAD"],
                    Moment = scores["Moment"],
                    GMoment = scores["GMoment"],
                    RMSE = scores["RMSE"],
                    Bias = scores["Bias"],
                    stringsAsFactors = FALSE,
                    row.names = NULL
                )
            )
        }
        fits_last <- fits
        if (!is.null(mask_mat)) mask_last <- mask_mat
    }

    summary <- aggregate(
        cbind(Hidden_N, MSPE, WMSPE, GMSPE, WGMSPE, MAD, Moment, GMoment, RMSE, Bias) ~ Model,
        data = records, FUN = mean
    )
    list(
        summary = summary,
        records = records,
        hide_mask = mask_last,
        fits = fits_last,
        criterion = criterion,
        scores = scores_last
    )
}
