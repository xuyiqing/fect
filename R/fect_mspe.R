fect_mspe <- function(out.fect, hide_mask = NULL, hide_n = 20, seed = NULL, n_rep = 1) {
    if (!is.null(seed)) set.seed(seed)
    hide_n_given <- !missing(hide_n)
    .is_fect_output <- function(obj) {
        is.list(obj) && !is.null(obj$Y.ct.full)
    }

    .get_last_matrix <- function(out_obj, nm, TT, N) {
        idx <- which(names(out_obj) == nm)
        if (length(idx) == 0) return(NULL)
        for (k in rev(idx)) {
            obj <- out_obj[[k]]
            if (is.matrix(obj) && all(dim(obj) == c(TT, N))) return(obj)
        }
        NULL
    }

    .get_last_X_array <- function(out_obj, TT, N) {
        idx <- which(names(out_obj) == "X")
        if (length(idx) == 0) return(NULL)
        for (k in rev(idx)) {
            obj <- out_obj[[k]]
            d <- dim(obj)
            if (!is.null(d) && length(d) == 3 && d[1] == TT && d[2] == N) return(obj)
        }
        NULL
    }

    .as_mask <- function(mask, TT, N) {
        if (is.null(mask)) return(NULL)
        if (is.matrix(mask)) {
            if (all(dim(mask) == c(TT, N))) return(as.logical(mask))
            if (all(dim(mask) == c(N, TT))) return(as.logical(t(mask)))
            stop("`hide_mask` must be TT x N or N x TT.")
        }
        if (is.vector(mask) && length(mask) == TT * N) {
            return(matrix(as.logical(mask), nrow = TT, ncol = N))
        }
        stop("`hide_mask` must be NULL, a matrix, or a vector of length TT*N.")
    }

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
    d_ref <- .get_last_matrix(ref, "D", TT, N)
    if (is.null(y_ref) || is.null(d_ref)) stop("`out.fect` must contain matrix Y and D.")

    mask_base <- .as_mask(hide_mask, TT, N)
    records <- data.frame(
        Rep = integer(0),
        Model = character(0),
        Hidden_N = integer(0),
        RMSE = numeric(0),
        Bias = numeric(0),
        stringsAsFactors = FALSE
    )
    fits_last <- NULL
    mask_last <- NULL

    for (rep_id in seq_len(n_rep)) {
        mask_mat <- mask_base
        if (is.null(mask_mat)) {
            cand <- which(d_ref == 0 & !is.na(y_ref))
            pick <- sample(cand, min(hide_n, length(cand)))
            mask_mat <- matrix(FALSE, nrow = TT, ncol = N)
            mask_mat[pick] <- TRUE
        } else {
            cand <- which(mask_mat & (d_ref == 0))
            if (hide_n_given) {
                pick <- sample(cand, min(hide_n, length(cand)))
                mask_mat <- matrix(FALSE, nrow = TT, ncol = N)
                mask_mat[pick] <- TRUE
            } else {
                mask_mat <- matrix(FALSE, nrow = TT, ncol = N)
                mask_mat[cand] <- TRUE
            }
        }

        fits <- vector("list", length(out_list))
        names(fits) <- names(out_list)

        for (i in seq_along(out_list)) {
            out_i <- out_list[[i]]
            y_true <- .get_last_matrix(out_i, "Y", TT, N)
            d_mat <- .get_last_matrix(out_i, "D", TT, N)
            x_arr <- .get_last_X_array(out_i, TT, N)
            if (is.null(y_true) || is.null(d_mat)) stop("Each out.fect must provide Y and D matrices.")

            y_col <- ".__Y__"
            d_col <- ".__D__"
            x_cols <- if (!is.null(x_arr) && dim(x_arr)[3] > 0) {
                paste0(".__X__", seq_len(dim(x_arr)[3]))
            } else {
                character(0)
            }

            data_long <- data.frame(
                id = rep(out_i$id, each = TT),
                time = rep(out_i$rawtime, N),
                stringsAsFactors = FALSE
            )
            data_long[[y_col]] <- as.numeric(c(y_true))
            data_long[[d_col]] <- as.numeric(c(d_mat))
            if (length(x_cols) > 0) {
                if (is.null(x_arr) || length(dim(x_arr)) != 3 || dim(x_arr)[3] < length(x_cols)) {
                    stop("X array in out.fect is not compatible with formula covariates.")
                }
                for (j in seq_along(x_cols)) {
                    data_long[[x_cols[j]]] <- as.numeric(c(x_arr[, , j]))
                }
            }

            long_mask <- as.vector(mask_mat)
            data_hidden <- data_long
            data_hidden[[y_col]][long_mask] <- NA_real_

            formula_obj <- stats::reformulate(c(d_col, x_cols), response = y_col)
            rerun_args <- list(
                formula = formula_obj,
                data = data_hidden,
                index = c("id", "time")
            )
            if (!is.null(out_i$call$method)) rerun_args$method <- eval(out_i$call$method, parent.frame())
            if (!is.null(out_i$call$r)) rerun_args$r <- eval(out_i$call$r, parent.frame())

            out_new <- do.call(fect::fect, rerun_args)
            fits[[i]] <- out_new

            pred <- out_new$Y.ct.full[mask_mat]
            actual <- y_true[mask_mat]
            err <- pred - actual
            records <- rbind(
                records,
                data.frame(
                    Rep = rep_id,
                    Model = names(out_list)[i],
                    Hidden_N = length(err),
                    RMSE = sqrt(mean(err^2, na.rm = TRUE)),
                    Bias = mean(err, na.rm = TRUE),
                    stringsAsFactors = FALSE
                )
            )
        }
        fits_last <- fits
        mask_last <- mask_mat
    }

    summary <- aggregate(cbind(Hidden_N, RMSE, Bias) ~ Model, data = records, FUN = mean)
    list(summary = summary, records = records, hide_mask = mask_last, fits = fits_last)
}


fect_mspe_sim <- function(out.fect, hide_mask = NULL, hide_mask_y0 = NULL, hide_n = 20, seed = NULL, n_rep = 1) {
    if (!is.null(seed)) set.seed(seed)
    hide_n_given <- !missing(hide_n)
    .is_fect_output <- function(obj) {
        is.list(obj) && !is.null(obj$Y.ct.full) && !is.null(obj$call)
    }
    .get_last_matrix <- function(out_obj, nm, TT, N) {
        idx <- which(names(out_obj) == nm)
        if (length(idx) == 0) return(NULL)
        for (k in rev(idx)) {
            obj <- out_obj[[k]]
            if (is.matrix(obj) && all(dim(obj) == c(TT, N))) return(obj)
        }
        NULL
    }
    .as_mask <- function(mask, TT, N) {
        if (is.null(mask)) return(NULL)
        if (is.matrix(mask)) {
            if (all(dim(mask) == c(TT, N))) return(as.logical(mask))
            if (all(dim(mask) == c(N, TT))) return(as.logical(t(mask)))
            stop("`hide_mask` must be TT x N or N x TT.")
        }
        if (is.vector(mask) && length(mask) == TT * N) {
            return(matrix(as.logical(mask), nrow = TT, ncol = N))
        }
        stop("`hide_mask` must be NULL, matrix, or vector length TT*N.")
    }
    .as_matrix <- function(mask, TT, N) {
        if (is.null(mask)) return(NULL)
        if (is.matrix(mask)) {
            if (all(dim(mask) == c(TT, N))) return(mask)
            if (all(dim(mask) == c(N, TT))) return(t(mask))
            stop("`hide_mask` must be TT x N or N x TT.")
        }
        if (is.vector(mask) && length(mask) == TT * N) {
            return(matrix(mask, nrow = TT, ncol = N))
        }
        stop("`hide_mask` must be NULL, matrix, or vector length TT*N.")
    }

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
    formula_env <- environment(ref$call$formula)
    data_ref <- eval(ref$call$data, envir = formula_env, enclos = parent.frame())
    idx_ref <- eval(ref$call$index, envir = formula_env, enclos = parent.frame())
    if (!"Y0" %in% colnames(data_ref)) stop("Source data must contain Y0.")

    rr <- match(data_ref[[idx_ref[2]]], ref$rawtime)
    cc <- match(data_ref[[idx_ref[1]]], ref$id)
    valid_map <- !is.na(rr) & !is.na(cc)
    d_ref_mat <- .get_last_matrix(ref, "D", TT, N)
    if (is.null(d_ref_mat)) stop("Cannot find TT x N D matrix in out.fect.")
    d_long <- rep(NA_real_, nrow(data_ref))
    d_long[valid_map] <- as.numeric(d_ref_mat[cbind(rr[valid_map], cc[valid_map])])

    mask_mat <- .as_mask(hide_mask, TT, N)
    if (is.null(mask_mat)) {
        stop("For fect_mspe_sim, provide `hide_mask` as a 0/1 matrix.")
    }
    y0_template <- .as_matrix(hide_mask_y0, TT, N)
    if (is.null(y0_template)) {
        y0_template <- matrix(NA_real_, nrow = TT, ncol = N)
        y0_template[cbind(rr[valid_map], cc[valid_map])] <- as.numeric(data_ref[valid_map, "Y0"])
    }
    long_mask <- as.logical(mask_mat[cbind(rr, cc)])
    long_val <- y0_template[cbind(rr, cc)]
    cand_rows <- which(valid_map & long_mask & !is.na(long_val))
    if (length(cand_rows) == 0) {
        stop("No valid hide positions with non-missing Y0 values.")
    }
    hide_rows <- if (hide_n_given && hide_n < length(cand_rows)) {
        sample(cand_rows, hide_n)
    } else {
        cand_rows
    }
    mask_mat <- matrix(FALSE, nrow = TT, ncol = N)
    mask_mat[cbind(rr[hide_rows], cc[hide_rows])] <- TRUE
    actual_mat <- matrix(NA_real_, nrow = TT, ncol = N)
    actual_mat[cbind(rr[hide_rows], cc[hide_rows])] <- as.numeric(long_val[hide_rows])

    records <- data.frame(
        Rep = integer(0),
        Model = character(0),
        Hidden_N = integer(0),
        RMSE = numeric(0),
        Bias = numeric(0),
        stringsAsFactors = FALSE
    )
    fits_last <- NULL
    mask_last <- NULL

    for (rep_id in seq_len(n_rep)) {
        hide_rows_rep <- if (hide_n_given && hide_n < length(cand_rows)) {
            sample(cand_rows, hide_n)
        } else {
            cand_rows
        }
        mask_mat_rep <- matrix(FALSE, nrow = TT, ncol = N)
        mask_mat_rep[cbind(rr[hide_rows_rep], cc[hide_rows_rep])] <- TRUE
        actual_mat_rep <- matrix(NA_real_, nrow = TT, ncol = N)
        actual_mat_rep[cbind(rr[hide_rows_rep], cc[hide_rows_rep])] <- as.numeric(long_val[hide_rows_rep])

        fits <- vector("list", length(out_list))
        names(fits) <- names(out_list)

        for (i in seq_along(out_list)) {
            out_i <- out_list[[i]]
            formula_env_i <- environment(out_i$call$formula)
            data_i <- eval(out_i$call$data, envir = formula_env_i, enclos = parent.frame())
            idx_i <- eval(out_i$call$index, envir = formula_env_i, enclos = parent.frame())
            f_vars <- all.vars(out_i$call$formula)
            y_col_i <- f_vars[1]
            d_col_i <- f_vars[2]
            x_cols_i <- if (length(f_vars) > 2) f_vars[3:length(f_vars)] else character(0)

            rr_i <- match(data_i[[idx_i[2]]], out_i$rawtime)
            cc_i <- match(data_i[[idx_i[1]]], out_i$id)
            long_mask_i <- as.logical(mask_mat_rep[cbind(rr_i, cc_i)])
            hide_rows_i <- which(!is.na(long_mask_i) & long_mask_i)

            data_hidden <- data_i
            data_hidden[hide_rows_i, y_col_i] <- NA_real_

            formula_obj <- stats::reformulate(c(d_col_i, x_cols_i), response = y_col_i)
            rerun_args <- list(formula = formula_obj, data = data_hidden, index = idx_i)
            if (!is.null(out_i$call$method)) rerun_args$method <- eval(out_i$call$method, parent.frame())
            if (!is.null(out_i$call$r)) rerun_args$r <- eval(out_i$call$r, parent.frame())

            out_new <- do.call(fect::fect, rerun_args)
            fits[[i]] <- out_new

            pred <- out_new$Y.ct.full[cbind(rr_i[hide_rows_i], cc_i[hide_rows_i])]
            actual <- actual_mat_rep[cbind(rr_i[hide_rows_i], cc_i[hide_rows_i])]
            err <- pred - actual
            records <- rbind(
                records,
                data.frame(
                    Rep = rep_id,
                    Model = names(out_list)[i],
                    Hidden_N = length(err),
                    RMSE = sqrt(mean(err^2)),
                    Bias = mean(err),
                    stringsAsFactors = FALSE
                )
            )
        }

        fits_last <- fits
        mask_last <- mask_mat_rep
    }

    summary <- aggregate(cbind(Hidden_N, RMSE, Bias) ~ Model, data = records, FUN = mean)
    list(summary = summary, records = records, hide_mask = mask_last, fits = fits_last)
}
