fect_mspe <- function(
    x,
    data = NULL,
    index = c("id", "time"),
    models = list(fe = list(method = "fe")),
    hide_ids = NULL,
    hide_times = NULL,
    hide_mask = NULL,
    hide_n = 20,
    hide_last_periods = 5,
    hide_pick = "head",
    id_col = "id",
    time_col = "time",
    y_col = NULL,
    d_col = NULL,
    print_table = TRUE
) {
    hide_mask_provided <- !is.null(hide_mask)
    if (hide_mask_provided && (!missing(hide_n) || !missing(hide_last_periods))) {
        stop("When `hide_mask` is provided, do not pass `hide_n` or `hide_last_periods`.")
    }

    .is_fect_output <- function(obj) {
        is.list(obj) && !is.null(obj$Y.ct.full) && !is.null(obj$call)
    }

    .infer_formula_cols <- function(formula_obj, out_obj, y_col, d_col) {
        if (!is.null(formula_obj)) {
            vars <- all.vars(formula_obj)
            if (length(vars) >= 2) {
                if (is.null(y_col)) y_col <- vars[1]
                if (is.null(d_col)) d_col <- vars[2]
            }
        }
        if (is.null(y_col)) y_col <- out_obj$Y
        if (is.null(d_col)) d_col <- out_obj$D
        list(y_col = y_col, d_col = d_col)
    }

    .run_call_with_hidden_data <- function(call_obj, data_hidden) {
        call_args <- as.list(call_obj)[-1]
        call_args$data <- data_hidden
        do.call(fect::fect, call_args)
    }

    .build_hide_mask <- function(data, id_col, time_col, d_col,
        hide_ids, hide_times, hide_mask,
        hide_n, hide_last_periods, hide_pick) {
        if (!is.null(hide_mask)) {
            if (is.matrix(hide_mask)) {
                ids <- unique(data[[id_col]])
                times <- unique(data[[time_col]])
                m <- hide_mask

                if (all(dim(m) == c(length(times), length(ids)))) {
                    mask_mat <- m
                } else if (all(dim(m) == c(length(ids), length(times)))) {
                    mask_mat <- t(m)
                } else {
                    stop(
                        "`hide_mask` matrix must be time x id or id x time, matching data dimensions."
                    )
                }

                vals <- as.vector(mask_mat)
                if (!all(vals %in% c(0, 1, FALSE, TRUE))) {
                    stop("`hide_mask` matrix must contain only 0/1 (or FALSE/TRUE).")
                }

                row_idx <- match(data[[time_col]], times)
                col_idx <- match(data[[id_col]], ids)
                if (any(is.na(row_idx)) || any(is.na(col_idx))) {
                    stop("Failed to align `hide_mask` matrix with data id/time indices.")
                }

                return(as.logical(mask_mat[cbind(row_idx, col_idx)]))
            }

            if (length(hide_mask) != nrow(data)) {
                stop("`hide_mask` vector must have length nrow(data).")
            }
            vals <- as.vector(hide_mask)
            if (!all(vals %in% c(0, 1, FALSE, TRUE))) {
                stop("`hide_mask` vector must contain only 0/1 (or FALSE/TRUE).")
            }
            return(as.logical(hide_mask))
        }

        custom_hide <- !is.null(hide_ids) || !is.null(hide_times)

        if (custom_hide) {
            mask <- rep(TRUE, nrow(data))
            if (!is.null(hide_ids)) {
                mask <- mask & data[[id_col]] %in% hide_ids
            }
            if (!is.null(hide_times)) {
                mask <- mask & data[[time_col]] %in% hide_times
            }
            treated_ids <- unique(data[data[[d_col]] == 1, id_col])
            mask <- mask & data[[id_col]] %in% treated_ids
            if (!is.null(hide_n)) {
                cand_ids <- unique(data[mask, id_col])
                if (length(cand_ids) > hide_n) {
                    if (hide_pick == "tail") {
                        chosen_ids <- tail(cand_ids, hide_n)
                    } else if (hide_pick == "head") {
                        chosen_ids <- head(cand_ids, hide_n)
                    } else if (hide_pick == "random") {
                        chosen_ids <- sample(cand_ids, hide_n)
                    } else {
                        stop("`hide_pick` must be one of: \"tail\", \"head\", \"random\".")
                    }
                    mask <- mask & data[[id_col]] %in% chosen_ids
                }
            }
        } else {
            # Default: select a subset of treated units, then hide them in last periods.
            treated_ids <- unique(data[data[[d_col]] == 1, id_col])
            if (length(treated_ids) == 0) {
                return(rep(FALSE, nrow(data)))
            }
            uniq_time <- sort(unique(data[[time_col]]))
            last_times <- tail(uniq_time, hide_last_periods)
            if (is.null(hide_n)) {
                chosen_ids <- treated_ids
            } else {
                units_to_hide <- max(1, hide_n)
                if (units_to_hide >= length(treated_ids)) {
                    chosen_ids <- treated_ids
                } else {
                    if (hide_pick == "tail") {
                        chosen_ids <- tail(treated_ids, units_to_hide)
                    } else if (hide_pick == "head") {
                        chosen_ids <- head(treated_ids, units_to_hide)
                    } else if (hide_pick == "random") {
                        chosen_ids <- sample(treated_ids, units_to_hide)
                    } else {
                        stop("`hide_pick` must be one of: \"tail\", \"head\", \"random\".")
                    }
                }
            }
            mask <- data[[id_col]] %in% chosen_ids & data[[time_col]] %in% last_times
        }

        mask
    }

    .build_eval_mask <- function(data, d_col, time_col, hide_mask) {
        hide_times_used <- sort(unique(data[hide_mask, time_col]))
        if (length(hide_times_used) == 0) {
            stop("Cannot infer hide periods from `hide_mask`.")
        }
        data[[d_col]] == 1 & data[[time_col]] %in% hide_times_used
    }

    # Mode A: x is out.fect or list of out.fect -> auto hide and rerun same model call(s)
    if (!inherits(x, "formula")) {
        out_list_in <- NULL
        if (.is_fect_output(x)) {
            out_list_in <- list(x)
        } else if (is.list(x) && length(x) > 0 &&
            all(vapply(x, .is_fect_output, logical(1)))) {
            out_list_in <- x
        }
        if (is.null(out_list_in)) {
            stop("Input must be either a formula or a `fect` output object.")
        }

        if (is.null(names(out_list_in)) || any(names(out_list_in) == "")) {
            names(out_list_in) <- paste0("out.fect_", seq_along(out_list_in))
        }

        ref <- out_list_in[[1]]
        if (is.null(data)) {
            if (is.null(ref$call$data)) {
                stop("Cannot infer original `data` from `out.fect$call`. Please pass `data` explicitly.")
            }
            data <- eval(ref$call$data, envir = parent.frame())
        }

        formula_obj <- ref$call$formula
        cols <- .infer_formula_cols(formula_obj, ref, y_col, d_col)
        y_col <- cols$y_col
        d_col <- cols$d_col

        required_cols <- c(id_col, time_col, y_col, d_col)
        if (!all(required_cols %in% colnames(data))) {
            stop("`data` must include columns: ", paste(required_cols, collapse = ", "))
        }

        if (!"Y0" %in% colnames(data)) {
            stop("`data` must contain column `Y0` for ground truth.")
        }
        target_col <- "Y0"

        mask <- .build_hide_mask(
            data = data,
            id_col = id_col,
            time_col = time_col,
            d_col = d_col,
            hide_ids = hide_ids,
            hide_times = hide_times,
            hide_mask = hide_mask,
            hide_n = hide_n,
            hide_last_periods = hide_last_periods,
            hide_pick = hide_pick
        )
        if (sum(mask) == 0) {
            stop("No observations selected to hide. Adjust hide settings.")
        }

        data_hidden <- data
        data_hidden[mask, y_col] <- NA
        eval_mask <- .build_eval_mask(data, d_col, time_col, mask)
        actual <- as.numeric(data[eval_mask, target_col])
        eval_id <- data[eval_mask, id_col]
        eval_time <- data[eval_mask, time_col]

        results <- data.frame(
            Model = names(out_list_in),
            Hidden_N = as.integer(sum(mask)),
            Eval_N = as.integer(sum(eval_mask)),
            RMSE = NA_real_,
            Bias = NA_real_,
            stringsAsFactors = FALSE
        )
        fits <- vector("list", length(out_list_in))
        names(fits) <- names(out_list_in)

        for (i in seq_along(out_list_in)) {
            out_new <- .run_call_with_hidden_data(out_list_in[[i]]$call, data_hidden)
            fits[[i]] <- out_new

            col_idx <- match(eval_id, out_new$id)
            row_idx <- match(eval_time, out_new$rawtime)
            if (any(is.na(col_idx)) || any(is.na(row_idx))) {
                stop("Failed to align hidden id/time with prediction matrix for model: ", names(out_list_in)[i])
            }

            pred <- out_new$Y.ct.full[cbind(row_idx, col_idx)]
            err <- pred - actual
            if (any(is.na(err))) {
                stop("NA found in prediction errors for model: ", names(out_list_in)[i], ".")
            }
            results$RMSE[i] <- sqrt(mean(err^2))
            results$Bias[i] <- mean(err)
        }

        if (print_table) print(results, row.names = FALSE)
        return(list(
            summary = results,
            hide_mask = mask,
            hidden_data = data_hidden,
            actual = actual,
            fits = fits
        ))
    }

    # Mode B: x is formula -> hide and run model list
    formula <- x
    vars <- all.vars(formula)
    if (length(vars) < 2) {
        stop("`formula` must include outcome and treatment variables.")
    }
    if (is.null(y_col)) y_col <- vars[1]
    if (is.null(d_col)) d_col <- vars[2]

    required_cols <- c(id_col, time_col, y_col, d_col)
    if (!all(required_cols %in% colnames(data))) {
        stop("`data` must include columns: ", paste(required_cols, collapse = ", "))
    }

    if (!"Y0" %in% colnames(data)) {
        stop("`data` must contain column `Y0` for ground truth.")
    }
    target_col <- "Y0"
    if (!is.list(models) || length(models) == 0) {
        stop("`models` must be a non-empty named list of fect() argument lists.")
    }
    if (is.null(names(models)) || any(names(models) == "")) {
        names(models) <- paste0("model_", seq_along(models))
    }

    mask <- .build_hide_mask(
        data = data,
        id_col = id_col,
        time_col = time_col,
        d_col = d_col,
        hide_ids = hide_ids,
        hide_times = hide_times,
        hide_mask = hide_mask,
        hide_n = hide_n,
        hide_last_periods = hide_last_periods,
        hide_pick = hide_pick
    )
    if (sum(mask) == 0) {
        stop("No observations selected to hide. Adjust hide settings.")
    }

    data_hidden <- data
    data_hidden[mask, y_col] <- NA
    eval_mask <- .build_eval_mask(data, d_col, time_col, mask)
    actual <- as.numeric(data[eval_mask, target_col])
    eval_id <- data[eval_mask, id_col]
    eval_time <- data[eval_mask, time_col]

    results <- data.frame(
        Model = names(models),
        Hidden_N = as.integer(sum(mask)),
        Eval_N = as.integer(sum(eval_mask)),
        RMSE = NA_real_,
        Bias = NA_real_,
        stringsAsFactors = FALSE
    )
    fits <- vector("list", length(models))
    names(fits) <- names(models)

    for (i in seq_along(models)) {
        fit_call <- c(list(formula = formula, data = data_hidden, index = index), models[[i]])
        out <- do.call(fect, fit_call)
        fits[[i]] <- out

        col_idx <- match(eval_id, out$id)
        row_idx <- match(eval_time, out$rawtime)
        if (any(is.na(col_idx)) || any(is.na(row_idx))) {
            stop("Failed to align hidden id/time with prediction matrix for model: ", names(models)[i])
        }

        pred <- out$Y.ct.full[cbind(row_idx, col_idx)]
        err <- pred - actual
        if (any(is.na(err))) {
            stop("NA found in prediction errors for model: ", names(models)[i], ".")
        }
        results$RMSE[i] <- sqrt(mean(err^2))
        results$Bias[i] <- mean(err)
    }

    if (print_table) print(results, row.names = FALSE)
    list(
        summary = results,
        hide_mask = mask,
        hidden_data = data_hidden,
        actual = actual,
        fits = fits
    )
}
