esplot <- function(data,  # time, ATT, CI.lower, CI.upper, count, ...
                   Period = NULL,
                   Estimate = "ATT",
                   SE = NULL,
                   CI.lower = "CI.lower", # Default name to look for
                   CI.upper = "CI.upper", # Default name to look for
                   Count = NULL,

                   proportion = 0.3,
                   est.lwidth = NULL,  # line thickness
                   est.pointsize = NULL,    # point size if show.points=TRUE
                   show.points = FALSE,  # for connected=TRUE, show points only at integer times
                   fill.gap = TRUE,      # fill missing times with 0
                   start0 = FALSE,       # if TRUE => vertical dashed line at x=-0.5
                   only.pre = FALSE,
                   only.post = FALSE,
                   show.count = NULL,    # whether to show count bars
                   stats = NULL,         # numeric p-values
                   stats.labs = NULL,    # labels for each p-value

                   highlight.periods = NULL,  # numeric vector of times to highlight
                   highlight.colors  = NULL,    # color for each highlight time
                   lcolor = NULL,
                   lwidth = NULL,
                   ltype = c("solid", "solid"),
                   connected = FALSE,    # if TRUE => line + ribbon
                   ci.outline = FALSE,
                   main      = NULL,
                   xlim      = NULL,
                   ylim      = NULL,
                   xlab      = NULL,
                   ylab      = NULL,
                   gridOff   = FALSE,
                   stats.pos = NULL,
                   theme.bw  = TRUE,
                   cex.main  = NULL,
                   cex.axis  = NULL,
                   cex.lab   = NULL,
                   cex.text  = NULL,
                   axis.adjust = FALSE,
                   color       = "#000000",
                   count.color = "gray70",
                   count.alpha = 0.4,
                   count.outline.color = "grey69"
)
{
  if (is.null(est.lwidth) || is.null(est.pointsize)) {
    default_est_lwidth <- .8
    default_est_pointsize <- 3

    if (!connected) {
      default_est_lwidth <- 0.6
      default_est_pointsize <- 2
    } else {
      if (show.points) {
        default_est_lwidth <- 0.7
        default_est_pointsize <- 1.2
      } else {
        default_est_lwidth <- 1.2
        default_est_pointsize <- 3
      }
    }
    if(is.null(est.lwidth)) est.lwidth <- default_est_lwidth
    if(is.null(est.pointsize)) est.pointsize <- default_est_pointsize
  }

  all_integer_like <- function(x) {
    x_chr <- if (is.factor(x)) as.character(x) else x
    all(grepl("^[-+]?[0-9]+$", x_chr))
  }
  # If input is a did_wrapper object, extract the event-study data
  if (inherits(data, "did_wrapper")) {
    data <- data$est.att
  }
  data <- as.data.frame(data)

  # Identify time/period column
  if (is.null(Period)) {
    common_period_names <- c("time", "Time", "period", "Period", "event.time", "event_time", "rel_time")
    found_period <- common_period_names[common_period_names %in% names(data)]
    if (length(found_period) > 0) {
      Period <- found_period[1]
    } else if (all_integer_like(rownames(data))) {
      data$Period_generated <- as.numeric(rownames(data))
      Period <- "Period_generated"
    } else {
      stop("Period column not specified and could not be inferred. Please specify the 'Period' argument.")
    }
  }
  if (!(Period %in% names(data))) stop(paste0("Period column '", Period, "' not found in data."))


  # Ensure Estimate column exists and is numeric
  if (!(Estimate %in% names(data))) stop(paste("Estimate column '", Estimate, "' not found in data."))
  data[[Estimate]] <- as.numeric(data[[Estimate]])

  # Ensure SE column is numeric if it exists and is specified
  if (!is.null(SE) && SE %in% names(data)) {
    data[[SE]] <- as.numeric(data[[SE]])
  }

  # --- START: Centralized CI calculation from SE ---
  # Check for CI.lower. If not present, try to calculate from SE.
  # CI.lower is the parameter name (e.g., "CI.lower" by default).
  if (!(CI.lower %in% names(data))) {
    if (!is.null(SE) && SE %in% names(data)) {
      calculated_ci_lower <- data[[Estimate]] - 1.96 * data[[SE]]
      if (any(!is.na(calculated_ci_lower))) {
        message(paste0("Column '", CI.lower, "' not found in input data. Calculating values for '", CI.lower, "' using '", Estimate, "' and '", SE, "'."))
      }
      data[[CI.lower]] <- calculated_ci_lower
    } else {
      data[[CI.lower]] <- NA_real_ # Create the column with NAs
    }
  } else { # CI.lower column exists in input data
    data[[CI.lower]] <- as.numeric(data[[CI.lower]]) # Ensure it's numeric
  }

  # Check for CI.upper. If not present, try to calculate from SE.
  if (!(CI.upper %in% names(data))) {
    if (!is.null(SE) && SE %in% names(data)) {
      calculated_ci_upper <- data[[Estimate]] + 1.96 * data[[SE]]
      if (any(!is.na(calculated_ci_upper))) {
        message(paste0("Column '", CI.upper, "' not found in input data. Calculating values for '", CI.upper, "' using '", Estimate, "' and '", SE, "'."))
      }
      data[[CI.upper]] <- calculated_ci_upper
    } else {
      data[[CI.upper]] <- NA_real_
    }
  } else { # CI.upper column exists in input data
    data[[CI.upper]] <- as.numeric(data[[CI.upper]]) # Ensure it's numeric
  }
  # --- END: Centralized CI calculation from SE ---
  # Now, 'data' is guaranteed to have columns named by CI.lower and CI.upper parameters,
  # containing original CIs, CIs from SE, or NAs. These columns are also numeric.

  # Possibly build an interpolated dataset if connected=TRUE
  if (connected) {
    min_period_val <- min(data[[Period]], na.rm = TRUE)
    max_period_val <- max(data[[Period]], na.rm = TRUE)

    if (!is.finite(min_period_val) || !is.finite(max_period_val)) {
      stop("Cannot interpolate with non-finite Period values.")
    }

    new_periods <- seq(
      from = min_period_val,
      to   = max_period_val,
      by   = 0.5
    )

    # SE for interpolation: use data[[SE]] if SE column exists, otherwise NA
    # data[[SE]] is already numeric or NA due to upfront logic.
    se_for_interp <- if (!is.null(SE) && SE %in% names(data)) {
      data[[SE]]
    } else {
      rep(NA_real_, nrow(data))
    }

    # CI columns are guaranteed to exist in 'data' and be numeric (or NA_real_).
    ci_lower_for_interp <- data[[CI.lower]]
    ci_upper_for_interp <- data[[CI.upper]]

    count_for_interp <- if (!is.null(Count) && Count %in% names(data)) {
      data[[Count]]
    } else {
      rep(NA_real_, nrow(data))
    }

    temp_data_for_interp <- data.frame(
      Period_val = data[[Period]],
      Estimate_val = data[[Estimate]],      # Already numeric
      SE_val = se_for_interp,              # Numeric or NA
      CI.lower_val = ci_lower_for_interp,  # Numeric or NA
      CI.upper_val = ci_upper_for_interp,  # Numeric or NA
      Count_val = count_for_interp
    )

    interp_col_safe <- function(y_vals) {
      if (all(is.na(y_vals))) return(rep(NA_real_, length(new_periods)))
      idx_valid_approx <- !is.na(temp_data_for_interp$Period_val) & !is.na(y_vals)
      if(sum(idx_valid_approx) < 2) return(rep(NA_real_, length(new_periods)))
      stats::approx(x = temp_data_for_interp$Period_val[idx_valid_approx],
                    y = y_vals[idx_valid_approx],
                    xout = new_periods, rule = 2)$y
    }

    df_interp <- data.frame(Period_placeholder = new_periods)
    df_interp[[Estimate]] <- interp_col_safe(temp_data_for_interp$Estimate_val)

    # Interpolate SE if SE parameter was provided
    if (!is.null(SE)) {
      # temp_data_for_interp$SE_val contains original SEs or NAs if SE column didn't exist/was all NA.
      # interp_col_safe will handle all-NA y_vals correctly.
      df_interp[[SE]] <- interp_col_safe(temp_data_for_interp$SE_val)
    }

    # Interpolate CIs. CI.lower and CI.upper are the names of the columns.
    df_interp[[CI.lower]] <- interp_col_safe(temp_data_for_interp$CI.lower_val)
    df_interp[[CI.upper]] <- interp_col_safe(temp_data_for_interp$CI.upper_val)


    if (!is.null(Count)) {
      df_interp[[Count]] <- NA_real_
      orig_t <- temp_data_for_interp$Period_val
      orig_c <- temp_data_for_interp$Count_val
      for (i in seq_len(nrow(df_interp))) {
        pval <- df_interp$Period_placeholder[i]
        if (abs(pval - round(pval)) < 1e-8) {
          idx_match <- which(abs(orig_t - pval) < 1e-8 & !is.na(orig_c))
          if (length(idx_match) >= 1) {
            df_interp[[Count]][i] <- orig_c[idx_match[1]]
          }
        }
      }
    }
    names(df_interp)[names(df_interp) == "Period_placeholder"] <- Period
    data <- df_interp
  } # End of if(connected)

  # Ensure relevant columns are numeric after potential interpolation.
  # Upfront logic and interpolation should handle this, but this is a safeguard.
  if (!is.null(SE) && SE %in% names(data)) data[[SE]] <- as.numeric(data[[SE]])
  # CI.lower and CI.upper columns are guaranteed to exist by upfront logic.
  data[[CI.lower]] <- as.numeric(data[[CI.lower]])
  data[[CI.upper]] <- as.numeric(data[[CI.upper]])


  # Check "show.count"
  if (is.null(show.count)) {
    show.count <- FALSE
  }
  if (!is.logical(show.count) && !is.numeric(show.count)) {
    stop("\"show.count\" must be TRUE/FALSE or NULL.")
  }
  if (show.count && is.null(Count)) {
    stop("\"Count\" column name is not specified via 'Count' parameter, but show.count is TRUE.")
  }
  if (show.count && !(Count %in% names(data))) {
    stop(paste0("Specified \"Count\" column '", Count, "' is not in the data, but show.count is TRUE."))
  }

  # Check other logical flags
  if (!is.logical(gridOff) && !gridOff %in% c(0, 1)) {
    stop("\"gridOff\" must be TRUE/FALSE.")
  }
  if (!is.logical(fill.gap) && !fill.gap %in% c(0, 1)) {
    stop("\"fill.gap\" must be TRUE/FALSE.")
  }
  if (!is.logical(axis.adjust) && !axis.adjust %in% c(0, 1)) {
    stop("\"axis.adjust\" must be TRUE/FALSE.")
  }

  # Title & label sizing
  if (!is.null(main)) {
    if (!is.character(main)) stop("\"main\" is not a string.")
    main <- main[1]
  }
  if (!is.null(cex.main)) {
    if (!is.numeric(cex.main)) stop("\"cex.main\" must be numeric.")
    cex.main <- 16 * cex.main # ggplot title size is different from base R
  } else {
    cex.main <- 16
  }

  if (!is.null(xlab) && !is.character(xlab)) stop("\"xlab\" is not a string.")
  if (!is.null(ylab) && !is.character(ylab)) stop("\"ylab\" is not a string.")

  if (!is.null(cex.lab)) {
    if (!is.numeric(cex.lab)) stop("\"cex.lab\" must be numeric.")
    cex.lab <- 15 * cex.lab
  } else {
    cex.lab <- 15
  }

  if (!is.null(cex.axis)) {
    if (!is.numeric(cex.axis)) stop("\"cex.axis\" must be numeric.")
    cex.axis <- 15 * cex.axis
  } else {
    cex.axis <- 15
  }

  if (!is.null(cex.text)) { # For annotate text
    if (!is.numeric(cex.text)) stop("\"cex.text\" must be numeric.")
    cex.text <- 5 * cex.text # ggplot text size for annotate
  } else {
    cex.text <- 5
  }

  # Stats
  if (!is.null(stats)) {
    stats <- c(stats)
    if (!is.numeric(stats)) {
      stop("The \"stats\" option must be numeric. Calling function should extract numeric values.")
    }
    stats_display <- formatC(stats, digits = 3, format = "f", flag = "#") # Keep as is

    n.stats <- length(stats_display)
    if (is.null(stats.labs)) {
      stats.labs <- rep("", n.stats)
    } else if (length(stats.labs) != n.stats) {
      stop("The \"stats.labs\" option should have the same length as \"stats\".")
    }
  }
  # Validation for stats.pos (already present and correct)
  if (!is.null(stats.pos)) {
    if (length(stats.pos) != 2) {
      stop("\"stats.pos\" must be of length 2.")
    }
    if (!is.numeric(stats.pos[1]) || !is.numeric(stats.pos[2])) {
      stop("Elements of \"stats.pos\" must be numeric.")
    }
  }


  # Axis rotation
  if (axis.adjust) {
    angle <- 45
    x.v <- 1 # vjust for x-axis text
    x.h <- 1 # hjust for x-axis text
  } else {
    angle <- 0
    x.v <- 0.5
    x.h <- 0.5
  }

  # xlim, ylim checks (user-provided)
  user_xlim <- xlim
  user_ylim <- ylim
  if (!is.null(user_xlim)) {
    if (!is.numeric(user_xlim) || length(user_xlim) != 2) {
      stop("\"xlim\" must be a numeric vector of length 2.")
    }
    if(user_xlim[1] > user_xlim[2]) user_xlim <- rev(user_xlim) # Ensure min < max
  }
  if (!is.null(user_ylim)) {
    if (!is.numeric(user_ylim) || length(user_ylim) != 2) {
      stop("\"ylim\" must be a numeric vector of length 2.")
    }
    if(user_ylim[1] > user_ylim[2]) user_ylim <- rev(user_ylim) # Ensure min < max
  }

  # Default axis labels
  if (is.null(xlab)) {
    xlab <- "Time Relative to Treatment"
  } else if (xlab == "") {
    xlab <- NULL # ggplot prefers NULL for no label
  }
  if (is.null(ylab)) {
    ylab <- "Effect on Y"
  } else if (ylab == "") {
    ylab <- NULL
  }

  if (is.null(lcolor)) {
    if (theme.bw) {
      lcolor <- "#AAAAAA70"
    } else {
      lcolor <- "white" # This might be invisible on a white background if not theme_bw
    }
  }

  if (length(as.vector(lcolor)) == 1) {
    lcolor <- rep(lcolor, 2)
  } else if (length(as.vector(lcolor)) != 2) {
    stop("\"lcolor\" must be a numeric vector of length 1 or 2.")
  }
  if (is.null(lwidth)) {
    if (theme.bw) {
      lwidth <- 1.5
    } else {
      lwidth <- 2
    }
  }
  if (length(as.vector(lwidth)) == 1) {
    lwidth <- rep(lwidth, 2)
  } else if (length(as.vector(lwidth)) != 2) {
    stop("\"lwidth\" must be a numeric vector of length 1 or 2.")
  }

  # Prepare data vectors for plotting logic
  current_time_vec <- data[[Period]]
  current_att_vec  <- data[[Estimate]]

  # CI.lower and CI.upper columns are now guaranteed to exist in 'data'
  # and be numeric (or NA_real_). Their names are given by the CI.lower and CI.upper parameters.
  CI.lower.val <- data[[CI.lower]]
  CI.upper.val <- data[[CI.upper]]

  current_count_vec <- if (!is.null(Count) && Count %in% names(data)) {
    data[[Count]]
  } else {
    rep(NA_real_, length(current_time_vec))
  }

  # Determine initial final_xlim for plotting
  final_xlim <- user_xlim
  if (is.null(final_xlim)) {
    if (show.count && !is.null(Count) && Count %in% names(data) && any(!is.na(data[[Count]]))) {
      finite_counts <- data[[Count]][is.finite(data[[Count]])]
      if (length(finite_counts) > 0) {
        maxC <- max(finite_counts, na.rm = TRUE)
        if (is.finite(maxC) && maxC > 0) {
          threshold <- maxC * proportion
          idx_valid_count_times <- data[[Count]] >= threshold & !is.na(data[[Count]]) & is.finite(data[[Period]])
          valid_times <- data[[Period]][idx_valid_count_times]

          if (length(valid_times) > 0) {
            final_xlim <- range(valid_times, na.rm = TRUE)
          } else {
            idx_finite_period <- is.finite(data[[Period]])
            final_xlim <- range(data[[Period]][idx_finite_period], na.rm = TRUE)
          }
        } else {
          idx_finite_period <- is.finite(data[[Period]])
          final_xlim <- range(data[[Period]][idx_finite_period], na.rm = TRUE)
        }
      } else {
        idx_finite_period <- is.finite(data[[Period]])
        final_xlim <- range(data[[Period]][idx_finite_period], na.rm = TRUE)
      }
    } else {
      idx_finite_period <- is.finite(data[[Period]])
      final_xlim <- range(data[[Period]][idx_finite_period], na.rm = TRUE)
    }
    if (any(!is.finite(final_xlim))) {
      idx_finite_period <- is.finite(data[[Period]])
      valid_times <- data[[Period]][idx_finite_period]
      if(length(valid_times) > 0) final_xlim <- range(valid_times) else final_xlim <- c(0,1)
    }
    if (final_xlim[1] == final_xlim[2]) final_xlim <- final_xlim + c(-0.5, 0.5)
  }


  # fill.gap logic
  time_vec_for_fillgap <- current_time_vec
  att_vec_for_fillgap  <- current_att_vec
  ci_l_vec_for_fillgap <- CI.lower.val # These are from data[[CI.lower]]
  ci_u_vec_for_fillgap <- CI.upper.val # These are from data[[CI.upper]]
  count_vec_for_fillgap<- current_count_vec

  if (fill.gap && !connected) {
    idx_finite_time <- is.finite(time_vec_for_fillgap)
    min_time_fill <- min(time_vec_for_fillgap[idx_finite_time], na.rm = TRUE)
    max_time_fill <- max(time_vec_for_fillgap[idx_finite_time], na.rm = TRUE)

    if (is.finite(min_time_fill) && is.finite(max_time_fill) && min_time_fill <= max_time_fill) {
      idx_finite_unique_times <- !is.na(time_vec_for_fillgap) & is.finite(time_vec_for_fillgap)
      unique_times_present <- unique(time_vec_for_fillgap[idx_finite_unique_times])


      diff_times <- diff(sort(unique_times_present))
      idx_valid_diffs <- is.finite(diff_times) & diff_times > 1e-9
      valid_diffs <- diff_times[idx_valid_diffs]

      step_by <- if(length(valid_diffs) > 0 && all(abs(valid_diffs - round(valid_diffs)) < 1e-9)) {
        1
      } else if (length(valid_diffs) > 0) {
        min(valid_diffs, na.rm = TRUE)
      } else {
        1
      }

      if(!is.finite(step_by) || step_by <= 1e-9) step_by <- 1

      if ( (max_time_fill - min_time_fill)/step_by + 1 > length(unique_times_present) ) {
        full_seq <- seq(min_time_fill, max_time_fill, by = step_by)
        time.add <- setdiff(full_seq, unique_times_present)
        if (length(time.add) > 0) {
          att_vec_for_fillgap  <- c(att_vec_for_fillgap, rep(0, length(time.add)))
          ci_l_vec_for_fillgap <- c(ci_l_vec_for_fillgap, rep(0, length(time.add)))
          ci_u_vec_for_fillgap <- c(ci_u_vec_for_fillgap, rep(0, length(time.add)))
          count_vec_for_fillgap<- c(count_vec_for_fillgap, rep(NA, length(time.add)))
          time_vec_for_fillgap <- c(time_vec_for_fillgap, time.add)
        }
      }
    }
  }

  plot_data <- data.frame(
    time = time_vec_for_fillgap
  )
  plot_data[[Estimate]] <- att_vec_for_fillgap
  plot_data[[CI.lower]] <- ci_l_vec_for_fillgap # Use the parameter name as col name
  plot_data[[CI.upper]] <- ci_u_vec_for_fillgap # Use the parameter name as col name
  if(!is.null(Count)) plot_data[[Count]] <- count_vec_for_fillgap
  else plot_data[["count_placeholder"]] <- count_vec_for_fillgap

  names(plot_data)[names(plot_data) == "time"] <- Period


  # Apply only.pre/only.post modifications to final_xlim
  if (!is.null(final_xlim) && all(is.finite(final_xlim))){
    if (only.pre) {
      final_xlim[2] <- min(final_xlim[2], 0, na.rm = TRUE)
      if (start0) final_xlim[2] <- min(final_xlim[2], -0.5, na.rm=TRUE)
    } else if (only.post) {
      final_xlim[1] <- max(final_xlim[1], 0, na.rm = TRUE)
      if (!start0) final_xlim[1] <- max(final_xlim[1], 0.5, na.rm=TRUE)
    }
    if(final_xlim[1] > final_xlim[2]) final_xlim <- rev(final_xlim)
  }

  if (!is.null(final_xlim)) {
    needs_resolve_min <- is.na(final_xlim[1]) || !is.finite(final_xlim[1])
    needs_resolve_max <- is.na(final_xlim[2]) || !is.finite(final_xlim[2])

    idx_finite_times_plot_data <- is.finite(plot_data[[Period]])
    finite_times_in_plot_data <- plot_data[[Period]][idx_finite_times_plot_data]
    data_derived_min_time <- if(length(finite_times_in_plot_data) > 0) min(finite_times_in_plot_data, na.rm=TRUE) else 0
    data_derived_max_time <- if(length(finite_times_in_plot_data) > 0) max(finite_times_in_plot_data, na.rm=TRUE) else 1
    if (data_derived_min_time > data_derived_max_time) { data_derived_min_time <- 0; data_derived_max_time <- 1; }


    if (needs_resolve_min) final_xlim[1] <- data_derived_min_time
    if (needs_resolve_max) final_xlim[2] <- data_derived_max_time


    if (final_xlim[1] > final_xlim[2]) {
      final_xlim <- range(c(data_derived_min_time, data_derived_max_time, final_xlim))
    }
    if (final_xlim[1] == final_xlim[2]) {
      final_xlim <- final_xlim + c(-0.5, 0.5)
    }
    idx <- plot_data[[Period]] >= final_xlim[1] &
      plot_data[[Period]] <= final_xlim[2]
    plot_data <- plot_data[idx, ]

  }


  if(nrow(plot_data) == 0 && !(is.null(final_xlim) && is.null(user_ylim))) {
    warning("No data points remaining after applying xlim/filters. Plot will be empty or may error.")
    empty_df_structure <- data.frame(temp_period = numeric(0), temp_est = numeric(0), temp_cil = numeric(0), temp_ciu = numeric(0))
    names(empty_df_structure) <- c(Period, Estimate, CI.lower, CI.upper)
    if(!is.null(Count)) empty_df_structure[[Count]] <- numeric(0)
    plot_data <- empty_df_structure
  }

  # Determine final ylim for plotting
  final_ylim <- user_ylim
  if (is.null(final_ylim)) {
    y_min_series <- plot_data[[Estimate]]
    y_max_series <- plot_data[[Estimate]]
    # CI.lower and CI.upper columns are guaranteed to be in plot_data
    idx_finite_cil <- is.finite(plot_data[[CI.lower]])
    y_min_series <- c(y_min_series, plot_data[[CI.lower]][idx_finite_cil])

    idx_finite_ciu <- is.finite(plot_data[[CI.upper]])
    y_max_series <- c(y_max_series, plot_data[[CI.upper]][idx_finite_ciu])


    idx_finite_ymin <- is.finite(y_min_series)
    y_min_data <- min(y_min_series[idx_finite_ymin], na.rm = TRUE)
    idx_finite_ymax <- is.finite(y_max_series)
    y_max_data <- max(y_max_series[idx_finite_ymax], na.rm = TRUE)


    if (!is.finite(y_min_data) && !is.finite(y_max_data) && nrow(plot_data)==0) {
      y_min_data <- 0; y_max_data <- 1;
    } else {
      if (!is.finite(y_min_data)) y_min_data <- if(is.finite(y_max_data)) y_max_data -1 else 0
      if (!is.finite(y_max_data)) y_max_data <- if(is.finite(y_min_data)) y_min_data +1 else 1
    }
    if (y_min_data == y_max_data) {
      y_min_data <- y_min_data - 0.5
      y_max_data <- y_max_data + 0.5
    }

    range_val <- y_max_data - y_min_data
    if (!is.finite(range_val) || range_val == 0) range_val <- 1

    top_expand_factor <- 0
    bot_expand_factor <- 0

    p_label_text_for_ylim_calc <- function(stats_vals, labs_vals) {
      text <- ""
      s_display <- formatC(stats_vals, digits = 3, format = "f", flag = "#")
      for (k in seq_along(s_display)) {
        text <- paste0(text, labs_vals[k], ifelse(labs_vals[k]=="", "", ": "), s_display[k], "\n")
      }
      sub("\n$", "", text)
    }
    if (!is.null(stats)) {
      num_stat_lines <- ceiling(length(unlist(strsplit(p_label_text_for_ylim_calc(stats, stats.labs), "\n"))))
      # top_expand_factor <- top_expand_factor # No change, already 0
    }
    if (show.count && !is.null(Count) && Count %in% names(plot_data)) {
      idx_valid_count_plot_data <- !is.na(plot_data[[Count]]) & plot_data[[Count]] > 0
      if (any(idx_valid_count_plot_data)) {
        bot_expand_factor <- bot_expand_factor + 0.15
      }
    }

    final_ylim <- c(
      y_min_data - bot_expand_factor * range_val,
      y_max_data + top_expand_factor * range_val
    )
    if(any(!is.finite(final_ylim))) {
      final_ylim <- c(y_min_data, y_max_data)
      if(final_ylim[1] == final_ylim[2]) final_ylim <- final_ylim + c(-0.5, 0.5)
    }
  }


  p <- ggplot(data = plot_data)


  if (theme.bw) {
    p <- p + theme_bw()
  }
  if (gridOff) {
    p <- p + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }

  p <- p + geom_hline(yintercept = 0, colour = lcolor[1], linewidth = lwidth[1], linetype = ltype[1])

  vline_pos <- if (start0) -0.5 else 0.5
  show_vline <- TRUE
  if (only.pre || only.post) show_vline <- FALSE
  if (!is.null(final_xlim) && all(is.finite(final_xlim))) {
    if (vline_pos < final_xlim[1] || vline_pos > final_xlim[2]) {
      show_vline <- FALSE
    }
  } else if (is.null(final_xlim)) {
    if (nrow(plot_data) > 0 && Period %in% names(plot_data)) {
      idx_finite_period_plot_data <- is.finite(plot_data[[Period]])
      if (any(idx_finite_period_plot_data)) {
        data_period_range <- range(plot_data[[Period]][idx_finite_period_plot_data], na.rm = TRUE)
        if (vline_pos < data_period_range[1] || vline_pos > data_period_range[2]) show_vline <- FALSE
      } else {
        show_vline <- FALSE
      }
    } else {
      show_vline <- FALSE
    }
  }


  if (show_vline) {
    p <- p + geom_vline(
      xintercept = vline_pos, colour = lcolor[2],
      linewidth = lwidth[2], linetype = ltype[2]
    )
  }


  intervals <- NULL
  if (!is.null(highlight.periods) && length(highlight.periods) > 0) {
    if (is.null(highlight.colors) || length(highlight.colors) != length(highlight.periods)) {
      warning("highlight.colors not provided or mismatched length; using default rainbow colors.")
      highlight.colors <- grDevices::rainbow(length(highlight.periods))
    }
    intervals <- data.frame(
      start = highlight.periods - 0.5,
      end   = highlight.periods + 0.5,
      color = highlight.colors
    )
  }

  remove_strictly_inside <- function(df, intervals_df, estimate_col, ci_l_col, ci_u_col, period_col) {
    if (is.null(intervals_df) || nrow(intervals_df) == 0) return(df)
    out_df <- df
    for (i in seq_len(nrow(intervals_df))) {
      s <- intervals_df$start[i]
      e <- intervals_df$end[i]
      idx_to_na <- which(out_df[[period_col]] > s & out_df[[period_col]] < e)
      if(length(idx_to_na) > 0) {
        out_df[[estimate_col]][idx_to_na] <- NA
        if(ci_l_col %in% names(out_df)) out_df[[ci_l_col]][idx_to_na] <- NA
        if(ci_u_col %in% names(out_df)) out_df[[ci_u_col]][idx_to_na] <- NA
      }
    }
    return(out_df)
  }

  if (connected) {
    data_base <- plot_data
    if (!is.null(intervals) && nrow(intervals) > 0) {
      data_base <- remove_strictly_inside(data_base, intervals, Estimate, CI.lower, CI.upper, Period)
    }
    p <- p +
      geom_ribbon(
        data = data_base, aes(x = .data[[Period]], ymin = .data[[CI.lower]], ymax = .data[[CI.upper]]),
        fill = color, alpha = 0.2, inherit.aes = FALSE, na.rm = TRUE,
        color =  ifelse(ci.outline, grDevices::adjustcolor(color, offset = c(0.3, 0.3, 0.3, 0)), NA)) +
      geom_line(
        data = data_base, aes(x = .data[[Period]], y = .data[[Estimate]]),
        linewidth = est.lwidth, color = color, inherit.aes = FALSE, na.rm = TRUE)
    if (show.points) {
      idx_int_points <- abs(as.numeric(data_base[[Period]]) - round(as.numeric(data_base[[Period]]))) < 1e-8 & !is.na(data_base[[Estimate]])
      data_base_int <- data_base[idx_int_points, ]

      if(nrow(data_base_int) > 0) {
        p <- p + geom_point(
          data = data_base_int, aes(x = .data[[Period]], y = .data[[Estimate]]),
          size = est.pointsize, color = color, inherit.aes = FALSE, na.rm = TRUE)
      }
    }
    if (!is.null(intervals) && nrow(intervals) > 0) {
      for (i in seq_len(nrow(intervals))) {
        sub_color <- intervals$color[i]; range_start <- intervals$start[i]; range_end <- intervals$end[i]
        idx_sub_data <- plot_data[[Period]] >= range_start & plot_data[[Period]] <= range_end
        sub.data <- plot_data[idx_sub_data, ]

        if(nrow(sub.data) > 0){
          p <- p + geom_ribbon(data = sub.data, aes(x = .data[[Period]], ymin = .data[[CI.lower]], ymax = .data[[CI.upper]]), inherit.aes = FALSE, na.rm = TRUE, fill = sub_color, alpha = 0.2, color =  ifelse(ci.outline, grDevices::adjustcolor(sub_color, offset = c(0.3,0.3,0.3,0)), NA)) +
            geom_line(data = sub.data, aes(x = .data[[Period]], y = .data[[Estimate]]), inherit.aes = FALSE, na.rm = TRUE, color = sub_color, linewidth = est.lwidth)
          if(show.points){
            idx_sub_data_int <- abs(as.numeric(sub.data[[Period]]) - round(as.numeric(sub.data[[Period]]))) < 1e-8 & !is.na(sub.data[[Estimate]])
            sub.data.int <- sub.data[idx_sub_data_int, ]

            if(nrow(sub.data.int) > 0) p <- p + geom_point(data = sub.data.int, aes(x = .data[[Period]], y = .data[[Estimate]]), size = est.pointsize, color = sub_color, inherit.aes = FALSE, na.rm = TRUE)
          }
        }
      }
    }
  } else { # Not connected
    p <- p + geom_pointrange(
      aes(x = .data[[Period]], y = .data[[Estimate]], ymin = .data[[CI.lower]], ymax = .data[[CI.upper]]),
      lwd = est.lwidth, color = color, fill = color, fatten = est.pointsize, na.rm = TRUE)

    if (!is.null(highlight.periods) && length(highlight.periods) > 0) {
      highlight_points_df <- data.frame(period_val = highlight.periods, color_val = highlight.colors)
      for (i in seq_len(nrow(highlight_points_df))) {
        hp_period <- highlight_points_df$period_val[i]; hp_color  <- highlight_points_df$color_val[i]
        idx_sub_data_point <- as.numeric(plot_data[[Period]]) == hp_period
        sub_data_point <- plot_data[idx_sub_data_point, ]

        if(nrow(sub_data_point) > 0) {
          p <- p + geom_pointrange(
            data = sub_data_point, aes(x = .data[[Period]], y = .data[[Estimate]], ymin = .data[[CI.lower]], ymax = .data[[CI.upper]]),
            lwd = est.lwidth, color = hp_color, fill = hp_color, fatten = est.pointsize, inherit.aes = FALSE, na.rm = TRUE)
        }
      }
    }
  }

  p <- p +
    xlab(xlab) + ylab(ylab) +
    theme(
      axis.title  = element_text(size = cex.lab),
      axis.text   = element_text(color = "black", size = cex.axis),
      axis.text.x = element_text(angle = angle, vjust = x.v, hjust = x.h),
      plot.title  = element_text(size = cex.main, hjust = 0.5, face = "bold")
    )



  if (show.count && !is.null(Count) && Count %in% names(plot_data)) {
    idx_rect_data <- !is.na(plot_data[[Count]]) & plot_data[[Count]] > 0
    rect_data <- plot_data[idx_rect_data, ]

    if(nrow(rect_data) > 0) {
      rect_data[[Period]] <- as.numeric(rect_data[[Period]])

      rect_data[,"xmin"] <- rect_data[[Period]] - 0.2
      rect_data[,"xmax"] <- rect_data[[Period]] + 0.2

      plot_range_y <- final_ylim[2] - final_ylim[1]
      if(!is.finite(plot_range_y) || plot_range_y <=0) plot_range_y <- 1

      rect_max_height_prop <- 0.1
      rect_bar_max_h <- plot_range_y * rect_max_height_prop
      rect.min_yval  <- final_ylim[1]

      max_rect_count <- max(rect_data[[Count]], na.rm = TRUE)

      if (is.finite(max_rect_count) && max_rect_count > 0) {
        rect_data[,"ymin"] <- rect.min_yval
        rect_data[,"ymax"] <- rect.min_yval + (rect_data[[Count]] / max_rect_count) * rect_bar_max_h
        p <- p + geom_rect(
          data = rect_data, aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
          fill = count.color, linewidth = 0.2, inherit.aes = FALSE, alpha = count.alpha, color = count.outline.color)

        max_count_val <- max_rect_count
        idx_max_count <- which(rect_data[[Count]] == max_count_val)
        max_count_time_pos <- rect_data[[Period]][idx_max_count[1]]

        if (length(idx_max_count) > 1) {
          preferred_pos <- if (start0) -1 else 0
          idx_preferred_pos_in_max <- rect_data[[Period]][idx_max_count] == preferred_pos
          if (any(idx_preferred_pos_in_max)) {
            max_count_time_pos <- preferred_pos
          } else {
            if(!is.null(final_xlim) && all(is.finite(final_xlim))){
              xlim_center <- mean(final_xlim)
              max_count_time_pos <- rect_data[[Period]][idx_max_count[which.min(abs(rect_data[[Period]][idx_max_count] - xlim_center))]]
            }
          }
        }
        p <- p + annotate("text",
                          x = max_count_time_pos,
                          y = rect.min_yval + rect_bar_max_h + 0.02 * plot_range_y,
                          label = format(max_count_val, big.mark=","),
                          size = cex.text * 0.7,
                          hjust = 0.5)
      }
    }
  }

  if (!is.null(stats)) {
    actual_stats_pos <- stats.pos
    stats_hjust_val <- 0
    stats_vjust_val <- 1

    if (is.null(actual_stats_pos)) {
      x_min_for_stats <- if(!is.null(final_xlim) && is.finite(final_xlim[1])) final_xlim[1] else {
        if(nrow(plot_data) > 0 && Period %in% names(plot_data)) {
          idx_finite_period_plot_data <- is.finite(plot_data[[Period]])
          if(any(idx_finite_period_plot_data)) min(plot_data[[Period]][idx_finite_period_plot_data], na.rm = TRUE) else 0
        } else 0
      }
      x_max_for_stats <- if(!is.null(final_xlim) && is.finite(final_xlim[2])) final_xlim[2] else {
        if(nrow(plot_data) > 0 && Period %in% names(plot_data)) {
          idx_finite_period_plot_data <- is.finite(plot_data[[Period]])
          if(any(idx_finite_period_plot_data)) max(plot_data[[Period]][idx_finite_period_plot_data], na.rm = TRUE) else 1
        } else 1
      }
      if (x_min_for_stats >= x_max_for_stats) {
        x_max_for_stats <- x_min_for_stats + 1
      }


      y_min_for_stats <- if(is.finite(final_ylim[1])) final_ylim[1] else {
        y_coords <- c(plot_data[[Estimate]], plot_data[[CI.lower]], plot_data[[CI.upper]])
        idx_finite_y_coords <- is.finite(y_coords)
        y_coords <- y_coords[idx_finite_y_coords]
        if(length(y_coords) > 0) min(y_coords, na.rm = TRUE) else 0
      }
      y_max_for_stats <- if(is.finite(final_ylim[2])) final_ylim[2] else {
        y_coords <- c(plot_data[[Estimate]], plot_data[[CI.lower]], plot_data[[CI.upper]])
        idx_finite_y_coords <- is.finite(y_coords)
        y_coords <- y_coords[idx_finite_y_coords]
        if(length(y_coords) > 0) max(y_coords, na.rm = TRUE) else 1
      }
      if (y_min_for_stats >= y_max_for_stats) {
        y_max_for_stats <- y_min_for_stats + 1
      }

      current_x_range <- x_max_for_stats - x_min_for_stats
      if (!is.finite(current_x_range) || current_x_range <= 1e-9) current_x_range <- 1

      current_y_range <- y_max_for_stats - y_min_for_stats
      if (!is.finite(current_y_range) || current_y_range <= 1e-9) current_y_range <- 1

      padding_x_factor <- 0.02
      padding_y_factor <- 0.02

      actual_stats_pos <- c(
        x_min_for_stats + padding_x_factor * current_x_range,
        y_max_for_stats - padding_y_factor * current_y_range
      )

      if(!is.finite(actual_stats_pos[1])) actual_stats_pos[1] <- x_min_for_stats
      if(!is.finite(actual_stats_pos[2])) actual_stats_pos[2] <- y_max_for_stats
    }

    p_label_text <- ""
    for (ii in seq_along(stats_display)) {
      p_label_text <- paste0(p_label_text, stats.labs[ii], ifelse(stats.labs[ii]=="", "", ": "), stats_display[ii], "\n")
    }
    p_label_text <- sub("\n$", "", p_label_text)

    p <- p + annotate("text",
                      x = actual_stats_pos[1], y = actual_stats_pos[2],
                      label = p_label_text,
                      size = cex.text,
                      hjust = stats_hjust_val,
                      vjust = stats_vjust_val,
                      lineheight = .8)
  }


  if (is.null(main)) {
    p <- p + ggtitle("Estimated Dynamic Treatment Effects")
  } else if (main != "") {
    p <- p + ggtitle(main)
  }
  p <- p + coord_cartesian(xlim = c(final_xlim[1]-0.2,final_xlim[2]+0.2), ylim = final_ylim)
  p <- p + scale_x_continuous(
    breaks = function(lims) {
      ## 1. ggplot’s own suggestion (needs {scales})
      br <- scales::pretty_breaks()(lims)

      ## 2. keep only integers
      br <- br[abs(br - round(br)) < 1e-8]

      ## 3. if none left (very narrow panel or weird limits), make 5 integers
      if (length(br) == 0) {
        br <- round(seq(lims[1], lims[2], length.out = 5))
      }
      br
    }
  )


  return(p)
}
