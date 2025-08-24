## an auxiliary function for estimating cumulative or sub-group treatment effect based on fect

## 1. cumulative effect for a specified event window
## 2. averaging treatment effect in a sub-group

effect <- function(x, ## a fect object
                   cumu = TRUE, ## whether to calculate cumulative effect
                   id = NULL, ## units to be averaged on
                   period = NULL, ## event window
                   plot = FALSE,
                   count = TRUE,
                   xlab = NULL,
                   ylab = NULL,
                   main = NULL) {
  if (is.null(x$eff.boot)){
    stop("No bootstrap/jackknife results available. Choose keep.sims = TRUE in fect().")
  }
  if (x$hasRevs) {
    warning("Cumulative effects are not well-defined for panels with treatment reversals.")
    return() # Return NA
  }

  # Select units for analysis
  if (is.null(id)) {
    # If no specific units provided, select all treated units
    mask <- (colSums(x$D.dat) > 0)
  } else {
    # Otherwise, select specified units
    mask <- (colnames(x$eff) %in% id)
  }

  # Extract relevant matrices for selected units
  eff <- x$eff[, mask]      # Treatment effects
  D <- x$D.dat[, mask]      # Treatment indicators
  I <- x$I.dat[, mask]      # Inclusion (non-missing) indicators
  inference <- x$call$vartype  # Inference type
  method <- x$method        # Method

  # Get dimensions of data
  TT <- dim(eff)[1]  # Number of time periods
  N <- dim(eff)[2]   # Number of units

  # Replace NA values in treatment indicators with zeros
  if (sum(is.na(c(D))) > 0) {
    D[which(is.na(D))] <- 0
  }

  # Store original treatment matrix
  D.old <- D

  # Calculate cumulative treatment for determining time periods
  D <- apply(D, 2, function(vec) {
    cumsum(vec)
  })

  # Find minimum pre-treatment period across all units
  minT0 <- min(apply(D, 2, function(vec) {
    sum(vec == 0)
  }))

  # Restore original treatment matrix
  D <- D.old

  # Set default period range if not provided
  period.raw <- c(1, TT - minT0)
  if (!is.null(period)) {
    # Validate user-provided period
    if (period[2] > period.raw[2]) {
      stop(paste("Ending period should not be greater than ", period.raw[2], sep = ""))
    }
  } else {
    period <- period.raw
  }

  # Calculate cumulative average treatment effect
  catt <- getEffect(D, I, eff, cumu, period)

  # Initialize bootstrap results
  catt.boot <- NULL

  # Check if uncertainty estimates are available
  if (is.null(x$est.avg)) {
    cat("No uncertainty estimates.")
  } else {
    # Perform bootstrap analysis
    nboots <- length(x$att.avg.boot)
    catt.boot <- matrix(NA, period[2] - period[1] + 1, nboots)

    # Calculate treatment effect for each bootstrap sample
    for (i in 1:nboots) {
      # Extract bootstrap matrices
      D.boot <- x$D.boot[, , i]
      I.boot <- x$I.boot[, , i]
      eff.boot <- x$eff.boot[, , i]

      # Select treated units in bootstrap sample
      if (is.null(id)){
        mask.boot <- (colSums(D.boot) > 0)
      } else {
        mask.boot <- (x$id[unlist(x$colnames.boot[i])] %in% id)
      }

      # Extract relevant matrices for selected units
      Itr.boot <- I.boot[, mask.boot]
      Dtr.boot <- D.boot[, mask.boot]
      eff.tr.boot <- eff.boot[, mask.boot]

      # Calculate treatment effect for this bootstrap sample
      catt.boot[, i] <- getEffect(
        as.matrix(Dtr.boot),
        as.matrix(Itr.boot),
        as.matrix(eff.tr.boot),
        cumu, period)
    }
  }

  # Function to calculate p-values for non-parametric test
  get.pvalue <- function(vec) {
    if (NaN %in% vec | NA %in% vec) {
      # Handle NaN and NA values
      nan.pos <- is.nan(vec)
      na.pos <- is.na(vec)
      pos <- c(which(nan.pos), which(na.pos))
      vec.a <- vec[-pos]
      # Calculate two-sided p-values
      a <- sum(vec.a >= 0) / (length(vec) - sum(nan.pos | na.pos)) * 2
      b <- sum(vec.a <= 0) / (length(vec) - sum(nan.pos | na.pos)) * 2
    } else {
      # Calculate two-sided p-values for complete data
      a <- sum(vec >= 0) / length(vec) * 2
      b <- sum(vec <= 0) / length(vec) * 2
    }
    # Return minimum p-value, capped at 1
    return(min(as.numeric(min(a, b)), 1))
  }

  # Calculate confidence intervals and statistics if bootstrap results exist
  if (!is.null(catt.boot)) {
    # Check if inference method is jackknife
    is_jackknife <- !is.null(inference) && inference == "jackknife"
    is_parametric <- !is.null(inference) && inference == "parametric"

    # Calculate standard errors with proper scaling for jackknife
    if (is_jackknife) {
      # For jackknife, scale by sqrt(N-1)
      N_samples <- ncol(catt.boot)
      jackknife_scale <- sqrt(N_samples - 1)
      se.att <- apply(catt.boot, 1, function(vec) sd(vec, na.rm = TRUE) * jackknife_scale)
    } else {
      # Standard calculation for bootstrap
      se.att <- apply(catt.boot, 1, function(vec) sd(vec, na.rm = TRUE))
    }

    # Calculate 95% confidence intervals
    if (is_jackknife || is_parametric) {
      # For jackknife, use t-distribution with N-1 degrees of freedom
      N_samples <- ncol(catt.boot)
      t_critical <- stats::qt(0.975, df = N_samples - 1)
      CI.att <- t(apply(cbind(catt, se.att), 1, function(row) {
        c(row[1] - t_critical * row[2], row[1] + t_critical * row[2])
      }))
    } else {
      # For bootstrap, use empirical quantiles
      CI.att <- t(apply(catt.boot, 1, function(vec) {
        quantile(vec, c(0.025, 0.975), na.rm = TRUE)
      }))
    }

    # Calculate p-values
    if (is_jackknife || is_parametric) {
      # For jackknife, use t-distribution for p-values
      N_samples <- ncol(catt.boot)
      pvalue.att <- sapply(1:nrow(catt.boot), function(i) {
        t_stat <- catt[i] / se.att[i]
        2 * stats::pt(-abs(t_stat), df = N_samples - 1)
      })
    } else {
      # For bootstrap, use empirical distribution
      pvalue.att <- apply(catt.boot, 1, get.pvalue)
    }
    # Combine results into a matrix
    est.catt <- cbind(catt, se.att, CI.att, pvalue.att)
    colnames(est.catt) <- c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value")
    rownames(est.catt) <- period[1]:period[2]
  }

  # Plot
  if (plot) {
    # Create a data frame for plotting
    plot_data <- data.frame(
      time = as.numeric(rownames(est.catt)),
      catt = est.catt[, "ATT"],
      ci_lower = est.catt[, "CI.lower"],
      ci_upper = est.catt[, "CI.upper"]
    )

    # Add count data for the bar chart at the bottom
    time_range <- unique(plot_data$time)
    # Calculate treated units at each relative time point
    # First re-create the relative time matrix (similar to what getEffect does)
    D_rel <- D.old  # Start with original treatment matrix
    for (i in 1:ncol(D_rel)) {
      cumD <- cumsum(D_rel[, i])
      t0 <- sum(cumD == 0)  # Find time of first treatment
      if (t0 > 0) {  # Only process treated units
        D_rel[, i] <- 1:nrow(D_rel) - t0  # Convert to relative time
      } else {
        D_rel[, i] <- NA  # Not treated
      }
    }

    # Count treated units at each time point
    count_data <- data.frame(
      time = time_range,
      count = sapply(time_range, function(t) {
        sum(D_rel == t, na.rm = TRUE)
      })
    )

    # Calculate rectangle dimensions for count display with more space between bars and main plot
    y_range <- diff(range(plot_data$ci_lower, plot_data$ci_upper, na.rm = TRUE))
    rect.min <- min(plot_data$ci_lower, na.rm = TRUE) - 0.3 * y_range
    rect.length <- 0.15 * y_range

    T.start <- time_range - 0.25
    T.end <- time_range + 0.25
    ymin <- rep(rect.min, length(time_range))
    ymax <- rect.min + rect.length * count_data$count / max(count_data$count)

    data.toplot <- data.frame(
      xmin = T.start,
      xmax = T.end,
      ymin = ymin,
      ymax = ymax
    )

    # xlab, ylab, and main
    if (is.null(xlab) == TRUE) {xlab <- "Time"}
    if (is.null(ylab) == TRUE) {
      if (cumu == TRUE) {
        ylab <- "Cumulative ATT"
      } else {
        ylab <- "ATT"
      }
    }
    if (is.null(main) == TRUE) {
      if (cumu == TRUE) {
        main <- "Average Cumulative Treatment Effects on the Treated"
      } else {
        main <- "Average Treatment Effects on the Treated"
      }
    }

    # Create the plot
    p <- ggplot() +
      geom_point(data = plot_data, aes(x = .data$time, y = .data$catt), size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_linerange(data = plot_data,
                     aes(x = .data$time, ymin = .data$ci_lower, ymax = .data$ci_upper),
                     size = 0.5) +
      labs(x = xlab, y = ylab, title = main) +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 11),
        axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
        axis.text = element_text(color = "black", size = 10),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold", margin = margin(10, 0, 10, 0))
      ) +
      # Force integer breaks on x-axis
      scale_x_continuous(breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 1))

    if (count == TRUE) {
      p <- p + geom_rect(data = data.toplot,
                         aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
                         fill = "gray50", alpha = 0.3, size = 0.3, color = "black") +
        annotate("text",
                 x = time_range[which.max(count_data$count)],
                 y = max(data.toplot$ymax) + 0.05 * y_range,
                 label = max(count_data$count),
                 size = 3, hjust = 0.5)
    }
    # Print the plot
    print(p)
  }
  # Prepare output
  out = x
  out$effect.est.avg <- catt
  if (!is.null(catt.boot)) {
    # Add estimation results to output (commented out line would also include bootstrap samples)
    out$effect.est.att = est.catt
  }
  class(out) <- "fect"
  return(out)
}


getEffect <- function(D,           # Treatment indicator matrix
                      I,           # Inclusion indicator matrix
                      eff,         # Effect matrix
                      cumu,        # Logical: whether to calculate cumulative effect
                      period) {    # Event window range: c(start, end)
  # Initialize output vector with NAs
  aeff <- rep(NA, period[2] - period[1] + 1)

  # Return empty result if all D values are NA
  if (sum(is.na(D)) == length(c(D))) {
    return(aeff)
  }

  # Replace NA values in D with 0
  if (sum(is.na(c(D))) > 0) {
    D[which(is.na(D))] <- 0
  }

  # Calculate cumulative sum of treatment for each unit
  D <- apply(D, 2, function(vec) {
    cumsum(vec)
  })

  TT <- dim(D)[1]  # Number of time periods
  N <- dim(D)[2]   # Number of units

  # Calculate relative time to treatment for each unit
  for (i in 1:N) {
    subd <- D[, i]
    t0 <- sum(subd == 0)  # Pre-treatment periods
    D[, i] <- 1:TT - t0   # Convert to relative time
  }

  # Set treatment indicators to NA where units are not included
  if (sum(c(I) == 0) > 0) {
    D[which(I == 0)] <- NA
  }

  # Flatten matrices to vectors for processing
  vd <- c(D)       # Vector of relative times
  veff <- c(eff)   # Vector of effects

  # Remove NA entries
  if (sum(is.na(vd)) > 0) {
    vd.rm <- which(is.na(vd))
    vd <- vd[-vd.rm]
    veff <- veff[-vd.rm]
  }

  # Get unique relative time periods
  uniT <- sort(unique(vd))

  # Define effective time range
  ts <- period[1]  # Start period
  te <- min(period[2], max(vd))  # End period (capped by available data)
  effT <- ts:te    # Effective time range

  if (cumu == TRUE) {
    # Calculate cumulative treatment effect
    if (sum(!(effT %in% uniT)) == 0) {
      pos <- c()
      for (i in 1:length(effT)) {
        # Accumulate positions for all periods up to current one
        pos <- c(pos, which(vd == effT[i]))
        # Calculate cumulative effect as mean effect times number of periods
        aeff[i] <- mean(veff[pos]) * i
      }
    }
  } else {
    # Calculate average treatment effect for each period
    ave <- as.numeric(tapply(veff, vd, mean))
    effT2 <- period[1]:period[2]
    for (i in 1:length(effT2)) {
      if (effT2[i] %in% uniT) {
        aeff[i] <- ave[which(uniT == effT2[i])]
      }
    }
  }
  return(aeff)  # Return vector of treatment effects
}
