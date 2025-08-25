.has_honest <- function() {
  requireNamespace("HonestDiDFEct", quietly = TRUE)
}

.honest <- function(fun) {
  if (!.has_honest()) {
    stop("Optional dependency ``HonestDiDFEct'' not available. Install from https://github.com/lzy318/HonestDiDFEct .", call. = FALSE)
  }
  getExportedValue("HonestDiDFEct", fun)
}


fect_sens <- function(
    fect.out,
    post.periods = NA, # Post-treatment periods
    l_vec = NA, # Optional custom weighting vector
    Mbarvec = seq(0, 1, by = 0.1), # Vector of Mbar for RM analysis
    Mvec = seq(0, 0.25, 0.05), # Vector of M for Smoothness analysis
    periodMbarvec = c(0, 0.5), # Vector of Mbar for period-by-period RM analysis
    periodMvec = c(0, 0.1), # Vector of M for Smoothness analysis
    parallel = FALSE,
    cores = NULL) {
  if (parallel && is.null(cores)) {
    cores <- min(parallel::detectCores() - 2, 8) # default to 8 cores if not specified
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }

  # -------------------------------------------------------------------
  # 1) Identify relevant periods and extract fect estimates + vcov
  # -------------------------------------------------------------------
  # We want event-time strings for pre and post
  pre.periods <- fect.out$placebo.period[1]:fect.out$placebo.period[2]
  if (any(is.na(post.periods))) {
    post.periods <- 1:length(fect.out$time) - length(fect.out$pre.periods) - (fect.out$placebo.period[2] - fect.out$placebo.period[1] + 1)
  }
  all.periods.char <- as.character(c(pre.periods, post.periods))

  # Make sure that:
  #   - rownames(fect.out$est.att) match these event times (e.g. "-2","-1","0","1", etc.)
  #   - fect.out$att.vcov has dimension names that match rownames(fect.out$est.att)
  # If 'att.vcov' lacks dimnames, you can do:
  #   dimnames(fect.out$att.vcov) <- list(rownames(fect.out$est.att),
  #                                       rownames(fect.out$est.att))

  # Numeric indices corresponding to those event-time strings
  idx <- match(all.periods.char, rownames(fect.out$est.att))
  idx <- idx[!is.na(idx)] # remove anything not found

  # Extract DTE estimates (beta.hat) and var-cov (vcov.hat)
  beta.hat <- fect.out$est.att[idx, 1]
  vcov.hat <- fect.out$att.vcov[idx, idx]

  # Counts of pre and post periods
  numPrePeriods <- length(pre.periods)
  numPostPeriods <- length(post.periods)

  # -------------------------------------------------------------------
  # 2) Construct weights for an "overall" post-treatment ATT
  # -------------------------------------------------------------------
  # If the user did NOT provide a custom weighting vector l_vec,
  # default to count-based weighting (counts from each post period).
  if (any(is.na(l_vec))) {
    # Indices for post-periods in the rownames
    post.idx <- match(as.character(post.periods), rownames(fect.out$est.att))
    post.idx <- post.idx[!is.na(post.idx)]

    count <- fect.out$count[post.idx]
    w.att <- count / sum(count)
  } else {
    w.att <- l_vec
  }

  # We will attach everything to the original fect.out in new sub-lists.
  # For each approach (RM or Smoothness), we will store both:
  #   (a) Weighted-average results
  #   (b) Period-by-period results

  # Initialize empty placeholders
  fect.out$RM_Sensitivity <- NULL
  fect.out$Smooth_Sensitivity <- NULL

  # -------------------------------------------------------------------
  # 3) Relative Magnitude Analysis (RM), if Mbarvec is non-empty
  # -------------------------------------------------------------------
  if (!is.null(Mbarvec) && length(Mbarvec) > 0) {
    # 3a) Weighted-average, across the entire post-treatment window
    rm_sens_results <- .honest("createSensitivityResults_relativeMagnitudes")(
      betahat = beta.hat,
      sigma = vcov.hat,
      numPrePeriods = numPrePeriods,
      numPostPeriods = numPostPeriods,
      l_vec = w.att,
      Mbarvec = Mbarvec,
      parallel = parallel
    )


    rm_original_cs <- .honest("constructOriginalCS")(
      betahat        = beta.hat,
      sigma          = vcov.hat,
      numPrePeriods  = numPrePeriods,
      numPostPeriods = numPostPeriods,
      l_vec          = w.att
    )
  }
  if (!is.null(periodMbarvec) && length(periodMbarvec) > 0) {
    # 3b) Period-by-period robust confidence sets
    #     We'll loop over each post-treatment period and set a vector of 0s except 1 for that period
    rm_period_output <- cbind.data.frame() # Will accumulate results

    for (t_i in seq_len(numPostPeriods)) {
      # l_vec for the "t_i-th" post period (in 1..numPostPeriods)
      # We need a vector of length numPrePeriods + numPostPeriods
      # Set that post period's position to 1
      dte_l <- rep(0, numPostPeriods)
      dte_l[t_i] <- 1

      # For each t_i, we run createSensitivityResults_relativeMagnitudes
      # across all Mbar in Mbarvec
      honest.dte <- .honest("createSensitivityResults_relativeMagnitudes")(
        betahat        = beta.hat,
        sigma          = vcov.hat,
        numPrePeriods  = numPrePeriods,
        numPostPeriods = numPostPeriods,
        l_vec          = dte_l,
        Mbarvec        = periodMbarvec,
        parallel       = parallel
      )

      # Convert to data.frame
      # The returned object typically has columns lb, ub, Mbar, etc.
      # We add a column for the actual "post period"
      honest.dte <- as.data.frame(honest.dte)
      # Actual event time is post.periods[t_i]
      honest.dte$postPeriod <- post.periods[t_i]

      # Append
      rm_period_output <- rbind(rm_period_output, honest.dte)
    }

    # Store results in fect.out$sensitivity.rm
    fect.out$sensitivity.rm <- list(
      results = rm_sens_results,
      original = rm_original_cs,
      periods = rm_period_output
    )
  }

  # -------------------------------------------------------------------
  # 4) Smoothness Analysis (C-LF), if Mvec is non-empty
  # -------------------------------------------------------------------
  if (!is.null(Mvec) && length(Mvec) > 0) {
    # 4a) Weighted-average analysis

    smooth_sens_results <- .honest("createSensitivityResults")(
      betahat = beta.hat,
      sigma = vcov.hat,
      numPrePeriods = numPrePeriods,
      numPostPeriods = numPostPeriods,
      method = "C-LF",
      l_vec = w.att,
      Mvec = Mvec,
      parallel = parallel
    )

    sm_original_cs <- .honest("constructOriginalCS")(
      betahat        = beta.hat,
      sigma          = vcov.hat,
      numPrePeriods  = numPrePeriods,
      numPostPeriods = numPostPeriods,
      l_vec          = w.att
    )
  }
  if (!is.null(periodMvec) && length(periodMvec) > 0) {
    # 4b) Period-by-period robust confidence sets
    smooth_period_output <- cbind.data.frame()

    for (t_i in seq_len(numPostPeriods)) {
      # l_vec for the "t_i-th" post period
      dte_l <- rep(0, numPostPeriods)
      dte_l[t_i] <- 1

      honest.dte <- .honest("createSensitivityResults")(
        betahat = beta.hat,
        sigma = vcov.hat,
        numPrePeriods = numPrePeriods,
        numPostPeriods = numPostPeriods,
        method = "C-LF",
        l_vec = dte_l,
        Mvec = periodMvec,
        parallel = parallel
      )

      honest.dte <- as.data.frame(honest.dte)
      honest.dte$postPeriod <- post.periods[t_i]

      smooth_period_output <- rbind(smooth_period_output, honest.dte)
    }

    # Store results in fect.out$sensitivity.smooth
    fect.out$sensitivity.smooth <- list(
      results = smooth_sens_results,
      original = sm_original_cs,
      periods = smooth_period_output
    )
  }
  if (parallel) {
    stopCluster(cl)
  }
  # -------------------------------------------------------------------
  # 5) Return the updated fect.out
  # -------------------------------------------------------------------
  return(fect.out)
}
