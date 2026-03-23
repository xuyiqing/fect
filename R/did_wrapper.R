# Treated-unit counts by event period
get_event_counts <- function(data, time_to_treat_var, unit_id) {
  time_sym <- rlang::sym(time_to_treat_var)
  unit_sym <- rlang::sym(unit_id)

  data %>%
    dplyr::filter(!is.na(!!time_sym)) %>%
    dplyr::group_by(!!time_sym, .drop = FALSE) %>%
    dplyr::summarise(
      count = dplyr::n_distinct(!!unit_sym),
      .groups = "drop"
    ) %>%
    dplyr::rename(Period = !!time_sym)
}


# Build event-study data frame
make_es_df <- function(Period_vec, est_vec, se_vec = NULL,
                       count_df = NULL, method = method) {
  df <- data.frame(ATT = est_vec, check.names = FALSE)

  if (!is.null(se_vec)) {
    df[["S.E."]] <- se_vec
    df[["CI.lower"]] <- df$ATT - 1.96 * se_vec
    df[["CI.upper"]] <- df$ATT + 1.96 * se_vec
    df[["p.value"]] <- 2 * (1 - stats::pnorm(abs(df$ATT / se_vec)))
  } else {
    df[["S.E."]] <- NA
    df[["CI.lower"]] <- NA
    df[["CI.upper"]] <- NA
    df[["p.value"]] <- NA
  }

  ## attach period counts if supplied
  if (!is.null(count_df)) {
    df$Period <- Period_vec + as.numeric(method != "didm")
    df <- df |>
      dplyr::left_join(
        count_df |>
          dplyr::mutate(Period = .data$Period + 1), # Assuming Period in count_df is 0-indexed relative to event
        by = "Period"
      )
    # If Period_vec can have NAs, this might cause issues with rownames
    # Ensure Period_vec used for rownames is clean or handle NAs appropriately
    valid_periods <- !is.na(df$Period) # Or use original Period_vec if df$Period might be modified by join
    if (any(valid_periods)) {
      rownames(df)[valid_periods] <- df$Period[valid_periods]
    }
    df$Period <- NULL
  } else {
    # Ensure Period_vec used for rownames is clean
    valid_rownames <- Period_vec + as.numeric(method != "didm")
    if (length(valid_rownames) == nrow(df) && !any(is.na(valid_rownames))) {
      rownames(df) <- valid_rownames
    }
    df$count <- NA
  }

  df <- df[, c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value", "count")]
  df
}

# ── Main wrapper ───────────────────────────────────────────────────────────

did_wrapper <- function(
    data, Y, D, X = NULL,
    index, # c("unit_id","time_id")
    method = c("twfe", "st", "iw", "cs_never", "cs_notyet", "didm"),
    se = c("default", "boot", "bootstrap", "jackknife"),
    nboots = 200,
    parallel = TRUE,
    core = NULL,
    ## TWFE / stacked-DID options
    time_to_treat_var = "Time_to_Treatment",
    treat_indicator = "treat",
    ## Callaway-Sant'Anna options
    csdid.base_period = "universal",
    ## DIDmultiplegtDYN options
    didm.effects = NA,
    didm.placebo = NA) {
  if (method == "bootstrap") method <- "boot"
  method <- match.arg(method)
  se <- match.arg(se)
  unit_id <- index[1]
  time_id <- index[2]

  # Define symbols for unit_id and D here for use throughout the function, including run_estimator_once
  unit_sym <- rlang::sym(unit_id)
  D_sym <- rlang::sym(D)
  # Y_sym    <- rlang::sym(Y) # If Y is used with dplyr NSE

  if (method == "didm" && is.na(didm.effects)) {
    stop('For "didm", you must supply "didm.effects" and "didm.placebo".')
  }

  ## ── Pre-processing ──────────────────────────────────────────────────────
  original_units <- length(unique(data[[unit_id]])) # Base R subsetting, fine

  data <- data %>%
    dplyr::group_by(!!unit_sym) %>% # Uses pre-defined unit_sym
    dplyr::mutate(treatment_mean = mean(!!D_sym, na.rm = TRUE)) %>% # Uses pre-defined D_sym
    dplyr::ungroup()
  data <- data[data$treatment_mean < 1, ] # Base R subsetting, fine
  message(
    "Dropped ", original_units - length(unique(data[[unit_id]])), # Base R subsetting, fine
    " always-treated units."
  )

  ## define cohorts & event time
  # Ensure fect::get.cohort doesn't have issues with column names if they contain special characters
  data <- fect::get.cohort(data, D = D, index = c(unit_id, time_id), start0 = TRUE)
  count_df <- get_event_counts(data, time_to_treat_var, unit_id)

  ## ── Core estimation routine (one sample) ───────────────────────────────
  run_estimator_once <- function(dd) { # dd is the data for this run

    out_ATT <- out_SE <- out_lower <- out_upper <- NA
    es_df <- data.frame()

    # ── 1) TWFE ───────────────────────────────────────────────────────────
    if (method == "twfe") {
      if (!treat_indicator %in% names(dd)) {
        dd <- dd |>
          dplyr::group_by(!!unit_sym) |> # MODIFIED: Use unit_sym from parent scope
          dplyr::mutate(tmp = as.numeric(mean(!!D_sym, na.rm = TRUE) > 0)) |> # MODIFIED: Use D_sym from parent scope
          dplyr::ungroup()
        names(dd)[names(dd) == "tmp"] <- treat_indicator
      }

      X_part <- if (!is.null(X)) paste(X, collapse = " + ") else "1"

      fit_att <- fixest::feols(
        stats::as.formula(
          paste0(Y, " ~ ", D, " + ", X_part, " | ", unit_id, " + ", time_id)
        ),
        data = dd,
        cluster = unit_id, # fixest can take string for cluster
        notes = FALSE
      )
      co <- fixest::coeftable(fit_att)
      est <- co[D, "Estimate"]
      se_ <- co[D, "Std. Error"]
      out_ATT <- est
      out_SE <- se_
      out_lower <- est - 1.96 * se_
      out_upper <- est + 1.96 * se_

      if (!time_to_treat_var %in% names(dd)) {
        dd <- fect::get.cohort(dd, D = D, index = c(unit_id, time_id), start0 = TRUE)
      }

      # Ensure time_to_treat_var column exists and handle NAs before using in i()
      if (time_to_treat_var %in% names(dd)) {
        dd[[time_to_treat_var]][is.na(dd[[time_to_treat_var]])] <- 0 # Or other suitable placeholder for fixest::i()
      } else {
        stop(paste("Column", time_to_treat_var, "not found in data for TWFE event study."))
      }
      # Ensure treat_indicator column exists
      if (!treat_indicator %in% names(dd)) {
        stop(paste("Column", treat_indicator, "not found in data for TWFE event study."))
      }


      fit_es <- fixest::feols(
        stats::as.formula(
          paste0(
            Y, " ~ i(", time_to_treat_var, ", ", treat_indicator, ", ref=-1) | ",
            unit_id, " + ", time_id
          )
        ),
        data = dd,
        cluster = unit_id, # fixest can take string for cluster
        notes = FALSE
      )

      ctab <- as.data.frame(fit_es$coeftable)
      # Robustly extract period numbers, handle cases where parsing fails
      period_nums <- vapply(
        rownames(ctab),
        function(x) {
          sp <- strsplit(x, "::")[[1]]
          if (length(sp) < 2) {
            return(NA_real_)
          }
          val_str <- strsplit(sp[2], ":")[[1]][1]
          suppressWarnings(as.numeric(val_str)) # Suppress warnings for non-numeric if any
        },
        numeric(1)
      )
      # Filter out NAs that might result from parsing non-coefficient rows or failed parsing
      valid_coeffs <- !is.na(period_nums)
      es_df <- make_es_df(
        period_nums[valid_coeffs], ctab$Estimate[valid_coeffs],
        ctab$`Std. Error`[valid_coeffs], count_df, method
      )


      # ── 2) Stacked DID ("st") ─────────────────────────────────────────────
    } else if (method == "st") {
      if (!"Cohort" %in% names(dd)) {
        dd <- fect::get.cohort(dd, D = D, index = c(unit_id, time_id), start0 = TRUE)
      }

      target.cohorts <- setdiff(unique(dd$Cohort), "Control")
      df.st_list <- vector("list", length(target.cohorts)) # Pre-allocate list
      for (i in seq_along(target.cohorts)) {
        c <- target.cohorts[i]
        # Use a temporary data frame to avoid modifying `dd` in a loop if it's large
        temp_df <- dd[dd$Cohort %in% c(c, "Control"), ]
        temp_df$stack_id <- i # Use loop index for stack_id
        df.st_list[[i]] <- temp_df
      }
      df.st <- do.call(rbind, df.st_list)


      df.st$st_unit <- as.numeric(factor(paste0(df.st$stack_id, "-", df.st[[unit_id]])))
      df.st$st_time <- as.numeric(factor(paste0(df.st$stack_id, "-", df.st[[time_id]])))

      fit_st <- fixest::feols(
        stats::as.formula(paste0(Y, " ~ ", D, " | st_unit + st_time")),
        data = df.st,
        cluster = "st_unit", # fixest can take string for cluster
        notes = FALSE,
        warn = FALSE
      )
      co <- fixest::coeftable(fit_st)
      est <- co[D, "Estimate"]
      se_ <- co[D, "Std. Error"]
      out_ATT <- est
      out_SE <- se_
      out_lower <- est - 1.96 * se_
      out_upper <- est + 1.96 * se_

      ## event study (same routine as above)
      if (!treat_indicator %in% names(df.st)) {
        df.st <- df.st |>
          dplyr::group_by(!!unit_sym) |> # MODIFIED: Use unit_sym from parent scope
          dplyr::mutate(tmp = as.numeric(mean(!!D_sym, na.rm = TRUE) > 0)) |> # MODIFIED: Use D_sym from parent scope
          dplyr::ungroup()
        names(df.st)[names(df.st) == "tmp"] <- treat_indicator
      }

      if (!time_to_treat_var %in% names(df.st)) {
        df.st <- fect::get.cohort(df.st, D = D, index = c(unit_id, time_id), start0 = TRUE)
      }

      # Ensure time_to_treat_var column exists and handle NAs
      if (time_to_treat_var %in% names(df.st)) {
        df.st[[time_to_treat_var]][is.na(df.st[[time_to_treat_var]])] <- 999999 # Or other placeholder
      } else {
        stop(paste("Column", time_to_treat_var, "not found in data for Stacked DID event study."))
      }
      if (!treat_indicator %in% names(df.st)) {
        stop(paste("Column", treat_indicator, "not found in data for Stacked DID event study."))
      }


      fit_es <- fixest::feols(
        stats::as.formula(
          paste0(
            Y, " ~ i(", time_to_treat_var, ", ", treat_indicator,
            ", ref=-1) | st_unit + st_time"
          )
        ),
        data = df.st,
        cluster = "st_unit" # fixest can take string for cluster
      )
      ctab <- as.data.frame(fit_es$coeftable)
      per <- vapply(
        rownames(ctab),
        function(x) {
          sp <- strsplit(x, "::")[[1]]
          if (length(sp) < 2) {
            return(NA_real_)
          }
          val_str <- strsplit(sp[2], ":")[[1]][1]
          suppressWarnings(as.numeric(val_str))
        },
        numeric(1)
      )
      valid_coeffs <- !is.na(per)
      es_df <- make_es_df(
        per[valid_coeffs], ctab$Estimate[valid_coeffs],
        ctab$`Std. Error`[valid_coeffs], count_df, method
      )


      # ── 3) IW (Sun-&-Abraham)  ────────────────────────────────────────────
    } else if (method == "iw") {
      if (!"FirstTreat" %in% names(dd)) { # 'FirstTreat' is typically generated by get.cohort or similar
        stop("Method 'iw' requires a FirstTreat column. Ensure get.cohort() was run or the column exists.")
      }

      dd$FirstTreat[is.na(dd$FirstTreat)] <- 1000 # Or other large value indicating never/not-yet treated for sunab
      X_part <- if (!is.null(X)) paste(X, collapse = " + ") else "1"

      fit_iw <- fixest::feols(
        stats::as.formula(
          # Ensure time_id is a valid column name for sunab
          paste0(
            Y, " ~ sunab(FirstTreat,", time_id, ") + ", X_part,
            " | ", unit_id, " + ", time_id
          )
        ),
        data = dd,
        cluster = unit_id # fixest can take string for cluster
      )
      att_sum <- summary(fit_iw, agg = "ATT")
      out_ATT <- att_sum$coeftable["ATT", "Estimate"]
      out_SE <- att_sum$coeftable["ATT", "Std. Error"]
      out_lower <- out_ATT - 1.96 * out_SE
      out_upper <- out_ATT + 1.96 * out_SE

      ctab <- as.data.frame(fixest::coeftable(fit_iw))
      # Robustly extract offsets for sunab
      offsets <- suppressWarnings(
        as.numeric(sub("^.*::", "", sub(":cohort::.*", "", rownames(ctab))))
      )
      valid <- !is.na(offsets)
      es_df <- make_es_df(
        offsets[valid], ctab$Estimate[valid],
        ctab$`Std. Error`[valid], count_df, method
      )


      # ── 4) Callaway-Sant'Anna (never-treated) ─────────────────────────────
    } else if (method == "cs_never") {
      if (!"FirstTreat" %in% names(dd)) { # 'FirstTreat' is gname for did::att_gt
        dd <- fect::get.cohort(dd, D = D, index = c(unit_id, time_id), start0 = TRUE)
      }

      # For 'nevertreated', control units usually have FirstTreat = 0 or Inf
      # Check did package documentation for exact requirement for gname with nevertreated
      dd$FirstTreat[is.na(dd$FirstTreat)] <- 0 # Assuming NA means never treated and coded as 0

      cs.out <- did::att_gt(
        yname = Y, gname = "FirstTreat", # Ensure 'FirstTreat' is the correct group identifier
        idname = unit_id, tname = time_id,
        xformla = if (!is.null(X)) stats::as.formula(paste("~", paste(X, collapse = "+"))) else ~1,
        control_group = "nevertreated",
        allow_unbalanced_panel = TRUE, # Consider implications
        data = as.data.frame(dd), # did package often prefers data.frame
        est_method = "reg", # Default, consider "dr" for doubly robust
        base_period = csdid.base_period
      )

      simple <- tryCatch(did::aggte(cs.out, type = "simple", na.rm = TRUE),
        error = function(e) {
          warning("Error in did::aggte (simple, cs_never): ", e$message)
          NULL
        }
      )
      if (!is.null(simple)) {
        out_ATT <- simple$overall.att
        out_SE <- simple$overall.se
        out_lower <- out_ATT - 1.96 * out_SE
        out_upper <- out_ATT + 1.96 * out_SE
      }

      dyn <- tryCatch(
        did::aggte(cs.out, type = "dynamic", na.rm = TRUE, cband = FALSE), # cband=FALSE for point CIs
        error = function(e) {
          warning("Error in did::aggte (dynamic, cs_never): ", e$message)
          NULL
        }
      )
      es_df <- if (is.null(dyn) || length(dyn$egt) == 0) { # Check if dyn$egt is empty
        make_es_df(numeric(0), numeric(0), numeric(0), count_df, method)
      } else {
        make_es_df(dyn$egt, dyn$att.egt, dyn$se.egt, count_df, method)
      }


      # ── 5) Callaway-Sant'Anna (not-yet-treated) ───────────────────────────
    } else if (method == "cs_notyet") {
      if (!"FirstTreat" %in% names(dd)) {
        dd <- fect::get.cohort(dd, D = D, index = c(unit_id, time_id), start0 = TRUE)
      }

      # For 'notyettreated', FirstTreat should be year of treatment or Inf/0 for never treated
      # Ensure NA handling for FirstTreat is appropriate for did::att_gt
      # dd$FirstTreat[is.na(dd$FirstTreat)] <- Inf # Or appropriate coding for never treated
      
      dd$FirstTreat[is.na(dd$FirstTreat)] <- 0

      cs.out <- did::att_gt(
        yname = Y, gname = "FirstTreat",
        idname = unit_id, tname = time_id,
        xformla = if (!is.null(X)) stats::as.formula(paste("~", paste(X, collapse = "+"))) else ~1,
        control_group = "notyettreated",
        allow_unbalanced_panel = TRUE,
        data = as.data.frame(dd), # did package often prefers data.frame
        est_method = "reg",
        base_period = csdid.base_period
      )

      simple <- tryCatch(did::aggte(cs.out, type = "simple", na.rm = TRUE),
        error = function(e) {
          warning("Error in did::aggte (simple, cs_notyet): ", e$message)
          NULL
        }
      )
      if (!is.null(simple)) {
        out_ATT <- simple$overall.att
        out_SE <- simple$overall.se
        out_lower <- out_ATT - 1.96 * out_SE
        out_upper <- out_ATT + 1.96 * out_SE
      }

      dyn <- tryCatch(
        did::aggte(cs.out, type = "dynamic", na.rm = TRUE, cband = FALSE),
        error = function(e) {
          warning("Error in did::aggte (dynamic, cs_notyet): ", e$message)
          NULL
        }
      )
      es_df <- if (is.null(dyn) || length(dyn$egt) == 0) {
        make_es_df(numeric(0), numeric(0), numeric(0), count_df, method)
      } else {
        make_es_df(dyn$egt, dyn$att.egt, dyn$se.egt, count_df, method)
      }


      # ── 6) DIDmultiplegtDYN Dyn ──────────────────────────────────────────────
    } else if (method == "didm") {
      res <- DIDmultiplegtDYN::did_multiplegt_dyn(
        df        = as.data.frame(dd), # DIDmultiplegtDYN prefers data.frame
        outcome   = Y,
        group     = unit_id,
        time      = time_id,
        treatment = D,
        controls  = X, # Pass X directly if it's a character vector of control names
        effects   = didm.effects,
        placebo   = didm.placebo,
        cluster   = unit_id, # Ensure this column exists and is appropriate for clustering
        graph_off = TRUE
      )

      # Check structure of res$results for safety
      if (!is.null(res$results) && "ATE" %in% names(res$results) && length(res$results$ATE) >= 2) {
        out_ATT <- res$results$ATE[1]
        out_SE <- res$results$ATE[2]
        if (length(res$results$ATE) >= 4) { # Check if CI is provided
          out_lower <- res$results$ATE[3]
          out_upper <- res$results$ATE[4]
        } else {
          out_lower <- out_ATT - 1.96 * out_SE
          out_upper <- out_ATT + 1.96 * out_SE
        }
      } else {
        warning("Results from did_multiplegt_dyn ATE not found or in unexpected format.")
      }


      Placebos <- res$results$Placebos
      Effects <- res$results$Effects
      T.pre <- if (!is.null(Placebos) && inherits(Placebos, "matrix")) nrow(Placebos) else 0
      T.post <- if (!is.null(Effects) && inherits(Effects, "matrix")) nrow(Effects) else 0


      if (T.pre + T.post == 0) {
        es_df <- make_es_df(numeric(0), numeric(0), numeric(0), count_df, method)
      } else {
        est_placebo <- if (T.pre > 0 && "Estimate" %in% colnames(Placebos)) Placebos[, "Estimate"] else numeric(0)
        est_effect <- if (T.post > 0 && "Estimate" %in% colnames(Effects)) Effects[, "Estimate"] else numeric(0)
        est <- c(est_placebo, est_effect)

        se_placebo <- if (T.pre > 0 && "SE" %in% colnames(Placebos)) Placebos[, "SE"] else numeric(0)
        se_effect <- if (T.post > 0 && "SE" %in% colnames(Effects)) Effects[, "SE"] else numeric(0)
        ses <- c(se_placebo, se_effect)

        # Ensure lengths match before creating data frame
        if (length(est) != length(ses)) {
          warning("Mismatch in length of estimates and standard errors from did_multiplegt_dyn.")
          # Fallback to empty or handle error appropriately
          es_df <- make_es_df(numeric(0), numeric(0), numeric(0), count_df, method)
        } else {
          peri <- c(
            if (T.pre > 0) seq(-T.pre, -1) else numeric(0),
            if (T.post > 0) seq(0, T.post - 1) else numeric(0)
          ) # didm effects often 0-indexed

          peri <- c(
            if (T.pre > 0) seq(-T.pre, -1) else numeric(0),
            if (T.post > 0) seq(1, T.post) else numeric(0)
          )

          es_df <- make_es_df(peri, est, ses, count_df, method)
        }
      }
    }

    list(
      ATT = out_ATT, ATT_se = out_SE,
      CI_lower = out_lower, CI_upper = out_upper,
      es = es_df
    )
  }

  ## ── SE / bootstrap handling ────────────────────────────────────────────
  if (se == "default") {
    out_full <- run_estimator_once(data)
    # Ensure ATT_se is not zero or NA before calculating p-value
    overall_p <- if (!is.na(out_full$ATT_se) && out_full$ATT_se != 0) {
      2 * (1 - stats::pnorm(abs(out_full$ATT / out_full$ATT_se)))
    } else {
      NA_real_
    }
  } else { # bootstrap or jackknife

    cluster_ids <- unique(data[[unit_id]]) # Base R, fine
    n_clusters <- length(cluster_ids)

    # Ensure nboots is sensible for jackknife
    if (se == "jackknife" && nboots != n_clusters) {
      message(paste0("For jackknife, nboots should be equal to the number of clusters (", n_clusters, "). Setting nboots to ", n_clusters, "."))
      nboots <- n_clusters
    }


    boot_fun <- function(b) { # b is the iteration number
      if (se %in% c("boot", "bootstrap")) {
        # Stratified bootstrap if applicable, or simple cluster bootstrap
        samp_indices <- sample(seq_along(cluster_ids), n_clusters, replace = TRUE)
        samp_cluster_names <- cluster_ids[samp_indices]
        # Efficiently create bootstrap sample by selecting rows
        # This assumes unit_id column is a factor or character
        d_b <- data[data[[unit_id]] %in% samp_cluster_names, ]
        # To handle duplicated clusters correctly if units within clusters are distinct in d_b:
        # This requires more care if unit IDs are not unique across original clusters
        # A common way is to replicate data for chosen clusters and assign new unique IDs if needed by downstream.
        # However, fixest handles repeated clusters in cluster argument if data is just rbind-ed.
        # The current approach data[data[[unit_id]] %in% samp, ] is standard for cluster bootstrap.
      } else { # jackknife
        omit_cluster_name <- cluster_ids[b] # b is 1 to n_clusters
        d_b <- data[data[[unit_id]] != omit_cluster_name, ]
      }
      run_estimator_once(d_b)
    }

    if (parallel) {
      if (is.null(core)) core <- parallelly::availableCores(omit = 1)
      future::plan(future::multisession, workers = core) # Use future::multisession
      # Ensure seed is handled correctly for parallel processing
      rep_list <- future.apply::future_lapply(seq_len(nboots), boot_fun,
        future.seed = TRUE
      ) # future.seed=TRUE is good
      future::plan(future::sequential) # Reset plan
    } else {
      rep_list <- lapply(seq_len(nboots), boot_fun)
    }

    att_reps <- vapply(rep_list, `[[`, numeric(1), "ATT")
    # Remove NAs from bootstrap replications if any estimator failed
    att_reps_clean <- att_reps[!is.na(att_reps)]
    if (length(att_reps_clean) < length(att_reps)) {
      warning(paste(length(att_reps) - length(att_reps_clean), "bootstrap/jackknife replications resulted in NA ATT values and were dropped."))
    }
    if (length(att_reps_clean) < 2) { # Need at least 2 for sd()
      warning("Too few successful bootstrap/jackknife replications to calculate SE/CI.")
      out_full <- run_estimator_once(data) # Get point estimate
      out_full$ATT_se <- NA_real_
      out_full$CI_lower <- NA_real_
      out_full$CI_upper <- NA_real_
      overall_p <- NA_real_
    } else {
      out_full <- run_estimator_once(data) # Get point estimate from full sample
      theta_hat <- out_full$ATT
      out_full$CI_lower <- 2 * theta_hat - stats::quantile(att_reps_clean, 0.975, na.rm = TRUE)
      out_full$CI_upper <- 2 * theta_hat - stats::quantile(att_reps_clean, 0.025, na.rm = TRUE)
      out_full$ATT_se <- stats::sd(att_reps_clean, na.rm = TRUE)
      p1 <- if (theta_hat > 0) {
        mean(att_reps_clean <= 0, na.rm = TRUE)
      } else {
        mean(att_reps_clean >= 0, na.rm = TRUE)
      }
      overall_p <- min(2 * p1, 1)
    }

    if (nrow(out_full$es) > 0) {
      periods_from_es <- as.numeric(rownames(out_full$es))


      boot_mat_list <- lapply(rep_list, function(x) {
        if (is.null(x$es) || nrow(x$es) == 0) {
          return(rep(NA_real_, length(periods_from_es)))
        }
        # Match based on rownames (periods)
        es_rep_periods <- as.numeric(rownames(x$es))
        matched_att <- x$es$ATT[match(periods_from_es, es_rep_periods)]
        # If some periods are not in x$es, matched_att will have NAs, which is fine.
        return(matched_att)
      })
      boot_mat <- do.call(rbind, boot_mat_list) # Rows are bootstrap reps, columns are periods

      # Filter out columns (periods) that are all NA if any
      valid_cols <- apply(boot_mat, 2, function(col) !all(is.na(col)))
      boot_mat_filtered <- boot_mat[, valid_cols, drop = FALSE]
      periods_filtered <- periods_from_es[valid_cols]
      original_att_filtered <- out_full$es$ATT[valid_cols]

      if (ncol(boot_mat_filtered) > 0) {
        boot_se <- apply(boot_mat_filtered, 2, stats::sd, na.rm = TRUE)
        ci_quantiles_upper <- apply(boot_mat_filtered, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
        ci_quantiles_lower <- apply(boot_mat_filtered, 2, stats::quantile, probs = 0.025, na.rm = TRUE)

        pvals_es <- vapply(
          seq_len(ncol(boot_mat_filtered)),
          function(j) {
            est_j <- original_att_filtered[j]
            reps_j <- boot_mat_filtered[, j]
            reps_j_clean <- reps_j[!is.na(reps_j)]
            if (length(reps_j_clean) < 2) {
              return(NA_real_)
            }
            p1 <- if (est_j > 0) {
              mean(reps_j_clean <= 0, na.rm = TRUE)
            } else {
              mean(reps_j_clean >= 0, na.rm = TRUE)
            }
            min(2 * p1, 1)
          },
          numeric(1)
        )

        # Assign back to the original out_full$es structure, potentially for a subset of periods
        out_full$es$`S.E.` <- NA_real_
        out_full$es$CI.lower <- NA_real_
        out_full$es$CI.upper <- NA_real_
        out_full$es$p.value <- NA_real_

        out_full$es$`S.E.`[valid_cols] <- boot_se
        out_full$es$CI.lower[valid_cols] <- 2 * original_att_filtered - ci_quantiles_upper
        out_full$es$CI.upper[valid_cols] <- 2 * original_att_filtered - ci_quantiles_lower
        out_full$es$p.value[valid_cols] <- pvals_es
      } else {
        warning("No valid data for bootstrapping event study standard errors.")
      }
    }
  }

  ## ── Return object ───────────────────────────────────────────────────────
  res <- list(
    est.avg = data.frame(
      ATT.avg = out_full$ATT,
      S.E. = out_full$ATT_se,
      CI.lower = out_full$CI_lower,
      CI.upper = out_full$CI_upper,
      p.value = overall_p,
      row.names = NULL # Ensure no rownames for this single-row data frame
    ),
    est.att = out_full$es # This should have periods as rownames from make_es_df
  )
  class(res) <- "did_wrapper"
  res
}
