# Treated-unit counts by event period
get_event_counts <- function(data, time_to_treat_var, unit_id) {
  data |>
    dplyr::filter(!is.na(rlang::.data[[time_to_treat_var]])) |>
    dplyr::group_by(rlang::.data[[time_to_treat_var]], .drop = FALSE) |>
    dplyr::summarise(count = dplyr::n_distinct(rlang::.data[[unit_id]]),
                     .groups = "drop") |>
    dplyr::rename(Period = {{ time_to_treat_var }})
}

# Build event-study data frame
make_es_df <- function(Period_vec, est_vec, se_vec = NULL,
                       count_df = NULL, method = method) {

  df <- data.frame(ATT = est_vec, check.names = FALSE)

  if (!is.null(se_vec)) {
    df[["S.E."]]      <- se_vec
    df[["CI.lower"]]  <- df$ATT - 1.96 * se_vec
    df[["CI.upper"]]  <- df$ATT + 1.96 * se_vec
    df[["p.value"]]   <- 2 * (1 - stats::pnorm(abs(df$ATT / se_vec)))
  } else {
    df[["S.E."]]      <- NA
    df[["CI.lower"]]  <- NA
    df[["CI.upper"]]  <- NA
    df[["p.value"]]   <- NA
  }

  ## attach period counts if supplied
  if (!is.null(count_df)) {
    df$Period <- Period_vec + as.numeric(method != "didm")
    df <- df |>
      dplyr::left_join(
        count_df |>
          dplyr::mutate(Period = Period + 1),
        by = "Period"
      )
    rownames(df) <- df$Period
    df$Period    <- NULL
  } else {
    rownames(df) <- Period_vec + as.numeric(method != "didm")
    df$count     <- NA
  }

  df <- df[, c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value", "count")]
  df
}

# ── Main wrapper ───────────────────────────────────────────────────────────

did_wrapper <- function(
    data, Y, D, X = NULL,
    index,                 # c("unit_id","time_id")
    method = c("twfe","st","iw","cs_never","cs_notyet","pm","didm"),
    se     = c("default","boot","bootstrap","jackknife"),
    nboots = 200,
    parallel = TRUE,
    core     = NULL,

    ## TWFE / stacked-DID options
    time_to_treat_var = "Time_to_Treatment",
    treat_indicator   = "treat",

    ## Callaway–Sant’Anna options
    csdid.base_period = "universal",

    ## DIDmultiplegt options
    didm.effects = NA,
    didm.placebo = NA
) {

  if (method == "bootstrap") method <- "boot"
  method <- match.arg(method)
  se     <- match.arg(se)
  unit_id <- index[1];  time_id <- index[2]

  if (method == "didm" && is.na(didm.effects))
    stop('For "didm", you must supply "didm.effects" and "didm.placebo".')

  ## ── Pre-processing ──────────────────────────────────────────────────────
  original_units <- length(unique(data[[unit_id]]))

  data <- data |>
    dplyr::group_by(rlang::.data[[unit_id]]) |>
    dplyr::mutate(treatment_mean = mean(rlang::.data[[D]], na.rm = TRUE)) |>
    dplyr::ungroup()

  data <- data[data$treatment_mean < 1, ]
  message("Dropped ", original_units - length(unique(data[[unit_id]])),
          " always-treated units.")

  ## define cohorts & event time
  data <- fect::get.cohort(data, D = D, index = c(unit_id, time_id), start0 = TRUE)
  count_df <- get_event_counts(data, time_to_treat_var, unit_id)

  ## ── Core estimation routine (one sample) ───────────────────────────────
  run_estimator_once <- function(dd) {

    out_ATT   <- out_SE <- out_lower <- out_upper <- NA
    es_df     <- data.frame()

    # ── 1) TWFE ───────────────────────────────────────────────────────────
    if (method == "twfe") {

      if (!treat_indicator %in% names(dd)) {
        dd <- dd |>
          dplyr::group_by(rlang::.data[[unit_id]]) |>
          dplyr::mutate(tmp = as.numeric(mean(rlang::.data[[D]], na.rm = TRUE) > 0)) |>
          dplyr::ungroup()
        names(dd)[names(dd) == "tmp"] <- treat_indicator
      }

      X_part <- if (!is.null(X)) paste(X, collapse = " + ") else "1"

      fit_att <- fixest::feols(
        stats::as.formula(
          paste0(Y, " ~ ", D, " + ", X_part, " | ", unit_id, " + ", time_id)
        ),
        data    = dd,
        cluster = unit_id,
        notes   = FALSE
      )
      co        <- fixest::coeftable(fit_att)
      est       <- co[D, "Estimate"]
      se_       <- co[D, "Std. Error"]
      out_ATT   <- est
      out_SE    <- se_
      out_lower <- est - 1.96 * se_
      out_upper <- est + 1.96 * se_

      if (!time_to_treat_var %in% names(dd))
        dd <- fect::get.cohort(dd, D = D, index = c(unit_id, time_id), start0 = TRUE)

      dd[[time_to_treat_var]][is.na(dd[[time_to_treat_var]])] <- 0

      fit_es <- fixest::feols(
        stats::as.formula(
          paste0(
            Y, " ~ i(", time_to_treat_var, ", ", treat_indicator, ", ref=-1) | ",
            unit_id, " + ", time_id
          )
        ),
        data    = dd,
        cluster = unit_id,
        notes   = FALSE
      )

      ctab  <- as.data.frame(fit_es$coeftable)
      period_nums <- vapply(
        rownames(ctab),
        function(x) {
          sp <- strsplit(x, "::")[[1]]
          if (length(sp) < 2) return(NA_real_)
          as.numeric(strsplit(sp[2], ":")[[1]][1])
        },
        numeric(1)
      )
      es_df <- make_es_df(period_nums, ctab$Estimate, ctab$`Std. Error`,
                          count_df, method)

      # ── 2) Stacked DID ("st") ─────────────────────────────────────────────
    } else if (method == "st") {

      if (!"Cohort" %in% names(dd))
        dd <- fect::get.cohort(dd, D = D, index = c(unit_id, time_id), start0 = TRUE)

      target.cohorts <- setdiff(unique(dd$Cohort), "Control")
      df.st <- NULL; k <- 1
      for (c in target.cohorts) {
        df.st <- rbind(
          df.st,
          within(dd[dd$Cohort %in% c(c,"Control"), ], { stack_id <- k })
        )
        k <- k + 1
      }

      df.st$st_unit <- as.numeric(factor(paste0(df.st$stack_id,"-",df.st[[unit_id]])))
      df.st$st_time <- as.numeric(factor(paste0(df.st$stack_id,"-",df.st[[time_id]])))

      fit_st <- fixest::feols(
        stats::as.formula(paste0(Y," ~ ",D," | st_unit + st_time")),
        data    = df.st,
        cluster = "st_unit",
        notes = FALSE,
        warn = FALSE
      )
      co        <- fixest::coeftable(fit_st)
      est       <- co[D,"Estimate"];  se_ <- co[D,"Std. Error"]
      out_ATT   <- est;               out_SE <- se_
      out_lower <- est - 1.96 * se_;  out_upper <- est + 1.96 * se_

      ## event study (same routine as above)
      if (!treat_indicator %in% names(df.st)) {
        df.st <- df.st |>
          dplyr::group_by(rlang::.data[[unit_id]]) |>
          dplyr::mutate(tmp = as.numeric(mean(rlang::.data[[D]],na.rm=TRUE) > 0)) |>
          dplyr::ungroup()
        names(df.st)[names(df.st) == "tmp"] <- treat_indicator
      }

      if (!time_to_treat_var %in% names(df.st))
        df.st <- fect::get.cohort(df.st, D = D, index = c(unit_id,time_id), start0 = TRUE)

      df.st[[time_to_treat_var]][is.na(df.st[[time_to_treat_var]])] <- 999999

      fit_es <- fixest::feols(
        stats::as.formula(
          paste0(Y," ~ i(", time_to_treat_var, ", ", treat_indicator,
                 ", ref=-1) | st_unit + st_time")
        ),
        data    = df.st,
        cluster = "st_unit"
      )
      ctab <- as.data.frame(fit_es$coeftable)
      per  <- vapply(
        rownames(ctab),
        function(x){
          sp<-strsplit(x,"::")[[1]]
          if(length(sp)<2) return(NA_real_)
          as.numeric(strsplit(sp[2],":")[[1]][1])
        },
        numeric(1)
      )
      es_df <- make_es_df(per, ctab$Estimate, ctab$`Std. Error`,
                          count_df, method)

      # ── 3) IW (Sun-&-Abraham)  ────────────────────────────────────────────
    } else if (method == "iw") {

      if (!"FirstTreat" %in% names(dd))
        stop("Method 'iw' requires a FirstTreat column.")

      dd$FirstTreat[is.na(dd$FirstTreat)] <- 1000
      X_part <- if (!is.null(X)) paste(X, collapse = " + ") else "1"

      fit_iw <- fixest::feols(
        stats::as.formula(
          paste0(Y," ~ sunab(FirstTreat,",time_id,") + ", X_part,
                 " | ", unit_id," + ", time_id)
        ),
        data    = dd,
        cluster = unit_id
      )
      att_sum  <- summary(fit_iw, agg = "ATT")
      out_ATT  <- att_sum$coeftable["ATT","Estimate"]
      out_SE   <- att_sum$coeftable["ATT","Std. Error"]
      out_lower<- out_ATT - 1.96*out_SE
      out_upper<- out_ATT + 1.96*out_SE

      ctab <- as.data.frame(fixest::coeftable(fit_iw))
      offsets <- suppressWarnings(
        as.numeric(sub("^.*::","",sub(":cohort::.*","",rownames(ctab))))
      )
      valid   <- !is.na(offsets)
      es_df   <- make_es_df(offsets[valid], ctab$Estimate[valid],
                            ctab$`Std. Error`[valid], count_df, method)

      # ── 4) Callaway–Sant’Anna (never-treated) ─────────────────────────────
    } else if (method == "cs_never") {

      if (!"FirstTreat" %in% names(dd))
        dd <- fect::get.cohort(dd, D = D, index = c(unit_id,time_id), start0 = TRUE)

      dd$FirstTreat[is.na(dd$FirstTreat)] <- 0

      cs.out <- did::att_gt(
        yname   = Y, gname = "FirstTreat",
        idname  = unit_id, tname = time_id,
        xformla = ~1, control_group = "nevertreated",
        allow_unbalanced_panel = TRUE,
        data    = dd, est_method = "reg",
        base_period = csdid.base_period
      )

      simple <- tryCatch(did::aggte(cs.out, type = "simple", na.rm = TRUE),
                         error = function(e) NULL)
      if (!is.null(simple)) {
        out_ATT   <- simple$overall.att
        out_SE    <- simple$overall.se
        out_lower <- out_ATT - 1.96*out_SE
        out_upper <- out_ATT + 1.96*out_SE
      }

      dyn <- tryCatch(
        did::aggte(cs.out, type="dynamic", na.rm = TRUE, cband = FALSE),
        error = function(e) NULL)
      es_df <- if (is.null(dyn))
        make_es_df(numeric(0), numeric(0), numeric(0), count_df, method)
      else
        make_es_df(dyn$egt, dyn$att.egt, dyn$se.egt, count_df, method)

      # ── 5) Callaway–Sant’Anna (not-yet-treated) ───────────────────────────
    } else if (method == "cs_notyet") {

      if (!"FirstTreat" %in% names(dd))
        dd <- fect::get.cohort(dd, D = D, index = c(unit_id,time_id), start0 = TRUE)

      cs.out <- did::att_gt(
        yname   = Y, gname = "FirstTreat",
        idname  = unit_id, tname = time_id,
        xformla = ~1, control_group = "notyettreated",
        allow_unbalanced_panel = TRUE,
        data    = dd, est_method = "reg",
        base_period = csdid.base_period
      )

      simple <- tryCatch(did::aggte(cs.out, type="simple", na.rm=TRUE),
                         error=function(e) NULL)
      if (!is.null(simple)) {
        out_ATT <- simple$overall.att
        out_SE  <- simple$overall.se
        out_lower <- out_ATT - 1.96*out_SE
        out_upper <- out_ATT + 1.96*out_SE
      }

      dyn <- tryCatch(
        did::aggte(cs.out, type="dynamic", na.rm=TRUE, cband=FALSE),
        error=function(e) NULL)
      es_df <- if (is.null(dyn))
        make_es_df(numeric(0), numeric(0), numeric(0), count_df, method)
      else
        make_es_df(dyn$egt, dyn$att.egt, dyn$se.egt, count_df, method)

      # ── 6) DIDmultiplegt Dyn ──────────────────────────────────────────────
    } else if (method == "didm") {

      res <- didmultiplegt::did_multiplegt_dyn(
        df        = dd,
        outcome   = Y,
        group     = unit_id,
        time      = time_id,
        treatment = D,
        controls  = NULL,
        effects   = didm.effects,
        placebo   = didm.placebo,
        cluster   = unit_id,
        graph_off = TRUE
      )

      out_ATT <- res$results$ATE[1]
      out_SE  <- res$results$ATE[2]
      if (length(res$results$ATE) >= 4) {
        out_lower <- res$results$ATE[3]
        out_upper <- res$results$ATE[4]
      } else {
        out_lower <- out_ATT - 1.96*out_SE
        out_upper <- out_ATT + 1.96*out_SE
      }

      Placebos <- res$results$Placebos
      Effects  <- res$results$Effects
      T.pre  <- if (!is.null(Placebos)) nrow(Placebos) else 0
      T.post <- if (!is.null(Effects )) nrow(Effects ) else 0

      if (T.pre + T.post == 0) {
        es_df <- make_es_df(numeric(0), numeric(0), numeric(0), count_df, method)
      } else {
        est  <- c(if (T.pre>0) Placebos[,"Estimate"] else numeric(0),
                  if (T.post>0) Effects[,"Estimate"]  else numeric(0))
        ses  <- c(if (T.pre>0) Placebos[,"SE"]       else numeric(0),
                  if (T.post>0) Effects[,"SE"]        else numeric(0))
        peri <- c(if (T.pre>0) seq(-T.pre,-1) else numeric(0),
                  if (T.post>0) seq(1,T.post) else numeric(0))
        es_df <- make_es_df(peri, est, ses, count_df, method)
      }
    }

    list(ATT = out_ATT, ATT_se = out_SE,
         CI_lower = out_lower, CI_upper = out_upper,
         es = es_df)
  }

  ## ── SE / bootstrap handling ────────────────────────────────────────────
  if (se == "default") {

    out_full <- run_estimator_once(data)
    overall_p <- 2 * (1 - stats::pnorm(abs(out_full$ATT / out_full$ATT_se)))

  } else {

    cluster_ids <- unique(data[[unit_id]])
    n_clusters  <- length(cluster_ids)

    boot_fun <- function(b) {
      if (se %in% c("boot","bootstrap")) {
        samp <- sample(cluster_ids, n_clusters, replace = TRUE)
        d_b  <- data[data[[unit_id]] %in% samp, ]
      } else {                             # jackknife
        omit <- cluster_ids[b]
        d_b  <- data[data[[unit_id]] != omit, ]
      }
      run_estimator_once(d_b)
    }

    if (parallel) {
      if (is.null(core)) core <- future::availableCores() - 1
      future::plan(multisession, workers = core)
      rep_list <- future.apply::future_lapply(seq_len(nboots), boot_fun,
                                              future.seed = TRUE)
      future::plan(sequential)
    } else {
      rep_list <- lapply(seq_len(nboots), boot_fun)
    }

    att_reps  <- vapply(rep_list, `[[`, numeric(1), "ATT")
    out_full  <- run_estimator_once(data)
    theta_hat <- out_full$ATT

    out_full$CI_lower <- 2*theta_hat - stats::quantile(att_reps, 0.975, na.rm=TRUE)
    out_full$CI_upper <- 2*theta_hat - stats::quantile(att_reps, 0.025, na.rm=TRUE)
    out_full$ATT_se   <- stats::sd(att_reps, na.rm = TRUE)

    overall_p <- {
      p1 <- if (theta_hat > 0) mean(att_reps <= 0, na.rm = TRUE)
      else               mean(att_reps >= 0, na.rm = TRUE)
      min(2*p1, 1)
    }

    ## bootstrapped event-study SE / CIs / p-values
    if (nrow(out_full$es) > 0) {
      periods <- out_full$es$` `
      boot_mat <- vapply(rep_list, function(x) {
        sapply(periods, function(p) {
          idx <- match(p, x$es$` `)
          if (!is.na(idx)) x$es$ATT[idx] else NA_real_
        })
      }, numeric(length(periods)))
      boot_mat <- t(boot_mat)

      boot_se  <- apply(boot_mat, 2, stats::sd, na.rm = TRUE)
      ci_lo    <- apply(boot_mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
      ci_hi    <- apply(boot_mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE)

      pvals_es <- vapply(
        seq_along(periods),
        function(j) {
          est <- out_full$es$ATT[j]
          p1  <- if (est > 0) mean(boot_mat[,j] <= 0, na.rm = TRUE)
          else         mean(boot_mat[,j] >= 0, na.rm = TRUE)
          min(2*p1, 1)
        },
        numeric(1)
      )

      out_full$es$`S.E.`     <- boot_se
      out_full$es$CI.lower   <- 2*out_full$es$ATT - ci_lo
      out_full$es$CI.upper   <- 2*out_full$es$ATT - ci_hi
      out_full$es$p.value    <- pvals_es
    }
  }

  ## ── Return object ───────────────────────────────────────────────────────
  res <- list(
    est.avg = data.frame(
      ATT.avg  = out_full$ATT,
      S.E.     = out_full$ATT_se,
      CI.lower = out_full$CI_lower,
      CI.upper = out_full$CI_upper,
      p.value  = overall_p
    ),
    est.att = out_full$es
  )
  class(res) <- "did_wrapper"
  res
}
