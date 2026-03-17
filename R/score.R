###################################################################
## Shared scoring function for cross-validation and MSPE diagnostics
###################################################################

.score_residuals <- function(resid,
                             obs_weights = NULL,
                             time_index = NULL,
                             count_weights = NULL,
                             norm.para = NULL) {

    ## ---- input validation ---- ##
    n <- length(resid)
    if (n == 0L) stop("No residuals to score.")
    if (!is.null(obs_weights) && length(obs_weights) != n) {
        stop("obs_weights must have the same length as resid.")
    }
    if (!is.null(time_index) && length(time_index) != n) {
        stop("time_index must have the same length as resid.")
    }

    ## ---- Step 1: squared residuals and period weights ---- ##
    e2 <- resid^2

    if (!is.null(count_weights) && !is.null(time_index)) {
        w_period <- count_weights[time_index]
        names(w_period) <- NULL
    } else {
        w_period <- rep(1, n)
    }

    ## ---- Step 2: accumulations ---- ##
    if (is.null(obs_weights)) {
        ## unweighted path
        SSE    <- sum(e2)
        WSSE   <- sum(w_period * e2)
        GSSE   <- sum(log(e2))
        ll     <- w_period * e2
        ll_pos <- ll[which(ll > 0)]
        WGSSE  <- sum(log(ll_pos))
        ll_len <- length(ll_pos)
        mad_list <- e2
        denom  <- n
    } else {
        ## weighted path
        SSE    <- sum(obs_weights * e2)
        WSSE   <- sum(obs_weights * w_period * e2)
        GSSE   <- sum(obs_weights * log(e2))
        ll     <- obs_weights * w_period * e2
        ll_pos <- ll[which(ll > 0)]
        WGSSE  <- sum(log(ll_pos))
        ll_len <- length(ll_pos)
        mad_list <- obs_weights * e2
        denom  <- sum(obs_weights)
    }

    ## ---- Step 3: criteria ---- ##
    MSPE   <- SSE / denom
    WMSPE  <- WSSE / denom
    GMSPE  <- exp(GSSE / denom)
    WGMSPE <- if (ll_len > 0) exp(WGSSE / ll_len) else Inf
    MAD    <- median(abs(mad_list - median(mad_list)))

    ## ---- Step 4: Moment and GMoment ---- ##
    if (!is.null(time_index)) {
        resid_mean   <- tapply(resid, time_index, mean)
        resid_mean   <- abs(resid_mean)

        gm_mean <- function(x) exp(sum(log(x)) / length(x))
        resid_g_mean <- tapply(abs(resid), time_index, gm_mean)

        if (!is.null(count_weights)) {
            w_moment   <- count_weights[names(resid_mean)]
            names(w_moment) <- NULL
            w_gmoment  <- count_weights[names(resid_g_mean)]
            names(w_gmoment) <- NULL
            Moment     <- sum(w_moment * resid_mean) / sum(w_moment)
            GMoment    <- sum(w_gmoment * resid_g_mean) / sum(w_gmoment)
        } else {
            Moment  <- mean(resid_mean)
            GMoment <- mean(resid_g_mean)
        }
    } else {
        Moment  <- NA_real_
        GMoment <- NA_real_
    }

    ## ---- Step 5: normalization ---- ##
    if (!is.null(norm.para)) {
        scale <- norm.para[1]^2
        MSPE    <- MSPE * scale
        WMSPE   <- WMSPE * scale
        GMSPE   <- GMSPE * scale
        WGMSPE  <- WGMSPE * scale
        MAD     <- MAD * scale
        Moment  <- Moment * scale
        GMoment <- GMoment * scale
    }

    ## ---- Step 6: convenience scores and return ---- ##
    RMSE <- sqrt(MSPE)
    Bias <- mean(resid)
    c(MSPE = MSPE, WMSPE = WMSPE, GMSPE = GMSPE, WGMSPE = WGMSPE,
      MAD = MAD, Moment = Moment, GMoment = GMoment,
      RMSE = RMSE, Bias = Bias)
}
