####
## diagtest
# diagnostic tests for no pre-trend, placebo effect and carryover effect
diagtest <- function(
    x, # fect object
    proportion = 0.3,
    pre.periods = NULL,
    f.threshold = NULL,
    tost.threshold = NULL,
    N_bar = NULL) {
    # get equivalence p values for two-one-sided-t tests
    tost <- function(coef, se, range) {
        z <- coef / se
        p1 <- 1 - pnorm((coef - range[1]) / se, lower.tail = TRUE) # left bound
        p2 <- 1 - pnorm((range[2] - coef) / se, lower.tail = TRUE) # right bound
        tost.p <- max(c(p1, p2))
        return(tost.p)
    }

    if (is.null(tost.threshold) == TRUE) {
        tost.threshold <- 0.36 * sqrt(x$sigma2.fect)
    }
    if (is.null(f.threshold) == TRUE) {
        f.threshold <- 0.6
    }

    # placebo test
    if (x$placeboTest == TRUE) {
        est.out <- x$est.placebo
        placebo.equiv.p <- tost(est.out[1], est.out[2], c(-tost.threshold, tost.threshold))
        out <- list(placebo.p = est.out[5], placebo.equiv.p = placebo.equiv.p)
    } # end of placebo test

    # carryover test
    if (x$carryoverTest == TRUE) {
        est.out <- x$est.carryover
        carryover.equiv.p <- tost(est.out[1], est.out[2], c(-tost.threshold, tost.threshold))
        out <- list(carryover.p = est.out[5], carryover.equiv.p = carryover.equiv.p)
    } # end of carryover test

    if (is.null(proportion) == TRUE) {
        proportion <- 0
    }

    max.pre.periods <- pre.pos <- NULL

    if (x$placeboTest == FALSE && x$carryoverTest == FALSE) {
        max.count <- max(x$count)
        max.pre.periods <- x$time[which(x$count >= max.count * proportion & x$time <= 0)]
        if (is.null(pre.periods) == TRUE) {
            pre.periods <- max.pre.periods
        } else if(length(pre.periods)>0) {
            pre.periods <- intersect(pre.periods[1]:pre.periods[length(pre.periods)], max.pre.periods)
        } else{
            pre.periods <- NA
        }
        pre.term <- pre.periods
        N_bar <- max(x$count[which(x$time %in% pre.periods)])
    }

    # testing no pre-trend
    if (x$placeboTest == FALSE & x$carryoverTest == FALSE & x$loo == TRUE) {
        max.pre.periods <- sum(x$time <= 0)
        pre.pos <- intersect(c(1:dim(x$pre.att.boot)[1]), which(x$time %in% pre.periods))
        res_boot <- x$pre.att.boot
        res_boot <- res_boot[, which(apply(!is.na(res_boot), 2, all))]
        if (length(pre.pos) == max.pre.periods) {
            pre.pos <- pre.pos[-1]
            #message("Cannot use full pre-treatment periods in F-test. The first period is removed.\n")
        }
        if (length(pre.pos) > 1) {
            res_boot <- res_boot[pre.pos, ]
        } else {
            res_boot <- t(as.matrix(res_boot[pre.pos, ]))
        }

        D <- as.matrix(x$pre.est.att[pre.pos, 1])
        coef_mat <- res_boot
        S <- cov(t(coef_mat)) ## * N_bar
        psi <- try(as.numeric(t(D) %*% solve(S) %*% D), silent = TRUE)
        
        if ("try-error" %in% class(psi)) {
            message("\n")
            #message("The estimated covariance matrix is irreversible.")
            message("F-test Failed. The estimated covariance matrix is singular.")
            message("\n")
            f.stat <- f.p <- f.equiv.p <- f.threshold <- NA
        } else {
            scale <- (N_bar - length(pre.pos)) / ((N_bar - 1) * length(pre.pos))
            ## F statistic

            if (scale <= 0) {
                message("Can't calculate the F statistic because of insufficient treated units.\n")
                f.stat <- NA
                f.p <- NA
                f.equiv.p <- NA
            } else {
                f.stat <- psi * scale
                f.p <- pf(f.stat,
                    df1 = length(pre.pos), df2 = N_bar - length(pre.pos),
                    lower.tail = FALSE
                )

                ## Equivalent F test
                f.equiv.p <- pf(f.stat,
                    df1 = length(pre.pos), df2 = N_bar - length(pre.pos),
                    ncp = N_bar * f.threshold
                )
            }
        }

        # TOST
        est.att <- x$pre.est.att[, c(1:2)]
        tost.equiv.p <- sapply(1:nrow(est.att), function(i) {
            return(tost(est.att[i, 1], est.att[i, 2], c(-tost.threshold, tost.threshold)))
        }) # keep the maximum p value

        tost.equiv.p <- max(tost.equiv.p)
        out <- list(
            f.stat = f.stat,
            f.p = f.p,
            f.threshold = f.threshold,
            f.equiv.p = f.equiv.p,
            df1 = length(pre.pos),
            df2 = N_bar - length(pre.pos),
            N_bar = N_bar,
            tost.equiv.p = tost.equiv.p,
            tost.threshold = tost.threshold
        )
    }

    if (x$placeboTest == FALSE & x$carryoverTest == FALSE & x$loo == FALSE) {
        max.pre.periods <- sum(x$time <= 0)
        pre.pos <- which(x$time %in% pre.periods)
        if (length(pre.pos) == max.pre.periods) {
            pre.pos <- pre.pos[-1]
            #message("Cannot use full pre-treatment periods in the F test. The first period is removed.\n")
        }

        res_boot <- x$att.boot
        res_boot <- res_boot[, which(apply(!is.na(res_boot), 2, all))]
        nboots <- ncol(res_boot)
        if (length(pre.pos) > 1) {
            res_boot <- res_boot[pre.pos, ]
        } else {
            res_boot <- t(as.matrix(res_boot[pre.pos, ]))
        }

        D <- as.matrix(x$est.att[pre.pos, 1])
        coef_mat <- res_boot

        S <- cov(t(coef_mat)) ## * N_bar

        psi <- try(as.numeric(t(D) %*% solve(S) %*% D), silent = TRUE)
        if ("try-error" %in% class(psi)) {
            message("\n")
            #message("The estimated covariance matrix is irreversible.")
            message("F-test Failed. The estimated covariance matrix is singular.")
            message("\n")
            f.stat <- f.p <- f.equiv.p <- f.threshold <- NA
        } else {
            scale <- (N_bar - length(pre.pos)) / ((N_bar - 1) * length(pre.pos))
            ## F statistic
            if (scale <= 0) {
                message("Can't calculate the F statistic because of insufficient treated units.\n")
                f.stat <- NA
                f.p <- NA
                f.equiv.p <- NA
            } else {
                f.stat <- psi * scale
                f.p <- pf(f.stat,
                    df1 = length(pre.pos), df2 = N_bar - length(pre.pos),
                    lower.tail = FALSE
                )

                ## Equivalent F test
                f.equiv.p <- pf(f.stat,
                    df1 = length(pre.pos), df2 = N_bar - length(pre.pos),
                    ncp = N_bar * f.threshold
                )
            }
        }

        # TOST

        est.att <- x$est.att[, c(1:2)]
        pos.zero <- which(x$time == 0)
        est.att <- est.att[pre.pos, , drop = FALSE]
        if (dim(est.att)[1] > 0) {
            tost.equiv.p <- sapply(1:nrow(est.att), function(i) {
                return(tost(est.att[i, 1], est.att[i, 2], c(-tost.threshold, tost.threshold)))
            }) # keep the maximum p value
            tost.equiv.p <- max(tost.equiv.p)
        } else {
            tost.equiv.p <- NA
        }


        out <- list(
            f.stat = f.stat,
            f.p = f.p,
            f.threshold = f.threshold,
            f.equiv.p = f.equiv.p,
            df1 = length(pre.pos),
            df2 = N_bar - length(pre.pos),
            N_bar = N_bar,
            tost.equiv.p = tost.equiv.p,
            tost.threshold = tost.threshold
        )
    }

    # end of testing no pre-trend

    return(out)
}
