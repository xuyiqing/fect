###################################################################
## Matrix Completion Function
###################################################################
fect_mc <- function(Y, # Outcome variable, (T*N) matrix
                    X, # Explanatory variables:  (T*N*p) array
                    D,
                    W,
                    I,
                    II,
                    T.on,
                    T.off = NULL,
                    T.on.carry = NULL,
                    T.on.balance = NULL,
                    balance.period = NULL,
                    lambda.cv = 1e10,
                    force,
                    hasF = 1,
                    hasRevs = 1,
                    tol, # tolerance level
                    max.iteration = 1000,
                    boot = FALSE, # bootstrapped sample
                    norm.para = NULL,
                    placeboTest = 0,
                    placebo.period = NULL,
                    carryoverTest = 0,
                    carryover.period = NULL,
                    calendar.enp.seq = NULL,
                    time.on.seq = NULL,
                    time.off.seq = NULL,
                    time.on.carry.seq = NULL,
                    time.on.balance.seq = NULL,
                    time.on.seq.W = NULL,
                    time.off.seq.W = NULL,
                    group.level = NULL,
                    group = NULL,
                    time.on.seq.group = NULL,
                    time.off.seq.group = NULL) {
    ## -------------------------------##
    ## Parsing data
    ## -------------------------------##
    carryover.pos <- placebo.pos <- na.pos <- NULL

    ## unit id and time
    TT <- dim(Y)[1]
    N <- dim(Y)[2]
    if (is.null(X) == FALSE) {
        p <- dim(X)[3]
    } else {
        p <- 0
    }

    ## replicate data
    YY <- Y
    YY[which(II == 0)] <- 0 ## reset to 0
    t.on <- c(T.on)

    ## initial fit using fastplm
    data.ini <- matrix(NA, (TT * N), (2 + 1 + p))
    data.ini[, 2] <- rep(1:N, each = TT) ## unit fe
    data.ini[, 3] <- rep(1:TT, N) ## time fe
    data.ini[, 1] <- c(Y) ## outcome
    if (p > 0) { ## covar
        for (i in 1:p) {
            data.ini[, (3 + i)] <- c(X[, , i])
        }
    }

    if (is.null(W)) {
        W.use <- as.matrix(0)
    } else {
        W.use <- W
        W.use[which(II == 0)] <- 0
    }

    ## observed Y0 indicator:
    oci <- which(c(II) == 1)
    if (!is.null(W)) {
        initialOut <- initialFit(data = data.ini, force = force, w = c(W), oci = oci)
    } else {
        initialOut <- initialFit(data = data.ini, force = force, w = NULL, oci = oci)
    }
    Y0 <- initialOut$Y0
    beta0 <- initialOut$beta0
    if (p > 0 && sum(is.na(beta0)) > 0) {
        beta0[which(is.na(beta0))] <- 0
    }

    ## -------------------------------##
    ## ----------- Main Algorithm ----------- ##
    ## -------------------------------##

    lambda.norm <- eigen.all <- NULL
    if (boot == FALSE) {
        Y.lambda <- YY - Y0
        Y.lambda[which(II == 0)] <- 0
        eigen.all <- svd(Y.lambda / (TT * N))$d
        lambda.norm <- lambda.cv / max(eigen.all)
    }

    validX <- 1 ## no multi-colinearity
    ## matrix completion
    est.best <- inter_fe_mc(YY, Y0, X, II, W.use, beta0, hasF, lambda.cv, force, tol, max.iteration)
    validX <- est.best$validX
    validF <- est.best$validF
    est.fect <- NULL
    if (boot == FALSE) {
        est.fect <- inter_fe_ub(YY, Y0, X, II, W.use, beta0, 0, force = force, tol, max.iteration)
    }

    ## ------------------------------##
    ## ----------- Summarize -------------- ##
    ## ------------------------------##

    ## -------------------------------##
    ##   ATT and Counterfactuals     ##
    ## -------------------------------##

    ## we first adjustment for normalization
    if (!is.null(norm.para)) {
        Y <- Y * norm.para[1]

        ## output of estimates
        est.best$mu <- est.best$mu * norm.para[1]
        if (force %in% c(1, 3)) {
            est.best$alpha <- est.best$alpha * norm.para[1]
        }
        if (force %in% c(2, 3)) {
            est.best$xi <- est.best$xi * norm.para[1]
        }
        # if (p>0) {
        #    est.best$beta <- est.best$beta * norm.para[1]
        # }
        est.best$residuals <- est.best$residuals * norm.para[1]
        est.best$fit <- est.best$fit * norm.para[1]
        if (boot == FALSE) {
            est.fect$fit <- est.fect$fit * norm.para[1]
        }
        est.fect$sigma2 <- est.fect$sigma2 * norm.para[1]
    }

    ## 0. revelant parameters
    if (p > 0) {
        na.pos <- is.nan(est.best$beta)
        beta <- est.best$beta
        if (sum(na.pos) > 0) {
            beta[na.pos] <- NA
        }
    } else {
        beta <- NA
    }

    ## 1. estimated att and counterfactuals
    eff <- Y - est.best$fit
    complete.index <- which(!is.na(eff))
    missing.index <- which(is.na(eff))
    if (length(missing.index) > 0) {
        I[missing.index] <- 0
        II[missing.index] <- 0
    }
    if (0 %in% I) {
        eff[which(I == 0)] <- NA
    }
    complete.index <- which(!is.na(eff))
    att.avg <- sum(eff[complete.index] * D[complete.index]) / (sum(D[complete.index]))

    att.avg.balance <- NA
    if (!is.null(balance.period)) {
        complete.index2 <- which(!is.na(T.on.balance))
        att.avg.balance <- sum(eff[complete.index2] * D[complete.index2]) / (sum(D[complete.index2]))
    }

    # weighted effect
    att.avg.W <- NA
    if (!is.null(W)) {
        att.avg.W <- sum(eff[complete.index] * D[complete.index] * W[complete.index]) / (sum(D[complete.index] * W[complete.index]))
    }

    ## att.avg.unit
    tr.pos <- which(apply(D, 2, sum) > 0)
    att.unit <- sapply(1:length(tr.pos), function(vec) {
        return(sum(eff[, tr.pos[vec]] * D[, tr.pos[vec]]) / sum(D[, tr.pos[vec]]))
    })
    att.avg.unit <- mean(att.unit, na.rm = TRUE)


    equiv.att.avg <- eff.equiv <- NULL
    if (boot == FALSE) {
        eff.equiv <- Y - est.fect$fit
        if (0 %in% I) {
            eff.equiv[which(I == 0)] <- NA
        }
        complete.index <- which(!is.na(eff.equiv))
        equiv.att.avg <- sum(eff.equiv[complete.index] * D[complete.index]) / (sum(D[complete.index]))
    }

    ## 2. rmse for treated units' observations under control
    tr <- which(apply(D, 2, sum) > 0)
    co <- which(apply(D, 2, sum) == 0)
    tr.co <- which((as.matrix(1 - D[, tr]) * as.matrix(II[, tr])) == 1)
    eff.tr <- as.matrix(eff[, tr])
    v.eff.tr <- eff.tr[tr.co]
    rmse <- sqrt(mean(v.eff.tr^2))

    ## 3. unbalanced output
    Y.ct.full <- est.best$fit
    if (0 %in% I) {
        eff[which(I == 0)] <- NA
        est.best$fit[which(I == 0)] <- NA
    }
    res.full <- est.best$residuals
    est.best$residuals[which(II == 0)] <- NA

    est.best$sigma2 <- mean(c(est.best$residuals[which(II == 1)])^2) ## mean squared error of residuals


    ## 4. dynamic effects
    t.on <- c(T.on)
    eff.v <- c(eff) ## a vector
    rm.pos1 <- which(is.na(eff.v))
    rm.pos2 <- which(is.na(t.on))
    eff.v.use1 <- eff.v
    t.on.use <- t.on
    n.on.use <- rep(1:N, each = TT)

    eff.equiv.v <- NULL
    if (boot == FALSE) {
        eff.equiv.v <- c(eff.equiv)
    }

    if (NA %in% eff.v | NA %in% t.on) {
        eff.v.use1 <- eff.v[-c(rm.pos1, rm.pos2)]
        t.on.use <- t.on[-c(rm.pos1, rm.pos2)]
        n.on.use <- n.on.use[-c(rm.pos1, rm.pos2)]
        if (boot == FALSE) {
            eff.equiv.v <- eff.equiv.v[-c(rm.pos1, rm.pos2)]
        }
    }

    pre.pos <- which(t.on.use <= 0)
    eff.pre <- cbind(eff.v.use1[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
    colnames(eff.pre) <- c("eff", "period", "unit")

    pre.sd <- eff.pre.equiv <- NULL
    if (boot == FALSE) {
        eff.pre.equiv <- cbind(eff.equiv.v[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
        colnames(eff.pre.equiv) <- c("eff.equiv", "period", "unit")

        pre.sd <- tapply(eff.pre.equiv[, 1], eff.pre.equiv[, 2], sd)
        pre.sd <- cbind(pre.sd, sort(unique(eff.pre.equiv[, 2])), table(eff.pre.equiv[, 2]))
        colnames(pre.sd) <- c("sd", "period", "count")
    }


    time.on <- sort(unique(t.on.use))
    att.on <- as.numeric(tapply(eff.v.use1, t.on.use, mean)) ## NA already removed
    count.on <- as.numeric(table(t.on.use))

    if (!is.null(time.on.seq)) {
        count.on.med <- att.on.med <- rep(NA, length(time.on.seq))
        att.on.med[which(time.on.seq %in% time.on)] <- att.on
        count.on.med[which(time.on.seq %in% time.on)] <- count.on
        att.on <- att.on.med
        count.on <- count.on.med
        time.on <- time.on.seq
    }

    ## weighted treatment effect

    if (!is.null(W)) {
        W.v <- c(W)
        rm.pos.W <- which(is.na(W))
        if (NA %in% eff.v | NA %in% t.on | NA %in% W.v) {
            eff.v.use.W <- eff.v[-c(rm.pos1, rm.pos2, rm.pos.W)]
            W.v.use <- W.v[-c(rm.pos1, rm.pos2, rm.pos.W)]
            t.on.use.W <- t.on[-c(rm.pos1, rm.pos2, rm.pos.W)]
            n.on.use.W <- n.on.use[-c(rm.pos1, rm.pos2, rm.pos.W)]
        } else {
            eff.v.use.W <- eff.v.use1
            t.on.use.W <- t.on.use
            n.on.use.W <- n.on.use
            W.v.use <- W.v
        }
        time.on.W <- sort(unique(t.on.use.W))
        att.on.sum.W <- as.numeric(tapply(eff.v.use.W * W.v.use, t.on.use.W, sum)) ## NA already removed
        W.on.sum <- as.numeric(tapply(W.v.use, t.on.use.W, sum))
        att.on.W <- att.on.sum.W / W.on.sum
        count.on.W <- as.numeric(table(t.on.use.W))

        if (!is.null(time.on.seq.W)) {
            att.on.sum.med.W <- W.on.sum.med <- count.on.med.W <- att.on.med.W <- rep(NA, length(time.on.seq.W))
            att.on.sum.med.W[which(time.on.seq.W %in% time.on.W)] <- att.on.sum.W
            att.on.med.W[which(time.on.seq.W %in% time.on.W)] <- att.on.W
            count.on.med.W[which(time.on.seq.W %in% time.on.W)] <- count.on.W
            W.on.sum.med[which(time.on.seq.W %in% time.on.W)] <- W.on.sum
            att.on.sum.W <- att.on.sum.med.W
            att.on.W <- att.on.med.W
            count.on.W <- count.on.med.W
            time.on.W <- time.on.seq.W
            W.on.sum <- W.on.sum.med
        }
    } else {
        att.on.sum.med.W <- att.on.sum.W <- count.on.med.W <- att.on.med.W <- W.on.sum.med <- att.on.W <- count.on.W <- time.on.W <- W.on.sum <- NULL
    }

    ## 4.1 carryover effect
    carry.att <- NULL
    if (!is.null(T.on.carry)) {
        t.on.carry <- c(T.on.carry)
        rm.pos4 <- which(is.na(t.on.carry))
        t.on.carry.use <- t.on.carry

        if (NA %in% eff.v | NA %in% t.on.carry) {
            eff.v.use3 <- eff.v[-c(rm.pos1, rm.pos4)]
            t.on.carry.use <- t.on.carry[-c(rm.pos1, rm.pos4)]
        }

        carry.time <- sort(unique(t.on.carry.use))
        carry.att <- as.numeric(tapply(eff.v.use3, t.on.carry.use, mean)) ## NA already removed

        if (!is.null(time.on.carry.seq)) {
            carry.att.med <- rep(NA, length(time.on.carry.seq))
            carry.att.med[which(time.on.carry.seq %in% carry.time)] <- carry.att
            carry.att <- carry.att.med
            carry.time <- time.on.carry.seq
        }
    }

    ## 4.2 balance effect
    balance.att <- NULL
    if (!is.null(balance.period)) {
        t.on.balance <- c(T.on.balance)
        rm.pos4 <- which(is.na(t.on.balance))
        t.on.balance.use <- t.on.balance

        if (NA %in% eff.v | NA %in% t.on.balance) {
            eff.v.use3 <- eff.v[-c(rm.pos1, rm.pos4)]
            t.on.balance.use <- t.on.balance[-c(rm.pos1, rm.pos4)]
        }

        balance.time <- sort(unique(t.on.balance.use))
        balance.att <- as.numeric(tapply(eff.v.use3, t.on.balance.use, mean)) ## NA already removed
        balance.count <- as.numeric(table(t.on.balance.use))

        if (!is.null(time.on.balance.seq)) {
            balance.att.med <- rep(NA, length(time.on.balance.seq))
            balance.count.med <- rep(0, length(time.on.balance.seq))
            balance.att.med[which(time.on.balance.seq %in% balance.time)] <- balance.att
            if (length(balance.count) > 0) {
                balance.count.med[which(time.on.balance.seq %in% balance.time)] <- balance.count
            }
            balance.count <- balance.count.med
            balance.att <- balance.att.med
            balance.time <- time.on.balance.seq
        }

        # placebo for balanced samples
        if (!is.null(placebo.period) && placeboTest == 1) {
            if (length(placebo.period) == 1) {
                balance.placebo.pos <- which(balance.time == placebo.period)
                balance.att.placebo <- balance.att[balance.placebo.pos]
            } else {
                balance.placebo.pos <- which(balance.time >= placebo.period[1] & balance.time <= placebo.period[2])
                balance.att.placebo <- sum(balance.att[balance.placebo.pos] * balance.count[balance.placebo.pos]) / sum(balance.count[balance.placebo.pos])
            }
        }
    }

    ## 5. placebo effect, if placeboTest == 1
    if (!is.null(placebo.period) && placeboTest == 1) {
        if (length(placebo.period) == 1) {
            placebo.pos <- which(time.on == placebo.period)
            att.placebo <- att.on[placebo.pos]
        } else {
            placebo.pos <- which(time.on >= placebo.period[1] & time.on <= placebo.period[2])
            att.placebo <- sum(att.on[placebo.pos] * count.on[placebo.pos]) / sum(count.on[placebo.pos])
        }

        if (!is.null(W)) {
            if (length(placebo.period) == 1) {
                placebo.pos.W <- which(time.on.W == placebo.period)
                att.placebo.W <- att.on.W[placebo.pos.W]
            } else {
                placebo.pos.W <- which(time.on.W >= placebo.period[1] & time.on.W <= placebo.period[2])
                att.placebo.W <- sum(att.on.sum.W[placebo.pos.W]) / sum(W.on.sum[placebo.pos.W])
            }
        }
    }

    eff.off.equiv <- off.sd <- eff.off <- NULL

    ## 6. switch-off effects
    if (hasRevs == 1) {
        t.off <- c(T.off)
        rm.pos3 <- which(is.na(t.off))
        eff.v.use2 <- eff.v
        t.off.use <- t.off

        if (NA %in% eff.v | NA %in% t.off) {
            eff.v.use2 <- eff.v[-c(rm.pos1, rm.pos3)]
            t.off.use <- t.off[-c(rm.pos1, rm.pos3)]
        }

        off.pos <- which(t.off.use > 0)
        eff.off <- cbind(eff.v.use2[off.pos], t.off.use[off.pos], n.on.use[off.pos])
        colnames(eff.off) <- c("eff", "period", "unit")

        if (boot == FALSE) {
            eff.off.equiv <- cbind(eff.equiv.v[off.pos], t.off.use[off.pos], n.on.use[off.pos])
            colnames(eff.off.equiv) <- c("off.equiv", "period", "unit")

            off.sd <- tapply(eff.off.equiv[, 1], eff.off.equiv[, 2], sd)
            off.sd <- cbind(off.sd, sort(unique(eff.off.equiv[, 2])), table(eff.off.equiv[, 2]))
            colnames(off.sd) <- c("sd", "period", "count")
        }

        time.off <- sort(unique(t.off.use))
        att.off <- as.numeric(tapply(eff.v.use2, t.off.use, mean)) ## NA already removed
        count.off <- as.numeric(table(t.off.use))

        if (!is.null(time.off.seq)) {
            count.off.med <- att.off.med <- rep(NA, length(time.off.seq))
            att.off.med[which(time.off.seq %in% time.off)] <- att.off
            count.off.med[which(time.off.seq %in% time.off)] <- count.off
            att.off <- att.off.med
            count.off <- count.off.med
            time.off <- time.off.seq
        }

        if (!is.null(W)) {
            if (NA %in% eff.v | NA %in% t.off | NA %in% W.v) {
                eff.v.use2.W <- eff.v[-c(rm.pos1, rm.pos3, rm.pos.W)]
                W.v.use2 <- W.v[-c(rm.pos1, rm.pos3, rm.pos.W)]
                t.off.use.W <- t.off[-c(rm.pos1, rm.pos3, rm.pos.W)]
            } else {
                eff.v.use2.W <- eff.v.use2
                t.off.use.W <- t.off.use
                W.v.use2 <- W.v
            }

            time.off.W <- sort(unique(t.off.use.W))
            att.off.sum.W <- as.numeric(tapply(eff.v.use2.W * W.v.use2, t.off.use.W, sum))
            W.off.sum <- as.numeric(tapply(W.v.use2, t.off.use.W, sum))
            att.off.W <- att.off.sum.W / W.off.sum ## NA already removed
            count.off.W <- as.numeric(table(t.off.use.W))

            if (!is.null(time.off.seq.W)) {
                att.off.sum.med.W <- W.off.sum.med <- count.off.med.W <- att.off.med.W <- rep(NA, length(time.off.seq.W))
                att.off.sum.med.W[which(time.off.seq.W %in% time.off.W)] <- att.off.sum.W
                att.off.med.W[which(time.off.seq.W %in% time.off.W)] <- att.off.W
                count.off.med.W[which(time.off.seq.W %in% time.off.W)] <- count.off.W
                W.off.sum.med[which(time.off.seq.W %in% time.off.W)] <- W.off.sum
                att.off.sum.W <- att.off.sum.med.W
                att.off.W <- att.off.med.W
                count.off.W <- count.off.med.W
                time.off.W <- time.off.seq.W
                W.off.sum <- W.off.sum.med
            }
        } else {
            W.off.sum.med <- W.off.sum <- att.off.sum.W <- att.off.sum.med.W <- count.off.med.W <- att.off.med.W <- count.off.med.W <- att.off.W <- count.off.W <- time.off.W <- NULL
        }
    }

    ## 7. carryover effects
    if (!is.null(carryover.period) && carryoverTest == 1 && hasRevs) {
        if (length(carryover.period) == 1) {
            carryover.pos <- which(time.off == carryover.period)
            att.carryover <- att.off[carryover.pos]
        } else {
            carryover.pos <- which(time.off >= carryover.period[1] & time.off <= carryover.period[2])
            att.carryover <- sum(att.off[carryover.pos] * count.off[carryover.pos]) / sum(count.off[carryover.pos])
        }
        if (!is.null(W)) {
            if (length(carryover.period) == 1) {
                carryover.pos.W <- which(time.off.W == carryover.period)
                att.carryover.W <- att.off.W[carryover.pos.W]
            } else {
                carryover.pos.W <- which(time.off.W >= carryover.period[1] & time.off.W <= carryover.period[2])
                att.carryover.W <- sum(att.off.sum.W[carryover.pos.W]) / sum(W.off.sum[carryover.pos.W])
            }
        }
    }

    ## 9. loess HTE by time
    D.missing <- D
    D.missing[which(D == 0)] <- NA
    eff.calendar <- apply(eff * D.missing, 1, mean, na.rm = TRUE)
    N.calendar <- apply(!is.na(eff * D.missing), 1, sum)
    T.calendar <- c(1:TT)
    if (sum(!is.na(eff.calendar)) > 1) {
        # loess fit
        if (!is.null(calendar.enp.seq)) {
            if (length(calendar.enp.seq) == 1 & is.na(calendar.enp.seq)) {
                calendar.enp.seq <- NULL
            }
        }
        if (is.null(calendar.enp.seq)) {
            loess.fit <- suppressWarnings(try(loess(eff.calendar ~ T.calendar, weights = N.calendar), silent = TRUE))
        } else {
            loess.fit <- suppressWarnings(try(loess(eff.calendar ~ T.calendar, weights = N.calendar, enp.target = calendar.enp.seq), silent = TRUE))
        }
        if ("try-error" %in% class(loess.fit)) {
            eff.calendar.fit <- eff.calendar
            calendar.enp <- NULL
        } else {
            eff.calendar.fit <- eff.calendar
            eff.calendar.fit[which(!is.na(eff.calendar))] <- loess.fit$fit
            calendar.enp <- loess.fit$enp
        }
    } else {
        eff.calendar.fit <- eff.calendar
        calendar.enp <- NULL
    }

    ## 8. cohort effects
    if (!is.null(group)) {
        cohort <- cbind(c(group), c(D), c(eff.v))
        rm.pos <- unique(c(rm.pos1, which(cohort[, 2] == 0)))
        cohort <- cohort[-rm.pos, ]

        g.level <- sort(unique(cohort[, 1]))
        raw.group.att <- as.numeric(tapply(cohort[, 3], cohort[, 1], mean))

        group.att <- rep(NA, length(group.level))
        group.att[which(group.level %in% g.level)] <- raw.group.att

        # by-group dynamic effects
        group.level.name <- names(group.level)

        group.output <- list()
        for (i in c(1:length(group.level))) {
            sub.group <- group.level[i]
            sub.group.name <- group.level.name[i]

            ## by-group dynamic effects
            t.on.sub <- c(T.on[which(group == sub.group)])
            eff.v.sub <- c(eff[which(group == sub.group)]) ## a vector
            rm.pos1.sub <- which(is.na(eff.v.sub))
            rm.pos2.sub <- which(is.na(t.on.sub))
            eff.v.use1.sub <- eff.v.sub
            t.on.use.sub <- t.on.sub
            if (NA %in% eff.v.sub | NA %in% t.on.sub) {
                eff.v.use1.sub <- eff.v.sub[-c(rm.pos1.sub, rm.pos2.sub)]
                t.on.use.sub <- t.on.sub[-c(rm.pos1.sub, rm.pos2.sub)]
            }
            if (length(t.on.use.sub) > 0) {
                time.on.sub <- sort(unique(t.on.use.sub))
                att.on.sub <- as.numeric(tapply(
                    eff.v.use1.sub,
                    t.on.use.sub,
                    mean
                )) ## NA already removed
                count.on.sub <- as.numeric(table(t.on.use.sub))
            } else {
                time.on.sub <- att.on.sub <- count.on.sub <- NULL
            }

            if (!is.null(time.on.seq.group)) {
                count.on.med.sub <- att.on.med.sub <- rep(NA, length(time.on.seq.group[[sub.group.name]]))
                time.on.seq.sub <- time.on.seq.group[[sub.group.name]]
                att.on.med.sub[which(time.on.seq.sub %in% time.on.sub)] <- att.on.sub
                count.on.med.sub[which(time.on.seq.sub %in% time.on.sub)] <- count.on.sub
                att.on.sub <- att.on.med.sub
                count.on.sub <- count.on.med.sub
                time.on.sub <- time.on.seq.sub
            }
            suboutput <- list(
                att.on = att.on.sub,
                time.on = time.on.sub,
                count.on = count.on.sub
            )

            ## placebo effect, if placeboTest == 1
            if (!is.null(placebo.period) && placeboTest == 1) {
                if (length(placebo.period) == 1) {
                    placebo.pos.sub <- which(time.on.sub == placebo.period)
                    if (length(placebo.pos.sub) > 0) {
                        att.placebo.sub <- att.on.sub[placebo.pos.sub]
                    } else {
                        att.placebo.sub <- NULL
                    }
                } else {
                    placebo.pos.sub <- which(time.on.sub >= placebo.period[1] & time.on.sub <= placebo.period[2])
                    if (length(placebo.pos.sub) > 0) {
                        att.placebo.sub <- sum(att.on.sub[placebo.pos.sub] * count.on.sub[placebo.pos.sub]) / sum(count.on.sub[placebo.pos.sub])
                    } else {
                        att.placebo.sub <- NULL
                    }
                }
                suboutput <- c(suboutput, list(att.placebo = att.placebo.sub))
            }

            ## T.off
            if (hasRevs == 1) {
                t.off.sub <- c(T.off[which(group == sub.group)])
                rm.pos3.sub <- which(is.na(t.off.sub))
                eff.v.use2.sub <- eff.v.sub
                t.off.use.sub <- t.off.sub
                if (NA %in% eff.v.sub | NA %in% t.off.sub) {
                    eff.v.use2.sub <- eff.v.sub[-c(rm.pos1.sub, rm.pos3.sub)]
                    t.off.use.sub <- t.off.sub[-c(rm.pos1.sub, rm.pos3.sub)]
                }
                if (length(t.off.use.sub) > 0) {
                    time.off.sub <- sort(unique(t.off.use.sub))
                    att.off.sub <- as.numeric(tapply(eff.v.use2.sub, t.off.use.sub, mean)) ## NA already removed
                    count.off.sub <- as.numeric(table(t.off.use.sub))
                } else {
                    time.off.sub <- att.off.sub <- count.off.sub <- NULL
                }

                if (!is.null(time.off.seq.group)) {
                    count.off.med.sub <- att.off.med.sub <- rep(NA, length(time.off.seq.group[[sub.group.name]]))
                    time.off.seq.sub <- time.off.seq.group[[sub.group.name]]
                    att.off.med.sub[which(time.off.seq.sub %in% time.off.sub)] <- att.off.sub
                    count.off.med.sub[which(time.off.seq.sub %in% time.off.sub)] <- count.off.sub
                    att.off.sub <- att.off.med.sub
                    count.off.sub <- count.off.med.sub
                    time.off.sub <- time.off.seq.sub
                }
                suboutput <- c(suboutput, list(
                    att.off = att.off.sub,
                    count.off = count.off.sub,
                    time.off = time.off.sub
                ))

                if (!is.null(carryover.period) && carryoverTest == 1) {
                    if (length(carryover.period) == 1) {
                        carryover.pos.sub <- which(time.off.sub == carryover.period)
                        if (length(carryover.pos.sub) > 0) {
                            att.carryover.sub <- att.off.sub[carryover.pos.sub]
                        } else {
                            att.carryover.sub <- NULL
                        }
                    } else {
                        carryover.pos.sub <- which(time.off.sub >= carryover.period[1] & time.off.sub <= carryover.period[2])
                        if (length(carryover.pos.sub) > 0) {
                            att.carryover.sub <- sum(att.off.sub[carryover.pos.sub] * count.off.sub[carryover.pos.sub]) / sum(count.off.sub[carryover.pos.sub])
                        } else {
                            att.carryover.sub <- NULL
                        }
                    }
                    suboutput <- c(suboutput, list(att.carryover = att.carryover.sub))
                }
            }
            group.output[[sub.group.name]] <- suboutput
        }
    }

    ## -------------------------------##
    ## Storage
    ## -------------------------------##
    out <- list(
        ## main results
        method = "mc",
        T.on = T.on,
        Y.ct = est.best$fit,
        Y.ct.full = Y.ct.full,
        eff = eff,
        I = I,
        II = II,
        D = D,
        Y = Y,
        X = X,
        att.avg = att.avg,
        att.avg.unit = att.avg.unit,
        ## supporting
        force = force,
        T = TT,
        N = N,
        tr = tr,
        co = co,
        p = p,
        lambda.cv = lambda.cv,
        lambda.norm = lambda.norm,
        eigen.all = eigen.all,
        beta = beta,
        est = est.best,
        mu = est.best$mu,
        validX = validX,
        validF = validF,
        niter = est.best$niter,
        time = time.on,
        att = att.on,
        count = count.on,
        eff.calendar = eff.calendar,
        N.calendar = N.calendar,
        eff.calendar.fit = eff.calendar.fit,
        calendar.enp = calendar.enp,
        eff.pre = eff.pre,
        eff.pre.equiv = eff.pre.equiv,
        pre.sd = pre.sd,
        rmse = rmse,
        res = est.best$residuals,
        sigma2 = est.best$sigma2,
        sigma2.fect = est.fect$sigma2
    )

    if (hasRevs == 1) {
        out <- c(out, list(
            time.off = time.off,
            att.off = att.off,
            count.off = count.off,
            eff.off = eff.off,
            eff.off.equiv = eff.off.equiv,
            off.sd = off.sd
        ))
    }

    if (!is.null(T.on.carry)) {
        out <- c(out, list(carry.att = carry.att, carry.time = carry.time))
    }

    if (!is.null(balance.period)) {
        out <- c(out, list(
            balance.att = balance.att, balance.time = balance.time,
            balance.count = balance.count, balance.avg.att = att.avg.balance
        ))
        if (!is.null(placebo.period) && placeboTest == 1) {
            out <- c(out, list(balance.att.placebo = balance.att.placebo))
        }
    }

    if (!is.null(W)) {
        out <- c(out, list(
            W = W,
            att.avg.W = att.avg.W,
            att.on.sum.W = att.on.sum.W,
            att.on.W = att.on.W,
            count.on.W = count.on.W,
            time.on.W = time.on.W,
            W.on.sum = W.on.sum
        ))
        if (hasRevs == 1) {
            out <- c(out, list(
                att.off.sum.W = att.off.sum.W,
                att.off.W = att.off.W,
                count.off.W = count.off.W,
                time.off.W = time.off.W,
                W.off.sum = W.off.sum
            ))
        }
        if (!is.null(placebo.period) && placeboTest == 1) {
            out <- c(out, list(att.placebo.W = att.placebo.W))
        }

        if (!is.null(carryover.period) && carryoverTest == 1) {
            out <- c(out, list(att.carryover.W = att.carryover.W))
        }
    }

    # if (boot == FALSE) {
    #    out <- c(out, list(equiv.att.avg = equiv.att.avg))
    # }

    if (force == 1) {
        out <- c(out, list(alpha = est.best$alpha))
    } else if (force == 2) {
        out <- c(out, list(xi = est.best$xi))
    } else if (force == 3) {
        out <- c(out, list(alpha = est.best$alpha, xi = est.best$xi))
    }

    if (!is.null(placebo.period) && placeboTest == 1) {
        out <- c(out, list(att.placebo = att.placebo))
    }

    if (!is.null(carryover.period) && carryoverTest == 1) {
        out <- c(out, list(att.carryover = att.carryover))
    }

    if (!is.null(group)) {
        out <- c(out, list(
            group.att = group.att,
            group.output = group.output
        ))
    }
    return(out)
} ## mc functions ends
