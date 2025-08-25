###################################################################
## Cross-validation
###################################################################
fect_cv <- function(Y, # Outcome variable, (T*N) matrix
                    X, # Explanatory variables:  (T*N*p) array
                    D, #  Indicator for treated unit (tr==1)
                    W,
                    I,
                    II,
                    T.on,
                    T.off = NULL,
                    T.on.carry = NULL,
                    T.on.balance = NULL,
                    balance.period = NULL,
                    method = "ife",
                    criterion = "mspe",
                    k = 5, # CV time
                    cv.prop = 0.1,
                    cv.treat = TRUE,
                    cv.nobs = 3,
                    cv.donut = 1,
                    min.T0 = 5,
                    r = 0, # initial number of factors considered if CV==1
                    r.end,
                    proportion = 0,
                    nlambda = 10,
                    lambda = NULL,
                    force,
                    hasRevs = 1,
                    tol, # tolerance level
                    max.iteration = 1000,
                    norm.para = NULL,
                    group.level = NULL,
                    group = NULL) {
    ## -------------------------------##
    ## Parsing data
    ## -------------------------------##
    placebo.pos <- na.pos <- NULL
    ## unit id and time
    TT <- dim(Y)[1]
    N <- dim(Y)[2]
    if (is.null(X) == FALSE) {
        p <- dim(X)[3]
    } else {
        p <- 0
        X <- array(0, dim = c(1, 1, 0))
    }

    if (is.null(W)) {
        W.use <- as.matrix(0)
        use_weight <- 0
    } else {
        use_weight <- 1
        W.use <- W
        W.use[which(II == 0)] <- 0
    }

    ## replicate data
    YY <- Y
    YY[which(II == 0)] <- 0 ## reset to 0
    if (use_weight) {
        WW <- W
        WW[which(II == 0)] <- 0 ## reset to 0
    }
    t.on <- c(T.on)
    T0.min <- min(apply(II, 2, sum))

    D.c <- apply(D, 2, function(vec) {
        cumsum(vec)
    })
    D.c <- ifelse(D.c > 0, 1, 0)
    D.sum <- colSums(D.c)
    tr <- which(D.sum >= 1)
    Ntr <- length(tr)
    co <- which(D.sum == 0)
    Nco <- length(co)

    ##  --------- initial fit using fastplm --------- ##
    data.ini <- matrix(NA, (TT * N), (2 + 1 + p))
    data.ini[, 2] <- rep(1:N, each = TT) ## unit fe
    data.ini[, 3] <- rep(1:TT, N) ## time fe
    data.ini[, 1] <- c(Y) ## outcome
    if (p > 0) { ## covar
        for (i in 1:p) {
            data.ini[, (3 + i)] <- c(X[, , i])
        }
    }
    ## observed Y0 indicator:
    oci <- which(c(II) == 1)
    if (use_weight == 1) {
        initialOut <- initialFit(
            data = data.ini,
            force = force,
            w = c(W),
            oci = oci
        )
    } else {
        initialOut <- initialFit(
            data = data.ini,
            force = force,
            oci = oci
        )
    }

    Y0 <- initialOut$Y0
    beta0 <- initialOut$beta0
    if (p > 0 && sum(is.na(beta0)) > 0) {
        beta0[which(is.na(beta0))] <- 0
    }

    ## ------------- restrictions on candidate hyper parameters ---------- ##
    obs.con <- (sum(II) - r.end * (N + TT) + r.end^2 - p) <= 0
    if (obs.con) {
        while ((sum(II) - r.end * (N + TT) + r.end^2 - p) <= 0) {
            r.end <- r.end - 1
        }
    }
    if (r.end >= T0.min) {
        if (method %in% c("both", "ife", "gsynth")) {
            message("Factor number should not be greater than ", T0.min - 1, "\n", sep = "")
        }
        r.end <- T0.min - 1
    } else {
        if (obs.con) {
            if (method %in% c("both", "ife", "gsynth")) {
                message("Factor number should not be greater than ", r.end, "\n", sep = "")
            }
        }
    }

    ## -------------------------------##
    ## ----------- Main Algorithm ----------- ##
    ## -------------------------------##

    validX <- 1 ## no multi-colinearity
    CV.out.ife <- CV.out.mc <- NULL

    ## ----------------------------------------------------##
    ##         Cross-validation of r and lambda           ##
    ## ----------------------------------------------------##

    r.max <- min(TT, r.end)
    r.cv <- 0 ## initial value

    if (method %in% c("ife", "both", "gsynth") && FALSE) {
        r.cv <- 0
        est.best <- inter_fe_ub(YY, Y0,
            X, II, W.use, beta0,
            0,
            force = force,
            tol, max.iteration
        )
        message("Cross validation cannot be performed since available pre-treatment records of treated units are too few. So set r.cv = 0.\n ")
    } else {
        r.old <- r ## save the minimal number of factors
        message("Cross-validating ...", "\n")
        if (criterion == "mspe") {
            message("Criterion: Mean Squared Prediction Error\n")
        } else if (criterion == "wmspe") {
            message("Criterion: Weighted Mean Squared Prediction Error\n")
        } else if (criterion == "gmspe") {
            message("Criterion: Geometric Mean Squared Prediction Error\n")
        } else if (criterion == "wgmspe") {
            message("Criterion: Weighted Geometric Mean Squared Prediction Error\n")
        } else if (criterion == "mad") {
            message("Criterion: Median Absolute Deviation\n")
        } else if (criterion == "moment") {
            message("Criterion: Moment Conditions\n")
        } else if (criterion == "gmoment") {
            message("Criterion: Geometric Moment Conditions\n")
        } else if (criterion == "pc") {
            message("Criterion: PC\n")
        }

        ## for gsynth, use the cross-validation function in fect_gsynth
        if (method == "gsynth") {
            message("Interactive fixed effects model...\n")
            out <- fect_gsynth(
                Y = Y, D = D, X = X, W = W, I = I, II = II,
                T.on = T.on, T.off = T.off,
                T.on.balance = T.on.balance,
                balance.period = balance.period,
                r = r, r.end = r.end, CV = TRUE,
                force = force, hasRevs = hasRevs,
                tol = tol, boot = 0,
                norm.para = norm.para,
                group.level = group.level, group = group
            )
            return(out)
        }

        ## ----- ##
        ## ------------- initialize ------------ ##
        ## ----- ##

        cv.pos <- which(t.on <= 0)
        t.on.cv <- t.on[cv.pos]
        count.on.cv <- as.numeric(table(t.on.cv))
        ## tot.id <- which(c(II)==1) ## observed control data
        ## cv.count <- ceiling((sum(II)*sum(II))/sum(I))
        rm.count <- floor(sum(II) * cv.prop)
        cv.count <- sum(II) - rm.count

        # ociCV <- matrix(NA, cv.count, k) ## store indicator
        # rmCV <- matrix(NA, rm.count, k) ## removed indicator
        ociCV <- list()
        rmCV <- list()
        estCV <- NULL ## used for mspe
        if (use_weight == 1) {
            W.rmCV <- list()
        }

        Y0CV <- array(NA, dim = c(TT, N, k)) ## store initial Y0
        if (p > 0) {
            beta0CV <- array(NA, dim = c(p, 1, k))
        } else {
            beta0CV <- array(0, dim = c(1, 0, k)) ## store initial beta0
        }

        ## cv.id.all <- c()
        flag <- 0
        for (i in 1:k) {
            cv.n <- 0
            repeat{
                cv.n <- cv.n + 1
                # cv.id <- cv.sample(II, as.integer(sum(II) - cv.count))
                get.cv <- cv.sample(II, D,
                    count = rm.count,
                    cv.count = cv.nobs,
                    cv.treat = cv.treat,
                    cv.donut = cv.donut
                )
                cv.id <- get.cv$cv.id
                ## cv.id <- sample(oci, as.integer(sum(II) - cv.count), replace = FALSE)
                II.cv.valid <- II.cv <- II
                II.cv[cv.id] <- 0
                II.cv.valid[cv.id] <- -1
                ## ziyi: if certain rows or columns doesn't satisfy con1 or con2,
                ## replace the row or column of II.cv using the corresponding rows or columns in II

                con1 <- sum(apply(II.cv, 1, sum) >= 1) == TT
                con2 <- sum(apply(II.cv, 2, sum) >= min.T0) == N

                if (con1 == TRUE & con2 == TRUE) {
                    break
                }

                if (cv.n >= 200) {
                    flag <- 1
                    # message("Some units have too few pre-treatment observations. Remove them automatically in Cross-Validation.")
                    keep.1 <- which(apply(II.cv, 1, sum) < 1)
                    keep.2 <- which(apply(II.cv, 2, sum) < min.T0)
                    II.cv[keep.1, ] <- II[keep.1, ]
                    II.cv[, keep.2] <- II[, keep.2]
                    II.cv.valid[keep.1, ] <- II[keep.1, ]
                    II.cv.valid[, keep.2] <- II[, keep.2]
                    cv.id <- which(II.cv.valid != II)
                    break
                }
            }

            if (length(cv.id) == 0) {
                stop("Some units have too few pre-treatment observations. Set a larger \"cv.prop\" or set \"cv.treat\" to FALSE.")
            }

            rmCV[[i]] <- cv.id
            ocicv <- setdiff(oci, cv.id)
            ociCV[[i]] <- ocicv
            if (use_weight) {
                W.estCV <- list()
            }

            if (cv.n < 200) {
                estCV <- c(estCV, list(get.cv$est.id))
            } else {
                cv.id.old <- get.cv$cv.id
                cv.diff <- setdiff(cv.id.old, cv.id)
                estCV <- c(estCV, list(setdiff(get.cv$est.id, cv.diff)))
            }

            if (use_weight == 0) {
                initialOutCv <- initialFit(
                    data = data.ini,
                    force = force,
                    oci = ocicv
                )
            } else {
                initialOutCv <- initialFit(
                    data = data.ini,
                    force = force,
                    w = c(W),
                    oci = ocicv
                )
            }

            Y0CV[, , i] <- initialOutCv$Y0

            if (p > 0) {
                beta0cv <- initialOutCv$beta0
                if (sum(is.na(beta0cv)) > 0) {
                    beta0cv[which(is.na(beta0cv))] <- 0
                }
                beta0CV[, , i] <- beta0cv
            }
        }

        if (flag == 1) {
            message("Some units have too few pre-treatment observations. Remove them automatically in Cross-Validation.\n")
        }

        ## get count matrix
        if (use_weight == 0) {
            count.T.cv <- count.T.cv.old <- table(T.on)
            count.T.cv.old <- count.T.cv <- count.T.cv[which(as.numeric(names(count.T.cv)) <= 0)]
            cv.prop.cut <- max(count.T.cv.old) * proportion
            cv.drop.index <- which(count.T.cv.old <= cv.prop.cut)

            count.T.cv <- count.T.cv / mean(count.T.cv)
            name.count.T.cv <- names(count.T.cv)
            count.T.cv <- c(count.T.cv, median(count.T.cv))
            names(count.T.cv) <- c(name.count.T.cv, "Control")
            count.T.cv[cv.drop.index] <- 0 # set weights to 0 for the periods when the number of treated observations is less than proportion
        }
        if (use_weight == 1) {
            count.T.cv <- count.T.cv.old <- aggregate(c(W), by = list(relative = c(T.on)), FUN = sum)
            count.T.cv.old <- count.T.cv <- count.T.cv[which(count.T.cv[, "relative"] <= 0), ]
            cv.prop.cut <- max(count.T.cv.old[, 2]) * proportion
            cv.drop.index <- which(count.T.cv.old[, 2] <= cv.prop.cut)

            name.count.T.cv <- count.T.cv[, 1]
            count.T.cv <- count.T.cv[, 2] / mean(count.T.cv[, 2])
            count.T.cv <- c(count.T.cv, median(count.T.cv))
            names(count.T.cv) <- c(name.count.T.cv, "Control")
            count.T.cv[cv.drop.index] <- 0
        }


        ##  --------------------------------------------- ##
        ##  ---------------- cross validation for ife model ------------------  ##
        ##  --------------------------------------------- ##

        if (method %in% c("ife", "both")) {
            message("Interactive fixed effects model...\n")

            r.pc <- est.pc.best <- MSPE.best <- WMSPE.best <- MSPE.pc.best <- NULL
            gmoment.best <- moment.best <- MAD.best <- GMSPE.best <- WGMSPE.best <- NULL

            if (criterion == "PC") {
                CV.out.ife <- matrix(NA, (r.max - r.old + 1), 6)
                colnames(CV.out.ife) <- c("r", "sigma2", "IC", "PC", "MSPTATT", "MSE")
            } else {
                CV.out.ife <- matrix(NA, (r.max - r.old + 1), 13)
                colnames(CV.out.ife) <- c(
                    "r", "sigma2", "IC", "PC",
                    "MSPE", "WMSPE", "GMSPE", "WGMSPE", "MAD", "Moment", "GMoment", "MSPTATT", "MSE"
                )
            }

            CV.out.ife[, "r"] <- c(r.old:r.max)
            CV.out.ife[, "PC"] <- CV.out.ife[, "GMoment"] <- CV.out.ife[, "Moment"] <- CV.out.ife[, "MAD"] <- CV.out.ife[, "MSPE"] <- CV.out.ife[, "WMSPE"] <- CV.out.ife[, "GMSPE"] <- CV.out.ife[, "WGMSPE"] <- 1e20

            for (i in 1:dim(CV.out.ife)[1]) { ## cross-validation loop starts
                ## inter FE based on control, before & after
                r <- CV.out.ife[i, "r"]
                ## k <- 5
                if (criterion %in% c("mspe", "wmspe", "gmspe", "wgmspe", "mad", "moment")) {
                    SSE <- 0
                    WSSE <- 0
                    GSSE <- 0
                    WGSSE <- 0
                    ll.length <- 0
                    moment.list <- c()
                    index.moment.list <- c()
                    MAD.list <- c()
                    for (ii in 1:k) {
                        II.cv <- II
                        II.cv[rmCV[[ii]]] <- 0
                        YY.cv <- YY
                        YY.cv[rmCV[[ii]]] <- 0
                        if (use_weight) {
                            W.use2 <- W.use
                            W.use2[rmCV[[ii]]] <- 0
                        } else {
                            W.use2 <- as.matrix(0)
                        }
                        est.cv.fit <- inter_fe_ub(
                            YY.cv, as.matrix(Y0CV[, , ii]), X, II.cv,
                            W.use2, as.matrix(beta0CV[, , ii]),
                            r, force, tol, max.iteration
                        )$fit
                        index.cv <- as.character(T.on[estCV[[ii]]])
                        index.cv[which(is.na(index.cv))] <- "Control"
                        weight.cv <- count.T.cv[index.cv]
                        names(weight.cv) <- NULL
                        if (use_weight == 0) {
                            SSE <- SSE + sum((YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2)
                            WSSE <- WSSE + sum(weight.cv * (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2)
                            GSSE <- GSSE + sum(log((YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2))
                            ll <- weight.cv * (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2
                            ll <- ll[which(ll > 0)]
                            WGSSE <- WGSSE + sum(log(ll))
                            ll.length <- ll.length + length(ll)
                            MAD.list <- c(MAD.list, (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2)
                        } else {
                            W.estCV[[ii]] <- WW[estCV[[ii]]]
                            # print(WW[estCV[[ii]]])
                            SSE <- SSE + sum(WW[estCV[[ii]]] * (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2)
                            WSSE <- WSSE + sum(WW[estCV[[ii]]] * weight.cv * (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2)
                            GSSE <- GSSE + sum(WW[estCV[[ii]]] * log((YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2))
                            ll <- WW[estCV[[ii]]] * weight.cv * (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2
                            ll <- ll[which(ll > 0)]
                            WGSSE <- WGSSE + sum(log(ll))
                            ll.length <- ll.length + length(ll)
                            MAD.list <- c(MAD.list, WW[estCV[[ii]]] * (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2)
                        }
                        # moment conditions
                        moment.list <- c(moment.list, (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]]))
                        index.moment.list <- c(index.moment.list, index.cv)

                        # resid.mean <- tapply((YY[estCV[[ii]]]-est.cv.fit[estCV[[ii]]]), index.cv, mean)
                        # resid.mean <- abs(resid.mean)
                        # weight.cv <- count.T.cv[names(resid.mean)]
                        # names(weight.cv) <- NULL
                        # moment <- c(moment, sum(weight.cv*resid.mean)/sum(weight.cv))
                    }

                    if (use_weight == 0) {
                        MSPE <- SSE / (length(unlist(estCV)))
                        WMSPE <- WSSE / (length(unlist(estCV)))
                        GMSPE <- exp(GSSE / (length(unlist(estCV))))
                        WGMSPE <- exp(WGSSE / ll.length)
                        MAD <- median(abs(MAD.list - median(MAD.list)))
                    } else {
                        MSPE <- SSE / (sum(unlist(W.estCV)))
                        WMSPE <- WSSE / (sum(unlist(W.estCV)))
                        GMSPE <- exp(GSSE / (sum(unlist(W.estCV))))
                        WGMSPE <- exp(WGSSE / ll.length)
                        MAD <- median(abs(MAD.list - median(MAD.list)))
                    }
                    # moment
                    resid.mean <- tapply(moment.list, index.moment.list, mean)
                    resid.mean <- abs(resid.mean)
                    # g-moment
                    gm_mean <- function(x) {
                        exp(sum(log(x)) / length(x))
                    }
                    resid.g.mean <- tapply(abs(moment.list), index.moment.list, gm_mean)
                    weight.cv.g <- count.T.cv[names(resid.g.mean)]
                    weight.cv <- count.T.cv[names(resid.mean)]
                    names(weight.cv) <- NULL
                    names(weight.cv.g) <- NULL
                    moment <- sum(weight.cv * resid.mean) / sum(weight.cv)
                    gmoment <- sum(weight.cv.g * resid.g.mean) / sum(weight.cv)
                }

                est.cv <- inter_fe_ub(
                    YY,
                    Y0,
                    X,
                    II,
                    W.use,
                    beta0,
                    r,
                    force,
                    tol,
                    max.iteration
                ) ## overall
                sigma2 <- est.cv$sigma2
                IC <- est.cv$IC
                PC <- est.cv$PC

                eff.v.cv <- c(Y - est.cv$fit)[cv.pos]
                meff <- as.numeric(tapply(eff.v.cv, t.on.cv, mean))
                MSPTATT <- sum(meff^2 * count.on.cv) / sum(count.on.cv)
                MSE <- sum(eff.v.cv^2) / length(eff.v.cv)

                if (!is.null(norm.para)) {
                    if (criterion %in% c("mspe", "wmspe", "gmspe", "wgmspe", "mad", "moment", "gmoment")) {
                        MSPE <- MSPE * (norm.para[1]^2)
                        WMSPE <- WMSPE * (norm.para[1]^2)
                        GMSPE <- GMSPE * (norm.para[1]^2)
                        WGMSPE <- WGMSPE * (norm.para[1]^2)
                        MAD <- MAD * (norm.para[1]^2)
                        moment <- moment * (norm.para[1]^2)
                        gmoment <- gmoment * (norm.para[1]^2)
                    }
                    sigma2 <- sigma2 * (norm.para[1]^2)
                    IC <- est.cv$IC - log(est.cv$sigma2) + log(sigma2)
                    PC <- PC * (norm.para[1]^2)
                }

                if (criterion == "mspe") {
                    if ((min(CV.out.ife[, "MSPE"]) - MSPE) > 0.01 * min(CV.out.ife[, "MSPE"])) {
                        ## at least 1% improvement for MPSE
                        MSPE.best <- MSPE
                        est.best <- est.cv
                        r.cv <- r
                    } else {
                        if (r == r.cv + 1) message("*")
                    }
                } else if (criterion == "wmspe") {
                    if ((min(CV.out.ife[, "WMSPE"]) - WMSPE) > 0.01 * min(CV.out.ife[, "WMSPE"])) {
                        ## at least 1% improvement for MPSE
                        WMSPE.best <- WMSPE
                        est.best <- est.cv
                        r.cv <- r
                    } else {
                        if (r == r.cv + 1) message("*")
                    }
                } else if (criterion == "gmspe") {
                    if ((min(CV.out.ife[, "GMSPE"]) - GMSPE) > 0.01 * min(CV.out.ife[, "GMSPE"])) {
                        ## at least 1% improvement for MPSE
                        GMSPE.best <- GMSPE
                        est.best <- est.cv
                        r.cv <- r
                    } else {
                        if (r == r.cv + 1) message("*")
                    }
                } else if (criterion == "wgmspe") {
                    if ((min(CV.out.ife[, "WGMSPE"]) - WGMSPE) > 0.01 * min(CV.out.ife[, "WGMSPE"])) {
                        ## at least 1% improvement for MPSE
                        WGMSPE.best <- WGMSPE
                        est.best <- est.cv
                        r.cv <- r
                    } else {
                        if (r == r.cv + 1) message("*")
                    }
                } else if (criterion == "mad") {
                    if ((min(CV.out.ife[, "MAD"]) - MAD) > 0.01 * min(CV.out.ife[, "MAD"])) {
                        MAD.best <- MAD
                        est.best <- est.cv
                        r.cv <- r
                    } else {
                        if (r == r.cv + 1) message("*")
                    }
                } else if (criterion == "moment") {
                    if ((min(CV.out.ife[, "Moment"]) - moment) > 0.01 * min(CV.out.ife[, "Moment"])) {
                        moment.best <- moment
                        est.best <- est.cv
                        r.cv <- r
                    } else {
                        if (r == r.cv + 1) message("*")
                    }
                } else if (criterion == "gmoment") {
                    if ((min(CV.out.ife[, "GMoment"]) - gmoment) > 0.01 * min(CV.out.ife[, "GMoment"])) {
                        gmoment.best <- gmoment
                        est.best <- est.cv
                        r.cv <- r
                    } else {
                        if (r == r.cv + 1) message("*")
                    }
                } else if (criterion == "pc") {
                    if (PC < min(CV.out.ife[, "PC"])) {
                        est.pc.best <- est.cv
                        r.pc <- r
                    }
                }

                if (criterion != "pc") {
                    CV.out.ife[i, 2:12] <- c(sigma2, IC, PC, MSPE, WMSPE, GMSPE, WGMSPE, MAD, moment, MSPTATT, MSE)
                } else {
                    CV.out.ife[i, 2:6] <- c(sigma2, IC, PC, MSPTATT, MSE)
                }

                if (criterion == "pc") {
                    message("r = ", r, "; sigma2 = ",
                        sprintf("%.5f", sigma2), "; IC = ",
                        sprintf("%.5f", IC), "; PC = ",
                        sprintf("%.5f", PC), "; MSPTATT = ",
                        sprintf("%.5f", MSPTATT), "; MSE = ",
                        sprintf("%.5f", MSE),
                        sep = ""
                    )
                } else {
                    # message("r = ",r, "; sigma2 = ",
                    #    sprintf("%.5f",sigma2), "; IC = ",
                    #    sprintf("%.5f",IC), "; PC = ",
                    #    sprintf("%.5f",PC), "; MSPE = ",
                    #    sprintf("%.5f",MSPE), "; GMSPE = ",
                    #    sprintf("%.5f",GMSPE), "; Moment = ",
                    #    sprintf("%.5f",moment), "; MSPTATT = ",
                    #    sprintf("%.5f",MSPTATT), "; MSE = ",
                    #    sprintf("%.5f",MSE), sep="")
                    message(
                        "r = ", r, "; sigma2 = ",
                        sprintf("%.5f", sigma2), "; IC = ",
                        sprintf("%.5f", IC), "; PC = ",
                        sprintf("%.5f", PC), "; MSPE = ",
                        sprintf("%.5f", MSPE)
                    )
                }
            } ## end of while: search for r_star over

            # MSPE.best <- min(CV.out[,"MSPE"])
            # PC.best <- min(CV.out[,"PC"])

            ## compare
            if (criterion == "both") {
                if (r.cv > r.pc) {
                    message("\n Factor number selected via cross validation may be larger than the true number. Using the PC criterion.\n\n ")
                    r.cv <- r.pc
                    est.best <- est.pc.best
                    MSPE.best <- MSPE.pc.best
                }
                est.best.ife <- est.best
                MSPE.best.ife <- MSPE.best
            } else if (criterion == "pc") {
                est.best.ife <- est.pc.best
                r.cv <- r.pc
            } else {
                est.best.ife <- est.best
                MSPE.best.ife <- MSPE.best
                WMSPE.best.ife <- WMSPE.best
                GMSPE.best.ife <- GMSPE.best
                WGMSPE.best.ife <- WGMSPE.best
                MAD.best.ife <- MAD.best
                moment.best.ife <- moment.best
                gmoment.best.ife <- gmoment.best
            }

            if (r > (TT - 1)) {
                message(" (r hits maximum)")
            }
            message("\n r* = ", r.cv, sep = "")
        }

        ##  ------------------------------------- ##
        ##  ---------------- cross validation for mc ---------------  ##
        ##  ------------------------------------- ##
        if (method %in% c("mc", "both")) {
            message("Matrix completion method...\n")
            eigen.all <- NULL
            if (is.null(lambda) || length(lambda) == 1) {
                ## create the hyper-parameter sequence
                ## biggest candidate lambda
                ## Y.lambda <- YY
                Y.lambda <- YY - Y0
                ## Y.lambda[which(II == 0)] <- Y0[which(II == 0)]
                Y.lambda[which(II == 0)] <- 0
                if (use_weight) {
                    Y.lambda <- Y.lambda * W
                }
                eigen.all <- svd(Y.lambda / (TT * N))$d
                lambda.max <- log10(max(eigen.all))
                lambda <- rep(NA, nlambda)
                lambda.by <- 3 / (nlambda - 2)
                for (i in 1:(nlambda - 1)) {
                    lambda[i] <- 10^(lambda.max - (i - 1) * lambda.by)
                }
                lambda[nlambda] <- 0
            } else {
                Y.lambda <- YY - Y0
                Y.lambda[which(II == 0)] <- 0
                if (use_weight) {
                    Y.lambda <- Y.lambda * W
                }
                eigen.all <- svd(Y.lambda / (TT * N))$d
            }

            ## store all MSPE
            MSPE.best <- WMSPE.best <- GMSPE.best <- WGMSPE.best <- MAD.best <- moment.best <- gmoment.best <- NULL
            CV.out.mc <- matrix(NA, length(lambda), 10)
            colnames(CV.out.mc) <- c("lambda.norm", "MSPE", "WMSPE", "GMSPE", "WGMSPE", "MAD", "Moment", "GMoment", "MSPTATT", "MSE")
            CV.out.mc[, "lambda.norm"] <- c(lambda / max(eigen.all))
            CV.out.mc[, "GMoment"] <- CV.out.mc[, "Moment"] <- CV.out.mc[, "MAD"] <- CV.out.mc[, "WGMSPE"] <- CV.out.mc[, "WGMSPE"] <- CV.out.mc[, "GMSPE"] <- CV.out.mc[, "WMSPE"] <- CV.out.mc[, "MSPE"] <- 1e20

            break_count <- 0
            break_check <- 0
            for (i in 1:length(lambda)) {
                ## k <- 5
                SSE <- 0
                WSSE <- 0
                GSSE <- 0
                WGSSE <- 0
                ll.length <- 0
                moment.list <- c()
                index.moment.list <- c()
                MAD.list <- c()

                for (ii in 1:k) {
                    II.cv <- II
                    II.cv[rmCV[[ii]]] <- 0
                    YY.cv <- YY
                    YY.cv[rmCV[[ii]]] <- 0
                    if (use_weight) {
                        W.use2 <- W.use
                        W.use2[rmCV[[ii]]] <- 0
                    } else {
                        W.use2 <- as.matrix(0)
                    }
                    est.cv.fit <- inter_fe_mc(
                        YY.cv, as.matrix(Y0CV[, , ii]),
                        X, II.cv, W.use2, as.matrix(beta0CV[, , ii]),
                        1, lambda[i], force, tol, max.iteration
                    )$fit
                    index.cv <- as.character(T.on[estCV[[ii]]])
                    index.cv[which(is.na(index.cv))] <- "Control"
                    weight.cv <- count.T.cv[index.cv]
                    names(weight.cv) <- NULL
                    if (use_weight == 0) {
                        SSE <- SSE + sum((YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2)
                        WSSE <- WSSE + sum(weight.cv * (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2)
                        GSSE <- GSSE + sum(log((YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2))
                        ll <- weight.cv * (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2
                        ll <- ll[which(ll > 0)]
                        WGSSE <- WGSSE + sum(log(ll))
                        ll.length <- ll.length + length(ll)
                        MAD.list <- c(MAD.list, (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2)
                    } else {
                        W.estCV[[ii]] <- WW[estCV[[ii]]]
                        SSE <- SSE + sum(WW[estCV[[ii]]] * (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2)
                        WSSE <- WSSE + sum(WW[estCV[[ii]]] * weight.cv * (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2)
                        GSSE <- GSSE + sum(WW[estCV[[ii]]] * log((YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2))
                        ll <- WW[estCV[[ii]]] * weight.cv * (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2
                        ll <- ll[which(ll > 0)]
                        WGSSE <- WGSSE + sum(log(ll))
                        ll.length <- ll.length + length(ll)
                        MAD.list <- c(MAD.list, WW[estCV[[ii]]] * (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]])^2)
                    }
                    # moment conditions
                    moment.list <- c(moment.list, (YY[estCV[[ii]]] - est.cv.fit[estCV[[ii]]]))
                    index.moment.list <- c(index.moment.list, index.cv)
                }
                if (use_weight == 0) {
                    MSPE <- SSE / (length(unlist(estCV)))
                    WMSPE <- WSSE / (length(unlist(estCV)))
                    GMSPE <- exp(GSSE / (length(unlist(estCV))))
                    WGMSPE <- exp(WGSSE / ll.length)
                    MAD <- median(abs(MAD.list - median(MAD.list)))
                } else {
                    MSPE <- SSE / (sum(unlist(W.estCV)))
                    WMSPE <- WSSE / (sum(unlist(W.estCV)))
                    GMSPE <- exp(GSSE / (sum(unlist(W.estCV))))
                    WGMSPE <- exp(WGSSE / ll.length)
                    MAD <- median(abs(MAD.list - median(MAD.list)))
                }

                # moment
                resid.mean <- tapply(moment.list, index.moment.list, mean)
                resid.mean <- abs(resid.mean)
                # g-moment
                gm_mean <- function(x) {
                    exp(sum(log(x)) / length(x))
                }
                resid.g.mean <- tapply(abs(moment.list), index.moment.list, gm_mean)
                weight.cv.g <- count.T.cv[names(resid.g.mean)]
                weight.cv <- count.T.cv[names(resid.mean)]
                names(weight.cv) <- NULL
                names(weight.cv.g) <- NULL
                moment <- sum(weight.cv * resid.mean) / sum(weight.cv)
                gmoment <- sum(weight.cv.g * resid.g.mean) / sum(weight.cv)

                est.cv <- inter_fe_mc(
                    YY, Y0, X, II, W.use, beta0,
                    1, lambda[i],
                    force, tol, max.iteration
                ) ## overall

                eff.v.cv <- c(Y - est.cv$fit)[cv.pos]
                meff <- as.numeric(tapply(eff.v.cv, t.on.cv, mean))
                MSPTATT <- sum(meff^2 * count.on.cv) / sum(count.on.cv)
                MSE <- sum(eff.v.cv^2) / length(eff.v.cv)

                if (!is.null(norm.para)) {
                    MSPE <- MSPE * (norm.para[1]^2)
                    WMSPE <- WMSPE * (norm.para[1]^2)
                    GMSPE <- GMSPE * (norm.para[1]^2)
                    WGMSPE <- WGMSPE * (norm.para[1]^2)
                    MAD <- MAD * (norm.para[1]^2)
                    moment <- moment * (norm.para[1]^2)
                    gmoment <- gmoment * (norm.para[1]^2)
                }

                if (criterion == "mspe") {
                    if ((min(CV.out.mc[, "MSPE"]) - MSPE) > 0.01 * min(CV.out.mc[, "MSPE"])) {
                        ## at least 1% improvement for MPSE
                        MSPE.best <- MSPE
                        est.best <- est.cv
                        lambda.cv <- lambda[i]
                        break_count <- 0
                        break_check <- 0
                    } else {
                        if (i > 1) {
                            if (lambda.cv == lambda[i - 1]) {
                                message("*")
                                break_check <- 1
                                break_count <- 0
                            }
                        }
                    }
                } else if (criterion == "wmspe") {
                    if ((min(CV.out.mc[, "WMSPE"]) - WMSPE) > 0.01 * min(CV.out.mc[, "WMSPE"])) {
                        ## at least 1% improvement for MPSE
                        WMSPE.best <- WMSPE
                        est.best <- est.cv
                        lambda.cv <- lambda[i]
                        break_check <- 0
                        break_count <- 0
                    } else {
                        if (i > 1) {
                            if (lambda.cv == lambda[i - 1]) {
                                message("*")
                                break_check <- 1
                                break_count <- 0
                            }
                        }
                    }
                } else if (criterion == "gmspe") {
                    if ((min(CV.out.mc[, "GMSPE"]) - GMSPE) > 0.01 * min(CV.out.mc[, "GMSPE"])) {
                        ## at least 1% improvement for MPSE
                        GMSPE.best <- GMSPE
                        est.best <- est.cv
                        lambda.cv <- lambda[i]
                        break_check <- 0
                        break_count <- 0
                    } else {
                        if (i > 1) {
                            if (lambda.cv == lambda[i - 1]) {
                                message("*")
                                break_check <- 1
                                break_count <- 0
                            }
                        }
                    }
                } else if (criterion == "wgmspe") {
                    if ((min(CV.out.mc[, "WGMSPE"]) - WGMSPE) > 0.01 * min(CV.out.mc[, "WGMSPE"])) {
                        ## at least 1% improvement for MPSE
                        WGMSPE.best <- WGMSPE
                        est.best <- est.cv
                        lambda.cv <- lambda[i]
                        break_check <- 0
                        break_count <- 0
                    } else {
                        if (i > 1) {
                            if (lambda.cv == lambda[i - 1]) {
                                message("*")
                                break_check <- 1
                                break_count <- 0
                            }
                        }
                    }
                } else if (criterion == "mad") {
                    if ((min(CV.out.mc[, "MAD"]) - MAD) > 0.01 * min(CV.out.mc[, "MAD"])) {
                        ## at least 1% improvement for MPSE
                        MAD.best <- MAD
                        est.best <- est.cv
                        lambda.cv <- lambda[i]
                        break_check <- 0
                        break_count <- 0
                    } else {
                        if (i > 1) {
                            if (lambda.cv == lambda[i - 1]) {
                                message("*")
                                break_check <- 1
                                break_count <- 0
                            }
                        }
                    }
                } else if (criterion == "moment") {
                    if ((min(CV.out.mc[, "Moment"]) - moment) > 0.01 * min(CV.out.mc[, "Moment"])) {
                        ## at least 1% improvement for MPSE
                        moment.best <- moment
                        est.best <- est.cv
                        lambda.cv <- lambda[i]
                        break_check <- 0
                        break_count <- 0
                    } else {
                        if (i > 1) {
                            if (lambda.cv == lambda[i - 1]) {
                                message("*")
                                break_check <- 1
                                break_count <- 0
                            }
                        }
                    }
                } else if (criterion == "gmoment") {
                    if ((min(CV.out.mc[, "GMoment"]) - gmoment) > 0.01 * min(CV.out.mc[, "GMoment"])) {
                        ## at least 1% improvement for MPSE
                        gmoment.best <- gmoment
                        est.best <- est.cv
                        lambda.cv <- lambda[i]
                        break_check <- 0
                        break_count <- 0
                    } else {
                        if (i > 1) {
                            if (lambda.cv == lambda[i - 1]) {
                                message("*")
                                break_check <- 1
                                break_count <- 0
                            }
                        }
                    }
                }

                if (break_check == 1) {
                    break_count <- break_count + 1
                }

                CV.out.mc[i, 2:10] <- c(MSPE, WMSPE, GMSPE, WGMSPE, MAD, moment, gmoment, MSPTATT, MSE)
                message("lambda.norm = ",
                    sprintf("%.5f", lambda[i] / max(eigen.all)), "; MSPE = ",
                    sprintf("%.5f", MSPE), "; MSPTATT = ",
                    sprintf("%.5f", MSPTATT), "; MSE = ",
                    sprintf("%.5f", MSE),
                    sep = ""
                )
                if (break_count == 3) {
                    break
                }
            }
            est.best.mc <- est.best
            MSPE.best.mc <- MSPE.best
            WMSPE.best.mc <- WMSPE.best
            GMSPE.best.mc <- GMSPE.best
            WGMSPE.best.mc <- WGMSPE.best
            MAD.best.mc <- MAD.best
            moment.best.mc <- moment.best
            gmoment.best.mc <- gmoment.best
            message("\n lambda.norm* = ", lambda.cv / max(eigen.all), sep = "")
            message("\n")
        }
    } ## End of Cross-Validation


    if (method == "ife") {
        est.best <- est.best.ife
        validF <- ifelse(r.cv > 0, 1, 0)
    } else if (method == "mc") {
        est.best <- est.best.mc
        validF <- est.best$validF
    } else {
        if (criterion == "mspe") {
            if (MSPE.best.ife <= MSPE.best.mc) {
                est.best <- est.best.ife
                validF <- ifelse(r.cv > 0, 1, 0)
                method <- "ife"
            } else {
                est.best <- est.best.mc
                validF <- est.best$validF
                method <- "mc"
            }
        }
        if (criterion == "wmspe") {
            if (WMSPE.best.ife <= WMSPE.best.mc) {
                est.best <- est.best.ife
                validF <- ifelse(r.cv > 0, 1, 0)
                method <- "ife"
            } else {
                est.best <- est.best.mc
                validF <- est.best$validF
                method <- "mc"
            }
        }
        if (criterion == "gmspe") {
            if (GMSPE.best.ife <= GMSPE.best.mc) {
                est.best <- est.best.ife
                validF <- ifelse(r.cv > 0, 1, 0)
                method <- "ife"
            } else {
                est.best <- est.best.mc
                validF <- est.best$validF
                method <- "mc"
            }
        }
        if (criterion == "wgmspe") {
            if (WGMSPE.best.ife <= WGMSPE.best.mc) {
                est.best <- est.best.ife
                validF <- ifelse(r.cv > 0, 1, 0)
                method <- "ife"
            } else {
                est.best <- est.best.mc
                validF <- est.best$validF
                method <- "mc"
            }
        }
        if (criterion == "mad") {
            if (MAD.best.ife <= MAD.best.mc) {
                est.best <- est.best.ife
                validF <- ifelse(r.cv > 0, 1, 0)
                method <- "ife"
            } else {
                est.best <- est.best.mc
                validF <- est.best$validF
                method <- "mc"
            }
        }
        if (criterion == "moment") {
            if (moment.best.ife <= moment.best.mc) {
                est.best <- est.best.ife
                validF <- ifelse(r.cv > 0, 1, 0)
                method <- "ife"
            } else {
                est.best <- est.best.mc
                validF <- est.best$validF
                method <- "mc"
            }
        }
        if (criterion == "gmoment") {
            if (gmoment.best.ife <= gmoment.best.mc) {
                est.best <- est.best.ife
                validF <- ifelse(r.cv > 0, 1, 0)
                method <- "ife"
            } else {
                est.best <- est.best.mc
                validF <- est.best$validF
                method <- "mc"
            }
        }
        message("\n Recommended method through cross-validation: ", method, sep = "")
        message("\n")
    }
    validX <- est.best$validX

    ## ------------------------------##
    ## ----------- Summarize -------------- ##
    ## ------------------------------##

    ## 00. run a fect to obtain residuals
    if (method == "ife") {
        if (r.cv == 0) {
            est.fect <- est.best
        } else {
            est.fect <- inter_fe_ub(YY, Y0, X, II, W.use, beta0, 0, force = force, tol, max.iteration)
        }
    } else {
        est.fect <- inter_fe_ub(YY, Y0, X, II, W.use, beta0, 0, force = force, tol, max.iteration)
    }


    ## -------------------------------##
    ##   ATT and Counterfactuals     ##
    ## -------------------------------##

    ## we first adjustment for normalization
    if (!is.null(norm.para)) {
        Y <- Y * norm.para[1]
        if (method == "ife") {
            ## variance of the error term
            sigma2 <- est.best$sigma2 * (norm.para[1]^2)
            IC <- est.best$IC - log(est.best$sigma2) + log(sigma2)
            PC <- est.best$PC * (norm.para[1]^2)
            est.best$sigma2 <- sigma2
            est.best$IC <- IC
            est.best$PC <- PC
        }

        ## output of estimates
        est.best$mu <- est.best$mu * norm.para[1]
        if (method == "ife" && r.cv > 0) {
            est.best$lambda <- est.best$lambda * norm.para[1]
            est.best$VNT <- est.best$VNT * norm.para[1]
        }
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
        est.fect$fit <- est.fect$fit * norm.para[1]
        est.fect$sigma2 <- est.fect$sigma2 * norm.para[1]
    }

    ## 0. revelant parameters
    sigma2 <- IC <- PC <- NULL
    if (method == "ife") {
        sigma2 <- est.best$sigma2
        IC <- est.best$IC
        PC <- est.best$PC
    }

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
    Y.ct.equiv <- Y.ct <- NULL
    Y.ct <- est.best$fit
    eff <- Y - Y.ct
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
    att.avg.unit <- mean(att.unit)

    eff.equiv <- Y - est.fect$fit
    equiv.att.avg <- sum(eff.equiv * D) / (sum(D))

    ## 2. rmse for treated units' observations under control
    tr <- which(apply(D, 2, sum) > 0)
    tr.co <- which((as.matrix(1 - D[, tr]) * as.matrix(II[, tr])) == 1)
    eff.tr <- as.matrix(eff[, tr])
    v.eff.tr <- eff.tr[tr.co]
    rmse <- sqrt(mean(v.eff.tr^2))

    ## 3. unbalanced output
    if (0 %in% I) {
        eff[which(I == 0)] <- NA
        est.best$fit[which(I == 0)] <- NA
        ## eff.equiv[which(I == 0)] <- NA
    }
    est.best$residuals[which(II == 0)] <- NA
    if (method == "mc") {
        est.best$sigma2 <- mean(c(est.best$residuals[which(II == 1)])^2) ## mean squared error of residuals
    }

    ## 4. dynamic effects
    t.on <- c(T.on)
    eff.v <- c(eff) ## a vector
    rm.pos1 <- which(is.na(eff.v))
    rm.pos2 <- which(is.na(t.on))
    eff.v.use1 <- eff.v
    t.on.use <- t.on
    n.on.use <- rep(1:N, each = TT)

    eff.equiv.v <- c(eff.equiv)

    if (NA %in% eff.v | NA %in% t.on) {
        eff.v.use1 <- eff.v[-c(rm.pos1, rm.pos2)]
        t.on.use <- t.on[-c(rm.pos1, rm.pos2)]
        n.on.use <- n.on.use[-c(rm.pos1, rm.pos2)]
        eff.equiv.v <- eff.equiv.v[-c(rm.pos1, rm.pos2)]
    }

    pre.pos <- which(t.on.use <= 0)
    eff.pre <- cbind(eff.v.use1[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
    colnames(eff.pre) <- c("eff", "period", "unit")

    eff.pre.equiv <- cbind(eff.equiv.v[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
    colnames(eff.pre.equiv) <- c("eff.equiv", "period", "unit")

    pre.sd <- tapply(eff.pre.equiv[, 1], eff.pre.equiv[, 2], sd)
    pre.sd <- cbind(pre.sd, sort(unique(eff.pre.equiv[, 2])), table(eff.pre.equiv[, 2]))
    colnames(pre.sd) <- c("sd", "period", "count")

    time.on <- sort(unique(t.on.use))
    att.on <- as.numeric(tapply(eff.v.use1, t.on.use, mean)) ## NA already removed
    count.on <- as.numeric(table(t.on.use))

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
    } else {
        att.on.sum.W <- att.on.W <- count.on.W <- time.on.W <- W.on.sum <- NULL
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
    }



    ## 5 carryover effect
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

        # if (!is.null(time.on.carry.seq)) {
        #    carry.att.med <- rep(NA, length(time.on.carry.seq))
        #    carry.att.med[which(time.on.carry.seq %in% carry.time)] <- carry.att
        #    carry.att <- carry.att.med

        # }
    }



    ## 6. switch-off effects
    eff.off <- eff.equiv <- off.sd <- NULL
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

        eff.off.equiv <- cbind(eff.equiv.v[off.pos], t.off.use[off.pos], n.on.use[off.pos])
        colnames(eff.off.equiv) <- c("off.equiv", "period", "unit")

        off.sd <- tapply(eff.off.equiv[, 1], eff.off.equiv[, 2], sd)
        off.sd <- cbind(off.sd, sort(unique(eff.off.equiv[, 2])), table(eff.off.equiv[, 2]))
        colnames(off.sd) <- c("sd", "period", "count")

        time.off <- sort(unique(t.off.use))
        att.off <- as.numeric(tapply(eff.v.use2, t.off.use, mean)) ## NA already removed
        count.off <- as.numeric(table(t.off.use))

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
        } else {
            W.off.sum <- att.off.sum.W <- att.off.W <- count.off.W <- time.off.W <- NULL
        }
    }

    ## 8. loess HTE by time
    D.missing <- D
    D.missing[which(D == 0)] <- NA
    eff.calendar <- apply(eff * D.missing, 1, mean, na.rm = TRUE)
    N.calendar <- apply(!is.na(eff * D.missing), 1, sum)
    T.calendar <- c(1:TT)
    if (sum(!is.na(eff.calendar)) > 1) {
        # loess fit
        loess.fit <- suppressWarnings(try(loess(eff.calendar ~ T.calendar, weights = N.calendar), silent = TRUE))
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

    ## 7. cohort effects
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

            suboutput <- list(
                att.on = att.on.sub,
                time.on = time.on.sub,
                count.on = count.on.sub
            )


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

                suboutput <- c(suboutput, list(
                    att.off = att.off.sub,
                    count.off = count.off.sub,
                    time.off = time.off.sub
                ))
            }
            group.output[[sub.group.name]] <- suboutput
        }
    }



    ## -------------------------------##
    ##           Storage
    ## -------------------------------##

    ## control group residuals
    out <- list(
        ## main results
        sigma2 = est.best$sigma2,
        sigma2.fect = est.fect$sigma2,
        T.on = T.on,
        Y.ct = est.best$fit,
        eff = eff,
        att.avg = att.avg,
        att.avg.unit = att.avg.unit,
        ## supporting
        force = force,
        T = TT,
        N = N,
        Ntr = Ntr,
        Nco = Nco,
        p = p,
        D = D,
        Y = Y,
        X = X,
        I = I,
        II = II,
        tr = tr,
        co = co,
        est = est.best,
        method = method,
        mu = est.best$mu,
        beta = beta,
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
        rmCV = rmCV,
        estCV = estCV,
        res = est.best$res
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
    }

    if (!is.null(T.on.carry)) {
        out <- c(out, list(carry.att = carry.att, carry.time = carry.time))
    }

    if (!is.null(balance.period)) {
        out <- c(out, list(balance.att = balance.att, balance.time = balance.time, balance.count = balance.count, balance.avg.att = att.avg.balance))
    }

    if (force == 1) {
        out <- c(out, list(
            alpha = est.best$alpha,
            alpha.tr = as.matrix(est.best$alpha[tr, ]),
            alpha.co = as.matrix(est.best$alpha[co, ])
        ))
    } else if (force == 2) {
        out <- c(out, list(xi = est.best$xi))
    } else if (force == 3) {
        out <- c(out, list(
            alpha = est.best$alpha, xi = est.best$xi,
            alpha.tr = as.matrix(est.best$alpha[tr, ]),
            alpha.co = as.matrix(est.best$alpha[co, ])
        ))
    }

    if (method == "ife") {
        out <- c(out, list(r.cv = r.cv, IC = IC, PC = PC))
        if (r.cv > 0) {
            out <- c(out, list(
                factor = as.matrix(est.best$factor),
                lambda = as.matrix(est.best$lambda),
                lambda.tr = as.matrix(est.best$lambda[tr, ]),
                lambda.co = as.matrix(est.best$lambda[co, ])
            ))
        }
    }

    if (method == "mc") {
        out <- c(out, list(
            lambda.cv = lambda.cv, lambda.seq = lambda,
            lambda.norm = lambda.cv / max(eigen.all),
            eigen.all = eigen.all
        ))
    }

    ## CV results
    if (!is.null(CV.out.ife)) {
        out <- c(out, list(CV.out.ife = CV.out.ife))
    }

    if (!is.null(CV.out.mc)) {
        out <- c(out, list(CV.out.mc = CV.out.mc))
    }

    if (!is.null(group)) {
        out <- c(out, list(
            group.att = group.att,
            group.output = group.output
        ))
    }

    return(out)
} ## cross-validation function ends
