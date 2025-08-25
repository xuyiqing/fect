fect_gsynth <- function(Y, # Outcome variable, (T*N) matrix
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
                        CV = TRUE,
                        criterion = "mspe",
                        r = 0, # r.end when CV==TRUE
                        r.end = 3,
                        binary = FALSE,
                        QR = FALSE,
                        force,
                        hasRevs = 1,
                        tol, # tolerance level
                        max.iteration = 1000,
                        boot = FALSE, # bootstrapped sample
                        placeboTest = 0,
                        placebo.period = NULL,
                        carryoverTest = 0,
                        carryover.period = NULL,
                        norm.para = NULL,
                        calendar.enp.seq = NULL,
                        time.on.seq = NULL,
                        time.off.seq = NULL,
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
    res.sd1 <- res.sd2 <- NULL
    ## unit id and time
    TT <- dim(Y)[1]
    N <- dim(Y)[2]
    id <- 1:N
    time <- 1:TT
    if (is.null(X) == FALSE) {
        p <- dim(X)[3]
    } else {
        p <- 0
        X <- array(0, dim = c(1, 1, 0))
    }

    ## replicate data
    YY <- Y
    YY[which(II == 0)] <- 0 ## reset to 0

    ## once treated, always treated
    ## careful unbalanced case: calculate treated units
    ## treat reversals as always treated
    D <- apply(D, 2, function(vec) {
        cumsum(vec)
    })
    D <- ifelse(D > 0, 1, 0)
    D.sum <- colSums(D)
    tr <- which(D.sum >= 1)
    Ntr <- length(tr)
    co <- which(D.sum == 0)
    Nco <- length(co)

    I.tr <- as.matrix(I[, tr]) ## maybe only 1 treated unit
    I.co <- I[, co]
    II.tr <- as.matrix(II[, tr])
    II.co <- II[, co]
    Y.tr <- as.matrix(Y[, tr])
    Y.co <- as.matrix(Y[, co])
    YY.tr <- as.matrix(YY[, tr])
    YY.co <- as.matrix(YY[, co])
    Ntr <- length(tr)
    Nco <- length(co)
    if (p == 0) {
        X.tr <- array(0, dim = c(TT, Ntr, 0))
        X.co <- array(0, dim = c(TT, Nco, 0))
    } else {
        X.tr <- array(NA, dim = c(TT, Ntr, p))
        X.co <- array(NA, dim = c(TT, Nco, p))
        for (j in 1:p) {
            X.tr[, , j] <- X[, tr, j]
            X.co[, , j] <- X[, co, j]
        }
    }

    if (is.null(W)) {
        W.use <- as.matrix(0)
    } 
    # else {
    #     stop("Gsynth doesn't suppoort weighted outcome model in this version.\n")
    # }

    if (!0 %in% I.tr) {
        ## a (TT*Ntr) matrix, time dimension: before treatment
        pre <- as.matrix(D[, tr] == 0 & II[, tr] == 1)
    } else {
        pre <- as.matrix(D[, tr] == 0 & I[, tr] == 1 & II[, tr] == 1)
    }

    T0 <- apply(pre, 2, sum)
    T0.min <- min(T0)
    pre.v <- as.vector(pre) ## vectorized "pre-treatment" indicator
    id.tr.pre.v <- rep(id, each = TT)[which(pre.v == 1)] ## vectorized pre-treatment grouping variable for the treated
    time.pre <- split(rep(time, Ntr)[which(pre.v == 1)], id.tr.pre.v) ## a list of pre-treatment periods
    sameT0 <- length(unique(T0)) == 1

    beta0 <- matrix(0, p, 1)
    # initial fit using Y.co
    if (0 %in% I.co) {
        data.ini <- matrix(NA, Nco * TT, (p + 3))
        data.ini[, 1] <- c(Y.co)
        data.ini[, 2] <- rep(1:Nco, each = TT)
        data.ini[, 3] <- rep(1:TT, Nco)
        if (p > 0) {
            for (i in 1:p) {
                data.ini[, (3 + i)] <- c(X.co[, , i])
            }
        }
        ## observed Y0 indicator:
        initialOut <- Y0.co <- beta0 <- FE0 <- xi0 <- factor0 <- NULL
        oci <- which(c(II.co) == 1)
        initialOut <- initialFit(data = data.ini, force = force, oci = oci)
        Y0.co <- initialOut$Y0
        beta0 <- initialOut$beta0
        if (p > 0 && sum(is.na(beta0)) > 0) {
            beta0[which(is.na(beta0))] <- 0
        }
    }


    validX <- 1 ## no multi-colinearity

    ## cross-validation only for gsynth
    if (CV == TRUE) {
        ## starting r
        if ((r > (T0.min - 1) & force %in% c(0, 2)) | (r > (T0.min - 2) & force %in% c(1, 3))) {
            message("r is too big compared with T0; reset to 0.")
            r <- 0
        }

        ## store all MSPE
        if (force %in% c(0, 2)) {
            r.max <- max(min((T0.min - 1), r.end), 0)
        } else {
            r.max <- max(min((T0.min - 2), r.end), 0)
        }

        if (r.max == 0) {
            r.cv <- 0
            message("Cross validation cannot be performed since available pre-treatment records of treated units are too few. So set r.cv = 0.")
            if (!0 %in% I.co) {
                est.co.best <- inter_fe(Y.co, X.co, 0, force = force, beta0 = beta0, tol, max.iteration)
            } else {
                est.co.best <- inter_fe_ub(Y.co, Y0.co, X.co, I.co, W = W.use, beta0, 0, force = force, tol, max.iteration)
            }
        } else {
            r.old <- r ## save the minimal number of factors

            message("Cross-validating ...", "\r")
            CV.out <- matrix(NA, (r.max - r.old + 1), 5)
            colnames(CV.out) <- c("r", "sigma2", "IC", "PC", "MSPE")
            CV.out[, "r"] <- c(r.old:r.max)
            CV.out[, "MSPE"] <- CV.out[, "PC"] <- 1e10
            r.pc <- est.co.pc.best <- NULL

            for (i in 1:dim(CV.out)[1]) {
                r <- CV.out[i, "r"]
                if (!0 %in% I.co) {
                    est.co <- inter_fe(
                        Y = Y.co, X = X.co, r,
                        force = force, beta0 = beta0, tol, max.iteration
                    )
                } else {
                    est.co <- inter_fe_ub(
                        Y = Y.co, Y0 = Y0.co, X = X.co, I = I.co, W = W.use,
                        beta0, r, force = force, tol, max.iteration
                    )
                }

                if (p > 0) {
                    na.pos <- is.nan(est.co$beta)
                    beta <- est.co$beta
                    beta[is.nan(est.co$beta)] <- 0
                }

                if (is.null(norm.para)) {
                    sigma2 <- est.co$sigma2
                    IC <- est.co$IC
                    PC <- est.co$PC
                } else {
                    sigma2 <- est.co$sigma2 * (norm.para[1]^2)
                    IC <- est.co$IC - log(est.co$sigma2) + log(sigma2)
                    PC <- est.co$PC * (norm.para[1]^2)
                }

                if (r != 0) {
                    F.hat <- as.matrix(est.co$factor)
                    if (force %in% c(1, 3)) {
                        F.hat <- cbind(F.hat, rep(1, TT))
                    }
                }
                U.tr <- Y.tr
                if (p > 0) {
                    for (j in 1:p) {
                        U.tr <- U.tr - X.tr[, , j] * beta[j]
                    }
                }
                U.tr <- U.tr - matrix(est.co$mu, TT, Ntr) ## grand mean
                if (force %in% c(2, 3)) {
                    U.tr <- U.tr - matrix(est.co$xi, TT, Ntr, byrow = FALSE)
                }

                if (0 %in% I.tr) {
                    U.tr[which(I.tr == 0)] <- 0
                }

                U.sav <- U.tr

                sum.e2 <- num.y <- 0
                for (lv in unique(unlist(time.pre))) {
                    U.tr <- U.sav
                    if (max(T0) == T0.min & (!0 %in% I.tr)) {
                        U.lv <- as.matrix(U.tr[setdiff(c(1:T0.min), lv), ]) ## setdiff : x
                    } else {
                        U.tr.pre.v <- as.vector(U.tr)[which(pre.v == 1)] ## pre-treatment residual in a vector
                        U.tr.pre <- split(U.tr.pre.v, id.tr.pre.v) ##  a list of pretreatment residuals
                        if (!0 %in% I.tr) {
                            U.lv <- lapply(U.tr.pre, function(vec) {
                                return(vec[-lv])
                            }) ## a list
                        } else {
                            ## U.tr.pre.sav <- U.tr.pre
                            for (i.tr in 1:Ntr) {
                                U.tmp <- U.tr.pre[[i.tr]]
                                U.tr.pre[[i.tr]] <- U.tmp[!time.pre[[i.tr]] == lv]
                            }
                            U.lv <- U.tr.pre
                        }
                    }

                    if (r == 0) {
                        if (force %in% c(1, 3)) { ## take out unit fixed effect
                            if (max(T0) == T0.min & (!0 %in% I.tr)) {
                                alpha.tr.lv <- colMeans(U.lv)
                                U.tr <- U.tr - matrix(alpha.tr.lv, TT, Ntr, byrow = TRUE)
                            } else {
                                alpha.tr.lv <- sapply(U.lv, mean)
                                U.tr <- U.tr - matrix(alpha.tr.lv, TT, Ntr, byrow = TRUE)
                            }
                        }
                        e <- U.tr[which(time == lv), ] ## that period
                    } else {
                        F.lv <- as.matrix(F.hat[which(time != lv), ])
                        if (max(T0) == T0.min & (!0 %in% I.tr)) {
                            F.lv.pre <- F.hat[setdiff(c(1:T0.min), lv), ]
                            lambda.lv <- try(
                                solve(t(F.lv.pre) %*% F.lv.pre) %*% t(F.lv.pre) %*% U.lv,
                                silent = TRUE
                            )
                            if ("try-error" %in% class(lambda.lv)) {
                                break
                            }
                        } else {
                            if (!0 %in% I.tr) {
                                lambda.lv <- try(as.matrix(sapply(U.lv, function(vec) {
                                    F.lv.pre <- as.matrix(F.lv[1:length(vec), ])
                                    l.lv.tr <- solve(t(F.lv.pre) %*% F.lv.pre) %*% t(F.lv.pre) %*% vec
                                    return(l.lv.tr)
                                })), silent = TRUE)
                                if ("try-error" %in% class(lambda.lv)) {
                                    break
                                } else {
                                    if ((r == 1) & (force %in% c(0, 2))) {
                                        lambda.lv <- t(lambda.lv)
                                    }
                                }
                            } else {
                                if (force %in% c(1, 3)) {
                                    lambda.lv <- matrix(NA, (r + 1), Ntr)
                                } else {
                                    lambda.lv <- matrix(NA, r, Ntr)
                                }
                                test <- try(
                                    for (i.tr in 1:Ntr) {
                                        F.lv.pre <- as.matrix(F.hat[setdiff(time.pre[[i.tr]], lv), ])
                                        lambda.lv[, i.tr] <- solve(t(F.lv.pre) %*% F.lv.pre) %*% t(F.lv.pre) %*% as.matrix(U.lv[[i.tr]])
                                    },
                                    silent = TRUE
                                )
                                if ("try-error" %in% class(test)) {
                                    break
                                }
                            }
                        }
                        lambda.lv <- t(lambda.lv) ## N*r
                        e <- U.tr[which(time == lv), ] - c(F.hat[which(time == lv), ] %*% t(lambda.lv))
                    }
                    if (sameT0 == FALSE | 0 %in% I.tr) {
                        e <- e[which(pre[which(time == lv), ] == TRUE)]
                    }
                    ## sum up
                    sum.e2 <- sum.e2 + t(e) %*% e
                    num.y <- num.y + length(e)
                } ## end of leave-one-out

                MSPE <- ifelse(num.y == 0, Inf, sum.e2 / num.y)
                if (!is.null(norm.para)) {
                    MSPE <- MSPE * (norm.para[1]^2)
                }

                if ((min(CV.out[, "MSPE"]) - MSPE) > tol * min(CV.out[, "MSPE"])) {
                    ## at least 5% improvement for MPSE
                    est.co.best <- est.co ## interFE result with the best r
                    r.cv <- r
                } else {
                    if (r == r.cv + 1) message("*")
                }

                if (PC < min(CV.out[, "PC"])) {
                    r.pc <- r
                    est.co.pc.best <- est.co
                }
                CV.out[i, 2:5] <- c(sigma2, IC, PC, MSPE)
                message("r = ", r, "; sigma2 = ",
                    sprintf("%.5f", sigma2), "; IC = ",
                    sprintf("%.5f", IC), "; PC = ",
                    sprintf("%.5f", PC), "; MSPE = ",
                    sprintf("%.5f", MSPE),
                    sep = ""
                )
            } ## end of while: search for r_star over

            MSPE.best <- min(CV.out[, "MSPE"])
            if (r > (T0.min - 1)) {
                message(" (r hits maximum)")
            }
            message("\n r* = ", r.cv, sep = "")
            message("\n")
        }
    } else {
        r.cv <- r
        r.min <- r.max <- r
    }

    est.co.fect <- NULL

    if (!0 %in% I.co) {
        est.co.best <- inter_fe(Y.co, X.co, r.cv,
            force = force, beta0 = beta0, tol, max.iteration
        )
    } else {
        est.co.best <- inter_fe_ub(YY.co, Y0.co, X.co, II.co, W = W.use, beta0, r.cv, force = force, tol, max.iteration)
    }

    if (boot == FALSE) {
        if (r.cv == 0) {
            est.co.fect <- est.co.best
        } else {
            if (!0 %in% I.co) {
                est.co.fect <- inter_fe(Y.co, X.co, 0,
                    force = force, beta0 = beta0, tol, max.iteration
                )
            } else {
                est.co.fect <- inter_fe_ub(YY.co, Y0.co, X.co, II.co, W = W.use, beta0, 0, force = force, tol, max.iteration)
            }
        }
    }
    validX <- est.co.best$validX
    validF <- ifelse(r.cv > 0, 1, 0)

    # get the counterfactual of X
    ## ## take out the effect of X
    U.tr.r0 <- U.tr <- Y.tr
    if (p > 0) {
        beta <- est.co.best$beta

        if (est.co.best$validX == 0) {
            beta <- matrix(0, p, 1)
        } else {
            beta <- est.co.best$beta
            beta[is.nan(est.co.best$beta)] <- 0
        }
        for (j in 1:p) {
            U.tr <- U.tr - X.tr[, , j] * beta[j]
        }

        if (boot == FALSE) {
            beta.r0 <- est.co.fect$beta
            if (est.co.fect$validX == 0) {
                beta.r0 <- matrix(0, p, 1)
            } else {
                beta.r0 <- est.co.fect$beta
                beta.r0[is.nan(est.co.fect$beta)] <- 0
            }
            for (j in 1:p) {
                U.tr.r0 <- U.tr.r0 - X.tr[, , j] * beta.r0[j]
            }
        }
    } else {
        beta <- NA
        beta.r0 <- NA
    }

    mu <- est.co.best$mu
    U.tr <- U.tr - matrix(mu, TT, Ntr) ## grand mean
    Y.fe.bar <- rep(mu, TT)

    if (boot == FALSE) {
        mu.r0 <- est.co.fect$mu
        U.tr.r0 <- U.tr.r0 - matrix(mu.r0, TT, Ntr)
        Y.fe.bar.r0 <- rep(mu.r0, TT)
    }

    if (force %in% c(2, 3)) {
        xi <- est.co.best$xi ## a (TT*1) matrix
        U.tr <- U.tr - matrix(c(xi), TT, Ntr, byrow = FALSE) ## will be adjusted at last
        Y.fe.bar <- Y.fe.bar + xi
        if (boot == FALSE) {
            xi.r0 <- est.co.fect$xi ## a (TT*1) matrix
            U.tr.r0 <- U.tr.r0 - matrix(c(xi.r0), TT, Ntr, byrow = FALSE)
            Y.fe.bar.r0 <- Y.fe.bar.r0 + xi.r0
        }
    }

    if (max(T0) == T0.min & (!0 %in% I.tr)) {
        U.tr.pre <- as.matrix(U.tr[1:T0.min, ])
        if (boot == FALSE) {
            U.tr.pre.r0 <- as.matrix(U.tr.r0[1:T0.min, ])
        }
    } else {
        ## not necessary to reset utr for ub data for pre.v doesn't include them
        U.tr.pre.v <- as.vector(U.tr)[which(pre.v == 1)] # pre-treatment residual in a vector
        U.tr.pre <- split(U.tr.pre.v, id.tr.pre.v) ##  a list of pretreatment residuals

        if (boot == FALSE) {
            U.tr.pre.v.r0 <- as.vector(U.tr.r0)[which(pre.v == 1)] # pre-treatment residual in a vector
            U.tr.pre.r0 <- split(U.tr.pre.v.r0, id.tr.pre.v) ##  a list of pretreatment residuals
        }
    }

    ## the error structure

    # for r=0
    if (force %in% c(1, 3)) { ## take out unit fixed effect
        if ((max(T0) == T0.min) & (!0 %in% I.tr)) {
            if (boot == FALSE) {
                alpha.tr.r0 <- as.matrix(colMeans(U.tr.pre.r0))
                U.tr.r0 <- U.tr.r0 - matrix(alpha.tr.r0, TT, Ntr, byrow = TRUE)
            }
        } else {
            if (boot == FALSE) {
                alpha.tr.r0 <- as.matrix(sapply(U.tr.pre.r0, mean))
                U.tr.r0 <- U.tr.r0 - matrix(alpha.tr.r0, TT, Ntr, byrow = TRUE)
            }
        }
    }
    if (boot == FALSE) {
        eff.r0 <- U.tr.r0
    }



    if (r.cv == 0) {
        if (force %in% c(1, 3)) { ## take out unit fixed effect
            if ((max(T0) == T0.min) & (!0 %in% I.tr)) {
                alpha.tr <- as.matrix(colMeans(U.tr.pre))
                U.tr <- U.tr - matrix(alpha.tr, TT, Ntr, byrow = TRUE)
            } else {
                alpha.tr <- as.matrix(sapply(U.tr.pre, mean))
                U.tr <- U.tr - matrix(alpha.tr, TT, Ntr, byrow = TRUE)
            }
        }
        eff <- U.tr
        lambda.tr <- NULL
        lambda.co <- NULL
    } else { ## Factors
        F.hat <- as.matrix(est.co.best$factor)
        if (force %in% c(1, 3)) {
            F.hat <- cbind(F.hat, rep(1, TT))
        }

        ## Lambda_tr (Ntr*r) or (Ntr*(r+1))
        if (max(T0) == T0.min & (!0 %in% I.tr)) {
            F.hat.pre <- F.hat[1:T0.min, ]
            lambda.tr <- try(solve(t(F.hat.pre) %*% F.hat.pre) %*% t(F.hat.pre) %*% U.tr.pre,
                silent = TRUE
            )
            if ("try-error" %in% class(lambda.tr)) {
                return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1)))
            }
        } else {
            if (!0 %in% I.tr) {
                lambda.tr <- try(as.matrix(sapply(U.tr.pre, function(vec) {
                    F.hat.pre <- as.matrix(F.hat[1:length(vec), ])
                    l.tr <- solve(t(F.hat.pre) %*% F.hat.pre) %*% t(F.hat.pre) %*% vec
                    return(l.tr) ## a vector of each individual lambdas
                })), silent = TRUE)
                if ("try-error" %in% class(lambda.tr)) {
                    return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1)))
                    ## stop("Error occurs. Please set a smaller value of factor number.")
                }
                if ((r.cv == 1) & (force %in% c(0, 2))) {
                    lambda.tr <- t(lambda.tr)
                }
            } else {
                if (force %in% c(1, 3)) {
                    lambda.tr <- matrix(NA, (r.cv + 1), Ntr)
                } else {
                    lambda.tr <- matrix(NA, r.cv, Ntr)
                }
                test <- try(
                    for (i.tr in 1:Ntr) {
                        F.hat.pre <- as.matrix(F.hat[time.pre[[i.tr]], ])
                        lambda.tr[, i.tr] <- solve(t(F.hat.pre) %*% F.hat.pre) %*% t(F.hat.pre) %*% as.matrix(U.tr.pre[[i.tr]])
                    },
                    silent = TRUE
                )
                if ("try-error" %in% class(test)) {
                    return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1), eff = matrix(NA, TT, Ntr)))
                    ## stop("Error occurs. Please set a smaller value of factor number.")
                }
            }
        }

        lambda.tr <- t(lambda.tr)
        eff <- U.tr - F.hat %*% t(lambda.tr)
        if (force %in% c(1, 3)) {
            alpha.tr <- as.matrix(lambda.tr[, (r.cv + 1), drop = FALSE])
            lambda.tr <- lambda.tr[, 1:r.cv, drop = FALSE]
        }
        if (boot == 0) {
            inv.tr <- try(
                ginv(t(as.matrix(lambda.tr))),
                silent = TRUE
            )

            if (!"try-error" %in% class(inv.tr)) {
                wgt.implied <- t(inv.tr %*% t(as.matrix(est.co.best$lambda)))
            }
        }
    } ## end of r!=0 case


    if (0 %in% I.tr) {
        eff[which(I.tr == 0)] <- 0 ## adjust
        if (boot == FALSE) {
            eff.r0[which(I.tr == 0)] <- 0
        }
    } ## missing data will be adjusted to NA finally

    ## -------------------------------##
    ## Summarize
    ## -------------------------------##

    ## counterfactuals for treated units
    Y.ct.tr <- as.matrix(Y.tr - eff)
    Y.ct.co <- Y.co - est.co.best$residuals
    # print(Y.ct.co)
    Y.ct <- Y
    Y.ct[, tr] <- Y.ct.tr
    Y.ct[, co] <- Y.ct.co

    if (boot == FALSE) {
        Y.ct.tr.r0 <- as.matrix(Y.tr - eff.r0)
        Y.ct.co.r0 <- Y.co - est.co.fect$residuals
        Y.ct.r0 <- Y
        Y.ct.r0[, tr] <- Y.ct.tr.r0
        Y.ct.r0[, co] <- Y.ct.co.r0
    }

    ## we first adjustment for normalization
    if (!is.null(norm.para)) {
        Y <- Y * norm.para[1]
        ## variance of the error term
        sigma2 <- est.co.best$sigma2 <- est.co.best$sigma2 * (norm.para[1]^2)
        IC <- est.co.best$IC <- est.co.best$IC - log(est.co.best$sigma2) + log(sigma2)
        PC <- est.co.best$PC <- est.co.best$PC * (norm.para[1]^2)

        ## output of estimates
        mu <- est.co.best$mu <- est.co.best$mu * norm.para[1]
        if (r.cv > 0) {
            est.co.best$lambda <- est.co.best$lambda * norm.para[1]
            lambda.tr <- lambda.tr * norm.para[1]
            est.co.best$VNT <- est.co.best$VNT * norm.para[1]
        }
        if (force %in% c(1, 3)) {
            est.co.best$alpha <- est.co.best$alpha * norm.para[1]
            alpha.tr <- alpha.tr * norm.para[1]
        }
        if (force %in% c(2, 3)) {
            xi <- est.co.best$xi <- est.co.best$xi * norm.para[1]
        }
        est.co.best$residuals <- est.co.best$residuals * norm.para[1]
        est.co.best$fit <- est.co.best$fit * norm.para[1]
        if (boot == FALSE) {
            est.co.fect$fit <- est.co.fect$fit * norm.para[1]
        }
        est.co.fect$sigma2 <- est.co.fect$sigma2 * norm.para[1]

        # estimated counterfactual
        Y.tr <- Y.tr * norm.para[1]
        Y.co <- Y.co * norm.para[1]

        Y.ct <- Y.ct * norm.para[1]
        Y.ct.tr <- Y.ct.tr * norm.para[1]
        Y.ct.co <- Y.ct.co * norm.para[1]
        eff <- eff * norm.para[1]

        if (boot == FALSE) {
            Y.ct.r0 <- Y.ct.r0 * norm.para[1]
            Y.ct.tr.r0 <- Y.ct.tr.r0 * norm.para[1]
            Y.ct.co.r0 <- Y.ct.co.r0 * norm.para[1]
            eff.r0 <- eff.r0 * norm.para[1]
        }
    }

    ## 0. relevant parameters
    IC <- est.co.best$IC
    sigma2 <- est.co.best$sigma2
    PC <- est.co.best$PC
    loglikelihood <- NULL

    if (p > 0) {
        na.pos <- is.nan(est.co.best$beta)
        beta <- est.co.best$beta
        if (sum(na.pos) > 0) {
            beta[na.pos] <- NA
        }
    } else {
        beta <- NA
    }

    ## 1. estimated att and counterfactuals
    if (boot == FALSE) {
        Y.ct.equiv <- Y.ct.r0
    } else {
        Y.ct.equiv <- NULL
    }

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
    marginal <- NULL

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
        eff.equiv <- Y - Y.ct.equiv
        if (0 %in% I) {
            eff.equiv[which(I == 0)] <- NA
        }
        complete.index <- which(!is.na(eff.equiv))
        equiv.att.avg <- sum(eff.equiv[complete.index] * D[complete.index]) / (sum(D[complete.index]))
    }

    ## 2. rmse for treated units' observations under control
    if (binary == 0) {
        tr <- which(apply(D, 2, sum) > 0)
        tr.co <- which((as.matrix(1 - D[, tr]) * as.matrix(II[, tr])) == 1)
        eff.tr <- as.matrix(eff[, tr])
        v.eff.tr <- eff.tr[tr.co]
        rmse <- sqrt(mean(v.eff.tr^2, na.rm = TRUE))
    }

    ## 3. unbalanced output
    Y.ct.full <- Y.ct
    res.full <- Y - Y.ct
    if (0 %in% I) {
        eff[which(I == 0)] <- NA
        Y.ct[which(I == 0)] <- NA
    }
    if (binary == FALSE) {
        res.full[which(II == 0)] <- NA
    }


    ## 4. dynamic effects
    t.on <- c(T.on)
    eff.v <- c(eff) ## a vector
    eff.equiv.v <- NULL
    if (binary == FALSE && boot == FALSE) {
        eff.equiv.v <- c(eff.equiv)
    }

    rm.pos1 <- which(is.na(eff.v))
    rm.pos2 <- which(is.na(t.on))

    eff.v.use1 <- eff.v
    t.on.use <- t.on
    n.on.use <- rep(1:N, each = TT)

    if (NA %in% eff.v | NA %in% t.on) {
        eff.v.use1 <- eff.v[-c(rm.pos1, rm.pos2)]
        t.on.use <- t.on[-c(rm.pos1, rm.pos2)]
        n.on.use <- n.on.use[-c(rm.pos1, rm.pos2)]
        if (binary == FALSE && boot == FALSE) {
            eff.equiv.v <- eff.equiv.v[-c(rm.pos1, rm.pos2)]
        }
    }

    pre.pos <- which(t.on.use <= 0)
    eff.pre <- cbind(eff.v.use1[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
    colnames(eff.pre) <- c("eff", "period", "unit")

    pre.sd <- eff.pre.equiv <- NULL
    if (binary == FALSE && boot == FALSE) {
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

    ## 6. switch-off effects
    eff.off.equiv <- off.sd <- eff.off <- NULL
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

        if (binary == FALSE && boot == FALSE) {
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
    if (!is.null(carryover.period) && carryoverTest == 1 && hasRevs == 1) {
        ## construct att.carryover
        ## eff is derived from eff.v
        ## period and Num.Units are derived from T.off
        if (length(carryover.period) == 1) {
            carryover.pos <- which(time.off == carryover.period)
            att.carryover <- att.off[carryover.pos]
        } else {
            carryover.pos <- which(time.off >= carryover.period[1] & time.off <= carryover.period[2])
            att.carryover <- sum(att.off[carryover.pos] * count.off[carryover.pos]) / sum(count.off[carryover.pos])
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
            if (length(att.on.sub) == 0) {
                att.on.sub <- NULL
            }
            if (length(time.on.sub) == 0) {
                time.on.sub <- NULL
            }
            if (length(count.on.sub) == 0) {
                count.on.sub <- NULL
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
                if (length(att.placebo.sub) == 0) {
                    att.placebo.sub <- NULL
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
                if (length(att.off.sub) == 0) {
                    att.off.sub <- NULL
                }
                if (length(time.off.sub) == 0) {
                    time.off.sub <- NULL
                }
                if (length(count.off.sub) == 0) {
                    count.off.sub <- NULL
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
                    if (length(att.carryover.sub) == 0) {
                        att.carryover.sub <- NULL
                    }
                    suboutput <- c(suboutput, list(att.carryover = att.carryover.sub))
                }
            }
            group.output[[sub.group.name]] <- suboutput
        }
    }

    method <- ifelse(r.cv > 0, "gsynth", "fe")

    ## -------------------------------##
    ##            Storage            ##
    ## -------------------------------##
    out <- list(
        ## main results
        method = method,
        Y.ct = Y.ct,
        Y.ct.full = Y.ct.full,
        Y = Y,
        D = D,
        tr = tr,
        co = co,
        eff = eff,
        eff.tr = eff[, tr],
        I = I,
        II = II,
        att.avg = att.avg,
        att.avg.unit = att.avg.unit,
        ## supporting
        force = force,
        T = TT,
        N = N,
        Ntr = Ntr,
        Nco = Nco,
        p = p,
        r.cv = r.cv,
        IC = IC,
        beta = beta,
        est = est.co.best,
        mu = est.co.best$mu,
        niter = est.co.best$niter,
        validX = validX,
        validF = validF,
        time = time.on,
        att = att.on,
        count = count.on,
        eff.calendar = eff.calendar,
        N.calendar = N.calendar,
        eff.calendar.fit = eff.calendar.fit,
        eff.pre = eff.pre,
        eff.pre.equiv = eff.pre.equiv,
        pre.sd = pre.sd
    )


    out <- c(out, list(
        PC = PC,
        sigma2 = sigma2,
        sigma2.fect = est.co.fect$sigma2,
        res = est.co.best$residuals,
        res.full = res.full,
        rmse = rmse
    ))


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
    if (r.cv > 0) {
        lambda.co <- as.matrix(est.co.best$lambda)
        rownames(lambda.co) <- co
        rownames(lambda.tr) <- tr
        out <- c(out, list(
            factor = as.matrix(est.co.best$factor),
            lambda.co = as.matrix(lambda.co),
            lambda.tr = as.matrix(lambda.tr)
        ))
        if (boot == 0) {
            if (!"try-error" %in% class(inv.tr)) {
                out <- c(out, list(wgt.implied = wgt.implied))
            }
        }
    }

    if (force == 1) {
        out <- c(out, list(
            alpha.tr = alpha.tr,
            alpha.co = est.co.best$alpha
        ))
    } else if (force == 2) {
        out <- c(out, list(xi = est.co.best$xi))
    } else if (force == 3) {
        out <- c(out, list(
            alpha.tr = alpha.tr,
            alpha.co = est.co.best$alpha,
            xi = est.co.best$xi
        ))
    }

    if (!is.null(placebo.period) && placeboTest == 1) {
        out <- c(out, list(att.placebo = att.placebo))
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
    }

    if (!is.null(balance.period)) {
        out <- c(out, list(balance.att = balance.att, balance.time = balance.time, balance.count = balance.count, balance.avg.att = att.avg.balance))
        if (!is.null(placebo.period) && placeboTest == 1) {
            out <- c(out, list(balance.att.placebo = balance.att.placebo))
        }
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
}
