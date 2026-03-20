#####################################
## goodness of fit test: wild
#####################################

fect_test <- function(
    out, # from fect
    Y,
    X,
    D, ## input
    I,
    II,
    T.on,
    T.off = NULL,
    method = "ife",
    degree = 2,
    knots = NULL, 
    cl = NULL,
    r = 0,
    lambda = Inf,
    force,
    hasRevs = 0,
    tol,
    norm.para,
    pre.period = NULL,
    off.period = NULL,
    nboots,
    parallel = TRUE,
    cores = NULL) {
    

    x <- y <- na.pos <- NULL
    TT <- dim(Y)[1]
    N <- dim(Y)[2]
    if (is.null(X) == FALSE) {
        p <- dim(X)[3]
    } else {
        p <- 0
    }
    
    ## output
    eff <- out$eff
    eff[which(is.na(eff))] <- 0

    ## fitted values Y0
    Y.f <- out$est$fit

    if (is.null(pre.period)) {
        pre.period <- c(min(c(T.on), na.rm = TRUE), 0)
    }

    if (is.null(off.period) & hasRevs == TRUE) {
        off.period <- c(0, max(c(T.off), na.rm = TRUE))
    }

    pre.pos <- which(T.on <= pre.period[2] & T.on >= pre.period[1] & I == 1)

    ## message("\n OK3 \n")

    ## demean pre-treatment tr eff under H0 for wild bootstrap
    #eff[pre.pos] <- eff[pre.pos] - mean(eff[pre.pos])

    ## observed F
    eff.pre <- as.data.frame(out$eff.pre)
    eff.pre <- eff.pre[which(eff.pre[,"period"] >= pre.period[1] & eff.pre[,"period"] <= pre.period[2]),]

    lm.fit <- lm(eff ~ factor(period), data = eff.pre)
    #res <- lm.fit$residuals
    #f <- ((sum(eff.pre[,"eff"]^2) - sum(res^2))/length(unique(eff.pre[,"period"])))/(sum(res^2)/(dim(eff.pre)[1] - length(unique(eff.pre[,"period"]))))
    f <- summary(lm.fit)$f[1]
    names(f) <- NULL

    ## message("\n OK4 \n")

    f2 <- NULL
    if (hasRevs) {
        ## switch-off 
        eff.off <- as.data.frame(out$eff.off)
        eff.off <- eff.off[which(eff.off[,"period"] >= off.period[1] & eff.off[,"period"] <= off.period[2]),]
        lm.fit2 <- lm(eff ~ factor(period), data = eff.off)

        f2 <- summary(lm.fit2)$f[1]
        names(f2) <- NULL
    }

    ## bootstrapped F under H0
    f.boot <- rep(NA, nboots)
    f2.boot <- NULL
    if (hasRevs) {
        f2.boot <- rep(NA, nboots)
    }

    ## message("\n OK1 \n")

    #if (method == "ife") {
        one.nonpara <- function() {
                
            if (is.null(cl)) {
                res.p <- matrix(sample(c(-1, 1), N*TT, replace = TRUE), TT, N)
            } else {
                res.p <- matrix(rep(sample(c(-1, 1), N, replace = TRUE), each = TT), TT, N)
            }
            Y.boot <- Y.f + res.p * eff
            Y.boot[which(I == 0)] <- 0
            if (method %in% c("ife", "fe")) {
                boot <- try(fect_fe(Y = Y.boot, X = X, D = D, I = I, II = II, 
                                T.on = T.on, T.off = NULL,
                                r.cv = r, force = force, hasRevs = 0, 
                                tol = tol, boot = 1,
                                norm.para = norm.para, 
                                placebo.period = NULL, placeboTest = 0), silent = TRUE)
            }
            else if (method == "mc") {
                boot <- try(fect_mc(Y = Y.boot, X = X, D = D, 
                                I = I, II = II, hasF = out$validF,
                                T.on = T.on, T.off = NULL, 
                                lambda.cv = lambda, force = force, hasRevs = 0, 
                                tol = tol, boot = 1,
                                norm.para = norm.para,
                                placebo.period = NULL, placeboTest = 0), silent = TRUE)
            }
            else if (method %in% c("polynomial", "bspline")) {
                boot <- try(fect_polynomial(Y = Y.boot, X = X, D = D, I = I, II = II, 
                                T.on = T.on, T.off = NULL, method = method, 
                                degree = degree, knots = knots, 
                                force = force, hasRevs = 0, tol = tol, boot = 1,
                                norm.para = norm.para, 
                                placebo.period = NULL, placeboTest = 0), silent = TRUE)
            }
            
            if ('try-error' %in% class(boot)) {
                ## message("NA")
                return(list(f.boot = NA, f.boot.2 = NA))
            } else {
                ## message(boot$niter)
                ## message("\n")
                data <- as.data.frame(boot$eff.pre)
                #colnames(data) <- c("period", "eff")
                data <- data[which(data[,"period"] <= pre.period[2] & data[,"period"] >= pre.period[1]),]
                lm.fit.boot <- lm(eff~factor(period), data = data)
                #f.boot <- ((sum(data[,"eff"]^2) - sum((lm.fit.boot$residuals)^2))/length(unique(data[,"period"])))/(sum((lm.fit.boot$residuals)^2)/(dim(data[1]-unique(data[,"period"]))))
                f.boot <- summary(lm.fit.boot)$f[1]
                ## if (r.cv == 0) {

                f.boot.2 <- NULL
                if (hasRevs) {
                    ## switch-off 
                    data2 <- as.data.frame(out$eff.off)
                    data2 <- data2[which(data2[,"period"] >= off.period[1] & data2[,"period"] <= off.period[2]),]
                    lm.fit.boot2 <- lm(eff ~ factor(period), data = data2)

                    f.boot.2 <- summary(lm.fit.boot2)$f[1]
                }
                if (hasRevs) {
                    return(list(f.boot = f.boot, f.boot.2 = f.boot.2))
                } else {
                    return(list(f.boot = f.boot))
                }
                
                ## } else {
                ##     return(list(f.boot = f.boot, factor.boot = boot$factor))
                ## } 
            }    
            
        } 
            
    #} else { ## mc
    #    one.nonpara <- function() {
            
    #        if (is.null(cl)) {
    #            res.p <- matrix(sample(c(-1, 1), N*TT, replace = TRUE), TT, N)
    #        } else {
    #            res.p <- matrix(rep(sample(c(-1, 1), N, replace = TRUE), each = TT), TT, N)
    #        }
    #        Y.boot <- Y.f + res.p * eff
    #        Y.boot[which(I == 0)] <- 0
    #        boot <- try(fect.mc(Y = Y.boot, X = X, D = D, 
    #                            I = I, II = II, hasF = out$validF,
    #                            T.on = T.on, T.off = NULL, 
    #                            lambda.cv = lambda, force = force, hasRevs = 0, 
    #                            tol = tol, boot = 1,
    #                            norm.para = norm.para,
    #                            placebo.period = NULL, placeboTest = 0), silent = TRUE)

    #        if ('try-error' %in% class(boot)) {
                ## message("NA")
    #            return(list(f.boot = NA, f.boot.2 = NA))
    #        } else {
                ## message(boot$niter)
                ## message("\n")
    #            data <- as.data.frame(boot$eff.pre)
                #colnames(data) <- c("period", "eff")
    #            data <- data[which(data[,"period"] <= pre.period[2] & data[,"period"] >= pre.period[1]),]
    #            lm.fit.boot <- lm(eff~factor(period), data = data)
                #f.boot <- ((sum(data[,"eff"]^2) - sum((lm.fit.boot$residuals)^2))/length(unique(data[,"period"])))/(sum((lm.fit.boot$residuals)^2)/(dim(data[1]-unique(data[,"period"]))))
    #            f.boot <- summary(lm.fit.boot)$f[1]
                ## if (r.cv == 0) {
                
    #            f.boot.2 <- NULL
    #            if (hasRevs) {
                    ## switch-off 
    #                data2 <- as.data.frame(out$eff.off)
    #                data2 <- data2[which(data2[,"period"] >= off.period[1] & data2[,"period"] <= off.period[2]),]
    #                lm.fit.boot2 <- lm(eff ~ factor(period), data = data2)

    #                f.boot.2 <- summary(lm.fit.boot2)$f[1]
    #            }
    #            if (hasRevs) {
    #                return(list(f.boot = f.boot, f.boot.2 = f.boot.2))
    #            } else {
    #                return(list(f.boot = f.boot))
    #            }

                ## } else {
                ##     return(list(f.boot = f.boot, factor.boot = boot$factor))
                ## } 
    #        }
    #    } 
    #}
    
    ## computing
    if (parallel == TRUE) { 
        boot.out <- foreach(j=1:nboots, 
                            .inorder = FALSE,
                            .export = c("fect_fe", "fect_mc", "fect_polynomial", "get_term"),
                            .packages = c("fect")
                            ) %dopar% {
                                return(one.nonpara())
                            }

        for (j in 1:nboots) { 
            f.boot[j] <- boot.out[[j]]$f.boot
        } 
    } else {
        for (j in 1:nboots) { 
            one.boot <- one.nonpara()
            f.boot[j] <- one.boot$f.boot
            ## report progress
            if (j%%100 == 0)  {
                message(".")   
            }  
        }  
    } 

    ## message("\n OK2 \n")

    ## end of bootstrapping
    if (sum(is.na(f.boot)) > 0) {
        f.boot <- f.boot[-which(is.na(f.boot))]
    }
    message(length(f.boot), " runs\n", sep = "")
    
    f.q <- quantile(f.boot, probs = c(0.025, 0.975))
    f.p <- sum(f.boot > f)/length(f.boot)
    wald.on <- list(stat = f, 
                    p = f.p,
                    quantile_lower = f.q[1], 
                    quantile_upper = f.q[2],
                    boot = f.boot)


    wald.off <- NULL
    if (hasRevs) {
        if (sum(is.na(f2.boot)) > 0) {
            f2.boot <- f2.boot[-which(is.na(f2.boot))]
        }
        f.q2 <- quantile(f2.boot, probs = c(0.025, 0.975))
        f.p2 <- sum(f2.boot > f2)/length(f2.boot)
        wald.off <- list(stat = f2, 
                         p = f.p2,
                         quantile_lower = f.q2[1], 
                         quantile_upper = f.q2[2],
                         boot = f2.boot)
    }

    if (hasRevs) {
        return(list(wald.on = wald.on, wald.off = wald.off))
    } else {
        return(list(wald.on = wald.on))
    }
    
} ## end of test

#####################################
## goodness of fit test: non para 
#####################################

#fect.test2 <- function(Y,
#                       X,
#                       D, ## input
#                       I,
#                       II,
#                       T.on,
#                       T.off = NULL,
#                       method = "ife",
#                       r = 0,
#                       lambda = Inf,
#                       force,
#                       hasRevs = 0,
#                       tol,
#                       norm.para,
#                       pre.period = NULL,
#                       nboots,
#                       parallel = TRUE,
#                       cores = NULL)  {
    
    
#    x <- y <- na.pos <- NULL
#    TT <- dim(Y)[1]
#    N <- dim(Y)[2]
#    if (is.null(X) == FALSE) {
#        p <- dim(X)[3]
#    } else {
#        p <- 0
#    }

#    if (hasRevs == 1) {
        ## D.fake : check reversals
#        D.fake <- apply(D, 2, function(vec){cumsum(vec)})
#        D.fake <- ifelse(D.fake > 0, 1, 0)
#        D.fake[which(I == 0)] <- 0

#        rev <- which(apply(D.fake == D, 2, sum) != TT)
#        co <- which(apply(D, 2, sum) == 0)
#        tr.all <- which(apply(D, 2, sum) > 0)
#        tr <- tr.all[which(!tr.all %in% rev)]

#        Nrev <- length(rev)
#        Ntr <- length(tr)
#        Nco <- length(co)
#    } else {
        ## treatement indicator
#        tr <- which(apply(D, 2, sum) > 0)
#        co <- which(apply(D, 2, sum) == 0)

#        Ntr <- length(tr)
#        Nco <- length(co)
#    }
    
    ## estimation
#    if (method == "ife") {
#        out <- fect.fe(Y = Y, X = X, D = D, I = I, II = II, 
#                       T.on = T.on, T.off = NULL,
#                       r.cv = r, force = force, hasRevs = 0, 
#                       tol = tol, boot = 1,
#                       norm.para = norm.para, 
#                       placebo.period = NULL, placeboTest = 0)
#    } else {
#        out <- fect.mc(Y = Y, X = X, D = D, I = I, II = II,
#                       T.on = T.on, T.off = NULL, 
#                       lambda.cv = lambda, force = force, hasRevs = 0, 
#                       tol = tol, boot = 1,
#                       norm.para = norm.para,
#                       placebo.period = NULL, placeboTest = 0)
#    }

    ## output
#    eff <- out$eff
#    eff[which(is.na(eff))] <- 0

    ## fitted values Y0
#    Y.f <- out$est$fit

#    if (is.null(pre.period)) {
#        pre.period <- c(min(c(T.on), na.rm = TRUE), 0)
#    }

#    pre.pos <- which(T.on <= pre.period[2] & T.on >= pre.period[1] & I == 1)

    ## demean pre-treatment tr eff under H0 for wild bootstrap
    #eff[pre.pos] <- eff[pre.pos] - mean(eff[pre.pos])

    ## observed F
#    eff.pre <- as.data.frame(out$eff.pre)
#    eff.pre <- eff.pre[which(eff.pre[,"period"] >= pre.period[1] & eff.pre[,"period"] <= pre.period[2]),]

#    eff.pre[,"period"] <- as.factor(eff.pre[,"period"])
#    lm.fit <- lm(eff ~ period, data = eff.pre)
    #res <- lm.fit$residuals
    #f <- ((sum(eff.pre[,"eff"]^2) - sum(res^2))/length(unique(eff.pre[,"period"])))/(sum(res^2)/(dim(eff.pre)[1] - length(unique(eff.pre[,"period"]))))
#    f <- summary(lm.fit)$f[1]


    ## bootstrapped F under H0
#    f.boot <- rep(NA, nboots)

#    message("\rBootstrapping ...\n")

#    if (method == "ife") {
#        one.nonpara <- function() {

#            if (hasRevs == 0) {
#                if (Nco > 0) {
#                    repeat{
#                        fake.co <- sample(co, Nco, replace=TRUE)
#                        fake.tr <- sample(tr, Ntr, replace=TRUE)
#                        boot.id <- c(fake.tr, fake.co)
#                        if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
#                            break
#                        }
#                    }
#                } else {
#                    repeat{
#                        boot.id <- sample(tr, Ntr, replace=TRUE)
#                        if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
#                            break
#                        }
#                    }
#                }
#            } else {
#                if (Ntr > 0) {
#                    if (Nco > 0) {
#                        repeat{
#                            fake.co <- sample(co, Nco, replace=TRUE)
#                            fake.tr <- sample(tr, Ntr, replace=TRUE)
#                            fake.rev <- sample(rev, Nrev, replace=TRUE)
#                            boot.id <- c(fake.rev, fake.tr, fake.co)
#                            if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
#                                break
#                            }
#                        }
#                    } else {
#                        repeat{
#                            fake.tr <- sample(tr, Ntr, replace=TRUE)
#                            fake.rev <- sample(rev, Nrev, replace=TRUE)
#                            boot.id <- c(fake.rev, fake.tr)
#                            if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
#                                break
#                            }
#                        }
#                    }
#                } else {
#                    if (Nco > 0) {
#                        repeat{
#                            fake.co <- sample(co, Nco, replace=TRUE)
#                            fake.rev <- sample(rev, Nrev, replace=TRUE)
#                            boot.id <- c(fake.rev, fake.co)
#                            if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
#                                break
#                            }
#                        }
#                    } else {
#                        repeat{
#                            boot.id <- sample(rev, Nrev, replace=TRUE)
#                            if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
#                                break
#                            }
#                        }
#                    }
#                }
#            }
                
#            X.boot <- X[,boot.id,,drop=FALSE]

                
#            boot <- try(fect.fe(Y = Y[, boot.id], X = X.boot, D = D[, boot.id],
#                                I = I[, boot.id], II = II[, boot.id], 
#                                T.on = T.on[, boot.id], T.off = NULL, 
#                                r.cv = r, force = force,
#                                hasRevs = 0, tol = tol, boot = 1,
#                                norm.para = norm.para,
#                                time.on.seq = NULL, time.off.seq = NULL,
#                                placebo.period = NULL, 
#                                placeboTest = 0), silent = TRUE)
#            if ('try-error' %in% class(boot)) {
#                return(NA)
#            } else {
#                data <- as.data.frame(boot$eff.pre)
                #colnames(data) <- c("period", "eff")
#                data <- data[which(data[,"period"] <= pre.period[2] & data[,"period"] >= pre.period[1]),]
#                data[,"period"] <- as.factor(data[,"period"])

#                eff.fit <- predict(lm.fit, data)
                ## remove pre-treatment effect
#                data[,"eff0"] <- data[,"eff"] - eff.fit

#                lm.fit.boot <- lm(eff0~period, data = data)
                #f.boot <- ((sum(data[,"eff"]^2) - sum((lm.fit.boot$residuals)^2))/length(unique(data[,"period"])))/(sum((lm.fit.boot$residuals)^2)/(dim(data[1]-unique(data[,"period"]))))
#                f.boot <- summary(lm.fit.boot)$f[1]
                ## if (r.cv == 0) {
#                return(list(f.boot = f.boot))
                ## } else {
                ##     return(list(f.boot = f.boot, factor.boot = boot$factor))
                ## } 
#            }
#        } 
            
#    } else { ## mc
#        one.nonpara <- function() {

#            if (hasRevs == 0) {
#                if (Nco > 0) {
#                    repeat{
#                        fake.co <- sample(co, Nco, replace=TRUE)
#                        fake.tr <- sample(tr, Ntr, replace=TRUE)
#                        boot.id <- c(fake.tr, fake.co)
#                        if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
#                            break
#                        }
#                    }
#                } else {
#                    repeat{
#                        boot.id <- sample(tr, Ntr, replace=TRUE)
#                        if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
#                            break
#                        }
#                    }
#                }
#            } else {
#                if (Ntr > 0) {
#                    if (Nco > 0) {
#                        repeat{
#                            fake.co <- sample(co, Nco, replace=TRUE)
#                            fake.tr <- sample(tr, Ntr, replace=TRUE)
#                            fake.rev <- sample(rev, Nrev, replace=TRUE)
#                            boot.id <- c(fake.rev, fake.tr, fake.co)
#                            if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
#                                break
#                            }
#                        }
#                    } else {
#                        repeat{
#                            fake.tr <- sample(tr, Ntr, replace=TRUE)
#                            fake.rev <- sample(rev, Nrev, replace=TRUE)
#                            boot.id <- c(fake.rev, fake.tr)
#                            if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
#                                break
#                            }
#                        }
#                    }
#                } else {
#                    if (Nco > 0) {
#                        repeat{
#                            fake.co <- sample(co, Nco, replace=TRUE)
#                            fake.rev <- sample(rev, Nrev, replace=TRUE)
#                            boot.id <- c(fake.rev, fake.co)
#                            if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
#                                break
#                            }
#                        }
#                    } else {
#                        repeat{
#                            boot.id <- sample(rev, Nrev, replace=TRUE)
#                            if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
#                                break
#                            }
#                        }
#                    }
#                }
#            }
                
#            X.boot <- X[,boot.id,,drop=FALSE]
            
#            boot <- try(fect.mc(Y = Y[,boot.id], X = X.boot, D = D[,boot.id],
#                                I = I[,boot.id], II = II[,boot.id],
#                                T.on = T.on[,boot.id], T.off = NULL, 
#                                lambda.cv = lambda, force = force, 
#                                hasF = out$validF, hasRevs = 0, 
#                                tol = tol, boot = 1,
#                                norm.para = norm.para,
#                                time.on.seq = NULL, time.off.seq = NULL,
#                                placebo.period = NULL, 
#                                placeboTest = 0), silent = TRUE)
            
#            if ('try-error' %in% class(boot)) {
#                return(NA)
#            } else {
#                data <- as.data.frame(boot$eff.pre)
                #colnames(data) <- c("period", "eff")
#                data <- data[which(data[,"period"] <= pre.period[2] & data[,"period"] >= pre.period[1]),]
#                data[,"period"] <- as.factor(data[,"period"])

#                eff.fit <- predict(lm.fit, data)
#                data[,"eff0"] <- data[,"eff"] - eff.fit

#                lm.fit.boot <- lm(eff0~period, data = data)
                #f.boot <- ((sum(data[,"eff"]^2) - sum((lm.fit.boot$residuals)^2))/length(unique(data[,"period"])))/(sum((lm.fit.boot$residuals)^2)/(dim(data[1]-unique(data[,"period"]))))
#                f.boot <- summary(lm.fit.boot)$f[1]
                ## if (r.cv == 0) {
#                return(list(f.boot = f.boot))
                ## } else {
                ##     return(list(f.boot = f.boot, factor.boot = boot$factor))
                ## }
#            }
            
#        } 
#    }

    ## computing
#    if (parallel == TRUE) { 
#        boot.out <- foreach(j=1:nboots, 
#                            .inorder = FALSE,
#                            .export = c("fect.fe", "fect.mc", "get_term"),
#                            .packages = c("fect")
#                            ) %dopar% {
#                                return(one.nonpara())
#                            }

#        for (j in 1:nboots) { 
#            f.boot[j] <- boot.out[[j]]$f.boot
#        } 
#    } else {
#        for (j in 1:nboots) { 
#            one.boot <- one.nonpara()
#            f.boot[j] <- one.boot$f.boot
#            ## report progress
#            if (j%%100 == 0)  {
#                message(".")   
#            }  
#        }  
#    } 
    ## end of bootstrapping
#    message("\r")

    
#    f.q <- quantile(f.boot, probs = c(0.025, 0.975))
#    f.p <- sum(f.boot > f)/nboots

#    title <- "Empirical distribution of F under H0"
    ## plot
#    f.data <- cbind.data.frame(f.boot = f.boot)
#    p <- ggplot(f.data, aes(f.boot)) + geom_density() + geom_vline(xintercept = f, colour="red",size = 0.5)

#    d <- ggplot_build(p)$data[[1]]

#    if (f < max(d[,"x"])) {
#        p <- p + geom_area(data = subset(d, x > f), aes(x=x, y=y), fill="red", alpha = 0.2)
#    }

#    p <- p + annotate("text", x = 0.9*max(d[,"x"]), y = 0.9*max(d[,"y"]), label = paste("P value: ", f.p, sep=""))

#    p <- p + ggtitle(title) +  theme(plot.title = element_text(size=20,
#                                                               hjust = 0.5,
#                                                               margin = margin(10, 0, 10, 0)))

    ## suppressWarnings(print(p))
  
    ## store results
#    Ftest <- matrix(NA, 1, 4)
#    Ftest[1, ] <- c(f, f.q, f.p)

#    colnames(Ftest) <- c("stat", "quantile_lower", "quantile_upper", "p_value")
    
#    return(list(stat = Ftest, boot = f.boot, graph = p))
    
#} ## end of test