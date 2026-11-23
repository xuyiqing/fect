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
            
    
    ## computing
    if (!identical(parallel, FALSE)) {
        boot.out <- foreach(j=1:nboots, 
                            .inorder = FALSE,
                            .export = c("fect_fe", "fect_mc", "get_term"),
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
