###################################################################
## IFE Model Function
###################################################################
fect.polynomial <- function(Y, # Outcome variable, (T*N) matrix
                            X, # Explanatory variables:  (T*N*p) array
                            D, #  Indicator for treated unit (tr==1) 
                            I,
                            II, 
                            T.on, 
                            T.off = NULL, 
                            T.on.carry = NULL,           
                            T.on.balance = NULL,
                            balance.period = NULL,
                            method,
                            degree = 1,
                            sfe = NULL,
                            cfe = NULL,
                            ind.matrix = NULL,
                            knots = NULL,
                            force, 
                            hasRevs = 1,
                            tol, # tolerance level
                            boot = FALSE, # bootstrapped sample
                            placeboTest = 0,
                            placebo.period = NULL,
                            carryoverTest = 0,
                            carryover.period = NULL,
                            norm.para = NULL,
                            calendar.enp.seq = NULL,
                            time.on.seq = NULL,
                            time.off.seq = NULL,
                            time.on.carry.seq = NULL,
                            time.on.balance.seq = NULL,
                            group.level = NULL,
                            group = NULL,
                            time.on.seq.group = NULL,
                            time.off.seq.group = NULL) {  
    
    ##-------------------------------##
    ## Parsing data
    ##-------------------------------##  
    carryover.pos <- placebo.pos <- na.pos <- NULL
    res.sd1 <- res.sd2 <- NULL

    ## unit id and time
    TT <- dim(Y)[1]
    N <- dim(Y)[2]
    if (is.null(X) == FALSE) {
        p <- dim(X)[3]
    } 
    else {
        p <- 0
        X <- array(0, dim = c(1, 1, 0))
    }

    ## replicate data
    YY <- Y
    YY[which(II == 0)] <- 0 ## reset to 0 

    ## initial fit using fastplm
    data.ini <- matrix(NA, (TT*N), (2 + 1 + p))
    data.ini[, 2] <- rep(1:N, each = TT)         ## unit fe
    data.ini[, 3] <- rep(1:TT, N)                ## time fe
    data.ini[, 1] <- c(Y)                       ## outcome
    if (p > 0) {                                ## covar
        for (i in 1:p) {
            data.ini[, (3 + i)] <- c(X[, , i])
        }
    }
    ## observed Y0 indicator:
    oci <- which(c(II) == 1)
    initialOut <- initialFit(data = data.ini, force = force, oci = oci)
    Y0 <- initialOut$Y0
    beta0 <- initialOut$beta0
    if (p > 0 && sum(is.na(beta0)) > 0) {
        beta0[which(is.na(beta0))] <- 0
    }
    est.fect <- NULL
    if (boot == FALSE) {
        est.fect <- inter_fe_ub(YY, Y0, X, II, beta0, 0, force = force, tol)
    }

    ## reshape 
    vy <- as.matrix(c(YY))
    vx.fit <- vx <- NULL
    if (p > 0) {
        vx <- matrix(NA, N*TT, p)
        for (i in 1:p) {
            vx[, i] <- c(X[,, i])
        }
        colnames(vx) <- paste0("x.",c(1:p))
        vx.fit <- as.matrix(vx[oci,])
    }
    
    vindex <- NULL
    sp <- NULL
    sf <- NULL
    cf <- NULL

    if(method == "cfe"){
        vindex <- cbind(rep(1:N, each = TT), rep(1:TT, N))  ## id time
        if (force == 1) {
            sf <- 1
        } 
        else if (force == 2) {
            sf <- 2
        } 
        else {
            sf <- c(1,2)
        }

        ## simple fixed effects
        for(ind.name in names(ind.matrix)){
            vindex <- cbind(vindex,matrix(ind.matrix[[ind.name]],ncol=1))
        }

        ind.name <- c("forceid","forcetime",names(ind.matrix))
        ind.index <- c(1:(2+length(names(ind.matrix))))
        colnames(vindex) <- names(ind.index) <- ind.name
        
        if(p>0){
            data.reg <- cbind.data.frame(vy,vx,vindex)
            formula.reg <- paste0("vy~",paste(paste0("x.",c(1:p)),collapse="+"),"|")    
        }
        else{
            data.reg <- cbind.data.frame(vy,vindex)
            formula.reg <- paste0("vy~1|")   
        }
        
        if(force==1){
            formula.reg <- paste0(formula.reg,"forceid")
        }
        else if(force==2){
            formula.reg <- paste0(formula.reg,"forcetime")
        }
        else if(force==3){
            formula.reg <- paste0(formula.reg,"forceid+forcetime")
        }


        if(!is.null(sfe)){
            formula.reg <- paste0(formula.reg,"+", paste(sfe,collapse="+"))
        }

        if(!is.null(cfe)){
            for(sub.cfe in cfe){
                sub.cf <- paste0(sub.cfe[1],"[",sub.cfe[2],"]")
                formula.reg <- paste0(formula.reg,"+",sub.cf)
            }
        }
        formula.reg <- as.formula(formula.reg)
        est.best <- suppressWarnings(invisible(feols(fml = formula.reg,
                                                   data = data.reg[oci,],
                                                   fixef.rm = "none")))
        
        yfit <- suppressWarnings(predict(est.best, newdata = data.reg))

    }
    else if (method == "polynomial") {
        vindex <- cbind(rep(1:N, each = TT), rep(1:TT, N))  ## id time
        for (i in 1:degree) {
            vindex <- cbind(vindex, rep((1:TT)^i, N))
        }
        
        colnames(vindex) <- c("forceid","forcetime",paste0("forcetime.",c(1:degree)))
        
        if(p>0){
            data.reg <- cbind.data.frame(vy,vx,vindex)
            formula.reg <- paste0("vy~",paste(paste0("x.",c(1:p)),collapse="+"),"|")    
        }
        else{
            data.reg <- cbind.data.frame(vy,vindex)
            formula.reg <- paste0("vy~1|")   
        }

        if(force==1){
            formula.reg <- paste0(formula.reg,"forceid")
        }
        else if(force==2){
            formula.reg <- paste0(formula.reg,"forcetime")
        }
        else if(force==3){
            formula.reg <- paste0(formula.reg,"forceid+forcetime")
        }

        for (i in 1:degree) {
            formula.reg <- paste0(formula.reg,paste0("+forceid","[",paste0("forcetime.",i),"]"))
        }
        formula.reg <- as.formula(formula.reg)
        est.best <- suppressWarnings(invisible(feols(fml = formula.reg,
                                                   data = data.reg[oci,],
                                                   fixef.rm = "none")))
        
        yfit <- suppressWarnings(predict(est.best, newdata = data.reg))
    } 


    Y.ct <- matrix(yfit, TT, N)
    if (p > 0) {
        beta <- as.matrix(c(est.best$coefficients)[1:p])
    } 
    else {
        beta <- matrix(0, 1, 0)
    }
    est.best$beta <- beta
    validX <- ifelse(p > 0, 1, 0) 

    ##-------------------------------##
    ##   ATT and Counterfactuals     ##
    ##-------------------------------##

    ## we first adjustment for normalization 
    if (!is.null(norm.para)) {
        Y <- Y * norm.para[1]
        ## variance of the error term 
        sigma2 <- est.best$sigma2 * (norm.para[1]^2)    
        est.best$sigma2 <- sigma2
        ## output of estimates
        est.best$mu <- est.best$mu * norm.para[1] 
        est.best$residuals <- est.best$residuals * norm.para[1]
        Y.ct <- Y.ct * norm.para[1]
        if (boot == FALSE) {
            est.fect$fit <- est.fect$fit * norm.para[1]
        }
        est.fect$sigma2 <- est.fect$sigma2 * norm.para[1]
    }

    ## 0. relevant parameters
    if (p>0) {
        na.pos <- is.nan(est.best$beta)
        beta <- est.best$beta
        if( sum(na.pos) > 0 ) {
            beta[na.pos] <- NA
        }
    } else {
        beta <- NA
    }
   
    ## 1. estimated att and counterfactuals
    eff <- Y - Y.ct   
    complete.index <- which(!is.na(eff))
    missing.index <- which(is.na(eff))
    if(length(missing.index)>0){
        I[missing.index] <- 0
        II[missing.index] <- 0
    }
    att.avg <- sum(eff[complete.index] * D[complete.index])/(sum(D[complete.index]))

    att.avg.balance <- NA
    if(!is.null(balance.period)){
        complete.index2 <- which(!is.na(T.on.balance))
        att.avg.balance <- sum(eff[complete.index2] * D[complete.index2])/(sum(D[complete.index2]))
    }


    ## att.avg.unit
    tr.pos <- which(apply(D, 2, sum) > 0)
    att.unit <- sapply(1:length(tr.pos), function(vec){return((eff[, tr.pos[vec]] * D[, tr.pos[vec]]) / sum(D[, tr.pos[vec]]))})
    att.avg.unit <- mean(att.unit,na.rm=TRUE)


    ## 2. rmse for treated units' observations under control
    tr <- which(apply(D, 2, sum) > 0)
    tr.co <- which((as.matrix(1 - D[,tr]) * as.matrix(II[,tr])) == 1)
    eff.tr <- as.matrix(eff[,tr])
    v.eff.tr <- eff.tr[tr.co]
    rmse <- sqrt(mean(v.eff.tr^2,na.rm=TRUE))
    

    ## 3. unbalanced output
    if (0 %in% I) {
        eff[which(I == 0)] <- NA
        Y.ct[which(I == 0)] <- NA
    }
      
    ## 4. dynamic effects
    t.on <- c(T.on)
    eff.v <- c(eff) ## a vector

    rm.pos1 <- which(is.na(eff.v))
    rm.pos2 <- which(is.na(t.on)) 

    eff.v.use1 <- eff.v
    t.on.use <- t.on
    n.on.use <- rep(1:N, each = TT)

    if (NA %in% eff.v | NA %in% t.on) {
        eff.v.use1 <- eff.v[-c(rm.pos1, rm.pos2)]
        t.on.use <- t.on[-c(rm.pos1, rm.pos2)]
        n.on.use <- n.on.use[-c(rm.pos1, rm.pos2)]
    }

    pre.pos <- which(t.on.use <= 0)
    eff.pre <- cbind(eff.v.use1[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
    colnames(eff.pre) <- c("eff", "period", "unit")

    sigma2.pre <- eff.pre.equiv <- NULL
    if (boot == FALSE) {
        eff.pre.equiv <- eff.pre

        sigma2.pre <- tapply(eff.pre.equiv[,1], eff.pre.equiv[,2], var)
        sigma2.pre <- cbind(sigma2.pre, sort(unique(eff.pre.equiv[, 2])), table(eff.pre.equiv[, 2]))
        colnames(sigma2.pre) <- c("sigma2", "period", "count")
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

    ## 4.1 carryover effect 
    carry.att <- NULL 
    if (!is.null(T.on.carry)) {
        t.on.carry <- c(T.on.carry)
        rm.pos4 <- which(is.na(t.on.carry)) 
        t.on.carry.use <- t.on.carry

        if (NA %in% eff.v | NA %in% t.on.carry) {
            eff.v.use3  <- eff.v[-c(rm.pos1, rm.pos4)]
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
            eff.v.use3  <- eff.v[-c(rm.pos1, rm.pos4)]
            t.on.balance.use <- t.on.balance[-c(rm.pos1, rm.pos4)]        
        }

        balance.time <- sort(unique(t.on.balance.use))
        balance.att <- as.numeric(tapply(eff.v.use3, t.on.balance.use, mean)) ## NA already removed
        balance.count <- as.numeric(table(t.on.balance.use))

        if (!is.null(time.on.balance.seq)) {
            balance.att.med <- rep(NA, length(time.on.balance.seq))
            balance.count.med <- rep(0, length(time.on.balance.seq))
            balance.att.med[which(time.on.balance.seq %in% balance.time)] <- balance.att
            if(length(balance.count)>0){
                balance.count.med[which(time.on.balance.seq %in% balance.time)] <- balance.count                
            }
            balance.count <- balance.count.med
            balance.att <- balance.att.med
            balance.time <- time.on.balance.seq
        }

        #placebo for balanced samples
        if(!is.null(placebo.period) && placeboTest == 1){
            if (length(placebo.period) == 1) {
                balance.placebo.pos <- which(balance.time == placebo.period)
                balance.att.placebo <- balance.att[balance.placebo.pos]
            } 
            else {
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
            eff.off.equiv <- eff.off

            off.sd <- tapply(eff.off.equiv[,1], eff.off.equiv[,2], sd)
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
        for(i in c(1:length(group.level))){
            sub.group <- group.level[i]
            sub.group.name <- group.level.name[i]

            ## by-group dynamic effects
            t.on.sub <- c(T.on[which(group==sub.group)])
            eff.v.sub <- c(eff[which(group==sub.group)]) ## a vector
            rm.pos1.sub <- which(is.na(eff.v.sub))
            rm.pos2.sub <- which(is.na(t.on.sub)) 
            eff.v.use1.sub <- eff.v.sub
            t.on.use.sub <- t.on.sub
            if (NA %in% eff.v.sub | NA %in% t.on.sub) {
                eff.v.use1.sub <- eff.v.sub[-c(rm.pos1.sub, rm.pos2.sub)]
                t.on.use.sub <- t.on.sub[-c(rm.pos1.sub, rm.pos2.sub)]
            }
            if(length(t.on.use.sub)>0){
                time.on.sub <- sort(unique(t.on.use.sub))
                att.on.sub <- as.numeric(tapply(eff.v.use1.sub, 
                                            t.on.use.sub, 
                                            mean)) ## NA already removed
                count.on.sub <- as.numeric(table(t.on.use.sub))
            }else{
                time.on.sub <- att.on.sub <- count.on.sub <- NULL
            }
            
            if (!is.null(time.on.seq.group)) {
                count.on.med.sub <- att.on.med.sub <- rep(NA, length(time.on.seq.group[[sub.group.name]]))
                time.on.seq.sub <- time.on.seq.group[[sub.group.name]]
                att.on.med.sub[which(time.on.seq.sub %in% time.on.sub)] <- att.on.sub
                count.on.med.sub[which(time.on.seq.sub %in% time.on.sub)] <- count.on.sub
                att.on.sub <- att.on.med.sub
                count.on.sub <- count.on.med.sub
                time.on.sub<- time.on.seq.sub
            }
            suboutput <- list(att.on=att.on.sub,
                              time.on=time.on.sub,
                              count.on=count.on.sub)

            ## placebo effect, if placeboTest == 1 
            if (!is.null(placebo.period) && placeboTest == 1) {              
                if (length(placebo.period) == 1) {
                    placebo.pos.sub <- which(time.on.sub == placebo.period)
                    if(length(placebo.pos.sub)>0){
                        att.placebo.sub <- att.on.sub[placebo.pos.sub]
                    }
                    else{att.placebo.sub <- NULL} 
                } 
                else {
                    placebo.pos.sub <- which(time.on.sub >= placebo.period[1] & time.on.sub <= placebo.period[2])
                    if(length(placebo.pos.sub)>0){
                        att.placebo.sub <- sum(att.on.sub[placebo.pos.sub] * count.on.sub[placebo.pos.sub]) / sum(count.on.sub[placebo.pos.sub])
                    }
                    else{att.placebo.sub <- NULL} 
                }
                suboutput <- c(suboutput, list(att.placebo = att.placebo.sub))
            }

            ## T.off
            if (hasRevs == 1) {    
                t.off.sub <- c(T.off[which(group==sub.group)])
                rm.pos3.sub <- which(is.na(t.off.sub))
                eff.v.use2.sub <- eff.v.sub
                t.off.use.sub <- t.off.sub
                if (NA %in% eff.v.sub | NA %in% t.off.sub) {
                    eff.v.use2.sub <- eff.v.sub[-c(rm.pos1.sub, rm.pos3.sub)]
                    t.off.use.sub <- t.off.sub[-c(rm.pos1.sub, rm.pos3.sub)]
                }
                if(length(t.off.use.sub)>0){
                    time.off.sub <- sort(unique(t.off.use.sub))
                    att.off.sub <- as.numeric(tapply(eff.v.use2.sub, t.off.use.sub, mean)) ## NA already removed
                    count.off.sub <- as.numeric(table(t.off.use.sub))
                }else{
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
                suboutput <- c(suboutput, list(att.off = att.off.sub,
                                               count.off = count.off.sub,
                                               time.off = time.off.sub))

                if (!is.null(carryover.period) && carryoverTest == 1) {
                    if (length(carryover.period) == 1) {
                        carryover.pos.sub <- which(time.off.sub == carryover.period)
                        if(length(carryover.pos.sub)>0){
                            att.carryover.sub <- att.off.sub[carryover.pos.sub]
                        } else{att.carryover.sub <- NULL}      
                    } else {
                        carryover.pos.sub <- which(time.off.sub >= carryover.period[1] & time.off.sub <= carryover.period[2])
                        if(length(carryover.pos.sub)>0){
                            att.carryover.sub <- sum(att.off.sub[carryover.pos.sub] * count.off.sub[carryover.pos.sub]) / sum(count.off.sub[carryover.pos.sub])
                        } else{att.carryover.sub <- NULL}  
                    }
                    suboutput <- c(suboutput,list(att.carryover = att.carryover.sub))
                }
            }
            group.output[[sub.group.name]] <- suboutput
        }
    }

    ## 9. loess HTE by time
    D.missing <- D
    D.missing[which(D==0)] <- NA
    eff.calendar <- apply(eff*D.missing,1,mean,na.rm=TRUE)
    N.calendar <- apply(!is.na(eff*D.missing),1,sum)
    T.calendar <- c(1:TT)
    if(sum(!is.na(eff.calendar))>1){
        #loess fit
        if(!is.null(calendar.enp.seq)){
            if(length(calendar.enp.seq)==1 & is.na(calendar.enp.seq)){
                calendar.enp.seq <- NULL
            }
        }
        if(is.null(calendar.enp.seq)){
            loess.fit <- suppressWarnings(try(loess(eff.calendar~T.calendar,weights = N.calendar),silent=TRUE))      
        }
        else{
            loess.fit <- suppressWarnings(try(loess(eff.calendar~T.calendar,weights = N.calendar,enp.target=calendar.enp.seq),silent=TRUE))
        }
        if('try-error' %in% class(loess.fit)){
            eff.calendar.fit <- eff.calendar
            calendar.enp <- NULL
        }
        else{
            eff.calendar.fit <- eff.calendar
            eff.calendar.fit[which(!is.na(eff.calendar))] <- loess.fit$fit
            calendar.enp <- loess.fit$enp              
        }
    }
    else{
        eff.calendar.fit <- eff.calendar
        calendar.enp <- NULL
    }

  
    ##-------------------------------##
    ##            Storage            ##
    ##-------------------------------##  
    out<-list(
        ## main results 
        method = method,
        Y.ct = Y.ct,
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
        p = p,
        beta = beta,
        est = est.best,
        sigma2 = est.best$sigma2,
        sigma2.fect = est.fect$sigma2,
        validX = validX,
        time = time.on,
        att = att.on,
        count = count.on,
        eff.calendar = eff.calendar,
        N.calendar = N.calendar,
        eff.calendar.fit = eff.calendar.fit,
        calendar.enp = calendar.enp,
        eff.pre = eff.pre,
        eff.pre.equiv = eff.pre.equiv,
        sigma2.pre = sigma2.pre)
    
    #print(att.on)

    if (hasRevs == 1) {
        out <- c(out, list(time.off = time.off, 
                           att.off = att.off,
                           count.off = count.off,
                           eff.off = eff.off,
                           eff.off.equiv = eff.off.equiv,
                           off.sd = off.sd))
    }

    if (!is.null(T.on.carry)) {
        out <- c(out, list(carry.att = carry.att, carry.time = carry.time))
    }

    if(!is.null(balance.period)){
        out <- c(out, list(balance.att = balance.att, balance.time = balance.time,balance.count = balance.count,balance.avg.att = att.avg.balance))        
        if (!is.null(placebo.period) && placeboTest == 1) {
            out <- c(out, list(balance.att.placebo = balance.att.placebo))
        }    
    }

    if (!is.null(placebo.period) && placeboTest == 1) {
        out <- c(out, list(att.placebo = att.placebo))
    }

    if (!is.null(carryover.period) && carryoverTest == 1) {
        out <- c(out, list(att.carryover = att.carryover))
    }

    if (!is.null(group)) {
        out <- c(out, list(group.att = group.att,
                           group.output = group.output))
    }
    return(out)
} 




