## Causal inference using counterfactual estimators 
## (fect: fixed effects counterfactuals)
## Version 1.4.4
## Author: Licheng Liu (Tsinghua), Ye Wang(NYU), Yiqing Xu(Stanford)
## Date: 2019.01.02

## MAIN FUNCTION
## fect.formula()
## fect.default()

## DEPENDENT FUNCTIONS
## fect.fe() ## interactive fixed effects model
## fect.mc() ## matrix completion
## fect.boot() ## bootstrap 

## fitness test
## fect.test ## wild bootstrap

## METHODS
## print.fect()
## plot.fect()

## several supporting functions


#####################################################################
## A Shell Function
#####################################################################

## generic function
fect <- function(formula = NULL, data, # a data frame (long-form)
                 Y, # outcome
                 D, # treatment 
                 X = NULL, # time-varying covariates
                 na.rm = FALSE, # remove missing values
                 index, # c(unit, time) indicators
                 force = "unit", # fixed effects demeaning
                 cl = "unit", 
                 r = 0, # nubmer of factors
                 lambda = NULL, # mc method: regularization parameter
                 nlambda = 10, ## mc method: regularization parameter
                 CV = TRUE, # cross-validation
                 k = 10, # times of CV
                 cv.prop = 0.1, ## proportion of CV counts
                 cv.treat = TRUE, 
                 cv.continue = 3,
                 binary = FALSE, # probit model
                 QR = FALSE, # QR or SVD for binary probit 
                 method = "fe", # method: e for fixed effects; ife for interactive fe; mc for matrix completion
                 criterion = "mspe", # for ife model: mspe, pc or both
                 se = FALSE, # report uncertainties
                 nboots = 200, # number of bootstraps
                 parallel = FALSE, # parallel computing
                 cores = NULL, # number of cores
                 tol = 0.001, # tolerance level
                 seed = NULL, # set seed
                 min.T0 = 5, # minimum T0
                 max.missing = NULL, # maximum missing
                 pre.period = NULL, # fit test period 
                 wald = FALSE, # fit test
                 placebo.period = NULL, # placebo test period
                 placeboTest = FALSE, # placebo test
                 normalize = FALSE # accelerate option
                ) {
    UseMethod("fect")
}

## formula method

fect.formula <- function(formula = NULL,data, # a data frame (long-form)
                         Y, # outcome
                         D, # treatment 
                         X = NULL, # time-varying covariates
                         na.rm = FALSE, # remove missing values
                         index, # c(unit, time) indicators
                         force = "unit", # fixed effects demeaning
                         cl = "unit", 
                         r = 0, # nubmer of factors
                         lambda = NULL, # mc method: regularization parameter
                         nlambda = 10, ## mc method: regularization parameter
                         CV = TRUE, # cross-validation
                         k = 10, # times of CV
                         cv.prop = 0.1, ## proportion of CV counts
                         cv.treat = TRUE, 
                         cv.continue = 3,
                         binary = FALSE, # probit model
                         QR = FALSE, # QR or SVD for binary probit 
                         method = "fe", # method: fe for fixed effects; ife for interactive fe; mc for matrix completion
                         criterion = "mspe", # for ife model: mspe, pc or both
                         se = FALSE, # report uncertainties
                         nboots = 200, # number of bootstraps
                         parallel = FALSE, # parallel computing
                         cores = NULL, # number of cores
                         tol = 0.001, # tolerance level
                         seed = NULL, # set seed
                         min.T0 = 5,
                         max.missing = NULL,
                         pre.period = NULL,
                         wald = FALSE,
                         placebo.period = NULL,
                         placeboTest = FALSE,
                         normalize = FALSE
                        ) {
    ## parsing
    varnames <- all.vars(formula)
    Yname <- varnames[1]
    Dname <- varnames[2]
    if (length(varnames) > 2) {
        Xname <- varnames[3:length(varnames)]
    } else {
        Xname <- NULL
    }

    namesData <- colnames(data)
    for (i in 1:length(varnames)) {
        if(!varnames[i] %in% namesData) {
            stop(paste("variable \"", varnames[i],"\" is not in the data set.", sep = ""))
        }
    }

    ## check binary outcome
    if (binary == TRUE) {
        unique_y <- sort(unique(data[,Yname]))
        if (length(unique_y) != 2) {
            stop("Outcome should only contain 0 and 1.")
        } else {
            if (sum(unique_y == c(0,1)) != 2) {
                stop("Outcome should only contain 0 and 1.")
            }
        } 
    }

    ## run the model
    out <- fect.default(formula = NULL, data = data, Y = Yname,
                        D = Dname, X = Xname,
                        na.rm, index, force, cl, r, lambda, nlambda, 
                        CV, k, cv.prop, cv.treat, cv.continue,
                        binary, QR, method, criterion, se, 
                        nboots, parallel, cores, tol, seed, min.T0,
                        max.missing, pre.period, wald,
                        placebo.period, placeboTest, normalize)
    
    out$call <- match.call()
    out$formula <- formula
    ## print(out)
    return(out)

}


## default function

fect.default <- function(formula = NULL, data, # a data frame (long-form)
                         Y, # outcome
                         D, # treatment 
                         X = NULL, # time-varying covariates
                         na.rm = FALSE, # remove missing values
                         index, # c(unit, time) indicators
                         force = "unit", # fixed effects demeaning
                         cl = "unit", 
                         r = 0, # nubmer of factors
                         lambda = NULL, ## mc method: regularization parameter
                         nlambda = 0, ## mc method: regularization parameter
                         CV = TRUE, # cross-validation
                         k = 10, # times of CV
                         cv.prop = 0.1,
                         cv.treat = TRUE, 
                         cv.continue = 3,
                         binary = FALSE, # probit model
                         QR = FALSE, # QR or SVD for binary probit 
                         method = "fe", # method: ife for interactive fe; mc for matrix completion
                         criterion = "mspe", 
                         se = FALSE, # report uncertainties
                         nboots = 200, # number of bootstraps
                         parallel = FALSE, # parallel computing
                         cores = NULL, # number of cores
                         tol = 0.001, # tolerance level
                         seed = NULL, # set seed
                         min.T0 = 5,
                         max.missing = NULL,
                         pre.period = NULL,
                         wald = FALSE,
                         placebo.period = NULL,
                         placeboTest = FALSE,
                         normalize = FALSE
                        ) {  
    
    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------## 
    ## read data 
    if (is.data.frame(data) == FALSE || length(class(data)) > 1) {
        data <- as.data.frame(data)
        ## warning("Not a data frame.")
    }
    ## index
    if (length(index) != 2 | sum(index %in% colnames(data)) != 2) {
        stop("\"index\" option misspecified. Try, for example, index = c(\"unit.id\", \"time\").")
    }
    
    ## check duplicated observations
    unique_label <- unique(paste(data[,index[1]],"_",data[,index[2]],sep=""))
    if (length(unique_label)!= dim(data)[1]) {
        stop("Some records may be duplicated or wrongly marked in the data set. Check the index.")
    }

    ## force
    if (force == "none") { # force = 0 "none": no additive fixed effects imposed
        force <- 0
    } else if (force == "unit") { # force = 1 "unit": unit fixed-effect (default)
        force <- 1
    } else if (force == "time") { # force = 2 "time": time fixed-effect
        force <- 2
    } else if (force == "two-way") { # force = 3 "two-way": two-way fixed-effect 
        force <- 3
    }
    if (!force %in% c(0, 1, 2, 3)) {
        stop("\"force\" option misspecified; choose from c(\"none\", \"unit\", \"time\", \"two-way\").")
    } 

    ## binary
    if (binary == 1) {
        method <- "ife"
        normalize <- FALSE
    }

    ## method
    if (!method %in% c("fe", "ife", "mc", "both")) {
        stop("\"method\" option misspecified; choose from c(\"fe\", \"ife\", \"mc\", \"both\").")
    }
    if (method == "fe") {
        r <- 0
        CV <- FALSE
        method <- "ife"
    }

    if (!criterion %in% c("mspe", "pc", "both")) {
        stop("\"method\" option misspecified; choose from c(\"mspe\", \"pc\", \"both\").")
    }

    if (method == "both" && criterion == "pc") {
        stop("\"ife\" and \"mc\" method cannot be compared; try using \"mspe\" as criteria.")
    }

    ## r
    if ( method %in% c("ife", "both") & r[1] < 0) {
        stop("\"r\" option misspecified. The number of factors must be non-negative.")
    }

    ## lambda
    if (method %in% c("mc", "both")) {
        if (!is.null(lambda)) {
            if (sum(lambda < 0) > 0) {
                stop("\"lambda\" option misspecified. It must be non-negative.")    
            }
        }
        if (CV == FALSE & is.null(lambda)) {
            cat("No lambda is supplied. FEct is applied.")
            method <- "ife"
            r <- 0
        }
    } 

    ## CV
    if (method == "both") {
        CV <- TRUE
    }
    if (CV == TRUE) {
        
        if (placeboTest == TRUE) {
            stop("Placebo test cannot be performed while doing cross-validation.")
        }

        if (method %in% c("ife", "both")) {
            if (length(r) == 2 & r[1] > r[2]) {
                stop("\"r\" option misspecified.")
            }
        } 
        if (method %in% c("mc", "both")) {
            if (nlambda <= 0) {
                stop("\"nlambda\" option misspecified.")
            }
        }
    } else {
        if (! method %in% c("ife", "mc")) {
            stop("\"method\" option misspecified; please choose from c(\"ife\", \"mc\").")
        }
    }

    if (placeboTest == TRUE) {
        wald <- FALSE
    }

    if (length(r) == 1) {
        if (r>=5) {
            r.end <- r
        } else {
            r.end <- 5
        }
    } else {
        r.end <- r[2]; r <- r[1]
    }

    ## uncertainty estimates
    if (is.logical(se) == FALSE & !se%in%c(0, 1)) {
        stop("\"se\" is not a logical flag.")
    } 

    ## normalize
    if (is.logical(normalize) == FALSE & !normalize%in%c(0, 1)) {
        stop("\"normalize\" is not a logical flag.")
    } 

    ## nboots
    if (se == TRUE & nboots <= 0) {
        stop("\"nboots\" option misspecified. Try, for example, nboots = 200.")
    }

    ## parallel & cores
    if (parallel == TRUE) {
        if (is.null(cores) == FALSE) {
            if (cores <= 0) {
                stop("\"cores\" option misspecified. Try, for example, cores = 2.")
            }
        }
    } 

    ## tol
    if (tol <= 0) {
        stop("\"tol\" option misspecified. Try using the default option.")
    }

    ## seed
    if (is.null(seed) == FALSE) {
        if (is.numeric(seed) == FALSE) {
            stop("\"seed\" should be a number.")
        }
    }

    ## remove missing values
    if (is.logical(na.rm) == FALSE & !na.rm%in%c(0, 1)) {
        stop("\"na.rm\" is not a logical flag.")
    } 

    ## select variable that are to be used 
    #if (!is.null(cl)) {
    #    if (cl %in% index) {
    #        data <- data[,c(index, Y, D, X)]
    #    } else {
    #        data <- data[,c(index, Y, D, X, cl)]
    #    }
    #} else {
        data <- data[,c(index, Y, D, X)] ## some variables may not be used
    #}
    
    if (na.rm == TRUE) {
        data <- na.omit(data)
    }

    ## placebo period
    if (placeboTest == TRUE && !is.null(placebo.period)) {
        if (sum(placebo.period > 0) > 0) {
            stop("\"placebo.period\" should not be greater than 0.")
        } else {
            if (length(placebo.period) > 2) {
                stop("\"placebo.period\" option misspecified. ")
            }
        }
    }

    ##-------------------------------##
    ## Parsing raw data
    ##-------------------------------##  

    ## store data and variable names
    data.old <- data
    Yname <- Y
    Dname <- D
    Xname <- X
    clname <- cl

    #if (!is.null(clname)) {
    #    if (!clname %in% index) {
    #        data[, clname] <- as.numeric(as.factor(data[, clname]))
    #    }
    #}

    ## normalize
    norm.para <- NULL
    if (normalize == TRUE) {
        sd.Y <- sd(as.matrix(data[,Yname]))
        data[,c(Yname, Xname)] <- data[,c(Yname, Xname)]/sd.Y
        norm.para <- sd.Y ## normalized parameter
    }

    ## check index and treatment indicator
    if (! class(data[, Dname]) %in% c("numeric", "integer")) {
        ## data[, Dname] <- as.numeric(as.character(data[, Dname]))
        stop("Treatment indicator should be a numeric value.")
    } 

    if (class(data[, index[1]]) == "factor") {
        data[, index[1]] <- as.character(data[, index[1]])
    } 

    if (class(data[, index[2]]) == "factor") {
        data[, index[2]] <- as.character(data[, index[2]])
    } 
    
    id <- index[1]
    time <- index[2]
    TT <- length(unique(data[,time]))
    N <- length(unique(data[,id]))
    p <- length(Xname)
    id.series <- unique(sort(data[,id])) ## unit id
    time.uni <- unique(sort(data[,time])) ## period

    ## sort data
    data <- data[order(data[,id], data[,time]), ]

    ## max.missing
    if (is.null(max.missing)) {
        max.missing <- TT
    }
    
    ## check missingness
    if (sum(is.na(data[, Yname])) > 0) {
        stop(paste("Missing values in variable \"", Yname,"\".", sep = ""))
    }
    if (sum(is.na(data[, Dname])) > 0) {
        stop(paste("Missing values in variable \"", Dname,"\".", sep = ""))
    }
    if (!(1%in%data[, Dname] & 0%in%data[,Dname] & length(unique(data[,Dname])) == 2)) {
        stop(paste("Error values in variable \"", Dname,"\".", sep = ""))
    }


    ## check variation in x
    if (p > 0) {
        for (i in 1:p) {
            if (sum(is.na(data[, Xname[i]])) > 0) {
                stop(paste("Missing values in variable \"", Xname[i],"\".", sep = ""))
            }

            if (sum(tapply(data[, Xname[i]], data[, id], var), na.rm = TRUE) == 0) {
                stop(paste("Variable \"",Xname[i], "\" is unit-invariant. Try to remove it.", sep = ""))
            }
            if (sum(tapply(data[, Xname[i]], data[, time], var), na.rm = TRUE) == 0) {
                stop(paste("Variable \"",Xname[i], "\" is time-invariant. Try to remove it.", sep = ""))
            }
        }
    }

    ## check index 
    if (sum(is.na(data[, id])) > 0) {
        stop(paste("Missing values in variable \"", id,"\".", sep = ""))
    }
    if (sum(is.na(data[, time])) > 0) {
        stop(paste("Missing values in variable \"", time,"\".", sep = ""))
    } 

    ## check balanced panel and fill unbalanced panel
    if (dim(data)[1] < TT*N) {
        
        data[,time] <- as.numeric(as.factor(data[,time]))
        ## ob <- "time_ob_ls"
        
        ## while (ob %in% colnames(data)) {
        ##     ob <- paste(ob, ob, sep = "_")
        ## }

        ## data[, ob] <- data[, time]
        ## for (i in 1:N) {
        ##     data[data[,id] == id.series[i], ob] <- data[data[,id] == id.series[i],time] + (i - 1) * TT  
        ## }

        ob.indicator <- data[,time]
        id.indicator <- table(data[, id])
        sub.start <- 1
        for (i in 1:(N - 1)) { 
            sub.start <- sub.start + id.indicator[i] 
            sub.end <- sub.start + id.indicator[i+1] - 1 
            ob.indicator[sub.start:sub.end] <- ob.indicator[sub.start:sub.end] + i * TT
        }

        #if (!is.null(clname)) {
        #    if (!clname %in% index) {
        #        variable <- c(Yname, Dname, Xname, clname)
        #    } else {
        #        variable <- c(Yname, Dname, Xname)
        #    }
        #} else {
            variable <- c(Yname, Dname, Xname)
        #}

        data_I <- matrix(0, N * TT, 1)
        data_I[ob.indicator, 1] <- 1
        data_ub <- as.matrix(data[, variable])
        data <- data_ub_adj(data_I, data_ub)
        colnames(data) <- variable
    }

    ## indicator matrix: index matrix that indicates if data is observed 
    I <- matrix(1, TT, N)
    Y.ind <- matrix(data[, Yname], TT, N)
    I[is.nan(Y.ind)] <- 0
    if (0%in%I) {
        data[is.nan(data)] <- 0
    }    
    ##treatment indicator: incorporates reversal treatments
    D <- matrix(data[, Dname], TT, N)
    ##outcome variable
    Y <- matrix(data[, Yname], TT, N)
    ## time-varying covariates
    X <- array(0, dim = c(TT, N, p))
    ## xp <- rep(0, p) ## label invariant x
    ## x.pos <- 0
    if (p > 0) {
        ## x.pos <- 1:p
        for (i in 1:p) {
            X[,,i] <- matrix(data[, Xname[i]], TT, N)
            ## if (force %in% c(1,3)) {
            ##     if (!0%in%I) {
            ##         tot.var.unit <- sum(apply(X[, , i], 2, var))
            ##     } else {
            ##         Xi <- X[,,i]
            ##         Xi[which(I == 0)] <- NA
            ##         tot.var.unit <- sum(apply(Xi, 2, var, na.rm = TRUE))
            ##     }
            ##     if(!is.na(tot.var.unit)) {
            ##         if (tot.var.unit == 0) {
                        ## time invariant covar can be removed
            ##             xp[i] <- 1
            ##             cat(paste("Variable \"", Xname[i],"\" is time-invariant.\n", sep = ""))   
            ##         }
            ##     }
            ## }
            ## if (force %in% c(2, 3)) {
            ##     if (!0%in%I) {
            ##         tot.var.time <- sum(apply(X[, , i], 1, var))
            ##     } else {
            ##         Xi <- X[,,i]
            ##         Xi[which(I == 0)] <- NA
            ##         tot.var.time <- sum(apply(Xi, 1, var, na.rm = TRUE))
            ##     } 
            ##     if (!is.na(tot.var.time)) {
            ##         if (tot.var.time == 0) {
                        ## can be removed in inter_fe
            ##             xp[i] <- 1
            ##             cat(paste("Variable \"", Xname[i],"\" has no cross-sectional variation.\n", sep = ""))
            ##         }
            ##     }
            ## } 
        } 
    }
    #if (!is.null(clname)) {
    #    if (clname %in% index) {
    #        cl <- 1:N
    #    } else {
    #        cl <- matrix(data[, clname], TT, N)
    #        v.cl <- c()
    #        for (i in 1:N) {
    #            if (sum(is.na(cl[,i])) > 0) {
    #                v.cl <- c(v.cl, na.omit(cl[,i])[1])
    #            } else {
    #                v.cl <- c(v.cl, cl[1, i])
    #            }
    #        }
    #        cl <- v.cl
    #    }
    #} else {
    #    cl <- NULL
    #}

    ## ----------------------------------------------------------- ##
    II <- I

    ## if (placeboTest == 1) {
    ##     T.on <- matrix(NA, TT, (N - length(rm.id)))
    ##     for (i in 1:(N - length(rm.id))) {
    ##         T.on[, i] <-  get_term(D[,i], I[,i], type = "on")
    ##     }
    ## } else {
        II[which(D==1)] <- 0 ## regard treated values as missing 
    ## }
    
    ## 1. remove units that have too few observations 
    T0 <- apply(II, 2, sum)
    T0.min <- min(T0)

    if (sum(T0[which(apply(D, 2, sum) > 0)] >= min.T0) == 0) {
        stop ("All treated units have been removed.\n")
    }   
    ## T0.min : minimum T0,  min.T0: manually set
    ## rm.tr.id: relative location of treated units (within all treated units) 
    ## that will be removed 
    if (T0.min < min.T0) {
        cat("Some treated units has too few pre-treatment periods; they are removed automatically.\n")
    }

    rm.id <- sort(unique(c(which((TT - apply(I, 2, sum)) > max.missing), which(T0 < min.T0))))
    ## rm.id <- which(T0 < min.T0) ## removed id
    ## rem.id <- which(T0 >= min.T0) ## remaining id  
    rem.id <- setdiff(1:N, rm.id)

    if (length(rm.id) == N) {
        stop("All units have been removed.\n")
    }

    if (length(rm.id) > 0) {
        X.old <- X
        if (p > 0) {
            X <- array(0,dim = c(TT, (N - length(rm.id)), p))
            for (i in 1:p) {
                subX <- X.old[, , i]
                X[, , i] <- as.matrix(subX[, -rm.id])
            }
        } else {
            X <- array(0,dim = c(TT, (N - length(rm.id)), 0))
        }

        # N <- N - length(rm.id)
        Y <- as.matrix(Y[,-rm.id])
        D <- as.matrix(D[,-rm.id])
        I <- as.matrix(I[,-rm.id]) ## after removing
        II <- as.matrix(II[,-rm.id])
        #if (!is.null(cl)) {
        #    cl <- cl[-rm.id]
        #}
    }  

    ## 2. check if some periods when all units are missing
    I.use <- apply(II, 1, sum) 
    if (0%in%I.use) {
        for (i in 1:TT) {
            if (I.use[i] == 0) {
                cat("There are not any observations under control at ",time.uni[i],", drop that period.\n")
            }
        }
        TT <- TT - sum(I.use == 0)
        time.uni <- time.uni[-which(I.use == 0)]
        
        I <- I[-which(I.use == 0),] ## remove that period
        II <- II[-which(I.use == 0),] ## remove that period        
        D <- D[-which(I.use == 0),] ## remove that period
        Y <- Y[-which(I.use == 0),] ## remove that period

        X.old <- X
        if (p > 0) {
            X <- array(0,dim = c(TT, (N - length(rm.id)), p))
            for (i in 1:p) {
                subX <- X.old[, , i]
                X[, , i] <- as.matrix(subX[-which(I.use == 0),])
            }
        } else {
            X <- array(0,dim = c(TT, (N - length(rm.id)), 0))
        }
    }

    ## 4. relative period 
    T.on <- matrix(NA, TT, (N - length(rm.id)))
    for (i in 1:(N - length(rm.id))) {
        T.on[, i] <-  get_term(D[, i], I[, i], type = "on")
    }

    ## 5. regard placebo period as under treatment
    if (placeboTest == TRUE) {
        if (length(placebo.period) == 1) {
            placebo.pos <- which(T.on == placebo.period)
            II[placebo.pos] <- 0
        } else {
            placebo.pos <- which(T.on >= placebo.period[1] & T.on <= placebo.period[2])
            II[placebo.pos] <- 0
        }

        ## remove treated units that have too few observations
        T0.2 <- apply(II, 2, sum)

        if (sum(T0.2[which(apply(D, 2, sum) > 0)] >= min.T0) == 0) {
            stop ("All treated units have been removed.\n")
        } 

        rm.id.2.pos <- sort(which(T0.2 < min.T0))
        rm.id.2 <- rem.id[rm.id.2.pos] 
        rem.id.2 <- setdiff(rem.id, rm.id.2)

        rem.id <- rem.id.2
        rm.id <- setdiff(1:N, rem.id)

        if (length(rm.id.2) > 0) {
            X.old <- X
            if (p > 0) {
                X <- array(0,dim = c(TT, (N - length(rm.id)), p))
                for (i in 1:p) {
                    subX <- X.old[, , i]
                    X[, , i] <- as.matrix(subX[, -rm.id.2.pos])
                }
            } else {
                X <- array(0,dim = c(TT, (N - length(rm.id)), 0))
            }

            # N <- N - length(rm.id)
            Y <- as.matrix(Y[,-rm.id.2.pos])
            D <- as.matrix(D[,-rm.id.2.pos])
            I <- as.matrix(I[,-rm.id.2.pos]) ## after removing
            II <- as.matrix(II[,-rm.id.2.pos])
            T.on <- as.matrix(T.on[,-rm.id.2.pos])
            #if (!is.null(cl)) {
            #    cl <- cl[-rm.id.2.pos]
            #}
        }  
    }

    ## check enough observations 
    if (min(apply(II, 1, sum)) == 0) {
        if (placeboTest == 0) {
            stop("Some periods do not have any observations.")
        } else {
            stop("Some periods do not have any observations. Please set a smaller range for placebo period.")
        }
    }

    if (min(apply(II, 2, sum)) == 0) {
        if (placeboTest == 0) {
            stop("Some units do not have any observations.")
        } else {
            stop("Some units do not have any observations. Please set a smaller range for placebo period.")
        }
    }

    ## 3. check reversals
    D.fake <- apply(D, 2, function(vec){cumsum(vec)})
    D.fake <- ifelse(D.fake > 0, 1, 0)
    D.fake[which(I==0)] <- 0
    Nrev <- sum(apply(D.fake == D, 2, sum) != TT)
    hasRevs <- ifelse(Nrev > 0, 1, 0)

    ## 6. switch-off periods  
    T.off <- NULL
    if (hasRevs == 1) {
        T.off <- matrix(NA, TT, (N - length(rm.id))) 
        for (i in 1:(N - length(rm.id))) {
            T.off[, i] <-  get_term(D[,i], I[,i], type = "off")
        }
    } 

    ## if (AR1 == TRUE) {
    ##     Y.first <- Y[1,]
    ##     Y.lag <- Y[1:(T-1),]
    ##     Y <- Y[2:T,]
    ##     D <- D[2:T,]
    ##     if (p == 0) {
    ##         X <- array(NA, dim=c((T-1),N,1))
    ##         X[,,1] <- Y.lag
    ##     } else {
    ##         X.first <- X[1,,]
    ##         X.sav <- X[2:T,,]
    ##         X <- array(NA,dim=c((T-1),N,(p+1)))
    ##         X[,,1] <- Y.lag
    ##         X[,,2:(p+1)] <- X.sav
    ##     }
    ##     T <- T-1
    ## } 

    ##-------------------------------##
    ## Register clusters
    ##-------------------------------##
    
    if ((se == TRUE | wald == TRUE) & parallel==TRUE) {

        ## set seed
        if (is.null(seed) == FALSE) {
            set.seed(seed)
        }
   
        if (is.null(cores) == TRUE) {
            cores <- detectCores()
        }
        para.clusters <- makeCluster(cores)
        registerDoParallel(para.clusters)
        cat("Parallel computing ...\n")
    }
    
    ##-------------------------------##
    ## run main program
    ##-------------------------------## 

    if (se == FALSE) {

        if (CV == TRUE) { 
            if (binary == FALSE) {
                out <- fect.cv(Y = Y, D = D, X = X, I = I, II = II, 
                               T.on = T.on, T.off = T.off, 
                               method = method,
                               criterion = criterion,
                               k = k, cv.prop = cv.prop,
                               cv.treat = cv.treat,
                               cv.continue = cv.continue, 
                               r = r, r.end = r.end, 
                               nlambda = nlambda, lambda = lambda,
                               force = force, hasRevs = hasRevs, 
                               tol = tol, norm.para = norm.para)
            } else {
                out <- fect.binary.cv(Y = Y, D = D, X = X,
                                      I = I, II = II, 
                                      T.on = T.on, T.off = T.off, 
                                      k = k, cv.prop = cv.prop,
                                      cv.treat = cv.treat, 
                                      cv.continue = cv.continue,
                                      r = r, r.end = r.end, 
                                      QR = QR, force = force, 
                                      hasRevs = hasRevs, tol = tol)
            }
            
        } else { ## non-binary case
            if (method == "ife") {
                out <- fect.fe(Y = Y, D = D, X = X, I = I, II = II,
                               T.on = T.on, T.off = T.off, r.cv = r,
                               binary = binary, QR = QR,
                               force = force, hasRevs = hasRevs, 
                               tol = tol, boot = 0,
                               norm.para = norm.para,
                               placeboTest = placeboTest, 
                               placebo.period = placebo.period)
            } else {
                out <- fect.mc(Y = Y, D = D, X = X, I = I, II = II,
                               T.on = T.on, T.off = T.off, 
                               lambda.cv = lambda,
                               force = force, hasRevs = hasRevs, 
                               tol = tol, boot = 0,
                               norm.para = norm.para,
                               placeboTest = placeboTest, 
                               placebo.period = placebo.period)
            }

        }

    } else { # SE == TRUE
        
        out <- fect.boot(Y = Y, D = D, X = X, I = I, II = II,
                         T.on = T.on, T.off = T.off, cl = NULL,
                         method = method, criterion = criterion,
                         CV = CV, k = k, r = r, r.end = r.end, 
                         nlambda = nlambda, lambda = lambda,
                         binary = binary, QR = QR,
                         force = force, hasRevs = hasRevs,
                         tol = tol, norm.para = norm.para,
                         placeboTest = placeboTest, 
                         placebo.period = placebo.period,
                         nboots = nboots, parallel = parallel,
                         cores = cores)

    }

    if ((out$validX == 0) & (p!=0) ) {
        warning("Multi-colinearity among covariates. Try removing some of them.\r")
    }     

    # Wald test
    if (wald == TRUE) {

        cat("Bootstrapping for Wald test ... ")

        out.wald <- fect.test(out = out, 
            Y = Y, D = D, X = X, I = I, II = II,
            T.on = T.on, T.off = NULL,
            method = method, cl = cl, r = out$r.cv, 
            lambda = out$lambda.cv,
            force = force, hasRevs = 0,
            tol = tol, norm.para = norm.para,
            pre.period = pre.period, nboots = nboots,
            parallel = parallel, cores = cores)
    } 

    
    if ((se == TRUE | wald == TRUE) & parallel == TRUE) {
        stopCluster(para.clusters)
        ##closeAllConnections()
    }

       
    
    ##-------------------------------##
    ## storage
    ##-------------------------------## 
    
    iname.old <- iname <- unique(sort(data.old[,id]))
    ## tname.old <- tname <- unique(sort(data.old[,time]))
    if (!0%in%I.use) {
        tname.old <- tname <- unique(sort(data.old[,time]))
    } else {
        tname.old <- tname <- unique(sort(data.old[,time]))[which(I.use != 0)]
    }

    if (length(rm.id)>0) {
        remove.id <- iname[rm.id]
        iname <- iname[-rm.id]
    }

    N.rem <- dim(D)[2]
    unit.type <- rep(NA, N.rem) ## 1 for control; 2 for treated; 3 for reversal;
    for (i in 1:N.rem) {
        di <- D[, i]
        ii <- I[, i]
        if (length(unique(di[which(ii==1)])) == 1) { ## treated or control
            if (0 %in% unique(di[which(ii==1)])) {
                unit.type[i] <- 1 ## control
            } else {
                unit.type[i] <- 2 ## treated
            }
        } else {
            unit.type[i] <- 3 ## reversal
        }
    }
    
    # 1 treated 2 control 3 missing 4 removed    
    obs.missing <- matrix(0, TT, N) ## not under treatment
    obs.missing[, rem.id] <- D + as.matrix(abs(I - 1)) * 3 ## under treatment
    obs.missing[which(obs.missing==0)] <- 2
    obs.missing[, rm.id] <- 4 ## removed

    colnames(obs.missing) <- unique(sort(data.old[,id]))
    rownames(obs.missing) <- tname
    
    if (p > 0) {
        Xname.tmp <- Xname
        rownames(out$beta) <- Xname.tmp
        if (se == TRUE) {
            rownames(out$est.beta) <- Xname.tmp
        }
    }  
    colnames(out$eff) <- iname
    rownames(out$eff) <- tname
   
    output <- c(list(Y.dat = Y,
                     D.dat = D,
                     I.dat = I,
                     Y = Yname,
                     D = Dname,
                     X = Xname,
                     index = index,
                     id = iname,
                     time = tname,
                     binary = binary,
                     placeboTest = placeboTest,
                     placebo.period = placebo.period,
                     unit.type = unit.type,
                     obs.missing = obs.missing), 
                     out)
                
    if (1 %in% rm.id) {
        output <- c(output,list(remove.id = remove.id))
        ## cat("list of removed units:",remove.id)
        ## cat("\n\n")
    }
    if (wald == TRUE) {
        output <- c(output,list(wald = out.wald))
    }
    output <- c(output, list(call = match.call()))
    class(output) <- "fect"
    return(output)
    
} ## Program fect ends 
