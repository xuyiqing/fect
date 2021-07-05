## Causal inference using counterfactual estimators 
## (fect: fixed effects counterfactuals)
## Version 0.1.0
## Author: Licheng Liu (Tsinghua), Ye Wang(NYU), Yiqing Xu(Stanford)
## Date: 2021.04.17

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
                 group = NULL, # cohort
                 na.rm = FALSE, # remove missing values
                 index, # c(unit, time) indicators
                 force = "two-way", # fixed effects demeaning
                 cl = "unit", 
                 r = 0, # number of factors
                 lambda = NULL, # mc method: regularization parameter
                 nlambda = 10, ## mc method: regularization parameter
                 CV = TRUE, # cross-validation
                 k = 10, # times of CV
                 cv.prop = 0.1, ## proportion of CV counts
                 cv.treat = TRUE, ## cv targeting treated units
                 cv.nobs = 3,  ## cv taking consecutive units
                 cv.donut = 0, ## cv mspe
                 binary = FALSE, # probit model
                 QR = FALSE, # QR or SVD for binary probit 
                 method = "fe", # method: e for fixed effects; ife for interactive fe; mc for matrix completion
                 criterion = "mspe", # for ife model: mspe, pc or both
                 alpha = 0.05, # significance level
                 se = FALSE, # report uncertainties
                 vartype = "bootstrap", # bootstrap or jackknife
                 nboots = 200, # number of bootstraps
                 parallel = TRUE, # parallel computing
                 cores = NULL, # number of cores
                 tol = 0.001, # tolerance level
                 seed = NULL, # set seed
                 min.T0 = 5, # minimum T0
                 max.missing = NULL, # maximum missing
                 proportion = 0.3, # use to fit the f test and equivalence test
                 pre.periods = NULL, # fit test period
                 f.threshold = 0.5, # equiv
                 tost.threshold = NULL, # equiv
                 knots = NULL,
                 degree = 2,  # wald = FALSE, # fit test
                 sfe = NULL,
                 cfe = NULL,
                 placeboTest = FALSE, # placebo test
                 placebo.period = NULL, # placebo test period
                 carryoverTest = FALSE, # carry-over test
                 carryover.period = NULL, # carry-over period
                 loo = FALSE, # leave one period out placebo  
                 permute = FALSE, ## permutation test
                 permu.dimension = 'time',
                 m = 2, ## block length
                 normalize = FALSE # accelerate option
                ) {
    UseMethod("fect")
}

## formula method

fect.formula <- function(formula = NULL,
                         data, # a data frame (long-form)
                         Y, # outcome
                         D, # treatment 
                         X = NULL, # time-varying covariates
                         group = NULL, # cohort
                         na.rm = FALSE, # remove missing values
                         index, # c(unit, time) indicators
                         force = "two-way", # fixed effects demeaning
                         cl = "unit", 
                         r = 0, # nubmer of factors
                         lambda = NULL, # mc method: regularization parameter
                         nlambda = 10, ## mc method: regularization parameter
                         CV = TRUE, # cross-validation
                         k = 10, # times of CV
                         cv.prop = 0.1, ## proportion of CV counts
                         cv.treat = TRUE, 
                         cv.nobs = 3,
                         cv.donut = 0, ## cv mspe
                         binary = FALSE, # probit model
                         QR = FALSE, # QR or SVD for binary probit 
                         method = "fe", # method: fe for fixed effects; ife for interactive fe; mc for matrix completion
                         criterion = "mspe", # for ife model: mspe, pc or both
                         alpha = 0.05, # significance level
                         se = FALSE, # report uncertainties
                         vartype = "bootstrap", # bootstrap or jackknife
                         nboots = 200, # number of bootstraps
                         parallel = TRUE, # parallel computing
                         cores = NULL, # number of cores
                         tol = 0.001, # tolerance level
                         seed = NULL, # set seed
                         min.T0 = 5,
                         max.missing = NULL,
                         proportion = 0.3,
                         pre.periods = NULL,
                         f.threshold = 0.5, # equiv
                         tost.threshold = NULL, 
                         knots = NULL,
                         degree = 2,   # wald = FALSE,
                         sfe = NULL,
                         cfe = NULL,
                         placeboTest = FALSE, # placebo test
                         placebo.period = NULL, # placebo test period
                         carryoverTest = FALSE, # carry-over test
                         carryover.period = NULL, # carry-over period
                         loo = FALSE, # leave one period out placebo
                         permute = FALSE, ## permutation test
                         permu.dimension = 'time',
                         m = 2, ## block length
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
    out <- fect.default(formula = NULL, 
                        data = data, 
                        Y = Yname,
                        D = Dname, 
                        X = Xname, 
                        group = group,
                        na.rm = na.rm, 
                        index = index, 
                        force = force, 
                        cl = cl, 
                        r = r, 
                        lambda = lambda, 
                        nlambda = nlambda, 
                        CV =CV, 
                        k = k, 
                        cv.prop = cv.prop, 
                        cv.treat = cv.treat, 
                        cv.nobs = cv.nobs, 
                        cv.donut = cv.donut,
                        binary = binary, 
                        QR = QR, 
                        method = method, 
                        criterion = criterion, 
                        alpha = alpha, 
                        se = se, 
                        vartype = vartype,
                        nboots = nboots, 
                        parallel = parallel, 
                        cores = cores, 
                        tol = tol, 
                        seed = seed, 
                        min.T0 = min.T0,
                        max.missing = max.missing, 
                        proportion = proportion, 
                        pre.periods = pre.periods, 
                        f.threshold = f.threshold, 
                        tost.threshold = tost.threshold,
                        knots = knots, 
                        degree = degree, 
                        sfe = sfe,
                        cfe = cfe,
                        placebo.period = placebo.period, 
                        placeboTest = placeboTest, 
                        carryoverTest = carryoverTest, 
                        carryover.period = carryover.period,
                        loo = loo,
                        permute = permute, 
                        permu.dimension = permu.dimension,
                        m = m, 
                        normalize = normalize)
    
    out$call <- match.call()
    out$formula <- formula
    print(out)
    return(out)

}


## default function

fect.default <- function(formula = NULL, data, # a data frame (long-form)
                         Y, # outcome
                         D, # treatment 
                         X = NULL, # time-varying covariates
                         group = NULL, # cohort
                         na.rm = FALSE, # remove missing values
                         index, # c(unit, time) indicators
                         force = "two-way", # fixed effects demeaning
                         cl = "unit", 
                         r = 0, # nubmer of factors
                         lambda = NULL, ## mc method: regularization parameter
                         nlambda = 0, ## mc method: regularization parameter
                         CV = TRUE, # cross-validation
                         k = 10, # times of CV
                         cv.prop = 0.1,
                         cv.treat = TRUE, 
                         cv.nobs = 3,
                         cv.donut = 1, ## cv mspe
                         binary = FALSE, # probit model
                         QR = FALSE, # QR or SVD for binary probit 
                         method = "fe", # method: ife for interactive fe; mc for matrix completion
                         criterion = "mspe", 
                         alpha = 0.05, # significance level
                         se = FALSE, # report uncertainties
                         vartype = "bootstrap", # bootstrap or jackknife
                         nboots = 200, # number of bootstraps
                         parallel = FALSE, # parallel computing
                         cores = NULL, # number of cores
                         tol = 0.001, # tolerance level
                         seed = NULL, # set seed
                         min.T0 = 5,
                         max.missing = NULL,
                         proportion = 0.3,
                         pre.periods = NULL,
                         f.threshold = 0.5, # equiv
                         tost.threshold = NULL, 
                         knots = NULL,
                         degree = 2,  # wald = FALSE,
                         sfe = NULL,
                         cfe = NULL,
                         placeboTest = FALSE, # placebo test
                         placebo.period = NULL, # placebo test period
                         carryoverTest = FALSE, # carry-over test
                         carryover.period = NULL, # carry-over period
                         loo = FALSE, # leave one period out placebo
                         permute = FALSE, ## permutation test
                         permu.dimension = "time", 
                         m = 2, ## block length
                         normalize = FALSE
                        ) {  
    
    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------## 
    placeboEquiv <- loo 
    ## read data 
    if (is.data.frame(data) == FALSE || length(class(data)) > 1) {
        data <- as.data.frame(data)
        ## warning("Not a data frame.")
    }
    ## index
    if (length(index) != 2 | sum(index %in% colnames(data)) != 2) {
        stop("\"index\" option misspecified. Try, for example, index = c(\"unit.id\", \"time\").")
    }

    if (se == 1) {
        if (! vartype %in% c("bootstrap", "jackknife", "parametric","wild")) {
            stop("\"vartype\" option misspecified.")
        }
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
    if (!method %in% c("fe", "ife", "mc", "both", "polynomial", "bspline", "cfe")) {
        stop("\"method\" option misspecified; choose from c(\"fe\", \"ife\", \"mc\", \"both\", \"polynomial\", \"bspline\",\"cfe\").")
    }
    if (method == "fe") {
        r <- 0
        CV <- FALSE
        method <- "ife"
    } else if (method %in% c("polynomial", "bspline","cfe")) {
        CV <- FALSE
    }

    if (!criterion %in% c("mspe", "wmspe", "gmspe", "moment", "gmoment", "mad", "pc")) {
        stop("\"criterion\" option misspecified; choose from c(\"mspe\", \"wmspe\", \"mad\", \"gmspe\",\"moment\",\"gmoment\", \"pc\").")
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

    ## leave one period out placebo 
    if (placeboEquiv == TRUE) {
        if(se!=TRUE){
            cat("For leave one period out placebo test, automatically set \"se\" to TRUE.")
            se <- TRUE
        }
        if(placeboTest == TRUE){
            cat("For leave one period out placebo test, automatically set \"placeboTest\" to FALSE.")
            placeboTest <- FALSE
        }
        if(carryoverTest == TRUE){
            cat("For leave one period out placebo test, automatically set \"carryoverTest\" to FALSE.")
            carryoverTest <- FALSE
        }
    }

    if(placeboTest==TRUE & carryoverTest==TRUE){
        stop("\"placeboTest\" and \"carryoverTest\" can't be performed simultaneously.")
    }

    ## CV
    if (method == "both") {
        CV <- TRUE
    }
    if (CV == TRUE) {
        
        if (placeboTest == TRUE) {
            stop("Placebo test cannot be performed while doing cross-validation.")
        }

        if (carryoverTest == TRUE) {
            stop("Carry-over test cannot be performed while doing cross-validation.")
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
    } 
    else {
        if (! method %in% c("ife", "mc", "polynomial", "bspline","cfe")) {
            stop("\"method\" option misspecified; please choose from c(\"ife\", \"mc\", \"polynomial\", \"bspline\").")
        }
    }

    if (method %in% c("polynomial", "bspline","cfe")) {
        if (permute == 1) {
            cat("Cannot do permutation test.\n")
            permute <- 0
        }
    }

    # wald <- FALSE 
    # if (se == TRUE && placeboTest == FALSE) {
    #     wald <- TRUE
    # }
    #if (placeboTest == TRUE) {
    #    wald <- FALSE
    #}
    if(permute == 1){
        if(placeboTest == TRUE){
            stop("\"permute\" can't be used with \"placeboTest\".")
        }
        if(carryoverTest == TRUE){
            stop("\"permute\" can't be used with \"carryoverTest\".")
        }
        if(loo == TRUE){
            stop("\"permute\" can't be used with \"loo\".")
        }
        if(!permu.dimension %in% c('time','unit')){
            stop("\"permu.dimension\" must be in c(\"unit\",\"time\") .")
        }
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

    # cohort 
    if (!is.null(group)) {
        if (! group %in% names(data)) {
            stop("\"group\" misspecified.\n")
        } 
    } 

    if(method == 'cfe'){
        if(is.null(sfe) & is.null(cfe)){
            cat("No additional sfe and cfe, use the \"fe\" estimator by default.\n")
            r <- 0
            CV <- FALSE
            method <- "ife"
        }
    }

    if(method == 'cfe'){
        if(!is.null(sfe)){
            for(sub.sfe in sfe){
                if (!sub.sfe %in% names(data)) {
                    stop("\"sfe\" misspecified.\n")
                }
                if(sub.sfe %in% index){
                    stop("\"sfe\" only contains additional fixed effects.\n")
                }    
            }            
        }

        if(!is.null(cfe)){
            if(!is.list(cfe)){
                stop("\"cfe\" should be a list.\n")
            }

            for(sub.cfe in cfe){
                for(sub.sub.cfe in sub.cfe){
                    if (!sub.sub.cfe %in% names(data)) {
                        stop("\"cfe\" misspecified.\n")
                    }
                }
                if(sub.cfe[1] == index[1] & sub.cfe[2] %in% X){
                    stop(paste0("Should remove ",sub.cfe[2]," from X.\n"))
                }
            }            
        }

    }

    if(method != 'cfe'){
        if (!is.null(group)) {
            data <- data[,c(index, Y, D, X, group)]
        } 
        else {
            data <- data[,c(index, Y, D, X)] ## some variables may not be used
        }        
    }
    else{
        all.var <- unique(c(index,sfe,unlist(cfe),Y,D,X,group))
        data <- data[,all.var]
    }

    
    if (na.rm == TRUE) {
        data <- na.omit(data)
    } 
    else{
        if(sum(is.na(data[,D]))>=1 | sum(is.na(data[,index[1]]))>=1  | sum(is.na(data[,index[2]]))>=1){
            stop("\"D\" or \"index\" should not have missing values when setting \"na.rm\" to FALSE.")
        }
        if(method=='cfe'){
            if(!is.null(sfe)){
                for(sub.sfe in sfe){
                    if(sum(is.na(data[,sub.sfe]))>=1){
                        stop("Variables in \"sfe\" should not have missing values when setting \"na.rm\" to FALSE.")
                    }
                }                
            }
        }
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

    ## carry-over period
    if (carryoverTest == TRUE && !is.null(carryover.period)) {
        if (sum(carryover.period <= 0) > 0) {
            stop("\"carryover.period\" should be larger than 0.")
        } else {
            if (length(carryover.period) > 2) {
                stop("\"carryover.period\" option misspecified. ")
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
        sd.Y <- sd(as.matrix(data[,Yname]),na.rm = TRUE)
        data[,c(Yname, Xname)] <- data[,c(Yname, Xname)]/sd.Y
        norm.para <- sd.Y ## normalized parameter
    }

    ## check index and treatment indicator
    if (! class(data[, Dname]) %in% c("numeric", "integer")) {
        ## data[, Dname] <- as.numeric(as.character(data[, Dname]))
        stop("Treatment indicator should be a numeric value.")
    } 

    if (class(data[, index[1]])[1] == "factor") {
        data[, index[1]] <- as.character(data[, index[1]])
    } 

    if (class(data[, index[2]])[1] == "factor") {
        data[, index[2]] <- as.character(data[, index[2]])
    } 
    
    id <- index[1]
    time <- index[2]
    TT.old <- TT <- length(unique(data[,time]))
    N.old <- N <- length(unique(data[,id]))
    p <- length(Xname)
    id.series <- unique(sort(data[,id])) ## unit id
    time.uni <- unique(sort(data[,time])) ## period

    if (!is.null(knots)) {
        for(i in 1:length(knots)) {
            knots[i] <- which(time.uni == knots[i])
        }
    }

    ## sort data
    data <- data[order(data[,id], data[,time]), ]
    if(na.rm==FALSE){
        # ziyi: if X and Y have missing values while D doesn't have missing values, save the value of D
        # check if some units or periods are completely missing after drop missing values of X or Y
        data.full <- data[,c(Dname,id,time)]
        data <- na.omit(data)
        data.old <- data
        TT <- length(unique(data[,time]))
        N <- length(unique(data[,id]))
        if(TT!=TT.old){
            cat("Some periods are totally removed after drop missing values of the outcome or covariates.\n")
        }
        if(N!=N.old){
            cat("Some units are totally removed after drop missing values of the outcome or covariates.\n")
        }
        id.series <- unique(sort(data[,id])) ## unit id
        time.uni <- unique(sort(data[,time])) ## period

        #remove these dropped units or periods in data.full
        #print(length(setdiff(data.full[,id],data[,id])))
        if(length(setdiff(data.full[,id],data[,id]))>0){
            rm.na.id <- setdiff(data.full[,id],data[,id])
            data.full <- data.full[which(!data.full[,id] %in% rm.na.id),]
        }else{
            rm.na.id <- NULL
        }

        if(length(setdiff(data.full[,time],data[,time]))>0){
            rm.na.time <- setdiff(data.full[,time],data[,time])
            data.full <- data.full[which(!data.full[,time] %in% rm.na.time),]
        }
        else{
            rm.na.time <- NULL
        }        
        # here the size of data.full should be smaller than TT*N, larger than length(data)
    }


    ## max.missing
    if (is.null(max.missing)) {
        max.missing <- TT
    }

    ## group.series <- NULL
    #group.ref <- NULL
    #if (!is.null(group)) {
    #    data.id <- data[, c(id, group)]
    #    rawgroup <- data.id[!duplicated(data.id[, 1]), 2]
    #    group <- as.numeric(as.factor(rawgroup))

        #group.ref <- cbind(group, rawgroup)
        #group.ref <- group.ref[!duplicated(group.ref[, 1]), ]
        #group.ref <- group.ref[order(group.ref[, 1]), ]
        ## rawgroup <- unique(rawgroup)
    #}

    ## gen group matrix
    if (!is.null(group)) {
        rawgroup <- data[, group]
        newgroup <- as.numeric(as.factor(rawgroup))
        data[, group] <- newgroup
        rawgroup <- cbind.data.frame(rawgroup, newgroup)
        rawgroup <- rawgroup[!duplicated(rawgroup[, 1]),]
    }

    if(method=='cfe'){
        if(!is.null(sfe)){
            for(sub.sfe in sfe){
                data[,sub.sfe] <- as.numeric(as.factor(data[,sub.sfe]))
            }        
        }

        if(!is.null(cfe)){
            for(sub.cfe in cfe){
                data[,sub.cfe[1]] <- as.numeric(as.factor(data[,sub.cfe[1]]))
            }        
        }        
    }

    ##cat("\nOK1\n")
    
    ## check missingness
    if (sum(is.na(data[, Yname])) > 0 & na.rm == TRUE) {
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
            if (sum(is.na(data[, Xname[i]])) > 0 & na.rm == TRUE) {
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
        if (!is.null(group)) {
            variable <- c(Yname, Dname, Xname, group)
        }
        if(method == 'cfe'){
            variable <- unique(c(sfe,unlist(cfe),variable))
        }

        data_I <- matrix(0, N * TT, 1)
        data_I[ob.indicator, 1] <- 1
        data_ub <- as.matrix(data[, variable])
        data <- data_ub_adj(data_I, data_ub)
        colnames(data) <- variable
        ## data is a TT*N matrix filled with observed pairs (Y/X, D).

        ## if these exists observations whose Y/X is missing but D is observed.
        if(na.rm==FALSE){
            data.full[,time] <- as.numeric(as.factor(data.full[,time]))
            ob.indicator.full <- data.full[,time]
            id.indicator.full <- table(data.full[, id])
            sub.start <- 1
            for (i in 1:(N - 1)) { 
                sub.start <- sub.start + id.indicator.full[i] 
                sub.end <- sub.start + id.indicator.full[i+1] - 1 
                ob.indicator.full[sub.start:sub.end] <- ob.indicator.full[sub.start:sub.end] + i * TT
            }
            data_I.full <- matrix(0, N * TT, 1)
            data_I.full[ob.indicator.full, 1] <- 1
            data_ub.full <- as.matrix(data.full[, Dname])
            data.D.full <- data_ub_adj(data_I.full, data_ub.full)
            colnames(data.D.full) <- Dname
            ## replace D in data with D in data.full
            data[,Dname] <- data.D.full[,Dname]
        }
    }

    ## indicator matrix: index matrix that indicates if data is observed 
    I.D <- I <- matrix(1, TT, N)
    Y.ind <- matrix(data[, Yname], TT, N)
    D.ind <- matrix(data[, Dname],TT,N)
    I[is.nan(Y.ind)] <- 0
    I.D[is.nan(D.ind)] <- 0
    ## I has more zeros than I.D 
    ## I.D is used in the function get_term

    if (0%in%I) {
        data[is.nan(data)] <- 0
    }

    ## group indicator 
    G.old <- G <- NULL
    if (!is.null(group)) {
        G <- matrix(data[, group], TT, N)
        ## replace group index 0(missing) for each unit
        if(!0 %in% G){
            if(sum(apply(G, 2, var))>0){
                stop("A unit in different periods should have the same group index.\n")
            }
        }
        else{
            if(max(apply(G,2,function(vec){return(length(table(vec)))}))>2){
                stop("A unit in different periods should have the same group index.\n")
            }
            G <- apply(G,2,function(vec){return(rep(max(vec),length(vec)))})
        }
        G.old <- G
    }
    
    ## each unit should have the same group index

    ## cat("\nOK2\n")

    ##treatment indicator: incorporates reversal treatments
    ## D==1 -> treatment
    D <- matrix(data[, Dname], TT, N)
    ##outcome variable
    Y <- matrix(data[, Yname], TT, N)
    ## time-varying covariates
    X <- array(0, dim = c(TT, N, p))
    if (p > 0) {
        for (i in 1:p) {
            X[,,i] <- matrix(data[, Xname[i]], TT, N)
        } 
    }

    index.matrix <- list()
    if(method=='cfe'){
        if(!is.null(sfe)){
            for(sub.sfe in sfe){
                data[,sub.sfe] <- as.numeric(as.factor(data[,sub.sfe]))
                sub.sfe.matrix <- matrix(data[,sub.sfe], TT, N)
                index.matrix[[sub.sfe]] <- sub.sfe.matrix
            }        
        }

        if(!is.null(cfe)){
            for(sub.cfe in cfe){
                sub.cfe.matrix <- matrix(as.numeric(as.factor(data[,sub.cfe[1]])), TT, N)
                index.matrix[[sub.cfe[1]]] <- sub.cfe.matrix
                sub.cfe.matrix <- matrix(data[,sub.cfe[2]], TT, N)
                index.matrix[[sub.cfe[2]]] <- sub.cfe.matrix
            }        
        }        
    }


    ## ----------------------------------------------------------- ##
    II <- I
    II[which(D==1)] <- 0 ## regard treated values as missing
    
    # Unbalance Check
    ## 1. remove units that have too control status
    T0 <- apply(II, 2, sum)
    T0.min <- min(T0)

    if (sum(T0[which(apply(D, 2, sum) > 0)] >= min.T0) == 0) {
        stop ("All treated units have been removed.\n")
    }   
    ## T0.min : minimum T0  
    ## min.T0: manually set
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
        } 
        else {
            X <- array(0,dim = c(TT, (N - length(rm.id)), 0))
        }

        # N <- N - length(rm.id)
        Y <- as.matrix(Y[,-rm.id])
        D <- as.matrix(D[,-rm.id])
        I <- as.matrix(I[,-rm.id]) ## after removing
        I.D <- as.matrix(I.D[,-rm.id])
        II <- as.matrix(II[,-rm.id])
        if (!is.null(group)) {
            G <- as.matrix(G[,-rm.id])
        }
        if(method == "cfe"){
            for(ind.name in names(index.matrix)){
                index.matrix[[ind.name]] <- as.matrix(index.matrix[[ind.name]][,-rm.id])
            }
        }
    }

    ## cat("\nOK1\n")  

    ## 2. check if some periods when all units are missing or treated
    I.use <- apply(II, 1, sum) 
    if (0%in%I.use) {
        for (i in 1:TT) {
            if (I.use[i] == 0) {
                cat("\nThere are not any observations under control at ",time.uni[i],", drop that period.\n")
            }
        }
        if (method %in% c("polynomial", "bspline")) {
            cat("\nThere are not any observations at some periods. Estimation results may not be reliable. Please use time fixed effects.\n")
        }
        TT <- TT - sum(I.use == 0)
        time.uni <- time.uni[-which(I.use == 0)]
        
        I <- I[-which(I.use == 0),] ## remove that period
        I.D <- I.D[-which(I.use == 0),]
        II <- II[-which(I.use == 0),] ## remove that period        
        D <- D[-which(I.use == 0),] ## remove that period
        Y <- Y[-which(I.use == 0),] ## remove that period

        if (!is.null(group)) {
            G <- G[-which(I.use == 0),]
        }

        if(method == "cfe"){
            for(ind.name in names(index.matrix)){
                index.matrix[[ind.name]] <- as.matrix(index.matrix[[ind.name]][-which(I.use == 0),])
            }
        }

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

    ## cat("\nOK2\n")  

    ## 3. relative period 
    T.on <- matrix(NA, TT, (N - length(rm.id)))
    for (i in 1:(N - length(rm.id))) {
        T.on[, i] <-  get_term(D[, i], I.D[, i], type = "on")
    }

    ## 4. check reversals
    D.fake <- apply(D, 2, function(vec){cumsum(vec)})
    D.fake <- ifelse(D.fake > 0, 1, 0)
    D.fake[which(I.D==0)] <- 0
    Nrev <- sum(apply(D.fake == D, 2, sum) != TT)
    hasRevs <- ifelse(Nrev > 0, 1, 0)
    if(hasRevs == FALSE & carryoverTest == TRUE){
        stop("Treatment status have no reversals. Cannot perform \"carryoverTest\" in this case.")
    }

    ## 5. switch-off periods
    T.off <- NULL
    if (hasRevs == 1) {
        T.off <- matrix(NA, TT, (N - length(rm.id))) 
        for (i in 1:(N - length(rm.id))) {
            T.off[, i] <-  get_term(D[,i], I.D[,i], type = "off")
        }
    }

    ## 6. regard placebo period as under treatment
    if (placeboTest == TRUE) {
        II.origin <- II
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
            stop ("All treated units have been removed in placebo test.\n")
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
            I.D <- as.matrix(I.D[,-rm.id.2.pos])
            II <- as.matrix(II[,-rm.id.2.pos])
            II.origin <- as.matrix(II.origin[,-rm.id.2.pos])
            T.on <- as.matrix(T.on[,-rm.id.2.pos])
            if(hasRevs){
                T.off <- as.matrix(T.off[,-rm.id.2.pos])                
            }
            if (!is.null(group)) {
                G <- as.matrix(G[,-rm.id.2.pos])
            }
            if(method == "cfe"){
                for(ind.name in names(index.matrix)){
                    index.matrix[[ind.name]] <- as.matrix(index.matrix[[ind.name]][,-rm.id.2.pos])
                }
            }
        }  
    }

    ## 7. Carryover Test 
    ## testcarryover.period = c(1,3)
    if(hasRevs == 1 & carryoverTest==TRUE & is.null(carryover.period)==FALSE){
        II.origin <- II
        if (length(carryover.period) == 1) {
            carryover.pos <- which(T.off == carryover.period)
            II[carryover.pos] <- 0
        } else {
            carryover.pos <- which(T.off >= carryover.period[1] & T.off <= carryover.period[2])
            II[carryover.pos] <- 0
        }
        ## remove treated units that have too few observations
        T0.3 <- apply(II, 2, sum)

        if (sum(T0.3[which(apply(D, 2, sum) > 0)] >= min.T0) == 0) {
            stop ("All treated units have been removed in carryover test.\n")
        } 

        rm.id.3.pos <- sort(which(T0.3 < min.T0))
        rm.id.3 <- rem.id[rm.id.3.pos] 
        rem.id.3 <- setdiff(rem.id, rm.id.3)

        rem.id <- rem.id.3
        rm.id <- setdiff(1:N, rem.id)

        if (length(rm.id.3) > 0) {
            X.old <- X
            if (p > 0) {
                X <- array(0,dim = c(TT, (N - length(rm.id)), p))
                for (i in 1:p) {
                    subX <- X.old[, , i]
                    X[, , i] <- as.matrix(subX[, -rm.id.3.pos])
                }
            } else {
                X <- array(0,dim = c(TT, (N - length(rm.id)), 0))
            }

            # N <- N - length(rm.id)
            Y <- as.matrix(Y[,-rm.id.3.pos])
            D <- as.matrix(D[,-rm.id.3.pos])
            I <- as.matrix(I[,-rm.id.3.pos]) ## after removing
            I.D <- as.matrix(I.D[,-rm.id.3.pos])
            II <- as.matrix(II[,-rm.id.3.pos])
            II.origin <- as.matrix(II.origin[,-rm.id.3.pos])
            T.on <- as.matrix(T.on[,-rm.id.3.pos])
            if(hasRevs){
                T.off <- as.matrix(T.off[,-rm.id.3.pos])                
            }
            if (!is.null(group)) {
                G <- as.matrix(G[,-rm.id.3.pos])
            }
            if(method == "cfe"){
                for(ind.name in names(index.matrix)){
                    index.matrix[[ind.name]] <- as.matrix(index.matrix[[ind.name]][,-rm.id.3.pos])
                }
            }
        }
    }

    ## 8. Finally, check enough observations 
    if (min(apply(II, 1, sum)) == 0) {
        if (placeboTest == 1) {
            stop("Some periods do not have any observations. Please set a smaller range for placebo period.")
        } 
        else if(carryoverTest == 1) {
            stop("Some periods do not have any observations. Please set a smaller range for carryover period.")
        }
        else{
            stop("Some periods do not have any observations.")
        }
    }

    if (min(apply(II, 2, sum)) == 0) {
       if (placeboTest == 1) {
            stop("Some units do not have any observations. Please set a smaller range for placebo period.")
        } 
        else if(carryoverTest == 1) {
            stop("Some units do not have any observations. Please set a smaller range for carryover period.")
        }
        else{
            stop("Some units do not have any observations.")
        }
    }
    ## cohort
    g.level <- NULL
    if (!is.null(group)) {
        #G[which(D == 0)] <- NA
        g.level <- unique(c(G))
        g.level <- g.level[!is.na(g.level)]
        rownames(rawgroup) <- rawgroup[,'newgroup']
        names(g.level) <- rawgroup[as.character(g.level),'rawgroup']
        #rawgroup <- rawgroup[order(rawgroup[, 2]),]
        #rawgroup <- rawgroup[which(rawgroup[, 2] %in% g.level), 1]
        ## tr.pos <- which(apply(D, 2, sum) > 0)
        ## rawgroup <- unique(rawgroup[tr.pos])
        ## ng <- length(group)
        ## group <- matrix(rep(group, each = TT), TT, ng)
        ## G <- as.matrix(G[, tr.pos])
    }

    ##cat("\nOK3\n")


    ## cat("\nOK3\n")

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

    if((se == TRUE | permute == TRUE) & parallel==FALSE){
        ## set seed
        if (is.null(seed) == FALSE) {
            set.seed(seed)
        } 
    }
    
    if ((se == TRUE | permute == TRUE) & parallel==TRUE) {
        ## set seed
        if (is.null(seed) == FALSE) {
            set.seed(seed)
        }
        if (is.null(cores) == TRUE) {
            cores <- detectCores()
        }
        para.clusters <- future::makeClusterPSOCK(cores)
        registerDoParallel(para.clusters)
        if (is.null(seed) == FALSE) {
            registerDoRNG(seed)
        }
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
                               cv.nobs = cv.nobs, 
                               cv.donut = cv.donut,
                               min.T0 = min.T0,
                               r = r, r.end = r.end,
                               proportion = proportion, 
                               nlambda = nlambda, lambda = lambda,
                               force = force, hasRevs = hasRevs, 
                               tol = tol, norm.para = norm.para, 
                               group.level = g.level, group = G)
            } 
            else {
                out <- fect.binary.cv(Y = Y, D = D, X = X,
                                      I = I, II = II, 
                                      T.on = T.on, T.off = T.off, 
                                      k = k, cv.prop = cv.prop,
                                      cv.treat = cv.treat, 
                                      cv.nobs = cv.nobs,
                                      r = r, r.end = r.end, 
                                      QR = QR, force = force, 
                                      hasRevs = hasRevs, tol = tol,
                                      group.level = g.level, group = G)
            }
            
        } 
        else { ## non-binary case
            if (method == "ife") {
                out <- try(fect.fe(Y = Y, D = D, X = X, I = I, II = II,
                               T.on = T.on, T.off = T.off, r.cv = r,
                               binary = binary, QR = QR,
                               force = force, hasRevs = hasRevs, 
                               tol = tol, boot = 0,
                               norm.para = norm.para,
                               placeboTest = placeboTest, 
                               placebo.period = placebo.period,
                               carryoverTest = carryoverTest,
                               carryover.period = carryover.period,
                               group.level = g.level, group = G), silent = TRUE)
            } 
            else if (method == "mc") {
                out <- try(fect.mc(Y = Y, D = D, X = X, I = I, II = II,
                               T.on = T.on, T.off = T.off, 
                               lambda.cv = lambda,
                               force = force, hasRevs = hasRevs, 
                               tol = tol, boot = 0,
                               norm.para = norm.para,
                               placeboTest = placeboTest, 
                               placebo.period = placebo.period,
                               carryoverTest = carryoverTest,
                               carryover.period = carryover.period,
                               group.level = g.level, group = G), silent = TRUE)
            } 
            else if (method %in% c("polynomial", "bspline", "cfe")) {
                out <- fect.polynomial(Y = Y, D = D, X = X, I = I, 
                                       II = II, T.on = T.on, 
                                       T.off = T.off, method = method,
                                       degree = degree,
                                       sfe = sfe, cfe = cfe,
                                       ind.matrix = index.matrix,
                                       knots = knots, force = force, 
                                       hasRevs = hasRevs, tol = tol, boot = 0, 
                                       placeboTest = placeboTest,
                                       placebo.period = placebo.period, 
                                       carryoverTest = carryoverTest,
                                       carryover.period = carryover.period,
                                       norm.para = norm.para,
                                       group.level = g.level, group = G)


            }

            if ('try-error' %in% class(out)) {
                stop("\nCannot estimate.\n")
            }
            # only for polynomial methods
            if(method %in% c("polynomial", "bspline", "cfe")){
                I <- out$I
                II <- out$II
            }
                
        }
    } 
    else { # SE == TRUE
        
        out <- fect.boot(Y = Y, D = D, X = X, I = I, II = II,
                         T.on = T.on, T.off = T.off, cl = NULL,
                         method = method, degree = degree,
                         sfe = sfe, cfe = cfe,
                         ind.matrix = index.matrix,
                         knots = knots, criterion = criterion,
                         CV = CV, k = k, cv.prop = cv.prop,
                         cv.treat = cv.treat, cv.nobs = cv.nobs,
                         r = r, r.end = r.end, 
                         nlambda = nlambda, lambda = lambda,
                         alpha = alpha, binary = binary, QR = QR,
                         force = force, hasRevs = hasRevs,
                         tol = tol, norm.para = norm.para,
                         placeboTest = placeboTest, 
                         placebo.period = placebo.period,
                         carryoverTest = carryoverTest,
                         carryover.period = carryover.period,
                         vartype = vartype,
                         nboots = nboots, parallel = parallel,
                         cores = cores, group.level = g.level, group = G)

        if(method %in% c("polynomial", "bspline", "cfe")){
            I <- out$I
            II <- out$II
        }
    
    }



    if ((out$validX == 0) & (p!=0) ) {
        warning("Multi-colinearity among covariates. Try removing some of them.\r")
    }

    pre.est.att <- pre.att.bound <- NULL
    pre.term <- NULL
    N_bar <- NULL
    ## leave one period out placebo test for pre-treatment periods

    if (is.null(proportion)==TRUE) {
        proportion <- 0    
    }
    max.count <- max(out$count)
    max.pre.periods <- out$time[which(out$count >= max.count * proportion & out$time <= 0)]
    if (is.null(pre.periods) == TRUE) {        
        pre.periods <- max.pre.periods        
    } 
    else {
        pre.periods <- intersect(pre.periods[1]:pre.periods[length(pre.periods)], max.pre.periods)
    }   
    pre.term <- pre.periods
    N_bar <- max(out$count[which(out$time %in% pre.periods)])
      
    if (placeboEquiv == TRUE) {

        r.cv <- out$r.cv 
        lambda.cv <- out$lambda.cv 
        method <- out$method
        if (method == "fe") {
            method <- "ife"
        }

        cat("\nOut-of-Sample Test...\n")
        #if (is.null(placebo.period)) {
        #    pre.term <- pre.term.min:0
        #} else {
        #    pre.term.bound <- sort(placebo.period)
        #    pre.term <- pre.term.bound[1]:pre.term.bound[length(pre.term.bound)]
        #} 

        placebo.period <- pre.term[1]:pre.term[length(pre.term)] 

        pre.est.att <- matrix(NA, length(pre.term), 6)
        pre.att.bound <- matrix(NA, length(pre.term), 2)
        pre.att.boot <- matrix(NA, length(pre.term), nboots)

        rownames(pre.est.att) <- rownames(pre.att.bound) <- pre.term
        colnames(pre.est.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                  "p.value", "count.on")
        colnames(pre.att.bound) <- c("CI.lower", "CI.upper")

        if(!is.null(group)){
            all_group_names <- names(g.level)
            pre.est.att.group <- list()
            pre.att.bound.group <- list()
            pre.att.boot.group <- list()

            #only use the selected pre.periods before
            for(group.name in all_group_names){
                pre.term.sub <- intersect(out$group.output[[group.name]]$time.on,pre.term)
                pre.est.att.sub <- matrix(NA, length(pre.term.sub), 6)
                pre.att.bound.sub <- matrix(NA, length(pre.term.sub), 2)
                pre.att.boot.sub <- matrix(NA, length(pre.term.sub), nboots)
                if(length(pre.term.sub)>0){
                    rownames(pre.att.boot.sub) <- rownames(pre.est.att.sub) <- rownames(pre.att.bound.sub) <- pre.term.sub
                }
                colnames(pre.est.att.sub) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                           "p.value", "count.on")
                colnames(pre.att.bound.sub) <- c("CI.lower", "CI.upper")

                pre.est.att.group[[group.name]] <- pre.est.att.sub
                pre.att.bound.group[[group.name]] <- pre.att.bound.sub
                pre.att.boot.group[[group.name]] <- pre.att.boot.sub
            }
        }

        jj <- length(pre.term)
        pre.term <- sort(pre.term, decreasing = TRUE)

        for (kk in pre.term) {
            
            placebo.pos <- which(T.on == kk)

            pX <- X 
            pY <- Y 
            pD <- D 
            pI <- I 
            pII <- II 
            pT.on <- T.on
            pT.off <- T.off
            pG <- G

            pII[placebo.pos] <- 0 ## placebo treatment


            ## remove treated units that have too few observations
            T0.2 <- apply(pII, 2, sum)

            if (sum(T0.2[which(apply(D, 2, sum) > 0)] >= min.T0) == 0) {
                cat("\n")
                cat(paste("All treated units have been removed for period ", kk, sep = ""))
                cat("\n")
                jj <- jj - 1
            
            } else {

                te <- paste("Pre-period ", kk, sep = "")
                if (kk == 0) {
                    te <- paste(te, "(one period before treatment)", sep = " ")
                }
                cat("\n")
                cat(te)
                cat("\n")

                rem.id.new <- rem.id
                rm.id.2.pos <- sort(which(T0.2 < min.T0))
                rm.id.2 <- rem.id.new[rm.id.2.pos] 
                rem.id.2 <- setdiff(rem.id.new, rm.id.2)

                rem.id.new <- rem.id.2
                rm.id.new <- setdiff(1:N, rem.id.new)
                 
                if (length(rm.id.2) > 0) {
                    X.old <- pX
                    if (p > 0) {
                        pX <- array(0,dim = c(TT, (N - length(rm.id.new)), p))
                        for (i in 1:p) {
                            subX <- X.old[, , i]
                            pX[, , i] <- as.matrix(subX[, -rm.id.2.pos])
                        }
                    } else {
                        pX <- array(0,dim = c(TT, (N - length(rm.id.new)), 0))
                    }

                    # N <- N - length(rm.id)
                    pY <- as.matrix(Y[,-rm.id.2.pos])
                    pD <- as.matrix(D[,-rm.id.2.pos])
                    pI <- as.matrix(I[,-rm.id.2.pos]) ## after removing
                    pII <- as.matrix(II[,-rm.id.2.pos])
                    pT.on <- as.matrix(T.on[,-rm.id.2.pos])
                    #if (!is.null(cl)) {
                    #    cl <- cl[-rm.id.2.pos]
                    #}
                    if(hasRevs){
                        pT.off <- as.matrix(T.off[,-rm.id.2.pos])
                    }
                    if (!is.null(group)) {
                        ## group <- group[-rm.id.2.pos]
                        ## rawgroup <- rawgroup[-rm.id.2.pos]
                        pG <- as.matrix(G[,-rm.id.2.pos])
                    }
                }

                p.out <- fect.boot(Y = pY, D = pD, X = pX, I = pI, II = pII,
                             T.on = pT.on, T.off = pT.off, cl = NULL,
                             method = method, degree = degree,
                             knots = knots, criterion = criterion,
                             CV = 0, k = k, cv.prop = cv.prop,
                             cv.treat = cv.treat, cv.nobs = cv.nobs,
                             r = r.cv, r.end = r.end, 
                             nlambda = nlambda, lambda = lambda.cv,
                             alpha = alpha, binary = binary, QR = QR,
                             force = force, hasRevs = 0,
                             tol = tol, norm.para = norm.para,
                             placeboTest = 0, 
                             placebo.period = NULL,
                             carryoverTest = 0,
                             carryover.period = NULL,
                             vartype = vartype,
                             nboots = nboots, parallel = parallel,
                             cores = cores, group.level = g.level, group = pG, 
                             dis = FALSE)

                p.est.att <- p.out$est.att 
                p.att.bound <- p.out$att.bound 
                p.pos <- which(as.numeric(rownames(p.est.att)) == kk)

                pre.est.att[jj, ] <- p.est.att[p.pos, ]
                pre.att.bound[jj, ] <- p.att.bound[p.pos, ]
                pre.att.boot[jj, ] <- p.out$att.boot[p.pos, ]
                pre.period.name <- rownames(pre.est.att)[jj]


                if(!is.null(group)){
                    for(group.name in all_group_names){
                        check.sub <- p.out$group.output[[group.name]]$att.on
                        if(!is.null(check.sub)){
                            p.est.att.sub <- p.out$est.group.output[[group.name]]$att.on
                            p.att.bound.sub <- p.out$est.group.output[[group.name]]$att.on.bound
                            p.att.boot.sub <- p.out$est.group.output[[group.name]]$att.on.boot
                            p.pos.sub <- which(as.numeric(rownames(p.est.att.sub)) == kk)
                            if(length(p.pos.sub)==1){
                                pre.est.att.group[[group.name]][pre.period.name,] <- p.est.att.sub[p.pos.sub,]
                                pre.att.bound.group[[group.name]][pre.period.name,] <- p.att.bound.sub[p.pos.sub,]
                                pre.att.boot.group[[group.name]][pre.period.name,] <- p.att.boot.sub[p.pos.sub,] 
                            }
                        }
                    }
                }
                jj <- jj - 1
            }
        }

        pre.est.group.output <- NULL
        if(!is.null(group)){
            pre.est.group.output <- list()
            for(group.name in all_group_names){
                sub.output <- list()
                sub.output$pre.est.att <- pre.est.att.group[[group.name]]
                sub.output$pre.att.bound <- pre.att.bound.group[[group.name]]
                sub.output$pre.att.boot <- pre.att.boot.group[[group.name]]
                pre.est.group.output[[group.name]] <- sub.output
            }
        }           
    }


    ## permutation test 
    if (permute == TRUE) {
        cat("Permuting under sharp null hypothesis ... ")

        out.permute <- fect.permu(Y = Y, X = X, D = D, I = I, r.cv = out$r.cv,
                                  lambda.cv = out$lambda.cv, m = m, 
                                  permu.dimension = permu.dimension,
                                  method = out$method, degree = degree, 
                                  knots = knots, force = force,                      
                                  tol = tol, norm.para = norm.para,
                                  nboots = nboots,
                                  parallel = parallel, cores = cores)

        permute.p <- sum(out.permute > abs(out$att.avg)) / length(out.permute)

        permute.result <- list(permute.att.avg = out.permute, p = permute.p)

    }

    
    if ((se == TRUE | permute) & parallel == TRUE) {
        stopCluster(para.clusters)
        ##closeAllConnections()
    }

    ## cat("\nOK4\n")

       
    
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
    
    # 1 treated 
    # 2 control 
    # 3 missing 
    # 4 removed    
    obs.missing <- matrix(0, TT, N) ## not under treatment
    obs.missing[, rem.id] <- D + as.matrix(abs(I - 1)) * 3 ## under treatment
    obs.missing[which(obs.missing==0)] <- 2
    if(placeboTest|carryoverTest){
        obs.missing[, rem.id] <- obs.missing[, rem.id] + 3*(II.origin!=II) ##placebo or carryover
    }
    obs.missing[which(obs.missing==4)] <- 3 # in case if I.D!=I
    obs.missing[, rm.id] <- 4 ## removed

    colnames(obs.missing) <- unique(sort(data.old[,id]))
    rownames(obs.missing) <- tname

    # if cross-validation:

    if (p > 0) {
        Xname.tmp <- Xname
        rownames(out$beta) <- Xname.tmp
        colnames(out$beta) <- c("Coef")
        if (binary == TRUE) {
            rownames(out$marginal) <- Xname.tmp
        }
        if (se == TRUE) {
            rownames(out$est.beta) <- Xname.tmp
            colnames(out$est.beta) <- c("Coef","S.E.","CI.lower","CI.upper","p.value")
            if (binary == TRUE) {
                rownames(out$est.marginal) <- Xname.tmp
            }
            if (placeboTest == TRUE) {
                colnames(out$est.placebo) <- c("Coef","S.E.","CI.lower","CI.upper","p.value","CI.lower(90%)","CI.upper(90%)")
                rownames(out$est.placebo) <- c("Placebo effect")
            }
            if (carryoverTest == TRUE) {
                colnames(out$est.carryover) <- c("Coef","S.E.","CI.lower","CI.upper","p.value","CI.lower(90%)","CI.upper(90%)")
                rownames(out$est.carryover) <- c("Carryover effect")
            }
        }
    }  
    colnames(out$eff) <- iname
    rownames(out$eff) <- tname

    ## cohort effect
    if (!is.null(group)) {
        out$group <- rawgroup
        #out$G <- group
        #out$group2 <- rawgroup2
        if (se == 1) {
            rownames(out$est.group.att) <- names(g.level)
        }
        out$g.level <- g.level
    }

    if (is.null(tost.threshold)==TRUE) {
        tost.threshold <- 0.36 * sqrt(out$sigma2.fect)
    }

    output <- c(list(Y.dat = Y,
                     D.dat = D,
                     I.dat = I,
                     Y = Yname,
                     D = Dname,
                     X = Xname,
                     T.on = T.on,
                     G = G.old,
                     hasRevs = hasRevs,
                     T.off = T.off,
                     index = index,
                     id = iname,
                     rawtime = tname,
                     binary = binary,
                     loo = loo,
                     proportion = proportion,
                     pre.periods = pre.periods,
                     tost.threshold = tost.threshold,
                     placeboTest = placeboTest,
                     placebo.period = placebo.period,
                     carryoverTest = carryoverTest,
                     carryover.period = carryover.period,
                     unit.type = unit.type,
                     obs.missing = obs.missing), 
                     out)
    
                
    if (1 %in% rm.id) {
        output <- c(output,list(remove.id = remove.id))
        ## cat("list of removed units:",remove.id)
        ## cat("\n\n")
    }    
    
    #if (se == TRUE) { 
    #    suppressWarnings(test.out <- diagtest(output, pre.periods = pre.periods, 
    #        f.threshold = f.threshold, tost.threshold = tost.threshold))
    #    output <- c(output, list(test.out = test.out))
    #}
    
    
    if (permute == TRUE) {
        output <- c(output,list(permute = permute.result))
    }

    if (placeboEquiv == TRUE) {
        output <- c(output, list(pre.est.att = pre.est.att, 
                                 pre.att.bound = pre.att.bound, 
                                 pre.att.boot = pre.att.boot,
                                 pre.est.group.output = pre.est.group.output))
    }

    # if (placeboEquiv || placeboTest || carryoverTest) {
    # classic equivalence test, placeboTest, and carryoverTest
    # this can also be used in placeboTest
    
    # 👇this is the classic equivalence test
    if(loo==TRUE){
        output$loo <- FALSE
    }

    if(se==1){
        suppressWarnings(
            test.out <- diagtest(output, pre.periods = pre.periods, 
                    f.threshold = f.threshold, 
                    tost.threshold = tost.threshold, 
                    N_bar = N_bar)
        )
        output <- c(output, list(test.out = test.out))
    }


    # 👇this is the loo equivalence test
    if(loo==TRUE){
        output$loo <- TRUE
    }
    if(loo==TRUE && se == 1){
        suppressWarnings(
        test.out <- diagtest(output, pre.periods = pre.periods, 
                    f.threshold = f.threshold, 
                    tost.threshold = tost.threshold, 
                    N_bar = N_bar)
        )
        output <- c(output, list(loo.test.out = test.out))
    }


    output <- c(output, list(call = match.call()))
    class(output) <- "fect"
    return(output)
} ## Program fect ends 





