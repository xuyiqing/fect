## Causal inference using counterfactual estimators
## (fect: fixed effects counterfactuals)

## MAIN FUNCTION
## fect.formula()
## fect.default()

## DEPENDENT FUNCTIONS
## fect_fe() ## interactive fixed effects model
## fect_mc() ## matrix completion
## fect_boot() ## bootstrap

## fitness test
## fect_test ## wild bootstrap

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
                 W = NULL, # weight
                 group = NULL, # cohort
                 na.rm = FALSE, # remove missing values
                 index, # c(unit, time) indicators
                 force = "two-way", # fixed effects demeaning
                 r = 0, # number of factors
                 lambda = NULL, # mc method: regularization parameter
                 nlambda = 10, ## mc method: regularization parameter
                 CV = NULL, # cross-validation
                 k = 10, # times of CV
                 cv.prop = 0.1, ## proportion of CV counts
                 cv.treat = FALSE, ## cv targeting treated units
                 cv.nobs = 3,  ## cv taking consecutive units
                 cv.donut = 0, ## cv mspe
                 criterion = "mspe", # for ife model: mspe, pc or both
                 binary = FALSE, # probit model
                 QR = FALSE, # QR or SVD for binary probit
                 method = "fe", # method: e for fixed effects; ife for interactive fe; mc for matrix completion
                 se = FALSE, # report uncertainties
                 vartype = "bootstrap", # bootstrap or jackknife
                 cl = NULL,
                 quantile.CI = FALSE,
                 nboots = 200, # number of bootstraps
                 alpha = 0.05, # significance level
                 parallel = TRUE, # parallel computing
                 cores = NULL, # number of cores
                 tol = 1e-3, # tolerance level
                 max.iteration = 1000,
                 seed = NULL, # set seed
                 min.T0 = NULL, # minimum T0
                 max.missing = NULL, # maximum missing
                 proportion = 0.3, # use to fit the f test and equivalence test
                 pre.periods = NULL, # fit test period
                 f.threshold = 0.5, # equiv
                 tost.threshold = NULL, # equiv
                 knots = NULL,
                 degree = 2,  # wald = FALSE, # fit test
                 sfe = NULL,
                 cfe = NULL,
                 balance.period = NULL, # the pre and post periods for balanced samples
                 fill.missing = FALSE, # whether to balance missing observations
                 placeboTest = FALSE, # placebo test
                 placebo.period = NULL, # placebo test period
                 carryoverTest = FALSE, # carry-over test
                 carryover.period = NULL, # carry-over period
                 carryover.rm = NULL,
                 loo = FALSE, # leave one period out placebo
                 permute = FALSE, ## permutation test
                 m = 2, ## block length
                 normalize = FALSE, # accelerate option
                 keep.sims = FALSE # keep individual bootstrap/jackknife simulations
                ) {
    UseMethod("fect")
}

## formula method

fect.formula <- function(formula = NULL,
                         data, # a data frame (long-form)
                         Y, # outcome
                         D, # treatment
                         X = NULL, # time-varying covariates
                         W = NULL, # weights
                         group = NULL, # cohort
                         na.rm = FALSE, # remove missing values
                         index, # c(unit, time) indicators
                         force = "two-way", # fixed effects demeaning
                         r = 0, # nubmer of factors
                         lambda = NULL, # mc method: regularization parameter
                         nlambda = 10, ## mc method: regularization parameter
                         CV = NULL, # cross-validation
                         k = 10, # times of CV
                         cv.prop = 0.1, ## proportion of CV counts
                         cv.treat = FALSE,
                         cv.nobs = 3,
                         cv.donut = 0, ## cv mspe
                         criterion = "mspe", # for ife model: mspe, pc or both
                         binary = FALSE, # probit model
                         QR = FALSE, # QR or SVD for binary probit
                         method = "fe", # method: fe for fixed effects; ife for interactive fe; mc for matrix completion
                         se = FALSE, # report uncertainties
                         vartype = "bootstrap", # bootstrap or jackknife
                         cl = NULL,
                         quantile.CI = FALSE,
                         nboots = 200, # number of bootstraps
                         alpha = 0.05, # significance level
                         parallel = TRUE, # parallel computing
                         cores = NULL, # number of cores
                         tol = 1e-3, # tolerance level
                         max.iteration = 1000,
                         seed = NULL, # set seed
                         min.T0 = NULL,
                         max.missing = NULL,
                         proportion = 0.3,
                         pre.periods = NULL,
                         f.threshold = 0.5, # equiv
                         tost.threshold = NULL,
                         knots = NULL,
                         degree = 2,   # wald = FALSE,
                         sfe = NULL,
                         cfe = NULL,
                         balance.period = NULL, # the pre and post periods for balanced samples
                         fill.missing = FALSE, # whether to balance missing observations
                         placeboTest = FALSE, # placebo test
                         placebo.period = NULL, # placebo test period
                         carryoverTest = FALSE, # carry-over test
                         carryover.period = NULL, # carry-over period
                         carryover.rm = NULL,
                         loo = FALSE, # leave one period out placebo
                         permute = FALSE, ## permutation test
                         m = 2, ## block length
                         normalize = FALSE,
                         keep.sims = FALSE
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
                        W = W,
                        group = group,
                        na.rm = na.rm,
                        balance.period = balance.period,
                        fill.missing = fill.missing,
                        index = index,
                        force = force,
                        r = r,
                        lambda = lambda,
                        nlambda = nlambda,
                        CV =CV,
                        k = k,
                        cv.prop = cv.prop,
                        cv.treat = cv.treat,
                        cv.nobs = cv.nobs,
                        cv.donut = cv.donut,
                        criterion = criterion,
                        binary = binary,
                        QR = QR,
                        method = method,
                        se = se,
                        vartype = vartype,
                        cl = cl,
                        quantile.CI = quantile.CI,
                        nboots = nboots,
                        alpha = alpha,
                        parallel = parallel,
                        cores = cores,
                        tol = tol,
                        max.iteration = max.iteration,
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
                        carryover.rm = carryover.rm,
                        loo = loo,
                        permute = permute,
                        m = m,
                        normalize = normalize,
                        keep.sims = keep.sims)

    out$call <- match.call()
    out$formula <- formula
    return(out)

}


## default function

fect.default <- function(formula = NULL, data, # a data frame (long-form)
                         Y, # outcome
                         D, # treatment
                         X = NULL, # time-varying covariates
                         W = NULL, # weights
                         group = NULL, # cohort
                         na.rm = FALSE, # remove missing values
                         index, # c(unit, time) indicators
                         force = "two-way", # fixed effects demeaning
                         r = 0, # nubmer of factors
                         lambda = NULL, ## mc method: regularization parameter
                         nlambda = 0,
                         CV = NULL, # cross-validation
                         k = 10, # times of CV
                         cv.prop = 0.1,
                         cv.treat = TRUE,
                         cv.nobs = 3,
                         cv.donut = 1, ## cv mspe
                         criterion = "mspe",
                         binary = FALSE, # probit model
                         QR = FALSE, # QR or SVD for binary probit
                         method = "fe", # method: ife for interactive fe; mc for matrix completion
                         se = FALSE, # report uncertainties
                         vartype = "bootstrap", # bootstrap or jackknife
                         cl = NULL,
                         quantile.CI = FALSE,
                         nboots = 200, # number of bootstraps
                         alpha = 0.05, # significance level
                         parallel = TRUE, # parallel computing
                         cores = NULL, # number of cores
                         tol = 1e-3, # tolerance level
                         max.iteration = 1000,
                         seed = NULL, # set seed
                         min.T0 = NULL,
                         max.missing = NULL,
                         proportion = 0.3,
                         pre.periods = NULL,
                         f.threshold = 0.5, # equiv
                         tost.threshold = NULL,
                         knots = NULL,
                         degree = 2,  # wald = FALSE,
                         sfe = NULL,
                         cfe = NULL,
                         balance.period = NULL, # the pre and post periods for balanced samples
                         fill.missing = FALSE, # whether to balance missing observations
                         placeboTest = FALSE, # placebo test
                         placebo.period = NULL, # placebo test period
                         carryoverTest = FALSE, # carry-over test
                         carryover.period = NULL, # carry-over period
                         carryover.rm = NULL,
                         loo = FALSE, # leave one period out placebo
                         permute = FALSE, ## permutation test
                         m = 2, ## block length
                         normalize = FALSE,
                         keep.sims = FALSE
                        ) {

    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------##
    placeboEquiv <- loo
    permu.dimension <- 'time'

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
        if (! vartype %in% c("bootstrap", "jackknife", "parametric")) {
            stop("\"vartype\" option misspecified.")
        }
      if (vartype  == "parametric" && method %in% c("fe","ife", "mc", "both","polynomial","cfe")) {
            stop("The \"parametric\" option is only available for the \"gsynth\" method.")
        }
    }

    if(!is.null(W)){
        if(length(W)!=1){
            stop("\"W\" should have only one element.")
        }
        if(!W %in% colnames(data)){
            stop("\"W\" is not in the dataset.")
        }
        if(is.numeric(data[,W])==FALSE){
            stop("\"W\" should be numeric.")
        }
        if(sum(data[,W]<0)>0){
            stop("\"W\" must be strictly positive.")
        }
        if(0 %in% data[,W]){
            data <- data[which(data[,W]>0),]
        }
    }

    ## check duplicated observations
    unique_label <- unique(paste(data[,index[1]],"_",data[,index[2]],sep=""))
    if (length(unique_label)!= dim(data)[1]) {
        stop("Observations are not uniquely defined by unit and time indicators.")
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
    if (!method %in% c("fe", "ife", "mc", "both", "polynomial", "cfe","gsynth")) {
        stop("\"method\" option misspecified; choose from c(\"fe\",\"gsynth\", \"ife\", \"mc\", \"both\", \"polynomial\",\"cfe\").")
    }

    if(is.null(min.T0)){
        if(method == "fe"){
            min.T0 <- 1
        }
        if(method %in% c("ife","mc","both","gsynth")){
            min.T0 <- 5
        }
        if(method %in% c("polynomial","cfe")){
            min.T0 <- 2
        }
    }
    else{
        if(min.T0 <= 0){
            stop("\"min.T0\" option should be larger than 0.\n")
        }
    }

    ## the default setting of CV
    if(is.null(CV)){
        if (method == "fe") {
            r <- 0
            CV <- FALSE
            method <- "ife"
        }
        else if (method %in% c("polynomial","cfe")) {
            CV <- FALSE
        }
        else if(method == "both"){
            CV <- TRUE
            if(length(r)==1 & r==0){
                r <- c(0:5)
            }
        }
        else if(method %in% c("ife","gsynth")){
            if(length(r)==1){
                CV <- FALSE
            }
            else if(length(r)>1){
                CV <- TRUE
            }
        }
        else if(method == "mc"){
            if(length(lambda)==1){
                CV <- FALSE
            }
            else if(length(lambda)>1 | is.null(lambda)){
                CV <- TRUE
            }
        }
    }
    else{
        if (method == "fe") {
            r <- 0
            CV <- FALSE
            method <- "ife"
        }
        else if (method %in% c("polynomial","cfe")) {
            CV <- FALSE
        }
        else if(method == "both"){
            CV <- TRUE
        }
    }



    if (!criterion %in% c("mspe", "gmspe", "moment", "pc")) {
        stop("\"criterion\" option misspecified; choose from c(\"mspe\", \"gmspe\",\"moment\", \"pc\").")
    }

    #if (!criterion %in% c("mspe","pc")) {
    #    stop("\"criterion\" option misspecified; choose from c(\"mspe\", \"pc\").")
    #}

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
            message("No lambda is supplied. FEct is applied.")
            method <- "ife"
            r <- 0
        }
    }

    ## leave one period out placebo
    if (placeboEquiv == TRUE) {
        if(se!=TRUE){
            message("For leave one period out placebo test, automatically set \"se\" to TRUE.")
            se <- TRUE
        }
        if(placeboTest == TRUE){
            message("For leave one period out placebo test, automatically set \"placeboTest\" to FALSE.")
            placeboTest <- FALSE
        }
        if(carryoverTest == TRUE){
            message("For leave one period out placebo test, automatically set \"carryoverTest\" to FALSE.")
            carryoverTest <- FALSE
        }
    }

    if(placeboTest==TRUE & carryoverTest==TRUE){
        stop("\"placeboTest\" and \"carryoverTest\" can't be performed simultaneously.")
    }

    ## CV

    if (CV == TRUE) {

        if (placeboTest == TRUE) {
            stop("Placebo test cannot be performed while doing cross-validation.")
        }

        if (carryoverTest == TRUE) {
            stop("Carry-over test cannot be performed while doing cross-validation.")
        }

        if (method %in% c("ife", "both","gsynth")) {
            if (length(r) == 2 & r[1] > r[2]) {
                stop("\"r\" option misspecified. The first element should be smaller than the second element in r().\n")
            }
        }
        if (method %in% c("mc", "both")) {
            if (nlambda <= 0) {
                stop("\"nlambda\" option misspecified.\n")
            }
        }
    }
    else {
        if (! method %in% c("gsynth","ife", "mc", "polynomial","cfe")) {
            stop("\"method\" option misspecified; please choose from c(\"gsynth\",\"ife\", \"mc\", \"polynomial\").")
        }
    }

    if (method %in% c("polynomial","cfe")) {
        if (permute == 1) {
            message("Cannot do permutation test.\n")
            permute <- 0
        }
    }


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
        #if (r>=5) {
        #    r.end <- r
        #}
        #else {
        #    r.end <- 5
        #}
        r.end <- r
    }
    else {
        r.end <- max(r)
        r <- min(r)
    }

    ## uncertainty estimates
    if (is.logical(se) == FALSE & !se%in%c(0, 1)) {
        stop("\"se\" is not a logical flag.")
    }

    if (is.logical(quantile.CI) == FALSE & !quantile.CI%in%c(0, 1)) {
        stop("\"quantile.CI\" is not a logical flag.")
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

    if(!is.null(balance.period)){
        if(length(balance.period)!=2){
            stop(" should be of length 2.")
        }
        if(balance.period[1]>0){
            stop("The first element in \"balance.period\" should be no larger than 0.")
        }
        if(balance.period[2]<1){
            stop("The second element in \"balance.period\" should be no smaller than 1.")
        }
        if (is.logical(fill.missing) == FALSE & !fill.missing%in%c(0, 1)) {
            stop("\"fill.missing\" is not a logical flag.")
        }
        if (!is.null(group)) {
            stop("\"group\" should not be used with \"balance.period\".\n")
        }
        if(!is.null(carryover.rm)){
            stop("\"balance.period\" should not be used with \"carryover.rm\".\n")
        }
        balance.periods <- c(balance.period[1]:balance.period[2])
        # treat the units with history balance.periods as a certain group
    }

    # cohort
    if (!is.null(group)) {
        if (! group %in% names(data)) {
            stop("\"group\" misspecified.\n")
        }
    }

    if (!is.null(cl)) {
        if (! cl %in% names(data)) {
            stop("\"cl\" misspecified.\n")
        }
        if(length(cl)!=1){
            stop("Length of \"cl\" must be 1.\n")
        }
        if(cl==index[1]){cl <- NULL}
    }

    if(method == 'cfe'){
        if(is.null(sfe) & is.null(cfe)){
            message("No additional sfe and cfe, use the \"fe\" estimator by default.\n")
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
            data <- data[,unique(c(index, Y, D, X, W, group,cl))]
        }
        else {
            data <- data[,unique(c(index, Y, D, X, W, cl))] ## some variables may not be used
        }
    }
    else{
        all.var <- unique(c(index,sfe,unlist(cfe),Y,D,X,W,group,cl))
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

    ## deal with Inf and -Inf
    if(Inf %in% data[,Y] | -Inf %in% data[,Y]){
        data[which(data[,Y]==Inf),Y] <- NA
        data[which(data[,Y]==-Inf),Y] <- NA
        warning("Detect infinite values in outcome, automatically replace them with NA.\n")
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
    Wname <- W

    if (!is.null(clname)) {
        if (!clname %in% index) {
            data[, clname] <- as.numeric(as.factor(data[, clname]))
        }
    }

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

    ## check missingness
    #if (sum(is.na(data[, Yname])) > 0 & na.rm == TRUE) {
    #    stop(paste("Missing values in variable \"", Yname,"\".", sep = ""))
    #}
    if (sum(is.na(data[, Dname])) > 0) {
        stop(paste("Missing values in variable \"", Dname,"\".", sep = ""))
    }
    if (!(1%in%data[, Dname] & 0%in%data[,Dname] & length(unique(data[,Dname])) == 2)) {
        stop(paste("Error values in variable \"", Dname,"\".", sep = ""))
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
        # if X and Y have missing values while D doesn't have missing values, save the value of D
        # check if some units or periods are completely missing after drop missing values of X or Y
        data.full <- data[,c(Dname,id,time)]
        data <- na.omit(data)
        data.old <- data
        TT <- length(unique(data[,time]))
        N <- length(unique(data[,id]))
        if(TT!=TT.old){
            message("Some periods are totally removed after drop missing values of the outcome or covariates.\n")
        }
        if(N!=N.old){
            message("Some units are totally removed after drop missing values of the outcome or covariates.\n")
        }
        id.series <- unique(sort(data[,id])) ## unit id
        time.uni <- unique(sort(data[,time])) ## period

        #remove these dropped units or periods in data.full
        if(length(setdiff(data.full[,id],data[,id]))>0){
            rm.na.id <- setdiff(data.full[,id],data[,id])
            data.full <- data.full[which(!data.full[,id] %in% rm.na.id),]
        }
        else{
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

    ## gen group matrix
    if (!is.null(group)) {
        rawgroup <- data[, group]
        newgroup <- as.numeric(as.factor(rawgroup))
        data[, group] <- newgroup
        rawgroup <- cbind.data.frame(rawgroup, newgroup)
        rawgroup <- rawgroup[!duplicated(rawgroup[, 1]),]
    }

    ##message("\nOK1\n")




    ## check variation in x
    if (p > 0) {
        for (i in 1:p) {
            if (sum(is.na(data[, Xname[i]])) > 0 & na.rm == TRUE) {
                stop(paste("Missing values in variable \"", Xname[i],"\".", sep = ""))
            }
            if(force %in% c(1,3)){
                if (sum(tapply(data[, Xname[i]], data[, id], var), na.rm = TRUE) == 0) {
                    stop(paste("Variable \"",Xname[i], "\" is unit-invariant. Try to remove it.", sep = ""))
                }
            }
            if(force %in% c(2,3)){
                if (sum(tapply(data[, Xname[i]], data[, time], var), na.rm = TRUE) == 0) {
                    stop(paste("Variable \"",Xname[i], "\" is time-invariant. Try to remove it.", sep = ""))
                }
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
        data[,id] <- as.numeric(as.factor(data[,id]))

        ob.indicator <- data[,time]
        id.indicator <- table(data[, id])
        sub.start <- 1
        for (i in 1:(N - 1)) {
            sub.start <- sub.start + id.indicator[i]
            sub.end <- sub.start + id.indicator[i+1] - 1
            ob.indicator[sub.start:sub.end] <- ob.indicator[sub.start:sub.end] + i * TT
        }

        variable <- c(Yname, Dname, Xname, id, time)

        if(!is.null(group)) {
            variable <- c(variable, group)
        }

        if(!is.null(W)){
            variable <- c(variable, Wname)
        }

        if(method == 'cfe'){
            variable <- unique(c(sfe,unlist(cfe),variable))
        }

        if(!is.null(cl)){
            variable <- unique(c(variable,cl))
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

    # recode treatment status for osbervations exposed to carryover effects
    hasCarryover <- 0
    if (!is.null(carryover.rm)) {
        if (length(carryover.rm) == 1 & class(carryover.rm)[1] == "numeric") {
            if (carryover.rm > 0) {
                #print(colnames(data))
                newT <- c(1:TT)
                data <- data[order(data[, id], data[, time]),]
                tempID <- unique(data[, id])
                for (i in tempID) {
                    subpos <- which(data[, id] == i)
                    subtime <- data[subpos, time]
                    subd <- data[subpos, Dname]
                    if (sum(subd) >= 1) {
                        tr.time <- subtime[which(subd == 1)]
                        cr.time <- c() # carryover period
                        for (k in 1:carryover.rm) {
                          cr.time <- c(cr.time, tr.time + k)
                        }
                        # note: if a period has both treatment effect and carryover effect,
                        # regard carryover effect as 0
                        cr.time <- unique(cr.time)
                        cr.time <- setdiff(cr.time, tr.time)

                        cr.pos <- subpos[which(subtime %in% cr.time)]

                        if (length(cr.pos) > 0) {
                          data[cr.pos, Dname] <- 2
                        }

                    }
                }

            }
        }
    }

    if (2 %in% data[, Dname]) {
        hasCarryover <- 1
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

    if(!is.null(W)){
        W <- matrix(data[, Wname], TT, N)
    }

    if(!is.null(cl)){
        cl <- matrix(data[, clname], TT, N)
        # for each column, replace 0 with the mean of non-zero values
        cl <- apply(cl, 2, function(column) {
                                    column[column == 0] <- mean(column[column != 0])
                                    return(column)})
    }

    ## each unit should have the same group index

    ## message("\nOK2\n")

    ##treatment indicator: incorporates reversal treatments
    ## D==1 -> treatment
    D.origin <- D <- matrix(data[, Dname], TT, N)
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
        stop ("All treated units have been removed. Please specify a smaller min.T0.\n")
    }
    ## T0.min : minimum T0
    ## min.T0: manually set
    ## rm.tr.id: relative location of treated units (within all treated units)
    ## that will be removed
    if (T0.min < min.T0) {
        message(paste0("For identification purposes, units whose number of untreated periods <",min.T0," are dropped automatically.\n"))
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
        if (!is.null(cl)) {
            cl <- as.matrix(cl[,-rm.id])
        }
        if (!is.null(W)) {
            W <- as.matrix(W[,-rm.id])
        }
        if(method == "cfe"){
            for(ind.name in names(index.matrix)){
                index.matrix[[ind.name]] <- as.matrix(index.matrix[[ind.name]][,-rm.id])
            }
        }
    }

    ## message("\nOK1\n")

    ## 2. check if some periods when all units are missing or treated
    I.use <- apply(II, 1, sum)
    if (0%in%I.use) {
        for (i in 1:TT) {
            if (I.use[i] == 0) {
                message("\nThere are not any observations under control at ",time.uni[i],", drop that period.\n")
            }
        }
        if (method %in% c("polynomial")) {
            message("\nThere are not any observations at some periods. Estimation results may not be reliable. Please use time fixed effects.\n")
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

        if (!is.null(W)) {
            W <- W[-which(I.use == 0),]
        }

        if (!is.null(cl)) {
            cl <- cl[-which(I.use == 0),]
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

    ## message("\nOK2\n")
    ## 3. relative period
    T.on <- matrix(NA, TT, (N - length(rm.id)))
    D1 <- D2 <- NULL
    T.on.carry <- NULL
    if (hasCarryover == 0) {
        for (i in 1:(N - length(rm.id))) {
            T.on[, i] <-  get_term(D[, i], I.D[, i], type = "on")
        }
    }
    else {
        # separate T.on and T.on.carry
        D1 <- D2 <- D
        D1[which(D1 == 2)] <- 0
        D2[which(D2 == 1)] <- 0
        D2[which(D2 == 2)] <- 1

        T.on.carry <- matrix(NA, TT, (N - length(rm.id)))

        for (i in 1:(N - length(rm.id))) {
            T.on[, i] <-  get_term(D1[, i], I.D[, i], type = "on")
            T.on.carry[, i] <-  get_term(D2[, i], I.D[, i], type = "on")
        }
        T.on[which(D == 2)] <- NA ## remove carryover effect
        T.on.carry[which(T.on.carry <= 0)] <- NA ## only keep carryover effect
    }
    rm(D1, D2)
    calendar.time <- as.matrix(replicate((N - length(rm.id)), c(time.uni)))

    ##3.1 balance samples
    ## for balance group, add group indicator
    T.on.balance <- matrix(NA, TT, (N - length(rm.id)))
    if(!is.null(balance.period)){
        if(fill.missing == FALSE){
            T.on.miss <- T.on
            T.on.miss[which(I==0)] <- NA
            T.on.balance <- apply(T.on.miss,2,function(x) v_replace(balance.periods,x))
        }
        if(fill.missing == TRUE){
            T.on.balance <- apply(T.on,2,function(x) v_replace(balance.periods,x))
        }
        if(sum(!is.na(T.on.balance))==0){
            stop("No Balanced Sample Found.\n")
        }
    }

    ## 4. check reversals
    D1 <- D
    if (hasCarryover == 1) {
        D1[which(D == 2)] <- 1
    }
    D.fake <- apply(D1, 2, function(vec){cumsum(vec)})
    D.fake <- ifelse(D.fake > 0, 1, 0)
    D.fake[which(I.D==0)] <- 0
    Nrev <- sum(apply(D.fake == D1, 2, sum) != TT)
    hasRevs <- ifelse(Nrev > 0, 1, 0)
    if(hasRevs == FALSE & carryoverTest == TRUE){
        stop("Treatment status have no reversals. Cannot perform \"carryoverTest\" in this case.")
    }
    if(hasRevs == TRUE & method == "gsynth"){
        stop("Gsynth can't be used when treatments have reversals.")
    }

    ## 5. switch-off periods
    T.off <- NULL
    if (hasRevs == 1) {
        T.off <- matrix(NA, TT, (N - length(rm.id)))
        for (i in 1:(N - length(rm.id))) {
            T.off[, i] <-  get_term(D1[,i], I.D[,i], type = "off")
        }
    }
    rm(D1)

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
            if (hasCarryover) {
                T.on.carry <- as.matrix(T.on.carry[, -rm.id.2.pos])
            }
            if(hasRevs){
                T.off <- as.matrix(T.off[,-rm.id.2.pos])
            }
            if (!is.null(group)) {
                G <- as.matrix(G[,-rm.id.2.pos])
            }
            if (!is.null(W)) {
                W <- as.matrix(W[,-rm.id.2.pos])
            }
            if (!is.null(cl)) {
                cl <- as.matrix(cl[,-rm.id.2.pos])
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
            if (hasCarryover) {
                T.on.carry <- as.matrix(T.on.carry[, -rm.id.3.pos])
            }
            if(hasRevs){
                T.off <- as.matrix(T.off[,-rm.id.3.pos])
            }
            if (!is.null(group)) {
                G <- as.matrix(G[,-rm.id.3.pos])
            }
            if (!is.null(W)) {
                W <- as.matrix(W[,-rm.id.3.pos])
            }
            if (!is.null(cl)) {
                cl <- as.matrix(cl[,-rm.id.3.pos])
            }
            if(method == "cfe"){
                for(ind.name in names(index.matrix)){
                    index.matrix[[ind.name]] <- as.matrix(index.matrix[[ind.name]][,-rm.id.3.pos])
                }
            }
        }
    }

    ## recover treatment indicators
    D.pos <- NULL
    ##DD <- D
    if (hasCarryover) {
        D.pos <- which(D == 2)
        if (length(D.pos) > 0) {
            D[D.pos] <- 0
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
            stop("Some units do not have any untreated observations.")
        }
    }
    ## cohort
    g.level <- NULL
    if (!is.null(group)) {
        g.level <- unique(c(G))
        g.level <- g.level[!is.na(g.level)]
        rownames(rawgroup) <- rawgroup[,'newgroup']
        names(g.level) <- rawgroup[as.character(g.level),'rawgroup']
        g.level <- sort(g.level)
    }

    ##-------------------------------##
    ## Register clusters
    ##-------------------------------##

    if((se == TRUE | permute == TRUE) & parallel==FALSE){
        ## set seed
        if (is.null(seed) == FALSE) {
            set.seed(seed+1)
        }
    }

    if ((se == TRUE | permute == TRUE) & parallel==TRUE) {
        ## set seed
        if (is.null(seed) == FALSE) {
            set.seed(seed)
        }
        if (is.null(cores) == TRUE) {
            cores <- min(detectCores() - 2, 8) # default to 8 cores if not specified
        }
        para.clusters <- future::makeClusterPSOCK(cores)
        registerDoParallel(para.clusters)
        if (is.null(seed) == FALSE) {
            registerDoRNG(seed)
        }
        message("Parallel computing ...\n")
    }

    ##-------------------------------##
    ## run main program
    ##-------------------------------##
    if (se == FALSE) {

        if (CV == TRUE) {
            if (binary == FALSE) {
                out <- fect_cv(Y = Y, D = D, X = X, W = W,
                               I = I, II = II,
                               T.on = T.on, T.off = T.off, T.on.carry = T.on.carry,
                               T.on.balance = T.on.balance,
                               balance.period = balance.period,
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
                               tol = tol, max.iteration = max.iteration, norm.para = norm.para,
                               group.level = g.level, group = G)
            }
            else {
                out <- fect_binary_cv(Y = Y, D = D, X = X,
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
                out <- fect_fe(Y = Y, D = D, X = X,
                               W = W, I = I, II = II,
                               T.on = T.on, T.off = T.off, r.cv = r, T.on.carry = T.on.carry,
                               T.on.balance = T.on.balance,
                               balance.period = balance.period,
                               binary = binary, QR = QR,
                               force = force, hasRevs = hasRevs,
                               tol = tol , max.iteration = max.iteration,  boot = 0,
                               norm.para = norm.para,
                               placeboTest = placeboTest,
                               placebo.period = placebo.period,
                               carryoverTest = carryoverTest,
                               carryover.period = carryover.period,
                               group.level = g.level, group = G)
            }
            else if(method == "gsynth"){
                out <- fect_gsynth(Y = Y, D = D, X = X,
                                   W = W, I = I, II = II,
                                   T.on = T.on, T.off = T.off, r = r, CV = 0,
                                   T.on.balance = T.on.balance,
                                   balance.period = balance.period,
                                   binary = binary, QR = QR,
                                   force = force, hasRevs = hasRevs,
                                   tol = tol , max.iteration = max.iteration, boot = 0,
                                   norm.para = norm.para,
                                   placeboTest = placeboTest,
                                   placebo.period = placebo.period,
                                   carryoverTest = carryoverTest,
                                   carryover.period = carryover.period,
                                   group.level = g.level, group = G)
            }
            else if (method == "mc") {
                out <- fect_mc(Y = Y, D = D, X = X,
                               W = W, I = I, II = II,
                               T.on = T.on, T.off = T.off, T.on.carry = T.on.carry,
                               T.on.balance = T.on.balance,
                               balance.period = balance.period,
                               lambda.cv = lambda,
                               force = force, hasRevs = hasRevs,
                               tol = tol , max.iteration = max.iteration, boot = 0,
                               norm.para = norm.para,
                               placeboTest = placeboTest,
                               placebo.period = placebo.period,
                               carryoverTest = carryoverTest,
                               carryover.period = carryover.period,
                               group.level = g.level, group = G)
            }
            else if (method %in% c("polynomial",  "cfe")) {
                out <- fect_polynomial(Y = Y, D = D, X = X,
                                       W = W, I = I, II = II,
                                       T.on = T.on, T.on.carry = T.on.carry,
                                       T.on.balance = T.on.balance,
                                       balance.period = balance.period,
                                       T.off = T.off, method = method,
                                       degree = degree,
                                       sfe = sfe, cfe = cfe,
                                       ind.matrix = index.matrix,
                                       knots = knots, force = force,
                                       hasRevs = hasRevs, tol = tol , max.iteration = max.iteration, boot = 0,
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
            if(method %in% c("polynomial",  "cfe")){
                I <- out$I
                II <- out$II
            }

        }
    }
    else { # SE == TRUE

        out <- fect_boot(Y = Y, D = D, X = X,
                         W = W, I = I, II = II,
                         T.on = T.on, T.off = T.off, T.on.carry = T.on.carry, cl = cl,
                         T.on.balance = T.on.balance, balance.period = balance.period,
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
                         tol = tol , max.iteration = max.iteration, norm.para = norm.para,
                         placeboTest = placeboTest,
                         placebo.period = placebo.period,
                         carryoverTest = carryoverTest,
                         carryover.period = carryover.period,
                         vartype = vartype,
                         quantile.CI = quantile.CI,
                         nboots = nboots, parallel = parallel,
                         cores = cores, group.level = g.level, group = G,
                         keep.sims=keep.sims)

        if(method %in% c("polynomial",  "cfe")){
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
    #print(pre.periods)
    max.pre.periods <- out$time[which(out$count >= max.count * proportion & out$time <= 0)]
    all.pre.periods <- out$time[which(out$time <= 0)]
    if (is.null(pre.periods) == TRUE) {
        pre.periods <- max.pre.periods
    }
    else if(length(pre.periods)>0) {
        pre.periods <- intersect(pre.periods[1]:pre.periods[length(pre.periods)], max.pre.periods)
    }
    else{
        pre.periods <- NULL
    }

    if(!is.null(pre.periods)){
        pre.term <- pre.periods
        N_bar <- max(out$count[which(out$time %in% pre.periods)])
    }

    if (placeboEquiv == TRUE) {
        pre.term <- all.pre.periods
        r.cv <- out$r.cv
        lambda.cv <- out$lambda.cv
        method <- out$method
        if (method == "fe") {
            method <- "ife"
        }

        message("\nOut-of-Sample Test...\n")
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
            pW <- W
            if (!is.null(cl)) {
                p.cl <- cl
            }else{p.cl <- NULL}

            pII[placebo.pos] <- 0 ## placebo treatment


            ## remove treated units that have too few observations
            T0.2 <- apply(pII, 2, sum)

            if (sum(T0.2[which(apply(D, 2, sum) > 0)] >= min.T0) == 0) {
                message("\n")
                message(paste("All treated units have been removed for period ", kk, sep = ""))
                message("\n")
                jj <- jj - 1

            } else {

                te <- paste("Pre-period ", kk, sep = "")
                if (kk == 0) {
                    te <- paste(te, "(one period before treatment)", sep = " ")
                }
                message("\n")
                message(te)
                message("\n")

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
                    if (!is.null(cl)) {
                        p.cl <- cl[-rm.id.2.pos]
                    }else{p.cl <- NULL}
                    if(hasRevs){
                        pT.off <- as.matrix(T.off[,-rm.id.2.pos])
                    }
                    if (!is.null(group)) {
                        ## group <- group[-rm.id.2.pos]
                        ## rawgroup <- rawgroup[-rm.id.2.pos]
                        pG <- as.matrix(G[,-rm.id.2.pos])
                    }
                    if (!is.null(W)) {
                        ## group <- group[-rm.id.2.pos]
                        ## rawgroup <- rawgroup[-rm.id.2.pos]
                        pW <- as.matrix(W[,-rm.id.2.pos])
                    }
                }


                p.out <- fect_boot(Y = pY, D = pD, X = pX,
                             W = pW, I = pI, II = pII,
                             T.on = pT.on, T.off = pT.off, cl = p.cl,T.on.carry = T.on.carry,
                             method = method, degree = degree,
                             knots = knots, criterion = criterion,
                             CV = 0, k = k, cv.prop = cv.prop,
                             cv.treat = cv.treat, cv.nobs = cv.nobs,
                             r = r.cv, r.end = r.end,
                             nlambda = nlambda, lambda = lambda.cv,
                             alpha = alpha, binary = binary, QR = QR,
                             force = force, hasRevs = 0,
                             tol = tol , max.iteration = max.iteration, norm.para = norm.para,
                             placeboTest = 0,
                             placebo.period = NULL,
                             carryoverTest = 0,
                             carryover.period = NULL,
                             vartype = vartype,
                             nboots = nboots, parallel = parallel,
                             quantile.CI = quantile.CI,
                             cores = cores, group.level = g.level, group = pG,
                             dis = FALSE)

                p.est.att <- p.out$est.att
                p.att.bound <- p.out$att.bound
                p.pos <- which(as.numeric(rownames(p.est.att)) == kk)
                pre.est.att[jj, ] <- p.est.att[p.pos, ]
                pre.att.bound[jj, ] <- p.att.bound[p.pos, ]
                pre.att.boot[jj, ] <- p.out$att.boot.original[p.pos, ]
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
        message("Permuting under sharp null hypothesis ... ")

        out.permute <- fect_permu(Y = Y, X = X, D = D, I = I, r.cv = out$r.cv,
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
        #closeAllConnections()
    }

    ## message("\nOK4\n")



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
    # 5 placebo or carryover
    obs.missing <- matrix(0, TT, N) ## not under treatment
    obs.missing[, rem.id] <- D + as.matrix(abs(I - 1)) * 3 ## under treatment
    obs.missing[which(obs.missing==0)] <- 2
    if(placeboTest|carryoverTest){
        obs.missing[, rem.id] <- obs.missing[, rem.id] + 3*(II.origin!=II) ##placebo or carryover
    }
    if(carryoverTest & !is.null(carryover.rm)){
        obs.missing[intersect(which(D.origin==2),which(I==1))] <- 6
    }
    obs.missing[which(obs.missing==4)] <- 3 # in case if I.D!=I
    obs.missing[, rm.id] <- 4 ## removed

    colnames(obs.missing) <- unique(sort(data.old[,id]))
    rownames(obs.missing) <- tname

    obs.missing.balance <- NA
    if(!is.null(balance.period)){
        obs.missing.balance <- obs.missing
        obs.missing.sub <- obs.missing[,rem.id]
        obs.missing.sub[which(T.on.balance>0)] <- 7
        obs.missing.sub[which(T.on.balance<=0)] <- 8
        obs.missing.sub[which(I==0)] <- 3
        obs.missing.balance[,rem.id] <- obs.missing.sub
    }

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
    out$eff.calendar <- cbind(matrix(out$eff.calendar,ncol=1),out$N.calendar)
    out$eff.calendar.fit <- cbind(matrix(out$eff.calendar.fit,ncol=1),out$N.calendar)
    rownames(out$eff.calendar) <- tname
    rownames(out$eff.calendar.fit) <- tname
    colnames(out$eff.calendar) <- c("ATT-calendar","count")
    colnames(out$eff.calendar.fit) <- c("ATT-calendar Fitted","count")

    if(se==TRUE){
        rownames(out$est.eff.calendar) <- tname
        rownames(out$est.eff.calendar.fit) <- tname
    }

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
                     W = Wname,
                     T.on = T.on,
                     G = G.old,
                     balance.period = balance.period,
                     hasRevs = hasRevs,
                     T.off = T.off,
                     T.on.balance = T.on.balance,
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
                     obs.missing = obs.missing,
                     obs.missing.balance = obs.missing.balance),
                     out)


    if (1 %in% rm.id) {
        output <- c(output,list(remove.id = remove.id))
        ## message("list of removed units:",remove.id)
        ## message("\n\n")
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

    # classic equivalence test
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


    # loo equivalence test
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





