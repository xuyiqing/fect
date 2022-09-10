get.cohort <- function(data,
                       D,
                       index,
                       varname = NULL, #length of 2 c("FirstTreated","Cohort")
                       entry.time = NULL){
    data.raw <- data
    if(length(D)>1){
        stop("Treatment indicator \"D\" should have length 1.\n")
    }
    if (is.data.frame(data) == FALSE || length(class(data)) > 1) {
        data <- as.data.frame(data)
        warning("Not a data frame.")
    }
    varnames <- colnames(data)
    if(!D%in%varnames){
        stop("\"D\" misspecified.\n")
    }

    if(is.null(varname)){
        varname1 <- 'FirstTreat'
        if(varname1 %in% colnames(data)){
            warnings("The variable \"FirstTreat\" will be replaced.\n")
        }
        varname2 <- 'Cohort'
        if(varname2 %in% colnames(data)){
            warnings("The variable \"Cohort\" will be replaced.\n")
        }
    }
    else{
        if(length(varname)!=2){
            stop("\"varname\" should have length 2.")
        }
        if(varname[1] %in% colnames(data)){
            stop(paste0("Variable ",varname[1]," is already in the data. Specify another group name."))
        }
        if(varname[2] %in% colnames(data)){
            stop(paste0("Variable ",varname[2]," is already in the data. Specify another group name."))
        }        
    }

    if(length(index)!=2){
        stop("\"index\" should have length 2.\n")
    }

    for(sub.index in index){
        if(!sub.index%in%index){
            stop(paste0(sub.index,"is not in the dataset.\n"))
        }
    }

    n.obs <- dim(data)[1]
    
    ## D and index should not have NA
    for(sub.var in c(D,index)){
        if(length(which(is.na(data[,sub.var])))>0){
            stop(paste0(sub.var,"contains missing values.\n"))
        }
    }

    ## check if uniquely identified
    unique_label <- unique(paste(data[,index[1]],"_",data[,index[2]],sep=""))
    if (length(unique_label)!= dim(data)[1]) {
        stop("Observations are not uniquely defined by unit and time indicators.")
    }

    if (!(class(data[, D]) %in% c("numeric", "integer"))) {
        stop("Treatment indicator should be a numeric value.")
    }

    if (!(1%in%data[, D] & 0%in%data[,D] & length(unique(data[,D])) == 2)) {
        stop(paste("Error values in variable \"", D,"\".", sep = ""))
    }


    ## index names
    index.id <- index[1]
    index.time <- index[2]
    data.raw[,index.id] <- as.character(data.raw[,index.id])

    ## raw id and time
    raw.id <- sort(unique(data[,index[1]]))
    raw.time <- sort(unique(data[,index[2]]))
    N <- length(raw.id)
    TT <- length(raw.time)


    ## entry time
    if(is.null(entry.time)){
        stagger <- 1
    }
    else{
        stagger <- 0
        if(!is.list(entry.time)){
            stop("\"entry.time\" should be a list.\n")
        }
        all.entry.time <- c()
        for(sub.time in entry.time){
            if(!is.numeric(sub.time)){
                stop("Elements in \"entry.time\" are misspecified.\n")
            }
            if(length(sub.time)!=2){
                stop("Elements in \"entry.time\" are misspecified.\n")
            }
            #if(!(sub.time[1]%in%raw.time && sub.time[2]%in%raw.time)){
            #    stop("Elements in \"entry.time\" are misspecified.\n")
            #}
            if(sub.time[1]>sub.time[2]){
                stop("Elements in \"entry.time\" are misspecified.\n")
            }
            all.entry.time <- c(all.entry.time,sub.time)
        }
        if(length(all.entry.time)>length(unique(all.entry.time))){
            stop("\"entry.time\" has overlapped periods.\n")
        }
    }

    ## sort data
    data <- data[order(data[,index.id], data[,index.time]), ]    
    data.raw <- data
    id.match <- as.numeric(as.factor(raw.id))
    names(id.match) <- as.character(raw.id)
    time.match <- as.numeric(as.factor(raw.time))
    names(time.match) <- as.character(raw.time)

    ## generate first treated
    ## check balanced panel and fill unbalanced panel
    Dname <- D
    T.on <- I <- D <- NULL
    if (dim(data)[1] != TT*N) {
        data[,index.id] <- as.numeric(as.factor(data[,index.id]))
        data[,index.time] <- as.numeric(as.factor(data[,index.time]))
        ob.indicator <- data[,index.time]
        id.indicator <- table(data[, index.id])
        sub.start <- 1
        for (i in 1:(N - 1)) { 
            sub.start <- sub.start + id.indicator[i] 
            sub.end <- sub.start + id.indicator[i+1] - 1 
            ob.indicator[sub.start:sub.end] <- ob.indicator[sub.start:sub.end] + i * TT
        }
        variable <- c(Dname)
        data_I <- matrix(0, N * TT, 1)
        data_I[ob.indicator, 1] <- 1
        data_D <- as.matrix(data[, variable])
        data_D <- data_ub_adj(data_I, data_D)
        colnames(data_D) <- variable
        I <- matrix(1, TT, N)
        D <- matrix(data_D[, Dname],TT,N)
        I[is.nan(D)] <- 0
        D[is.nan(D)] <- 0
    } 
    else {
        data[,index.id] <- as.numeric(as.factor(data[,index.id]))
        data[,index.time] <- as.numeric(as.factor(data[,index.time]))
        I <- matrix(1, TT, N)
        D <- matrix(data[,Dname], TT, N)
    }

    T.on <- matrix(NA, TT, N)
    for (i in 1:N) {
        T.on[, i] <-  get_term(D[, i], I[, i], type = "on")
    }
    
    t.on <- c(T.on)
    use.index <- (data[,index.id]-1)*TT + data[,index.time]
    t.on <- t.on[use.index]-1

    tr.pos <- which(apply(D, 2, sum) > 0)
    co.pos <- which(apply(D, 2, sum) == 0)
    tr.name <- names(id.match)[tr.pos]
    co.name <- names(id.match)[co.pos]
    
    D.cum1 <- apply(D,2,cumsum)
    D.cum2 <- apply(D.cum1,2,cumsum)
    T0.tr <- apply(D.cum2[,tr.pos],2,function(vec){which(vec==1)})
    T0.tr.origin <- as.numeric(sapply(T0.tr,function(x){names(time.match)[which(time.match==x)]}))
    T0.co.origin <- rep(NA,length(co.name))


    first.treat <- cbind.data.frame(index.id = c(tr.name,co.name),
                                    FirstTreat = c(T0.tr.origin,T0.co.origin))

    first.treat[,varname1] <- first.treat[,"FirstTreat"]
    if(varname1!='FirstTreat'){
        first.treat[,"FirstTreat"] <- NULL
    }
    data.raw <- merge(data.raw,first.treat,by.x = index.id, by.y = "index.id")
    data.raw[,varname2] <- 'Control'
    data.raw[,'Time_to_Treatment'] <- t.on
    if(stagger==1){
        all.first.treat <- sort(unique(data.raw[,varname1]))
        for(sub.first in all.first.treat){
            cohort.name <- paste0("Cohort:",sub.first)
            data.raw[which(data.raw[,varname1] == sub.first),varname2] <- cohort.name
        }
    }
    else{
        for(sub.time in entry.time){
            cohort.name <- paste0("Cohort:",sub.time[1],'-',sub.time[2])
            data.raw[which(data.raw[,varname1]>=sub.time[1] & data.raw[,varname1]<=sub.time[2]),varname2] <- cohort.name
        }
        
        cohort.name <- "Cohort:Other"
        data.raw[which(!is.na(data.raw[,varname1]) & data.raw[,varname2] == 'Control'),varname2] <- cohort.name     
    }

    return(data.raw)
}