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
        stop("Unit and time variables do not uniquely identify all observations. Some may be duplicated or Incorrectly marked in the dataset.")
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

    id.match <- as.numeric(as.factor(raw.id))
    names(id.match) <- as.character(raw.id)
    time.match <- as.numeric(as.factor(raw.time))
    names(time.match) <- as.character(raw.time)

    ## generate first treated
    ## check balanced panel and fill unbalanced panel
    Dname <- D
    I <- D <- NULL
    if (dim(data)[1] != TT*N) {
        data[,index.id] <- as.numeric(as.factor(data[,index.id]))
        data[,index.time] <- as.numeric(as.factor(data[,index.time]))
        I <- matrix(0, TT, N)
        D <- matrix(0, TT, N)
        for (i in 1:dim(data)[1]) {
            D[data[i,index.time],data[i,index.id]] <- data[i,Dname]
            I[data[i,index.time],data[i,index.id]] <- 1
        }
    } 
    else {
        data[,index.id] <- as.numeric(as.factor(data[,index.id]))
        data[,index.time] <- as.numeric(as.factor(data[,index.time]))
        I <- matrix(1, TT, N)
        D <- matrix(data[,Dname], TT, N)
    }

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