get.cohort <- function(data,
                       D,
                       index,
                       varname = NULL, # length of 4 c("FirstTreated","Cohort",'Time_to_Treatment',"Time_to_Exit")
                       start0 = FALSE, # default 0 as the last pre-treatment period, 1 as the first post-treatment period
                       entry.time = NULL,
                       drop.always.treat = FALSE) {
    data.raw <- data
    if (length(D) > 1) {
        stop("Treatment indicator \"D\" should have length 1.")
    }
    if (is.data.frame(data) == FALSE || length(class(data)) > 1) {
        data <- as.data.frame(data)
    }
    varnames <- colnames(data)
    if (!D %in% varnames) {
        stop("\"D\" misspecified.\n")
    }

    if (is.logical(start0) == FALSE & !start0 %in% c(0, 1)) {
        stop("\"start0\" is not a logical flag.")
    }

    if (is.logical(drop.always.treat) == FALSE & !drop.always.treat %in% c(0, 1)) {
        stop("\"drop.always.treat\" is not a logical flag.")
    }


    if (is.null(varname)) {
        varname1 <- "FirstTreat"
        if (varname1 %in% colnames(data)) {
            warnings("The variable \"FirstTreat\" will be replaced.\n")
            data[, varname1] <- NULL
        }
        varname2 <- "Cohort"
        if (varname2 %in% colnames(data)) {
            warnings("The variable \"Cohort\" will be replaced.\n")
            data[, varname2] <- NULL
        }
        varname3 <- "Time_to_Treatment"
        if (varname3 %in% colnames(data)) {
            warnings("The variable \"Time_to_Treatment\" will be replaced.\n")
            data[, varname3] <- NULL
        }
        varname4 <- "Time_to_Exit"
        if (varname4 %in% colnames(data)) {
            warnings("The variable \"Time_to_Exit\" will be replaced.\n")
            data[, varname4] <- NULL
        }
    } else {
        if (length(varname) != 4) {
            stop("\"varname\" should have three elements.")
        }
        if (varname[1] %in% colnames(data)) {
            warnings(paste0("Variable ", varname[1], " will be replaced.\n"))
            data[, varname1] <- NULL
        }
        if (varname[2] %in% colnames(data)) {
            warnings(paste0("Variable ", varname[2], " will be replaced.\n"))
            data[, varname2] <- NULL
        }
        if (varname[3] %in% colnames(data)) {
            warnings(paste0("Variable ", varname[3], " will be replaced.\n"))
            data[, varname3] <- NULL
        }
        if (varname[4] %in% colnames(data)) {
            warnings(paste0("Variable ", varname[4], " will be replaced.\n"))
            data[, varname4] <- NULL
        }
    }

    if (length(index) != 2) {
        stop("\"index\" should have length 2.\n")
    }

    for (sub.index in index) {
        if (!sub.index %in% index) {
            stop(paste0(sub.index, "is not in the dataset.\n"))
        }
    }

    n.obs <- dim(data)[1]

    ## D and index should not have NA
    for (sub.var in c(D, index)) {
        if (length(which(is.na(data[, sub.var]))) > 0) {
            stop(paste0(sub.var, "contains missing values.\n"))
        }
    }

    ## check if uniquely identified
    unique_label <- unique(paste(data[, index[1]], "_", data[, index[2]], sep = ""))
    if (length(unique_label) != dim(data)[1]) {
        stop("Observations are not uniquely defined by unit and time indicators.")
    }

    if (!(class(data[, D]) %in% c("numeric", "integer"))) {
        stop("Treatment indicator should be a numeric value.")
    }

    if (!(1 %in% data[, D] & 0 %in% data[, D] & length(unique(data[, D])) == 2)) {
        stop(paste("Error values in variable \"", D, "\".", sep = ""))
    }


    ## index names
    index.id <- index[1]
    index.time <- index[2]
    data.raw[, index.id] <- as.character(data.raw[, index.id])

    ## raw id and time
    raw.id <- sort(unique(data[, index[1]]))
    raw.time <- sort(unique(data[, index[2]]))
    N <- length(raw.id)
    TT <- length(raw.time)

    ## entry time
    if (is.null(entry.time)) {
        staggered <- 1
    } else {
        staggered <- 0
        if (!is.list(entry.time)) {
            stop("\"entry.time\" should be a list.\n")
        }
        all.entry.time <- c()
        for (sub.time in entry.time) {
            if (!is.numeric(sub.time)) {
                stop("Elements in \"entry.time\" are misspecified.\n")
            }
            if (length(sub.time) != 2) {
                stop("Elements in \"entry.time\" are misspecified. Each element should have length 2.\n")
            }
            if (sub.time[1] > sub.time[2]) {
                stop("Elements in \"entry.time\" are misspecified. For each element, the second should not be smaller than the first.\n")
            }
            all.entry.time <- c(all.entry.time, sub.time)
        }
        if (length(all.entry.time) > length(unique(all.entry.time))) {
            stop("\"entry.time\" has overlapped periods.\n")
        }
    }

    ## sort data
    data <- data[order(data[, index.id], data[, index.time]), ]
    data.raw <- data
    id.match <- as.numeric(as.factor(raw.id))
    names(id.match) <- as.character(raw.id)
    time.match <- as.numeric(as.factor(raw.time))
    names(time.match) <- as.character(raw.time)

    ## generate first treated
    ## check balanced panel and fill unbalanced panel
    Dname <- D
    T.on <- I <- D <- NULL
    if (dim(data)[1] != TT * N) {
        data[, index.id] <- as.numeric(as.factor(data[, index.id]))
        data[, index.time] <- as.numeric(as.factor(data[, index.time]))
        ob.indicator <- data[, index.time]
        id.indicator <- table(data[, index.id])
        sub.start <- 1
        for (i in 1:(N - 1)) {
            sub.start <- sub.start + id.indicator[i]
            sub.end <- sub.start + id.indicator[i + 1] - 1
            ob.indicator[sub.start:sub.end] <- ob.indicator[sub.start:sub.end] + i * TT
        }
        variable <- c(Dname)
        data_I <- matrix(0, N * TT, 1)
        data_I[ob.indicator, 1] <- 1

        data_D <- matrix(NaN, N * TT, 1)
        data_D[ob.indicator, 1] <- as.matrix(data[, variable])
        colnames(data_D) <- variable
        I <- matrix(1, TT, N)
        D <- matrix(data_D[, Dname], TT, N)
        I[is.nan(D)] <- 0
        D[is.nan(D)] <- 0
    } else {
        data[, index.id] <- as.numeric(as.factor(data[, index.id]))
        data[, index.time] <- as.numeric(as.factor(data[, index.time]))
        I <- matrix(1, TT, N)
        D <- matrix(data[, Dname], TT, N)
    }

    T.on <- matrix(NA, TT, N)
    for (i in 1:N) {
        T.on[, i] <- get_term(D[, i], I[, i], type = "on")
    }
    T.off <- matrix(NA, TT, N)
    for (i in 1:N) {
        T.off[, i] <- get_term(D[, i], I[, i], type = "off")
    }

    t.on <- c(T.on)
    t.off <- c(T.off)
    use.index <- (data[, index.id] - 1) * TT + data[, index.time]
    if (start0 == TRUE) {
        t.on <- t.on[use.index] - 1
        t.off <- t.off[use.index] - 1
    } else {
        t.on <- t.on[use.index]
        t.off <- t.off[use.index]
    }

    tr.pos <- which(apply(D, 2, sum) > 0)
    co.pos <- which(apply(D, 2, sum) == 0)
    tr.name <- names(id.match)[tr.pos]
    co.name <- names(id.match)[co.pos]

    D.cum1 <- apply(D, 2, cumsum)
    D.cum2 <- apply(D.cum1, 2, cumsum)
    T0.tr <- apply(matrix(D.cum2[, tr.pos], nrow = dim(D.cum2)[1]), 2, function(vec) {
        which(vec == 1)
    })
    T0.tr.origin <- as.numeric(sapply(T0.tr, function(x) {
        names(time.match)[which(time.match == x)]
    }))
    T0.co.origin <- rep(NA, length(co.name))


    first.treat <- cbind.data.frame(
        index.id = c(tr.name, co.name),
        FirstTreat = c(T0.tr.origin, T0.co.origin)
    )

    first.treat[, varname1] <- first.treat[, "FirstTreat"]
    if (varname1 != "FirstTreat") {
        first.treat[, "FirstTreat"] <- NULL
    }
    data.raw <- merge(data.raw, first.treat, by.x = index.id, by.y = "index.id")
    data.raw[, varname2] <- "Control"
    data.raw[, varname3] <- t.on
    if (sum(!is.na(t.off)) > 0) {
        data.raw[, varname4] <- t.off
    }

    if (staggered == 1) {
        all.first.treat <- sort(unique(data.raw[, varname1]))
        for (sub.first in all.first.treat) {
            cohort.name <- paste0("Cohort:", sub.first)
            data.raw[which(data.raw[, varname1] == sub.first), varname2] <- cohort.name
        }
    } else {
        for (sub.time in entry.time) {
            cohort.name <- paste0("Cohort:", sub.time[1], "-", sub.time[2])
            data.raw[intersect(which(data.raw[, varname1] >= sub.time[1]), which(data.raw[, varname1] <= sub.time[2])), varname2] <- cohort.name
        }
        cohort.name <- "Cohort:Other"
        data.raw[which(!is.na(data.raw[, varname1]) & data.raw[, varname2] == "Control"), varname2] <- cohort.name
    }


    if (drop.always.treat) {
        size.0 <- dim(data.raw)[1]
        unit.index <- index[1]
        data.raw <- as.data.frame(data.raw %>% group_by(get(unit.index)) %>% mutate(treatment_mean = mean(get(Dname), na.rm = TRUE)))
        data.raw <- data.raw[which(data.raw$treatment_mean < 1), ]
        size.1 <- dim(data.raw)[1]
        diff.1 <- size.0 - size.1
        message(paste0("Number of Always Treated Observations Removed:", diff.1, ".\n"))
        data.raw[, "treatment_mean"] <- NULL
    }

    return(data.raw)
}


get_term <- function(d,
                     ii,
                     type = "on") {
    dd <- d
    iii <- ii
    first.pos <- min(which(iii == 1))
    if (first.pos != 1) {
        dd <- dd[-(1:(first.pos - 1))]
        iii <- iii[-(1:(first.pos - 1))]
    }
    T <- length(dd)
    if (0 %in% iii) {
        if (T > 1) {
            for (i in 1:(T - 1)) {
                if (iii[i + 1] == 0) {
                    dd[i + 1] <- dd[i]
                }
            }
        }
    }

    if (type == "off") {
        dd <- abs(dd - 1)
    }
    d1 <- dd[1:(T - 1)]
    d2 <- dd[2:T]

    if (T == 1) {
        term <- rep(NA, 1)
    } else if (sum(d1 == d2) == (T - 1)) {
        term <- rep(NA, T)
    } else {
        change.pos <- which(d1 != d2) + 1
        change.length <- length(change.pos)
        term <- NULL
        if (dd[1] == 0) {
            for (i in 1:(change.length)) {
                if (i == 1) {
                    part.term <- (2 - change.pos[i]):0
                } else {
                    if (i %% 2 == 0) {
                        part.term <- 1:(change.pos[i] - change.pos[i - 1])
                    } else {
                        part.term <- (change.pos[i - 1] - change.pos[i] + 1):0
                    }
                }
                term <- c(term, part.term)
            }
        } else if (dd[1] == 1) {
            for (i in 1:(change.length)) {
                if (i == 1) {
                    part.term <- rep(NA, change.pos[i] - 1)
                } else {
                    if (i %% 2 == 0) {
                        part.term <- (change.pos[i - 1] - change.pos[i] + 1):0
                    } else {
                        part.term <- 1:(change.pos[i] - change.pos[i - 1])
                    }
                }
                term <- c(term, part.term)
            }
        }
        if (dd[change.pos[change.length]] == 0) {
            term <- c(term, rep(NA, (T - change.pos[change.length] + 1)))
        } else {
            term <- c(term, 1:(T - change.pos[change.length] + 1))
        }
    }
    if (first.pos != 1) {
        term <- c(rep(NA, (first.pos - 1)), term)
    }
    return(term)
}
