##---------------------------------------##
## Cumulative treatment effect           ##
##---------------------------------------##
att.cumu <- function(x, ## a fect object
                     period = NULL, ## range, length = 2
                     weighted = TRUE, ## weighted cumulative effect
                     alpha = 0.05, 
                     type = "on", ## switch on or switch off
                     plot = FALSE
                    ) {

    end <- catt <- CI.lower <- CI.upper <- NULL
    if (is.null(period)) {
        period <- c(1, max(x$time))
    }

    if (length(period) != 2) {
        stop(" Length of \"period\" should equal 2.\n")
    }

    if(period[1] > period[2]){
        stop("period[1] should be smaller than period[2].\n")
    }

    ## start point
    se <- 0
    time <- x$time
    att <- x$att
    est.att <- x$est.att

    if (!is.null(est.att) & sum(abs(c(x$att.boot)),na.rm = TRUE) != 0) {
        se <- 1
    }

    if (type == "off") {
        time <- x$time.off
        att <- x$att.off
        est.att <- x$est.att.off
        if (is.null(time)) {
            stop("No switch-off effect.")
        }
    }

    p.start <- max(period[1], min(time))
    p.end <- min(period[2], max(time))

    pl <- p.end - p.start + 1 

    if(se == 0){
        cumu.att <- matrix(NA, pl, 3)
        cumu.att[1, c(1,2)] <- p.start
        cumu.att[1, 3] <- att[which(time == p.start)]
        colnames(cumu.att) <- c("start", "end", "catt")
    }
    else if (se == 1) {
        cumu.att <- matrix(NA, pl, 7)
        cumu.att[1, c(1,2)] <- p.start
        cumu.att[1, 3] <- att[which(time == p.start)]
        cumu.att[1, 4:7] <- est.att[which(time == p.start), c("S.E.", "CI.lower", "CI.upper", "p.value")]
        colnames(cumu.att) <- c("start", "end", "catt", "S.E.", "CI.lower", "CI.upper", "p.value")
    }

    for(i in 2:pl) {
        cumu.att[i, ] <- att.cumu.sub(x, c(p.start, p.start + i - 1), weighted, alpha, type)
    }

    ## plot 
    if (plot == TRUE) {
        data <- as.data.frame(cumu.att)
        data <- data[, -1] ## remove the first column
        data <- rbind(0, data)

        lp.start <- p.start - 1

        data[1, 1] <- lp.start

        p <- ggplot(data)
        p <- p + geom_hline(yintercept = 0, colour="#FFFFFF", size = 1)
        p <- p + geom_line(data = data, aes(end, catt), size = 1.2)
        if (se == 1) {
            p <- p + geom_ribbon(data = data, aes(x = end, ymin = CI.lower, ymax = CI.upper), alpha = 0.2) 
        }

        p <- p + xlab("Relative Period") +  ylab("Estimated Cumulative Effects")

        suppressWarnings(print(p))
    }

    return(cumu.att)
}





att.cumu.sub <- function(x, ## a fect object
                         period, ## range, length = 2
                         weighted = TRUE, ## weighted cumulative effect
                         alpha = 0.05, 
                         type = "on" ## switch on or switch off
                        ) {
  
    if (length(period) != 2) {
        stop(" Length of \"period\" should equal 2.")
    }
    
    ## period length
    pl <- period[2] - period[1] + 1  
    se <- 0

    if (type == "on") {
        time <- x$time
        att <- x$att
        count <- x$count

        if (!is.null(x$est.att)) {
            se <- 1
            att.boot <- x$att.boot
            count.boot <- x$att.count.boot
        }
    } else {
        time <- x$time.off
        att <- x$att.off
        count <- x$count.off

        if (is.null(time)) {
            stop("No switch-off effect.")
        }

        if (!is.null(x$est.att.off)) {
            se <- 1
            att.boot <- x$att.off.boot
            count.boot <- x$att.off.count.boot
        }
    }


    att.pos <- which(time>=period[1]&time<=period[2])
    att <- att[att.pos]
    count <- count[att.pos]
    rm.pos1 <- which(is.na(att))
    rm.pos2 <- which(is.na(count))
    if (NA %in% att | NA %in% count) {
        att <- att[-c(rm.pos1, rm.pos2)]
        count <- count[-c(rm.pos1, rm.pos2)]
    }
    catt <- sum(att*count*(length(count)/sum(count)))


    if (se == 1) {
        att.boot <- as.matrix(att.boot[att.pos,])
        count.boot <- as.matrix(count.boot[att.pos,])
        
        nboots <- dim(att.boot)[2]
        catt.boot <- rep(NA, nboots)
        
        for (i in 1:nboots) {
            att.sub <- att.boot[,i]
            count.sub <- count.boot[,i]
            rm.pos1 <- which(is.na(att.sub))
            rm.pos2 <- which(is.na(count.sub))
            if (NA %in% att.sub | NA %in% count.sub) {
                att.sub <- att.sub[-c(rm.pos1, rm.pos2)]
                count.sub <- count.sub[-c(rm.pos1, rm.pos2)]
            }
            catt.boot[i] <- sum(att.sub*count.sub*(length(count.sub)/sum(count.sub)))*(pl)/(length(count.sub))
        }

        catt.se <- sd(catt.boot, na.rm = TRUE)
        catt.ci <- quantile(catt.boot, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
        catt.p <- get.pvalue(catt.boot)
    }
  
    if (se == 0) {
        result <- t(matrix(c(period, catt)))
        colnames(result) <- c("start", "end", "cumulative att")
    } else { ## with uncertainties
        result <- t(matrix(c(period, catt, catt.se, catt.ci, catt.p)))
        colnames(result) <- c("start", "end", "cumulative att", "S.E.", "CI.lower", "CI.upper", "p.value")
    }
  
    return(result)
} 

##############################
##   equivalence test       ##
##############################

## function to get two-sided p-values
get.pvalue <- function(vec) {
  if (NaN%in%vec|NA%in%vec) {
    nan.pos <- is.nan(vec)
    na.pos <- is.na(vec)
    pos <- c(which(nan.pos),which(na.pos))
    vec.a <- vec[-pos]
    a <- sum(vec.a >= 0)/(length(vec)-sum(nan.pos|na.pos)) * 2
    b <- sum(vec.a <= 0)/(length(vec)-sum(nan.pos|na.pos)) * 2  
  } else {
    a <- sum(vec >= 0)/length(vec) * 2
    b <- sum(vec <= 0)/length(vec) * 2  
  }
  return(min(as.numeric(min(a, b)),1))
}



