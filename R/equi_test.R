##---------------------------------------##
## Cumulative treatment effect           ##
##---------------------------------------##
att.cumu <- function(x, ## a fect object
                     period, ## range, length = 2
                     weighted = TRUE, ## weighted cumulative effect 
                     type = "on" ## switch on or switch off
                    ) {
  
  if (length(period) != 2) {
    stop(" Length of \"period\" should equal 2.")
  }
  
  ## period length
  pl <- period[2] - period[1] + 1
  
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
  
  se <- 0
  
  if (type == "on") { ## switch on effect
    time <- x$time.on
    if (weighted == TRUE) { ## weighted
      att.pos <- which(time>=period[1]&time<=period[2])
      att <- x$att.on[att.pos]
      count <- x$count.on[att.pos]
      rm.pos1 <- which(is.na(att))
      rm.pos2 <- which(is.na(count))
      if (NA %in% att | NA %in% count) {
        att <- att[-c(rm.pos1, rm.pos2)]
        count <- count[-c(rm.pos1, rm.pos2)]
      }
      catt <- sum(att*count*(length(count)/sum(count)))
      if (!is.null(x$est.att.on)) {
        se <- 1
        att.pos <- which(time>=period[1]&time<=period[2])
        att.boot <- as.matrix(x$att.on.boot[att.pos,])
        count.boot <- as.matrix(x$att.on.count.boot[att.pos,])
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
        catt.ci <- quantile(catt.boot, c(0.025, 0.975), na.rm = TRUE)
        catt.p <- get.pvalue(catt.boot)
      }
      
    } else { ## simple
      att <- x$att.on
      catt <- sum(att[which(time>=period[1]&time<=period[2])], na.rm = TRUE)
      if (!is.null(x$est.att.on)) {
        se <- 1
        att.boot <- x$att.on.boot
        catt.boot <- apply(as.matrix(att.boot[which(time>=period[1]&time<=period[2]),]), 2, sum, na.rm = TRUE)
        catt.se <- sd(catt.boot, na.rm = TRUE)
        catt.ci <- quantile(catt.boot, c(0.025, 0.975), na.rm = TRUE)
        catt.p <- get.pvalue(catt.boot)
      }
    }
  } else {
    time <- x$time.off
    if (is.null(time)) {
      stop("No switch-off effect.")
    }
    if (weighted == TRUE) {
      att.pos <- which(time>=period[1]&time<=period[2])
      att <- x$att.off[att.pos]
      count <- x$count.off[att.pos]
      rm.pos1 <- which(is.na(att))
      rm.pos2 <- which(is.na(count))
      if (NA %in% att | NA %in% count) {
        att <- att[-c(rm.pos1, rm.pos2)]
        count <- count[-c(rm.pos1, rm.pos2)]
      }
      catt <- sum(att*count*(length(count)/sum(count)))
      if (!is.null(x$est.att.off)) {
        se <- 1
        att.pos <- which(time>=period[1]&time<=period[2])
        att.boot <- as.matrix(x$att.off.boot[att.pos,])
        count.boot <- as.matrix(x$att.off.count.boot[att.pos,])
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
        catt.ci <- quantile(catt.boot, c(0.025, 0.975), na.rm = TRUE)
        catt.p <- get.pvalue(catt.boot)
      }
      
    } else {
      att <- x$att.off
      catt <- sum(att[which(time>=period[1]&time<=period[2])], na.rm = TRUE)
      if (!is.null(x$est.att.off)) {
        se <- 1
        att.boot <- x$att.off.boot
        catt.boot <- apply(as.matrix(att.boot[which(time>=period[1]&time<=period[2]),]), 2, sum, na.rm = TRUE)
        catt.se <- sd(catt.boot, na.rm = TRUE)
        catt.ci <- quantile(catt.boot, c(0.025, 0.975), na.rm = TRUE)
        catt.p <- get.pvalue(catt.boot)
      }
    }
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





