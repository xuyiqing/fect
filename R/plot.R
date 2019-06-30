##########
## Plot
##########
# x: a fect object
# type of the plot; axes limits; axes labels; 
# main: whether to show the title;
# id: plot a part of units
plot.fect <- function(x,  
                      type = "raw",
                      gap.type = "on",
                      bound = NULL,
                      vis = "connected",
                      count = TRUE,
                      proportion = 0.3,
                      equiv.bound = NULL,
                      equiv.bound.size = 1,
                      effect.bound.ratio = FALSE,
                      p.value = TRUE,
                      main = NULL,
                      xlim = NULL, 
                      ylim = NULL,
                      xlab = NULL, 
                      ylab = NULL,
                      legendOff = FALSE,
                      legend.pos = NULL,
                      legend.nrow = NULL,
                      legend.labs = NULL,
                      text.pos = NULL,
                      theme.bw = FALSE,
                      gridOff = FALSE,
                      id = NULL,
                      cex.main = NULL,
                      cex.main.sub = NULL,
                      cex.axis = NULL,
                      cex.lab = NULL, 
                      cex.legend = NULL,
                      cex.text = NULL,
                      axis.adjust = FALSE,
                      axis.lab = "both",
                      axis.lab.gap = c(0, 0),
                      ...){
    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------## 
    ATT.ON2 <- ATT.ON3 <- ATT.ON4 <- NULL
    CI.lower3 <-  CI.lower4 <- CI.upper3 <- CI.upper4 <- NULL

    p <- NULL
    outcome <- NULL ## global variable
    labels1 <- labels2 <- labels3 <- NULL
    ATT.OFF <- ATT.ON <- CI.lower <- CI.upper <- NULL
    xmax <- xmin <- ymax <- ymin <- NULL
    placeboTest <- x$placeboTest
    placebo.period <- x$placebo.period
    binary <- x$binary
    wald <- !is.null(x$wald$p)


    if (class(x) != "fect") {
        stop("Not a \"fect\" object.")
    }
    if (!type %in% c("missing", "raw", "gap","equiv")) {
        stop("\"type\" option misspecified.")
    }
    if (type == "gap" && gap.type == "off" && is.null(x$att.off)) {
        stop("No switch-off treatment effect to be plotted.")
    }
    if (is.null(bound) == FALSE) {
        if (!bound %in% c("none", "min", "equiv", "both")) {
            stop("\"bound\" option misspecified.")
        } 
    }    
    
    if (is.null(xlim)==FALSE) {
        if (is.numeric(xlim)==FALSE) {
            stop("Some element in \"xlim\" is not numeric.")
        } else {
            if (length(xlim)!=2) {
                stop("xlim must be of length 2.")
            }
        }
    }
    if (is.null(ylim)==FALSE) {
        if (is.numeric(ylim)==FALSE) {
            stop("Some element in \"ylim\" is not numeric.")
        } else {
            if (length(ylim)!=2) {
                stop("ylim must be of length 2.")
            }
        }
    }

    ## gap or equiv
    if (type == "gap") {
        if (is.null(bound) == TRUE) {
            bound <- "none"
        }
        maintext <- "Estimated Average Treatment Effect"
    } else if (type == "equiv") {
        if (is.null(bound) == TRUE) {
            bound <- "both"
        }
        if (length(xlim)==0) {
            xlim <- c(-1e5, 0.5)
        } else {
            if (xlim[2]>0) {
                xlim[2]<-0.5
            }
        }
        maintext <- "Equivalence Test" 
        type <- "gap"
    }       


    if (is.null(xlab)==FALSE) {
        if (is.character(xlab) == FALSE) {
            stop("\"xlab\" is not a string.")
        } else {
            xlab <- xlab[1]
        }   
    }
    if (is.null(ylab)==FALSE) {
        if (is.character(ylab) == FALSE) {
            stop("\"ylab\" is not a string.")
        } else {
            ylab <- ylab[1]
        }   
    }
    if (is.logical(legendOff) == FALSE & is.numeric(legendOff)==FALSE) {
        stop("\"legendOff\" is not a logical flag.")
    }
    if (is.logical(gridOff) == FALSE & is.numeric(gridOff)==FALSE) {
        stop("\"gridOff\" is not a logical flag.")
    }

    if (is.null(main)==FALSE) {
        if (is.character(main) == FALSE) {
            stop("\"main\" is not a string.")
        } else {
            main <- main[1]
        }   
    }
    if (axis.adjust == TRUE) {
        angle <- 45
        x.v <- 1
        x.h <- 1
    } else {
        angle <- 0
        x.v <- 0
        if (type == "missing") {
            x.h <- 0.5
        } else {
            x.h <- 0
        }
    }


    #### font size
    ## title
    if (is.null(cex.main)==FALSE) {
        if (is.numeric(cex.main)==FALSE) {
            stop("\"cex.main\" is not numeric.")
        }
        cex.main <- 18 * cex.main
    } else {
        cex.main <- 18
    }
    ## subtitle
    if (is.null(cex.main.sub)==FALSE) {
        if (is.numeric(cex.main.sub)==FALSE) {
            stop("\"cex.main.sub\" is not numeric.")
        }
        cex.main.sub <- 16 * cex.main.sub
    } else {
        cex.main.sub <- 16
    }
    ## axis label
    if (is.null(cex.lab)==FALSE) {
        if (is.numeric(cex.lab)==FALSE) {
            stop("\"cex.lab\" is not numeric.")
        }
        cex.lab <- 15 * cex.lab
    } else {
        cex.lab <- 15
    }
    ## axis number
    if (is.null(cex.axis)==FALSE) {
        if (is.numeric(cex.axis)==FALSE) {
            stop("\"cex.axis\" is not numeric.")
        }
        cex.axis <- 15 * cex.axis
    }  else {
        cex.axis <- 15
    }
    ## legend
    if (is.null(cex.legend)==FALSE) {
        if (is.numeric(cex.legend)==FALSE) {
            stop("\"cex.legend\" is not numeric.")
        }
        cex.legend <- 15 * cex.legend
    }  else {
        cex.legend <- 15
    }
    ## text
    if (is.null(cex.text)==FALSE) {
        if (is.numeric(cex.text)==FALSE) {
            stop("\"cex.text\" is not numeric.")
        }
        cex.text <- 6 * cex.text
    }  else {
        cex.text <- 6
    }

    ## text label position
    if (!is.null(text.pos)) {
       if (length(text.pos) != 2) {
           stop(" \"text.pos\" must be of length 2. ")
       }
    }

    ##-------------------------------##
    ## Plotting
    ##-------------------------------## 

    Y <- x$Y.dat
    D <- x$D.dat
    I <- x$I.dat
    Yname <- x$Y

    index <- x$index
    unit.type <- x$unit.type
    obs.missing <- x$obs.missing
    tname <- time <- x$time

    TT <- dim(Y)[1]
    N <- dim(Y)[2]
    
    if (type == "missing") {
        if (is.null(id) == TRUE) {
            ## if (is.null(show.id) == TRUE) {
                id <- colnames(obs.missing)
            ## } else {
            ##     id <- colnames(obs.missing)[show.id]
            ## }
        }
        m.l <- length(id)
        for (i in 1:m.l) {
            if (!id[i]%in%colnames(obs.missing)) {
                stop("Some specified units are not in the data.")
            }
        }
    } else { ## raw plot
        if (is.null(id) == TRUE) {
            ## if (is.null(show.id) == TRUE) {
                id <- x$id
            ## } else {
            ##     id <- colnames(obs.missing)[show.id]
            ## }
        }
        m.l <- length(id)
        id.pos <- rep(NA, m.l)
        for (i in 1:m.l) {
            if (!id[i] %in% x$id) {
                stop("Some specified units are not in the data.")
            } else {
                id.pos[i] <- which(x$id == id[i])
            }
        }
        Y <- as.matrix(Y[, id.pos])
        I <- as.matrix(I[, id.pos])
        D <- as.matrix(D[, id.pos])
        unit.type <- unit.type[id.pos]
    }
  
    ## type of plots
    if (type == "gap") {
        if (gap.type == "on") {
            time <- x$time.on
        } else {
            time <- x$time.off
        }
    } else {
        if (!is.numeric(time[1])) {
            time <- 1:TT
        }
    }

    ## periods to show 
    time.end <- length(time)
    show.time <- 1:time.end    
    if (length(xlim) != 0) {
        show.time <- which(time >= xlim[1] & time <= xlim[2])
    }

    ## count 
    show.count <- 1:time.end
    count.on <- x$count.on
    max.count <- max(count.on[which(time >= 0)])

    if (!is.null(proportion)) {
        ## max.count <- max(count.on[which(time > 0)])
        show.count <- which(count.on >= max.count * proportion)
    }

    show <- intersect(show.count, show.time) 

    if (length(show) <= 2) {
        stop("Cannot plot.\n")
    }    

    nT <- length(show)
    time.label <- tname[show]
    T.b <- 1:length(show)
 
    ## legend on/off
    if (legendOff == TRUE) {
        legend.pos <- "none"
    } else {
        if (is.null(legend.pos)) {
            legend.pos <- "bottom"
        }
    }

    ## adjust x-axis
    scaleFUN <- function(x) sprintf("%.f", x) ## integer value at x axis

    ############  START  ############### 
    if (type == "gap") {

        ## axes labels
        if (is.null(xlab) == TRUE) {
            xlab <- paste("Time relative to Treatment")
        } else if (xlab == "") {
            xlab <- NULL
        }
    
        if (is.null(ylab) == TRUE) {
            ylab <- "Coefficient"
        } else if (ylab == "") {
            ylab <- NULL
        }

        ## 0 lines
        lcolor <- "white"
        lwidth <- 2
        if (theme.bw == TRUE) {
            lcolor <- "#AAAAAA70"
            lwidth <- 1.5
        }           
        

        if (gap.type == "on") {

            data2 <- NULL
            if (bound != "none") {
                if (is.null(x$est.att.on)) {
                    cat("No uncertainty estimates.\n")
                    bound <- "none"
                }
                if (placeboTest == TRUE) {
                    bound <- "none"
                }
            } 

            if (bound != "none") {
                
                time0 <- NULL
                if (sum(time[show] <= 0) == 0) {
                    cat("No pre-treatment periods are to be plotted.\n")
                    time0 <- 1:length(time[show])
                } else {
                    time0 <- which(time[show] <= 0)
                }
                
                att.on.sub <- as.matrix(x$est.att.on[show, ])
                
                minBound <- max(abs(att.on.sub[time0, c("CI.lower", "CI.upper")]))
                minbound <- c(-minBound, minBound)

                if (is.null(equiv.bound)) {

                    if (is.null(equiv.bound.size)) {
                        equiv.bound.size <- 1
                    } 
                    equiv.bound <- c(-abs(equiv.bound.size * x$att.avg), abs(equiv.bound.size * x$att.avg))
                
                } else {
                
                    if (length(equiv.bound) == 1) {
                        equiv.bound <- c(-abs(equiv.bound), abs(equiv.bound))
                    }

                }

                bound.time <- time[show]
                bound.time <- bound.time[which(bound.time <= 0)]

                ## add legend for 95\% CI
                set.limits <- "ci"
                if (is.null(legend.labs)==TRUE) {
                    set.labels <- "ATT (w/ 95% CI)"                    
                } else {
                    set.labels <- legend.labs
                }
                set.colors <- "#000000FF"
                set.linetypes <- "solid"
                set.size <- 1

                ## create a dataset for bound
                if (bound == "equiv") {
                    data2 <- cbind.data.frame(c(rep(equiv.bound, each = length(bound.time))))
                    names(data2) <- "bound"
                    data2$time <- rep(bound.time, 2)
                    data2$type <- rep(c("equiv"), 2 * length(bound.time))
                    data2$id <- rep(1:2, each = length(bound.time))

                    set.limits <- c(set.limits, "equiv")
                    if (is.null(legend.labs)==TRUE) {
                        set.labels <- c(set.labels, "equiv. Bound")
                    } else {
                        set.labels <- legend.labs
                    }
                    set.colors <- c(set.colors, "red")
                    set.linetypes <- c(set.linetypes, "dashed")
                    set.size <- c(set.size, 0.7)
                } 
                else if (bound == "min") {
                    data2 <- cbind.data.frame(c(rep(minbound, each = length(bound.time))))
                    names(data2) <- "bound"
                    data2$time <- rep(bound.time, 2)
                    data2$type <- rep(c("min"), 2 * length(bound.time))
                    data2$id <- rep(1:2, each = length(bound.time))
                    set.limits <- c(set.limits, "min")
                    if (is.null(legend.labs)==TRUE) {
                        set.labels <- c(set.labels, "Min. Bound")
                    } else {
                        set.labels <- legend.labs
                    }
                    set.colors <- c(set.colors, "gray50")
                    set.linetypes <- c(set.linetypes, "dashed")
                    set.size <- c(set.size, 0.7)
                }
                else if (bound == "both") {
                    data2 <- cbind.data.frame(c(rep(minbound, each = length(bound.time)), rep(equiv.bound, each = length(bound.time))))
                    names(data2) <- "bound"
                    data2$time <- rep(bound.time, 4)
                    data2$type <- rep(c("min", "equiv"), each = 2 * length(bound.time))
                    data2$id <- rep(1:4, each = length(bound.time))
                    set.limits <- c(set.limits, "min", "equiv")
                    if (is.null(legend.labs)==TRUE) {
                        set.labels <- c(set.labels, "Min. Bound", "Equiv. Bound")
                    } else {
                        set.labels <- legend.labs
                    }
                    set.colors <- c(set.colors, "gray50", "red")
                    set.linetypes <- c(set.linetypes, "dashed", "dashed")
                    set.size <- c(set.size, 0.7, 0.7)
                }
            }

            ## switch-on effect
            if (is.null(x$est.att.on)==TRUE) { 
                
                cat("Uncertainty estimates not available.\n")
                data <- cbind.data.frame(time, ATT.ON = x$att.on, count = x$count.on)[show,]
                if (length(ylim) != 0) {
                    rect.length <- (ylim[2] - ylim[1]) / 5
                    rect.min <- ylim[1]
                } else {
                    rect.length <- (max(data[,"ATT.ON"]) - min(data[,"ATT.ON"]))/2
                    rect.min <- min(data[,"ATT.ON"]) - rect.length
                } 

            } else {
                
                tb <- x$est.att.on
                data <- cbind.data.frame(time, tb, count = x$count.on)[show,]
                ## rect.length <- (max(data[,"CI.upper"]) - min(data[,"CI.lower"]))/2
                if (length(ylim) != 0) {
                    rect.length <- (ylim[2] - ylim[1]) / 5
                    rect.min <- ylim[1]
                } else {
                    rect.length <- (max(data[,"CI.upper"]) - min(data[,"CI.lower"]))/2
                    rect.min <- min(data[,"CI.lower"]) - rect.length
                }  
                if (placeboTest == TRUE) {
                    data[,"ATT.ON4"] <- data[,"ATT.ON3"] <- data[,"ATT.ON2"] <- data[,"ATT.ON"]
                    data[,"CI.lower4"] <- data[,"CI.lower3"] <- data[,"CI.lower2"] <- data[,"CI.lower"]
                    data[,"CI.upper4"] <- data[,"CI.upper3"] <- data[,"CI.upper2"] <- data[,"CI.upper"]
                    
                    pos1 <- intersect(which(data[,"time"] >= (placebo.period[1] + 1)), which(data[,"time"] <= (placebo.period[2] - 1)))
                    pos2 <- c(which(data[,"time"] <= (placebo.period[1] - 1)), which(data[,"time"] >= (placebo.period[2] + 1)))

                    pos3 <- intersect(which(data[,"time"] >= (placebo.period[1])), which(data[,"time"] <= (placebo.period[2] - 1)))
                    pos4 <- c(which(data[,"time"] <= (placebo.period[1] - 2)), which(data[,"time"] >= (placebo.period[2] + 1)))
                    
                    data[pos1, c("ATT.ON","CI.lower","CI.upper")] <- NA
                    data[pos2, c("ATT.ON2","CI.lower2","CI.upper2")] <- NA
                    data[pos3, c("ATT.ON3","CI.lower3","CI.upper3")] <- NA
                    data[pos4, c("ATT.ON4","CI.lower4","CI.upper4")] <- NA
                }

            }
             
            ## plotting
            if (bound != "none") {
                p <- ggplot(data2) +
                    xlab(xlab) +  ylab(ylab) + 
                    scale_x_continuous(labels=scaleFUN) 

                p <- p + geom_line(aes(time, bound,
                                       colour = type,
                                       linetype = type,
                                       size = type,
                                       group = id)) 
                ## legend
                if (is.null(legend.nrow) == TRUE) {
                    legend.nrow <- ifelse(length(set.limits) <= 3, 1, 2)    
                } 
                
                p <- p + scale_colour_manual(limits = set.limits,
                                             labels = set.labels,
                                             values =set.colors) +
                         scale_size_manual(limits = set.limits,
                                           labels = set.labels,
                                           values = set.size) +
                         scale_linetype_manual(limits = set.limits,
                                               labels = set.labels,
                                               values = set.linetypes) +
                         guides(linetype = guide_legend(title=NULL, nrow=legend.nrow),
                                colour = guide_legend(title=NULL, nrow=legend.nrow),
                                size = guide_legend(title=NULL, nrow=legend.nrow)) 
                
                if (effect.bound.ratio == TRUE) {
                    if (is.null(text.pos)) {
                        text.pos[1] <- min(data[,"time"], na.rm = 1)
                        text.pos[2] <- ifelse(is.null(ylim), max(data[,"CI.upper"], na.rm = 1), ylim[1])
                    }
                    p.label <- paste("ATT / Min. Bound = ", sprintf("%.3f",x$att.avg / minBound), sep="")
                    p <- p + annotate("text", x = text.pos[1], y = text.pos[2], 
                        label = p.label, size = cex.text, hjust = 0)
                
                }
                
            } else {
                p <- ggplot(data) +
                    xlab(xlab) +  ylab(ylab) + 
                    scale_x_continuous(labels=scaleFUN) 
            }

            ## theme
            if (theme.bw == TRUE) {
                p <- p + theme_bw() 
            }

            # horizontal 0 line
            p <- p + geom_hline(yintercept = 0, colour=lcolor,size = lwidth)
            # vertical 0 line
            if (length(xlim)!=0) {
                if (xlim[2]>=1) {
                    p <- p + geom_vline(xintercept = 0, colour=lcolor,size = lwidth)
                }
            } else {
                p <- p + geom_vline(xintercept = 0, colour=lcolor,size = lwidth)
            }

            ## legend and axes
            p <- p + theme(legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = cex.legend),
               legend.position = legend.pos,
               legend.background = element_rect(fill="transparent",colour=NA),
               axis.title=element_text(size=cex.lab),
               axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
               axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
               axis.text = element_text(color="black", size=cex.axis),
               axis.text.x = element_text(size = cex.axis, angle = angle, hjust=x.h, vjust=x.v),
               axis.text.y = element_text(size = cex.axis),
               plot.title = element_text(size = cex.main, hjust = 0.5, face="bold", margin = margin(10, 0, 10, 0)))

            ## grid
            if (gridOff == TRUE) {
                p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
            }


           
            if (placeboTest == FALSE) {
                ## point estimates
                p <- p + geom_line(data = data, aes(time, ATT.ON), size = 1.2)
                ## CIs
                p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower, ymax=CI.upper),alpha=0.2)
                ## wald
                if (wald == TRUE) {
                    p.label <- paste0("Wald p value: ", sprintf("%.3f",x$wald$p))
                    wald <- TRUE
                } 
            } else {
                ## point estimates
                p <- p + geom_line(data = data, aes(time, ATT.ON3), size = 0.7)
                p <- p + geom_line(data = data, aes(time, ATT.ON4), size = 0.7, colour = "#4671D5")
                if (vis == "connected") {
                    p <- p + geom_point(data = data, aes(time, ATT.ON), size = 1.2)
                    p <- p + geom_point(data = data, aes(time, ATT.ON2), size = 1.2, colour = "#4671D5")
                }
                ## CIs
                p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower3, ymax=CI.upper3),alpha=0.2)
                p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower4, ymax=CI.upper4),alpha=0.2, fill = "#0000FF")
                ## p value
                p.label <- paste0("p value: ", sprintf("%.3f",x$est.placebo[5]))
            }       
                        
            
            ## mark p value
            if ((wald == TRUE | placeboTest == TRUE) & p.value == TRUE) {
                if (is.null(text.pos)) {
                    text.pos[1] <- min(data[,"time"], na.rm = 1)
                    text.pos[2] <- ifelse(is.null(ylim), max(data[,"CI.upper"], na.rm = 1), ylim[2])
                }
                p <- p + annotate("text", x = text.pos[1], y = text.pos[2], 
                    label = p.label, size = cex.text, hjust = 0)        
            }
            
            if (count == TRUE) {
                data[,"xmin"] <- data[,"time"] - 0.2
                data[,"xmax"] <- data[,"time"] + 0.2
                data[,"ymin"] <- rep(rect.min, dim(data)[1])
                data[,"ymax"] <- rect.min + (data[,"count"]/max.count) * 0.8 * rect.length
                p <- p + geom_rect(data = data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey70", colour = "grey69", alpha = 0.4, size = 0.2)
            }

        } else {
            ## switch-off effects
            if (is.null(x$est.att.off) == TRUE) { 
                
                cat("Uncertainty estimates not available.\n")
                data <- cbind.data.frame(time, ATT.OFF = x$att.off, count = x$count.off)[show,]
                ## rect.length <- (max(data[,"ATT.OFF"]) - min(data[,"ATT.OFF"]))/2 
                if (length(ylim) != 0) {
                    rect.length <- (ylim[2] - ylim[1]) / 5
                    rect.min <- ylim[1]
                } else {
                    rect.length <- (max(data[,"ATT.OFF"]) - min(data[,"ATT.OFF"]))/2 
                    rect.min <- min(data[,"ATT.OFF"]) - rect.length
                }  

            } else {                
                tb <- x$est.att.off
                data <- cbind.data.frame(time, tb, count = x$count.off)[show,]
                ## rect.length <- (max(data[,"CI.upper"]) - min(data[,"CI.lower"]))/2
                if (length(ylim) != 0) {
                    rect.length <- (ylim[2] - ylim[1]) / 5
                    rect.min <- ylim[1]
                } else {
                    rect.length <- (max(data[,"CI.upper"]) - min(data[,"CI.lower"]))/2
                    rect.min <- min(data[,"CI.lower"]) - rect.length
                }  
            }
             
            ## plotting
            p <- ggplot(data) +
                    geom_vline(xintercept = 0, colour=lcolor,size = 2) +
                    geom_hline(yintercept = 0, colour=lcolor,size = 2) +
                    ## annotate("rect", xmin= time.bf, xmax= Inf,
                    ##          ymin=-Inf, ymax=Inf, alpha = .1,
                    ##          fill = "yellow") +
                    xlab(xlab) +  ylab(ylab) + 
                    scale_x_continuous(labels=scaleFUN) 

            if (theme.bw == TRUE) {
                p <- p + theme_bw() 
            } 

            if (gridOff == TRUE) {
                p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
            }   

            p <- p + theme(legend.position = legend.pos,
             axis.title=element_text(size=cex.lab),
             axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
             axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
             axis.text = element_text(color="black", size=cex.axis),
             axis.text.x = element_text(size = cex.axis, angle = angle, hjust=x.h, vjust=x.v),
             axis.text.y = element_text(size = cex.axis),
             plot.title = element_text(size = cex.main,
               hjust = 0.5,
               face="bold",
               margin = margin(10, 0, 10, 0)))
                    
           
            ## point estimates
            p <- p + geom_line(aes(time, ATT.OFF), size = 1.2)
             
            ## confidence intervals
            if (is.null(x$est.att.off)==FALSE) {
                p <- p + geom_ribbon(aes(x = time, ymin=CI.lower, ymax=CI.upper),alpha=0.2)
            }

            if (count == TRUE) {
                data[,"xmin"] <- data[,"time"] - 0.2
                data[,"xmax"] <- data[,"time"] + 0.2
                data[,"ymin"] <- rep(rect.min, dim(data)[1])
                data[,"ymax"] <- rect.min + (data[,"count"]/max.count) * 0.8 * rect.length
                p <- p + geom_rect(data = data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey70", colour = "grey69", alpha = 0.4, size = 0.2)
            }

        }
            
        ## title
        if (is.null(main) == TRUE) {
            p <- p + ggtitle(maintext)
        } else if (main!=""){
            p <- p + ggtitle(main)
        }

        ## ylim
        if (is.null(ylim) == FALSE) {
            p <- p + coord_cartesian(ylim = ylim)
        }
        suppressWarnings(print(p))
    }
    else if (type == "raw") {

        ## axes labels
        if (is.null(xlab)==TRUE) {
            xlab <- index[2]
        } else if (xlab == "") {
            xlab <- NULL
        }
        if (is.null(ylab)==TRUE) {
            ylab <- Yname
        } else if (ylab == "") {
            ylab <- NULL
        }

        if (is.null(ylim) == TRUE && binary == FALSE) {
            ylim <- c(min(c(Y[show,]), na.rm = TRUE), max(c(Y[show,]), na.rm = TRUE))
        }

        if (is.null(main) == TRUE) {
            main <- "Raw Data"
        }

        set.labels = c("Control", "Treatment") 
      
        if (1 %in% unit.type) {
            co.pos <- which(unit.type == 1)
            Nco <- length(co.pos)
            data1 <- cbind.data.frame("time" = c(rep(time[show], Nco)),
                                      "outcome" = c(Y[show, co.pos]),
                                      "type" = c(rep("co", (Nco*nT))),
                                      "id" = c(rep(1:Nco, each = nT)))
            limits1 <- c("co", "tr")
            #limits1 <- "co"
            colors1 <- c("#99999950", "#FC8D6280")
            #colors1 <- "#99999950"
            ## main1 <- main
            main1 <- "Always Under Control"
        }

        if (2 %in% unit.type) {
            tr.pos <- which(unit.type == 2)
            Ntr <- length(tr.pos)
            data2 <- cbind.data.frame("time" = c(rep(time[show], Ntr)),
                                      "outcome" = c(Y[show, tr.pos]),
                                      "type" = c(rep("tr",(Ntr*nT))),
                                      "id" = c(rep(1:Ntr,each = nT)))
            limits2 <- c("co", "tr")
            #limits2 <- "tr"
            colors2 <- c("#99999950", "#FC8D6280") 
            #colors2 <- "#FC8D6280" 
            ## main2 <- ifelse(1%in%unit.type, "", main)
            main2 <- "Always Under Treatment"
        }

        if (3 %in% unit.type) {
            rv.pos <- which(unit.type == 3)
            Nrv <- length(rv.pos)

            D.plot <- D

            D.plot[which(D.plot == 0)] <- NA
            D.plot[which(I == 0)] <- NA

            D.rv <- as.matrix(D.plot[, rv.pos])
            Y.rv <- as.matrix(Y[, rv.pos])

            Y.trt <- Y.rv * D.rv
            Y.trt.show <- as.matrix(Y.trt[show,])
            time.trt.show <- time[show]
            ut.time <- ut.id <- NULL
            for (i in 1:Nrv) {
                if (sum(is.na(Y.trt.show[,i])) != nT) {
                    ut.id <- c(ut.id, rep(i, nT - sum(is.na(Y.trt.show[,i]))))
                    ut.time <- c(ut.time, time.trt.show[which(!is.na(Y.trt.show[,i]))])
                }
            }


            data3 <- cbind.data.frame("time" = c(rep(time[show], Nrv), ut.time),
                                      "outcome" = c(c(Y[show, rv.pos]),
                                                  c(Y.trt.show[which(!is.na(Y.trt.show))])),
                                      "type" = c(rep("co",(Nrv*nT)),
                                               rep("tr",length(ut.id))),
                                      "id" = c(rep(1:Nrv,each = nT), ut.id))
            ## remove NA
            if (sum(is.na(data3[,"outcome"])) > 0) {
                data3 <- data3[-which(is.na(data3[,"outcome"])),]
            }
            limits3 <- c("co", "tr")
            colors3 <- c("#99999950", "#FC8D6280")
            main3 <- "Treatment Status Changed"
        }
            
        subplot <- function (data, limits, labels, colors, main, binary, theme.bw, gridOff) {
            ## binary
            if (binary == TRUE) {
                data$outcome <- factor(data$outcome, ordered = 1)
            }

            ## theme
            p <- ggplot(data) + xlab(xlab) +  ylab(ylab) 

            ## legend
            set.limits = limits
            set.colors = colors
            set.linetypes = rep("solid", length(limits))
            set.linewidth = rep(0.5, length(limits))

            if (theme.bw == TRUE) {
                p <- p + theme_bw() +
                theme(legend.margin = margin(c(0, 5, 5, 0)),
                      legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = cex.legend),
                      axis.title=element_text(size=cex.lab),
                      axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
                      axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
                      axis.text = element_text(color="black", size=cex.axis),
                      axis.text.x = element_text(size = cex.axis, angle = angle, hjust=x.h, vjust=x.v),
                      axis.text.y = element_text(size = cex.axis),
                      plot.title = element_text(size = cex.main.sub,
                                                hjust = 0.5,
                                                face="bold",
                                                margin = margin(10, 0, 10, 0)))
                if (gridOff == TRUE) {
                    p <- p + theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank())
                }
            } else {
                p <- p +  theme(legend.margin = margin(c(0, 5, 5, 0)),
                                legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = cex.legend),
                                axis.title=element_text(size=cex.lab),
                                axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
                                axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
                                axis.text = element_text(color="black", size=cex.axis),
                                axis.text.x = element_text(size = cex.axis, angle = angle, hjust=x.h, vjust=x.v),
                                axis.text.y = element_text(size = cex.axis),
                                plot.title = element_text(size = cex.main.sub,
                                                          hjust = 0.5,
                                                          face="bold",
                                                          margin = margin(10, 0, 10, 0)))
            }

            ## main
            if (binary == FALSE) {
                p <- p + geom_line(aes(time, outcome,
                                   colour = type,
                                   size = type,
                                   linetype = type,
                                   group = id))

                p <- p + scale_colour_manual(limits = set.limits,
                                             labels = set.labels,
                                             values =set.colors) +
                        scale_linetype_manual(limits = set.limits,
                                              labels = set.labels,
                                              values = set.linetypes) +
                        scale_size_manual(limits = set.limits,
                                          labels = set.labels,
                                          values = set.linewidth) +
                        guides(linetype = guide_legend(title=NULL, nrow=1),
                               colour = guide_legend(title=NULL, nrow=1),
                               size = guide_legend(title=NULL, nrow=1)) 

            } else {
                p <- p + geom_jitter(width = 0.15, height = 0.15, shape = 1,
                                     aes(x = time, y = outcome, colour = type))

                p <- p + scale_colour_manual(limits = set.limits,
                                             labels = set.labels,
                                             values =set.colors) +
                         guides(linetype = guide_legend(title=NULL, nrow=1),
                                colour = guide_legend(title=NULL, nrow=1),
                                size = guide_legend(title=NULL, nrow=1))
            }            
        
            if (!is.numeric(time.label)) {
                p <- p + 
                    scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
            }

            ## title
            if (is.null(main) == TRUE) {
                p <- p + ggtitle("Raw Data")
            } else if (main!="") {
                p <- p + ggtitle(main)
            }

            ## ylim
            if (is.null(ylim) == FALSE) {
                p <- p + coord_cartesian(ylim = ylim)
            }
            return(p)
        }

        cex.main.top <- cex.main

        if (length(unique(unit.type))==1) {
            if (1%in%unit.type) {
                p1 <- subplot(data1, limits1, labels1, colors1, main1, binary, theme.bw, gridOff)
                if (legend.pos != "none") {
                    suppressWarnings(g <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs)
                    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                    suppressWarnings(grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                     legend, nrow = 2, heights = c (1, 1/5)),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2))))
                } else {
                    suppressWarnings(grid.arrange(p1 + theme(legend.position="none"),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2))))
                }   
            }
            else if (2%in%unit.type) {
                p2 <- subplot(data2, limits2, labels2, colors2, main2, binary, theme.bw, gridOff)
                if (legend.pos != "none") {
                    suppressWarnings(g <- ggplotGrob(p2 + theme(legend.position="bottom"))$grobs)
                    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                    suppressWarnings(grid.arrange(arrangeGrob(p2 + theme(legend.position="none"),
                                     legend, nrow = 2, heights = c (1, 1/5)),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2)))) 
                } else {
                    suppressWarnings(grid.arrange(p2 + theme(legend.position="none"),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2))))
                }   
            }
            else if (3%in%unit.type) {
                p3 <- subplot(data3, limits3, labels3, colors3, main3, binary, theme.bw, gridOff)
                if (legend.pos != "none") {
                    suppressWarnings(g <- ggplotGrob(p3 + theme(legend.position="bottom"))$grobs)
                    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                    suppressWarnings(grid.arrange(arrangeGrob(p3 + theme(legend.position="none"),
                                     legend, nrow = 2, heights = c (1, 1/5)),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2)))) 
                } else {
                    suppressWarnings(grid.arrange(p3 + theme(legend.position="none"),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2))))
                }  
            }
        }
        else if (length(unique(unit.type)) == 2) {
            if (!1%in%unit.type) {
                p2 <- subplot(data2, limits2, labels2, colors2, main2, binary, theme.bw, gridOff)
                p3 <- subplot(data3, limits3, labels3, colors3, main3, binary, theme.bw, gridOff)
                if (legend.pos != "none") {
                    suppressWarnings(g <- ggplotGrob(p2 + theme(legend.position="bottom"))$grobs)
                    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                    suppressWarnings(grid.arrange(arrangeGrob(p2 + theme(legend.position="none"), p3 + theme(legend.position="none"),
                                     legend, nrow = 3, heights = c (1, 1, 1/5)),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2))))  
                } else {
                    suppressWarnings(grid.arrange(p2 + theme(legend.position="none"),
                                     p3 + theme(legend.position="none"),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2))))
                }  
            }
            else if (!2%in%unit.type) {
                p1 <- subplot(data1, limits1, labels1, colors1, main1, binary, theme.bw, gridOff)
                p3 <- subplot(data3, limits3, labels3, colors3, main3, binary, theme.bw, gridOff)
                if (legend.pos != "none") {
                    suppressWarnings(g <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs)
                    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                    suppressWarnings(grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p3 + theme(legend.position="none"),
                                     legend, nrow = 3, heights = c (1, 1, 1/5)),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2))))  
                } else {
                    suppressWarnings(grid.arrange(p1 + theme(legend.position="none"),
                                     p3 + theme(legend.position="none"),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2))))
                }  
            }
            else if (!3%in%unit.type) {
                p1 <- subplot(data1, limits1, labels1, colors1, main1, binary, theme.bw, gridOff)
                p2 <- subplot(data2, limits2, labels2, colors2, main2, binary, theme.bw, gridOff)
                if (legend.pos != "none") {
                    suppressWarnings(g <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs)
                    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                    suppressWarnings(grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"),
                                     legend, nrow = 3, heights = c (1, 1, 1/5)),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2))))  
                } else {
                    suppressWarnings(grid.arrange(p1 + theme(legend.position="none"),
                                     p2 + theme(legend.position="none"),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2))))
                }   
            }
        }
        else {
            p1 <- subplot(data1, limits1, labels1, colors1, main1, binary, theme.bw, gridOff)
            p2 <- subplot(data2, limits2, labels2, colors2, main2, binary, theme.bw, gridOff)
            p3 <- subplot(data3, limits3, labels3, colors3, main3, binary, theme.bw, gridOff)
            if (legend.pos != "none") {
                suppressWarnings(g <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs)
                legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                suppressWarnings(grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"),
                                 p3 + theme(legend.position="none"), legend, nrow = 4, heights = c (1, 1, 1, 1/5)),
                                 top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2))))
            } else {
                suppressWarnings(grid.arrange(p1 + theme(legend.position="none"),
                                     p2 + theme(legend.position="none"),
                                     p3 + theme(legend.position="none"),
                                     top = textGrob(main, gp = gpar(fontsize=cex.main.top,font=2))))
            }
        }
        ## end of raw plot
    }
    else if (type=="missing") {
        
        if (is.null(xlab)==TRUE) {
            xlab <- index[2]
        } else if (xlab == "") {
            xlab <- NULL
        }
        if (is.null(ylab)==TRUE) {
            ylab <- index[1]
        } else if (ylab == "") {
            ylab <- NULL
        }
        if (is.null(main)==TRUE) {
            main <- "Treatment Status"
        } else if (main == "") {
            main <- NULL
        }

        m <- obs.missing
        m <- as.matrix(m[show,which(colnames(m)%in%id)])

        all <- unique(c(m))
        col <- col2 <- breaks <- label <- NULL

        if (1%in%all) {
            col <- c(col,"#06266F")
            col2 <- c(col2, "1"=NA)
            breaks <- c(breaks,1)
            label <- c(label,"Under Treatment")
        }
        if (2%in%all) {
            col <- c(col,"#B0C4DE")
            col2 <- c(col2, "2"=NA)
            breaks <- c(breaks,2)
            label <- c(label,"Under Control")
        }
        if (3%in%all) {
            col <- c(col,"#FFFFFF")
            col2 <- c(col2, "3"=NA)
            breaks <- c(breaks,3)
            label <- c(label,"Missing")
        }
        if (4%in%all) {
            col <- c(col,"#A9A9A9")
            col2 <- c(col2, "4"="red")
            breaks <- c(breaks,4)
            label <- c(label,"Treated (Removed)")
        }
        
        TT <- dim(m)[1]
        N <- dim(m)[2]
        units <- rep(rev(1:N), each = TT)
        period <- rep(1:TT, N)
        res <- c(m)
        data <- cbind.data.frame(units=units, period=period, res=res)
        data[,"res"] <- as.factor(data[,"res"])

        ## check if N >= 200
        if (dim(m)[2] >= 200) {
            if (axis.lab == "both") {
                axis.lab <- "time"
            }
            else if (axis.lab == "unit") {
                axis.lab <- "off"
            }
        }

        ## labels
        N.b <- 1:N
        if (axis.lab == "both") {
            if (length(axis.lab.gap)==2) {
                x.gap <- axis.lab.gap[1]
                y.gap <- axis.lab.gap[2] 
            } else {
                x.gap <- y.gap <- axis.lab.gap[1]
            }
        } else {
            x.gap <- y.gap <- axis.lab.gap[1]
        }

        if (x.gap != 0) {
            T.b <- seq(from = 1, to = length(show), by = (x.gap + 1))
        }
        if (y.gap != 0) {
            N.b <- seq(from = N, to = 1, by = -(y.gap + 1))
        }
        id <- rev(id)
        
        p <- ggplot(data, aes(x = period, y = units,
                              fill = res), position = "identity") 
        p <- p + geom_tile(colour="gray90", size=0.1, stat="identity") 
  
        p <- p +
            labs(x = xlab, y = ylab,  
                title=main) +
            theme_bw() + 
            scale_fill_manual(NA, breaks = breaks, values = col, labels=label)

        if(4%in%all) {
            p <- p + geom_point(aes(colour=res),size=0.5)
            p <- p + scale_color_manual(NA, breaks=breaks,
                                        values=col2, labels=label)
        }

        p <- p +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(fill=NA,color="gray90", size=0.5, linetype="solid"),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.title=element_text(size=cex.lab),
              axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
              axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
              axis.text = element_text(color="black", size=cex.axis),
              axis.text.x = element_text(size = cex.axis, angle = angle, hjust=x.h, vjust=x.v),
              axis.text.y = element_text(size = cex.axis),
              plot.background = element_rect(fill = "gray90"),
              legend.background = element_rect(fill = "gray90"),
              legend.position = legend.pos,
              legend.margin = margin(c(0, 5, 5, 0)),
              legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = cex.legend),
              legend.title=element_blank(),
              plot.title = element_text(size=cex.main, hjust = 0.5,face="bold",margin = margin(8, 0, 8, 0)))

        if (axis.lab == "both") {
            p <- p + scale_x_continuous(expand = c(0, 0), breaks = T.b, labels = time.label[T.b]) +
            scale_y_continuous(expand = c(0, 0), breaks = N.b, labels = id[N.b])
        }
        else if (axis.lab == "unit") {
            p <- p + scale_x_continuous(expand = c(0, 0), breaks = T.b, labels = NULL) +
            scale_y_continuous(expand = c(0, 0), breaks = N.b, labels = id[N.b])            
        }
        else if (axis.lab == "time") {
            p <- p + scale_x_continuous(expand = c(0, 0), breaks = T.b, labels = time.label[T.b]) +
            scale_y_continuous(expand = c(0, 0), breaks = N.b, labels = NULL)
        }
        else if (axis.lab == "off") {
            p <- p + scale_x_continuous(expand = c(0, 0), breaks = 1:length(show), labels = NULL) +
            scale_y_continuous(expand = c(0, 0), breaks = 1:N, labels = NULL)
        }
        
        if(length(all)>=4) {
            p <- p + guides(fill=guide_legend(nrow=2,byrow=TRUE))
        }
        suppressWarnings(print(p))
        ## end of missing plot
    }
}

