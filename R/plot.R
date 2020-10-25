##########
## Plot
##########
# x: a fect object
# type of the plot; axes limits; axes labels; 
# main: whether to show the title;
# id: plot a part of units
plot.fect <- function(x,  
  type = "gap", # gap, equiv, missing
  switch.on = TRUE,
  bound = NULL,
  vis = "connected",
  count = TRUE,
  proportion = 0.3, # control the xlim
  pre.periods = NULL, # for testing
  f.threshold = NULL, # equiv f
  tost.threshold = NULL, # pre-trend and placebo
  effect.bound.ratio = FALSE,
  stats = NULL,       ## "none", "F.p", "F.equiv.p", "placebo.p", "equiv.p"
  stats.labs = NULL,
  main = NULL,
  xlim = NULL, 
  ylim = NULL,
  xlab = NULL, 
  ylab = NULL,
  legendOff = FALSE,
  legend.pos = NULL,
  legend.nrow = NULL,
  legend.labs = NULL,
  stats.pos = NULL,
  theme.bw = TRUE,
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
    
    equiv.p <- NULL
    type.old <- type
    ATT <- ATT2 <- ATT3 <- ATT4 <- NULL
    CI.lower3 <-  CI.lower4 <- CI.upper3 <- CI.upper4 <- NULL

    p <- NULL
    outcome <- NULL ## global variable
    labels1 <- labels2 <- labels3 <- NULL
    ATT.OFF <- ATT.ON <- CI.lower <- CI.upper <- NULL
    xmax <- xmin <- ymax <- ymin <- NULL
    placeboTest <- x$placeboTest
    placebo.period <- x$placebo.period
    binary <- x$binary
    Ftest <- !is.null(x$pre.test)
     

    if (class(x) != "fect") {
        stop("Not a \"fect\" object.")
    }
    if (!type %in% c("missing", "gap","equiv")) {
        stop("\"type\" option misspecified.")
    }
    if (type == "gap" && switch.on == FALSE && is.null(x$att.off)) {
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

    if (!is.null(stats)) {
        if (!placeboTest) {
            for (i in 1:length(stats)) {
                if (!stats[i] %in% c("none", "F.p", "F.equiv.p", "F.stat","equiv.p")) {
                    stop ("Choose \"stats\" from c(\"none\", \"F.stat\", \"F.p\", \"F.equiv.p\", \"equiv.p\").")
                }

            }    
        } else {
            for (i in 1:length(stats)) {
                 if (!stats[i] %in% c("none", "placebo.p", "equiv.p")) {
                    stop ("Choose \"stats\" from c(\"none\", \"placebo.p\", \"equiv.p\").")
                }
            }          
        }
        if ("none" %in% stats) {
            stats <- "none"
        }
    }


    ## show statistics
    if (is.null(stats)== TRUE) {
        if (type == "gap") {
            if (placeboTest == TRUE) {
                stats <- "placebo.p"
            } else {
                stats <- c("none")
            }
        } else {
            stats <- c("F.p","equiv.p")
        }
    }

    if (!("none" %in% stats)) {
        if (is.null(stats.labs)==FALSE) {
            if (length(stats.labs)!=length(stats)) {
                stop("\"stats.lab\" should have the same length as \"stats\".")
            }               
        } else {
            stats.labs <- rep(NA, length(stats)) 
            for (i in 1:length(stats)) {
                if (stats[i] == "F.p") {
                    stats.labs[i] <- "F p-value"
                }
                if (stats[i] == "F.equiv.p") {
                    stats.labs[i] <- "F equivalence p-value"
                }
                if (stats[i] == "F.stat") {
                    stats.labs[i] <- "F statistics"
                }
                if (stats[i] == "placebo.p") {
                    stats.labs[i] <- "Placebo p-value"
                }
                if (stats[i] == "equiv.p") {
                    stats.labs[i] <- ifelse(x$placeboTest, "Placebo equivalence p-value", "Equivalence p-value")
                }
            }
        }
    }
    
    ## gap or equiv
    ytitle <- NULL
    bound.old <- bound
    if (type == "gap") {
        if (is.null(bound) == TRUE) {
            bound.old <- bound <- "none"
        }
        maintext <- "Estimated ATT"
        ytitle <- paste("Effect on",x$Y)

    } else if (type == "equiv") { # equiv plot is a gap plot (with some options changes)
        
        if (is.null(x$est.att)) {
            stop("No uncertainty estimates.\n")
        }

        ## change to 90% CI
        if (placeboTest == FALSE) {
            est.att <- x$est.att
            est.bound <- x$att.bound

            est.att[, "CI.lower"] <- est.bound[, "CI.lower"]
            est.att[, "CI.upper"] <- est.bound[, "CI.upper"]

            x$est.att <- est.att

            if (switch.on == TRUE) {
                if (length(xlim)==0) {
                    xlim <- c(-1e5, 0)
                } else {
                    if (xlim[2]>0) {
                        xlim[2]<-0
                    }
                }
            } else {
                if (length(xlim)==0) {
                    xlim <- c(1, 1e5)
                } else {
                    if (xlim[1]<=0) {
                        xlim[1]<-1
                    }
                }
            }
            maintext <- "Equivalence Test" 
            ytitle <- paste("Residual Average of",x$Y)        

        } else {
            maintext <- "Estimated ATT"
            ytitle <- paste("Effect on",x$Y)
        }

        bound <- "both"
        type <- "gap"
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
    if (!is.null(stats.pos)) {
       if (length(stats.pos) != 2) {
           stop(" \"stats.pos\" must be of length 2. ")
       }
    }

    ##-------------------------------##
    ## Plotting
    ##-------------------------------## 
    show.T0 <- which(x$time == 0)
    if (switch.on == FALSE) {
        show.T0 <- which(x$time.off == 0)
    }

    Y <- x$Y.dat
    D <- x$D.dat
    I <- x$I.dat
    Yname <- x$Y

    index <- x$index
    unit.type <- x$unit.type
    obs.missing <- x$obs.missing
    tname <- time <- x$rawtime

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
        if (switch.on == TRUE) {
            time <- x$time
            count.num <- x$count
            best.pos <- 1
        } else if (switch.on == FALSE) {
            time <- x$time.off
            count.num <- x$count.off
            best.pos <- 0
        }
        max.count <- max(count.num)
        ## bound <- "none"
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
    if (type == "gap") {
        
        if (is.null(proportion) == TRUE) {
            show.count <- 1:time.end
        } else {    
            show.count <- which(count.num >= max.count * proportion)
        }
        # which periods to be shown
        show <- intersect(show.count, show.time) 

        # maximum number of cases to be shown
        max.count <- max(count.num[show])
        
        # where on x-axis to show the number
        max.count.pos <- time[intersect(show,which(count.num == max.count))]
        
        if (length(max.count.pos)>1) {
            if (best.pos %in% max.count.pos) {
                max.count.pos <- best.pos
            } else if ((1-best.pos) %in% max.count.pos) {
                max.count.pos <- 1-best.pos
            } else {
                max.count.pos <- max.count.pos[1]
            }
        }

    } else {
        show <- show.time
    }
    

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

    # get equivalence p values for two-one-sided-t tests
    tost <- function(coef, se, range) {
        z <- coef/se
        p1 <- 1 - pnorm((-range[1]+coef)/se) # left bound
        p2 <- 1 - pnorm((range[2]-coef)/se) # right bound
        tost.p <- max(p1,p2)
        return(tost.p)
    }

    ############  START  ############### 
    if (type == "gap") {

        ## axes labels
        if (is.null(xlab) == TRUE) {
            if (switch.on == TRUE) {
                xlab <- paste("Time Relative to the Treatment")
            } else {
                xlab <- paste("Time Relative to Exiting the Treatment")
            }            
        } else if (xlab == "") {
            xlab <- NULL
        }
    
        if (is.null(ylab) == TRUE) {
            ylab <- ytitle
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

        ## bound
        data2 <- NULL
        if (bound != "none") {
            if (is.null(x$est.att)) {
                cat("No uncertainty estimates.\n")
                bound <- "none"
            }
            #if (placeboTest == TRUE) {
            #    bound <- "none"
            #}
        } 

        if (bound != "none"|| "equiv.p" %in% stats) {

            time0 <- NULL

            if (switch.on == TRUE) {
                if (sum(time[show] <= 0) == 0) {
                    cat("No pre-treatment periods are to be plotted.\n")
                    time0 <- 1:length(time[show])
                } else {
                    time0 <- which(time[show] <= 0)
                }
                att.sub <- as.matrix(x$att.bound[show, ])
                ## att.sub <- as.matrix(x$est.att[show, ])

            } else {
                if (sum(time[show] > 0) == 0) {
                    cat("No non-treatment periods are to be plotted.\n")
                    time0 <- 1:length(time[show])
                } else {
                    time0 <- which(time[show] >= 1)
                }
                att.sub <- as.matrix(x$att.off.bound[show, ])
            }

            minBound <- max(abs(att.sub[time0, c("CI.lower", "CI.upper")]))
            minbound <- c(-minBound, minBound)

            ## equivalence range
            if (is.null(f.threshold)==TRUE) {
                f.threshold <- x$test.out$tost.threshold
            }
            if (is.null(tost.threshold)==TRUE) {
                tost.threshold <- x$test.out$tost.threshold
            }
            equiv.range <- c(-1 * tost.threshold, tost.threshold)

            ## recover bound type
            if (placeboTest == TRUE) {
                if (is.null(bound.old)) {
                    bound.old <- "none"
                }
            } else {
                if (is.null(bound.old)) {
                    bound.old <- "both"
                }
            }

            ## change ci if plot bounds
            if (bound.old != "none") {
                if (!is.null(x$est.att)) {
                    c.est.att <- x$est.att
                    c.est.att[, "CI.lower"] <- x$att.bound[, "CI.lower"]
                    c.est.att[, "CI.upper"] <- x$att.bound[, "CI.upper"]
                    x$est.att <- c.est.att
                }
                if (!is.null(x$est.att.off)) {
                    c.est.att <- x$est.att.off
                    c.est.att[, "CI.lower"] <- x$att.off.bound[, "CI.lower"]
                    c.est.att[, "CI.upper"] <- x$att.off.bound[, "CI.upper"]
                    x$est.att.off <- c.est.att
                }
            }

            # bound period
            bound.time <- time[show]
            if (switch.on == TRUE) {
                bound.time <- bound.time[which(bound.time <= 0)]
            } else {
                bound.time <- bound.time[which(bound.time >= 1)]
            }
            
            ## add legend for 95\% CI
            set.limits <- "ci"
            if (is.null(legend.labs)==TRUE) {
                if (bound.old != "none") {
                    set.labels <- "Residual Average (w/ 90% CI)" 
                } else {
                    set.labels <- "ATT (w/ 95% CI)" 
                }
                                   
            } else {
                set.labels <- legend.labs
            }
            set.colors <- "#000000FF"
            set.linetypes <- "solid"
            set.size <- 1

            ## create a dataset for bound
            if (bound.old == "equiv") {
                data2 <- cbind.data.frame(c(rep(equiv.range, each = length(bound.time))))
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
            else if (bound.old == "min") {
                data2 <- cbind.data.frame(c(rep(minbound, each = length(bound.time))))
                names(data2) <- "bound"
                data2$time <- rep(bound.time, 2)
                data2$type <- rep(c("min"), 2 * length(bound.time))
                data2$id <- rep(1:2, each = length(bound.time))
                set.limits <- c(set.limits, "min")
                if (is.null(legend.labs)==TRUE) {
                    set.labels <- c(set.labels, "Min. Range")
                } else {
                    set.labels <- legend.labs
                }
                set.colors <- c(set.colors, "gray50")
                set.linetypes <- c(set.linetypes, "dashed")
                set.size <- c(set.size, 0.7)
            }
            else if (bound.old == "both") {
                data2 <- cbind.data.frame(c(rep(minbound, each = length(bound.time)), rep(equiv.range, each = length(bound.time))))
                names(data2) <- "bound"
                data2$time <- rep(bound.time, 4)
                data2$type <- rep(c("min", "equiv"), each = 2 * length(bound.time))
                data2$id <- rep(1:4, each = length(bound.time))
                set.limits <- c(set.limits, "min", "equiv")
                if (is.null(legend.labs)==TRUE) {
                    set.labels <- c(set.labels, "Min. Range", "Equiv. Range")
                } else {
                    set.labels <- legend.labs
                }
                set.colors <- c(set.colors, "gray50", "red")
                set.linetypes <- c(set.linetypes, "dashed", "dashed")
                set.size <- c(set.size, 0.7, 0.7)
            }
        } 

        CI <- NULL
        if (switch.on == TRUE) {
            if (is.null(x$est.att)==TRUE) {
                CI <- FALSE
            } else {
                CI <- TRUE
            }
        } else if (switch.on == FALSE) {
            if (is.null(x$est.att.off)==TRUE) {
                CI <- FALSE
            } else {
                CI <- TRUE
            }           
        }
        

        ## data frame for main estimates
        if (switch.on == TRUE) {            

            ## switch-on effect
            if (CI==FALSE) {               

                data <- cbind.data.frame(time, ATT = x$att, count = count.num)[show,]                

            } else {

                tb <- x$est.att
                data <- cbind.data.frame(time, tb)[show,]
                colnames(data)[2] <- "ATT"

                if (placeboTest == TRUE) {
                    data[,"ATT4"] <- data[,"ATT3"] <- data[,"ATT2"] <- data[,"ATT"]
                    data[,"CI.lower4"] <- data[,"CI.lower3"] <- data[,"CI.lower2"] <- data[,"CI.lower"]
                    data[,"CI.upper4"] <- data[,"CI.upper3"] <- data[,"CI.upper2"] <- data[,"CI.upper"]                    
                    pos1 <- intersect(which(data[,"time"] >= (placebo.period[1] + 1)), which(data[,"time"] <= (placebo.period[2] - 1)))
                    pos2 <- c(which(data[,"time"] <= (placebo.period[1] - 1)), which(data[,"time"] >= (placebo.period[2] + 1)))
                    pos3 <- intersect(which(data[,"time"] >= (placebo.period[1])), which(data[,"time"] <= (placebo.period[2] - 1)))
                    pos4 <- c(which(data[,"time"] <= (placebo.period[1] - 2)), which(data[,"time"] >= (placebo.period[2] + 1)))
                    data[pos1, c("ATT","CI.lower","CI.upper")] <- NA
                    data[pos2, c("ATT2","CI.lower2","CI.upper2")] <- NA
                    data[pos3, c("ATT3","CI.lower3","CI.upper3")] <- NA
                    data[pos4, c("ATT4","CI.lower4","CI.upper4")] <- NA
                }
            }

        }        

        # height of the histogram
        if (CI == FALSE) {
            cat("Uncertainty estimates not available.\n")
            if (length(ylim) != 0) {
                rect.length <- (ylim[2] - ylim[1]) / 5
                rect.min <- ylim[1]
            } else {
                rect.length <- (max(data[,"ATT"]) - min(data[,"ATT"]))/2
                rect.min <- min(data[,"ATT"]) - rect.length
            } 
        } else {
            if (length(ylim) != 0) {
                rect.length <- (ylim[2] - ylim[1]) / 5
                rect.min <- ylim[1]
            } else {
                rect.length <- (max(data[,"CI.upper"]) - min(data[,"CI.lower"]))/2
                rect.min <- min(data[,"CI.lower"]) - rect.length
            }  
        } 

        ## plotting
        if (bound.old == "none") {

            p <- ggplot(data)             

        } else  { ## with bounds
            
            p <- ggplot(data2) 

            p <- p + geom_line(aes(time, bound, colour = type, linetype = type, size = type, group = id)) 
            
            ## legends for bounds
            if (is.null(legend.nrow) == TRUE) {
                legend.nrow <- ifelse(length(set.limits) <= 3, 1, 2)    
            } 

            p <- p + scale_colour_manual(limits = set.limits, labels = set.labels, values =set.colors) +
                scale_size_manual(limits = set.limits, labels = set.labels, values = set.size) +  
                scale_linetype_manual(limits = set.limits, labels = set.labels, values = set.linetypes) +
                guides(linetype = guide_legend(title=NULL, nrow=legend.nrow), colour = guide_legend(title=NULL, nrow=legend.nrow),
                size = guide_legend(title=NULL, nrow=legend.nrow)) 

            if (effect.bound.ratio == TRUE) {
                if (is.null(stats.pos)) {
                    stats.pos[1] <- min(data[,"time"], na.rm = 1)
                    stats.pos[2] <- ifelse(is.null(ylim), max(data[,"CI.upper"], na.rm = 1), ylim[1])
                }
                p.label <- paste("ATT / Min. Range = ", sprintf("%.3f",x$att.avg / minBound), sep="")
                p <- p + annotate("text", x = stats.pos[1], y = stats.pos[2], 
                    label = p.label, size = cex.text, hjust = 0)                
            }
        } 

        ## xlab and ylab 
        p <- p + xlab(xlab) +  ylab(ylab) + scale_x_continuous(labels=scaleFUN)

        ## theme
        if (theme.bw == TRUE) {
            p <- p + theme_bw() 
        }


        ## grid
        if (gridOff == TRUE) {
            p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }

        # horizontal 0 line
        p <- p + geom_hline(yintercept = 0, colour=lcolor,size = lwidth)
            # vertical 0 line
        if (length(xlim)!=0) {
            if ((xlim[2]>=1 & switch.on == TRUE) | (xlim[1]<=0 & switch.on == FALSE)) {
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


        ## point estimates
        if (placeboTest == FALSE) {
            ## point estimates
            p <- p + geom_line(data = data, aes(time, ATT), size = 1.2)
            ## CIs
            if (CI == TRUE) {
                p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower, ymax=CI.upper),alpha=0.2)                
            }
            ## stats
            p.label <- NULL
            for (i in 1:length(stats)) {
                if ("F.p" %in% stats[i]) {
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",x$test.out$f.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("F.stat" %in% stats[i]) {
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",x$test.out$f.stat))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("F.equiv.p" %in% stats[i]) {
                    # calculate new p value
                    if (is.null(f.threshold)==FALSE | is.null(proportion) == FALSE | is.null(pre.periods) == FALSE) {
                        test.out <- diagtest(x, proportion = proportion, pre.periods = pre.periods, f.threshold = f.threshold)
                        f.equiv.p <- test.out$f.equiv.p
                    } else {
                        f.equiv.p <- x$test.out$f.equiv.p
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",f.equiv.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("equiv.p" %in% stats[i]) {
                    # calculate new p value
                    if (is.null(tost.threshold)==FALSE | is.null(proportion) == FALSE | is.null(pre.periods) == FALSE) {
                        test.out <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
                        tost.equiv.p <- test.out$tost.equiv.p
                    } else {
                        tost.equiv.p <- x$test.out$tost.equiv.p
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",tost.equiv.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
            }                
               

        } else { # placebo plot
            ## point estimates
            p.label <- NULL
            p <- p + geom_line(data = data, aes(time, ATT3), size = 0.7)
            p <- p + geom_line(data = data, aes(time, ATT4), size = 0.7, colour = "#4671D5")
            if (vis == "connected") {
                p <- p + geom_point(data = data, aes(time, ATT), size = 1.2)
                p <- p + geom_point(data = data, aes(time, ATT2), size = 1.2, colour = "#4671D5")
            }
            ## CIs
            p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower3, ymax=CI.upper3),alpha=0.2)
            p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower4, ymax=CI.upper4),alpha=0.2, fill = "#0000FF")
            
            ## stats
            for (i in 1:length(stats)) {
                if ("placebo.p" %in% stats[i]) {
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",x$test.out$placebo.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("equiv.p" %in% stats[i]) {
                    p.label1 <- NULL
                    if (is.null(tost.threshold)==FALSE | is.null(proportion) == FALSE| is.null(pre.periods) == FALSE) {
                        test.out <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
                        placebo.equiv.p <- test.out$placebo.equiv.p
                    } else {
                        placebo.equiv.p <- x$test.out$placebo.equiv.p
                    }
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f", placebo.equiv.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }                
            }
        }        

        ## mark stats
        hpos <- ifelse(switch.on == TRUE, 0, 1)
        if ("none" %in% stats == FALSE) {
            if (is.null(stats.pos)) {                
                if (switch.on == TRUE) {
                    stats.pos[1] <- min(data[,"time"], na.rm = 1)
                    ## hpos <- 0
                } else {
                    stats.pos[1] <- max(data[,"time"], na.rm = 1)
                    ## hpos <- 1
                }
                ci.top <- max(data[,"CI.upper"], na.rm = 1)
                stats.pos[2] <- ifelse(is.null(ylim), ci.top, ylim[2]) 
            }
            if (!is.null(p.label)) {
                p <- p + annotate("text", x = stats.pos[1], y = stats.pos[2], 
                              label = p.label, size = cex.text, hjust = hpos, vjust = "top")
            }            
        }    

        ## histogram
        if (count == TRUE) {
            data[,"xmin"] <- data[,"time"] - 0.2
            data[,"xmax"] <- data[,"time"] + 0.2
            data[,"ymin"] <- rep(rect.min, dim(data)[1])
            data[,"ymax"] <- rect.min + (data[,"count"]/max.count) * 0.8 * rect.length
            xx <- range(data$time)
            p <- p + geom_rect(data = data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                fill = "grey70", colour = "grey69", alpha = 0.4, size = 0.2)
            p <- p + annotate("text", x = max.count.pos - 0.02 * (xx[2]-xx[1]), 
                y = max(data$ymax) + 0.2 * rect.length, 
                label = max.count, size = cex.text * 0.8, hjust = 0.5)
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
    }


    ## missing/treatment plot
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
        ## end of missing plot
    }

    # if (!is.null(equiv.p)) {
    #     return(list(p = p, equiv.p = equiv.p))
    # } else {
    #     return(list(p = p))
    # }
    return(suppressWarnings(print(p)))
    
}

