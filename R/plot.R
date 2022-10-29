## new plot
# x: a fect object
# type of the plot; axes limits; axes labels; 
# main: whether to show the title;
# id: plot a part of units

plot.fect <- function(x,  
  type = NULL, # gap, equiv, status, exit, factors, loadings, calendar
  loo = FALSE,
  highlight = NULL, ## for carryover test and placebo test
  plot.ci = NULL, ## "0.9", "0.95", "none"
  show.points = NULL,
  show.group = NULL,
  bound = NULL, # "none", "min", "equiv", "both"
  vis = NULL,
  count = TRUE,
  proportion = 0.3, # control the xlim
  pre.periods = NULL, # for testing
  f.threshold = NULL, # equiv f
  tost.threshold = NULL, # pre-trend placebo carryover
  effect.bound.ratio = FALSE,
  stats = NULL,       ## "none", "F.p", "F.equiv.p", "placebo.p", "carryover.p", "equiv.p"
  stats.labs = NULL,
  main = NULL,
  xlim = NULL, 
  ylim = NULL,
  xlab = NULL, 
  ylab = NULL,
  gridOff = FALSE,
  legendOff = FALSE,
  legend.pos = NULL,
  legend.nrow = NULL,
  legend.labs = NULL,
  stats.pos = NULL,
  theme.bw = TRUE,
  nfactors = NULL,
  include.FE = TRUE,
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
  start0 = FALSE,
  return.test = FALSE,
  balance = NULL,
  ...){

    group <- ATT5 <- ATT6 <- CI.lower.90 <- CI.lower6 <- CI.upper.90 <- CI.upper6 <- L1 <- eff <- NULL

    # come needed variables
    equiv.p <- NULL
    ATT <- ATT2 <- ATT3 <- ATT4 <- NULL
    CI.lower3 <-  CI.lower4 <- CI.upper3 <- CI.upper4 <- NULL
    p <- NULL
    outcome <- NULL ## global variable
    labels1 <- labels2 <- labels3 <- NULL
    ATT.OFF <- ATT.ON <- CI.lower <- CI.upper <- NULL
    xmax <- xmin <- ymax <- ymin <- NULL

    # check the class
    if (class(x)[1] != "fect") {
        stop("Not a \"fect\" object.")
    }
    loo.test.out <- test.out <- x$test.out

    if(is.null(balance)){
        if(!is.null(x$balance.period)){
            balance <- TRUE
        }
        else{
            balance <- FALSE
        }
    }


    # check if the input has the loo results
    pequiv <- !is.null(x$pre.est.att) ## if have leave one out pre-treatment results
    if (loo == TRUE && pequiv == FALSE) {
        stop("No leave one out results for pre-treatment periods.")
    }
    else if(loo == TRUE && pequiv == TRUE){
        loo <- 1 
    }
    else{
        loo <- 0
    }

    placeboTest <- x$placeboTest
    placebo.period <- x$placebo.period
    carryoverTest <- x$carryoverTest

    carryover.period <- x$carryover.period
    binary <- x$binary
    Ftest <- !is.null(x$pre.test)

    if(!is.null(show.group)){
        if(length(show.group)>1){
            stop("\"show.group\" should contain only one group.\n")
        }
        all.group.name <- names(x$g.level)
        if(!show.group%in%all.group.name){
            message("The specified group does not exist or its treatment effects cannot be estimated.\n")
            return(0)
        } 
    }

    # check the key option type
    if(!is.null(type)){
        if (!type %in% c("status", "gap","equiv","exit","factors","loadings","calendar","box")) {
            stop("\"type\" option misspecified. Must be one of followings:\"status\",\"gap\",\"equiv\",\"exit\",\"calendar\",\"box\".")
        }
        if (type == "exit" && is.null(x$att.off)) {
            stop("No exiting treatment effect to be plotted.")
        }
    }
    else{ #default setting
        if(loo == 1){
            type <- 'equiv'
        }
        else if(placeboTest == 1){
            type <- "gap"
        }
        else if(carryoverTest == 1){
            type <- "exit"
        }
        else{
            type <- "gap"
        }
    }
    type.old <- type

    # factors, loadings, fe
    if(type %in% c("loadings","factors")){
        if(type == "loadings"){
            if(!x$method %in% c("gsynth","ife","fe")){
                stop("Can't Visualize the Loadings.\n")
            }
            if (x$r.cv==0) {
                stop("No factors are included in the model.\n") 
            } 
            else {
                ## number of loadings to be plotted
                if (is.null(nfactors)==TRUE) {
                    nfactors<-min(x$r.cv,4) 
                } 
                else if (nfactors>x$r.cv) {
                    message("Too many factors specified. ")
                    nfactors<-min(x$r.cv,4) 
                }
                if (nfactors == 1) {
                    message("Loadings for the first factor are shown...\n")
                } 
                else if (nfactors < x$r.cv) {
                    message(paste("Loadings for the first",nfactors,"factors are shown...\n"))
                }

                ## title
                if (is.null(main) == TRUE) {
                    main <- "Factor Loadings"
                } 
                else if (main=="") {
                    main <- NULL
                }

                ## prepare data
                L.hat <- rbind(x$lambda.tr, x$lambda.co)
                Lname <- Llabel <- c()
                r <- x$r.cv
                for (i in 1:r) {
                    Lname<-c(Lname,paste("L",i,sep=""))
                    Llabel<-c(Llabel,paste("Factor",i))
                }
                colnames(L.hat) <- Lname
                rownames(L.hat) <- c()

                if(x$force %in% c(1,3) & include.FE == TRUE){
                    L.hat <- cbind(c(x$alpha.tr,x$alpha.co),L.hat)
                    colnames(L.hat) <- c(paste("Factor",0),Lname)
                    rownames(L.hat) <- c()
                    nfactors <- nfactors + 1
                    Llabel<-c("FE",Llabel)
                }

                data <- cbind.data.frame(L.hat,
                                         "id"=c(x$tr, x$co),
                                         "group"=as.factor(c(rep("Treated",x$Ntr),
                                         rep("Control",x$Nco))))
                
                if (nfactors == 1) {
                    p <- ggplot(data, aes(x=group, y=L1, fill = group)) +
                    geom_boxplot(alpha = 0.7) +
                    coord_flip() + guides(fill=FALSE) +
                    xlab("") + ylab("Factor Loading")  
                } 
                else {
                    if (x$Ntr >= 5) {
                        my_dens <- function(data, mapping, ...) {
                        ggplot(data = data, mapping = mapping) +
                        geom_density(..., alpha = 0.7, color = NA)
                        }
                        p <- GGally::ggpairs(data, mapping = aes(color = group, fill = group),
                        columns = 1:nfactors,
                        columnLabels = Llabel[1:nfactors],
                        diag = list(continuous = my_dens),
                        title = main) +
                        theme(plot.title = element_text(hjust = 0.5))
                    } 
                    else if(x$Ntr > 1) {
                        my_dens <- function(data, mapping, ...) {
                        ggplot(data = data, mapping = mapping) +
                        geom_density(..., fill = "gray", alpha = 0.7, color = "gray50")
                        }
                        p <- GGally::ggpairs(data, mapping = aes(color = group),
                        columns = 1:nfactors,
                        columnLabels = Llabel[1:nfactors],
                        diag = list(continuous = my_dens),
                        title = main)
                    }
                    else{
                        my_dens <- function(data, mapping, ...) {
                        ggplot(data = data, mapping = mapping) +
                        geom_density(..., fill = "gray", alpha = 0.7, color = "gray50")
                        }
                        p <- GGally::ggpairs(data, mapping = aes(color = group),
                        columns = 1:nfactors, upper = 'blank',
                        columnLabels = Llabel[1:nfactors],
                        diag = list(continuous = my_dens),
                        title = main)                     
                    }
                }
                #suppressWarnings(print(p))
                return(p)
            }
        }
        if(type == "factors"){
            if (theme.bw == TRUE) {
                line.color <- "#AAAAAA70"
            } 
            else {
                line.color <- "white"
            }
            if (axis.adjust==TRUE) {
                angle <- 45
                x.v <- 1
                x.h <- 1
            } 
            else {
                angle <- 0
                x.v <- 0
                if (type=="missing") {
                    x.h <- 0.5
                } else {
                    x.h <- 0
                }
            }
            r <- x$r.cv
            if(!x$method %in% c("gsynth","ife")){
                stop("Can't Visualize the Loadings.\n")
            }
            if (x$r.cv==0) {
                stop("No factors are included in the model.\n") 
            }

            time <- x$rawtime
            if (!is.numeric(time[1])) {
                time <- 1:x$T
            }
            if (length(xlim) != 0) {
                show <- which(time>=xlim[1]& time<=xlim[2])
            } 
            else {
                show <- 1:length(time)
            }
            nT <- length(show)
            time.label <- x$rawtime[show]
            F.hat <- x$factor

            if(x$force %in% c(2,3) & include.FE == TRUE){
                F.hat <- cbind(x$xi,F.hat)
                r <- r + 1
            }

            if (x$r.cv==0) {
                message("No factors included in the model.\n")
            } 
            else {
                ## axes labels
                if (is.null(xlab)==TRUE) {
                    xlab <- x$index[2]
                } else if (xlab == "") {
                    xlab <- NULL
                }
                if (is.null(ylab)==TRUE) {
                    ylab <- "Estimate"
                } else if (ylab == "") {
                    ylab <- NULL
                }
                ## title
                if (is.null(main) == TRUE) {
                    main <- "Latent Factors"
                } else if (main=="") {
                    main <- NULL
                }
                ## prepare data
                L.co <- x$lambda.co
                
                if(x$force %in% c(2:3) & include.FE == TRUE){
                    L.co <- cbind(rep(1,dim(L.co)[1]),L.co)
                    r.use <- c(0:(r-1))
                }
                else{
                    r.use <- c(1:r)
                }
                norm<-sqrt(diag(t(L.co)%*%L.co)/(x$N-x$Ntr))
                data <- cbind.data.frame("time" = rep(time[show],r),
                                        "factor" = c(F.hat[show,])*rep(norm,each=nT),
                                        "group" = as.factor(c(rep(r.use,each=nT))))
                ## theme
                p <- ggplot(data) 
                if (theme.bw == TRUE) {
                    p <- p + theme_bw()
                }
                p <- p + xlab(xlab) +  ylab(ylab) + ggtitle(main) +
                    geom_hline(yintercept=0,colour=line.color,size = 2) +
                    theme(legend.position = legend.pos,
                        axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                        plot.title = element_text(size=20,
                                                    hjust = 0.5,
                                                    face="bold",
                                                    margin = margin(10, 0, 10, 0)))  
                ## main plot
                p <- p + geom_line(aes(time, factor,
                                    colour = group,
                                    group = group), size = 1.2)


                brew.colors <- c("black","steelblue","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9")
                set.colors = brew.colors[1:r]
                p <- p + scale_colour_manual(values =set.colors) 

                ## legend
                p <- p + guides(colour = guide_legend(title="Factor(s)", ncol=4)) 

                if (!is.numeric(time.label)) {
                    p <- p + 
                        scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
                }
            
                ## ylim
                if (is.null(ylim) == FALSE) {
                    p <- p + coord_cartesian(ylim = ylim)
                }            
                #suppressWarnings(print(p))
                return(p)
            }

        }
    }


    if(!is.null(show.group)){
        if(is.null(x$group.output[[show.group]]$att.on) & type!='status'){
            stop(paste0(show.group, " doesn't contain treated units. Can't plot.\n"))
        }        
    }        

    if(type=='status' && !is.null(show.group)){
        if(!is.null(id)){
            stop("\"show.group\" can't be used with \"id\" in the status plot.\n")
        }
    }

    if(type == 'gap' | type == 'equiv'){
        carryoverTest <- FALSE
    }

    if(type == 'exit'){
        placeboTest <- FALSE
    }

    if (!is.null(plot.ci)) {
        if(!plot.ci  %in% c("0.9", "0.95", "none")){
            stop("\"plot.ci\" must be one of \"0.95\", \"0.9\" or \"none\".")
        }
        if (plot.ci  %in% c("0.90", "0.95") && is.null(x$est.att)) {
            stop("No uncertainty estimates.")
        }
        if(plot.ci == "0.90" && type%in%c("gap","exit")){
            warning("90% CI in gap plots or exiting treatment plots.\n")
        }
        if(plot.ci == "0.95" && type=="equiv"){
            warning("95% CI in equivalence test plots.\n")
        }
    } 
    else { # default settings for plot.ci
        if(is.null(x$est.att)){
            plot.ci <- "none"
        }
        else if(type=='equiv'){
            plot.ci <- "0.9"
        } 
        else { #gap plot or exiting plot
            plot.ci <- "0.95"
        }
    }

    if(plot.ci == "0.95"){
        plot.ci <- "95"
    }
    if(plot.ci == "0.9"){
        plot.ci <- "90"
    }
    
    if(type=='equiv' && plot.ci=='none'){
        stop("No uncertainty estimates. Can't perform equivalence tests.\n")
    }

    if((placeboTest|carryoverTest) && plot.ci=='none' && type!="status"){
        stop("No uncertainty estimates. Can't perform placebo test or carryover test.\n")
    }

    if (is.null(vis)) {
        if (placeboTest|carryoverTest) {
            vis <- "connected"
        } else {
            vis <- "none"
        }
    }
    else{
        if (!vis %in% c("connected","none")){
            stop("\"vis\" must be \"connected\" or \"none\".")
        }
    }

    if(is.null(show.points)){
        if(placeboTest==TRUE){
            if(length(placebo.period)==1){
                show.points <- TRUE
            }
            else {
                show.points <- TRUE
            }
        }
        else {
            show.points <- TRUE
        }
    }

    if(!show.points %in% c(TRUE,FALSE)){
        stop("\"show.points\" mis-specified.\n")
    }

    if(show.points){
        if(is.null(plot.ci)){
            plot.ci.point <- "both"
        }
        else{
            plot.ci.point <- plot.ci
        }   
        plot.ci <- "both"
    }

    if (is.null(highlight)) {
        if (placeboTest|carryoverTest) {
            highlight <- TRUE
        } else {
            highlight <- FALSE
        }
    }

    if (is.null(bound) == FALSE) {
        if (!bound %in% c("none", "min", "equiv", "both")) {
            stop("\"bound\" option misspecified.")
        } 
    }
    else{ # default settings for bound
        if(type=="equiv"){
            bound <- "both"
        }
        else{
            bound <- "none"
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
        } 
        else {
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
        if (!placeboTest && !carryoverTest) {
            for (i in 1:length(stats)) {
                if (!stats[i] %in% c("none", "F.p", "F.equiv.p", "F.stat","equiv.p")) {
                    stop ("Choose \"stats\" from c(\"none\", \"F.stat\", \"F.p\", \"F.equiv.p\", \"equiv.p\").")
                }
            }    
        } 
        else if(placeboTest) {
            for (i in 1:length(stats)) {
                 if (!stats[i] %in% c("none", "placebo.p", "equiv.p")) {
                    stop ("Choose \"stats\" from c(\"none\", \"placebo.p\", \"equiv.p\").")
                }
            }          
        }
        else if(carryoverTest){ # carry over test
            for (i in 1:length(stats)) {
                 if (!stats[i] %in% c("none", "carryover.p", "equiv.p")) {
                    stop ("Choose \"stats\" from c(\"none\", \"carryover.p\", \"equiv.p\").")
                }
            }  
        }
        if ("none" %in% stats) {
            stats <- "none"
        }
    }
    else{ # default settings for the option stats
        if (type == "gap") {
            if (placeboTest == TRUE) {
                stats <- c("placebo.p", "equiv.p")
            }
            else {
                stats <- c("none")
            }             
        } 
        else if(type == 'equiv') {
            stats <- c("F.p","equiv.p")
            if (placeboTest == TRUE) {
                stats <- c("placebo.p", "equiv.p")
            }
        }
        else if(type=='exit'){
            if(carryoverTest==TRUE){
                stats <- c("carryover.p", "equiv.p")
            }else{
                stats <- "none"
            }
        }
        else{
            stats <- "none"
        }
    }

    if(type=='calendar'){
        stats <- "none"
    }

    # names for all statistics
    if (!("none" %in% stats)) {
        if (is.null(stats.labs)==FALSE) {
            if (length(stats.labs)!=length(stats)) {
                stop("\"stats.lab\" should have the same length as \"stats\".")
            }               
        } 
        else {
            stats.labs <- rep(NA, length(stats)) 
            for (i in 1:length(stats)) {
                if (stats[i] == "F.p") {
                    stats.labs[i] <- "F test p-value"
                }
                if (stats[i] == "F.equiv.p") {
                    stats.labs[i] <- "F equivalence test p-value"
                }
                if (stats[i] == "F.stat") {
                    stats.labs[i] <- "F statistics"
                }
                if (stats[i] == "placebo.p") {
                    stats.labs[i] <- "Placebo test p-value"
                }
                if (stats[i] == "carryover.p") {
                    stats.labs[i] <- "Carryover effect test p-value"
                }
                if (stats[i] == "equiv.p") {
                    if(placeboTest){
                        stats.labs[i] <- "Placebo equivalence test p-value"
                    }
                    else if(carryoverTest){
                        stats.labs[i] <- "Carryover effect equivalence test p-value"
                    }
                    else{
                        stats.labs[i] <- "Equivalence test p-value"
                    }
                }
            }
        }
    }

    # titles; xlim and ylim
    ytitle <- NULL
    bound.old <- bound
    if (type == "gap") { # classic plots or placebo tests
        maintext <- "Estimated ATT"
        ytitle <- paste("Effect on",x$Y)
        if(placeboTest==1){
            maintext <- "Placebo Test"
        }
    }
    else if (type == "equiv") { # equiv plot is a gap plot (with some options changes)
        if (is.null(x$est.att)) {
            stop("No uncertainty estimates.\n")
        }
        
        # classic equivalence test
        if (length(xlim)==0) {
            xlim <- c(-1e5, 0)
        } 
        else {
            if (xlim[2]>0) {
                xlim[2] <- 0
            }
        }
        if (loo==0) { 
            maintext <- "Equivalence Test" 
            #ytitle <- paste("Residual Average of",x$Y)
            ytitle <- paste("Effect on",x$Y)        
        } 
        else {
            # loo equivalence test
            maintext <- "Leave-one-out Equivalence Test"
            ytitle <- paste("Effect on",x$Y)
        }
    }
    else if (type=='exit') {
        maintext <- "Estimated ATT"
        ytitle <- paste("Effect on",x$Y)
        if(carryoverTest==1){
            maintext <- "Carryover Effects"
        }
    }
    else if (type=='calendar'){
        maintext <- "ATT by Calendar Time"
        ytitle <- paste("Effect on",x$Y)
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
        if (type == "status") {
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
        cex.main <- 16 * cex.main
    } else {
        cex.main <- 16
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
        cex.text <- 5 * cex.text
    }  else {
        cex.text <- 5
    }

    ## text label position
    if (!is.null(stats.pos)) {
       if (length(stats.pos) != 2) {
           stop(" \"stats.pos\" must be of length 2. ")
       }
    }

    # key function
    # generate a data frame contains the results for the plots
    # colnames of est.att: ATT S.E. CI.lower CI.upper p.value count CI.lower.90 CI.upper.90
    est.bound <- est.att <- NULL
    if(!is.null(show.group)){
        target.group <- x$est.group.output[[show.group]]
        info.group <- x$group.output[[show.group]]

        x$att <- info.group$att.on
        x$time <- info.group$time.on
        x$count <- info.group$count.on
        x$time.off <- info.group$time.off
        x$count.off <- info.group$count.off

        x$est.att <- target.group$att.on
        x$att.bound <- target.group$att.on.bound
        x$att.boot <- target.group$att.on.boot

        x$est.att.off <- target.group$att.off
        x$att.off.bound <- target.group$att.off.bound

        x$est.placebo <- target.group$att.placebo
        x$est.carryover <- target.group$att.carryover

        if(loo == 1){
            loo.group <- x$pre.est.group.output[[show.group]]
            x$pre.est.att <- loo.group$pre.est.att
            x$pre.att.bound <- loo.group$pre.att.bound
            x$pre.att.boot <- loo.group$pre.att.boot
        }

        if(type == 'status'){
            NN <- dim(x$obs.missing)[2]
            TT <- dim(x$obs.missing)[1]
            T.name <- rownames(x$obs.missing)
            N.name <- colnames(x$obs.missing)
            use.obs.missing <- x$obs.missing[which(x$G==x$g.level[show.group])]
            use.id <- colnames(x$obs.missing)
            use.index <- apply(x$G,2,mean)
            use.id <- use.id[which(use.index==x$g.level[show.group])]
            use.obs.missing <- matrix(use.obs.missing, nrow = TT)
            rownames(use.obs.missing) <- T.name
            colnames(use.obs.missing) <- use.id
            x$obs.missing <- use.obs.missing
        }
    }

    use.balance <- FALSE
    if(!is.null(x$balance.period) & balance==TRUE){
        x$att <- x$balance.att
        x$time <- x$balance.time
        x$count <- x$balance.count
        x$att.avg <- x$balance.att.avg
        x$est.att <- x$est.balance.att
        x$att.bound <- x$balance.att.bound
        x$att.boot <- x$balance.att.boot
        x$est.placebo <- x$est.balance.placebo
        x$obs.missing <- x$obs.missing.balance
        use.balance <- TRUE
    }


    if (!is.null(x$est.att)) { # have uncertainty estimation
        est.att <- x$est.att 
        est.bound <- x$att.bound
        
        colnames(est.bound) <- c("CI.lower.90", "CI.upper.90") ## 90% ci 
        est.att <- cbind(est.att, est.bound)
        pre.est.att <- pre.att.bound <- NULL
        if (loo==1) { ## replace pre-treatment period with loo results 
            pre.est.att <- x$pre.est.att
            pre.att.bound <- x$pre.att.bound 
            colnames(pre.att.bound) <- c("CI.lower.90", "CI.upper.90") ## 90% ci
            pre.est.att <- cbind(pre.est.att, pre.att.bound)
            t0 <- t1 <- NULL 
            t.s <- t.e <- NULL
            t0 <- rownames(pre.est.att)
            t1 <- rownames(est.att)
            t.s <- which(t1 == t0[1])
            t.e <- which(t1 == t0[length(t0)])
            est.att[t.s:t.e ,] <- pre.est.att 
        }
    }

    # est.att.off for the exit plot
    est.bound.off <- est.att.off <- NULL
    if(!is.null(x$est.att.off)){
        est.att.off <- x$est.att.off
        est.bound.off <- x$att.off.bound
        colnames(est.bound.off) <- c("CI.lower.90", "CI.upper.90") ## 90% ci 
        est.att.off <- cbind(est.att.off, est.bound.off)
    }

    ##-------------------------------##
    ## Plotting
    ##-------------------------------##






    show.T0 <- which(x$time == 0)
    if (type == 'exit') {
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

    if (type == "status") {
        if (is.null(id) == TRUE) {
            id <- colnames(obs.missing)
        }
        m.l <- length(id)
        for (i in 1:m.l) {
            if (!id[i]%in%colnames(obs.missing)) {
                stop("Some specified units are not in the data.")
            }
        }
    }
    else { ## raw plot
        if (is.null(id) == TRUE) {
            id <- x$id
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
    if (type == "gap" | type == "equiv") {
        time <- x$time
        count.num <- x$count
        best.pos <- 1
        max.count <- max(count.num)
    }
    else if (type == "exit"){
        time <- x$time.off
        count.num <- x$count.off
        best.pos <- 0
        max.count <- max(count.num)
    }
    else {
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
    if (type %in% c("gap","equiv","exit")) {
        if (is.null(proportion) == TRUE) {
            show.count <- 1:time.end
        } 
        else {    
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
    } 
    else {
        show <- show.time
    }

    if (length(show) < 2 & type %in% c("gap","equiv","exit")) {
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
        p2 <- 1 - pnorm((range[2]-coef)/se)  # right bound
        tost.p <- max(p1,p2)
        return(tost.p)
    }


    if(type %in% c("gap","equiv","exit")){
        if(type=="exit"){
            switch.on <- FALSE
        }
        else{
            switch.on <- TRUE
        }
        
        ## axes labels
        if (is.null(xlab) == TRUE) {
            if (switch.on == TRUE) {
                xlab <- paste("Time since the Treatment Began")
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

        ## y=0 line type
        lcolor <- "white"
        lwidth <- 2
        if (theme.bw == TRUE) {
            lcolor <- "#AAAAAA70"
            lwidth <- 1.5
        }

        ## equivalence range

        if (is.null(f.threshold)==TRUE) {
            f.threshold <- x$test.out$f.threshold
            change.f.threshold <- 0
        }
        else{
            change.f.threshold <- 1
        }

        if (is.null(tost.threshold)==TRUE) {
            if(placeboTest==1 | x$carryoverTest==1){
                tost.threshold <- x$tost.threshold
            }
            else{
                tost.threshold <- x$test.out$tost.threshold
            }
            change.tost.threshold <- 0
        }
        else{
            change.tost.threshold <- 1
        }

        if(proportion == x$proportion){
            change.proportion <- 0
        }
        else{
            change.proportion <- 1
        }

        if(is.null(pre.periods)){
            pre.periods <- x$pre.periods
        }
        else{
            max.count.test <- max(x$count)
            max.pre.periods <- x$time[which(x$count >= max.count.test * proportion & x$time <= 0)]
            pre.periods <- intersect(pre.periods[1]:pre.periods[length(pre.periods)], max.pre.periods)
        }
        
        if(length(pre.periods) != length(x$pre.periods)){
            change.pre.periods <- 1
        }
        else if(all(pre.periods == x$pre.periods)){
            change.pre.periods <- 0
        }
        else{
            change.pre.periods <- 1
        }

        ## bound
        data2 <- NULL
        if (bound != "none") {
            if (is.null(x$est.att)) {
                message("No uncertainty estimates.\n")
                bound <- "none"
            }
        }
        if (bound != "none"|| "equiv.p" %in% stats) {
            time0 <- NULL
            if (switch.on == TRUE) {
                if (sum(time[show] <= 0) == 0) {
                    message("No pretreatment periods are to be plotted.\n")
                    time0 <- 1:length(time[show])
                } else {
                    time0 <- which(time[show] <= 0)
                }
                att.sub <- as.matrix(est.att[show, c("CI.lower.90", "CI.upper.90")])
                minBound <- max(abs(att.sub[time0, c("CI.lower.90", "CI.upper.90")]), na.rm = TRUE)
            } 
            else {
                if (sum(time[show] > 0) == 0) {
                    message("No non-treatment periods are to be plotted.\n")
                    time0 <- 1:length(time[show])
                } else {
                    time0 <- which(time[show] >= 1)
                }
                att.sub <- as.matrix(x$att.off.bound[show, ])
                minBound <- max(abs(att.sub[time0, c("CI.lower", "CI.upper")]),na.rm = TRUE)
            }
            # Min range
            minbound <- c(-minBound, minBound)
            equiv.range <- c(-1 * tost.threshold, tost.threshold)

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
                if (plot.ci == "90") {
                    set.labels <- "Residual Average (w/ 90% CI)" 
                } 
                else if(plot.ci == "95") {
                    set.labels <- "ATT (w/ 95% CI)" 
                }
                else{
                    set.labels <- "ATT"
                }                      
            } 
            else {
                set.labels <- legend.labs
            }
            set.colors <- "#000000FF"
            set.linetypes <- "solid"
            set.size <- 1

            ## create a dataset for bound

            if (bound.old == "equiv") {
                use2 <- c(rep(equiv.range, each = length(bound.time)))
                data2 <- cbind.data.frame(bound = use2)
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
        } 
        else if (switch.on == FALSE) {
            if (is.null(x$est.att.off)==TRUE) {
                CI <- FALSE
            } else {
                CI <- TRUE
            }           
        }

        if(plot.ci=="none"){
            CI <- FALSE
        }

        ## data frame for main estimates
        if (switch.on == TRUE) {            
            ## switch-on effect
            if (CI == FALSE) {             
                data <- cbind.data.frame(time, ATT = x$att, count = count.num)[show,]                
            } 
            else {
                tb <- est.att
                data <- cbind.data.frame(time, tb)[show,]
                colnames(data)[2] <- "ATT"
                if (plot.ci %in% c("90", "95")) {
                    if (placeboTest == TRUE && length(placebo.period) != 1) {
                        data[,"ATT4"] <- data[,"ATT3"] <- data[,"ATT2"] <- data[,"ATT"]
                        ci.name <- NULL
                        if (plot.ci == "95") {
                            ci.name <- c("CI.lower", "CI.upper")
                        } else if (plot.ci == "90") {
                            ci.name <- c("CI.lower.90", "CI.upper.90")
                        }
                        data[,"CI.lower4"] <- data[,"CI.lower3"] <- data[,"CI.lower2"] <- data[,ci.name[1]]
                        data[,"CI.upper4"] <- data[,"CI.upper3"] <- data[,"CI.upper2"] <- data[,ci.name[2]]                    
                        pos1 <- intersect(which(data[,"time"] >= (placebo.period[1] + 1)), which(data[,"time"] <= (placebo.period[2] - 1)))
                        pos2 <- c(which(data[,"time"] <= (placebo.period[1] - 1)), which(data[,"time"] >= (placebo.period[2] + 1)))
                        pos3 <- intersect(which(data[,"time"] >= (placebo.period[1])), which(data[,"time"] <= (placebo.period[2] - 1)))
                        pos4 <- c(which(data[,"time"] <= (placebo.period[1] - 2)), which(data[,"time"] >= (placebo.period[2] + 1)))
                        data[pos1, c("ATT", ci.name[1], ci.name[2])] <- NA
                        data[pos2, c("ATT2","CI.lower2","CI.upper2")] <- NA
                        data[pos3, c("ATT3","CI.lower3","CI.upper3")] <- NA
                        data[pos4, c("ATT4","CI.lower4","CI.upper4")] <- NA
                    }
                }
            }
        }
        else{
            ## exit treatment plot
            if (CI == FALSE) {               
                data <- cbind.data.frame(time, ATT = x$att.off, count = count.num)[show,]                
            }
            else {
                tb <- est.att.off
                data <- cbind.data.frame(time, tb)[show,]
                colnames(data)[2] <- "ATT"
                colnames(data)[7] <- 'count'
                if (plot.ci %in% c("90", "95")) {
                    if (carryoverTest == TRUE && length(carryover.period) != 1) {
                        data[,"ATT6"] <- data[,"ATT5"] <- data[,"ATT4"] <- data[,"ATT3"] <- data[,"ATT2"] <- data[,"ATT"]
                        ci.name <- NULL
                        if (plot.ci == "95") {
                            ci.name <- c("CI.lower", "CI.upper")
                        } else if (plot.ci == "90") {
                            ci.name <- c("CI.lower.90", "CI.upper.90")
                        }
                        if (is.null(x$est.carry.att)) {
                            data[,"CI.lower4"] <- data[,"CI.lower3"] <- data[,"CI.lower2"] <- data[,ci.name[1]]
                            data[,"CI.upper4"] <- data[,"CI.upper3"] <- data[,"CI.upper2"] <- data[,ci.name[2]]                    
                            pos1 <- intersect(which(data[,"time"] >= (carryover.period[1] + 1)), which(data[,"time"] <= (carryover.period[2] - 1)))
                            pos2 <- c(which(data[,"time"] <= (carryover.period[1] - 1)), which(data[,"time"] >= (carryover.period[2] + 1)))
                            pos3 <- intersect(which(data[,"time"] >= (carryover.period[1])), which(data[,"time"] <= (carryover.period[2] - 1)))
                            pos4 <- c(which(data[,"time"] <= (carryover.period[1] - 2)), which(data[,"time"] >= (carryover.period[2] + 1)))
                            data[pos1, c("ATT", ci.name[1], ci.name[2])] <- NA
                            data[pos2, c("ATT2","CI.lower2","CI.upper2")] <- NA
                            data[pos3, c("ATT3","CI.lower3","CI.upper3")] <- NA
                            data[pos4, c("ATT4","CI.lower4","CI.upper4")] <- NA
                        } 
                        else {
                            data[,"CI.lower6"] <- data[,"CI.lower5"] <- data[,"CI.lower4"] <- data[,"CI.lower3"] <- data[,"CI.lower2"] <- data[,ci.name[1]]
                            data[,"CI.upper6"] <- data[,"CI.upper5"] <- data[,"CI.upper4"] <- data[,"CI.upper3"] <- data[,"CI.upper2"] <- data[,ci.name[2]]                    
                            
                            pos1 <- intersect(which(data[,"time"] >= (carryover.period[1])), which(data[,"time"] <= (carryover.period[2])))
                            pos2 <- c(which(data[,"time"] <= min(-dim(x$est.carry.att)[1] -1 , -2)), which(data[,"time"] >= (carryover.period[2]) + 1))
                            pos5 <- intersect(which(data[,"time"] >= min(-dim(x$est.carry.att)[1], -1)), which(data[,"time"] <= -1))
                            
                            pos3 <- intersect(which(data[,"time"] >= (carryover.period[1])), which(data[,"time"] <= (carryover.period[2] - 1)))
                            pos4 <- c(which(data[,"time"] <= min(-dim(x$est.carry.att)[1]-1, -3)), which(data[,"time"] >= (carryover.period[2] + 1)))
                            pos6 <- intersect(which(data[,"time"] >= min(-dim(x$est.carry.att)[1] + 1, -1)), which(data[,"time"] <=  -1))

                            data[unique(c(pos1, pos5)), c("ATT", ci.name[1], ci.name[2])] <- NA
                            data[unique(c(pos2, pos5)), c("ATT2","CI.lower2","CI.upper2")] <- NA
                            data[unique(c(pos3, pos6)), c("ATT3","CI.lower3","CI.upper3")] <- NA
                            data[unique(c(pos4, pos6)), c("ATT4","CI.lower4","CI.upper4")] <- NA
                            data[unique(c(pos1, pos2)), c("ATT5","CI.lower5","CI.upper5")] <- NA
                            data[unique(c(pos3, pos4)), c("ATT6","CI.lower6","CI.upper6")] <- NA
                            
                        }
                    }
                }
            } 
        }


        # height of the histogram
        if (CI == FALSE) {
                message("Uncertainty estimates not available.\n")
                if (length(ylim) != 0) {
                    rect.length <- (ylim[2] - ylim[1]) / 5
                    rect.min <- ylim[1]
                } else {
                    rect.length <- (max(data[,"ATT"], na.rm = TRUE) - min(data[,"ATT"], na.rm = TRUE))/2
                    rect.min <- min(data[,"ATT"], na.rm = TRUE) - rect.length
                } 
        } 
        else {
                if (length(ylim) != 0) {
                    rect.length <- (ylim[2] - ylim[1]) / 5
                    rect.min <- ylim[1]
                } else {
                    rect.length <- (max(data[,"CI.upper"], na.rm = TRUE) - min(data[,"CI.lower"], na.rm = TRUE))/2
                    rect.min <- min(data[,"CI.lower"], na.rm = TRUE) - rect.length 
                }  
        }

        ## plotting
        ## line

        if(start0 == TRUE){
            data$time <- data$time - 1
            data2$time <- data2$time - 1
            if(!is.null(placebo.period)){
                placebo.period <- placebo.period - 1
            }
            if(!is.null(carryover.period)){
                carryover.period <- carryover.period - 1
            }
        }

        
            # plot bound
            if (bound.old == "none") {
                p <- ggplot(data)             
            } 
            else { ## with bounds
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
        p <- p + xlab(xlab) +  ylab(ylab) 

        ## theme
        if (theme.bw == TRUE) {
            p <- p + theme_bw() 
        }

        ## grid
        if (gridOff == TRUE) {
            p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }

        # horizontal 0 line
        p <- p + geom_hline(yintercept = 0, colour = lcolor,size = lwidth)
        # vertical 0 line

            if (plot.ci == "both") {
                lwidth <- lwidth * 0.5
            }
        
            if (length(xlim)!=0) {
                if ((xlim[2]>=1 & switch.on == TRUE) | (xlim[1]<=0 & switch.on == FALSE)) {
                    if(start0 == FALSE){
                        if(plot.ci == 'both'){                                
                                p <- p + geom_vline(xintercept = 0.5, colour=lcolor,size = lwidth)
                            }
                            else{                                
                                p <- p + geom_vline(xintercept = 0, colour=lcolor,size = lwidth)                             
                            }
                        }
                        else{
                            if(plot.ci == 'both'){                                
                                p <- p + geom_vline(xintercept = -0.5, colour=lcolor,size = lwidth)
                            }
                            else{
                                p <- p + geom_vline(xintercept = -1, colour=lcolor,size = lwidth)                             
                            }
                        }
                }
            } 
            else {
                    if(start0 == FALSE){
                        if(plot.ci == 'both'){
                            p <- p + geom_vline(xintercept = 0.5, colour=lcolor,size = lwidth)
                        }
                        else{
                            p <- p + geom_vline(xintercept = 0, colour=lcolor,size = lwidth)                             
                        }
                    }
                    else{
                        if(plot.ci == 'both'){
                            p <- p + geom_vline(xintercept = -0.5, colour=lcolor,size = lwidth)
                        }
                        else{
                            p <- p + geom_vline(xintercept = -1, colour=lcolor,size = lwidth)                             
                        }
                    }
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

        ## add ATT point estimates
        classic <- 0 
        if(highlight==FALSE){
            classic <- 1
        }
        if(placeboTest==TRUE && length(placebo.period) == 1 && plot.ci %in% c("90", "95")){
            classic <- 1
        }
        if(carryoverTest==TRUE && length(carryover.period) == 1 && plot.ci %in% c("90", "95")){
            classic <- 1
        }
        if(carryoverTest == FALSE && placeboTest==FALSE){
            classic <- 1
        }



        if (classic==1) {
            ## point estimates
            if (plot.ci %in% c("90", "95","none")) {
                p <- p + geom_line(data = data, aes(x=time, y=ATT), size = 1.2)
                if (vis == "connected") {
                    p <- p + geom_point(data = data, aes(x=time, y=ATT), size = 1.2, na.rm=TRUE)
                }
                ## CIs
                if (CI == TRUE) {
                    if (plot.ci == "95") {
                        p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower, ymax=CI.upper),alpha=0.2)                
                    } 
                    else {
                        p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower.90, ymax=CI.upper.90),alpha=0.2)                
                    } 
                }
            } 
            else if(plot.ci == 'both') {
                if(CI==TRUE){
                    if(plot.ci.point %in% c("both","95")){
                        p <- p + geom_pointrange(data = data, aes(x = time, y = ATT, ymin=CI.lower, ymax=CI.upper), lwd=0.6,fatten = 2)
                    }
                    if(plot.ci.point %in% c("both","90")){
                        p <- p + geom_pointrange(data = data, aes(x = time, y = ATT, ymin=CI.lower.90, ymax=CI.upper.90), lwd=0.6,fatten = 2)
                    }                    
                }
                else{
                    p <- p + geom_point(data = data, aes(x = time, y = ATT), size=1.2)                    
                }

            }
        } 
        else if(classic==0 && switch.on==TRUE) {
            ## point estimates
            ## placebo tests
            p.label <- NULL
            if (plot.ci %in% c("90", "95","none")) {
                if(vis=='none'){
                    p <- p + geom_line(data = data, aes(time, ATT3), size = 1.2)
                    p <- p + geom_line(data = data, aes(time, ATT4), size = 1.2, colour = "#4671D5")
                }
                if (vis == "connected") {
                    p <- p + geom_line(data = data, aes(time, ATT3), size = 0.7)
                    p <- p + geom_line(data = data, aes(time, ATT4), size = 0.7, colour = "#4671D5")
                    p <- p + geom_point(data = data, aes(time, ATT), size = 1.2, na.rm=TRUE)
                    p <- p + geom_point(data = data, aes(time, ATT2), size = 1.2, colour = "#4671D5", na.rm = TRUE)
                }
                ## CIs
                p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower3, ymax=CI.upper3),alpha=0.2, na.rm = FALSE)
                p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower4, ymax=CI.upper4),alpha=0.2, fill = "#0000FF", na.rm = FALSE)
            } 
            else if(plot.ci == "both") {
                pos.ci <- intersect(which(data[,"time"] >= (placebo.period[1])), which(data[,"time"] <= (placebo.period[length(placebo.period)])))
                pos.ci2 <- setdiff(1:dim(data)[1], pos.ci)
                if(CI==TRUE){
                    if(plot.ci.point %in% c("both","95")){
                        p <- p + geom_pointrange(data = data[pos.ci,], aes(x = time, y = ATT, ymin=CI.lower, ymax=CI.upper), lwd=0.6, color="blue", fill="blue",fatten = 2)
                        p <- p + geom_pointrange(data = data[pos.ci2,], aes(x = time, y = ATT, ymin=CI.lower, ymax=CI.upper), lwd=0.6,fatten = 2)
                    }
                    if(plot.ci.point %in% c("both","90")){
                        p <- p + geom_pointrange(data = data[pos.ci,], aes(x = time, y = ATT, ymin=CI.lower.90, ymax=CI.upper.90), lwd=0.6, color="blue", fill="blue",fatten = 2)
                        p <- p + geom_pointrange(data = data[pos.ci2,], aes(x = time, y = ATT, ymin=CI.lower.90, ymax=CI.upper.90), lwd=0.6,fatten = 2)
                    }                    
                }
                else{
                    p <- p + geom_point(data = data[pos.ci,], aes(x = time, y = ATT), lwd=0.6, color="blue", fill="blue",size=1.2)
                    p <- p + geom_point(data = data[pos.ci2,], aes(x = time, y = ATT), lwd=0.6,size=1.2)                    
                }
            }
        }
        else if(classic==0 && switch.on==FALSE){
            ## point estimates
            ## carryover tests
            p.label <- NULL
            if (plot.ci %in% c("90", "95","none")) {
                if (plot.ci %in% c("90", "95","none")) {
                    if (is.null(x$est.carry.att)) {
                        if(vis == "none"){
                            p <- p + geom_line(data = data, aes(time, ATT3), size = 1.2)
                            p <- p + geom_line(data = data, aes(time, ATT4), size = 1.2, colour = "red")
                        }
                        if (vis == "connected") {
                            p <- p + geom_line(data = data, aes(time, ATT3), size = 0.7)
                            p <- p + geom_line(data = data, aes(time, ATT4), size = 0.7, colour = "red")
                            p <- p + geom_point(data = data, aes(time, ATT), size = 1.2 , na.rm=TRUE)
                            p <- p + geom_point(data = data, aes(time, ATT2), size = 1.2, colour = "red", na.rm=TRUE)
                        }
                        ## CIs
                        p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower3, ymax=CI.upper3),alpha=0.2,na.rm = FALSE)
                        p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower4, ymax=CI.upper4),alpha=0.2, fill = "pink",na.rm = FALSE)
                    } 
                    else {
                        data[,'time'] <- length(x$carry.att) + data[,'time']
                        if(vis == "none"){
                            p <- p + geom_line(data = data, aes(time, ATT3), size = 1.2)
                            p <- p + geom_line(data = data, aes(time, ATT4), size = 1.2, colour = "red")
                            p <- p + geom_line(data = data, aes(time, ATT6), size = 1.2, colour = "blue")
                        }
                        if (vis == "connected") {
                            p <- p + geom_line(data = data, aes(time, ATT3), size = 0.7)
                            p <- p + geom_line(data = data, aes(time, ATT4), size = 0.7, colour = "red")
                            p <- p + geom_line(data = data, aes(time, ATT6), size = 0.7, colour = "blue")
                            p <- p + geom_point(data = data, aes(time, ATT), size = 1.2 , na.rm=TRUE)
                            p <- p + geom_point(data = data, aes(time, ATT2), size = 1.2, colour = "red", na.rm=TRUE)
                            p <- p + geom_point(data = data, aes(time, ATT5), size = 1.2, colour = "blue", na.rm=TRUE)
                            p <- p + geom_point(aes(x=0, y=data[which(data$time==0),'ATT5']), size = 1.2 , na.rm=TRUE)
                        }
                        ## CIs
                        p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower3, ymax=CI.upper3),alpha=0.2,na.rm = FALSE)
                        p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower4, ymax=CI.upper4),alpha=0.2, fill = "pink",na.rm = FALSE)
                        p <- p + geom_ribbon(data = data, aes(x = time, ymin=CI.lower6, ymax=CI.upper6),alpha=0.2, fill = "blue",na.rm = FALSE)
                    }
                } 
            } 
            else if(plot.ci == "both") {
                if (is.null(x$est.carry.att)) {
                    pos.ci <- intersect(which(data[,"time"] >= (carryover.period[1])), which(data[,"time"] <= (carryover.period[length(carryover.period)])))
                    pos.ci2 <- setdiff(1:dim(data)[1], pos.ci)
                    if(CI==TRUE){
                        if(plot.ci.point %in% c("both","95")){
                            p <- p + geom_pointrange(data = data[pos.ci,], aes(x = time, y = ATT, ymin=CI.lower, ymax=CI.upper), lwd=0.6, color="red", fill="red",fatten = 2)
                            p <- p + geom_pointrange(data = data[pos.ci2,], aes(x = time, y = ATT, ymin=CI.lower, ymax=CI.upper), lwd=0.6,fatten = 2)
                        }
                        if(plot.ci.point %in% c("both","90")){
                            p <- p + geom_pointrange(data = data[pos.ci,], aes(x = time, y = ATT, ymin=CI.lower.90, ymax=CI.upper.90), lwd=0.6, color="red", fill="red",fatten = 2)
                            p <- p + geom_pointrange(data = data[pos.ci2,], aes(x = time, y = ATT, ymin=CI.lower.90, ymax=CI.upper.90), lwd=0.6,fatten = 2)
                        }
                    }
                    else{
                        p <- p + geom_point(data = data[pos.ci,], aes(x = time, y = ATT), lwd=0.6, color="red", fill="red",size=1.2)
                        p <- p + geom_point(data = data[pos.ci2,], aes(x = time, y = ATT), lwd=0.6,size=1.2)                    
                    }                    
                }
                else{
                    T.carry <- -dim(x$est.carry.att)[1] + 1
                    pos.ci <- intersect(which(data[,"time"] >= (carryover.period[1])), which(data[,"time"] <= (carryover.period[length(carryover.period)])))
                    pos.ci2 <- intersect(which(data[,"time"] >= T.carry), which(data[,"time"] <= 0))
                    pos.ci3 <- setdiff(1:dim(data)[1], c(pos.ci, pos.ci2))
                    data[,'time'] <- length(x$carry.att) + data[,'time']
                    if(CI==TRUE){
                        if(plot.ci.point %in% c("both","95")){
                            p <- p + geom_pointrange(data = data[pos.ci,], aes(x = time, y = ATT, ymin=CI.lower, ymax=CI.upper), lwd=0.6, color="red", fill="red",fatten = 2)
                            p <- p + geom_pointrange(data = data[pos.ci2,], aes(x = time, y = ATT, ymin=CI.lower, ymax=CI.upper), lwd=0.6, color="blue", fill="blue",fatten = 2)
                            p <- p + geom_pointrange(data = data[pos.ci3,], aes(x = time, y = ATT, ymin=CI.lower, ymax=CI.upper), lwd=0.6 ,fatten = 2)
                        }
                        if(plot.ci.point %in% c("both","90")){
                            p <- p + geom_pointrange(data = data[pos.ci,], aes(x = time, y = ATT, ymin=CI.lower.90, ymax=CI.upper.90), lwd=0.6, color="red", fill="red",fatten = 2)
                            p <- p + geom_pointrange(data = data[pos.ci2,], aes(x = time, y = ATT, ymin=CI.lower.90, ymax=CI.upper.90), lwd=0.6, color="blue", fill="blue",fatten = 2)
                            p <- p + geom_pointrange(data = data[pos.ci3,], aes(x = time, y = ATT, ymin=CI.lower.90, ymax=CI.upper.90), lwd=0.6,fatten = 2)
                        }
                    }
                    else{
                        p <- p + geom_point(data = data[pos.ci,], aes(x = time, y = ATT), lwd=0.6, color="red", fill="red",size=1.2)
                        p <- p + geom_point(data = data[pos.ci2,], aes(x = time, y = ATT), lwd=0.6, color="blue", fill="red",size=1.2)
                        p <- p + geom_point(data = data[pos.ci3,], aes(x = time, y = ATT), lwd=0.6, size=1.2)                       
                    }
  
                }

            }
        }

        ## print stats
        p.label <- NULL

        if(x$loo==TRUE && loo==TRUE){
            #recalculate p value and f value
            loo.equiv <- 1
        }
        else{
            loo.equiv <- 0
        }
        
        if (type %in% c('equiv','gap') && loo.equiv == 0) { 
            for (i in 1:length(stats)) {
                if ("F.p" %in% stats[i]) {
                    if (change.proportion | change.pre.periods | !is.null(show.group) | use.balance) {
                        x$loo <- FALSE
                        test.out <- diagtest(x, 
                                             proportion = proportion, 
                                             pre.periods = pre.periods, 
                                             f.threshold = f.threshold)
                        f.p <- test.out$f.p
                    } 
                    else {
                        f.p <- x$test.out$f.p
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",f.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("F.stat" %in% stats[i]) {
                    if (change.proportion | change.pre.periods | !is.null(show.group) | use.balance) {
                        x$loo <- FALSE
                        test.out <- diagtest(x, 
                                             proportion = proportion, 
                                             pre.periods = pre.periods, 
                                             f.threshold = f.threshold)
                        f.stat <- test.out$f.stat
                    } 
                    else {
                        f.stat <- x$test.out$f.stat
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",f.stat))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("F.equiv.p" %in% stats[i]) {
                    # calculate new p value (ziyi: re-add this)
                    if (change.f.threshold | change.proportion | change.pre.periods | !is.null(show.group)|use.balance) {
                        x$loo <- FALSE
                        # some problems here, should change to change.f.threshold; change.proportion; change.pre.periods
                        test.out <- diagtest(x, 
                                             proportion = proportion, 
                                             pre.periods = pre.periods, 
                                             f.threshold = f.threshold)
                        f.equiv.p <- test.out$f.equiv.p
                    } 
                    else {
                        f.equiv.p <- x$test.out$f.equiv.p
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",f.equiv.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("equiv.p" %in% stats[i] && placeboTest == 0) {
                    # calculate new p value (ziyi: re-add this)
                    if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)|use.balance) {
                        x$loo <- FALSE
                        test.out <- diagtest(x, 
                                             proportion = proportion, 
                                             pre.periods = pre.periods, 
                                             tost.threshold = tost.threshold)
                        tost.equiv.p <- test.out$tost.equiv.p
                    } 
                    else {
                        tost.equiv.p <- x$test.out$tost.equiv.p
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",tost.equiv.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("placebo.p" %in% stats[i]) {
                    if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
                        test.out <- diagtest(x, 
                                             proportion = proportion, 
                                             pre.periods = pre.periods, 
                                             tost.threshold = tost.threshold)
                        placebo.p <- test.out$placebo.p
                    } 
                    else {
                        placebo.p <- x$test.out$placebo.p
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",placebo.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("equiv.p" %in% stats[i] && placeboTest==1) {
                    p.label1 <- NULL
                    # calculate new p value (ziyi: re-add this)
                    if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
                        test.out <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
                        placebo.equiv.p <- test.out$placebo.equiv.p
                    } 
                    else {
                        placebo.equiv.p <- x$test.out$placebo.equiv.p
                    }
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f", placebo.equiv.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                } 
            } 
        }
        else if(type %in% c('equiv','gap') && loo.equiv == 1){ #loo
            for (i in 1:length(stats)) {
                if ("F.p" %in% stats[i]) {
                    if (change.proportion | change.pre.periods | !is.null(show.group)) {
                        test.out <- diagtest(x, 
                                             proportion = proportion, 
                                             pre.periods = pre.periods, 
                                             f.threshold = f.threshold)
                        f.p <- test.out$f.p
                    } 
                    else {
                        f.p <- x$test.out$f.p
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",f.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("F.stat" %in% stats[i]) {
                    if (change.proportion | change.pre.periods | !is.null(show.group)) {
                        test.out <- diagtest(x, 
                                             proportion = proportion, 
                                             pre.periods = pre.periods, 
                                             f.threshold = f.threshold)
                        f.stat <- test.out$f.stat
                    } 
                    else {
                        f.stat <- x$test.out$f.stat
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",f.stat))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("F.equiv.p" %in% stats[i]) {
                    # calculate new p value (ziyi: re-add this)
                    if (change.f.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
                        loo.test.out <- diagtest(x, 
                                             proportion = proportion, 
                                             pre.periods = pre.periods, 
                                             f.threshold = f.threshold)
                        f.equiv.p <- loo.test.out$f.equiv.p
                    } 
                    else {
                        f.equiv.p <- x$loo.test.out$f.equiv.p
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",f.equiv.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("equiv.p" %in% stats[i]) {
                    # calculate new p value (ziyi: re-add this)
                    if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
                        loo.test.out <- diagtest(x, 
                                             proportion = proportion, 
                                             pre.periods = pre.periods, 
                                             tost.threshold = tost.threshold)
                        tost.equiv.p <- loo.test.out$tost.equiv.p
                    } 
                    else {
                        tost.equiv.p <- x$loo.test.out$tost.equiv.p
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",tost.equiv.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
            }
        } 
        else if(type=='gap' && placeboTest==TRUE){               
            ## stats
            for (i in 1:length(stats)) {
                if ("placebo.p" %in% stats[i]) {
                    if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
                        test.out <- diagtest(x, 
                                             proportion = proportion, 
                                             pre.periods = pre.periods, 
                                             tost.threshold = tost.threshold)
                        placebo.p <- test.out$placebo.p
                    } 
                    else {
                        placebo.p <- x$test.out$placebo.p
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",placebo.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("equiv.p" %in% stats[i]) {
                    p.label1 <- NULL
                    # calculate new p value (ziyi: re-add this)
                    if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
                        test.out <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
                        placebo.equiv.p <- test.out$placebo.equiv.p
                    } 
                    else {
                        placebo.equiv.p <- x$test.out$placebo.equiv.p
                    }
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f", placebo.equiv.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }                
            } 
        }
        else if(type=='exit' && carryoverTest==TRUE){
            ## stats
            for (i in 1:length(stats)) {
                if ("carryover.p" %in% stats[i]) {
                    if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
                        test.out <- diagtest(x, 
                                             proportion = proportion, 
                                             pre.periods = pre.periods, 
                                             tost.threshold = tost.threshold)
                        carryover.p <- test.out$carryover.p
                    } 
                    else {
                        carryover.p <- x$test.out$carryover.p
                    }
                    p.label1 <- NULL
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f",carryover.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }
                if ("equiv.p" %in% stats[i]) {
                    p.label1 <- NULL
                    # calculate new p value (ziyi: re-add this)
                    if (change.tost.threshold | change.proportion | change.pre.periods| !is.null(show.group)) {
                        test.out <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
                        carryover.equiv.p <- test.out$carryover.equiv.p
                    } 
                    else {
                        carryover.equiv.p <- x$test.out$carryover.equiv.p
                    }
                    p.label1 <- paste0(stats.labs[i],": ", sprintf("%.3f", carryover.equiv.p))
                    p.label <- paste0(p.label, p.label1, "\n")
                }                
            } 
        }
        



        ## mark stats
        ##hpos <- ifelse(switch.on == TRUE, 0, 1)
        hpos <- 0
        if ("none" %in% stats == FALSE) {
            if (is.null(stats.pos)) {                
                if (switch.on == TRUE) {
                    stats.pos[1] <- min(data[,"time"], na.rm = 1)
                } 
                else {
                    stats.pos[1] <- min(data[,"time"], na.rm = 1)
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

        if(dim(data)[1]>4){
            p <- p + scale_x_continuous(labels=scaleFUN)
        }
        else{
            p <- p + scale_x_continuous(breaks=c(data[,'time']))
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

        ##xlim

        
    }

    if(type == "calendar"){

        CI <- NULL
        if (is.null(x$est.eff.calendar)==TRUE) {
            CI <- FALSE
        } 
        else {
            CI <- TRUE
        }
        if(plot.ci=="none"){
            CI <- FALSE
        }
        ## axes labels
        if (is.null(xlab) == TRUE) {
            xlab <- "Calendar Time"         
        } 
        else if (xlab == "") {
            xlab <- NULL
        }

        if (is.null(ylab) == TRUE) {
            ylab <- ytitle
        } 
        else if (ylab == "") {
            ylab <- NULL
        }

        ## y=0 line type
        lcolor <- "white"
        lwidth <- 2
        if (theme.bw == TRUE) {
            lcolor <- "#AAAAAA70"
            lwidth <- 1.5
        }

        if (CI == FALSE) {
            message("Uncertainty estimates not available.\n")
            data.1 <- x$eff.calendar
                data.2 <- x$eff.calendar.fit
                if (length(ylim) != 0) {
                    rect.length <- (ylim[2] - ylim[1]) / 5
                    rect.min <- ylim[1]
                } else {
                    rect.length <- (max(c(data.1,data.2), na.rm = TRUE) - min(c(data.1,data.2), na.rm = TRUE))/2
                    rect.min <- min(c(data.1,data.2), na.rm = TRUE) - rect.length
                }
                d1 <- data.1 <- as.matrix(x$eff.calendar[which(!is.na(x$eff.calendar[,1])),])
                d2 <- data.2 <- as.matrix(x$eff.calendar.fit[which(!is.na(x$eff.calendar.fit[,1])),])
                if(dim(d1)[2]==1){
                    d1 <- data.1 <- t(d1)
                    rownames(d1) <- rownames(data.1) <- rownames(x$eff.calendar)[which(!is.na(x$est.eff.calendar[,1]))]
                }
                if(dim(d2)[2]==1){
                    d2 <- data.2 <- t(d2)
                    rownames(d2) <- rownames(data.2) <- rownames(x$eff.calendar.fit)[which(!is.na(x$est.eff.calendar.fit[,1]))]
                }
        } 
        else {
                if(is.null(x$est.eff.calendar)){
                    stop("Uncertainty estimates not available.\n")
                }
                d1 <- data.1 <- as.matrix(x$est.eff.calendar[which(!is.na(x$est.eff.calendar[,1])),])
                d2 <- data.2 <- as.matrix(x$est.eff.calendar.fit[which(!is.na(x$est.eff.calendar.fit[,1])),])
                if(dim(d1)[2]==1){
                    d1 <- data.1 <- t(d1)
                    rownames(d1) <- rownames(data.1) <- rownames(x$est.eff.calendar)[which(!is.na(x$est.eff.calendar[,1]))]
                }
                if(dim(d2)[2]==1){
                    d2 <- data.2 <- t(d2)
                    rownames(d2) <- rownames(data.2) <- rownames(x$est.eff.calendar.fit)[which(!is.na(x$est.eff.calendar.fit[,1]))]
                }
                
                if (length(ylim) != 0) {
                    rect.length <- (ylim[2] - ylim[1]) / 5
                    rect.min <- ylim[1]
                } else {
                    rect.length <- (max(c(data.1[,4],data.2[,4]), na.rm = TRUE) - min(c(data.1[,3],data.2[,3]), na.rm = TRUE))/2
                    rect.min <- min(c(data.1[,3],data.2[,3]), na.rm = TRUE) - rect.length 
                }  
        }
        p <- ggplot()
        ## xlab and ylab 
        p <- p + xlab(xlab) +  ylab(ylab)

        ## theme
        if (theme.bw == TRUE) {
            p <- p + theme_bw() 
        }

        ## grid
        if (gridOff == TRUE) {
            p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }

        # horizontal 0 line
        p <- p + geom_hline(yintercept = 0, colour = lcolor,size = lwidth)

        TTT <- as.numeric(rownames(data.1))
        TTT.2 <- as.numeric(rownames(data.2))

        if(CI==FALSE){
            p <- p + geom_point(aes(x=TTT,y=d1[,1]),color='gray50',fill='gray50',alpha=1,size=1.2)
            p <- p + geom_line(aes(x=TTT.2,y=d2[,1]),color='skyblue',size=1.1)
        }
        else{
            p <- p + geom_line(aes(x=TTT.2,y=d2[,1]),color='skyblue',size=1.1)
            p <- p + geom_ribbon(aes(x=TTT.2,ymin=d2[,3],ymax=d2[,4]),color='skyblue',fill='skyblue',alpha=0.5,size=0)
            p <- p + geom_pointrange(aes(x=TTT,y=d1[,1],ymin=d1[,3],ymax=d1[,4]),color='gray50',fill='gray50',alpha=1,size=0.6)
        }

        if(count==TRUE){
                T.start <- c()
                T.end <- c()
                ymin <- c()
                ymax <- c()
                T.gap <- (max(TTT)-min(TTT))/length(TTT)
                for(i in c(1:dim(d1)[1])){
                    T.start <- c(T.start,TTT[i]-0.25*T.gap)
                    T.end <- c(T.end,TTT[i]+0.25*T.gap)
                    ymin <- c(ymin, rect.min)
                    ymax <- c(ymax, rect.min+rect.length*d1[i,'count']/max(d1[,'count']))
                }
                data.toplot <- cbind.data.frame(xmin=T.start,
                                                xmax=T.end,
                                                ymin=ymin,
                                                ymax=ymax)
                max.count.pos <- mean(TTT[which.max(d1[,'count'])])
                p <- p + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),data=data.toplot,fill='gray50',alpha=0.3,size=0.3,color='black')
                p <- p + annotate("text", x = max.count.pos - 0.02 * T.gap, 
                    y = max(data.toplot$ymax) + 0.2 * rect.length, 
                    label = max(x$N.calendar), size = cex.text * 0.8, hjust = 0.5)                  
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

        if(length(TTT)<=4){
            p <- p + scale_x_continuous(breaks=TTT)
        }
        else{
            p <- p + scale_x_continuous(labels=scaleFUN)
        }

        # ## xlim
        # if(is.null(xlim)){
        #     if(is.na(d1[1,1])){
        #         ## drop all periods before first non-missing
        #         for(j in c(2:dim(d1)[1])){
        #             if(!is.na(d1[j,1])){
        #                 xlim <- c(TTT[j],max(TTT))
        #                 break
        #             }
        #         }
        #     }            
        # }
        ## xlim
        if (is.null(xlim) == FALSE) {
            p <- p + coord_cartesian(xlim = xlim)
        }

        p <- p + geom_hline(yintercept = x$att.avg,color='red',size=0.8,linetype='dashed')

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

    }

    if(type == "box"){
        if (is.null(xlab)==TRUE) {
            xlab <- index[2]
        } else if (xlab == "") {
            xlab <- NULL
        }
        if (is.null(ylab)==TRUE) {
            ylab <- "Estimated Treatment Effects"
        } else if (ylab == "") {
            ylab <- NULL
        }
        if (is.null(main)==TRUE) {
            main <- "Treatment Effects"
        } else if (main == "") {
            main <- NULL
        }

        ## y=0 line type
        lcolor <- "white"
        lwidth <- 2
        if (theme.bw == TRUE) {
            lcolor <- "#AAAAAA70"
            lwidth <- 1.5
        }

        p <- ggplot()
        ## xlab and ylab 
        p <- p + xlab(xlab) +  ylab(ylab) 

        ## theme
        if (theme.bw == TRUE) {
            p <- p + theme_bw() 
        }

        ## title
        if (is.null(main) == TRUE) {
            p <- p + ggtitle(maintext)
        } 
        else if (main!=""){
            p <- p + ggtitle(main)
        }

        ## grid
        if (gridOff == TRUE) {
            p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }

        # horizontal 0 line
        p <- p + geom_hline(yintercept = 0, colour = lcolor,size = lwidth)

        complete.index.eff <- which(!is.na(x$eff))
        complete.index.time <- which(!is.na(x$T.on))
        complete.index <- intersect(complete.index.eff,complete.index.time)
        eff.use <- x$eff[complete.index]
        time.use <- x$T.on[complete.index]
        id.mat <- rep(colnames(x$eff),each=dim(x$eff)[1])
        id.use <- id.mat[complete.index]
        data.toplot <- cbind.data.frame(time=time.use,id=id.use,eff=eff.use)
        data.count <- cbind.data.frame(time=x$time,count=x$count)

        if(start0==TRUE){
            data.toplot$time <- data.toplot$time - 1
            data.count$time <- data.count$time - 1
        }

        if(!is.null(xlim)){
            data.count <- data.count[which(data.count[,'time']>=min(xlim) & data.count[,'time']<=max(xlim)),]
            data.toplot <- data.toplot[which(data.toplot[,'time']>=min(xlim) & data.toplot[,'time']<=max(xlim)),]
        }
        

        data.use <- merge(data.toplot,data.count,by = "time")
        #print(data.use)

        if (length(ylim) != 0) {
            rect.length <- (ylim[2] - ylim[1]) / 5
            rect.min <- ylim[1]
        } 
        else {
            rect.length <- (max(data.use[,'eff'], na.rm = TRUE) - min(data.use[,'eff'], na.rm = TRUE))/3
            rect.min <- min(data.use[,'eff'], na.rm = TRUE) - rect.length
        }

        if(start0==FALSE){
            data.pre.1 <- data.use[which(data.use$time<=0 & data.use$count>=10),]
            data.pre.2 <- data.use[which(data.use$time<=0 & data.use$count<10),]
            data.post.1 <- data.use[which(data.use$time>0 & data.use$count>=10),]
            data.post.2 <- data.use[which(data.use$time>0 & data.use$count<10),]
        }
        else{
            data.pre.1 <- data.use[which(data.use$time<0 & data.use$count>=10),]
            data.pre.2 <- data.use[which(data.use$time<0 & data.use$count<10),]
            data.post.1 <- data.use[which(data.use$time>=0 & data.use$count>=10),]
            data.post.2 <- data.use[which(data.use$time>=0 & data.use$count<10),]
        }


        levels <- as.factor(as.character(data.count[,1]))
        data.pre.1$time <- as.character(data.pre.1$time)
        data.pre.2$time <- as.character(data.pre.2$time)
        data.post.1$time <- as.character(data.post.1$time)
        data.post.2$time <- as.character(data.post.2$time)
        data.pre.1$time <- factor(data.pre.1$time,levels=levels)
        data.pre.2$time <- factor(data.pre.2$time,levels=levels)
        data.post.1$time <- factor(data.post.1$time,levels=levels)
        data.post.2$time <- factor(data.post.2$time,levels=levels)

        p <- p + geom_boxplot(aes(x=time,y=eff),position="dodge", alpha=0.5, 
                        data = data.pre.1,fill='skyblue',
                        outlier.fill='skyblue',outlier.size = 1.25,
                        outlier.color='skyblue',
                        outlier.alpha = 0.5)
        p <- p + geom_boxplot(aes(x=time,y=eff),position="dodge", alpha=0.5, 
                            data = data.post.1,fill='pink',outlier.fill = 'pink',
                            outlier.size = 1.25,outlier.color = 'pink',
                            outlier.alpha = 0.5)

        p <- p + geom_point(aes(x=time,y=eff),data = data.post.2,
                            color="pink", size=1.25, alpha=0.8)
        p <- p + geom_point(aes(x=time,y=eff),data = data.pre.2,
                            color="skyblue", size=1.25, alpha=0.8)
        p <- p + scale_x_discrete(limits =levels)
        
        if(count==TRUE){
                T.start <- c()
                T.end <- c()
                ymin <- c()
                ymax <- c()
                T.gap <- 1
                for(i in c(1:dim(data.count)[1])){
                    T.start <- c(T.start,data.count[i,1]-0.25*T.gap- min(data.count[,1]) + 1) 
                    T.end <- c(T.end,data.count[i,1]+0.25*T.gap- min(data.count[,1])+1) 
                    ymin <- c(ymin, rect.min)
                    ymax <- c(ymax, rect.min+rect.length*data.count[i,2]/max(data.count[,2]))
                }
                data.toplot <- cbind.data.frame(xmin=T.start,
                                                xmax=T.end,
                                                ymin=ymin,
                                                ymax=ymax)
                max.count.pos <- data.count[which.max(data.count[,2]),1][1]- min(data.count[,1])+1
                p <- p + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),data=data.toplot,fill='gray50',alpha=0.3,size=0.3,color='black')
                p <- p + annotate("text", x = max.count.pos - 0.02 * T.gap, 
                    y = max(data.toplot$ymax) + 0.1 * rect.length, 
                    label = max(data.count[,2]), size = cex.text * 0.7, hjust = 0.5)       
        }

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


    }


    if(type == 'status'){
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
            col2 <- c(col2, "4"=NA)
            breaks <- c(breaks,4)
            label <- c(label,"Removed")
        }
        if(5%in%all){
            col2 <- c(col2, "5"=NA)
            breaks <- c(breaks,5)
            if(placeboTest==TRUE){
                col <- c(col,'#66C2A5')
                label <- c(label,"Placebo Tests")
            }
            else if(carryoverTest==TRUE){
                col <- c(col,"#E78AC3")
                label <- c(label,"Carryover Tests")
            }
        }
        if(6%in%all){
            col <- c(col,"#ffc425")
            col2 <- c(col2, "6"=NA)
            breaks <- c(breaks,6)
            label <- c(label,"Carryover Removed")
        }
        if(7%in%all){
            col <- c(col,"#00852B")
            col2 <- c(col2, "7"=NA)
            breaks <- c(breaks,7)
            label <- c(label,"Balanced Sample: Post")
        }
        if(8%in%all){
            col <- c(col,"#A5CA18")
            col2 <- c(col2, "8"=NA)
            breaks <- c(breaks,8)
            label <- c(label,"Balanced Sample: Pre")
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
        
        if (gridOff == FALSE) {
            p <- p + geom_tile(colour="gray90", size=0.05, stat="identity") 
        } else {
            p <- p + geom_tile(stat="identity")
        }
        
        p <- p +
            labs(x = xlab, y = ylab,  
                title=main) +
            theme_bw() + 
            scale_fill_manual(NA, breaks = breaks, values = col, labels=label)

        #if(4%in%all) {
        #    p <- p + geom_point(aes(colour=res),size=0.5)
        #    p <- p + scale_color_manual(NA, breaks=breaks,
        #                                values=col2, labels=label)
        #}

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
              legend.title = element_blank(),
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
        
        if(length(all)>=3) {
            p <- p + guides(fill=guide_legend(nrow=2,byrow=TRUE))
        }
    }

    #suppressWarnings(print(p))
    if(return.test==TRUE){
        return(list(p=p,test.out=test.out))        
    }
    else{
        return(p)
    }


}