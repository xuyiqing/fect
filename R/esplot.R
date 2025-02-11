## new coefplot
# x: time ATT CI.lower CI.upper count CI.lower.90 CI.upper.90
esplot <- function(data,# time ATT CI.lower CI.upper count
                   Period, 
                   Estimate,
                   SE,
                   CI.lower = NULL,
                   CI.upper = NULL,
                   Count = NULL,
                   fill.gap = TRUE, # use values 0 to fill the gaps of the dynamic treatment effects
                   start0 = FALSE, # default: 1 as the first post-treatment period
                   show.count = NULL, # whether show the bar of number of observations
                   stats = NULL, # a list of p-values
                   stats.labs = NULL, # a list of p-value labels
                   highlight.periods = NULL, # a list of periods
                   highlight.colors = NULL, # a list of colors
                   main = NULL,
                   xlim = NULL,
                   ylim = NULL,
                   xlab = NULL, 
                   ylab = NULL,
                   gridOff = FALSE,
                   stats.pos = NULL,
                   theme.bw = TRUE,
                   cex.main = NULL,
                   cex.axis = NULL,
                   cex.lab = NULL, 
                   cex.text = NULL,
                   axis.adjust = FALSE){
  
  scaleFUN <- function(x) sprintf("%.f", x)
  xmax<-xmin<-ymax<-ymin<- NULL

  data <- as.data.frame(data)

  # count
  if (is.logical(show.count) == FALSE & is.numeric(show.count)==FALSE & is.null(show.count)==FALSE) {
    stop("\"show.count\" is not a logical flag.")
  }
  if (is.null(show.count)==TRUE){
    show.count <- FALSE
  }
  if(show.count == TRUE & is.null(Count)){
    stop("\"Count\" is not specified.")
  }
  if(show.count == TRUE){
    if(!(Count %in% colnames(data))){
      stop("\"Count\" is not in the data.")      
    }
  }
  
  # gridOff
  if (is.logical(gridOff) == FALSE & !gridOff%in%c(0, 1)) {
    stop("\"gridOff\" is not a logical flag.")
  } 

  if (is.logical(fill.gap) == FALSE & !fill.gap%in%c(0, 1)) {
    stop("\"fill.gap\" is not a logical flag.")
  } 

  if (is.logical(axis.adjust) == FALSE & !axis.adjust%in%c(0, 1)) {
    stop("\"axis.adjust\" is not a logical flag.")
  } 
  
  # title
  if (is.null(main)==FALSE) {
    if (is.character(main) == FALSE) {
      stop("\"main\" is not a string.")
    } else {
      main <- main[1]
    }   
  }
  if (is.null(cex.main)==FALSE) {
    if (is.numeric(cex.main)==FALSE) {
      stop("\"cex.main\" is not numeric.")
    }
    cex.main <- 16 * cex.main
  } else {
    cex.main <- 16
  }
  
  # axis label
  if (is.null(xlab) == FALSE) {
    if (is.character(xlab) == FALSE) {
      stop("\"xlab\" is not a string.")
    }  
  }
  if (is.null(ylab) == FALSE) {
    if (is.character(ylab) == FALSE) {
      stop("\"ylab\" is not a string.")
    }  
  }
  if (is.null(cex.lab)==FALSE) {
    if (is.numeric(cex.lab)==FALSE) {
      stop("\"cex.lab\" is not numeric.")
    }
    cex.lab <- 14 * cex.lab
  } else {
    cex.lab <- 14
  }
  
  # axis number
  if (is.null(cex.axis)==FALSE) {
    if (is.numeric(cex.axis)==FALSE) {
      stop("\"cex.axis\" is not numeric.")
    }
    cex.axis <- 15 * cex.axis
  }  else {
    cex.axis <- 15
  }
  
  # text
  if (is.null(cex.text)==FALSE) {
    if (is.numeric(cex.text)==FALSE) {
      stop("\"cex.text\" is not numeric.")
    }
    cex.text <- 5 * cex.text
  }  else {
    cex.text <- 5
  }
  
  # statistics values
  if (is.null(stats)==FALSE){
    stats <- c(stats)
    if(!is.numeric(stats)){
      stop("The \"stats\" option must be numeric.")
    }
    stats <- ifelse(stats %% 1 != 0, round(stats,3), stats)
    stats <- as.character(stats)
    n.stats <- length(stats)

    if(length(stats.labs)!=n.stats){
      stop("The \"stats.labs\" option should have the same length as the \"stats\" option.")
    }
  }
  
  # stats positions
  if (!is.null(stats.pos)) {
    if (length(stats.pos) != 2) {
      stop(" \"stats.pos\" must be of length 2. ")
    }
    if (is.numeric(stats.pos[0])==FALSE) {
      stop("Elements of \"stats.pos\" are not numeric.")
    }else{
      if (is.numeric(stats.pos[1])==FALSE) {
        stop("Elements of \"stats.pos\" are not numeric.")    
      }
    }
  }
  
  # axis.adjust
  if (is.logical(axis.adjust) == FALSE & is.numeric(axis.adjust)==FALSE) {
    stop("\"axis.adjust\" is not a logical flag.")
  }
  if (axis.adjust == TRUE) {
    angle <- 45
    x.v <- 1
    x.h <- 1
  } else {
    angle <- 0
    x.v <- 0
    x.h <- 0
  }

  
  # xlim&ylim
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
  
  # axes labels
  if (is.null(xlab) == TRUE) {
    xlab <- paste("Time Relative to Treatment")   
  } 
  else if (xlab == "") {
    xlab <- NULL
  }
    
  if (is.null(ylab) == TRUE) {
    ylab <- paste("Effect on Y")
  } 
  else if (ylab == "") {
    ylab <- NULL
  }
  
  # y=0 line type
  lcolor <- "white"
  lwidth <- 2
  if (theme.bw == TRUE) {
    lcolor <- "#AAAAAA70"
    lwidth <- 1.5
  }
  
  # data
  data <- as.data.frame(data)
  if(!is.null(xlim)){
    data <- data[which(data[,Period]>=xlim[1] & data[,Period]<=xlim[2]),]
  }
  time <- data[,Period]
  ATT <- data[,Estimate]
  se <- data[,SE]
  if(!is.null(CI.lower)){
    CI.lower <- data[,CI.lower]
  }
  if(!is.null(CI.upper)){
    CI.upper <- data[,CI.upper]
  }
  if(!is.null(Count)){
    count.num <- data[,Count]
  }
  else{
    count.num <- rep(0,length(time))
  }
  if(!is.null(se)){
    if(is.null(CI.lower)){
      CI.lower <- ATT - 1.96*se
    }
    if(is.null(CI.upper)){
      CI.upper <- ATT + 1.96*se
    }
  }

  

  # add default zero
  time_lag <- max(time)-min(time)+1
  if(fill.gap == TRUE){
    if(time_lag>dim(data)[1]){
      time.add <- setdiff(c(min(time):max(time)),time)
      time <- c(time,time.add)
      ATT <- c(ATT,rep(0,length(time.add)))
      se <- c(se,rep(0,length(time.add)))
      CI.lower <- c(CI.lower,rep(0,length(time.add)))
      CI.upper <- c(CI.upper,rep(0,length(time.add)))
      count.num <- c(count.num,rep(0,length(time.add)))
    }    
  }

  if(!is.null(highlight.periods)){
    for(sub.period in highlight.periods){
      if(!sub.period %in% time){
        stop(paste0("Period ",sub.period," is out of the range of periods in the graph."))
      }
    }

    n.highlight <- length(highlight.periods)
    if(length(highlight.colors)!=n.highlight){
      stop("The \"highlight.colors\" option should have the same length as the \"highlight.periods\" option.")
    }
  }

  # Length of variables
  if (!(length(ATT)==length(time) & length(CI.lower)==length(time) & length(CI.upper)==length(time))){
    stop("The length of time, ATT, and uncertainty estimations must be the same.")
  }
  
  if (show.count == TRUE & !(length(count.num)==length(time))){
    stop("The length of time and observation counts must be the same.")
  }
  
  # count
  if (is.null(count.num) == TRUE){
    max.count <- NULL
    if (show.count == TRUE){
      stop("Input \"data\" must include variable \"Count\" if option \"show.count\" is TRUE.")
    }
  }
  else{
    max.count <- max(count.num)
  }
  
  
  data <- cbind.data.frame(time = time, ATT = ATT, CI.lower = CI.lower,
                           CI.upper=CI.upper,count=count.num)

  # set ylim if NULL
  if (is.null(ylim) == TRUE) {
    ylim <- c(min(data[,"CI.lower"], na.rm = TRUE) * 1.3, max(data[,"CI.upper"], na.rm = TRUE) * 1.3)
  } 
  ymin <- ylim[1]
  ymax <- ylim[2]
                         
#             $$\            $$\     
#           $$ |           $$ |    
#  $$$$$$\  $$ | $$$$$$\ $$$$$$\   
# $$  __$$\ $$ |$$  __$$\\_$$  _|  
# $$ /  $$ |$$ |$$ /  $$ | $$ |    
# $$ |  $$ |$$ |$$ |  $$ | $$ |$$\ 
# $$$$$$$  |$$ |\$$$$$$  | \$$$$  |
# $$  ____/ \__| \______/   \____/ 
# $$ |                             
# $$ |                             
# \__|                               

  p <- ggplot(data)
  
  max.count.pos <- time[which(count.num == max.count)]
  

  if(start0 == TRUE){
    best.pos <- 0
  }
  else{best.pos <- 1}
  
  if (length(max.count.pos)>1) {
    if (best.pos %in% max.count.pos) {
      max.count.pos <- best.pos
    } 
    else if ((1-best.pos) %in% max.count.pos) {
      max.count.pos <- 1-best.pos
    } 
    else {
      max.count.pos <- max.count.pos[1]
    }
  }
  
  # height of the histogram
  rect.length <- (ymax - ymin)/6
  rect.min <- ymin 
  
  # xlab and ylab
  p <- p + xlab(xlab) +  ylab(ylab) 
  
  # theme
  if (theme.bw == TRUE) {
    p <- p + theme_bw() 
  }
  
  ## grid
  if (gridOff == TRUE) {
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  
  # horizontal 0 line
  p <- p + geom_hline(yintercept = 0, 
                      colour = lcolor, 
                      linewidth = lwidth)
  
  # vertical 0 line
  if(start0 == FALSE){
    p <- p + geom_vline(xintercept = 0.5, colour=lcolor,linewidth = lwidth*0.7) 
  }
  else{
    p <- p + geom_vline(xintercept = -0.5, colour=lcolor,linewidth = lwidth*0.7) 
  }
     
  p <- p + geom_pointrange(data = data, aes(x = time, y = ATT, ymin=CI.lower, ymax=CI.upper), lwd=0.6, color="black", fill="black",fatten = 2)  
  T0 <- which(data[,"time"] == 0)

  
  
  maintext <- "Estimated Dynamic Treatment Effects"

  if(!is.null(highlight.periods)){
    kk <- 1
    data.highlight.plot <- list()
    for(sub.period in highlight.periods){
      sub.time <- sub.period
      sub.index <- which(time==sub.time)
      sub.ATT <- ATT[sub.index]
      sub.CI.lower <- CI.lower[sub.index]
      sub.CI.upper <- CI.upper[sub.index]
      sub.color <- highlight.colors[kk]
      data.highlight.plot[[kk]] <- cbind.data.frame(time=sub.time,
                                                    ATT=sub.ATT,
                                                    CI.lower=sub.CI.lower,
                                                    CI.upper=sub.CI.upper)
      
      p <- p + geom_pointrange(data=data.highlight.plot[[kk]],
                               aes(x = time, y = ATT, ymin=CI.lower, ymax=CI.upper), lwd=0.6, 
                               color=sub.color, fill=sub.color,fatten = 2)  
      kk <- kk + 1
    }
  }
  

  # axes
  p <- p + theme(axis.title=element_text(size=cex.lab),
                 axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
                 axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
                 axis.text = element_text(color="black", size=cex.axis),
                 axis.text.x = element_text(size = cex.axis, angle = angle, hjust=x.h, vjust=x.v),
                 axis.text.y = element_text(size = cex.axis),
                 plot.title = element_text(size = cex.main, hjust = 0.5, face="bold", margin = margin(10, 0, 10, 0)))
  
  # histogram
  if (show.count == TRUE) {
    data[,"xmin"] <- data[,"time"] - 0.2
    data[,"xmax"] <- data[,"time"] + 0.2
    data[,"ymin"] <- rep(rect.min, dim(data)[1])
    data[,"ymax"] <- rect.min + (data[,"count"]/max.count) * 0.6 * rect.length
    xx <- range(data$time)
    p <- p + geom_rect(data = data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                       fill = "grey70", colour = "grey69", alpha = 0.4, size = 0.2)
    p <- p + annotate("text", x = max.count.pos, #- 0.02 * (xx[2]-xx[1]), 
                      y = max(data$ymax) + 0.2 * rect.length, 
                      label = max.count, size = cex.text * 0.8, hjust = 0.5)                
  }
  
  # stats.pos
  if (!is.null(stats)){
    if (is.null(stats.pos)) {
      stats.pos[1] <- min(data[,"time"], na.rm = 1)
      stats.pos[2] <- ifelse(is.null(ylim), max(data[,"CI.upper"], na.rm = 1)*1.16, ylim[1])
    }
  }
  
  # p.label
  p.label <- NULL
  if(!is.null(stats)){
    for (ii in 1:length(stats)) {
      p.label <- paste0(p.label,stats.labs[ii],":",stats[ii],"\n")
    }
    p <- p + annotate("text", x = stats.pos[1], y = stats.pos[2], 
                      label = p.label, size = cex.text * 0.8, hjust = 0)      
  }


  
  # xlim & ylim  
  p <- p + ylim(ylim = ylim)
  

  # title
  if (is.null(main) == TRUE) {
    p <- p + ggtitle(maintext)
  } 
  else if (main!=""){
    p <- p + ggtitle(main)
  }
  
  if(dim(data)[1]>4){
    p <- p + scale_x_continuous(labels=scaleFUN)
  }
  else{
    p <- p + scale_x_continuous(breaks=c(data[,'time']))
  }

  
  
  return(p)
}

