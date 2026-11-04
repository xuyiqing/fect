#!/usr/bin/env Rscript
# Plotting-only revision: preserve the saved estimates and interval definitions.
# Usage: Rscript redraw_figures.R BUNDLE_DIR OUTPUT_DIR
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
args <- commandArgs(TRUE)
if(length(args)!=2L)stop("Usage: Rscript redraw_figures.R BUNDLE_DIR OUTPUT_DIR")
bundle <- normalizePath(args[1], mustWork=TRUE)
out <- normalizePath(args[2], mustWork=TRUE)
results <- file.path(bundle,"results")
for(d in c("figures","panels","tables"))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
read <- function(path)read.csv(file.path(results,path),check.names=FALSE)
pub_theme <- function(){
  theme_classic(base_size=13,base_family="sans") +
    theme(axis.title=element_text(size=17,color="black"),
      axis.text=element_text(size=12,color="black"),
      axis.title.x=element_text(margin=margin(t=10)),
      axis.title.y=element_text(margin=margin(r=10)),
      axis.line=element_line(linewidth=.35,color="#555555"),
      axis.ticks=element_line(linewidth=.35,color="#555555"),
      plot.margin=margin(14,18,12,15),
      legend.position="inside",legend.position.inside=c(.02,.98),
      legend.justification=c(0,1),legend.direction="horizontal",
      legend.background=element_blank(),legend.key=element_blank(),
      legend.title=element_blank(),legend.text=element_text(size=11),
      legend.key.width=grid::unit(.9,"lines"),legend.key.height=grid::unit(.8,"lines"),
      legend.margin=margin(0,0,0,0),plot.title=element_blank(),plot.subtitle=element_blank(),
      plot.caption=element_blank())
}
yrange <- function(d){r<-range(c(d$ci_lower,d$ci_upper,0),na.rm=TRUE);r+c(-1,1)*diff(r)*.13}
save <- function(p,stem,width=7,height=5,folder="figures"){
  base<-file.path(out,folder,stem)
  ggsave(paste0(base,".pdf"),p,width=width,height=height,device=grDevices::cairo_pdf,bg="white")
  ggsave(paste0(base,".png"),p,width=width,height=height,dpi=600,device=ragg::agg_png,bg="white",limitsize=FALSE)
}
support <- function(d,xlim,breaks,labels=NULL,xlab=NULL,constant=FALSE){
  d<-d[!duplicated(d$event_time),]
  d$height<-d$n/max(d$n)
  p<-ggplot(d,aes(event_time,height)) + geom_col(width=.44,fill="#dddddd",color=NA)+
    scale_x_continuous(limits=xlim,breaks=breaks,labels=if(is.null(labels))waiver() else labels,expand=expansion(mult=0))+
    scale_y_continuous(limits=c(0,1.6),breaks=NULL,expand=expansion(mult=0))+
    labs(x=xlab,y=NULL)+pub_theme()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),
      plot.margin=margin(0,18,9,15))
  if(constant)p<-p+annotate("text",x=mean(range(d$event_time)),y=1.28,label=format(unique(d$n)[1],big.mark=","),size=3.7)
  else p<-p+geom_text(aes(label=format(n,big.mark=",")),vjust=0,nudge_y=.10,size=3.5)
  p
}
discrete_plot <- function(d,labels,ylim,xlab){
  d$event_time<-match(d$level,sort(unique(d$level)))
  d$n<-d$n_cells
  xlim<-c(.5,2.5) # each estimate sits one-quarter of the panel width from the edge
  top<-ggplot(d,aes(event_time,estimate))+
    geom_hline(yintercept=0,color="black",linewidth=.65)+
    geom_linerange(aes(ymin=ci_lower,ymax=ci_upper),linewidth=.65,color="black")+
    geom_point(shape=21,fill="white",size=3.2,stroke=.9,color="black")+
    scale_x_continuous(limits=xlim,breaks=1:2,expand=expansion(mult=0))+
    coord_cartesian(ylim=ylim)+labs(x=NULL,y="Estimate")+pub_theme()+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),plot.margin=margin(14,18,0,15))
  (top/support(d,xlim,1:2,labels,xlab)) + plot_layout(heights=c(4,.9))
}
combine <- function(left,right) (left|right)+plot_layout(widths=c(1,1))

# Fig. 2: same scale across day +7 and held-out day -1; separate count strip.
cmend<-read("case1-lawall/cm/tables/cm_hte_day7_and_day_minus1.csv")
lim2<-yrange(cmend)
f2a<-discrete_plot(subset(cmend,event_time==7),c("Party cue","Issue cue"),lim2,"Cue type")
f2b<-discrete_plot(subset(cmend,event_time==-1),c("Party cue","Issue cue"),lim2,"Cue type")
save(f2a,"lawall-cm-day7",folder="panels");save(f2b,"lawall-cm-heldout-day-minus1",folder="panels")
save(combine(f2a,f2b),"Fig02",13,4.9)

# Fig. 3: one recipient-count strip per cue, EM/CM distinguished by point shape.
em<-read("case1-lawall/cm/tables/em_dynamic_heldout.csv");em$method<-"EM"
cm<-read("case1-lawall/cm/tables/cm_dynamic_with_ci.csv");cm$method<-"CM"
fields<-c("event_time","level","estimate","ci_lower","ci_upper","n_cells","method")
dyn<-rbind(em[,fields],cm[,fields]);lim3<-yrange(dyn)
dynamic_cue<-function(level){
  d<-dyn[as.character(dyn$level)==as.character(level),]
  d$method<-factor(d$method,levels=c("EM","CM"));d$xx<-d$event_time+ifelse(d$method=="EM",-.11,.11)
  top<-ggplot(d,aes(xx,estimate,shape=method,color=method))+
    geom_hline(yintercept=0,color="black",linewidth=.65)+
    geom_vline(xintercept=.5,color="#999999",linetype="dashed",linewidth=.35)+
    geom_linerange(aes(ymin=ci_lower,ymax=ci_upper),linewidth=.45)+
    geom_point(size=2.4,stroke=.7,fill="white")+
    annotate("point",x=0,y=0,shape=23,size=2.8,fill="white",color="black",stroke=.7)+
    scale_color_manual(values=c(EM="#777777",CM="black"))+
    scale_shape_manual(values=c(EM=16,CM=21))+
    scale_x_continuous(limits=c(-7.75,14.75),breaks=c(-7,-4,-1,0,1,4,7,10,14),expand=expansion(mult=0))+
    coord_cartesian(ylim=lim3)+labs(x=NULL,y="Estimate")+pub_theme()+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),plot.margin=margin(14,18,0,15))
  day7<-d[d$event_time==7,]
  label<-sprintf("Day 7: EM %.5f; CM %.5f",day7$estimate[day7$method=="EM"],day7$estimate[day7$method=="CM"])
  top<-top+annotate("text",x=-7,y=lim3[2]-.11*diff(lim3),label=label,hjust=0,size=3.3,color="black")
  count<-d[d$method=="CM",];count$n<-count$n_cells
  (top/support(count,c(-7.75,14.75),c(-7,-4,-1,0,1,4,7,10,14),xlab="Event day",constant=TRUE))+
    plot_layout(heights=c(4,.85))
}
f3a<-dynamic_cue(0);f3b<-dynamic_cue(1)
save(f3a,"lawall-em-cm-dynamic-party",folder="panels");save(f3b,"lawall-em-cm-dynamic-issue",folder="panels")
save(combine(f3a,f3b),"Fig03",14,5.1)

# Fig. 4: saved GPS coefficients applied to recovered original covariates.
gps_path<-file.path(out,"inputs","fig4_gps.csv")
if(file.exists(gps_path)){
  gps<-read.csv(gps_path)
  gps$level<-gps$moderator
  if(!all(c("probability","level")%in%names(gps)))stop("Fig4 requires probability, level")
  if(any(gps$probability<.4|gps$probability>.6))stop("GPS observations outside requested display window: do not hide them.")
  bins<-seq(.4,.6,length.out=41)
  histdata<-do.call(rbind,lapply(0:1,function(k){
    v<-gps$probability[gps$level==k];h<-hist(v,breaks=bins,plot=FALSE)
    data.frame(x=h$mids,y=h$counts/sum(h$counts)*ifelse(k==1,1,-1),group=ifelse(k==1,"Issue cue","Party cue"))
  }))
  f4<-ggplot(histdata,aes(x,y,fill=group))+
    geom_col(width=diff(bins)[1]*.94)+geom_hline(yintercept=0,color="black",linewidth=.6)+
    geom_vline(xintercept=.5,color="#777777",linetype="dashed",linewidth=.4)+
    scale_fill_manual(values=c("Issue cue"="#333333","Party cue"="#bcbcbc"))+
    scale_x_continuous(limits=c(.4,.6),breaks=seq(.4,.6,.05),expand=expansion(mult=0))+
    scale_y_continuous(labels=function(x)paste0(abs(round(x*100)),"%"),expand=expansion(mult=.15))+
    labs(x="Probability of issue cue",y="Share within cue group")+pub_theme()
  save(f4,"Fig04",7,4.8);save(f4,"lawall-cm-positivity-overlap",folder="panels")
}

# Fig. 5: retain the original native fixed-gap t intervals (not paired-bootstrap intervals).
d5<-read("case2-distelhorst-locke/grouped-em/tables/profile_effect_summary.csv")
d5$level<-ifelse(d5$profile_id=="clothing",1,0);d5$n_cells<-d5$unit_period_count
lim5<-yrange(d5)
f5a<-discrete_plot(subset(d5,period=="posttreatment"),c("Non-clothing","Clothing"),lim5,"Industry")
f5b<-discrete_plot(subset(d5,period=="pretreatment"),c("Non-clothing","Clothing"),lim5,"Industry")
save(f5a,"distelhorst-em-post-groups",folder="panels");save(f5b,"distelhorst-em-pre-groups",folder="panels")
save(combine(f5a,f5b),"Fig05",13,4.9)

# General paired event-time plot: no capped error bars, half-period edge padding.
event_plot<-function(d,groups,shapes=c(16,21),colors=c("black","#777777"),ylim=NULL){
  d$group<-factor(d$group,levels=groups)
  d$xx<-d$event_time+ifelse(d$group==groups[1],-.11,.11)
  if(is.null(ylim))ylim<-yrange(d)
  ggplot(d,aes(xx,estimate,color=group,shape=group))+
    geom_hline(yintercept=0,color="black",linewidth=.65)+
    geom_vline(xintercept=.5,color="#999999",linetype="dashed",linewidth=.35)+
    geom_linerange(aes(ymin=ci_lower,ymax=ci_upper),linewidth=.6)+
    geom_point(size=2.8,fill="white",stroke=.8)+
    scale_color_manual(values=setNames(colors,groups))+scale_shape_manual(values=setNames(shapes,groups))+
    scale_x_continuous(limits=range(d$event_time)+c(-.65,.65),breaks=sort(unique(d$event_time)),expand=expansion(mult=0))+
    coord_cartesian(ylim=ylim)+labs(x="Event period",y="Estimate")+pub_theme()
}
d6<-read("case2-distelhorst-locke/paired-em/minT0-2/tables/dynamic_industry_profiles.csv")
d6<-d6[d6$profile_id %in% c("clothing","non_clothing"),]
d6$group<-ifelse(d6$profile_id=="clothing","Clothing","Non-clothing")
f6<-event_plot(d6,c("Non-clothing","Clothing"))
save(f6,"Fig06",7.5,5);save(f6,"distelhorst-dynamic-em",folder="panels")

d7<-read("case2-distelhorst-locke/factory-trend/tables/dynamic_effects.csv")
d7$group<-ifelse(grepl("^ordinary",d7$spec),"Ordinary FE","Factory trends")
lim7<-yrange(d7)
f7a<-event_plot(subset(d7,min_T0==3),c("Ordinary FE","Factory trends"),ylim=lim7)
f7b<-event_plot(subset(d7,min_T0==2),c("Ordinary FE","Factory trends"),ylim=lim7)
save(combine(f7a,f7b),"Fig07",13,4.9)
save(combine(f7a,f7b),"distelhorst-factory-trend-dynamic",13,4.9,folder="panels")

# Fig. 8: original fixed-gap loess/predict pointwise interval, unchanged.
gap_path<-file.path(out,"inputs","fig8_kramon_gaps.csv")
if(file.exists(gap_path)){
  gaps<-read.csv(gap_path)
  continuous<-function(phase){
    d<-gaps[gaps$phase==phase&is.finite(gaps$moderator)&is.finite(gaps$gap),]
    d<-d[order(d$moderator),]
    fit<-loess(gap~moderator,data=d)
    pr<-predict(fit,newdata=d$moderator,se=TRUE)
    dd<-data.frame(x=d$moderator,estimate=as.numeric(pr$fit),ci_lower=as.numeric(pr$fit-qt(.975,pr$df)*pr$se.fit),ci_upper=as.numeric(pr$fit+qt(.975,pr$df)*pr$se.fit))
    dd<-dd[!duplicated(dd$x),]
    p<-ggplot(dd,aes(x,estimate))+geom_hline(yintercept=0,color="black",linewidth=.65)+
      geom_ribbon(aes(ymin=ci_lower,ymax=ci_upper),fill="#d0d0d0",alpha=.65)+
      geom_line(color="black",linewidth=.75)+scale_x_continuous(limits=c(-8.7,11.7),breaks=c(-8,-4,0,4,8,11),expand=expansion(mult=0))+
      labs(x=NULL,y="Estimate")+pub_theme()+
      theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),plot.margin=margin(14,18,0,15))
    # Preserve the original pretty(..., n=50) support bins, in a separate strip.
    hh<-hist(d$moderator,breaks=pretty(range(d$moderator),n=50),plot=FALSE)
    hd<-data.frame(xmin=head(hh$breaks,-1),xmax=tail(hh$breaks,-1),count=hh$counts)
    peak<-which.max(hd$count)
    histp<-ggplot(hd)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=0,ymax=count),fill="#d8d8d8",color=NA)+
      annotate("text",x=mean(c(hd$xmin[peak],hd$xmax[peak])),y=hd$count[peak]*1.15,label=hd$count[peak],size=3.5)+
      scale_x_continuous(limits=c(-8.7,11.7),breaks=c(-8,-4,0,4,8,11),expand=expansion(mult=0))+
      scale_y_continuous(breaks=NULL,expand=expansion(mult=c(0,.13)))+
      labs(x="Parental education (centered years)",y=NULL)+pub_theme()+
      theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(0,18,9,15))
    write.csv(dd,file.path(out,"tables",paste0("fig8_",phase,"_loess.csv")),row.names=FALSE)
    (p/histp)+plot_layout(heights=c(4,.85))
  }
  f8a<-continuous("posttreatment");f8b<-continuous("pretreatment")
  save(f8a,"kramon-em-posttreatment",folder="panels");save(f8b,"kramon-em-pretreatment",folder="panels")
  save(combine(f8a,f8b),"Fig08",14,5.3)
}

# Fig. 9: remove unobserved period -2; same y-scale, no pooled-average line.
d9<-read("case3-kramon/em/tables/event_time_support.csv")
d9$event_time<-d9$relative_time;d9$n<-d9$effect_unit_periods
lim9<-yrange(d9)
parent_dynamic<-function(id){
  d<-d9[d9$profile_id==id,]
  top<-ggplot(d,aes(event_time,estimate))+
    geom_hline(yintercept=0,color="black",linewidth=.65)+
    geom_vline(xintercept=.5,color="#999999",linetype="dashed",linewidth=.35)+
    geom_linerange(aes(ymin=ci_lower,ymax=ci_upper),color="black",linewidth=.65)+geom_point(color="black",size=2.8)+
    scale_x_continuous(limits=c(-1.65,2.65),breaks=-1:2,expand=expansion(mult=0))+
    coord_cartesian(ylim=lim9)+labs(x=NULL,y="Estimate")+pub_theme()+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),plot.margin=margin(14,18,0,15))
  (top/support(d,c(-1.65,2.65),-1:2,xlab="Event period"))+plot_layout(heights=c(4,.9))
}
f9a<-parent_dynamic("lower_parent_education");f9b<-parent_dynamic("higher_parent_education")
save(f9a,"kramon-em-dynamic-lower",folder="panels");save(f9b,"kramon-em-dynamic-higher",folder="panels")
save(combine(f9a,f9b),"Fig09",13,4.9)

# Caption metadata lives outside images (not in-plot title/subtitle/note).
captions<-data.frame(figure=2:9,
 left=c("Day +7 CM","Party cue","Issue above; party below","Posttreatment","Non-clothing versus clothing","At least 3 untreated periods","Posttreatment","At or below mean parental education"),
 right=c("Held-out day -1 CM","Issue cue","","Pretreatment","","At least 2 untreated periods","Pretreatment","Above mean parental education"),
 note=c("Original 200-draw recipient bootstrap percentile intervals.","EM and CM preserve their original interval algorithms. Diamond at day 0 is a zero reference with no CI, not an observation. N=17,988 at observed dates; no day-0 count. The original day-7 reference is shown numerically, not as an average line.","Probability scale 0.4-0.6; no recovered score falls outside this interval. Same original target recipients and moderator GPS.","Original native fixed-gap t intervals, not the separate paired-bootstrap contrast intervals. Counts are factory-year cells.","Original 1000 paired-factory-bootstrap normal intervals; min.T0=2.","Original 50-draw factory-bootstrap pilot. Left/right samples differ as before. Same vertical scale.","Original fixed-gap loess pointwise intervals; neither full-pipeline bootstrap nor simultaneous bands. Support histogram below the estimate panel.","Original native fixed-gap t intervals and counts. Period -2 is absent in the input. Same vertical scale across groups."))
write.csv(captions,file.path(out,"captions.csv"),row.names=FALSE)
cat("Rendered:",paste(list.files(file.path(out,"figures"),pattern="pdf$"),collapse=", "),"\n")
