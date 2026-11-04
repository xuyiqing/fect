#!/usr/bin/env Rscript
# Reconstruct plotting inputs from checksum-verified public releases.
# No bootstrap is rerun. Original saved bootstrap summaries are used elsewhere.
suppressPackageStartupMessages({
  library(haven)
  library(digest)
  library(jsonlite)
  library(fect)
})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("Usage: Rscript recover_inputs.R BUNDLE_DIR OUTPUT_DIR")
bundle <- normalizePath(args[1],mustWork=TRUE)
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
out <- file.path(normalizePath(args[2],mustWork=TRUE),"inputs")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

fetch_source <- function(id, filename, expected_md5, expected_sha256=NULL) {
  path <- file.path(out, filename)
  if (!file.exists(path)) {
    old <- options(HTTPUserAgent="Mozilla/5.0", timeout=120)
    on.exit(options(old), add=TRUE)
    download.file(paste0("https://dataverse.harvard.edu/api/access/datafile/", id,
                         "?format=original"), path, mode="wb", quiet=TRUE)
  }
  md5 <- unname(tools::md5sum(path))
  sha <- digest(path, algo="sha256", file=TRUE, serialize=FALSE)
  stopifnot(identical(md5,expected_md5))
  if (!is.null(expected_sha256)) stopifnot(identical(sha,expected_sha256))
  list(path=path, file_id=id, md5=md5, sha256=sha)
}

# Figure 4: use the exact saved pooled GPS coefficients. Their tiny numeric-time
# terms are included at the displayed event day +7. Individual covariates are
# copied exactly as in the original panel builder (Euro is the omitted region).
lawall_source <- fetch_source(7416888, "fieldexpdata.csv",
  "e6b489ac0f26221425f0f2ee8d65e9ca")
lawall <- read.csv(lawall_source$path, stringsAsFactors=FALSE)
lawall <- lawall[lawall$treat_email==1 & lawall$treat_negative==1,]
stopifnot(nrow(lawall)==35976L)
gps_coef <- read.csv(file.path(bundle,"results/case1-lawall/cm/diagnostics/positivity/pooled_gps_coefficients.csv"), check.names=FALSE)
linear_predictor <- rep(0,nrow(lawall))
for (i in seq_len(nrow(gps_coef))) {
  term <- gsub("`","",gps_coef$term[i],fixed=TRUE)
  value <- if (term=="(Intercept)") rep(1,nrow(lawall)) else
    if (term==".cm_time_num") rep(7,nrow(lawall)) else
    if (term=="I(.cm_time_num^2)") rep(49,nrow(lawall)) else
    if (startsWith(term,"region_")) as.integer(lawall$dem_region==substring(term,8)) else
    lawall[[term]]
  stopifnot(length(value)==nrow(lawall),!anyNA(value))
  linear_predictor <- linear_predictor+gps_coef$estimate[i]*value
}
probability <- plogis(linear_predictor)
gps <- data.frame(subject_id=lawall$X, moderator=lawall$treat_issue,
                  event_time=7L, probability=probability, log_odds=linear_predictor)
gps_validation <- do.call(rbind,lapply(0:1,function(m) {
  p <- if(m==1) probability[gps$moderator==m] else 1-probability[gps$moderator==m]
  w <- 1/p
  data.frame(level=as.character(m),n=length(p),observed_probability_min=min(p),
             observed_probability_p01=unname(quantile(p,.01)),
             observed_probability_p05=unname(quantile(p,.05)),
             observed_probability_median=median(p),
             observed_probability_p95=unname(quantile(p,.95)),
             observed_probability_max=max(p),max_weight=max(w),
             ess=sum(w)^2/sum(w^2),ess_ratio=(sum(w)^2/sum(w^2))/length(w))
}))
gps_reference <- read.csv(file.path(bundle,"results/case1-lawall/cm/diagnostics/positivity/positivity_by_level_event_time.csv"))
gps_reference <- gps_reference[gps_reference$event_time==7,]
gps_reference <- gps_reference[match(gps_validation$level,gps_reference$level),]
gps_cols <- setdiff(names(gps_validation),"level")
gps_error <- max(abs(as.matrix(gps_validation[gps_cols])-as.matrix(gps_reference[gps_cols])))
stopifnot(gps_error<1e-7)
write.csv(gps,file.path(out,"fig4_gps.csv"),row.names=FALSE)

# Figures 8/9: reconstruct the author-model complete-case sample, then the
# original untreated-period gate and one common FE fit, with SE disabled.
kramon_source <- fetch_source(4101843,"Kramon_etal_analysisdata.dta",
 "70a75c7c4d22b421567140115a660265",
 "aefc584ec55aad12ba6bf241844726947e5488fa2e6565594badde8f827d462e")
k <- as.data.frame(read_dta(kramon_source$path))
k[] <- lapply(k,function(v)if(inherits(v,"haven_labelled"))as.numeric(v)else v)
k <- k[order(k$pupid,k$round),]
k$top2_female <- k$top2*k$female
k$top2_ageDM <- k$top2*k$ageDM
k$top2_parent_eduDM <- k$top2*k$parent_eduDM
columns <- c("pupid","round","ethnicimp","ethnicimp_urband","top2","female","age","ageDM",
 "parent_eduDM","top2_female","top2_ageDM","top2_parent_eduDM","voced","voced_treat_voucher",
 "coethnic_enum","KLPSg_w2","KLPS1_w1","KLPS1_w2","KLPS2_w1","KLPS2_w2","KLPS3_w1","KLPS3_w2","sample","weight","schid")
k <- k[k$sample==1 & complete.cases(k[columns]),columns]
stopifnot(nrow(k)==18571L,length(unique(k$pupid))==8709L)
has_control <- tapply(k$top2==0,k$pupid,any)
k <- k[unname(has_control[as.character(k$pupid)]),]
stopifnot(nrow(k)==17807L,length(unique(k$pupid))==8204L)
fit <- fect(ethnicimp_urband ~ top2+ageDM+voced+voced_treat_voucher+coethnic_enum,
 data=k,index=c("pupid","round"),method="fe",force="two-way",se=FALSE,
 parallel=FALSE,na.rm=TRUE,min.T0=1)
saveRDS(fit,file.path(out,"kramon_common_fit_no_bootstrap.rds"))
mod <- tapply(k$parent_eduDM,k$pupid,function(x)unique(x)[1])
moderator_matrix <- matrix(rep(unname(mod[as.character(fit$id)]),each=length(fit$rawtime)),
 nrow=length(fit$rawtime))
unit_matrix <- matrix(rep(fit$id,each=length(fit$rawtime)),nrow=length(fit$rawtime))
time_matrix <- matrix(rep(fit$rawtime,length(fit$id)),nrow=length(fit$rawtime))
keep <- !is.na(fit$eff)
gap_input <- data.frame(person_id=as.vector(unit_matrix)[keep],
 time=as.vector(time_matrix)[keep],moderator=as.vector(moderator_matrix)[keep],
 gap=as.vector(fit$eff)[keep],treatment=as.vector(fit$D.dat)[keep],
 event_time=as.vector(fit$T.on)[keep])
gap_input$phase <- ifelse(gap_input$treatment==1,"posttreatment",
 ifelse(!is.na(gap_input$event_time)&gap_input$event_time>=-1&gap_input$event_time<=0,
        "pretreatment","other"))
profile_ranges <- list(lower_parent_education=c(-8,0),higher_parent_education=c(0,11))
validation <- do.call(rbind,lapply(names(profile_ranges),function(profile){
 r <- profile_ranges[[profile]]
 do.call(rbind,lapply(c("posttreatment","pretreatment"),function(phase){
  vals <- gap_input$gap[gap_input$phase==phase&gap_input$moderator>=r[1]&gap_input$moderator<=r[2]]
  n <- length(vals);se <- sd(vals)/sqrt(n);ci <- qt(.975,n-1)*se
  data.frame(profile_id=profile,period=phase,estimate=mean(vals),
   dispersion_standard_error=se,ci_lower=mean(vals)-ci,ci_upper=mean(vals)+ci,unit_period_count=n)
 }))
}))
reference <- read.csv(file.path(bundle,"results/case3-kramon/em/tables/profile_effect_summary.csv"))
validation_key <- paste(validation$profile_id,validation$period)
ref_key <- paste(reference$profile_id,reference$period)
reference <- reference[match(validation_key,ref_key),]
numeric_cols <- c("estimate","dispersion_standard_error","ci_lower","ci_upper","unit_period_count")
profile_error <- max(abs(as.matrix(validation[numeric_cols])-as.matrix(reference[numeric_cols])))
write.csv(validation,file.path(out,"kramon_profile_validation.csv"),row.names=FALSE)
message("Maximum Kramon profile discrepancy: ",profile_error)
if(profile_error>1e-6)stop("Recovered fit does not match saved manuscript profiles; do not use inputs.")
write.csv(gap_input,file.path(out,"fig8_kramon_gaps.csv"),row.names=FALSE)
validation_manifest <- list(lawall_source=lawall_source,kramon_source=kramon_source,
 gps_max_absolute_reference_error=gps_error,kramon_max_absolute_profile_error=profile_error,
 kramon_fect_version=as.character(packageVersion("fect")),
 kramon_sample_rows=nrow(k),kramon_sample_units=length(unique(k$pupid)),
 kramon_att=mean(gap_input$gap[gap_input$treatment==1]),
 bootstrap_rerun=FALSE)
write_json(validation_manifest,file.path(out,"recovered_input_validation.json"),pretty=TRUE,auto_unbox=TRUE)
print(validation_manifest)
