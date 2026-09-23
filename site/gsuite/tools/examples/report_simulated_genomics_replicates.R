# Summarize all ten planned attempts; no new simulations or model fits.
root <- getOption("gsuite.qualify.root", "build/examples/simulated-genomics/replicates")
plan <- read.csv(file.path(root,"plan.csv"))
paths <- file.path(root,sprintf("replicate-%02d.rds",plan$replicate))
stopifnot(nrow(plan)==10L,all(file.exists(paths)))
runs <- lapply(paths,readRDS)
stopifnot(identical(vapply(runs,`[[`,integer(1),"seed"),plan$seed))
status <- do.call(rbind,lapply(runs,function(x)data.frame(replicate=x$replicate,seed=x$seed,
  status=x$status,error=if(is.null(x$error)) "" else x$error,
  seconds=as.numeric(difftime(x$finished,x$started,units="secs")))))
write.csv(status,file.path(root,"status.csv"),row.names=FALSE)
completed <- runs[vapply(runs,function(x)x$status=="completed",logical(1))]
if(!length(completed)) stop("No completed replicates; inspect status.csv")
stopifnot(all(vapply(completed,function(x)identical(x$scripts,completed[[1]]$scripts),logical(1))))
stopifnot(all(vapply(completed,function(x)all(vapply(x$results$provenance$stages,
  function(s)identical(s$seed,x$seed),logical(1))),logical(1))))
native_reference <- completed[[1]]$results$provenance$stages$data$native_md5
stopifnot(all(vapply(completed,function(x)all(vapply(x$results$provenance$stages,
  function(s)identical(s$native_md5,native_reference) &&
    identical(unname(s$script_md5),unname(x$scripts[[1L]])),logical(1))),logical(1))))
ld_identities <- vapply(completed,function(x)
  paste(unname(x$results$provenance$LD_checksums),collapse=":"),character(1))
stopifnot(!anyDuplicated(ld_identities))
checks <- do.call(rbind,lapply(completed,function(x)data.frame(replicate=x$replicate,seed=x$seed,
  max_ols_error=unname(x$results$checks["max_ols_error"]),
  truth_score_error=unname(x$results$checks["truth_score_error"]))))
stopifnot(all(is.finite(as.matrix(checks[3:4]))),all(as.matrix(checks[3:4])<1e-9))
write.csv(checks,file.path(root,"checks.csv"),row.names=FALSE)
write.csv(data.frame(script=names(completed[[1]]$scripts),md5=unname(completed[[1]]$scripts)),
  file.path(root,"source-hashes.csv"),row.names=FALSE)
combine <- function(name) do.call(rbind,lapply(completed,function(x)
  cbind(replicate=x$replicate,seed=x$seed,x$results[[name]])))
tables <- setNames(lapply(c("estimation","prediction","association","marginal_effects",
  "bayesian","posterior","credible_sets","coloc","pathways","availability","tails","timings"),combine),
  c("estimation","prediction","association","marginal_effects","bayesian","posterior",
    "credible_sets","coloc","pathways","availability","tails","timings"))
for(k in names(tables)) write.csv(tables[[k]],file.path(root,paste0(k,".csv")),row.names=FALSE)
e <- tables$estimation
e$error <- e$estimate-e$truth
e$interval_available <- is.finite(e$lower)&is.finite(e$upper)&is.finite(e$truth)
e$covered <- ifelse(e$interval_available,e$lower<=e$truth & e$upper>=e$truth,NA)
# Trait/region rows within a replicate are dependent; retain these labels.
# Only global trait-specific estimates get an across-replicate coverage summary.
global <- e[e$scope=="reference",]
groups <- split(global,interaction(global$analysis,global$quantity,global$trait1,
  ifelse(is.na(global$trait2),"none",global$trait2),drop=TRUE))
accuracy <- do.call(rbind,lapply(groups,function(d) data.frame(
  analysis=d$analysis[1],quantity=d$quantity[1],trait1=d$trait1[1],trait2=d$trait2[1],
  planned=nrow(plan),completed=nrow(d),finite=sum(is.finite(d$error)),
  mean_truth=mean(d$truth),mean_estimate=mean(d$estimate,na.rm=TRUE),
  bias=mean(d$error,na.rm=TRUE),rmse=sqrt(mean(d$error^2,na.rm=TRUE)),
  intervals=sum(d$interval_available),covered=sum(d$covered,na.rm=TRUE),
  interval_kind=paste(unique(d$interval_kind),collapse=";"))))
write.csv(accuracy,file.path(root,"accuracy.csv"),row.names=FALSE)
saveRDS(list(plan=plan,status=status,tables=tables,accuracy=accuracy,checks=checks),file.path(root,"summary.rds"))
draw <- function(name,code,width=1800,height=1250) {
  png(file.path(root,paste0(name,".png")),width=width,height=height,res=150)
  tryCatch(force(code),finally=dev.off())
}
cols <- c("#28678b","#bf6b2d","#38836c","#8661a1","#ad4262","#686868")
method_label <- function(x) {
  key <- sub("_(h2|rg|cov)$","",x)
  lookup <- c(ldsc="LDSC",gnova="GNOVA",sumher="SumHer",
    ldsc_annotation="LDSC +\nannotations",bayesc="BayesC",bayesr="BayesR",
    hess="HESS",supergnova="SUPERGNOVA",lava="LAVA",ridge="Ridge")
  unname(lookup[key])
}
strip <- function(values,groups,main,ylab,reference=NULL,labels=NULL,ylim=NULL) {
  groups<-factor(groups,levels=unique(groups)); lev<-levels(groups)
  if(is.null(labels)) labels<-lev
  if(is.null(ylim)) ylim<-extendrange(range(c(values,reference),finite=TRUE))
  plot(1,1,type="n",xlim=c(.6,length(lev)+.4),ylim=ylim,xaxt="n",
    xlab="",ylab=ylab,main=main)
  axis(1,seq_along(lev),labels,las=2,cex.axis=.72)
  if(!is.null(reference)) abline(h=reference,lty=2,col="#777777")
  for(j in seq_along(lev)) {
    v<-values[groups==lev[j]]; offsets<-seq(-.13,.13,length.out=length(v))
    points(j+offsets,v,pch=16,cex=.7,col=adjustcolor(cols[(j-1)%%length(cols)+1],.65))
    segments(j-.22,mean(v,na.rm=TRUE),j+.22,mean(v,na.rm=TRUE),lwd=2)
  }
}
draw("replicate-estimation", {
  par(mfrow=c(2,2),mar=c(7,4.8,3,1),las=1)
  for(t in c("A","B","C")) {
    d<-global[global$quantity=="h2" & global$trait1==t,]
    strip(d$error,method_label(d$analysis),paste("Heritability error: trait",t),
      "Estimate - population truth",0)
  }
  d<-global[global$quantity=="rg",]
  # One mean error per replicate/method across its three trait pairs.
  a<-aggregate(error~replicate+analysis,d,mean)
  strip(a$error,method_label(a$analysis),"Genetic correlation: mean pair error",
    "Mean error across three trait pairs",0)
})
draw("replicate-prediction-mapping", {
  par(mfrow=c(2,2),mar=c(5.8,4.8,3,1),las=1)
  d<-tables$prediction
  strip(d$cor_genetic,method_label(d$method),"Held-out genetic prediction","Correlation",ylim=c(0,1))
  strip(d$calibration_slope,method_label(d$method),"Prediction calibration","Truth-on-prediction slope",1)
  d<-tables$credible_sets
  labels<-paste(d$region,d$trait)
  strip(as.numeric(d$causal_covered),labels,"Component sets: causal coverage","Covered (0 / 1)",
    ylim=c(-.05,1.05))
  d<-tables$coloc
  probability<-ifelse(d$truth=="H4",d$H4,d$H3)
  strip(probability,paste(d$region,d$truth),"Colocalisation: true hypothesis",
    "Posterior probability",ylim=c(0,1.05))
})
draw("replicate-uncertainty", {
  par(mfrow=c(2,2),mar=c(7,4.8,3,1),las=1)
  d<-accuracy[accuracy$quantity=="h2" & accuracy$trait1=="A",]
  fraction<-ifelse(d$intervals>0,d$covered/d$intervals,NA_real_)
  x<-barplot(fraction,names.arg=method_label(d$analysis),
    ylim=c(0,1.2),las=2,col=cols,ylab="Fraction containing population truth",
    main="Trait A intervals: observed coverage",cex.names=.7)
  text(x,ifelse(is.finite(fraction),fraction+.06,.06),
    ifelse(d$intervals>0,paste0(d$covered,"/",d$intervals),"Unavailable"),cex=.8)
  abline(h=.95,lty=2);mtext("Different interval constructions; ten simulations",side=3,cex=.7,line=.1)
  d<-tables$availability[grepl("_cov$",tables$availability$analysis),]
  strip(d$available_se/d$estimates,method_label(d$analysis),
    "Regional uncertainty availability","Fraction of 60 SEs available",ylim=c(0,1.08))
  d<-aggregate(cbind(controlled_available,approximate_available,requested,fallback)~replicate,tables$tails,sum)
  strip(c(d$requested-d$controlled_available,d$fallback),
    rep(c("Controlled\nunavailable","Approximate\nmoment fallback"),each=nrow(d)),
    "Gene-tail diagnostics","Gene / trait count (of 3,000)",0)
  d<-global[global$quantity=="rg",]
  d$outside<-as.numeric(d$estimate < -1 | d$estimate > 1)
  a<-aggregate(outside~replicate+analysis,d,sum)
  strip(a$outside,method_label(a$analysis),"Out-of-range genetic correlations",
    "Count / three trait pairs",0,ylim=c(-.1,3.1))
})
print(status);print(accuracy)
cat("REPLICATE_REPORT|completed\n")
