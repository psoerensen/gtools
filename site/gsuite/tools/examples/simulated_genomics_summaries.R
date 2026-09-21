# Sourced by simulated_genomics_workflow.R after fitting. No additional fits.
# Population LD gives marginal GWAS effects; direct effects are used for BLR.
marginal_truth <- do.call(rbind,lapply(seq_len(nr),function(r) {
  R <- rho[r]^abs(outer(seq_len(size),seq_len(size),"-"))
  R %*% data$simulation$B[regions[[r]],,drop=FALSE]
}))
association_effects <- do.call(rbind,lapply(traits,function(t) {
  truth <- marginal_truth[,t]; estimate <- gwas$stat[[t]]$beta
  data.frame(trait=t,correlation=cor(estimate,truth),
    rmse=sqrt(mean((estimate-truth)^2)),
    slope=cov(estimate,truth)/var(truth))
}))

comparison <- estimates[estimates$scope %in% c("reference","annotation_category","region") &
  estimates$quantity %in% c("h2","cov","rg"),]
comparison$truth <- vapply(seq_len(nrow(comparison)),function(i) {
  e <- comparison[i,]; a <- e$trait1; b <- e$trait2
  G <- if(e$scope=="region") data$local_G[[match(e$region,names(regions))]] else
    if(e$scope=="annotation_category") {
      selected <- switch(e$component,enriched=1:80,background=81:nr,
        stop("Unknown annotation category"))
      Reduce(`+`,data$local_G[selected])
    } else data$total_G
  if(e$quantity=="h2") return(G[a,a]/data$total_V[a])
  if(e$quantity=="cov") return(G[a,b]/sqrt(data$total_V[a]*data$total_V[b]))
  if(G[a,a]<=0 || G[b,b]<=0) return(NA_real_)
  G[a,b]/sqrt(G[a,a]*G[b,b])
},numeric(1))
# Preserve method-specific uncertainty; do not manufacture intervals for LAVA.
comparison$lower <- ifelse(comparison$se_available,
  comparison$estimate-1.96*comparison$se,NA_real_)
comparison$upper <- ifelse(comparison$se_available,
  comparison$estimate+1.96*comparison$se,NA_real_)
comparison$interval_kind <- ifelse(comparison$se_available,
  paste0("normal_approximation:",comparison$se_method),"unavailable")
estimation <- comparison[c("analysis","scope","component","region","quantity","trait1",
  "trait2","estimate","truth","lower","upper","interval_kind","uncertainty_reason")]
for(method in names(bayesian)) {
  p <- posterior_table[posterior_table$method==method & posterior_table$parameter=="h2[A]",]
  estimation <- rbind(estimation,data.frame(analysis=method,scope="reference",component="all",
    region=NA_character_,quantity="h2",trait1="A",trait2=NA_character_,estimate=p$mean,
    truth=data$truth_h2["A"],lower=p$lower,upper=p$upper,
    interval_kind="conditional_equal_tailed_posterior",uncertainty_reason=""))
}

credible_sets <- do.call(rbind,lapply(names(regions)[1:2],function(r)
  do.call(rbind,lapply(c("A","B"),function(t) {
    x <- mapping$fit$sets[[r]][[t]]
    causal <- which(data$simulation$B[regions[[r]],t]!=0)
    stopifnot(length(causal)==1L)
    sets <- x$sets
    data.frame(region=r,trait=t,causal_marker=regions[[r]][causal],
      sets=length(sets),total_members=sum(vapply(sets,function(s)length(s$indices),integer(1))),
      causal_covered=any(vapply(sets,function(s)causal %in% s$indices,logical(1))),
      all_complete=length(sets)>0 && all(vapply(sets,`[[`,logical(1),"complete")),
      all_pass_ld=length(sets)>0 && all(vapply(sets,`[[`,logical(1),"passes_ld_threshold")),
      converged=mapping$fit$diagnostic_summary$converged[
        mapping$fit$diagnostic_summary$region==r & mapping$fit$diagnostic_summary$trait==t])
  }))))
bayesian_recovery <- do.call(rbind,lapply(names(bayesian),function(method) {
  e <- bayesian[[method]]$estimates
  truth <- data$simulation$B[,"A"]*data$dosage_sd[,"A"]/sd(data$simulation$Y[rows$A,"A"])
  causal <- truth!=0
  data.frame(method=method,true_active=sum(causal),expected_active=sum(e$pip),
    effect_correlation=cor(e$mean,truth),effect_rmse=sqrt(mean((e$mean-truth)^2)),
    mean_pip_causal=mean(e$pip[causal]),mean_pip_zero_effect=mean(e$pip[!causal]),
    pip_above_half=sum(e$pip>.5),causal_pip_above_half=sum(e$pip>.5 & causal),
    max_chain_pip_range=max(e$chain_pip_range))
}))
coloc_results <- do.call(rbind,lapply(names(mapping$coloc),function(r) {
  e <- mapping$coloc[[r]]$estimates
  data.frame(region=r,truth=if(r==names(regions)[1]) "H4" else "H3",e[paste0("H",0:4)])
}))
tail_summary <- do.call(rbind,lapply(traits,function(t) {
  e <- tail_comparison[tail_comparison$trait==t,]
  data.frame(trait=t,requested=nrow(e),controlled_available=sum(e$controlled_available),
    approximate_available=sum(e$approximate_available),fallback=sum(e$moment_fallback))
}))
stage_metadata <- lapply(stage_names,function(k) {
  s <- readRDS(file.path(out,paste0(k,".rds"))); s$value <- NULL; s
}); names(stage_metadata) <- stage_names
compact <- list(seed=seed,design=c(individuals=n,markers=m,traits=length(traits)),
  estimation=estimation,prediction=scoring$metrics,association=association_summary,
  marginal_effects=association_effects,bayesian=bayesian_recovery,
  posterior=posterior_table,credible_sets=credible_sets,coloc=coloc_results,
  pathways=pathway_results,availability=availability,tails=tail_summary,timings=timings,
  checks=c(max_ols_error=max(gwas$oracle$max_error),truth_score_error=scoring$truth_score_error),
  provenance=list(stages=stage_metadata,LD_checksums=ld_hash,
    prior=lapply(bayesian,`[[`,"prior"),control=lapply(bayesian,`[[`,"control")))
saveRDS(compact,file.path(out,"compact-results.rds"))
for(k in c("estimation","marginal_effects","bayesian","credible_sets","coloc","tails"))
  write.csv(compact[[k]],file.path(out,paste0(gsub("_","-",k),"-summary.csv")),row.names=FALSE)

draw("association-effects", {
  par(mfrow=c(1,3),mar=c(4.8,4.8,3,1),las=1)
  for(t in traits) {
    x<-marginal_truth[,t]; y<-gwas$stat[[t]]$beta
    lim<-extendrange(range(c(x,y)))
    plot(x,y,pch=16,cex=.3,col=adjustcolor(palette[1],.2),xlim=lim,ylim=lim,
      xlab="True marginal effect (population LD)",ylab="GWAS effect (raw dosage)",
      main=paste("Trait",t)); abline(0,1,lty=2)
  }
},height=700)
draw("regional-annotation", {
  par(mfrow=c(2,2),mar=c(4.8,4.8,3,1),las=1)
  for(method in c("hess_cov","supergnova_cov","lava_cov")) {
    e<-estimation[estimation$analysis==method,]
    lim<-extendrange(range(c(e$truth,e$estimate),finite=TRUE))
    plot(e$truth,e$estimate,pch=19,col=palette[1],xlim=lim,ylim=lim,
      xlab="True regional genetic covariance",ylab="Estimated covariance",
      main=sub("_cov","",toupper(method),ignore.case=TRUE)); abline(0,1,lty=2)
  }
  e<-estimation[estimation$scope=="annotation_category" & estimation$quantity=="h2",]
  labels<-paste(e$trait1,e$component)
  lim<-extendrange(range(c(e$estimate,e$truth,e$lower,e$upper),finite=TRUE))
  plot(seq_len(nrow(e)),e$estimate,xaxt="n",xlab="",ylab="Annotation heritability",
    ylim=lim,pch=19,col=palette[1],main="LDSC annotation categories")
  axis(1,seq_len(nrow(e)),labels,las=2,cex.axis=.65)
  points(seq_len(nrow(e)),e$truth,pch=4,col=palette[2],cex=1.2)
  arrows(seq_len(nrow(e)),e$lower,seq_len(nrow(e)),e$upper,angle=90,code=3,length=.025)
  legend("topright",c("Estimate / approximate 95% interval","Population truth"),
    pch=c(19,4),col=palette[1:2],bty="n",cex=.7)
},height=1250)
draw("pathway-comparison", {
  par(mfrow=c(1,2),mar=c(5,5,3,1),las=1)
  p<-pathway_results[pathway_results$trait=="A" & pathway_results$method!="bayesc",]
  methods<-unique(p$method); sets<-names(pathways$sets)
  plot(1,1,type="n",xlim=c(.7,length(methods)+.3),ylim=extendrange(range(-log10(pmax(p$p_value,1e-300)))),
    xaxt="n",xlab="",ylab="Nominal -log10(p)",main="Pathway tests: trait A")
  axis(1,seq_along(methods),methods,cex.axis=.8)
  for(j in seq_along(sets)) {
    e<-p[p$pathway==sets[j],]; points(match(e$method,methods)+(j-2)*.09,
      -log10(pmax(e$p_value,1e-300)),pch=14+j,col=palette[j])
  }
  legend("topright",sets,pch=15:17,col=palette,bty="n",cex=.75)
  p<-pathway_results[pathway_results$method=="bayesc",]
  barplot(p$pip,names.arg=p$pathway,col=palette,ylim=c(0,1),
    ylab="Posterior inclusion probability",main="Bayesian pathways: trait A")
},height=800)
