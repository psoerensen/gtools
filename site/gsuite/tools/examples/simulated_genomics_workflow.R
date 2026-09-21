# One connected qualification example; run from the gsuite repository root.
# Default: regenerate the fixed outputs. Pass --resume to reuse completed stages.
# No package installation, source edits, reference download or website publishing.
.libPaths(c("build/r-library", .libPaths()))
library(gsuite)
stopifnot(requireNamespace("gsim", quietly=TRUE))
out <- "build/examples/simulated-genomics"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
out <- normalizePath(out, winslash="/", mustWork=TRUE)
resume <- "--resume" %in% commandArgs(TRUE)
only_arg <- grep("^--only=",commandArgs(TRUE),value=TRUE)
only <- if(length(only_arg)) strsplit(sub("^--only=","",only_arg),",",fixed=TRUE)[[1]] else NULL
stage_names <- c("data","ld","gwas","gcorr","annotation","gbayes","gscore","gmap","gsea","gbayes-summary")
if(!is.null(only) && any(!only %in% stage_names)) stop("Unknown --only stage")
seed <- 20260921L
n <- 10000L; m <- 50000L; size <- 100L; nr <- m %/% size
traits <- c("A", "B", "C")
ids <- sprintf("m%05d", seq_len(m))
individuals <- sprintf("id%05d", seq_len(n))
region <- rep(sprintf("region%03d", seq_len(nr)), each=size)
regions <- split(ids, region)
rho <- seq(.2, .8, length.out=nr)
rows <- setNames(lapply(0:2, function(k) k*3000L+seq_len(3000L)), traits)
holdout <- 9001:10000
stage <- function(name, expression) {
  path <- file.path(out, paste0(name, ".rds"))
  if(!is.null(only) && !name %in% only) {
    if(file.exists(path)) return(readRDS(path)$value)
    return(NULL)
  }
  if (resume && file.exists(path)) {
    cat("REUSE|",name,"\n",sep=""); flush.console()
    return(readRDS(path)$value)
  }
  cat("START|",name,"\n",sep=""); flush.console()
  start <- Sys.time(); warning_text <- character()
  result <- tryCatch(withCallingHandlers(force(expression), warning=function(w) {
    warning_text <<- c(warning_text, conditionMessage(w))
  }), error=function(e) {
    writeLines(c(name,conditionMessage(e)),file.path(out,"failure.txt"))
    stop(e)
  })
  seconds <- as.numeric(difftime(Sys.time(),start,units="secs"))
  saveRDS(list(value=result, started=start, seconds=seconds,
    warnings=unique(warning_text), script_md5=unname(tools::md5sum(
      "tools/examples/simulated_genomics_workflow.R")), session=sessionInfo(),
    versions=gsuite_versions(),
    native_md5=tools::md5sum(list.files(system.file("libs",package="gsuite"),
      pattern="[.]dll$",recursive=TRUE,full.names=TRUE))),path)
  cat("DONE|",name,"|",round(seconds,2)," seconds\n",sep=""); flush.console()
  result
}

data <- stage("data", {
  set.seed(seed)
  prefix <- file.path(out,"genotypes")
  con <- file(paste0(prefix,".bed"),"wb")
  writeBin(as.raw(c(0x6c,0x1b,0x01)),con)
  # Only one genotype column is generated at a time. The first region is kept
  # for an independent OLS check; never allocate the full N-by-M matrix.
  first_region <- matrix(0L,n,size,dimnames=list(individuals,ids[seq_len(size)]))
  frequencies <- matrix(0,m,3,dimnames=list(ids,traits))
  dosage_sd <- frequencies
  for(j in seq_len(m)) {
    w <- rbinom(n,2,.3)
    r <- (j-1L)%/%size+1L
    if((j-1L)%%size) { keep<-runif(n)<rho[r]; w[keep]<-previous[keep] }
    previous <- w
    if(j<=size) first_region[,j] <- w
    for(t in traits) {
      frequencies[j,t] <- mean(w[rows[[t]]])/2
      dosage_sd[j,t] <- sd(w[rows[[t]]])
    }
    codes <- c(3L,2L,0L)[w+1L]
    writeBin(as.raw(colSums(matrix(codes,4L)*c(1,4,16,64))),con)
    if(j%%10000L==0L) {cat("GENERATED|",j," markers\n",sep="");flush.console()}
  }
  close(con)
  write.table(data.frame(rep(seq_len(nr),each=size),ids,0,
    as.character(rep(seq_len(size)*1000L,nr)),"A","G"),paste0(prefix,".bim"),
    quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table(data.frame(individuals,individuals,0,0,0,-9),paste0(prefix,".fam"),
    quote=FALSE,row.names=FALSE,col.names=FALSE)
  Glist <- gprep(bedfiles=paste0(prefix,".bed"))
  # Enriched first 80 regions; the first two are reserved for simple mapping.
  active <- runif(m)<ifelse(rep(seq_len(nr),each=size)<=80,.08,.008)
  B <- matrix(0,m,3,dimnames=list(ids,traits))
  correlation <- matrix(.6,3,3); diag(correlation)<-1
  B[active,] <- matrix(rnorm(sum(active)*3,sd=.035),sum(active),3)%*%chol(correlation)
  B[1:200,] <- 0
  B[50,] <- .5                 # shared variant in region001
  B[125,1] <- .5; B[175,2:3] <- .5 # distinct A/B signals in region002
  simulation <- gsim::gsim(Glist=Glist,nt=3,architecture="fixed",beta=B,
    h2=c(.3,.4,.35),vg=c(.3,.4,.35),re=0,standardize_W=FALSE,
    scale_effects=TRUE,seed=seed+1L,compute_sumstats=FALSE,chunk_size=64)
  colnames(simulation$Y) <- colnames(simulation$G) <- colnames(simulation$B) <- traits
  # Independent population truth: no estimated reference LD enters this calculation.
  local_G <- lapply(seq_len(nr),function(r) {
    b <- simulation$B[region==names(regions)[r],,drop=FALSE]
    crossprod(b,(.42*rho[r]^abs(outer(seq_len(size),seq_len(size),"-")))%*%b)
  })
  total_G <- Reduce(`+`,local_G)
  total_V <- diag(total_G)+diag(simulation$Sigma_e)
  list(Glist=Glist,simulation=simulation,first_region=first_region,
    frequencies=frequencies,dosage_sd=dosage_sd,local_G=local_G,
    total_G=total_G,total_V=total_V,
    truth_h2=diag(total_G)/total_V,truth_rg=cov2cor(total_G))
})

LD <- stage("ld", {
  reference <- gprep(bedfiles=data$Glist$bedfiles,ids=individuals[1:9000])
  x <- ldprep(reference,reference="simulated-training-pool",assembly="artificial",
    task="sparseld",out_prefix=file.path(out,"LD"),max_distance_bp=0,
    max_distance_variants=size,r2=0,nthreads=1,block_size=128,overwrite=TRUE)
  x <- ldprep(LDlist=x,task="scores")
  for(r in names(regions)[1:20])
    x <- ldprep(LDlist=x,task="eigen",markers=regions[[r]],region_name=r)
  validate_LDlist(x)
  x
})
ld_hash <- tools::md5sum(unlist(LD$resource$paths))

gwas <- stage("gwas", {
  fits <- stat <- setNames(vector("list",3),traits)
  oracle <- list()
  for(t in traits) {
    g <- gprep(bedfiles=data$Glist$bedfiles,ids=individuals[rows[[t]]])
    y <- data$simulation$Y[,t]; names(y)<-individuals
    fits[[t]] <- glma(y,g,method="linear",threads=1,block_size=128)
    a <- fits[[t]]$associations
    stopifnot(nrow(a)==m,all(a$status=="ok"),identical(a$marker,ids))
    stat[[t]] <- data.frame(marker=ids,allele1="A",allele2="G",
      chromosome=LD$resource$markers$chromosome,
      position_bp=LD$resource$markers$base_pair_position,
      beta=a$beta,se=a$se,n=3000,p_value=a$p,z=a$beta/a$se,block=region)
    # Check both an effect marker and background markers against ordinary R OLS.
    for(j in c(1,25,50,75,100)) {
      ref <- summary(lm(y[rows[[t]]]~data$first_region[rows[[t]],j]))$coefficients[2,]
      err <- max(abs(c(a$beta[j],a$se[j],a$p[j])-ref[c(1,2,4)]))
      stopifnot(err<1e-9)
      oracle[[length(oracle)+1]] <- data.frame(trait=t,marker=ids[j],max_error=err)
    }
  }
  list(stat=stat,oracle=do.call(rbind,oracle),
    prepared=sumstat(stat,LD,task="standardize")$stat,
    diagnostics=lapply(fits,`[[`,"diagnostics"))
})

correlation <- stage("gcorr", {
  fits <- list()
  baseline <- matrix(1,m,1,dimnames=list(ids,"baseline"))
  baseline_LD <- ldprep(LDlist=LD,task="scores",annotation=baseline)
  for(method in c("ldsc","gnova","sumher")) for(task in c("h2","rg")) {
    ctrl <- if(task=="rg") list(overlap_intercept=0) else list()
    if(method=="ldsc") ctrl$intercept <- "fixed_one"
    if(method=="sumher") ctrl$tagging <- ld_scores(LD)
    fits[[paste(method,task,sep="_")]] <- gcorr(gwas$stat,baseline_LD,method=method,
      task=task,annotation=if(method=="sumher") baseline else NULL,control=ctrl)
  }
  for(method in c("hess","supergnova","lava")) {
    ctrl <- switch(method,hess=list(shared_sample_size=0,phenotypic_correlation=0,
      minimum_eigenvalue=0,maximum_components=0),supergnova=list(overlap_intercept=0),
      lava=list(sampling_correlation=diag(3),explained_ld_fraction=1))
    if(method=="lava") dimnames(ctrl$sampling_correlation)<-list(traits,traits)
    regional_stat <- lapply(gwas$stat,function(s)s[c("marker","z","n")])
    fits[[paste0(method,"_cov")]] <- gcorr(regional_stat,LD,method=method,task="cov",
      regions=regions[1:20],control=ctrl)
  }
  fits
})

annotation_fit <- stage("annotation", {
  enriched <- as.numeric(rep(seq_len(nr),each=size)<=80)
  annotation <- cbind(enriched=enriched,background=1-enriched)
  rownames(annotation) <- ids
  annotated_LD <- ldprep(LDlist=LD,task="scores",annotation=annotation)
  list(
    ldsc_annotation_h2=gcorr(gwas$stat,annotated_LD,annotation=annotation,
      method="ldsc",task="h2",control=list(intercept="fixed_one")),
    ldsc_annotation_rg=gcorr(gwas$stat,annotated_LD,annotation=annotation,
      method="ldsc",task="rg",control=list(intercept="fixed_one",overlap_intercept=0)))
})

bayesian <- stage("gbayes", {
  fits <- list()
  elapsed <- numeric()
  for(method in c("bayesc","bayesr")) {
    prior <- list(residual_variance=.7,effect_variance=.3/(m*.02))
    if(method=="bayesc") prior$inclusion_probability <- .02 else {
      prior$weights <- c(.98,.006,.007,.007)
      prior$variance_multipliers <- c(0,.01,.1,1)
      prior$effect_variance <- .3/(m*sum(prior$weights*prior$variance_multipliers))
    }
    fit_start <- proc.time()[["elapsed"]]
    fits[[method]] <- gbayes(gwas$prepared$A,LD,method=method,trait="A",prior=prior,
      control=list(burnin=200,sampling_sweeps=1000,seeds=c(11,29),threads=1))
    elapsed[method] <- proc.time()[["elapsed"]]-fit_start
    cat("FIT|",method,"|",round(elapsed[method],2)," seconds\n",sep=""); flush.console()
    stopifnot(all(is.finite(fits[[method]]$estimates$mean)),
      all(fits[[method]]$estimates$pip>=0 & fits[[method]]$estimates$pip<=1))
  }
  attr(fits,"elapsed_seconds") <- elapsed
  fits
})

scoring <- stage("gscore", {
  ridge <- gscore(gwas$prepared$A,LD,method="ridge",control=list(penalty=1))
  weights <- ridge$weights[c("marker","allele1","allele2")]
  # Convert effects per empirical SD to gscore's sqrt(2p(1-p)) dosage coding.
  af <- setNames(data$frequencies[,"A"],ids)
  multiplier <- sqrt(2*af*(1-af))/data$dosage_sd[,"A"]
  weights$ridge <- ridge$weights$score*multiplier
  for(method in names(bayesian)) weights[[method]]<-bayesian[[method]]$estimates$mean*multiplier
  # A truth score is only an execution check, never a fitted candidate.
  weights$truth <- data$simulation$B[,"A"]*sqrt(2*af*(1-af))
  test <- gprep(bedfiles=data$Glist$bedfiles,ids=individuals[holdout])
  prediction <- gscore(task="apply",Glist=test,weights=weights,allele_frequencies=af)
  genetic <- data$simulation$G[holdout,"A"]
  phenotype <- data$simulation$Y[holdout,"A"]
  centered_error <- max(abs(scale(prediction$scores[,"truth"],scale=FALSE)-scale(genetic,scale=FALSE)))
  stopifnot(centered_error<1e-9)
  metrics <- do.call(rbind,lapply(c("ridge","bayesc","bayesr"),function(method) {
    p<-prediction$scores[,method]
    data.frame(method=method,cor_genetic=cor(p,genetic),cor_phenotype=cor(p,phenotype),
      rmse_standardized_genetic=sqrt(mean((scale(p,scale=FALSE)-
        scale(genetic,scale=FALSE)/sd(data$simulation$Y[rows$A,"A"]))^2)))
  }))
  list(ridge=ridge,prediction=prediction,metrics=metrics,truth_score_error=centered_error)
})

mapping <- stage("gmap", {
  vy <- setNames(vapply(traits,function(t)var(data$simulation$Y[rows[[t]],t]),numeric(1)),traits)
  stat <- gmap_stat(gwas$stat,LD,phenotype_variance=vy)
  # Raw dosage priors match the recovered OLS sufficient-statistic scale.
  fit <- gmap(stat[1:2],LD,regions=regions[1:2],method="multi_effect",
    prior=list(residual_variance=1,effect_variance=.04),
    control=list(effects=1,max_sweeps=200))
  stopifnot(all(fit$diagnostic_summary$converged))
  coloc <- setNames(lapply(names(regions)[1:2],function(r) {
    s <- lapply(gwas$stat[1:2],function(s)s[s$marker %in% regions[[r]],])
    gmap(s,task="coloc",method="abf",trait1="A",trait2="B",regions=r,
      phenotype_sd=sqrt(vy[1:2]),
      control=list(complete_region_coverage=TRUE,non_overlapping_samples=TRUE))
  }),names(regions)[1:2])
  list(fit=fit,coloc=coloc)
})

pathways <- stage("gsea", {
  # Two artificial genes per independent region; varying sizes, no real annotation.
  genes <- blocks <- list()
  for(r in seq_len(nr)) {
    cut <- 20L+r%%60L
    gn <- sprintf("gene%04d",2*r-c(1,0))
    genes[[gn[1]]] <- regions[[r]][seq_len(cut)]
    genes[[gn[2]]] <- regions[[r]][(cut+1L):size]
    blocks[[names(regions)[r]]] <- gn
  }
  sets <- list(enriched=names(genes)[1:160],background=names(genes)[401:560],
    mixed=names(genes)[seq(1,999,by=5)])
  metadata <- setNames(lapply(traits,function(t)data.frame(gene=names(genes),sample_size=3000,
    mean_minor_allele_count=vapply(genes,function(g) {
      p<-data$frequencies[match(g,ids),t];mean(6000*pmin(p,1-p))
    },numeric(1)))),traits)
  controlled <- gstat(gwas$stat,LD,sets=genes,blocks=blocks,metadata=metadata,
    control=list(independent_blocks=TRUE,tail_method="controlled_series"))
  saveRDS(controlled,file.path(out,"gene-evidence.rds"))
  # Separate, explicitly approximate route for every gene; never splice p-values
  # into the controlled-series result or remove genes whose bound failed.
  evidence <- gstat(gwas$stat,LD,sets=genes,blocks=blocks,metadata=metadata,
    control=list(independent_blocks=TRUE,tail_method="saddlepoint"))
  saveRDS(evidence,file.path(out,"gene-evidence-saddlepoint.rds"))
  stopifnot(all(vapply(evidence$stat,function(s)all(s$p_available),logical(1))))
  fits <- list(
    ora=gsea(lapply(evidence$stat,function(s)transform(s,selected=p_value<.05)),sets,"ora"),
    preranked=gsea(evidence$stat,sets,"preranked",control=list(seed=31,replicates=999)),
    competitive=gsea(evidence$stat,sets,"competitive",sampling=evidence$sampling),
    magma=gsea(evidence$stat,sets,"magma",sampling=evidence$sampling),
    bayesc=gsea(evidence$stat["A"],sets,"bayesc",sampling=evidence$sampling,
      prior=list(residual_variance=1,effect_variance=1),
      control=list(burnin=200,sampling_sweeps=1000,seeds=c(11,29))))
  list(evidence=evidence,controlled=controlled,fits=fits,genes=genes,sets=sets)
})

# Capture joint-draw genetic variance without changing the original samplers or
# overwriting the original timing evidence. Fixed hyperparameters are not learned.
bayesian_summary <- stage("gbayes-summary", {
  fits <- lapply(names(bayesian),function(method) {
    original <- bayesian[[method]]
    fit <- gbayes(gwas$prepared$A,LD,method=method,trait="A",
      prior=original$prior,control=original$control,
      posterior=list(reference=LD,phenotype_variance=c(A=1)))
    stopifnot(identical(fit$estimates,original$estimates),
      identical(fit$parameter_mean,original$parameter_mean))
    fit
  })
  setNames(fits,names(bayesian))
})

if(!is.null(only)) {
  cat("SELECTED_STAGES|completed; full report requires all stages\n")
} else {
stopifnot(identical(ld_hash,tools::md5sum(unlist(LD$resource$paths))))
correlation <- c(correlation,annotation_fit)
estimates <- do.call(rbind,lapply(names(correlation),function(k)
  cbind(analysis=k,correlation[[k]]$estimates)))
write.csv(estimates,file.path(out,"correlation.csv"),row.names=FALSE)
write.csv(scoring$metrics,file.path(out,"prediction.csv"),row.names=FALSE)
write.csv(gwas$oracle,file.path(out,"ols-check.csv"),row.names=FALSE)
write.csv(mapping$fit$estimates,file.path(out,"mapping.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(mapping$coloc,`[[`,"estimates")),file.path(out,"coloc.csv"),row.names=FALSE)
write.csv(data.frame(trait=traits,h2=data$truth_h2),file.path(out,"truth.csv"),row.names=FALSE)
write.csv(data$truth_rg,file.path(out,"truth-rg.csv"))
pathway_results <- do.call(rbind,lapply(names(pathways$fits),function(method) {
  p<-pathways$fits[[method]]$pathways
  data.frame(method=method,trait=p$trait,pathway=p$pathway,
    p_value=if("p_value" %in% names(p)) p$p_value else NA_real_,
    pip=if("pip" %in% names(p)) p$pip else NA_real_)
}))
frequentist <- pathway_results$method!="bayesc"
stopifnot(all(is.finite(pathway_results$p_value[frequentist])),
  all(pathway_results$p_value[frequentist]>=0 & pathway_results$p_value[frequentist]<=1),
  all(is.finite(pathway_results$pip[!frequentist])),
  all(pathway_results$pip[!frequentist]>=0 & pathway_results$pip[!frequentist]<=1),
  all(pathways$fits$preranked$pathways$inference_available))
write.csv(pathway_results,file.path(out,"pathways.csv"),row.names=FALSE)
tail_comparison <- do.call(rbind,lapply(traits,function(t) {
  exact<-pathways$controlled$stat[[t]]; approx<-pathways$evidence$stat[[t]]
  stopifnot(identical(exact$gene,approx$gene),
    identical(exact$statistic,approx$statistic))
  data.frame(trait=t,gene=exact$gene,statistic=exact$statistic,
    controlled_available=exact$p_available,controlled_p=exact$p_value,
    controlled_error_bound=exact$tail_error_bound,controlled_terms=exact$tail_terms,
    controlled_reason=exact$unavailable_reason,
    approximate_available=approx$p_available,approximate_p=approx$p_value,
    approximate_method=approx$tail_method,moment_fallback=approx$moment_fallback)
}))
write.csv(tail_comparison,file.path(out,"gene-tail-comparison.csv"),row.names=FALSE)
stages <- stage_names
timings <- do.call(rbind,lapply(stages,function(k) {
  s<-readRDS(file.path(out,paste0(k,".rds")))
  data.frame(stage=k,seconds=s$seconds,started=as.character(s$started),script_md5=s$script_md5,
    warnings=paste(s$warnings,collapse="; "))
}))
write.csv(timings,file.path(out,"timings.csv"),row.names=FALSE)
availability <- do.call(rbind,lapply(names(correlation),function(k) {
  e<-correlation[[k]]$estimates
  data.frame(analysis=k,estimates=nrow(e),unavailable=sum(!is.finite(e$estimate)),
    out_of_range=sum(e$out_of_range,na.rm=TRUE),available_se=sum(e$se_available),
    uncertainty=paste(unique(e$uncertainty_reason[nzchar(e$uncertainty_reason)]),collapse="; "))
}))
write.csv(availability,file.path(out,"availability.csv"),row.names=FALSE)
chain_summary <- do.call(rbind,lapply(names(bayesian),function(k) {
  e<-bayesian[[k]]$estimates
  data.frame(method=k,max_pip_range=max(e$chain_pip_range),
    markers_pip_range_over_point1=sum(e$chain_pip_range>.1),
    effect_correlation=cor(e$mean,data$simulation$B[,"A"]*data$dosage_sd[,"A"]))
}))
write.csv(chain_summary,file.path(out,"chains.csv"),row.names=FALSE)

# LD can make a marker with zero direct effect significant as well.
association_summary <- do.call(rbind,lapply(traits,function(t) {
  s <- gwas$stat[[t]]; causal <- data$simulation$B[,t]!=0
  detected <- s$p_value < .05/m
  stopifnot(identical(s$marker,rownames(data$simulation$B)),
    all(is.finite(s$p_value)),all(s$p_value>=0 & s$p_value<=1))
  data.frame(trait=t,markers=m,causal_markers=sum(causal),
    threshold=.05/m,significant_markers=sum(detected),
    significant_causal=sum(detected & causal),
    significant_zero_effect=sum(detected & !causal))
}))
write.csv(association_summary,file.path(out,"association-summary.csv"),row.names=FALSE)
posterior_table <- do.call(rbind,lapply(names(bayesian_summary),function(method) {
  fit <- bayesian_summary[[method]]
  stopifnot(identical(fit$estimates,bayesian[[method]]$estimates))
  tab <- fit$posterior; tab$method <- method; tab$truth <- NA_real_
  tab$truth[tab$parameter=="h2[A]"] <- data$truth_h2["A"]
  tab$truth[tab$parameter=="active_markers"] <- sum(data$simulation$B[,"A"]!=0)
  tab
}))
write.csv(posterior_table,file.path(out,"bayesian-posterior.csv"),row.names=FALSE)
capture.output(for(method in names(bayesian_summary)) {
  cat("\n",toupper(method),"; fixed hyperparameters; conditional posterior summaries\n")
  print(summary(bayesian_summary[[method]]))
},file=file.path(out,"bayesian-fit-summary.txt"))
draw <- function(name, code, width=1800,height=1100) {
  png(file.path(out,paste0(name,".png")),width=width,height=height,res=150)
  tryCatch(force(code),finally=dev.off())
}
palette <- c("#28678b","#bf6b2d","#38836c")
causal_colour <- "#c55a11"
region_index <- rep(seq_len(nr),each=size)
marker_x <- (seq_len(m)-.5)/size+.5
draw("association-manhattan", {
  par(mfrow=c(3,1),mar=c(3.8,4.5,2.8,1),oma=c(0,0,1,0),las=1)
  for(t in traits) {
    s <- gwas$stat[[t]]; causal <- data$simulation$B[,t]!=0
    y <- -log10(pmax(s$p_value,.Machine$double.xmin))
    plot(marker_x,y,type="n",xlim=c(.5,nr+.5),ylim=c(0,max(y)*1.13),
      xlab="Artificial region (100 markers each)",ylab="-log10(p)",
      main=paste("Linear GWAS: trait",t))
    rect(.5,-1,80.5,max(y)*1.2,col="#edf4fa",border=NA)
    points(marker_x[!causal],y[!causal],pch=16,cex=.32,
      col=ifelse(region_index[!causal]%%2,"#64748b","#aab5c3"))
    points(marker_x[causal],y[causal],pch=17,cex=.48,col=causal_colour)
    abline(h=-log10(.05/m),lty=2,col="#333333")
    legend("topright",c("True causal marker","Zero direct effect (may tag a causal marker)",
      "0.05 / 50,000 per trait"),pch=c(17,16,NA),lty=c(NA,NA,2),
      col=c(causal_colour,"#64748b","#333333"),bty="n",cex=.72)
  }
},height=1500)

draw("bayesian-markers", {
  par(mfrow=c(2,2),mar=c(4.5,4.8,3,1),las=1)
  truth_effect <- data$simulation$B[,"A"]*data$dosage_sd[,"A"]/
    sd(data$simulation$Y[rows$A,"A"])
  causal <- truth_effect!=0
  for(method in names(bayesian_summary)) {
    e <- bayesian_summary[[method]]$estimates
    plot(marker_x,e$pip,type="n",ylim=c(0,1.12),
      xlab="Artificial region",ylab="Posterior inclusion probability",
      main=paste(toupper(method),"- trait A"))
    points(marker_x[!causal],e$pip[!causal],pch=16,cex=.3,col="#aab5c3")
    points(marker_x[causal],e$pip[causal],pch=17,cex=.5,col=causal_colour)
    legend("topright",c("True causal","Zero direct effect"),pch=c(17,16),
      col=c(causal_colour,"#aab5c3"),bty="n",cex=.8)
    limits <- extendrange(range(c(truth_effect,e$mean)))
    plot(truth_effect,e$mean,type="n",xlim=limits,ylim=limits,
      xlab="True standardized effect",ylab="Posterior mean effect",
      main="Marker-effect recovery")
    abline(0,1,lty=2,col="#555555")
    points(truth_effect[!causal],e$mean[!causal],pch=16,cex=.35,col="#aab5c3")
    points(truth_effect[causal],e$mean[causal],pch=17,cex=.55,col=causal_colour)
  }
},height=1250)

draw("bayesian-heritability", {
  par(mfrow=c(2,2),mar=c(4.5,4.8,3,1),las=1)
  for(method in names(bayesian_summary)) {
    fit <- bayesian_summary[[method]]
    traces <- lapply(fit$quantities$chains,function(chain)
      as.numeric(chain$genetic$covariance_trace[,1]))
    h <- fit$posterior[fit$posterior$parameter=="h2[A]",]
    stopifnot(nrow(h)==1L,all(is.finite(unlist(traces))),
      abs(mean(unlist(traces))-h$mean)<1e-10)
    yrange <- extendrange(range(c(unlist(traces),data$truth_h2["A"])))
    plot(seq_along(traces[[1]]),traces[[1]],type="n",ylim=yrange,
      xlab="Retained draw within chain",ylab="Reference heritability",
      main=paste(toupper(method),"- two chains"))
    for(k in seq_along(traces)) lines(traces[[k]],col=adjustcolor(palette[k],.65))
    abline(h=data$truth_h2["A"],lty=2)
    legend("topright",c("Chain 1","Chain 2","Population truth"),
      col=c(palette[1:2],"black"),lty=c(1,1,2),bty="n",cex=.8)
    hist(unlist(traces),breaks=35,col="#c6dbea",border="white",
      xlim=yrange,xlab="Reference heritability",main="Conditional posterior",ylab="Draw count")
    abline(v=data$truth_h2["A"],lty=2,lwd=2)
    abline(v=h$mean,col=palette[1],lwd=2)
    abline(v=c(h$lower,h$upper),col=palette[1],lty=3)
    legend("topright",c("Population truth","Posterior mean","95% interval"),
      col=c("black",palette[1],palette[1]),lty=c(2,1,3),bty="n",cex=.8)
  }
},height=1250)

draw("bayesian-prediction", {
  par(mfrow=c(1,2),mar=c(4.8,4.8,3,1),las=1)
  truth_score <- as.numeric(scale(data$simulation$G[holdout,"A"],scale=FALSE))/
    sd(data$simulation$Y[rows$A,"A"])
  for(method in names(bayesian_summary)) {
    p <- as.numeric(scale(scoring$prediction$scores[,method],scale=FALSE))
    lim <- extendrange(range(c(truth_score,p)))
    plot(truth_score,p,pch=16,cex=.55,col=adjustcolor(palette[1],.45),
      xlim=lim,ylim=lim,xlab="True genetic value / training phenotype SD",
      ylab="Predicted genetic value (centered)",
      main=sprintf("%s: r = %.3f",toupper(method),cor(p,truth_score)))
    abline(0,1,lty=2,col="#555555")
  }
},height=800)

draw("heritability", {
  par(mfrow=c(2,2),mar=c(5,5,4,2),las=1)
  for(t in traits) {
    e<-estimates[estimates$scope=="reference" & estimates$quantity=="h2" & estimates$trait1==t,]
    e<-e[match(c("ldsc_h2","gnova_h2","sumher_h2","ldsc_annotation_h2"),e$analysis),]
    lo<-e$estimate-1.96*e$se;hi<-e$estimate+1.96*e$se
    truth<-data$truth_h2[match(t,traits)]
    cols<-c(palette,"#8661a1")
    plot(1:4,e$estimate,xaxt="n",xlab="",ylab="Heritability",pch=19,col=cols,
      xlim=c(.7,4.3),ylim=extendrange(range(c(lo,hi,truth),finite=TRUE)),main=paste("Trait",t))
    axis(1,1:4,c("LDSC","GNOVA","SumHer","LDSC+ann"),cex.axis=.8);abline(h=truth,lty=2)
    ok<-is.finite(lo)&is.finite(hi)
    arrows(which(ok),lo[ok],which(ok),hi[ok],angle=90,code=3,length=.04,col=cols[ok])
  }
  e<-estimates[estimates$scope=="reference" & estimates$quantity=="rg",]
  target<-mapply(function(a,b)data$truth_rg[a,b],e$trait1,e$trait2)
  cols<-ifelse(grepl("annotation",e$analysis),"#8661a1",palette[match(e$method,c("ldsc","gnova","sumher"))])
  plot(target,e$estimate,pch=19,col=cols,
    xlab="True genetic correlation",ylab="Estimated genetic correlation",
    xlim=extendrange(range(c(target,e$estimate),finite=TRUE)),
    ylim=extendrange(range(c(target,e$estimate),finite=TRUE)),main="Three trait pairs")
  abline(0,1,lty=2);legend("topleft",c("LDSC","GNOVA","SumHer","LDSC+ann"),pch=19,col=c(palette,"#8661a1"),bty="n",cex=.75)
})
draw("mapping-prediction", {
  par(mfrow=c(2,2),mar=c(5,5,4,2),las=1)
  for(r in names(regions)[1:2]) {
    e<-mapping$fit$estimates
    a<-e[e$region==r & e$trait=="A",];b<-e[e$region==r & e$trait=="B",]
    plot(seq_len(size),a$pip,type="h",col=palette[1],lwd=2,ylim=c(0,1),xaxt="n",
      xlab="Marker within region",ylab="Posterior inclusion probability",main=r)
    axis(1,c(1,25,50,75,100))
    lines(seq_len(size)+.2,b$pip,type="h",col=palette[2],lwd=2)
    causal<-which(rowSums(abs(data$simulation$B[regions[[r]],1:2,drop=FALSE]))>0)
    abline(v=causal,lty=3,col="#777777")
    legend("topright",c("Trait A","Trait B","Causal position"),col=c(palette[1:2],"#777777"),lty=c(1,1,3),bty="n",cex=.8)
  }
  probabilities<-sapply(mapping$coloc,function(f)as.numeric(f$estimates[1,paste0("H",0:4)]))
  barplot(probabilities,beside=TRUE,names.arg=c("Shared signal","Distinct signals"),
    col=c("#cccccc","#95b8d1","#e3b793",palette[1],palette[2]),ylim=c(0,1.18),yaxt="n",
    ylab="Posterior probability",main="Regional ABF colocalisation")
  axis(2,seq(0,1,.2))
  legend("top",paste0("H",0:4),fill=c("#cccccc","#95b8d1","#e3b793",palette[1:2]),bty="n",cex=.7,horiz=TRUE)
  barplot(scoring$metrics$cor_genetic,names.arg=scoring$metrics$method,col=palette,
    ylim=c(0,1),ylab="Correlation with true genetic value",main="1,000 held-out individuals")
})
draw("pathways", {
  par(mfrow=c(2,2),mar=c(5,5,4,2),las=1)
  for(t in traits) {
    p<-pathway_results[pathway_results$trait==t & pathway_results$method=="magma",]
    barplot(-log10(pmax(p$p_value,.Machine$double.xmin)),names.arg=p$pathway,
      col=palette,ylab="Pathway -log10(p)",main=paste("Bounded MAGMA: trait",t))
  }
  p<-pathway_results[pathway_results$method=="bayesc",]
  barplot(p$pip,names.arg=p$pathway,ylim=c(0,1),col=palette,yaxt="n",
    ylab="Posterior inclusion probability",main="Bayesian pathway model: trait A")
  axis(2,seq(0,1,.2))
})
print(timings[c("stage","seconds")]);print(scoring$metrics)
print(estimates[estimates$scope=="reference",c("method","quantity","trait1","trait2","estimate","se")])
print(pathway_results)
if(file.exists(file.path(out,"failure.txt"))) unlink(file.path(out,"failure.txt"))
cat("SIMULATED_GENOMICS|completed\n")
}
