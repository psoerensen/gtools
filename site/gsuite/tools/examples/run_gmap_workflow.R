# Run with installed gsuite and gsim. All outputs stay in one directory.
library(gsuite)
stopifnot(requireNamespace("gsim",quietly=TRUE))
out <- getOption("gsuite.mapping.out","build/examples/gmap-workflow")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,mustWork=TRUE)
started <- Sys.time()
set.seed(20260915)
n <- 10000L; nref <- 4000L; size <- 20L; m <- 2L*size
ids <- sprintf("m%02d",seq_len(m))
regions <- list(shared=ids[1:20],distinct=ids[21:40])
genotypes <- function(count) {
  W <- matrix(0L,count,m,dimnames=list(NULL,ids))
  for(j in seq_len(m)) {
    W[,j] <- rbinom(count,2,.3)
    if ((j-1L) %% size) {
      copy <- runif(count)<.8
      W[copy,j] <- W[copy,j-1L]
    }
  }
  W
}
ref <- genotypes(nref); W1 <- genotypes(n); W2 <- genotypes(n)
# A small artificial reference; A is the counted allele in these BED files.
prefix <- file.path(out,"reference")
con <- file(paste0(prefix,".bed"),"wb")
writeBin(as.raw(c(0x6c,0x1b,0x01)),con)
for(j in seq_len(m)) {
  codes <- c(3L,2L,0L)[ref[,j]+1L]
  writeBin(as.raw(colSums(matrix(codes,4L)*c(1,4,16,64))),con)
}
close(con)
write.table(data.frame(rep(1:2,each=size),ids,0,rep(seq_len(size)*1000,2),"A","G"),
  paste0(prefix,".bim"),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(data.frame(1:nref,1:nref,0,0,0,-9),paste0(prefix,".fam"),
  quote=FALSE,row.names=FALSE,col.names=FALSE)
LD <- ldprep(gs_gprep(bedfiles=paste0(prefix,".bed")),reference="artificial-panel",
  assembly="artificial",task="sparseld",out_prefix=file.path(out,"LD"),
  max_distance_bp=0,max_distance_variants=size,r2=0,nthreads=1,overwrite=TRUE)
for(r in names(regions)) LD <- ldprep(LDlist=LD,task="eigen",markers=regions[[r]],region_name=r)
saveRDS(LD,file.path(out,"LDlist.rds")); LD <- readRDS(file.path(out,"LDlist.rds"))
checksums_before <- tools::md5sum(unlist(LD$resource$paths))
make_trait <- function(W,causal,seed) {
  b <- setNames(rep(0,m),ids); b[causal] <- c(.12,.14)
  gv <- var(drop(scale(W)%*%b))
  sim <- gsim::gsim(W=W,architecture="fixed",nt=1,beta=b,h2=gv/(gv+1),
    standardize_W=TRUE,scale_effects=FALSE,seed=seed,compute_sumstats=FALSE)
  stopifnot(max(abs(drop(scale(W)%*%b)-drop(sim$G)))<1e-10,
            abs(sim$Sigma_e[1,1]-1)<1e-10)
  y <- drop(sim$Y); X <- scale(W,center=TRUE,scale=FALSE)
  yc <- y-mean(y); d <- colSums(X^2); score <- drop(crossprod(X,yc))
  beta <- score/d; se <- sqrt((sum(yc^2)-score^2/d)/((n-2)*d))
  md <- LD$resource$markers
  stat <- data.frame(marker=md$marker_id,allele1=md$allele1,allele2=md$allele2,beta=beta,se=se,n=n)
  exact <- stat[1:3]; exact$score <- score; exact$marker_ss <- d
  list(stat=stat,variance=var(y),exact=exact,truth=b/apply(W,2,sd))
}
A <- make_trait(W1,c(8,25),20260916)
B <- make_trait(W2,c(8,36),20260917)
gwas <- list(A=A$stat,B=B$stat)
stat <- gmap_stat(gwas,LD,phenotype_variance=c(A=A$variance,B=B$variance))
exact <- list(A=A$exact,B=B$exact)
conversion <- do.call(rbind,lapply(names(stat),function(t) data.frame(trait=t,
  max_score_error=max(abs(stat[[t]]$score-exact[[t]]$score)),
  max_marker_ss_error=max(abs(stat[[t]]$marker_ss-exact[[t]]$marker_ss)))))
stopifnot(max(conversion$max_score_error)<1e-8,max(conversion$max_marker_ss_error)<1e-8)
prior <- list(residual_variance=1,effect_variance=.04)
mc <- list(burnin=1000,sampling_sweeps=5000,seeds=c(31,97),retain_diagnostics=TRUE)
fits <- list(); timing <- list(); equivalence <- list()
for(method in c("bayesc","bayesr","multi_effect")) {
  p <- prior
  if(method=="bayesc") p$inclusion_probability <- .05
  if(method=="bayesr") p$weights <- c(.95,.01,.01,.03)
  ctl <- if(method=="multi_effect") list(effects=3,max_sweeps=500,
    estimate_effect_variance=TRUE) else mc
  elapsed <- system.time(fits[[method]] <- gmap(stat,LD,regions,method,p,ctl))[["elapsed"]]
  timing[[method]] <- data.frame(method=method,seconds=elapsed)
  direct <- gmap(exact,LD,regions,method,p,ctl)
  equivalence[[method]] <- data.frame(method=method,
    max_pip_error=max(abs(fits[[method]]$estimates$pip-direct$estimates$pip)),
    max_mean_error=max(abs(fits[[method]]$estimates$mean-direct$estimates$mean)))
}
equivalence <- do.call(rbind,equivalence)
stopifnot(max(equivalence$max_pip_error)<1e-8,max(equivalence$max_mean_error)<1e-8,
  all(fits$multi_effect$diagnostic_summary$converged))
shared <- coloc(fits$multi_effect,fits$multi_effect,trait1="A",trait2="B",
  control=list(complete_region_coverage=TRUE,non_overlapping_samples=TRUE))
# Single-causal ABF route uses the regional GWAS summaries directly, without LD.
abf <- lapply(names(regions),function(region) {
  ix <- match(regions[[region]],gwas$A$marker)
  coloc(gwas$A[ix,],gwas$B[ix,],method="abf",trait1="A",trait2="B",regions=region,
    phenotype_sd=sqrt(c(A=A$variance,B=B$variance)),
    control=list(complete_region_coverage=TRUE,non_overlapping_samples=TRUE),
    sensitivity=list(lower_shared=list(p1=1e-4,p2=1e-4,p12=1e-6)))
})
names(abf) <- names(regions)
write.csv(do.call(rbind,lapply(abf,`[[`,"estimates")),file.path(out,"abf.csv"),row.names=FALSE)
stopifnot(identical(checksums_before,tools::md5sum(unlist(LD$resource$paths))))
truth <- data.frame(marker=ids,A=A$truth,B=B$truth)
saveRDS(list(gwas=gwas,phenotype_variance=c(A=A$variance,B=B$variance),stat=stat,
  regions=regions,truth=truth,fits=fits,coloc=shared,abf=abf,conversion=conversion,
  equivalence=equivalence,timing=do.call(rbind,timing)),file.path(out,"workflow.rds"))
write.csv(do.call(rbind,timing),file.path(out,"timing.csv"),row.names=FALSE)
write.csv(conversion,file.path(out,"conversion.csv"),row.names=FALSE)
write.csv(equivalence,file.path(out,"equivalence.csv"),row.names=FALSE)
write.csv(shared$estimates,file.path(out,"coloc.csv"),row.names=FALSE)
summaries <- do.call(rbind,lapply(names(fits),function(method) {
  z <- summary(fits[[method]])$estimates; z$method <- method; z
}))
write.csv(summaries,file.path(out,"summaries.csv"),row.names=FALSE)
png(file.path(out,"mapping.png"),width=1800,height=1200,res=150)
par(mfrow=c(2,2),mar=c(5,4.5,3,1))
colors <- c("#235789","#d97922","#27826a")
for(trait in c("A","B")) for(r in names(regions)) {
  pips <- sapply(fits,function(f) f$estimates$pip[f$estimates$region==r & f$estimates$trait==trait])
  matplot(seq_len(size),pips,type="b",pch=c(16,17,15),lty=1,col=colors,
    ylim=c(0,1),xlab="Marker within region",ylab="Posterior inclusion probability",
    main=paste(r,"region / trait",trait))
  causal <- which(truth[[trait]][match(regions[[r]],ids)]!=0)
  abline(v=causal,lty=3,col="grey35")
  legend("topright",c("BayesC","BayesR","Multi-effect"),col=colors,pch=c(16,17,15),lty=1,bty="n",cex=.75)
}
dev.off()
png(file.path(out,"coloc.png"),width=1600,height=700,res=140)
par(mfrow=c(1,2),mar=c(5,4.5,3,1))
for(r in names(regions)) plot(shared,region=r)
dev.off()
capture.output({print(conversion);print(equivalence);print(do.call(rbind,timing));
  for(f in fits) print(summary(f)); print(summary(shared)); print(sessionInfo())},
  file=file.path(out,"report.txt"))
cat("Workflow completed in",round(as.numeric(difftime(Sys.time(),started,units="secs")),2),"seconds.\n")
cat("Results:",file.path(out,"workflow.rds"),"\n")
