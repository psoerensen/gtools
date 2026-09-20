# One reproducible summary-statistic regression example; installed gsuite + gsim.
library(gsuite)
stopifnot(requireNamespace("gsim",quietly=TRUE))
out <- getOption("gsuite.bayes.out","build/examples/gbayes-workflow")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,mustWork=TRUE)
set.seed(20260920)
n <- 10000L; m <- 40L
ids <- sprintf("m%02d",seq_len(m))
W <- matrix(rbinom(n*m,2,.3),n,m,dimnames=list(NULL,ids))
for(j in 2:m) {
  copy <- runif(n)<.5
  W[copy,j] <- W[copy,j-1L]
}
b <- setNames(rep(0,m),ids); b[c(8,25)] <- c(.15,.2)
gv <- var(drop(scale(W)%*%b))
sim <- gsim::gsim(W=W,architecture="fixed",nt=1,beta=b,h2=gv/(gv+1),
  standardize_W=TRUE,scale_effects=FALSE,seed=20260921,compute_sumstats=FALSE)
y <- drop(sim$Y)
# Population SD convention gives diagonal X'X=N, as used by the native adapter.
center <- colMeans(W)
scale_x <- sqrt(colMeans(sweep(W,2,center)^2))
Z <- sweep(sweep(W,2,center),2,scale_x,"/")
scale_y <- sqrt(mean((y-mean(y))^2)); zy <- (y-mean(y))/scale_y

# A small full-LD reference from these same individuals; A is the counted allele.
prefix <- file.path(out,"reference")
con <- file(paste0(prefix,".bed"),"wb")
writeBin(as.raw(c(0x6c,0x1b,0x01)),con)
for(j in seq_len(m)) {
  codes <- c(3L,2L,0L)[W[,j]+1L]
  writeBin(as.raw(colSums(matrix(codes,4L)*c(1,4,16,64))),con)
}
close(con)
write.table(data.frame(1,ids,0,seq_len(m)*1000,"A","G"),
  paste0(prefix,".bim"),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(data.frame(1:n,1:n,0,0,0,-9),paste0(prefix,".fam"),
  quote=FALSE,row.names=FALSE,col.names=FALSE)
LD <- ldprep(gs_gprep(bedfiles=paste0(prefix,".bed")),reference="artificial-training-panel",
  assembly="artificial",task="sparseld",out_prefix=file.path(out,"LD"),
  max_distance_bp=0,max_distance_variants=m,r2=0,nthreads=1,overwrite=TRUE)
saveRDS(LD,file.path(out,"LDlist.rds")); LD <- readRDS(file.path(out,"LDlist.rds"))
stat <- data.frame(marker=ids,allele1="A",allele2="G",
  beta_std=drop(crossprod(Z,zy))/n,n=n)
ctl <- list(burnin=1000,sampling_sweeps=4000,seeds=c(31,97),threads=1)
fixed <- gbayes(stat,LD,method="bayesr",trait="simulated",
  prior=list(residual_variance=1/scale_y^2,effect_variance=.04,
    weights=c(.9,.04,.04,.02)),control=ctl)
learned <- gbayes(stat,LD,method="bayesc",trait="simulated",
  prior=list(residual_variance=1/scale_y^2,effect_variance=.04,inclusion_probability=.05,
    effect_variance_prior=list(df=4,scale=.04),weight_prior=c(38,2)),
  control=c(ctl,list(estimate_effect_variance=TRUE,estimate_weights=TRUE)))
print(fixed); print(summary(learned))

# Apply posterior-mean effects to another dosage matrix on the SAME scale.
# This illustrates weight reuse; it is not an assessment of prediction accuracy.
target <- matrix(rbinom(10*m,2,.3),10,m,dimnames=list(NULL,ids))
weights <- fixed$estimates[,c("marker","allele1","allele2","mean")]
stopifnot(identical(colnames(target),weights$marker)) # counted allele A throughout
target_z <- sweep(sweep(target,2,center),2,scale_x,"/")
score <- drop(target_z %*% weights$mean) # standardized phenotype units
stopifnot(all(is.finite(score)),all(fixed$estimates$pip>=0 & fixed$estimates$pip<=1))
bundle <- list(stat=stat,LDlist=LD,fixed=fixed,learned=learned,weights=weights,
  coding=list(center=center,scale=scale_x,phenotype_scale=scale_y),scores=score)
saveRDS(bundle,file.path(out,"workflow.rds"))
stopifnot(identical(readRDS(file.path(out,"workflow.rds")),bundle))
cat("Saved results in",file.path(out,"workflow.rds"),"\n")
