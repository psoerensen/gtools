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
# Three continuous traits measured in the same 10K individuals.
traits <- c("height","weight","lipid")
B <- matrix(0,m,3,dimnames=list(ids,traits))
B[8,] <- c(.15,.12,-.1); B[25,1:2] <- c(.2,.15); B[32,3] <- .2
gv3 <- apply(scale(W)%*%B,2,var)
re <- matrix(c(1,.25,-.1,.25,1,.15,-.1,.15,1),3)
sim3 <- gsim::gsim(W=W,architecture="fixed",nt=3,beta=B,h2=gv3/(gv3+1),re=re,
  standardize_W=TRUE,scale_effects=FALSE,seed=20260922,compute_sumstats=FALSE)
Y <- sim3$Y; colnames(Y) <- traits
scale_y3 <- sqrt(colMeans(sweep(Y,2,colMeans(Y))^2))
Ystd <- sweep(sweep(Y,2,colMeans(Y)),2,scale_y3,"/")
stats3 <- setNames(lapply(seq_along(traits),function(t) data.frame(
  marker=ids,allele1="A",allele2="G",beta_std=drop(crossprod(Z,Ystd[,t]))/n,n=n)),traits)
# The simulation supplies the residual covariance. A real analysis needs a
# justified sampling-error covariance; overlap counts alone do not supply it.
Omega <- sim3$Sigma_e/outer(scale_y3,scale_y3)
dimnames(Omega) <- list(traits,traits)
V <- matrix(.01,3,3); diag(V) <- .04; dimnames(V) <- list(traits,traits)
pattern_weights <- c(.9,rep(.1/7,7)) # native order 000,100,010,110,001,101,011,111
joint <- gbayes(stats3,LD,method="bayesr",
  prior=list(effect_covariance=V,pattern_weights=pattern_weights,
    covariance_prior=list(df=6,scale=V*2),pattern_prior=pattern_weights*20,component_prior=c(1,1,1)),
  sampling=list(dependence="shared_ld",covariance=Omega),
  control=list(burnin=1000,sampling_sweeps=2000,seeds=c(31,97),threads=1,
    estimate_covariance=TRUE,estimate_pattern_weights=TRUE,estimate_component_weights=TRUE))
print(joint); print(summary(joint))
joint_scores <- target_z %*% joint$estimates$mean
stopifnot(identical(colnames(joint_scores),traits),all(is.finite(joint_scores)))
# Artificial annotations demonstrate the interface, not biological enrichment.
A <- cbind(category=as.numeric(seq_len(m) %% 3L==0),position=seq_len(m)/m)
rownames(A) <- ids
annotated <- gbayes(stat,LD,method="bayesc",
  prior=list(residual_variance=1/scale_y^2,effect_variance=.04,inclusion_probability=.05),
  annotation=list(inclusion=list(design=A),variance=list(design=A)),control=ctl)
annotated_joint <- gbayes(stats3,LD,method="bayesr",
  prior=list(effect_covariance=V,pattern_weights=pattern_weights),
  sampling=list(dependence="shared_ld",covariance=Omega),
  annotation=list(pattern=list(design=A),component=list(design=A),variance=list(design=A)),
  control=list(burnin=1000,sampling_sweeps=2000,seeds=c(31,97),threads=1))
print(annotated$annotations$inclusion$coefficient_mean)
print(annotated_joint$posterior[grepl("^(pattern|component|variance)\\[",
  annotated_joint$posterior$parameter),])
stopifnot(all(is.finite(annotated_joint$annotations$variance$multiplier_mean)))
# One combined workflow: group priors, reference summaries and selected scores.
groups <- setNames(rep(c("block_1","block_2"),each=m/2),ids)
G <- gs_gprep(bedfiles=paste0(prefix,".bed"))
prediction <- list(Glist=G,ids=G$ids[1:3],target_set="three-example-individuals",
  allele_frequencies=setNames(center/2,ids),
  effect_multipliers=setNames(sqrt(2*(center/2)*(1-center/2))/scale_x,ids),
  retain_draws=TRUE,batch_size=2)
complete <- gbayes(stats3,LD,method="bayesr",
  prior=list(group_covariance=list(block_1=V,block_2=V),covariance_groups=groups,
    group_covariance_prior=list(block_1=list(df=6,scale=V*2),block_2=list(df=6,scale=V*2)),
    pattern_weights=pattern_weights),
  sampling=list(dependence="shared_ld",covariance=Omega),
  annotation=list(variance=list(design=A)),
  posterior=list(reference=LD,partition=groups,phenotype_variance=setNames(rep(1,3),traits),
    prediction=prediction),
  control=list(burnin=1000,sampling_sweeps=2000,seeds=c(31,97),estimate_group_covariance=TRUE))
print(summary(complete))
stopifnot(isTRUE(all.equal(unname(complete$quantities$prediction$mean),
  unname(Z[1:3,]%*%complete$estimates$mean),tolerance=1e-10)))
bundle <- list(complete=complete,annotation=A,annotated=annotated,annotated_joint=annotated_joint,stat=stat,LDlist=LD,fixed=fixed,learned=learned,weights=weights,
  coding=list(center=center,scale=scale_x,phenotype_scale=scale_y),scores=score,
  multivariate=list(stat=stats3,fit=joint,scores=joint_scores,phenotype_scale=scale_y3,
    sampling_error_covariance=Omega,true_effects=sweep(B*sqrt((n-1)/n),2,scale_y3,"/")))
saveRDS(bundle,file.path(out,"workflow.rds"))
stopifnot(identical(readRDS(file.path(out,"workflow.rds")),bundle))
cat("Saved results in",file.path(out,"workflow.rds"),"\n")
