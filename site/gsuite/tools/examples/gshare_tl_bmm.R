# Run from the gsuite checkout after installing its development package.
library(gsuite)
set.seed(2409)
make_site <- function(n,effects) {
  X <- cbind(m1=rnorm(n),m2=rnorm(n),m3=rnorm(n))
  X <- sweep(X,2,colMeans(X))
  y <- drop(X%*%effects+rnorm(n,sd=.7))
  list(X=X,y=y-mean(y))
}
A <- make_site(60,c(.35,-.2,0))
B <- make_site(100,c(.28,-.1,.1))
markers <- data.frame(marker=colnames(A$X),allele="A",center=0,scale=1)
prior <- function(ids) {
  k <- length(ids)
  list(residual_variance=.49,
    weights=setNames(rep(1/k,k),ids),
    variances=setNames(rep(.12,k),ids),
    weight_counts=setNames(rep(1,k),ids),
    variance_priors=list(df=setNames(rep(4,k),ids),
      scale=setNames(rep(.12,k),ids)),
    learn_weights=k>1,learn_variances=TRUE,
    learn_residual=FALSE,residual_prior=NULL)
}
common <- list(markers=markers,trait="simulated trait",units="trait units",
  covariates="centred phenotype and predictors; no fitted covariates",
  control=list(burnin=500,sampling_sweeps=2000,seeds=c(71L,73L)))
source_fit <- do.call(gshare_bmm,c(common,A,
  list(site="A",prior=prior("local"))))
message <- gshare_bmm_source(source_fit)
message_file <- file.path("build","examples","gshare-workflow","tl-bmm-source.txt")
dir.create(dirname(message_file),recursive=TRUE,showWarnings=FALSE)
write_gshare_bmm(message,message_file,overwrite=TRUE)
receiver_fit <- do.call(gshare_bmm,c(common,B,
  list(site="B",prior=prior(c("local","A")),
    sources=list(read_gshare_bmm(message_file)),disjoint=TRUE)))
separate_fit <- do.call(gshare_bmm,c(common,B,
  list(site="B",prior=prior("local"))))
print(source_fit)
print(receiver_fit)
print(data.frame(marker=markers$marker,truth=c(.28,-.1,.1),
  separate=separate_fit$effects$mean,transfer=receiver_fit$effects$mean))
print(receiver_fit$allocation_probability)
# A deliberately conflicting smaller receiver illustrates that transfer
# can be unhelpful. It is a scenario check, not a calibration study.
B_conflict <- make_site(30,c(-.4,.25,0))
conflict_transfer <- do.call(gshare_bmm,c(common,B_conflict,
  list(site="C",prior=prior(c("local","A")),
    sources=list(message),disjoint=TRUE)))
conflict_separate <- do.call(gshare_bmm,c(common,B_conflict,
  list(site="C",prior=prior("local"))))
truth_conflict <- c(-.4,.25,0)
effect_rmse <- function(fit) sqrt(mean((fit$effects$mean-truth_conflict)^2))
print(data.frame(scenario="conflicting effects",
  method=c("separate","transferred"),
  effect_rmse=c(effect_rmse(conflict_separate),
    effect_rmse(conflict_transfer))))
