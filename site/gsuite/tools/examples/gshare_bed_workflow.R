# Run from gsuite after the focused native and package builds.
# Uses the existing 100-marker gsim example; all outputs stay under build/.
.libPaths(c("build/r-library",.libPaths()))
library(gsuite)
args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args)) args[[1]] else "fit"
stopifnot(mode %in% c("fit","wrong","report"))
base <- "build/examples/gshare-workflow"
if(!file.exists(file.path(base,"simulation.rds")))
  stop("Run tools/examples/gshare_workflow.R once to prepare the gsim dataset")
out <- "build/examples/gshare-bed"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
sim <- readRDS(file.path(base,"simulation.rds"))
p <- ncol(sim$W)
stopifnot(p==100L,all(colnames(sim$W)==paste0("m",seq_len(p))))
bed_path <- file.path(base,"genotypes.bed")
stopifnot(file.exists(bed_path))
frequency <- setNames(sim$reference_frequency,colnames(sim$W))
markers <- data.frame(marker=names(frequency),allele="G",center=2*unname(frequency),
  scale=sqrt(2*unname(frequency)*(1-unname(frequency))))
ids <- rownames(sim$W)
site_ids <- lapply(sim$rows[c("A","B")],function(i) ids[i])
settings <- list(trait="site phenotype",units="simulation units",
  covariates="known zero intercept; donor-reference centering/scaling",
  control=list(burnin=2000L,sampling_sweeps=8000L,seeds=c(2027L,2029L),threads=1L))
source_prior <- list(residual_variance=1,effect_variance=.04,inclusion_probability=.1)
relationship <- list(heterogeneous=c(.7,.0204),similar=c(.95,.0039),
  wrong=c(-1,.0001),receiver_only=c(.7,.0204))
errors <- sim$simulation$Y - sim$X%*%sim$B
effects <- list(heterogeneous=sim$B,similar=sim$B,wrong=sim$B,receiver_only=sim$B)
effects$similar[,2] <- .95*sim$B[,1]+.15*(sim$B[,2]-.7*sim$B[,1])
extra <- setdiff(seq_len(p),sim$causal)[1:2]
effects$receiver_only[extra,2] <- c(.2,-.2)
Y <- lapply(effects,function(B) sim$X%*%B+errors)
bed <- lapply(site_ids,function(selected) gprep(bedfiles=bed_path,ids=selected))
fit_one <- function(site,response,prior,transfer=NULL,packed=TRUE) {
  i <- sim$rows[[site]]
  core <- list(site=site,prior=prior,transfer=transfer,disjoint=!is.null(transfer))
  if(packed) do.call(gshare,c(settings,core,list(Glist=bed[[site]],
    data=data.frame(id=site_ids[[site]],`site phenotype`=response[i],check.names=FALSE),
    reference_frequency=frequency))) else do.call(gshare,c(settings,core,list(
      X=sim$X[i,,drop=FALSE],y=response[i],markers=markers)))
}
if(mode %in% c("fit","wrong")) {
  saved <- if(mode=="wrong") readRDS(file.path(out,"fits.rds")) else list()
  source <- if(mode=="wrong") saved$source else fit_one("A",Y$heterogeneous[,1],source_prior)
  saved$source <- source
  metrics <- list()
  for(name in if(mode=="wrong") "wrong" else names(relationship)) {
    r <- relationship[[name]][1];tau2 <- relationship[[name]][2]
    receiver_prior <- list(residual_variance=1,
      effect_variance=.04*r*r+tau2,inclusion_probability=.1)
    separate <- fit_one("B",Y[[name]][,2],receiver_prior)
    message <- gshare_transfer(source,r,tau2,probability_floor=1e-6)
    file <- file.path(out,paste0(name,".gshare"))
    write_gshare(message,file,overwrite=TRUE)
    loaded <- read_gshare(file)
    stopifnot(identical(message,loaded))
    transferred <- fit_one("B",Y[[name]][,2],list(residual_variance=1),loaded)
    if(name=="heterogeneous")
      stopifnot(identical(transferred$effects,
        fit_one("B",Y[[name]][,2],list(residual_variance=1),message)$effects))
    truth <- drop(sim$X[sim$rows$test,,drop=FALSE]%*%effects[[name]][,2])
    score <- function(fit) sqrt(mean((drop(sim$X[sim$rows$test,,drop=FALSE]%*%
      fit$effects$mean)-truth)^2))
    metrics[[name]] <- data.frame(scenario=name,method=c("separate_bed","transfer_bed"),
      genetic_rmse=c(score(separate),score(transferred)),
      seconds=c(separate$seconds,transferred$seconds))
    saved[[name]] <- list(separate=separate,transfer=transferred)
    input <- file.path(out,paste0(name,"-statistics.txt"))
    con <- file(input,"wt")
    writeLines(as.character(p),con)
    write.table(markers[c("marker","center","scale")],con,
      row.names=FALSE,col.names=FALSE,quote=FALSE)
    for(site in c("A","B")) {
      i <- sim$rows[[site]];x <- sim$X[i,,drop=FALSE]
      response <- Y[[name]][i,match(site,c("A","B"))]
      writeLines(paste(site,length(i),sprintf("%.17g",sum(response^2))),con)
      writeLines(paste(sprintf("%.17g",c(crossprod(x,response),crossprod(x))),
        collapse=" "),con)
    }
    close(con)
  }
  if(mode=="fit") {
    dense_source <- fit_one("A",Y$heterogeneous[,1],source_prior,packed=FALSE)
    dense_msg <- gshare_transfer(dense_source,.7,.0204,probability_floor=1e-6)
    dense_receiver <- fit_one("B",Y$heterogeneous[,2],list(residual_variance=1),
      dense_msg,packed=FALSE)
    stopifnot(max(abs(source$effects$mean-dense_source$effects$mean))<.02,
      max(abs(source$effects$conditional_inclusion-
        dense_source$effects$conditional_inclusion))<.02,
      max(abs(saved$heterogeneous$transfer$effects$mean-
        dense_receiver$effects$mean))<.02)
    saved$dense <- list(source=dense_source,receiver=dense_receiver)
  }
  saveRDS(saved,file.path(out,"fits.rds"))
  table <- do.call(rbind,metrics)
  if(mode=="wrong") {
    previous <- read.csv(file.path(out,"bed-metrics.csv"))
    table <- rbind(subset(previous,scenario!="wrong"),table)
  }
  write.csv(table,file.path(out,"bed-metrics.csv"),row.names=FALSE)
  cat("Packed fits and dense reference inputs written to",out,"\n")
  if(mode=="fit") cat("Dense-versus-BED maximum source mean difference:",
    max(abs(source$effects$mean-dense_source$effects$mean)),"\n")
  cat("Run genomic_two_sites for each *-statistics.txt, then run this script with report.\n")
} else {
  rows <- read.csv(file.path(out,"bed-metrics.csv"))
  for(name in names(relationship)) {
    native <- file.path(out,paste0(name,"-native-effects.csv"))
    if(!file.exists(native)) stop("Missing native joint reference: ",native)
    effects_native <- read.csv(native)
    joint <- subset(effects_native,method=="joint_fixed" & site=="B")
    joint <- joint[match(colnames(sim$X),joint$marker),]
    stopifnot(nrow(joint)==p,!anyNA(joint$mean))
    truth <- drop(sim$X[sim$rows$test,,drop=FALSE]%*%effects[[name]][,2])
    predicted <- drop(sim$X[sim$rows$test,,drop=FALSE]%*%joint$mean)
    rows <- rbind(rows,data.frame(scenario=name,method="joint_fixed_dense",
      genetic_rmse=sqrt(mean((predicted-truth)^2)),
      seconds=read.csv(file.path(out,paste0(name,"-native-timing.csv")))$seconds[4]))
  }
  write.csv(rows,file.path(out,"comparison.csv"),row.names=FALSE)
  print(rows,row.names=FALSE)
}
