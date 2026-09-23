# Self-contained gsim preparation and native gshare R example; run from gsuite.
# Simulation follows the gshare genomic demonstration; no sibling source dependency.
.libPaths(c("build/r-library",.libPaths()))
library(gsuite)
out<-"build/examples/gshare-workflow"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out<-normalizePath(out,winslash="/",mustWork=TRUE)
existing<-Sys.getenv("GSHARE_EXAMPLE_SIMULATION",unset="")
if(nzchar(existing)) {
  stopifnot(file.copy(existing,file.path(out,"simulation.rds"),overwrite=TRUE))
} else {

  stopifnot(requireNamespace("gsim",quietly=TRUE))
  started <- proc.time()[[3]]
  m <- 100L; n <- 3000L; donors <- 100L
  markers <- paste0("m",seq_len(m)); donor_ids <- paste0("donor",seq_len(donors))
  set.seed(20260923)
  # Synthetic phased reference, with donor haplotypes correlated within ten-marker blocks.
  H <- matrix(0L,donors*2L,m)
  for(j in seq_len(m)) {
    H[,j] <- rbinom(2L*donors,1L,.3)
    if((j-1L)%%10L) {keep<-runif(2L*donors)<.8; H[keep,j]<-H[keep,j-1L]}
  }
  freq <- colMeans(H); stopifnot(all(freq>.05 & freq<.95))
  vcf <- file.path(out,"reference.vcf")
  header <- paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",donor_ids),collapse="\t")
  records <- vapply(seq_len(m),function(j) paste(c("1",j*1000L,markers[j],"A","G",".","PASS",".","GT",
      paste0(H[seq_len(donors),j],"|",H[donors+seq_len(donors),j])),collapse="\t"),character(1))
  writeLines(c("##fileformat=VCFv4.2",header,records),vcf)
  map <- data.frame(chromosome="1",variant_id=markers,genetic_position_cm=seq(0,1,length.out=m))
  reference <- gsim::gsim_import_vcf(vcf,map,file.path(out,"reference"),unsupported="error",overwrite=TRUE)
  pedigree <- gsim::gsim_pedigree(n_generations=2L,animals_per_generation=(n+2L)/2L,founder_generations=1L,
      new_founder_probability=1,unknown_sire_probability=0,unknown_dam_probability=0,
      unphenotyped_probability=0,permute_external_order=FALSE,seed=7010)
  # The pedigree constructor retains two descendants. Export them but use only
  # the 3,000 original founders in every analysis and holdout assignment.
  founder_ids<-pedigree$pedigree$animal[is.na(pedigree$pedigree$sire)&is.na(pedigree$pedigree$dam)]
  stopifnot(length(founder_ids)==n)
  founders <- gsim::gsim_simulate_founders(reference=gsim::gsim_reference(reference$prefix),
      founder_ids=founder_ids,populations=setNames(rep("P",donors),donor_ids),
      ancestry_weights=c(P=1),mutation_age=setNames(rep(1e9,m),markers),N=c(P=100),Ne=c(P=1000),
      rho=c(P=.02),seed=7011,output=file.path(out,"founders"),batch_size=512L,threads=1L,overwrite=TRUE)
  bed <- gsim::gsim_simulate_pedigree(founders,pedigree,seed=7012,output=file.path(out,"genotypes"),format="bed",overwrite=TRUE)
  # Bounded example decoder of standard SNP-major PLINK BED (G is counted A1).
  con <- file(bed$paths[["bed"]],"rb")
  stopifnot(identical(readBin(con,"raw",3L),as.raw(c(0x6c,0x1b,0x01))))
  physical_n<-length(pedigree$canonical_order)
  W <- matrix(0,physical_n,m,dimnames=list(pedigree$canonical_order,markers))
  for(j in seq_len(m)) {
    bytes <- as.integer(readBin(con,"raw",ceiling(physical_n/4)))
    stopifnot(length(bytes)==ceiling(physical_n/4))
    # vapply is byte-columns; flatten in byte order directly.
    code <- as.vector(vapply(bytes,function(x) bitwAnd(bitwShiftR(x,c(0L,2L,4L,6L)),3L),integer(4)))[seq_len(physical_n)]
    W[,j] <- c(2,NA,1,0)[code+1L]
  }
  close(con)
  bim<-read.table(bed$paths[["bim"]]); fam<-read.table(bed$paths[["fam"]])
  stopifnot(!anyNA(W),identical(as.character(bim[[2]]),markers),all(bim[[5]]=="G"),
            identical(as.character(fam[[2]]),rownames(W)))
  W<-W[founder_ids,,drop=FALSE]
  center<-2*freq; sx<-sqrt(2*freq*(1-freq)); X<-sweep(sweep(W,2,center),2,sx,"/")
  # Training/test assignments precede phenotype generation; no holdout-derived preprocessing.
  set.seed(7013); order<-sample.int(n); rows<-list(A=order[1:500],B=order[501:2000],test=order[2001:3000])
  set.seed(7014); causal<-sort(sample.int(m,10)); B<-matrix(0,m,2,dimnames=list(markers,c("A","B")))
  B[causal,1]<-rnorm(10,sd=.2); B[causal,2]<-.7*B[causal,1]+rnorm(10,sd=sqrt(.0204))
  vg<-apply(X%*%B,2,var)
  simulation<-gsim::gsim(W=X,nt=2,architecture="fixed",beta=B,h2=vg/(vg+1),re=0,
      standardize_W=FALSE,scale_effects=FALSE,seed=7015,compute_sumstats=FALSE)
  stopifnot(max(abs(simulation$G-X%*%B))<1e-10,max(abs(simulation$Sigma_e-diag(2)))<1e-10)
  path<-file.path(out,"statistics.txt"); con<-file(path,"wt")
  writeLines(as.character(m),con)
  write.table(data.frame(markers,center,sx),con,row.names=FALSE,col.names=FALSE,quote=FALSE)
  for(s in 1:2) {
    idx<-rows[[s]]; y<-simulation$Y[idx,s]; xx<-X[idx,,drop=FALSE]
    writeLines(paste(c(names(rows)[s],length(idx),sprintf("%.17g",sum(y*y))),collapse=" "),con)
    writeLines(paste(sprintf("%.17g",c(crossprod(xx,y),crossprod(xx))),collapse=" "),con)
  }
  close(con)
  saveRDS(list(W=W,X=X,B=B,rows=rows,simulation=simulation,causal=causal,center=center,scale=sx,
      reference_frequency=freq,seconds=proc.time()[[3]]-started,session=sessionInfo(),
      gsim_version=as.character(packageVersion("gsim"))),file.path(out,"simulation.rds"))
  cat("Prepared",n,"gsim founder genotypes at",m,"markers; elapsed",proc.time()[[3]]-started,"seconds\n")

}
sim<-readRDS(file.path(out,"simulation.rds"))
markers<-data.frame(marker=colnames(sim$X),allele="G",center=sim$center,scale=sim$scale)
control<-list(burnin=2000,sampling_sweeps=8000,seeds=c(2027,2029))
prior<-list(residual_variance=1,effect_variance=.04,inclusion_probability=.1)
local_fit<-function(site,transfer=NULL) {
  s<-match(site,c("A","B"));idx<-sim$rows[[site]]
  gshare(X=sim$X[idx,,drop=FALSE],y=sim$simulation$Y[idx,s],site=site,
    markers=markers,trait="site phenotype",units="simulation units",
    covariates="known zero intercept; donor-reference centering/scaling",
    prior=if(is.null(transfer)) prior else list(residual_variance=1),
    transfer=transfer,disjoint=!is.null(transfer),control=control)
}
source<-local_fit("A")
message<-gshare_transfer(source,multiplier=.7,heterogeneity_variance=.0204,probability_floor=1e-6)
file<-file.path(out,"site-A.gshare")
write_gshare(message,file,overwrite=TRUE)
received<-read_gshare(file)
stopifnot(identical(received,message))
separate<-local_fit("B");transfer<-local_fit("B",received)
test<-sim$rows$test;truth<-sim$simulation$G[test,2]
metrics<-do.call(rbind,lapply(list(separate=separate,transfer=transfer),function(f) {
  prediction<-drop(sim$X[test,]%*%f$effects$mean)
  data.frame(genetic_rmse=sqrt(mean((prediction-truth)^2)),genetic_correlation=cor(prediction,truth),
    causal_pip=mean(f$effects$pip[sim$causal]),seconds=f$seconds)
}))
print(metrics);print(message);print(transfer)
write.csv(metrics,file.path(out,"comparison.csv"))
saveRDS(list(source=source,separate=separate,transfer=transfer,message=message),file.path(out,"fits.rds"))
png(file.path(out,"prediction.png"),width=1000,height=600,res=120)
barplot(metrics$genetic_rmse,names.arg=rownames(metrics),col=c("#617487","#d78028"),
  ylab="Site-B held-out genetic-value RMSE",ylim=c(0,max(metrics$genetic_rmse)*1.15))
dev.off()
png(file.path(out,"pips.png"),width=1200,height=600,res=120)
par(mfrow=c(1,2),mar=c(4,4,5,1))
for(name in c("separate","transfer")) {
  f<-get(name);active<-seq_len(nrow(f$effects))%in%sim$causal
  plot(f$effects$pip,ylim=c(0,1.25),yaxt="n",col=ifelse(active,"#c44332","#617487"),pch=ifelse(active,17,16),
    xlab="Marker",ylab="PIP",main=name)
  axis(2,at=seq(0,1,.2))
  legend("topright",c("True causal","Other"),col=c("#c44332","#617487"),pch=c(17,16),bty="n")
}
dev.off()
