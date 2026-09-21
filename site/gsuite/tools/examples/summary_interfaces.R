# Small end-to-end interface demonstration, not a calibration benchmark.
library(gsuite)
stopifnot(requireNamespace("gsim",quietly=TRUE))
out<-getOption("gsuite.summary.out","build/examples/summary-interfaces")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(2109);n<-10000L;m<-24L;ids<-paste0("m",seq_len(m))
W<-matrix(rbinom(n*m,2,.3),n,m,dimnames=list(NULL,ids))
prefix<-file.path(out,"reference")
con<-file(paste0(prefix,".bed"),"wb");writeBin(as.raw(c(0x6c,0x1b,0x01)),con)
for(j in seq_len(m)) {
  codes<-c(3L,2L,0L)[W[,j]+1L]
  writeBin(as.raw(colSums(matrix(codes,4L)*c(1,4,16,64))),con)
}
close(con)
write.table(data.frame("1",ids,seq_len(m)/10,seq_len(m)*1000,"A","G"),paste0(prefix,".bim"),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(data.frame(seq_len(n),seq_len(n),0,0,0,-9),paste0(prefix,".fam"),quote=FALSE,row.names=FALSE,col.names=FALSE)
Glist<-gs_gprep(bedfiles=paste0(prefix,".bed"))
LD<-ldprep(Glist,reference="artificial-full",assembly="artificial",task="sparseld",out_prefix=file.path(out,"LD"),max_distance_bp=0,max_distance_variants=m,r2=0,nthreads=1,overwrite=TRUE)
beta<-c(.2,-.15,.1,rep(0,m-3))
simulation<-gsim::gsim(W=W,architecture="fixed",nt=1,beta=beta,h2=.3,standardize_W=TRUE,scale_effects=TRUE,seed=7,compute_sumstats=FALSE)
y<-drop(simulation$Y)
summaries<-function(rows) {
  X<-scale(W[rows,,drop=FALSE],center=TRUE,scale=FALSE);yc<-y[rows]-mean(y[rows]);d<-colSums(X^2);xy<-drop(crossprod(X,yc))
  b<-xy/d;se<-sqrt((sum(yc^2)-xy^2/d)/((length(rows)-2)*d))
  raw<-data.frame(marker=ids,allele1="A",allele2="G",chromosome="1",position_bp=seq_len(m)*1000,beta=b,se=se,n=length(rows),p_value=2*pnorm(-abs(b/se)))
  gprep_stat(raw,LD,task="standardize")$stat
}
full<-summaries(seq_len(n));training<-summaries(1:6000);validation<-summaries(6001:n)
ridge<-gscore(full,LD,method="ridge",control=list(penalty=.1))
# Disjoint GWAS samples; full-sample LD is an explicit shared reference here.
selected<-gscore(full,LD,task="select",training=training,validation=validation,validation_LDlist=LD,
  positions_cm=setNames(seq_len(m)/10,ids),candidates=list(ridge=list(method="ridge",penalty=.1),sparse=list(method="l1",lambda=.01,delta=.1)),control=list(independent=TRUE))
scores<-gscore(task="apply",Glist=Glist,weights=selected,allele_frequencies=setNames(colMeans(W)/2,ids))
vb<-gbayes(full,LD,"variational",prior=list(p=.1,h2=.3))
mixture<-gbayes(full,LD,"mixture_variational",prior=list(variances=c(0,.01),weights=c(.9,.1)))
gibbs<-gbayes(full,LD,"point_normal",prior=list(p=.1,h2=.3),control=list(seed=7,burnin=100,sampling_sweeps=200))
shrinkage<-gbayes(full,LD,"continuous_shrinkage",control=list(block_offsets=c(0,m),seed=7,burnin=100,sampling_sweeps=200))
saveRDS(list(reference=LD,stat=full,ridge=ridge,selected=selected,scores=scores,vb=vb,mixture=mixture,gibbs=gibbs,shrinkage=shrinkage),file.path(out,"workflow.rds"))
print(summary(selected));print(summary(vb))
cat("Saved interface demonstration to",normalizePath(out),"\n")
