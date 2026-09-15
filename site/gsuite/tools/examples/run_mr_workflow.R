# A constructed MR example, not a simulation or statistical calibration study.
library(gsuite)
out <- getOption("gsuite.mr.out","build/examples/mr-workflow")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,mustWork=TRUE)
h <- matrix(1,1,1)
for(i in 1:4) h <- kronecker(h,matrix(c(1,1,1,-1),2))
W <- 1+h[rep(1:16,64),2:13] # 1,024 artificial reference individuals; exactly orthogonal markers.
ids <- paste0("mr",1:12);prefix <- file.path(out,"reference")
con <- file(paste0(prefix,".bed"),"wb");writeBin(as.raw(c(0x6c,0x1b,0x01)),con)
for(j in seq_len(ncol(W))) {
  codes <- c(3L,2L,0L)[W[,j]+1L]
  writeBin(as.raw(colSums(matrix(codes,4L)*c(1,4,16,64))),con)
}
close(con)
write.table(data.frame(1,ids,0,seq_along(ids)*1000,"A","G"),paste0(prefix,".bim"),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(data.frame(1:nrow(W),1:nrow(W),0,0,0,-9),paste0(prefix,".fam"),quote=FALSE,row.names=FALSE,col.names=FALSE)
LD <- ldprep(gs_gprep(bedfiles=paste0(prefix,".bed")),reference="artificial-MR",ancestry="artificial",
  assembly="artificial",task="sparseld",out_prefix=file.path(out,"LD"),
  max_distance_bp=0,max_distance_variants=20,r2=0,nthreads=1,overwrite=TRUE)
i <- 0:11;sign <- ifelse(i%%2,-1,1)
X <- data.frame(marker=ids,allele1="A",allele2="G",beta=sign*(.06+.01*i),se=.004,n=100000)
Y <- X;Y$beta <- sign*(.3*abs(X$beta)+.001*sin(i));Y$se <- .006;Y$n <- 90000
stat <- list(X=X,Y=Y)
prepared <- gcorr(stat,LD,method="gsmr",task="prepare_mr",exposure="X",outcome="Y",
  control=list(non_overlapping_samples=TRUE,effect_units=c(X="X units",Y="Y units"),
               ancestry=c(X="artificial",Y="artificial")))
stopifnot(prepared$available)
fits <- lapply(c("gsmr","gsmr2","ivw","egger","weighted_median"),function(method)
  gcorr(prepared,method=method,task="mr",control=if(method=="weighted_median")
    list(seed=19,bootstrap_replicates=1000) else list()))
names(fits) <- vapply(fits,`[[`,character(1),"method")
stopifnot(all(vapply(fits,`[[`,logical(1),"available")))
results <- do.call(rbind,lapply(fits,`[[`,"estimates"));print(results)
saveRDS(list(stat=stat,LD=LD,prepared=prepared,fits=fits),file.path(out,"workflow.rds"))
write.csv(results,file.path(out,"estimates.csv"),row.names=FALSE)
png(file.path(out,"mr.png"),width=1200,height=850,res=140)
plot(fits$gsmr,main="Constructed MR example: generating slope 0.3")
dev.off()
cat("Saved MR workflow to",out,"\n")
