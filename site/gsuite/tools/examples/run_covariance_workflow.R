# Constructed covariance example: workflow illustration, not a calibration study.
library(gsuite)
out <- getOption("gsuite.covariance.out","build/examples/covariance-workflow")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
tt <- c("A","B","C","D")
loading <- c(.8,.7,.6,.5)
G <- tcrossprod(loading)+diag(c(.4,.5,.6,.7))
dimnames(G) <- list(tt,tt)
labels <- unlist(lapply(seq_along(tt),function(j) paste(tt[j:4],tt[j],sep="|")))
V <- .0004*.3^abs(outer(1:10,1:10,"-"));dimnames(V) <- list(labels,labels)
input <- gcorr_covariance(G,V,continuous=TRUE,units=setNames(rep("standardized units",4),tt),
  provenance="Constructed one-factor covariance; illustrative full sampling covariance")
model <- list(variables=c(tt,"F"),
  parameters=data.frame(name=c(paste0("l",tt),paste0("v",tt)),start=c(rep(.5,4),rep(.6,4)),lower=c(rep(-Inf,4),rep(1e-8,4))),
  entries=data.frame(matrix=c(rep("directed",4),rep("disturbance",5)),
    row=c(tt,tt,"F"),column=c(rep("F",4),tt,"F"),
    parameter=c(paste0("l",tt),paste0("v",tt),""),value=c(rep(0,8),1)))
sem <- gcorr(input,method="sem",task="fit",model=model)
network <- gcorr(input,method="ggm",task="network")
graph <- list(edges=matrix(c("A","B","B","C","C","D"),ncol=2,byrow=TRUE))
ggm <- gcorr(input,method="ggm",task="fit",model=graph)
stopifnot(sem$available,sem$converged,network$available,ggm$available,ggm$converged)
saveRDS(list(input=input,sem=sem,network=network,ggm=ggm),file.path(out,"workflow.rds"))
write.csv(sem$estimates,file.path(out,"sem-estimates.csv"),row.names=FALSE)
write.csv(ggm$edges,file.path(out,"ggm-edges.csv"),row.names=FALSE)
png(file.path(out,"covariance-models.png"),width=1400,height=650)
par(mfrow=c(1,2),mar=c(4,4,4,2));plot(network);plot(ggm)
dev.off()
print(sem);print(ggm)
