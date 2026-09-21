# Runnable GWAS -> native gene evidence -> pathway analysis demonstration.
# Outputs stay in one directory; genotypes are used only to create this example.
library(gsuite)
stopifnot(requireNamespace("gsim",quietly=TRUE))
out<-getOption("gsuite.gsea.out","build/examples/gwas-pathways")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out<-normalizePath(out,mustWork=TRUE)
set.seed(20260921)
n<-10000L;ngenes<-80L;m<-3L*ngenes
ids<-paste0("m",seq_len(m));gene_ids<-paste0("gene",seq_len(ngenes))
frequencies<-runif(m,.15,.45)
W<-matrix(rbinom(n*m,2,rep(frequencies,each=n)),n,m,dimnames=list(NULL,ids))
for(j in seq_len(m)) if((j-1L)%%60L) {
  copy<-runif(n)<.25;W[copy,j]<-W[copy,j-1L]
}
genes<-stats::setNames(lapply(seq_len(ngenes),function(g)ids[(3*g-2):(3*g)]),gene_ids)
blocks<-split(gene_ids,rep(paste0("chr",1:4),each=20))
sets<-list(pathway_A=gene_ids[1:15],pathway_B=gene_ids[31:45],
           pathway_C=gene_ids[61:75])
# Save a counted-allele A reference and retain all within-chromosome correlations.
prefix<-file.path(out,"reference")
con<-file(paste0(prefix,".bed"),"wb")
writeBin(as.raw(c(0x6c,0x1b,0x01)),con)
for(j in seq_len(m)) {
  codes<-c(3L,2L,0L)[W[,j]+1L]
  writeBin(as.raw(colSums(matrix(codes,4L)*c(1,4,16,64))),con)
}
close(con)
write.table(data.frame(rep(1:4,each=60),ids,0,rep(seq_len(60)*1000,4),"A","G"),
  paste0(prefix,".bim"),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(data.frame(seq_len(n),seq_len(n),0,0,0,-9),paste0(prefix,".fam"),
  quote=FALSE,row.names=FALSE,col.names=FALSE)
LD<-ldprep(gs_gprep(bedfiles=paste0(prefix,".bed")),reference="artificial-in-sample",
  assembly="artificial",task="sparseld",out_prefix=file.path(out,"LD"),
  max_distance_bp=0,max_distance_variants=60,r2=0,nthreads=1,overwrite=TRUE)
md<-LD$resource$markers
X<-scale(W,center=TRUE,scale=FALSE);d<-colSums(X^2)
make_trait<-function(active,seed) {
  b<-rep(0,m);b[ids %in% unlist(genes[active],use.names=FALSE)]<-.018
  gv<-var(drop(scale(W)%*%b))
  simulation<-gsim::gsim(W=W,architecture="fixed",nt=1,beta=b,h2=gv/(gv+1),
    standardize_W=TRUE,scale_effects=FALSE,seed=seed,compute_sumstats=FALSE)
  y<-drop(simulation$Y);yc<-y-mean(y);xy<-drop(crossprod(X,yc))
  beta<-xy/d;se<-sqrt((sum(yc^2)-xy^2/d)/((n-2)*d))
  data.frame(marker=md$marker_id,allele1=md$allele1,allele2=md$allele2,z=beta/se)
}
stat<-list(trait_A=make_trait(sets$pathway_A,111),
           trait_B=make_trait(c(sets$pathway_A,sets$pathway_B),222))
# Counts come from these GWAS individuals, not from reference-N labels.
minor_count<-pmin(colSums(W),2*n-colSums(W))
metadata<-data.frame(gene=gene_ids,sample_size=n,
  mean_minor_allele_count=vapply(genes,function(g)mean(minor_count[match(g,ids)]),numeric(1)))
started<-proc.time()[["elapsed"]]
evidence<-glma_genes(stat,LD,genes,blocks,metadata,
  control=list(independent_blocks=TRUE,tail_method="controlled_series"))
stopifnot(all(vapply(evidence$stat,function(s)all(s$p_available&s$p_value>0&s$p_value<1),logical(1))))
fits<-list(
  ora=gsea(lapply(evidence$stat,function(s)transform(s,selected=p_value<.05)),sets,"ora",control=list(adjustment="BH")),
  preranked=gsea(evidence$stat,sets,"preranked",control=list(seed=31,replicates=999,adjustment="BH")),
  competitive=gsea(evidence$stat,sets,"competitive",sampling=evidence$sampling,control=list(adjustment="BH")),
  magma=gsea(evidence$stat,sets,"magma",sampling=evidence$sampling,control=list(adjustment="BH")),
  bayesc=gsea(evidence$stat["trait_A"],sets,"bayesc",sampling=evidence$sampling,
    prior=list(residual_variance=1,effect_variance=1),control=list(burnin=100,sampling_sweeps=1000,seeds=c(11,29)))
)
elapsed<-proc.time()[["elapsed"]]-started
results<-do.call(rbind,lapply(names(fits),function(method) {
  p<-fits[[method]]$pathways
  data.frame(method=method,trait=p$trait,pathway=p$pathway,
    p_value=if("p_value"%in%names(p))p$p_value else NA_real_,
    pip=if("pip"%in%names(p))p$pip else NA_real_)
}))
write.csv(results,file.path(out,"pathways.csv"),row.names=FALSE)
saveRDS(list(evidence=evidence,fits=fits,sets=sets,stat=stat,LDlist=LD,
            metadata=metadata,gene_to_marker=genes,elapsed_seconds=elapsed),file.path(out,"workflow.rds"))
png(file.path(out,"gwas-pathways.png"),width=1800,height=1100,res=150)
par(mfrow=c(2,2),mar=c(4,4,3,1))
for(trait in names(stat)) {
  s<-evidence$stat[[trait]];active<-if(trait=="trait_A")sets$pathway_A else c(sets$pathway_A,sets$pathway_B)
  plot(seq_len(ngenes),-log10(s$p_value),pch=16,col=ifelse(s$gene%in%active,"#c05a22","#38678f"),
    xlab="Artificial gene index",ylab="Gene -log10(p)",main=trait)
  legend("topright",c("Contains simulated effects","Other genes"),pch=16,col=c("#c05a22","#38678f"),bty="n",cex=.8)
}
a<-fits$magma$pathways
values<-sapply(names(stat),function(t)-log10(pmax(a$p_value[a$trait==t],.Machine$double.xmin)))
barplot(t(values),beside=TRUE,names.arg=names(sets),col=c("#38678f","#c05a22"),
        ylab="Pathway -log10(p)",main="Bounded MAGMA workflow")
legend("topright",names(stat),fill=c("#38678f","#c05a22"),bty="n",cex=.8)
barplot(fits$bayesc$pathways$pip,names.arg=names(sets),ylim=c(0,1),col="#3b8775",
        ylab="Posterior inclusion probability",main="BayesC: trait_A (working score model)")
dev.off()
print(evidence)
print(results,row.names=FALSE)
cat("Gene preparation and pathway fits:",elapsed,"seconds\n")
