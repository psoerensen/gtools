library(gsuite)
set.seed(2401)
genes <- paste0("gene", seq_len(120))
sets <- list(pathway_A=genes[1:25], pathway_B=genes[21:45],
             pathway_C=genes[71:95])
X <- vapply(sets, function(s) as.numeric(genes %in% s), numeric(120))
rownames(X) <- genes
noise <- matrix(rnorm(240),120,2)
stat <- list(
  trait_A=data.frame(gene=genes, score=2*X[,1]+noise[,1]),
  trait_B=data.frame(gene=genes, score=1.5*X[,1]+X[,2]+.3*noise[,1]+sqrt(.91)*noise[,2])
)
# Independence across genes is part of this artificial data-generating model.
sampling <- list(gene_covariance="independent")
ora <- gsea(transform(stat$trait_A, selected=score>2), sets, "ora",
            control=list(adjustment="BH"))
ranked <- gsea(stat$trait_A, sets, "preranked",
               control=list(replicates=999, seed=2401, adjustment="BH"))
competitive <- gsea(stat, sets, "competitive", sampling=sampling,
                    control=list(adjustment="BH"))
# Artificial p-values illustrate the API, not a calibrated GWAS-to-gene pipeline.
gene_p <- transform(stat$trait_A, p_value=pnorm(score,lower.tail=FALSE),
                    marker_count=10, effective_parameters=8,
                    sample_size=10000, mean_minor_allele_count=2000)
magma <- gsea(gene_p, sets, "magma", sampling=sampling,
              control=list(adjustment="BH"))
scalar <- gsea(stat$trait_A, sets, "bayesc", sampling=sampling,
               prior=list(residual_variance=1,effect_variance=1),
               control=list(burnin=100,sampling_sweeps=1000,seeds=c(11,29)))
Omega <- matrix(c(1,.3,.3,1),2,dimnames=list(names(stat),names(stat)))
G <- diag(2);dimnames(G)<-dimnames(Omega)
joint <- gsea(stat, sets, "bayesr",
               sampling=list(gene_covariance="independent",trait_covariance=Omega),
               prior=list(effect_covariance=G),
               control=list(burnin=100,sampling_sweeps=1000,seeds=c(11,29)))
print(competitive)
print(joint)
summary(joint)
# Optional interactive graphics and persistence:
# plot(joint, trait="trait_A")
# plot(ranked)
# saveRDS(joint, "pathway-fit.rds")
