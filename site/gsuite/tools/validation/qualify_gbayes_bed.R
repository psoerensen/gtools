.libPaths(c(file.path('build','r-library'),.libPaths()))
library(gsuite)
library(sblr)
x <- readRDS('build/examples/simulated-genomics/data.rds')$value
n <- as.integer(Sys.getenv('COMPARE_N','3000'))
m <- as.integer(Sys.getenv('COMPARE_M','1000'))
burn <- as.integer(Sys.getenv('COMPARE_BURN','300'))
draws <- as.integer(Sys.getenv('COMPARE_DRAWS','1000'))
ids <- x$Glist$ids[seq_len(n)]
markers <- x$Glist$rsids[[1]][seq_len(m)]
g <- gprep(bedfiles=x$Glist$bedfiles,ids=ids,rsids=markers)
stopifnot(identical(g$rsids[[1]],markers),identical(g$study_ids,ids),
          all(g$source_marker_indices_0based[[1]]==seq_len(m)-1L))
sg <- x$Glist
sg$af[[1]][seq_len(m)] <- g$af[[1]]
y <- x$simulation$Y[seq_len(n),'A']
y <- y-mean(y)
names(y) <- ids
dat <- data.frame(id=ids,A=as.double(y))
results <- list()
for (method in c('bayesc','bayesr')) {
  cat('Starting',method,'n',n,'m',m,'burn',burn,'draws',draws,'\n')
  s_args <- list(y=y,Glist=sg,cls=list(seq_len(m)),method=method,
    nit=draws,nburn=burn,nchains=2L,chain_seeds=c(11L,29L),
    convergence='none',h2=.3,updateE=FALSE,updateB=FALSE,updatePi=FALSE,
    full_sweep_every=1L,null_skip_base=1L,null_skip_max=1L,
    verbose=FALSE)
  if (method=='bayesc')
    s_args <- c(s_args,list(pi_init=.05,pi_vb_init=.05,pi_prior_mean=.05))
  else s_args$pi <- c(.95,rep(.05/3,3))
  s_time <- system.time(s <- do.call(stblr_bed,s_args))
  print(s$input[c('B','E','pi','mixture_var')])
  prior <- list(residual_variance=as.double(s$input$E),
                effect_variance=as.double(s$input$B))
  if (method=='bayesc') prior$inclusion_probability <- .05 else {
    prior$weights <- as.double(s$input$pi)
    prior$variance_multipliers <- c(0,.01,.1,1)
  }
  g_time <- system.time(z <- gbayes(data=dat,Glist=g,trait='A',method=method,
    prior=prior,control=list(burnin=burn,sampling_sweeps=draws,
      seeds=c(11L,29L),threads=1L)))
  sbm <- as.double(s$bm[,1]); spip <- as.double(s$dm[,1])
  gbm <- z$estimates$mean; gpip <- z$estimates$pip
  stopifnot(length(sbm)==m,length(gbm)==m,
            identical(rownames(s$bm),z$estimates$marker))
  compare <- list(method=method,n=n,m=m,burn=burn,draws=draws,
    sblr_seconds=unname(s_time['elapsed']),gsuite_seconds=unname(g_time['elapsed']),
    effects_cor=cor(sbm,gbm),effects_rmse=sqrt(mean((sbm-gbm)^2)),
    pip_cor=cor(spip,gpip),pip_mae=mean(abs(spip-gpip)),
    sblr_pip_sum=sum(spip),gsuite_pip_sum=sum(gpip),
    sblr_nonzero=sum(spip>.5),gsuite_nonzero=sum(gpip>.5),
    sblr_input=s$input[c('B','E','pi')],gsuite_parameters=z$parameter_mean,
    per_marker=data.frame(marker=markers,sblr_mean=sbm,gsuite_mean=gbm,
                          sblr_pip=spip,gsuite_pip=gpip))
  print(compare[c('method','n','m','burn','draws','sblr_seconds',
                  'gsuite_seconds','effects_cor','effects_rmse','pip_cor',
                  'pip_mae','sblr_pip_sum','gsuite_pip_sum')])
  results[[method]] <- compare
  saveRDS(results,'build/examples/gbayes-bed-comparison/results.rds')
}
cat('Parameter learning on 1000 people and 200 markers\n')
x <- readRDS('build/examples/simulated-genomics/data.rds')$value
n <- 1000L;m <- 200L
ids <- x$Glist$ids[seq_len(n)]
g <- gprep(bedfiles=x$Glist$bedfiles,ids=ids,
           rsids=x$Glist$rsids[[1]][seq_len(m)])
sg <- x$Glist;sg$af[[1]][seq_len(m)] <- g$af[[1]]
y <- x$simulation$Y[seq_len(n),'A'];y <- y-mean(y);names(y) <- ids
dat <- data.frame(id=ids,A=as.double(y))
out <- list()
for(method in c('bayesc','bayesr')) {
  s_args <- list(y=y,Glist=sg,cls=list(seq_len(m)),method=method,
    nit=600,nburn=200,nchains=2L,chain_seeds=c(11L,29L),
    convergence='none',h2=.3,full_sweep_every=1L,
    null_skip_base=1L,null_skip_max=1L,verbose=FALSE)
  if(method=='bayesc')
    s_args <- c(s_args,list(pi_init=.05,pi_vb_init=.05,pi_prior_mean=.05))
  else s_args$pi <- c(.95,rep(.05/3,3))
  s_time <- system.time(s <- do.call(stblr_bed,s_args))
  prior <- list(residual_variance=as.double(s$input$E),
    effect_variance=as.double(s$input$B),
    effect_variance_prior=list(df=4,scale=as.double(s$input$ssb_prior)),
    residual_variance_prior=list(df=4,scale=as.double(s$input$sse_prior)))
  if(method=='bayesc') {
    prior$inclusion_probability <- .05
    prior$weight_prior <- c(.95,.05)*5e5
  } else {
    prior$weights <- as.double(s$input$pi)
    prior$variance_multipliers <- c(0,.01,.1,1)
    prior$weight_prior <- as.double(s$input$alpha)
  }
  g_time <- system.time(z <- gbayes(data=dat,Glist=g,trait='A',method=method,
    prior=prior,control=list(burnin=200,sampling_sweeps=600,
      seeds=c(11,29),estimate_effect_variance=TRUE,
      estimate_residual_variance=TRUE,estimate_weights=TRUE)))
  out[[method]] <- list(sblr_sec=s_time['elapsed'],gsuite_sec=g_time['elapsed'],
    sblr_vb=s$vb,sblr_ve=s$ve,sblr_pi=s$pi_mean,
    gsuite=z$parameter_mean,effects_cor=cor(s$bm[,1],z$estimates$mean),
    pip_cor=cor(s$dm[,1],z$estimates$pip))
  cat(method,'sblr seconds',unname(s_time['elapsed']),
      'gsuite seconds',unname(g_time['elapsed']),
      'effect correlation',out[[method]]$effects_cor,
      'PIP correlation',out[[method]]$pip_cor,'\n')
}
saveRDS(out,'build/examples/gbayes-bed-comparison/learning-results.rds')
