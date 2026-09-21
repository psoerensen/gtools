# Ten fixed-design repetitions of the connected example, from the repo root.
# Rscript --vanilla tools/examples/simulated_genomics_replicates.R
# Optional --worker=1 or --worker=2 splits odd/even replicates across two R jobs.
# Existing compact results are retained. Failed seeds are recorded, never replaced.
args <- commandArgs(TRUE)
worker_arg <- grep("^--worker=",args,value=TRUE)
worker <- if(length(worker_arg)) as.integer(sub("^--worker=","",worker_arg)) else 0L
stopifnot(length(worker)==1L,!is.na(worker),worker %in% 0:2)
root <- "build/examples/simulated-genomics/replicates"
dir.create(root,recursive=TRUE,showWarnings=FALSE)
# Space seeds so one replicate's gsim seed (seed + 1) is never another
# replicate's genotype/effect seed.
plan <- data.frame(replicate=1:10,seed=20261000L+100L*(1:10))
plan_path <- file.path(root,"plan.csv")
if(file.exists(plan_path)) stopifnot(identical(read.csv(plan_path),plan)) else
  write.csv(plan,plan_path,row.names=FALSE)
scripts <- file.path("tools/examples",c("simulated_genomics_workflow.R",
  "simulated_genomics_summaries.R","simulated_genomics_replicates.R"))
hashes <- tools::md5sum(scripts)
# A single reusable directory per worker bounds large genotype/LD/checkpoint files.
scratch <- file.path(root,paste0("worker-",worker))
selected <- if(worker==0L) 1:10 else seq(worker,10L,by=2L)
for(i in selected) {
  result_path <- file.path(root,sprintf("replicate-%02d.rds",i))
  if(file.exists(result_path)) {
    old <- readRDS(result_path)
    stopifnot(identical(old$seed,plan$seed[i]),identical(old$scripts,hashes))
    cat("REPLICATE_RETAINED|",i,"|",old$status,"\n",sep=""); next
  }
  stopifnot(identical(tools::md5sum(scripts),hashes))
  cat("REPLICATE_START|",i,"|seed=",plan$seed[i],"\n",sep=""); flush.console()
  options(gsuite.example.out=scratch,gsuite.example.seed=plan$seed[i],
    gsuite.example.resume=FALSE,gsuite.example.figures=FALSE)
  began <- Sys.time()
  e <- new.env(parent=globalenv())
  error <- tryCatch({sys.source(scripts[1],envir=e);NULL},error=function(e)conditionMessage(e))
  result <- list(replicate=i,seed=plan$seed[i],status=if(is.null(error)) "completed" else "failed",
    error=error,started=began,finished=Sys.time(),scripts=hashes,
    results=if(is.null(error)) e$compact else NULL)
  # On failure retain the evidence already computed, without retaining large fits.
  if(!is.null(error)) {
    result$stage_metadata <- lapply(list.files(scratch,pattern="[.]rds$",full.names=TRUE),function(p) {
      x <- readRDS(p)
      if(is.list(x) && identical(x$seed,plan$seed[i])) {x$value<-NULL;x} else NULL
    })
  }
  saveRDS(result,result_path)
  cat("REPLICATE_END|",i,"|",result$status,"|",error,"\n",sep=""); flush.console()
  rm(e,result);gc()
}
cat("REPLICATE_WORKER|completed\n")
