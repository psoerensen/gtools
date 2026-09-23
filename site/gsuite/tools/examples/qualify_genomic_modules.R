# Ten new fixed-seed attempts at the connected genomic qualification workflow.
# Run from the gsuite root. --replicate=1 is the one-seed preflight;
# --worker=1 and --worker=2 process odd/even seeds without replacing failures.
.libPaths(c("build/r-library", .libPaths()))
args <- commandArgs(TRUE)
worker_arg <- grep("^--worker=", args, value = TRUE)
replicate_arg <- grep("^--replicate=", args, value = TRUE)
stopifnot(length(worker_arg) <= 1L, length(replicate_arg) <= 1L)
worker <- if(length(worker_arg)) as.integer(sub("^--worker=", "", worker_arg)) else 0L
one <- if(length(replicate_arg)) as.integer(sub("^--replicate=", "", replicate_arg)) else NA_integer_
stopifnot(!is.na(worker), worker %in% 0:2,
          is.na(one) || one %in% 1:10)
root <- "build/examples/simulated-genomics/qualification"
dir.create(root, recursive = TRUE, showWarnings = FALSE)
plan <- data.frame(replicate = 1:10,
                   seed = 20262000L + 100L * (1:10))
plan_path <- file.path(root, "plan.csv")
if(file.exists(plan_path)) stopifnot(identical(read.csv(plan_path), plan)) else
  write.csv(plan, plan_path, row.names = FALSE)
scripts <- file.path("tools/examples", c("simulated_genomics_workflow.R",
  "simulated_genomics_summaries.R", "qualify_genomic_modules.R"))
hashes <- tools::md5sum(scripts)
selected <- if(!is.na(one)) one else if(worker == 0L) 1:10 else
  seq.int(worker, 10L, by = 2L)
scratch <- file.path(root, paste0("worker-", if(is.na(one)) worker else one))
for(i in selected) {
  path <- file.path(root, sprintf("replicate-%02d.rds", i))
  if(file.exists(path)) {
    old <- readRDS(path)
    stopifnot(identical(old$seed, plan$seed[i]), identical(old$scripts, hashes))
    cat("RETAINED|", i, "|", old$status, "\n", sep = "")
    next
  }
  stopifnot(identical(tools::md5sum(scripts), hashes))
  cat("REPLICATE_START|", i, "|seed=", plan$seed[i], "\n", sep = "")
  data_path <- file.path(scratch, "data.rds")
  resume <- !"--fresh" %in% args && file.exists(data_path) &&
    identical(readRDS(data_path)$seed, plan$seed[i])
  options(gsuite.example.out = scratch, gsuite.example.seed = plan$seed[i],
          gsuite.example.resume = resume, gsuite.example.figures = FALSE,
          gsuite.qualify.loco.root = NULL)
  began <- Sys.time()
  e <- new.env(parent = globalenv())
  error <- tryCatch({
    sys.source(scripts[[1L]], envir = e)
    NULL
  }, error = function(condition) conditionMessage(condition))
  result <- list(replicate = i, seed = plan$seed[i],
                 status = if(is.null(error)) "completed" else "failed",
                 error = error, started = began, finished = Sys.time(),
                 scripts = hashes,
                 results = if(is.null(error)) e$compact else NULL)
  if(!is.null(error)) result$stage_files <- list.files(scratch, pattern = "[.]rds$")
  saveRDS(result, path)
  cat("REPLICATE_END|", i, "|", result$status, "|",
      if(is.null(error)) "" else error, "\n", sep = "")
  rm(e, result)
  gc()
}
cat("WORKER_COMPLETED|", worker, "\n", sep = "")
