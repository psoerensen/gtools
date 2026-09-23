# One quantitative-trait LOCO scan on the connected simulation's genotype panel.
# Usage: Rscript --vanilla tools/examples/qualify_glma_loco.R WORKFLOW_DIR
# The main workflow treats each 100-marker region as a separate LD block. For
# LOCO, group those same markers into 20 artificial chromosomes; no genotypes or
# phenotypes are regenerated. This keeps the chromosome exclusion meaningful.
.libPaths(c("build/r-library", .libPaths()))
library(gsuite)
args <- commandArgs(TRUE)
root_arg <- getOption("gsuite.qualify.loco.root", if(length(args) == 1L) args[[1L]] else NULL)
stopifnot(is.character(root_arg), length(root_arg) == 1L)
run <- function(root) {
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
record <- readRDS(file.path(root, "data.rds"))
data <- record$value
seed <- record$seed
original <- sub("[.]bed$", "", data$Glist$bedfiles[[1L]])
prefix <- file.path(root, "loco-chromosomes")
paths <- paste0(prefix, c(".bed", ".bim", ".fam"))
on.exit(unlink(paths), add = TRUE)
stopifnot(file.copy(paste0(original, ".bed"), paths[[1L]], overwrite = TRUE),
          file.copy(paste0(original, ".fam"), paths[[3L]], overwrite = TRUE))
bim <- read.table(paste0(original, ".bim"), header = FALSE,
                  colClasses = c("integer", "character", "numeric", "integer",
                                 "character", "character"))
stopifnot(nrow(bim) == 50000L, identical(sort(unique(bim[[1L]])), 1:500))
region <- bim[[1L]]
bim[[1L]] <- (region - 1L) %/% 25L + 1L
bim[[4L]] <- ((region - 1L) %% 25L) * 1000000L + bim[[4L]]
write.table(bim, paths[[2L]], quote = FALSE, col.names = FALSE, row.names = FALSE)
ids <- sprintf("id%05d", 1:3000)
panel <- gprep(bedfiles = paths[[1L]], ids = ids)
y <- data$simulation$Y[, "A"]
names(y) <- sprintf("id%05d", seq_len(nrow(data$simulation$Y)))
started <- Sys.time()
fit <- glma(y, panel, method = "infinitesimal_loco",
            algorithm = "observation_pcg", threads = 2L,
            block_size = 256L, background_markers = seq.int(1L, 50000L, by = 100L),
            controls = list(seed = seed))
seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
a <- fit$associations
stopifnot(nrow(a) == 50000L, identical(a$marker, rownames(data$simulation$B)))
status <- table(a$status)
causal <- data$simulation$B[, "A"] != 0
ok <- a$status == "ok" & is.finite(a$p)
result <- list(seed = seed, started = started, seconds = seconds,
               method = "infinitesimal_loco", algorithm = "observation_pcg",
               traits = "A", individuals = 3000L,
               markers_tested = 50000L, background_markers = 500L,
               artificial_chromosomes = 20L, threads = 2L,
               completed = sum(ok), status = status,
               significant = sum(ok & a$p < .05 / 50000),
               significant_causal = sum(ok & a$p < .05 / 50000 & causal),
               significant_zero_effect = sum(ok & a$p < .05 / 50000 & !causal),
               diagnostics = fit$diagnostics, background = fit$fit,
               versions = gsuite_versions(), script_md5 = tools::md5sum(
                 "tools/examples/qualify_glma_loco.R"))
saveRDS(result, file.path(root, "loco-qualification.rds"))
cat("LOCO_COMPLETED|", sum(ok), "|", round(seconds, 2), " seconds\n", sep = "")
}
run(root_arg)
