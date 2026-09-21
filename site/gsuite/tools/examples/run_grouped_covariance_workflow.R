# Constructed four-trait example for annotation-specific SEM and GGM.
# This demonstrates the interface; one simulation cannot calibrate inference.
library(gsuite)
library(gsim)

out <- getOption("gsuite.grouped.covariance.out",
                 "build/examples/grouped-covariance-workflow")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20260916)

n <- 10000L
nref <- 2000L
nblocks <- 40L
size <- 20L
m <- nblocks * size
markers <- sprintf("snp%04d", seq_len(m))
blocks <- rep(sprintf("block%02d", seq_len(nblocks)), each = size)
groups <- rep(c("coding", "regulatory"), each = m / 2L)
annotation <- cbind(coding = as.numeric(groups == "coding"),
                    regulatory = as.numeric(groups == "regulatory"))
rownames(annotation) <- markers

# Independent LD blocks with a copy-or-redraw dosage generator.
W <- matrix(0L, n + nref, m, dimnames = list(NULL, markers))
for (j in seq_len(m)) {
  W[, j] <- rbinom(n + nref, 2, .3)
  if ((j - 1L) %% size) {
    reuse <- runif(n + nref) < .4
    W[reuse, j] <- W[reuse, j - 1L]
  }
}
sim <- gsim::gsim(W = W, nt = 4L, h2 = rep(.35, 4),
                  vg = rep(.35, 4), rg = .4, re = 0,
                  architecture = "bayesc", pi = c(0, 1),
                  standardize_W = TRUE, scale_effects = TRUE,
                  seed = 20260917, compute_sumstats = FALSE)
traits <- c("A", "B", "C", "D")
phenotypes <- sim$Y[seq_len(n), , drop = FALSE]
association <- cor(W[seq_len(n), , drop = FALSE], phenotypes)
Z <- association * sqrt((n - 2) / (1 - association^2))
stat <- setNames(lapply(seq_along(traits), function(k) {
  data.frame(marker = markers, z = Z[, k], n = n, block = blocks)
}), traits)
rm(sim, phenotypes, association, Z)

# Write a small PLINK reference, then retain only its prepared LD resource.
prefix <- file.path(out, "reference")
con <- file(paste0(prefix, ".bed"), "wb")
writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
for (j in seq_len(m)) {
  codes <- c(3L, 2L, 0L)[W[n + seq_len(nref), j] + 1L]
  writeBin(as.raw(colSums(matrix(codes, 4L) * c(1, 4, 16, 64))), con)
}
close(con)
write.table(data.frame(rep(1L, m), markers, 0L, seq_len(m) * 1000L,
                       "A", "G"), paste0(prefix, ".bim"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)
write.table(data.frame(seq_len(nref), seq_len(nref), 0L, 0L, 0L, -9L),
            paste0(prefix, ".fam"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)
rm(W)

LD <- ldprep(gprep(bedfiles = paste0(prefix, ".bed")),
             reference = "simulation", task = "sparseld",
             out_prefix = file.path(out, "LD"), max_distance_bp = 0,
             max_distance_variants = size, r2 = 0, nthreads = 1,
             overwrite = TRUE)
LD <- ldprep(LDlist = LD, task = "scores", annotation = annotation)
saveRDS(LD, file.path(out, "LDlist.rds"))
stopifnot(all(file.remove(paste0(prefix, c(".bed", ".bim", ".fam")))))
source_fit <- gcorr(stat, LD, annotation = annotation, method = "gnova",
                    task = "rg", control = list(overlap_intercept = 0))
saveRDS(source_fit, file.path(out, "source-fit.rds"))

units <- setNames(rep("standardized trait units", 4), traits)
sem_model <- list(
  variables = c(traits, "F"),
  parameters = data.frame(name = c(paste0("l", traits), paste0("v", traits)),
                          start = c(rep(.5, 4), rep(.2, 4)),
                          lower = c(rep(-Inf, 4), rep(1e-8, 4))),
  entries = data.frame(matrix = c(rep("directed", 4), rep("disturbance", 5)),
                       row = c(traits, traits, "F"),
                       column = c(rep("F", 4), traits, "F"),
                       parameter = c(paste0("l", traits), paste0("v", traits), ""),
                       value = c(rep(0, 8), 1)))
sem <- gcorr(source_fit, method = "sem", task = "fit",
             annotation = c("coding", "regulatory"), model = sem_model,
             units = units)
# An illustrative fixed graph with A-D absent, specified independently of
# the estimated edges. The missing edge is constrained to zero in each group.
graph <- list(edges = matrix(c("A", "B", "A", "C", "B", "C",
                               "B", "D", "C", "D"),
                             ncol = 2, byrow = TRUE))
ggm <- gcorr(source_fit, method = "ggm", task = "fit",
             annotation = c("coding", "regulatory"), model = graph,
             units = units)

print(sem$diagnostics)
print(ggm$diagnostics)
stopifnot(sem$inference_available, ggm$inference_available,
          ggm$partial_correlations$inference_available,
          all(sem$diagnostics$available), all(ggm$diagnostics$available))
R <- matrix(0, 1, nrow(sem$estimates),
            dimnames = list("coding minus regulatory lA", sem$estimates$label))
R[1, c("coding::lA", "regulatory::lA")] <- c(1, -1)
sem_test <- gcorr(sem, method = "sem", task = "test", constraints = R)
P <- ggm$partial_correlations$estimates
R <- matrix(0, 1, nrow(P),
            dimnames = list("coding minus regulatory partial A-B", P$label))
R[1, c("coding::partial[A,B]", "regulatory::partial[A,B]")] <- c(1, -1)
ggm_test <- gcorr(ggm, method = "ggm", task = "test", constraints = R)
stopifnot(sem_test$available, ggm_test$available)
print(sem_test$results)
print(ggm_test$results)

saveRDS(list(sem = sem, ggm = ggm, sem_test = sem_test,
             ggm_test = ggm_test), file.path(out, "models.rds"))
write.csv(sem$estimates, file.path(out, "sem-estimates.csv"), row.names = FALSE)
write.csv(ggm$partial_correlations$estimates,
          file.path(out, "ggm-partial-correlations.csv"), row.names = FALSE)
write.csv(rbind(sem_test$results, ggm_test$results),
          file.path(out, "contrasts.csv"), row.names = FALSE)
write.csv(rbind(cbind(method = "sem", sem$diagnostics),
                cbind(method = "ggm", ggm$diagnostics)),
          file.path(out, "diagnostics.csv"), row.names = FALSE)

# Show point estimates without implying that one simulated draw calibrates SEs.
draw <- function(M, title) {
  cols <- grDevices::colorRampPalette(c("#2C6B9A", "white", "#B8582B"))(201)
  plot.new()
  plot.window(xlim = c(.5, 4.5), ylim = c(.5, 4.5), asp = 1)
  for (i in seq_len(4)) for (j in seq_len(4)) {
    y <- 5 - i
    color <- if (i == j) "#E8EBEF" else
      cols[pmax(1L, pmin(201L, round((M[i, j] + 1) * 100) + 1L))]
    rect(j - .49, y - .49, j + .49, y + .49,
         col = color, border = "white", lwd = 2)
    label <- if (i == j) "" else if (abs(M[i, j]) < .005) "0.00" else
      sprintf("%.2f", M[i, j])
    text(j, y, label, cex = 1.15)
  }
  axis(1, at = 1:4, labels = traits, tick = FALSE, line = -.3)
  axis(2, at = 4:1, labels = traits, tick = FALSE, las = 1, line = -.3)
  title(main = title, cex.main = 1.05)
  box()
}
grDevices::png(file.path(out, "grouped-covariance.png"),
               width = 1450, height = 1150, res = 145)
old <- par(mfrow = c(2, 2), mar = c(2.4, 2.8, 3.2, 1.1),
           oma = c(2, 0, 0, 0))
for (group in c("coding", "regulatory")) {
  G <- source_fit$matrices[[LD$reference]]$annotation_category[[group]]$cov
  draw(cov2cor(G), paste(group, "genetic correlation"))
}
for (group in c("coding", "regulatory")) {
  precision <- solve(ggm$fits[[group]]$fitted_covariance)
  partial <- -cov2cor(precision)
  draw(partial, paste(group, "GGM partial correlation"))
}
mtext("Artificial data | point estimates; uncertainty is in the saved tables",
      side = 1, outer = TRUE, line = .5, cex = .85)
par(old)
grDevices::dev.off()
