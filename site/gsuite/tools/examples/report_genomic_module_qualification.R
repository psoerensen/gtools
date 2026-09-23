# Report the ten fixed new attempts; do not simulate or fit anything here.
root <- "build/examples/simulated-genomics/qualification"
options(gsuite.qualify.root = root)
source("tools/examples/report_simulated_genomics_replicates.R", local = TRUE)
summary <- readRDS(file.path(root, "summary.rds"))
runs <- lapply(paths, readRDS)
stopifnot(length(runs) == 10L,
          all(vapply(runs, function(x) x$status == "completed", logical(1))))
loco_record <- readRDS(file.path(root, "loco-seed2.rds"))
stopifnot(identical(loco_record$seed, runs[[2L]]$seed),
          loco_record$completed == loco_record$markers_tested,
          identical(names(loco_record$status), "ok"),
          all(loco_record$background$loco$status == ""))
loco <- data.frame(replicate = 2L, seed = loco_record$seed,
                   markers_tested = loco_record$markers_tested,
                   background_markers = loco_record$background_markers,
                   completed = loco_record$completed,
                   seconds = loco_record$seconds,
                   significant = loco_record$significant,
                   significant_causal = loco_record$significant_causal,
                   significant_zero_effect = loco_record$significant_zero_effect,
                   calibration_factor = loco_record$diagnostics$calibration_factor,
                   fitted_genetic_variance = loco_record$background$genetic_variance,
                   script_md5 = unname(loco_record$script_md5))
write.csv(loco, file.path(root, "loco.csv"), row.names = FALSE)
public_assets <- "website/assets/genomic-qualification"
dir.create(public_assets, recursive = TRUE, showWarnings = FALSE)
asset_names <- c("status.csv", "source-hashes.csv", "checks.csv", "accuracy.csv",
  "association.csv", "estimation.csv", "availability.csv", "bayesian.csv", "posterior.csv",
  "prediction.csv", "credible_sets.csv", "coloc.csv", "pathways.csv",
  "tails.csv", "loco.csv", "replicate-estimation.png",
  "replicate-prediction-mapping.png")
stopifnot(all(file.copy(file.path(root, asset_names),
                        file.path(public_assets, asset_names), overwrite = TRUE)))
fmt <- function(x, digits = 3L) sprintf(paste0("%.", digits, "f"), x)
tab <- function(headings, rows) c(
  paste0("| ", paste(headings, collapse = " | "), " |"),
  paste0("| ", paste(rep("---", length(headings)), collapse = " | "), " |"),
  vapply(rows, function(x) paste0("| ", paste(x, collapse = " | "), " |"),
         character(1)))
save_page <- function(name, lines) {
  while (length(lines) > 0L && !nzchar(tail(lines, 1L))) {
    lines <- head(lines, -1L)
  }
  writeLines(lines, file.path("docs", paste0(name, "-qualification.md")),
             useBytes = TRUE)
}
opening <- function(lib, guide) c(
  paste0("# ", lib, " qualification"), "",
  paste0("Ten fixed simulations use the [shared design](genomic-qualification-design.md). ",
         "The [", lib, " guide](", guide, ") describes inputs and options. ",
         "The [runnable workflow](../tools/examples/simulated_genomics_workflow.R) ",
         "contains the exact settings."), "")

association <- tables$association
checks <- summary$checks
gwas_time <- tables$timings$seconds[tables$timings$stage == "gwas"]
stopifnot(all(status$status == "completed"),
          all(checks$max_ols_error < 1e-9),
          all(checks$truth_score_error < 1e-9), length(gwas_time) == 10L)
save_page("glma", c(opening("glma", "methods.md#association-mapping"),
  "## Use", "",
  "```r", "linear <- glma(y, Glist, method = \"linear\", threads = 1)",
  "loco <- glma(y, chromosome_Glist, method = \"infinitesimal_loco\",",
  "             algorithm = \"observation_pcg\",",
  "             background_markers = seq(1, 50000, by = 100), threads = 2)",
  "```", "",
  "## Results", "",
  tab(c("Check", "Observed result"), list(
    c("Linear marker tests", format(sum(association$markers), big.mark = ",", scientific = FALSE)),
    c("Maximum selected OLS difference", format(max(checks$max_ols_error), scientific = TRUE, digits = 3)),
    c("Median time, three linear scans", paste0(fmt(median(gwas_time), 1), " seconds")),
    c("LOCO tests with successful status (seed 2)", paste0(loco$completed, "/", loco$markers_tested)),
    c("LOCO elapsed time (seed 2)", paste0(fmt(loco$seconds, 1), " seconds")),
    c("LOCO calibration factor (seed 2)", fmt(loco$calibration_factor))
  )), "",
  "The linear checks compare selected markers with R ordinary least squares.",
  "The separate LOCO timing run uses the second simulated dataset, groups",
  "artificial regions into 20 chromosomes and fits a 500-marker background.",
  "LOCO was timed once, so repeatability was not assessed. This unrelated-individual",
  "design does not establish mixed-model null-tail calibration, parameter recovery",
  "or a benefit over ordinary regression.", "",
  "[Association counts](../website/assets/genomic-qualification/association.csv) ·",
  "[LOCO timing record](../website/assets/genomic-qualification/loco.csv)", ""))

acc <- summary$accuracy
pick <- function(method, quantity, trait1, trait2 = NA_character_) {
  x <- acc[acc$analysis == method & acc$quantity == quantity &
    acc$trait1 == trait1 &
    (if(is.na(trait2)) is.na(acc$trait2) else acc$trait2 == trait2), ]
  stopifnot(nrow(x) == 1L)
  x
}
gc_methods <- c("ldsc", "gnova", "sumher", "ldsc_annotation")
gc_labels <- c("LDSC", "GNOVA", "SumHer", "Annotation LDSC")
gc_rows <- lapply(seq_along(gc_methods), function(j) {
  h <- pick(paste0(gc_methods[j], "_h2"), "h2", "A")
  r <- pick(paste0(gc_methods[j], "_rg"), "rg", "A", "B")
  c(gc_labels[j], fmt(h$mean_estimate), fmt(h$mean_truth), fmt(h$bias),
    fmt(r$mean_estimate), fmt(r$mean_truth), fmt(r$bias))
})
save_page("gcorr", c(opening("gcorr", "summary-statistic-analysis.md"),
  "## Use", "", "```r",
  "genome <- gcorr(stat, LDlist, method = \"ldsc\", task = \"rg\")",
  "partition <- gcorr(stat, annotated_LD, annotation = annotation,",
  "                   method = \"ldsc\", task = \"rg\")",
  "```", "", "## Results", "",
  "Trait A heritability and A–B genetic correlation are shown below. Other traits,",
  "regions and uncertainty availability remain in the compact result tables.", "",
  tab(c("Method", "h2 estimate", "h2 truth", "h2 bias", "rg estimate", "rg truth", "rg bias"), gc_rows), "",
  "![Heritability and genetic-correlation error across the fixed simulations.](../website/assets/genomic-qualification/replicate-estimation.png)", "",
  "These methods have different assumptions and uncertainty constructions. Ten",
  "replicates do not establish interval coverage or rank methods. Estimates that",
  "were unavailable or outside the correlation range are retained in the results.", "",
  "[Seed-level estimates](../website/assets/genomic-qualification/estimation.csv) ·",
  "[Accuracy summary](../website/assets/genomic-qualification/accuracy.csv) ·",
  "[Uncertainty availability](../website/assets/genomic-qualification/availability.csv)", ""))

b <- tables$bayesian
p <- tables$posterior
bayes_rows <- lapply(c("bayesc", "bayesr"), function(method) {
  h <- pick(method, "h2", "A")
  e <- b[b$method == method, ]
  q <- p[p$method == method & p$estimated & is.finite(p$rhat), ]
  c(method, fmt(h$mean_estimate), fmt(h$mean_truth), fmt(h$bias),
    fmt(h$rmse), paste0(h$covered, "/", h$intervals),
    fmt(mean(e$effect_correlation)), fmt(max(q$rhat)))
})
save_page("gbayes", c(opening("gbayes", "gbayes-workflow.qmd"),
  "## Use", "", "```r",
  "fit <- gbayes(prepared$A, LDlist, method = \"bayesc\", trait = \"A\",",
  "              prior = prior, control = list(residual_policy = \"fixed\",",
  "                                          burnin = 500, sampling_sweeps = 700,",
  "                                          seeds = c(11, 29, 47, 71)))",
  "```", "", "## Results", "",
  tab(c("Method", "h2 estimate", "h2 truth", "bias", "RMSE", "intervals", "effect r", "max R-hat"), bayes_rows), "",
  "The 50K-marker connected fits hold residual variance, effect scale and mixture",
  "weights fixed at declared values. They estimate marker effects and their",
  "posterior genetic summaries. [Separate evidence](simulated-genomics-qualification.md#bayesian-linear-regression-qualification)",
  "covers parameter learning with a matched independent LD reference. The",
  "interval counts above describe conditional posterior intervals, not calibrated",
  "frequentist coverage.", "",
  "[Effect recovery](../website/assets/genomic-qualification/bayesian.csv) ·",
  "[Posterior diagnostics](../website/assets/genomic-qualification/posterior.csv)", ""))

s <- tables$prediction
score_rows <- lapply(c("ridge", "bayesc", "bayesr"), function(method) {
  x <- s[s$method == method, ]
  c(method, fmt(mean(x$cor_genetic)), fmt(mean(x$cor_phenotype)),
    fmt(mean(x$calibration_slope)))
})
save_page("gscore", c(opening("gscore", "summary-interfaces.qmd"),
  "## Use", "", "```r",
  "weights <- gscore(prepared$A, LDlist, method = \"ridge\",",
  "                  control = list(penalty = 1))",
  "prediction <- gscore(task = \"apply\", Glist = held_out_Glist,",
  "                     weights = weights$weights,",
  "                     allele_frequencies = frequencies)",
  "```", "", "## Results", "",
  tab(c("Weights", "Mean r with genetic value", "Mean r with phenotype", "Mean calibration slope"), score_rows), "",
  paste0("Applying the true effects reproduced centered simulated genetic values ",
         "within ", format(max(checks$truth_score_error), scientific = TRUE, digits = 3),
         " at worst. This checks marker order and coding. Prediction quality ",
         "depends on the fitted weights and the deliberately difficult polygenic architecture."), "",
  "The ridge fit uses a fixed penalty of 1 without tuning. Its low held-out",
  "correlation and calibration slope are reported as observed; this setting is",
  "an execution check, not a recommended prediction recipe.", "",
  "![Held-out prediction, regional mapping and colocalisation across the fixed simulations.](../website/assets/genomic-qualification/replicate-prediction-mapping.png)", "",
  "[Held-out prediction results](../website/assets/genomic-qualification/prediction.csv)", ""))

sets <- tables$credible_sets
coloc <- tables$coloc
shared <- coloc[coloc$truth == "H4", ]
distinct <- coloc[coloc$truth == "H3", ]
save_page("gmap", c(opening("gmap", "gmap-workflow.qmd"),
  "## Use", "", "```r",
  "mapping_stat <- gmap_stat(stat, LDlist,",
  "                          phenotype_variance = phenotype_variance)",
  "fit <- gmap(mapping_stat[c(\"A\", \"B\")], LDlist,",
  "            regions = regions, method = \"multi_effect\",",
  "            prior = prior, control = control)",
  "shared <- gmap(stat[c(\"A\", \"B\")], task = \"coloc\",",
  "               method = \"abf\", trait1 = \"A\", trait2 = \"B\",",
  "               regions = regions[1])",
  "```", "", "## Results", "",
  tab(c("Check", "Ten-replicate result"), list(
    c("Component sets containing the causal marker", paste0(sum(sets$causal_covered), "/", nrow(sets))),
    c("Regional fits converged", paste0(sum(sets$converged), "/", nrow(sets))),
    c("Mean H4 in shared-signal region", fmt(mean(shared$H4))),
    c("Mean H3 in distinct-signal region", fmt(mean(distinct$H3)))
  )), "",
  "![Held-out prediction, regional mapping and colocalisation across the fixed simulations.](../website/assets/genomic-qualification/replicate-prediction-mapping.png)", "",
  "The first two regions were constructed for these hypotheses and contain one",
  "causal variant per trait. This is a focused regional check, not broad credible-set",
  "coverage or colocalisation calibration across realistic loci.", "",
  "[Component credible sets](../website/assets/genomic-qualification/credible_sets.csv) ·",
  "[Colocalisation probabilities](../website/assets/genomic-qualification/coloc.csv)", ""))

pathways <- tables$pathways
tails <- tables$tails
method_rows <- lapply(c("ora", "preranked", "competitive", "magma"), function(method) {
  x <- pathways[pathways$method == method, ]
  a <- x[x$pathway == "enriched", ]
  z <- x[x$pathway == "background", ]
  c(method, paste0(sum(a$p_value < .05), "/", nrow(a)),
    paste0(sum(z$p_value < .05), "/", nrow(z)))
})
bayes_pathway_rows <- lapply(c("enriched", "background", "mixed"), function(pathway) {
  x <- pathways[pathways$method == "bayesc" & pathways$pathway == pathway, ]
  stopifnot(nrow(x) == 10L, all(is.finite(x$pip)))
  c(pathway, fmt(mean(x$pip)))
})
save_page("gsea", c(opening("gsea", "gsea-workflow.qmd"),
  "## Use", "", "```r",
  "evidence <- gstat(stat, LDlist, sets = genes, blocks = blocks,",
  "                  metadata = metadata,",
  "                  control = list(independent_blocks = TRUE,",
  "                                 tail_method = \"saddlepoint\"))",
  "fit <- gsea(evidence$stat, pathways, method = \"competitive\",",
  "            sampling = evidence$sampling)",
  "```", "", "## Results", "",
  tab(c("Method", "Enriched p < 0.05", "Background p < 0.05"), method_rows), "",
  "The scalar BayesC pathway analysis reports posterior inclusion probabilities",
  "on a separate scale:", "",
  tab(c("Pathway", "Mean PIP"), bayes_pathway_rows), "",
  paste0("The controlled gene-tail calculation was available for ",
         sum(tails$controlled_available), "/", sum(tails$requested),
         " requested tests; the documented saddlepoint route was available for ",
         sum(tails$approximate_available), "/", sum(tails$requested),
         ". Its moment fallback was used ", sum(tails$fallback), " times."), "",
  "The enriched and background sets are artificial and fixed before analysis.",
  "Counts are descriptive across dependent traits and pathways; they are not",
  "power or false-positive-rate estimates. PIPs and p-values have different",
  "interpretations.", "",
  "[Pathway results](../website/assets/genomic-qualification/pathways.csv) ·",
  "[Gene-tail availability](../website/assets/genomic-qualification/tails.csv)", ""))
cat("GENOMIC_MODULE_REPORT|completed=", nrow(status),
    "|loco=", nrow(loco), "\n", sep = "")
