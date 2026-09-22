# Qualify the supported scalar reference-LD BayesC/BayesR workflow.
#
# The public result is a supported-use contract, not a sensitivity study. Raw
# repetitions and temporary genotype/LD resources remain under build/.

.libPaths(c(file.path("build", "r-library"), .libPaths()))
suppressPackageStartupMessages(library(gsuite))
stopifnot(requireNamespace("gsim", quietly = TRUE))

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else
  file.path("tools", "validation", "qualify_gbayes_reference_ld.R")
script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
resume <- "--resume" %in% args
reps_arg <- grep("^--reps=", args, value = TRUE)
repetitions <- if (length(reps_arg)) as.integer(sub("^--reps=", "", reps_arg[[1L]])) else 5L
if (!is.finite(repetitions) || repetitions < 1L)
  stop("--reps must be a positive integer.", call. = FALSE)

out <- getOption("gsuite.gbayes.qualification.out",
  file.path("build", "qualification", "gbayes-reference-ld"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, winslash = "/", mustWork = TRUE)

design <- list(
  seed = 20260922L,
  n_study = 3000L,
  n_reference = 9000L,
  markers = 2000L,
  causal_markers = 40L,
  target_h2 = 0.30,
  copy_probability = 0.65,
  ld_window_markers = 1000L,
  ld_r2_threshold = 0.001,
  residual_adjustment = 0.9,
  weight_prior_strength = 5000,
  burnin = 500L,
  bayesc_sweeps = 700L,
  bayesr_sweeps = 5000L,
  chains = 4L,
  threads = 4L
)

marker_ids <- sprintf("m%04d", seq_len(design$markers))
active_probability <- design$causal_markers / design$markers
bayesr_active_counts <- c(24L, 12L, 4L)
bayesr_multipliers <- c(0, 0.01, 0.1, 1)
bayesr_weights <- c(1 - active_probability,
  bayesr_active_counts / design$markers)

simulate_genotypes <- function(n, seed, row_prefix) {
  set.seed(seed)
  W <- matrix(0L, n, design$markers,
    dimnames = list(sprintf("%s%05d", row_prefix, seq_len(n)), marker_ids))
  W[, 1L] <- rbinom(n, 2L, 0.30)
  for (j in 2:design$markers) {
    fresh <- rbinom(n, 2L, 0.30)
    copy <- runif(n) < design$copy_probability
    W[, j] <- ifelse(copy, W[, j - 1L], fresh)
  }
  W
}

standardize <- function(W) {
  centered <- sweep(W, 2L, colMeans(W), "-")
  scale <- sqrt(colSums(centered^2) / nrow(W))
  if (any(!is.finite(scale)) || any(scale <= 0))
    stop("All simulated markers must be polymorphic.", call. = FALSE)
  list(value = sweep(centered, 2L, scale, "/"), scale = scale)
}

write_bed <- function(W, prefix) {
  con <- file(paste0(prefix, ".bed"), "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  for (j in seq_len(ncol(W))) {
    codes <- c(3L, 2L, 0L)[W[, j] + 1L]
    codes <- c(codes, rep(1L, (-length(codes)) %% 4L))
    writeBin(as.raw(colSums(matrix(codes, 4L) * c(1, 4, 16, 64))), con)
  }
  write.table(data.frame(1L, marker_ids, 0,
    seq_len(ncol(W)) * 1000L, "A", "G"), paste0(prefix, ".bim"),
    quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(data.frame(rownames(W), rownames(W), 0, 0, 0, -9),
    paste0(prefix, ".fam"), quote = FALSE, row.names = FALSE,
    col.names = FALSE)
}

make_effects <- function(seed) {
  set.seed(seed)
  B <- matrix(0, design$markers, 2L,
    dimnames = list(marker_ids, c("BayesC", "BayesR")))
  causal_c <- sort(sample(seq_len(design$markers), design$causal_markers))
  remaining <- setdiff(seq_len(design$markers), causal_c)
  causal_r <- sort(sample(remaining, design$causal_markers))
  B[causal_c, "BayesC"] <- rnorm(design$causal_markers)
  component <- rep(seq_along(bayesr_active_counts), bayesr_active_counts)
  component <- sample(component)
  B[causal_r, "BayesR"] <- rnorm(design$causal_markers,
    sd = sqrt(bayesr_multipliers[component + 1L]))
  list(beta = B, causal = list(BayesC = causal_c, BayesR = causal_r),
    component = component)
}

method_prior <- function(method) {
  expected_multiplier <- if (method == "bayesc") active_probability else
    sum(bayesr_weights * bayesr_multipliers)
  vb <- design$target_h2 / (design$markers * expected_multiplier)
  prior <- list(
    residual_variance = 1 - design$target_h2,
    effect_variance = vb,
    effect_variance_prior = list(df = 4, scale = 0.5 * vb),
    residual_variance_prior = list(df = 4,
      scale = 0.5 * (1 - design$target_h2))
  )
  if (method == "bayesc") {
    prior$inclusion_probability <- active_probability
    prior$weight_prior <- design$weight_prior_strength * c(1 - active_probability,
      active_probability)
  } else {
    prior$weights <- bayesr_weights
    prior$variance_multipliers <- bayesr_multipliers
    prior$weight_prior <- design$weight_prior_strength * bayesr_weights
  }
  prior
}

posterior_row <- function(fit, pattern) {
  hit <- grep(pattern, fit$posterior$parameter)
  if (length(hit) != 1L)
    stop("Expected one posterior row matching: ", pattern, call. = FALSE)
  fit$posterior[hit, , drop = FALSE]
}

fit_summary <- function(fit, method, repetition, truth_effect,
                        truth_score, truth_reference_h2, elapsed) {
  h2 <- posterior_row(fit, "^h2\\[")
  estimated <- fit$posterior$estimated & fit$posterior$status == "available"
  rhat <- fit$posterior$rhat[estimated]
  ess <- fit$posterior$ess_bulk[estimated]
  effect <- setNames(fit$estimates$mean, fit$estimates$marker)[marker_ids]
  predicted <- drop(truth_score$Z %*% effect)
  weights <- as.numeric(fit$parameter_mean$weights)
  active <- posterior_row(fit, "^active_markers$")$mean[[1L]]
  data.frame(
    repetition = repetition,
    method = if (method == "bayesc") "BayesC" else "BayesR",
    target_h2 = design$target_h2,
    reference_truth_h2 = truth_reference_h2,
    h2_mean = h2$mean,
    h2_lower = h2$lower,
    h2_upper = h2$upper,
    interval_contains_truth = h2$lower <= truth_reference_h2 &&
      h2$upper >= truth_reference_h2,
    residual_variance = as.numeric(fit$parameter_mean$residual_variance),
    effect_variance = as.numeric(fit$parameter_mean$effect_variance),
    active_probability = 1 - weights[[1L]],
    active_markers = active,
    effect_correlation = cor(effect, truth_effect),
    score_correlation = cor(predicted, truth_score$value),
    score_calibration = unname(coef(lm(truth_score$value ~ 0 + predicted))[[1L]]),
    max_rhat = if (any(is.finite(rhat))) max(rhat, na.rm = TRUE) else NA_real_,
    min_ess_bulk = if (any(is.finite(ess))) min(ess, na.rm = TRUE) else NA_real_,
    elapsed_seconds = elapsed,
    stringsAsFactors = FALSE
  )
}

run_repetition <- function(repetition) {
  rep_dir <- file.path(out, sprintf("rep%02d", repetition))
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(rep_dir, "result.rds")
  if (resume && file.exists(result_path)) {
    saved <- readRDS(result_path)
    if (!identical(saved$design, design))
      stop("Saved repetition does not match the current qualification design: ",
        result_path, call. = FALSE)
    return(saved)
  }

  rep_seed <- design$seed + 10000L * repetition
  W <- simulate_genotypes(design$n_study, rep_seed + 1L, "study")
  W_reference <- simulate_genotypes(design$n_reference, rep_seed + 2L, "ref")
  effects <- make_effects(rep_seed + 3L)
  simulation <- gsim::gsim(W = W, nt = 2L, architecture = "fixed",
    beta = effects$beta, h2 = rep(design$target_h2, 2L),
    vg = rep(design$target_h2, 2L), re = 0, standardize_W = FALSE,
    scale_effects = TRUE, seed = rep_seed + 4L, compute_sumstats = FALSE)
  colnames(simulation$Y) <- colnames(simulation$G) <-
    colnames(simulation$B) <- c("BayesC", "BayesR")

  study_prefix <- file.path(rep_dir, "study-genotypes")
  reference_prefix <- file.path(rep_dir, "reference-genotypes")
  write_bed(W, study_prefix)
  write_bed(W_reference, reference_prefix)
  study_glist <- gprep(bedfiles = paste0(study_prefix, ".bed"))
  reference_glist <- gprep(bedfiles = paste0(reference_prefix, ".bed"))
  LD <- ldprep(reference_glist, reference = "matched-independent-reference",
    assembly = "artificial", task = "sparseld",
    out_prefix = file.path(rep_dir, "LD"), max_distance_bp = 0,
    max_distance_variants = design$ld_window_markers,
    r2 = design$ld_r2_threshold, block_size = 128L,
    nthreads = 1L, overwrite = TRUE)
  validate_LDlist(LD)

  study_standardized <- standardize(W)
  reference_standardized <- standardize(W_reference)
  summaries <- vector("list", 2L)
  fits <- vector("list", 2L)
  names(fits) <- c("BayesC", "BayesR")

  for (trait in c("BayesC", "BayesR")) {
    method <- tolower(trait)
    y <- simulation$Y[, trait]
    names(y) <- rownames(W)
    scan <- glma(y, study_glist, method = "linear", threads = 1L,
      block_size = 128L)
    associations <- scan$associations
    stat <- data.frame(marker = associations$marker, allele1 = "A",
      allele2 = "G", chromosome = "1",
      position_bp = seq_len(design$markers) * 1000L,
      beta = associations$beta, se = associations$se,
      n = design$n_study, p_value = associations$p)
    prepared <- sumstat(setNames(list(stat), trait), LD,
      task = "standardize")$stat[[trait]]

    fit_seed <- rep_seed + if (method == "bayesc") 100L else 200L
    control <- list(
      burnin = design$burnin,
      sampling_sweeps = if (method == "bayesc") design$bayesc_sweeps else
        design$bayesr_sweeps,
      seeds = fit_seed + seq_len(design$chains),
      threads = design$threads,
      residual_policy = "reference_ld",
      residual_adjustment = design$residual_adjustment,
      estimate_effect_variance = TRUE,
      estimate_weights = TRUE
    )
    started <- proc.time()[["elapsed"]]
    fit <- gbayes(prepared, LD, method = method, trait = trait,
      prior = method_prior(method), control = control,
      posterior = list(reference = LD,
        phenotype_variance = setNames(1, trait)))
    elapsed <- proc.time()[["elapsed"]] - started
    fits[[trait]] <- fit

    y_scale <- sqrt(sum((y - mean(y))^2) / length(y))
    truth_effect <- simulation$B[, trait] * study_standardized$scale / y_scale
    truth_score <- list(Z = study_standardized$value,
      value = drop(study_standardized$value %*% truth_effect))
    reference_score <- drop(reference_standardized$value %*% truth_effect)
    truth_reference_h2 <- mean(reference_score^2)
    summaries[[match(trait, c("BayesC", "BayesR"))]] <- fit_summary(
      fit, method, repetition, truth_effect, truth_score,
      truth_reference_h2, elapsed)
  }

  result <- list(summary = do.call(rbind, summaries), fits = fits,
    truth = list(beta = simulation$B, genetic = simulation$G,
      causal = effects$causal), design = design,
    ld = list(stored_entries = LD$resource$stored_entry_count,
      marker_count = LD$resource$marker_count),
    versions = c(gsuite = as.character(packageVersion("gsuite")),
      gsim = as.character(packageVersion("gsim"))), session = sessionInfo())
  saveRDS(result, result_path, compress = FALSE)
  result
}

results <- vector("list", repetitions)
for (rep in seq_len(repetitions)) {
  cat(sprintf("START|rep%02d\n", rep)); flush.console()
  results[[rep]] <- run_repetition(rep)
  cat(sprintf("DONE|rep%02d\n", rep)); flush.console()
}

replicate_results <- do.call(rbind, lapply(results, `[[`, "summary"))
rownames(replicate_results) <- NULL
write.csv(replicate_results, file.path(out, "replicate-results.csv"),
  row.names = FALSE)

aggregate_one <- function(x) {
  error <- x$h2_mean - x$reference_truth_h2
  data.frame(
    method = x$method[[1L]],
    repetitions = nrow(x),
    mean_reference_truth_h2 = mean(x$reference_truth_h2),
    mean_h2 = mean(x$h2_mean),
    h2_bias = mean(error),
    h2_rmse = sqrt(mean(error^2)),
    interval_coverage = mean(x$interval_contains_truth),
    mean_effect_correlation = mean(x$effect_correlation),
    mean_score_correlation = mean(x$score_correlation),
    mean_score_calibration = mean(x$score_calibration),
    max_rhat = max(x$max_rhat),
    min_ess_bulk = min(x$min_ess_bulk),
    mean_elapsed_seconds = mean(x$elapsed_seconds),
    stringsAsFactors = FALSE
  )
}
qualification_summary <- do.call(rbind,
  lapply(split(replicate_results, replicate_results$method), aggregate_one))
rownames(qualification_summary) <- NULL

criteria <- data.frame(
  criterion = c("absolute h2 bias", "h2 RMSE", "interval coverage",
    "mean effect correlation", "mean score correlation", "maximum R-hat",
    "minimum bulk ESS"),
  rule = c("<= 0.05", "<= 0.075", ">= 0.80", ">= 0.80", ">= 0.90",
    "<= 1.05", ">= 100"), stringsAsFactors = FALSE)
qualification_summary$pass_h2_bias <- abs(qualification_summary$h2_bias) <= 0.05
qualification_summary$pass_h2_rmse <- qualification_summary$h2_rmse <= 0.075
qualification_summary$pass_coverage <- qualification_summary$interval_coverage >= 0.80
qualification_summary$pass_effect_correlation <-
  qualification_summary$mean_effect_correlation >= 0.80
qualification_summary$pass_score_correlation <-
  qualification_summary$mean_score_correlation >= 0.90
qualification_summary$pass_rhat <- qualification_summary$max_rhat <= 1.05
qualification_summary$pass_ess <- qualification_summary$min_ess_bulk >= 100
pass_columns <- grep("^pass_", names(qualification_summary), value = TRUE)
qualification_summary$passed <- apply(qualification_summary[pass_columns], 1L, all)

write.csv(qualification_summary, file.path(out, "qualification-summary.csv"),
  row.names = FALSE)
write.csv(criteria, file.path(out, "acceptance-criteria.csv"), row.names = FALSE)

png(file.path(out, "qualification.png"), width = 1800, height = 760,
  res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.8, 2.2, 1.0),
  las = 1, bty = "l")
palette <- c(BayesC = "#2f6f99", BayesR = "#c46a2d")
offset <- c(BayesC = -0.08, BayesR = 0.08)
plot(NA, xlim = c(0.65, repetitions + 0.35), ylim = range(c(
  replicate_results$h2_lower, replicate_results$h2_upper,
  replicate_results$reference_truth_h2)), xlab = "Simulation",
  ylab = expression("Reference " * h^2), xaxt = "n",
  main = "Heritability recovery")
axis(1, seq_len(repetitions))
for (method in names(palette)) {
  z <- replicate_results[replicate_results$method == method, ]
  x <- z$repetition + offset[[method]]
  arrows(x, z$h2_lower, x, z$h2_upper, angle = 90, code = 3,
    length = 0.035, col = palette[[method]], lwd = 1.5)
  points(x, z$h2_mean, pch = 19, col = palette[[method]], cex = 1.05)
  points(x, z$reference_truth_h2, pch = 1, col = palette[[method]],
    cex = 1.05, lwd = 1.4)
}
legend("topright", c("BayesC estimate", "BayesR estimate", "Reference truth"),
  col = c(palette, "#333333"), pch = c(19, 19, 1), bty = "n", cex = 0.86)

plot(NA, xlim = c(0.65, repetitions + 0.35), ylim = c(0.85, 1),
  xlab = "Simulation", ylab = "Correlation with truth", xaxt = "n",
  main = "Effects and genetic scores")
axis(1, seq_len(repetitions))
abline(h = c(0.8, 0.9), lty = 3, col = "#aaaaaa")
for (method in names(palette)) {
  z <- replicate_results[replicate_results$method == method, ]
  x <- z$repetition + offset[[method]]
  points(x, z$effect_correlation, pch = 1, col = palette[[method]],
    cex = 1.05, lwd = 1.4)
  points(x, z$score_correlation, pch = 19, col = palette[[method]],
    cex = 1.05)
}
legend("bottomright", c("Effect correlation", "Score correlation"),
  pch = c(1, 19), col = "#333333", bty = "n", cex = 0.86)
par(op)
dev.off()

saveRDS(list(design = design, criteria = criteria,
  qualification_summary = qualification_summary,
  replicate_results = replicate_results,
  script_md5 = unname(tools::md5sum(script_path)), session = sessionInfo()),
  file.path(out, "qualification-record.rds"), compress = FALSE)

print(qualification_summary)
if (repetitions >= 5L && !all(qualification_summary$passed))
  stop("The supported reference-LD qualification did not meet every acceptance criterion.",
    call. = FALSE)
