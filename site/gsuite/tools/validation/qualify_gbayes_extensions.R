# Qualify advanced native gbayes transitions under their supported use cases.
#
# Raw fits and temporary genotype/LD resources remain under build/. The public
# record contains only compact summaries and one figure.

.libPaths(c(file.path("build", "r-library"), .libPaths()))
suppressPackageStartupMessages(library(gsuite))

args <- commandArgs(trailingOnly = TRUE)
reps_arg <- grep("^--reps=", args, value = TRUE)
repetitions <- if (length(reps_arg)) as.integer(sub("^--reps=", "", reps_arg[[1L]])) else 5L
if (!is.finite(repetitions) || repetitions < 1L)
  stop("--reps must be a positive integer.", call. = FALSE)

out <- getOption("gsuite.gbayes.extensions.out",
  file.path("build", "qualification", "gbayes-extensions"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, winslash = "/", mustWork = TRUE)

design <- list(seed = 20260922L, chains = 4L, threads = 4L,
  burnin = 300L, sweeps = 900L)

standardize <- function(x) {
  x <- sweep(x, 2L, colMeans(x), "-")
  sweep(x, 2L, sqrt(colSums(x^2) / nrow(x)), "/")
}

rmvn <- function(n, covariance) {
  matrix(rnorm(n * nrow(covariance)), n, nrow(covariance)) %*%
    chol(covariance)
}

write_bed <- function(W, prefix, marker, chromosome, position) {
  con <- file(paste0(prefix, ".bed"), "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  for (j in seq_len(ncol(W))) {
    codes <- c(3L, 2L, 0L)[W[, j] + 1L]
    codes <- c(codes, rep(1L, (-length(codes)) %% 4L))
    writeBin(as.raw(colSums(matrix(codes, 4L) * c(1, 4, 16, 64))), con)
  }
  write.table(data.frame(chromosome, marker, 0, position, "A", "G"),
    paste0(prefix, ".bim"), quote = FALSE, row.names = FALSE,
    col.names = FALSE)
  ids <- sprintf("i%05d", seq_len(nrow(W)))
  write.table(data.frame(ids, ids, 0, 0, 0, -9), paste0(prefix, ".fam"),
    quote = FALSE, row.names = FALSE, col.names = FALSE)
}

prepare_ld <- function(W, directory, marker, chromosome, position,
                       window = ncol(W), r2 = 0) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(directory, "genotypes")
  write_bed(W, prefix, marker, chromosome, position)
  gl <- gprep(bedfiles = paste0(prefix, ".bed"))
  ldprep(gl, reference = "simulated qualification reference",
    assembly = "artificial", task = "sparseld",
    out_prefix = file.path(directory, "LD"), max_distance_bp = 0,
    max_distance_variants = window, r2 = r2, block_size = 64L,
    nthreads = 1L, overwrite = TRUE)
}

stat_frame <- function(beta, n, LD) {
  md <- LD$resource$markers
  data.frame(marker = md$marker_id, allele1 = md$allele1,
    allele2 = md$allele2, beta_std = as.numeric(beta), n = n)
}

posterior_value <- function(fit, pattern) {
  i <- grep(pattern, fit$posterior$parameter)
  if (length(i) != 1L) stop("Expected one posterior value for ", pattern)
  fit$posterior$mean[[i]]
}

run_block <- function(rep, directory) {
  set.seed(design$seed + rep * 1000L + 11L)
  n <- 1400L; block_size <- 24L; p <- 2L * block_size
  marker <- sprintf("block%02d", seq_len(p))
  chromosome <- rep(c("1", "2"), each = block_size)
  position <- rep(seq_len(block_size) * 1000L, 2L)
  af <- runif(p, .15, .45)
  W <- sapply(af, function(q) rbinom(n, 2L, q))
  Z <- standardize(W)
  beta_hat <- numeric(p)
  truth <- c(left = .55, right = .82)
  response_ss <- numeric(2L)
  for (b in seq_len(2L)) {
    ix <- ((b - 1L) * block_size + 1L):(b * block_size)
    causal <- ix[seq(2L, block_size, length.out = 6L)]
    raw <- rnorm(length(causal))
    g <- drop(Z[, causal, drop = FALSE] %*% raw)
    g <- g * sqrt((1 - truth[[b]]) / mean(g^2))
    e <- rnorm(n, sd = sqrt(truth[[b]]))
    y <- g + e
    y <- (y - mean(y)) / sqrt(mean((y - mean(y))^2))
    truth[[b]] <- mean((e / sqrt(mean((g + e - mean(g + e))^2)))^2)
    beta_hat[ix] <- drop(crossprod(Z[, ix, drop = FALSE], y)) / n
    response_ss[[b]] <- sum(y^2)
  }
  LD <- prepare_ld(W, file.path(directory, "block"), marker, chromosome,
    position, window = block_size, r2 = 0)
  blocks <- setNames(rep(c("left", "right"), each = block_size), marker)
  prior <- list(residual_variance = 1, effect_variance = .015,
    inclusion_probability = .25)
  fit <- gbayes(stat_frame(beta_hat, n, LD), LD, "bayesc", prior = prior,
    control = list(burnin = design$burnin, sampling_sweeps = design$sweeps,
      seeds = design$seed + rep * 100L + seq_len(design$chains),
      threads = design$threads,
      block_residual = list(block = blocks, response_ss = response_ss,
        likelihood_dimension = rep(n, 2L), initial_variance = rep(.7, 2L),
        genetic_variance_denominator = rep(n, 2L), prior_df = 4,
        prior_scale = rep(.35, 2L), mode = "sample",
        minimum_variance_ratio = 1e-6)))
  estimate <- fit$block_residual_variance_mean
  data.frame(repetition = rep, scenario = "Block residual",
    metric = names(truth), truth = as.numeric(truth), estimate = estimate,
    error = estimate - truth, auxiliary = NA_real_)
}

run_bayess <- function(rep, directory) {
  set.seed(design$seed + rep * 1000L + 21L)
  n <- 4000L; p <- 600L; causal_n <- 300L; true_s <- -.50
  marker <- sprintf("maf%03d", seq_len(p))
  af <- seq(.05, .48, length.out = p)[sample.int(p)]
  W <- sapply(af, function(q) rbinom(n, 2L, q))
  Z <- standardize(W)
  causal <- sample.int(p, causal_n)
  q <- 2 * af * (1 - af)
  beta <- numeric(p)
  beta[causal] <- rnorm(causal_n, sd = sqrt(q[causal]^true_s))
  g <- drop(Z %*% beta); beta <- beta * sqrt(.60 / mean(g^2)); g <- drop(Z %*% beta)
  y <- g + rnorm(n, sd = sqrt(.40)); y <- as.numeric(scale(y))
  LD <- prepare_ld(W, file.path(directory, "bayess"), marker, "1",
    seq_len(p) * 1000L, window = p, r2 = .02)
  stat <- stat_frame(drop(crossprod(Z, y)) / n, n, LD)
  af_named <- setNames(af, marker)
  fit <- gbayes(stat, LD, "bayesc",
    prior = list(residual_variance = .40, effect_variance = .004,
      inclusion_probability = causal_n / p,
      effect_variance_prior = list(df = 4, scale = .002)),
    control = list(burnin = design$burnin, sampling_sweeps = design$sweeps,
      seeds = design$seed + rep * 100L + 20L + seq_len(design$chains),
      threads = design$threads, estimate_effect_variance = TRUE),
    annotation = list(bayess = list(allele_frequency = af_named,
      coefficient_prior_sd = 1.5)))
  estimate <- unname(fit$annotations$variance$coefficient_mean[["maf_exponent_s"]])
  data.frame(repetition = rep, scenario = "BayesS exponent",
    metric = "S", truth = true_s, estimate = estimate,
    error = estimate - true_s, auxiliary = NA_real_)
}

run_multivariate <- function(rep, directory) {
  set.seed(design$seed + rep * 1000L + 31L)
  n <- 1800L; p <- 90L; traits <- c("trait1", "trait2", "trait3")
  marker <- sprintf("multi%03d", seq_len(p))
  af <- runif(p, .12, .45)
  W <- sapply(af, function(q) rbinom(n, 2L, q)); Z <- standardize(W)
  active <- sample.int(p, 72L)
  base_cor <- matrix(c(1,.45,.20,.45,1,.35,.20,.35,1), 3L,
    dimnames = list(traits, traits))
  base_sd <- c(.055, .050, .045)
  V <- diag(base_sd) %*% base_cor %*% diag(base_sd)
  B <- matrix(0, p, 3L); B[active, ] <- rmvn(length(active), V)
  G <- Z %*% B
  omega <- matrix(c(.62,.12,.05,.12,.68,.09,.05,.09,.72), 3L,
    dimnames = list(traits, traits))
  Y0 <- G + rmvn(n, omega)
  sy <- sqrt(colMeans(sweep(Y0, 2L, colMeans(Y0), "-")^2))
  Y <- sweep(sweep(Y0, 2L, colMeans(Y0), "-"), 2L, sy, "/")
  B <- sweep(B, 2L, sy, "/")
  truth_V <- cov(B[active, , drop = FALSE])
  truth_omega <- diag(1 / sy) %*% omega %*% diag(1 / sy)
  dimnames(truth_V) <- dimnames(truth_omega) <- list(traits, traits)
  LD <- prepare_ld(W, file.path(directory, "multivariate"), marker, "1",
    seq_len(p) * 1000L, window = p, r2 = 0)
  bhat <- crossprod(Z, Y) / n
  stat <- setNames(lapply(seq_along(traits), function(k)
    stat_frame(bhat[, k], n, LD)), traits)
  patterns <- rbind(none = c(0,0,0), shared = c(1,1,1)); colnames(patterns) <- traits
  initial_V <- truth_V + diag(1e-6, 3L)
  fit <- gbayes(stat, LD, "bayesc",
    prior = list(effect_covariance = initial_V, patterns = patterns,
      pattern_weights = c(.2, .8),
      covariance_prior = list(df = 7, scale = 3 * initial_V),
      sampling_covariance_prior = list(df = 7, scale = 3 * truth_omega)),
    sampling = list(dependence = "shared_ld", covariance = truth_omega,
      response_crossproduct = crossprod(Y), likelihood_dimension = n),
    control = list(burnin = design$burnin, sampling_sweeps = design$sweeps,
      seeds = design$seed + rep * 100L + 40L + seq_len(design$chains),
      threads = design$threads, estimate_covariance = TRUE,
      estimate_sampling_covariance = TRUE))
  estimated_V <- fit$parameter_mean$effect_covariance
  estimated_omega <- fit$parameter_mean$sampling_covariance
  pick <- lower.tri(truth_V, diag = TRUE)
  rbind(
    data.frame(repetition = rep, scenario = "Effect covariance",
      metric = paste(row(truth_V)[pick], col(truth_V)[pick], sep = ","),
      truth = truth_V[pick], estimate = estimated_V[pick],
      error = estimated_V[pick] - truth_V[pick], auxiliary = NA_real_),
    data.frame(repetition = rep, scenario = "Sampling covariance",
      metric = paste(row(truth_omega)[pick], col(truth_omega)[pick], sep = ","),
      truth = truth_omega[pick], estimate = estimated_omega[pick],
      error = estimated_omega[pick] - truth_omega[pick], auxiliary = NA_real_))
}

run_swap <- function(rep, directory) {
  set.seed(design$seed + rep * 1000L + 41L)
  n <- 1800L; pairs <- 20L; p <- pairs * 2L
  marker <- sprintf("swap%02d", seq_len(p)); chromosome <- rep(seq_len(pairs), each = 2L)
  af <- runif(pairs, .15, .4)
  W <- matrix(0L, n, p)
  for (k in seq_len(pairs)) {
    lead <- rbinom(n, 2L, af[[k]])
    proxy <- ifelse(runif(n) < .97, lead, rbinom(n, 2L, af[[k]]))
    W[, 2L*k-c(1L,0L)] <- cbind(lead, proxy)
  }
  Z <- standardize(W); causal <- seq(1L, p, by = 2L)
  beta <- numeric(p); beta[causal] <- rnorm(pairs)
  g <- drop(Z %*% beta); beta <- beta * sqrt(.45 / mean(g^2)); g <- drop(Z %*% beta)
  y <- g + rnorm(n, sd = sqrt(.55)); y <- as.numeric(scale(y))
  LD <- prepare_ld(W, file.path(directory, "swap"), marker, chromosome,
    rep(c(1000L, 2000L), pairs), window = 2L, r2 = 0)
  stat <- stat_frame(drop(crossprod(Z, y)) / n, n, LD)
  prior <- list(residual_variance = .55, effect_variance = .45 / pairs,
    inclusion_probability = .5)
  common <- list(burnin = design$burnin, sampling_sweeps = design$sweeps,
    seeds = design$seed + rep * 100L + 60L + seq_len(design$chains),
    threads = design$threads)
  ordinary <- gbayes(stat, LD, "bayesc", prior = prior, control = common)
  moved <- gbayes(stat, LD, "bayesc", prior = prior,
    control = c(common, list(ld_swap = list(probability = .5,
      minimum_r2 = .8, maximum_friends = 4L, moves = 2L))))
  score <- function(fit) cor(drop(Z %*% fit$estimates$mean), g)
  data.frame(repetition = rep, scenario = "LD relocation",
    metric = "score correlation", truth = score(ordinary), estimate = score(moved),
    error = score(moved) - score(ordinary),
    auxiliary = moved$ld_swap$accepted / moved$ld_swap$attempted)
}

run_one <- function(rep) {
  directory <- file.path(out, sprintf("rep%02d", rep))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  rbind(run_block(rep, directory), run_bayess(rep, directory),
    run_multivariate(rep, directory), run_swap(rep, directory))
}

results <- vector("list", repetitions)
for (rep in seq_len(repetitions)) {
  cat(sprintf("START|rep%02d\n", rep)); flush.console()
  results[[rep]] <- run_one(rep)
  cat(sprintf("DONE|rep%02d\n", rep)); flush.console()
}
results <- do.call(rbind, results); rownames(results) <- NULL
write.csv(results, file.path(out, "replicate-results.csv"), row.names = FALSE)

summary <- do.call(rbind, lapply(split(results, results$scenario), function(x) {
  scale <- pmax(abs(x$truth), .05)
  data.frame(scenario = x$scenario[[1L]], repetitions = length(unique(x$repetition)),
    mean_truth = mean(x$truth), mean_estimate = mean(x$estimate),
    mean_absolute_error = mean(abs(x$error)),
    mean_scaled_absolute_error = mean(abs(x$error) / scale),
    minimum_auxiliary = if (all(is.na(x$auxiliary))) NA_real_ else min(x$auxiliary, na.rm = TRUE))
}))
rownames(summary) <- NULL
summary$passed <- with(summary,
  (scenario == "Block residual" & mean_absolute_error <= .05) |
  (scenario == "BayesS exponent" & mean_absolute_error <= .25) |
  (scenario == "Effect covariance" & mean_absolute_error <= 5e-4) |
  (scenario == "Sampling covariance" & mean_absolute_error <= .05) |
  (scenario == "LD relocation" & mean_absolute_error <= .02 & minimum_auxiliary > 0))
write.csv(summary, file.path(out, "qualification-summary.csv"), row.names = FALSE)
criteria <- data.frame(scenario = c("Block residual", "BayesS exponent",
    "Effect covariance", "Sampling covariance", "LD relocation"),
  criterion = c("mean absolute error <= 0.05",
    "mean absolute error <= 0.25",
    "mean absolute error <= 0.0005",
    "mean absolute error <= 0.05",
    "mean score-correlation change <= 0.02 and accepted fraction > 0"))
write.csv(criteria, file.path(out, "acceptance-criteria.csv"), row.names = FALSE)

png(file.path(out, "qualification.png"), width = 1800, height = 900, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.4, 4.7, 2.3, 1), bty = "l", las = 1)
parameter <- results[results$scenario %in% c("Block residual", "BayesS exponent"), ]
plot(parameter$truth, parameter$estimate, pch = ifelse(parameter$scenario == "Block residual", 19, 1),
  col = ifelse(parameter$scenario == "Block residual", "#2f6f99", "#c46a2d"),
  xlab = "Simulation truth", ylab = "Posterior mean", main = "Scale and BayesS recovery")
abline(0, 1, lty = 2, col = "#666666")
legend("topleft", c("Block residual", "BayesS exponent"), pch = c(19,1),
  col = c("#2f6f99", "#c46a2d"), bty = "n")
covariance <- results[results$scenario %in% c("Effect covariance", "Sampling covariance"), ]
plot(covariance$truth, covariance$estimate,
  pch = ifelse(covariance$scenario == "Effect covariance", 19, 1),
  col = ifelse(covariance$scenario == "Effect covariance", "#2f6f99", "#31845f"),
  xlab = "Simulation truth", ylab = "Posterior mean", main = "Multivariate covariance recovery")
abline(0, 1, lty = 2, col = "#666666")
legend("topleft", c("Effect covariance", "Sampling covariance"), pch = c(19,1),
  col = c("#2f6f99", "#31845f"), bty = "n")
par(op); dev.off()

saveRDS(list(design = design, criteria = criteria, results = results, summary = summary,
  session = sessionInfo()), file.path(out, "qualification-record.rds"),
  compress = FALSE)
print(summary)
if (repetitions >= 5L && !all(summary$passed))
  stop("The advanced gbayes qualification did not meet every acceptance criterion.",
    call. = FALSE)
