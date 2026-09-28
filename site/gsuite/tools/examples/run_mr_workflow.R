suppressPackageStartupMessages({
  library(gbase)
  library(gcorr)
})

output <- file.path("build", "examples", "mr-workflow")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
markers <- paste0("mr", 1:12)
index <- 0:11
orientation <- ifelse(index %% 2, -1, 1)
exposure_effect <- orientation * (.06 + .01 * index)
outcome_effect <- orientation * (.3 * abs(exposure_effect) + .001 * sin(index))

make_summary <- function(trait, effect, se, n) {
  data.frame(rsids = markers, chr = "1", pos = seq_along(markers) * 100,
    ea = "A", nea = "G", b = effect, seb = se, n = n)
}

statistics <- list(X = make_summary("X", exposure_effect, .004, 100000),
  Y = make_summary("Y", outcome_effect, .006, 90000))
signed_ld <- diag(12)
dimnames(signed_ld) <- list(markers, markers)
prepared <- gcorr_prepare_mr(statistics, signed_ld, exposure = "X",
  outcome = "Y", reference_sample_size = 1024,
  ancestry = c(X = "artificial", Y = "artificial"),
  reference_ancestry = "artificial",
  effect_units = c(X = "standard deviations", Y = "standard deviations"),
  non_overlapping_samples = TRUE,
  preparation = "gsmr2")
fits <- list(gsmr = gcorr_gsmr(prepared), gsmr2 = gcorr_gsmr2(prepared),
  ivw = gcorr_mr_ivw(prepared), egger = gcorr_mr_egger(prepared),
  weighted_median = gcorr_mr_weighted_median(prepared,
    bootstrap_replicates = 1000, seed = 19))
estimates <- do.call(rbind, lapply(fits, `[[`, "estimates"))
saveRDS(list(statistics = statistics, signed_ld = signed_ld,
  prepared = prepared, fits = fits), file.path(output, "workflow.rds"))
write.csv(estimates, file.path(output, "estimates.csv"), row.names = FALSE)
cat("GCORR_MR_EXAMPLE|passed|", normalizePath(output, winslash = "/"),
  "\n", sep = "")
