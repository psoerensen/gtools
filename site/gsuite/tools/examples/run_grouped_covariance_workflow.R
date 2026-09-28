# Deterministic standalone gcorr grouped SEM and GGM example.
# The source fit is deliberately small so this example tests the public
# grouped-covariance contract without rerunning a genome-wide estimator.
library(gcorr)

out <- getOption("gsuite.grouped.covariance.out",
  "build/examples/grouped-covariance-workflow")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

traits <- c("A", "B", "C")
groups <- c("coding", "regulatory")
make_covariance <- function(path) {
  directed <- matrix(0, 3, 3)
  directed[2, 1] <- path
  directed[3, 2] <- .25
  transform <- solve(diag(3) - directed)
  value <- transform %*% diag(c(1, .8, .7)) %*% t(transform)
  dimnames(value) <- list(traits, traits)
  value
}
covariance <- setNames(lapply(c(.35, .15), make_covariance), groups)
moment_labels <- c("A|A", "B|A", "C|A", "B|B", "C|B", "C|C")
labels <- unlist(lapply(groups, function(group) {
  paste(group, moment_labels, sep = "::")
}), use.names = FALSE)
sampling <- matrix(.00002, length(labels), length(labels),
  dimnames = list(labels, labels))
diag(sampling) <- .0004
parameters <- do.call(rbind, lapply(groups, function(group) {
  data.frame(
    quantity = c("category_heritability", "category_covariance",
      "category_covariance", "category_heritability",
      "category_covariance", "category_heritability"),
    trait1 = c("A", "B", "C", "B", "C", "C"),
    trait2 = c("", "A", "A", "", "B", ""), component = group,
    stringsAsFactors = FALSE)
}))
source_fit <- structure(list(method = "gnova", traits = traits,
  matrices = list(categories = lapply(covariance, function(value) {
    list(genetic_covariance = value)
  })), diagnostics = list(joint_jackknife = list(
    parameters = parameters, sampling_covariance = sampling))),
  class = "gcorr_gnova_fit")

input <- gcorr_grouped_covariance(source_fit, continuous = TRUE,
  units = setNames(rep("standardized genetic covariance", 3), traits))
model <- list(variables = traits,
  parameters = data.frame(name = c("ab", "bc", "vA", "vB", "vC"),
    start = c(.2, .2, 1, 1, 1)),
  entries = data.frame(matrix = c("directed", "directed",
    rep("disturbance", 3)), row = c("B", "C", traits),
    column = c("A", "B", traits),
    parameter = c("ab", "bc", "vA", "vB", "vC"), value = 0,
    stringsAsFactors = FALSE))
sem <- gcorr_grouped_sem(input, model)
ggm <- gcorr_grouped_ggm(input,
  matrix(c("A", "B", "B", "C"), 2, 2, byrow = TRUE))

sem_constraints <- matrix(0, 1, nrow(sem$estimates), dimnames = list(
  "coding minus regulatory ab", sem$estimates$label))
sem_constraints[1, c("coding::ab", "regulatory::ab")] <- c(1, -1)
sem_test <- gcorr_grouped_sem_test(sem, sem_constraints)
partial <- ggm$partial_correlations$estimates
ggm_constraints <- matrix(0, 1, nrow(partial), dimnames = list(
  "coding minus regulatory partial A-B", partial$label))
ggm_constraints[1, c("coding::partial[A,B]",
  "regulatory::partial[A,B]")] <- c(1, -1)
ggm_test <- gcorr_grouped_ggm_test(ggm, ggm_constraints)

stopifnot(sem$inference_available, ggm$inference_available,
  sem_test$available, ggm_test$available)
saveRDS(list(source = source_fit, covariance = input, sem = sem, ggm = ggm,
  sem_test = sem_test, ggm_test = ggm_test), file.path(out, "models.rds"))
write.csv(sem$estimates, file.path(out, "sem-estimates.csv"),
  row.names = FALSE)
write.csv(partial, file.path(out, "ggm-partial-correlations.csv"),
  row.names = FALSE)
print(sem_test)
print(ggm_test)
cat("GCORR_GROUPED_COVARIANCE_EXAMPLE|passed|",
  normalizePath(out, winslash = "/"), "\n", sep = "")
