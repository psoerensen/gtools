suppressPackageStartupMessages(library(gcorr))

output <- file.path("build", "examples", "covariance-workflow")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

traits <- LETTERS[1:4]
loadings <- c(.8, .7, .6, .5)
unique_variance <- c(.4, .5, .6, .7)
G <- tcrossprod(loadings) + diag(unique_variance)
dimnames(G) <- list(traits, traits)
labels <- unlist(lapply(seq_along(traits), function(column) {
  paste(traits[column:length(traits)], traits[column], sep = "|")
}), use.names = FALSE)
V <- .0004 * .3^abs(outer(seq_along(labels), seq_along(labels), "-"))
dimnames(V) <- list(labels, labels)
units <- setNames(rep("standardized units", 4), traits)
input <- gcorr_covariance(G, V, continuous = TRUE, units = units,
  provenance = "Constructed one-factor covariance")

model <- list(variables = c(traits, "F"),
  parameters = data.frame(name = c(paste0("l", traits), paste0("v", traits)),
    start = c(rep(.5, 4), rep(.6, 4)),
    lower = c(rep(-Inf, 4), rep(1e-8, 4)), upper = rep(Inf, 8)),
  entries = data.frame(matrix = c(rep("directed", 4), rep("disturbance", 5)),
    row = c(traits, traits, "F"), column = c(rep("F", 4), traits, "F"),
    parameter = c(paste0("l", traits), paste0("v", traits), ""),
    value = c(rep(0, 8), 1)))
sem <- gcorr_sem(input, model)
network <- gcorr_ggm_network(input)
edges <- t(combn(traits, 2))
ggm <- gcorr_ggm(input, edges)

constraints <- matrix(0, 2, nrow(sem$estimates),
  dimnames = list(c("equal first loadings", "third loading"),
    sem$estimates$parameter))
constraints[1, c("lA", "lB")] <- c(1, -1)
constraints[2, "lC"] <- 1
sem_test <- gcorr_sem_test(sem, constraints, target = c(0, .5))
edge_test <- gcorr_ggm_test(ggm,
  matrix(c("A", "B", "B", "C"), 2, 2, byrow = TRUE))

saveRDS(list(input = input, sem = sem, network = network, ggm = ggm,
  sem_test = sem_test, edge_test = edge_test),
  file.path(output, "workflow.rds"))
cat("GCORR_COVARIANCE_EXAMPLE|passed|", normalizePath(output,
  winslash = "/"), "\n", sep = "")
