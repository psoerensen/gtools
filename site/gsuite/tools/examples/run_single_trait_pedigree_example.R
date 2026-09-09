# Run with an installed gsuite: source("tools/examples/run_single_trait_pedigree_example.R")
library(gsuite)

# Reuse the BW illustration in tests/testthat/test-formula-reml.R.
# These deterministic toy values are NOT draws from a pedigree animal model.
# In particular, its artificial sex labels vary between records of an animal;
# real animal records should carry the animal's consistent sex.
pedigree <- data.frame(
  id = c("a4", "a1", "a7", "a2", "a5", "a3", "a6", "a8"),
  sire = c("a1", NA, "a3", NA, "a2", "a1", "a3", "a4"),
  dam = c("a2", NA, "a5", NA, "a1", "a2", "a5", "a6"))
row <- seq_len(35L)
id <- rep(c("a5", "a1", "a7", "a3", "a6", "a2", "a4"), 5L)
sex <- factor(rep(c("F", "M", "F", "M", "F"), each = 7L),
              levels = c("F", "M"))
x <- ((row * 7L) %% 19L - 9) / 7
animal_bw <- c(a1 = -.30, a2 = .15, a3 = .38, a4 = -.12,
               a5 = .27, a6 = -.20, a7 = .42)
BW <- 1.4 + .35 * (sex == "M") - .28 * x + animal_bw[id] +
  .32 * sin(row * 1.7) + .18 * cos(row * .43)
data <- data.frame(id, BW = as.double(BW), sex, x)
data <- data[c(9:35, 1:8), ]
rownames(data) <- NULL

# Replace the two tables above with your own data; IDs link them explicitly.
# NA means unknown parent. Known parents must have a row in pedigree.
relationship <- gs_pedigree(pedigree, id = "id", sire = "sire", dam = "dam")
formulas <- BW ~ sex + x + (1 | id)
components <- list(
  animal = gvc("id", "BW", kernel = relationship, start = 0.5),
  residual = gvc("residual", "BW", start = 0.3))

# Existing formula-fixture controls, unchanged; no tuning to this example.
control <- gs_fit_control(
  direct_backend = "cholmod", trace_route = "selected_inverse",
  maximum_iterations = 80L, score_tolerance = 5e-6)
fit <- gfit(formulas, data = data, vc = components, task = "reml",
              control = control)

print(summary(fit))                  # named covariance estimates and convergence
print(coef(fit))                      # fixed intercept, sex contrast and x slope
animals <- fit$random_effects$animal  # public ID-labelled animal BLUP table
animals$observed <- animals$id %in% data$id
print(animals)                        # includes unphenotyped a8
print(list(converged = fit$converged, reason = fit$convergence_reason))
