# ---- Variance Partitioning: Env vs Phylogeny ----
# Uses the original phenotype matrix (not PCs) as multivariate response,
# and uses phylogenetic eigenvectors from cophenetic distances.

# ---- Load Required Packages ----
library(install.load)
install_load("here")

# ---- Load configuration and setup ----
source(here("config.R"))
source(here("Figures_5th_Draft", "00_Setup.R"))
source(here("Figures_5th_Draft", "00_color_palette.R"))
source(here("Figures_5th_Draft", "00_functions.R"))

# ---- Load Processed Data ----
filtered_experimental_metadata <- readRDS(
  data_here("00_data", "Figure_Processed_Data", "filtered_experimental_metadata.RDS")
)

pca_res <- readRDS(
  data_here("00_data", "Figure_Processed_Data", "pca_res_experimental_phenotypes.RDS")
)

filtered_experimental_tree <- readRDS(
  data_here("00_data", "Figure_Processed_Data", "filtered_experimental_tree.RDS")
)

cols <- c("Elastase Production", "Biofilm Fraction", "Growth Rate")

# 1) Build analysis dataframe with complete cases
keep <- complete.cases(
  filtered_experimental_metadata[, cols],
  filtered_experimental_metadata$Env,
  filtered_experimental_metadata$treeNames
)

df <- filtered_experimental_metadata[keep, , drop = FALSE]
df$Env <- recode(df$Env, "Non-CF Human Infection" = "Acute Human Infection")

# 2) Response matrix (phenotypes)
Y <- as.matrix(df[, cols, drop = FALSE])

# Drop zero-variance columns before scaling to avoid NaN values
sdv <- apply(Y, 2, sd, na.rm = TRUE)
Y <- Y[, sdv > 0, drop = FALSE]
Y <- scale(Y)
rownames(Y) <- df$treeNames

# 3) Align to phylogenetic distance matrix
phylo_dist <- cophenetic(filtered_experimental_tree)
ids <- intersect(df$treeNames, rownames(phylo_dist))

df <- df[base::match(ids, df$treeNames), , drop = FALSE]
Y <- Y[base::match(ids, rownames(Y)), , drop = FALSE]
D <- phylo_dist[ids, ids]

# 4) Phylogenetic eigenvectors from distance matrix
pcoa <- cmdscale(as.dist(D), eig = TRUE, add = TRUE, k = nrow(D) - 1)
pos <- which(pcoa$eig > 0)

phy_axes_all <- as.data.frame(pcoa$points[, pos, drop = FALSE])
phy_prop <- pcoa$eig[pos] / sum(pcoa$eig[pos])

# Keep enough axes to explain 80% positive-eigenvalue variance (cap at 10)
k_phy <- min(10, which(cumsum(phy_prop) >= 0.80)[1])
phy_axes <- phy_axes_all[, seq_len(k_phy), drop = FALSE]
colnames(phy_axes) <- paste0("Phy", seq_len(k_phy))

# 5) Environment design matrix
env_mm <- model.matrix(~ Env, data = df)[, -1, drop = FALSE]

# 6) Variance partitioning
vp <- varpart(Y, env_mm, phy_axes)
print(vp)
plot(vp, Xnames = c("Environment", "Phylogeny"))

# 7) Significance tests for unique fractions
# Env unique effect controlling for phylogeny
rda_env_unique <- rda(Y, env_mm, phy_axes)
env_unique_test <- anova.cca(rda_env_unique, permutations = 999)

# Phylogeny unique effect controlling for env
rda_phy_unique <- rda(Y, phy_axes, env_mm)
phy_unique_test <- anova.cca(rda_phy_unique, permutations = 999)

print(env_unique_test)
print(phy_unique_test)

# 8) Optional: adonis2 with same predictor logic (marginal tests)
adon_df <- cbind(data.frame(Env = factor(df$Env)), phy_axes)
adon_formula <- as.formula(
  paste("Y ~ Env +", paste(colnames(phy_axes), collapse = " + "))
)

adon_both <- adonis2(adon_formula, data = adon_df, method = "euclidean", by = "margin")
print(adon_both)

# 9) PERMANOVA assumption check for Env dispersion
dY <- dist(Y, method = "euclidean")
bd_env <- betadisper(dY, group = factor(df$Env))
print(permutest(bd_env, permutations = 999))

# 10) Extended PERMDISP Analysis: Testing Dispersion for Both Marginal and Partial Effects
# PERMDISP (betadisper + permutest) tests whether average distance to group centroid
# differs among groups. A non-significant PERMDISP result indicates that significant
# RDA/PERMANOVA results reflect true centroid differences (phenotype shifts) rather than
# dispersion artifacts (differential within-group variability).

cat("\n=== PERMDISP: Heterogeneity of Dispersion ===\n")
cat("Tests whether phenotype variation within environment groups differs.\n")
cat("Non-significant results support interpretation of true phenotype shifts.\n\n")

# A) Overall Env effect (marginal): Does within-group phenotype variation differ by Env?
cat("A) Overall Env effect - within-group dispersion:\n")
bd_env_overall <- betadisper(dY, group = factor(df$Env))
test_env_overall <- permutest(bd_env_overall, permutations = 999)
print(test_env_overall)
cat("Interpretation: If non-significant (p > 0.05), Env groups have similar within-group\n")
cat("variability, strengthening confidence that Env RDA significance reflects true shifts.\n\n")

# B) Phylogeny effect: Does within-group phenotype variation differ by Phylogeny?
cat("B) Phylogeny effect - within-group dispersion:\n")
phylo_groups <- df$cluster  # Using clusters derived from phylogenetic distance
if (!is.null(phylo_groups) && !any(is.na(phylo_groups))) {
  bd_phy <- betadisper(dY, group = factor(phylo_groups))
  test_phy <- permutest(bd_phy, permutations = 999)
  print(test_phy)
  cat("Interpretation: If non-significant, phylogenetic groups show similar dispersion.\n\n")
} else {
  cat("Phylogenetic groups unavailable or contain NAs; skipping phylogeny dispersion test.\n\n")
}

# C) Env effect AFTER removing phylogeny: Test dispersion on residuals
# This is the parallel PERMDISP test for the partial RDA (rda_env_unique)
cat("C) Env effect AFTER controlling for phylogeny - residual dispersion:\n")
cat("Testing whether environment groups differ in within-group dispersion of\n")
cat("phenotypes after phylogenetic effects are removed.\n")

# Get residuals from a model with only phylogeny
rda_phy_only <- rda(Y, phy_axes)
Y_residuals <- residuals(rda_phy_only)

# Calculate distances on phylogeny-adjusted residuals
dY_residuals <- dist(Y_residuals, method = "euclidean")
bd_env_partial <- betadisper(dY_residuals, group = factor(df$Env))
test_env_partial <- permutest(bd_env_partial, permutations = 999)
print(test_env_partial)
cat("Interpretation: If non-significant, the partial RDA Env effect (rda_env_unique)\n")
cat("reflects true phenotype shifts among Env groups, not dispersion differences.\n\n")

# Summary interpretation
cat("=== Summary ===\n")
cat("If both tests A and C are non-significant, this strongly supports that the\n")
cat("detected environment effect represents real phenotype differentiation rather\n")
cat("than differential within-group variability (which can inflate PERMANOVA p-values).\n\n")