# ============================================================================
# ---- Experimental Data PCA ----
# ============================================================================

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

# Create PCA data frame with grouping variables
cols <- c("Elastase Production", "Biofilm Fraction", "Growth Rate")

na_mask <- !is.na(rowSums(filtered_experimental_metadata[, cols]))
pca_df <- data.frame(
  pca_res$x,
  Cluster = factor(filtered_experimental_metadata$cluster[na_mask]),
  Env = filtered_experimental_metadata$Env[na_mask]
) %>%
  mutate(Env = recode(Env, "Non-CF Human Infection" = "Acute Human Infection"))

# ---- Prepare PCA Loadings ----

loadings_df <- data.frame(
  Feature = c("Elastase Production", "Biofilm\nFraction", "Growth Rate"),
  PC1 = pca_res$rotation[, 1],
  PC2 = pca_res$rotation[, 2]
)

arrow_scale <- 3

env_colors <- c(
  "Cystic Fibrosis"                     = "#F3BA36",
  "Acute Human Infection"  = "#B6502B",
  "Environmental"            = "#9CA780"
)

var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
pc1_lab <- paste0("PC1 (", sprintf("%.1f", 100 * var_explained[1]), "%)")
pc2_lab <- paste0("PC2 (", sprintf("%.1f", 100 * var_explained[2]), "%)")

# ---- Plot PCA with Loadings ----

p.PCA_Env_loadings <- ggplot(pca_df, aes(PC1, PC2, color = Env)) +
  geom_point(size = 3) +
  scale_color_manual(values = env_colors) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-3, 4) +
  ylim(-3, 4) +
  labs(x = pc1_lab, y = pc2_lab) +
  geom_segment(
    data = loadings_df,
    aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
    arrow = arrow(length = unit(0.2, "cm")), 
    color = "black"
  ) +
    geom_text(
    data = loadings_df,
    aes(x = PC1 * arrow_scale * 1.1, y = PC2 * arrow_scale * 1.1, label = Feature),
    color = "black", 
    size = 4
  )

# ---- Save PCA Plot ----
ggsave(
  filename = data_here("00_data", "Figures", "Fig5_ExperimentalPCA_Env_Loadings_v3.png"),
  plot = p.PCA_Env_loadings,
  width = 5,
  height = 5,
  dpi = 300
)

# -----------------------------
# 0) Build clean feature table
# -----------------------------
cols <- c("ecr.n", "cont.max.rate", "bfrac")

feature_df <- filtered_experimental_metadata %>%
  mutate(Env = recode(Env, "Non-CF Human Infection" = "Acute Human Infection")) %>%
  select(Env, all_of(cols)) %>%
  tidyr::drop_na()
  
# z-score features so effects are comparable
feature_z <- feature_df %>%
  mutate(across(all_of(cols), ~ as.numeric(scale(.x))))

# ------------------------------------------------------
# 1) Direct interpretation: env means per feature (z)
# ------------------------------------------------------
env_feature_summary <- feature_z %>%
  pivot_longer(cols = all_of(cols), names_to = "Feature", values_to = "z_value") %>%
  group_by(Env, Feature) %>%
  summarise(
    mean_z = mean(z_value),
    sd_z = sd(z_value),
    n = n(),
    se_z = sd_z / sqrt(n),
    .groups = "drop"
  ) %>%
  arrange(Feature, desc(mean_z))

env_feature_summary
# Interpretation:
# mean_z > 0 = that env tends higher-than-overall for that feature
# mean_z < 0 = that env tends lower-than-overall

# -------------------------------------------------------------------
# 3) PCA-based env-feature association using loadings + env centroids
# -------------------------------------------------------------------
# Uses the first 2 PCs that your plot is based on.

# loadings from PCA
loadings_long <- as_tibble(pca_res$rotation[, 1:2], rownames = "Feature") %>%
  rename(load_PC1 = PC1, load_PC2 = PC2)

# env centroids in PCA space
pc_centroids <- pca_df %>%
  group_by(Env) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

# association score: projection of env centroid onto each feature loading
pca_env_feature_assoc <- tidyr::crossing(pc_centroids, loadings_long) %>%
  mutate(
    assoc_score_PC1_PC2 = PC1 * load_PC1 + PC2 * load_PC2
  ) %>%
  arrange(Env, desc(assoc_score_PC1_PC2))

pca_env_feature_assoc
# Higher positive score => feature aligns more strongly with that env in PCA space.
