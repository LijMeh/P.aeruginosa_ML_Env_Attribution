# ============================================================================
# Plot ROC curve for non-blocked Env XGBoost model
# ============================================================================

# ---- Load Required Packages ----
library(install.load)
install_load("here")

# ---- Load configuration and setup ----
source(here("config.R"))
source(here("00_R_configs", "00_Setup.R"))
source(here("00_R_configs", "00_color_palette.R"))
source(here("00_R_configs", "00_functions.R"))

# ---- Load Processed Data ----
roc_data <- readRDS(data_here("00_data", "Figure_Processed_Data", "477_roc_data_processed.rds"))

probs_long <- readRDS(data_here("00_data", "Figure_Processed_Data", "477_probs_processed.rds"))

metadata_unmapped <- read.csv(data_here("00_data", "blockedModelInputMetadata.csv")) %>%
  mutate(Genome_ID = str_replace_all(Genome_ID, "\\.", "_"))

niche_map <- c(
  "Bronchiectasis" = "Bronchiectasis",
  "CF" = "Adult CF",
  "early.CF" = "Pediatric CF",
  "Early CF" = "Pediatric CF",
  "Early.CF" = "Pediatric CF",
  "Urinary" = "Urinary",
  "Gastrointestinal" = "Gastrointestinal",
  "Wound" = "Wound",
  "Blood" = "Blood",
  "Aquatic" = "Aquatic",
  "Ocean" = "Aquatic",
  "Hospital Environment" = "Hospital",
  "Hospital.Environment" = "Hospital",
  "Hospital.Enviornment" = "Hospital",
  "Rectal/Feces" = "Fecal",
  "Rectal-Feces" = "Fecal",
  "Rectal.Feces" = "Fecal"
)

metadata_mapped <- metadata_unmapped %>%
  mutate(Niche = remap_niches(Niche)) %>%
  filter(!is.na(Niche))

# ---- Define Color Palette ----
niche_colors <- c(
  "Hospital"  = "#727272",
  "Aquatic"               = "#9CBED2",
  "Adult CF"                    = "#F3BA36",
  "Pediatric CF"              = "#F7D074",
  "Bronchiectasis"        = "#E36B2E",
  "Blood"                 = "#7F4348",
  "Urinary"               = "#A8441B",
  "Fecal"          = "#B89C8D",
  "Wound"                 = "#C57066"
)

# ---- Calculate AUC for each environment ----
env_order <- c("Aquatic", "Hospital", "Adult CF", "Pediatric CF", "Bronchiectasis", "Blood", "Urinary", "Fecal", "Wound")

auc_data <- roc_data %>%
  group_by(Env) %>%
  arrange(FPR) %>%
  summarise(
    AUC = sum(diff(FPR) * (head(TPR, -1) + tail(TPR, -1)) / 2),
    .groups = "drop"
  ) %>%
  mutate(Env = factor(Env, levels = env_order))

# ---- Plot ROC Curve ----
roc_data$Env <- factor(roc_data$Env, levels = env_order)

p.roc_curve <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Env)) +
  geom_line(size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = niche_colors, name = "Environment") +
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  )             

# ---- Save the plot ----
ggsave(
  filename = data_here("00_data", "Figures", "Fig2_ROC_Curve.png"),
  plot = p.roc_curve,
  width = 7,
  height = 5,
  dpi = 300
)       

# ---- Small AUC bar plot ----
auc_data$Env <- factor(auc_data$Env, levels = env_order)

p.auc_bar <- ggplot(auc_data, aes(x = Env, y = AUC, fill = Env)) +
  geom_col() +
  scale_fill_manual(values = niche_colors) +
  labs(x = "Environment", y = "AUROC", title = "AUROC by Environment") +
  theme_bw(base_size = 10) +
  theme(
    plot.background = element_rect(fill = NA, color = NA)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# ---- Save the plot ----
ggsave(
  filename = data_here("00_data", "Figures", "Fig2_AUC_Bar.png"),
  plot = p.auc_bar,
  width = 2.5,
  height = 3,
  dpi = 300
)

# ---- Process Data for Next Best Threshold Plot ----

# Identify focal environment per sample from REAL metadata
probs_summary <- probs_long %>%
  left_join(
    metadata_mapped %>% 
      select(Genome_ID, Niche) %>% 
      rename(Genome_IDs = Genome_ID, Focal_Env = Niche),
    by = "Genome_IDs"
  )

# For each Focal_Env, calculate mean probability for ALL environments  
plot_data <- probs_summary %>%
  group_by(Focal_Env, Env) %>%
  summarise(
    Mean_Prob = mean(Probability),
    .groups = "drop"
  ) %>%
  # Reorder so Focal appears first within each group
  mutate(
    Env_Order = case_when(
      Env == Focal_Env ~ 0,  # Focal first
      TRUE ~ 1
    )
  ) %>%
  arrange(Focal_Env, Env_Order, Env)

env_order <- c("Aquatic", "Hospital", "Adult CF", "Pediatric CF", "Bronchiectasis", "Blood", "Urinary", "Fecal", "Wound")
plot_data$Env <- factor(plot_data$Env, levels = env_order)
plot_data$Focal_Env <- factor(plot_data$Focal_Env, levels = env_order)

p.best_threshold <- ggplot(plot_data, aes(x = Focal_Env, y = Mean_Prob, fill = Env)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = niche_colors) +
  labs(x = "Actual Environment", y = "Mean Probability", fill = "Evaluated Environment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  )


# ---- Save the plot ----
ggsave(
  filename = data_here("00_data", "Figures", "Fig2_Best_Threshold.png"),
  plot = p.best_threshold,
  width = 8,
  height = 4,
  dpi = 300
)

# ---- Histogram of Probability Distributions by Focal_Env ----
probs_summary$Env <- factor(probs_summary$Env, levels = env_order)
probs_summary$Focal_Env <- factor(probs_summary$Focal_Env, levels = env_order)

p.prob_dist <- ggplot(probs_summary, aes(x = Probability, fill = Env)) +
  geom_histogram(position = "dodge", bins = 10, alpha = 1) +
  facet_wrap(~Focal_Env, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = niche_colors) +
  labs(
    x = "Probability",
    y = "Frequency",
    fill = "Environment",
    title = "Probability Distributions by Focal Environment"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold")
  )

# ---- Save the plot ----
ggsave(
  filename = data_here("00_data", "Figures", "Fig2_Prob_Distributions.png"),
  plot = p.prob_dist,
  width = 12,
  height = 10,
  dpi = 300
)




