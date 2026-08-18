# ============================================================================
# Plot importance and mean shap across genes
# ============================================================================

# ---- Load Required Packages ----
library(install.load)
install_load("here")

# ---- Load configuration and setup ----
source(here("config.R"))
source(here("00_R_configs", "00_setup.R"))
source(here("00_R_configs", "00_color_palette.R"))
source(here("00_R_configs", "00_functions.R"))

# ---- Load Processed Data ----
flat_df_annotated_gene_names <- readRDS(
  file = data_here("00_data", "Figure_Processed_Data", "flat_df_annotated_gene_names_SHAP_data.RDS")
)

gene_info <- readRDS(
  file = data_here("00_data", "Figure_Processed_Data", "gene_variant_PA_Annotations.rds")
)

gene_info %>% filter(Gene.Name %in% "tonB2")

# Define exact groupings
chronic_envs <- c("Adult CF", "Pediatric CF", "Bronchiectasis")
acute_envs <- c("Fecal", "Blood", "Urinary", "Wound")

# Define niche colors
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

# Calculate GLOBAL mean_abs_shap per gene (across all 9 envs)
global_importance <- flat_df_annotated_gene_names %>%
  group_by(gene) %>%
  summarise(global_mean_abs_shap = mean(mean_abs_shap, na.rm = TRUE), .groups = "drop")

# Define genes to include in plot (ALL)
genes_above_global_threshold <- global_importance$gene

# Merge importance with gene info and finalize gene names
global_importance_named <- global_importance %>%
  left_join(gene_info, by = c("gene" = "Prot_acc")) %>%
  mutate(
    gene_display_name = ifelse(
      !is.na(Gene.Name) & Gene.Name != "",
      Gene.Name,
      gene
    )
  )

# ---- Create data for lines connecting min/max NPV for each gene ----
gene_lines <- flat_df_annotated_gene_names %>%
  left_join(global_importance_named, by = "gene") %>%
  filter(gene %in% genes_above_global_threshold) %>%
  group_by(gene) %>%
  arrange(Net_Predictive_Value) %>%
  slice(c(1, n())) %>%
  ungroup()

# ---- Create Plot ----
# First, prepare all data with transformations
plot_data <- flat_df_annotated_gene_names %>%
  left_join(global_importance_named, by = "gene") %>%
  group_by(gene) %>%
  mutate(
    max_npv = max(Net_Predictive_Value, na.rm = TRUE),
    min_npv = min(Net_Predictive_Value, na.rm = TRUE),
    is_max = Net_Predictive_Value == max_npv,
    is_min = Net_Predictive_Value == min_npv,
    meets_threshold = gene %in% genes_above_global_threshold
  ) %>%
  ungroup()

# Create alternating positions for labeled genes based on their x-axis order
labeled_genes_temp <- plot_data %>%
  distinct(gene, global_mean_abs_shap) %>%
  arrange(desc(global_mean_abs_shap)) %>%
  mutate(
    label_position = rep(c("top", "bottom"), length.out = n())
  ) %>%
  select(gene, label_position)

plot_data <- plot_data %>%
  left_join(labeled_genes_temp, by = "gene") %>%
  group_by(gene) %>%
  mutate(
    label = case_when(
      meets_threshold & label_position == "top" & is_max ~ gene_display_name,
      meets_threshold & label_position == "bottom" & is_min ~ gene_display_name,
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

# Extract gene lines from plot data
gene_lines_data <- plot_data %>%
  group_by(gene) %>%
  arrange(Net_Predictive_Value) %>%
  slice(c(1, n())) %>%
  ungroup()

# Create plot (version 1: with alternating labels)
p_shap_variability_all <- plot_data %>%
  ggplot(aes(x = global_mean_abs_shap, y = Net_Predictive_Value, 
             color = Env, size = global_mean_abs_shap)) +
  geom_line(data = gene_lines_data, aes(group = gene, x = global_mean_abs_shap, y = Net_Predictive_Value), 
            color = "gray50", alpha = 0.3, linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(alpha = 0.65, stroke = 0) +
  ggrepel::geom_text_repel(aes(label = label), 
            size = 3, fontface = "bold", color = "black", 
            min.segment.length = 0, box.padding = 1, max.overlaps = Inf,
            force = 3, force_pull = 0) +
  scale_color_manual(values = niche_colors) +
  scale_size_continuous(range = c(0.5, 5), guide = "none") +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = scales::breaks_pretty(),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  labs(
    x = "Global Mean Absolute SHAP",
    y = "Net Predictive Value (per ENV)",
    color = "Environment",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14)
  )

# Create version 2: without alternating labels (all on max points)
plot_data_TopOnly <- plot_data %>%
  mutate(
    label = ifelse(meets_threshold & is_max, gene_display_name, NA_character_),
    label_fontface = case_when(
      Variant_Type == "WT" ~ "plain",
      Variant_Type == "Mut" ~ "bold",
      TRUE ~ "plain"
    )
  )

wave_df <- make_hist_wave_npv_units(
  x = plot_data$global_mean_abs_shap,
  y_npv = plot_data$Net_Predictive_Value,
  bins = 80,
  smooth_window = 6
)

p_shap_variability_all <- plot_data_TopOnly %>%
  ggplot(aes(
    x = global_mean_abs_shap,
    y = Net_Predictive_Value,
    color = Env,
    size = global_mean_abs_shap
  )) +
  geom_area(
    data = wave_df %>% dplyr::arrange(type, x),
    aes(x = x, y = y, fill = type, group = type),
    inherit.aes = FALSE,
    alpha = 0.20,
    color = NA
  ) +
  geom_hline(yintercept = 0, color = "gray70", linewidth = 0.3) +
  geom_line(
    data = gene_lines_data,
    aes(group = gene, x = global_mean_abs_shap, y = Net_Predictive_Value),
    color = "gray5", alpha = 0.7, linewidth = 0.7, inherit.aes = FALSE
  ) +
  geom_point(alpha = 0.8, stroke = 0) +
  ggrepel::geom_text_repel(
    data = ~ .x %>% dplyr::filter(global_mean_abs_shap > 0.01),
    aes(label = label, fontface = label_fontface),
    family = "sans",
    size = 2.5,
    color = "black",
    min.segment.length = 0,
    box.padding = 0.5,
    max.overlaps = Inf,
    force = 2,
    force_pull = 0
  ) +
  scale_fill_manual(
    values = c(npv_pos = "gray40", npv_neg = "gray40"),
    guide = "none"
  ) +
  scale_color_manual(values = niche_colors) +
  scale_size_continuous(range = c(0.5, 5), guide = "none") +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 20, sigma = 0.001),
    breaks = scales::breaks_pretty(),
    labels = scales::label_number(accuracy = 0.01)
  ) +
    scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 20, sigma = 0.05),
    breaks = c(-2, -1, - 0.5, -0.25, -0.1, 0,  0.1, 0.25, 0.5, 1, 2),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  labs(
    x = "Global Mean Absolute SHAP",
    y = "Net Predictive Value (per ENV)",
    color = "Environment"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none"
  )

# ---- Save Plot ----
ggsave(
  filename = data_here("00_data", "Figures", "Fig4_SHAP_Variability_AllEnvs.png"),
  plot = p_shap_variability_all,
  width = 8,
  height = 5,
  dpi = 300
)







