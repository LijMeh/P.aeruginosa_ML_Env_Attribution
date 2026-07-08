# ============================================================================
# Plot confusion matrix for non-blocked Env XGBoost model
# ============================================================================

# ---- Load Required Packages ----
library(install.load)
install_load("here")

# ---- Load configuration and setup ----
source(here("config.R"))
source(here("Figures_4th_Draft", "00_Setup.R"))
source(here("Figures_4th_Draft", "00_color_palette.R"))
source(here("Figures_4th_Draft", "00_functions.R"))

# ---- Load Processed Data ----
final_model_performance_data_processed <- readRDS(data_here("00_data", "Figure_Processed_Data", "Final_Model_Performance_Data_Processed.RDS"))

metadata_filtered <- readRDS(
  file = data_here("00_data", "Figure_Processed_Data", "Metadata_Filtered.rds")
)

# ---- Define environment order ----
env_labels <- c(
   "Aquatic", "Hospital",
  "Adult CF", "Pediatric CF", "Bronchiectasis", "Blood", "Urinary", "Fecal", "Wound"
)

cm_final <- final_model_performance_data_processed$confusion_matrix
accuracy_df_final <- final_model_performance_data_processed$accuracy_metrics
stats_final <- final_model_performance_data_processed$stats

# ---- Plot Confusion Matrix Heatmap ----
cm_plot_data_final <- cm_final$confusion_row_pct %>%
  filter(True %in% env_labels, Pred %in% env_labels) %>%
  mutate(
    True = factor(True, levels = env_labels),
    Pred = factor(Pred, levels = env_labels)
  )

p.confusion_matrix <- ggplot(cm_plot_data_final, aes(x = Pred, y = True, fill = row_pct)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = n), size = 3.5, color = "black") +
  scale_fill_gradient(
    low = "white", 
    high = "#883268",
    limits = c(0, 1), 
    breaks = seq(0, 0.8, 0.1)
  ) +
  labs(
    x = "Predicted labels",
    y = "True labels",
    fill = "Proportion of\nRow Total"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10, margin = margin(t = 20)),
    axis.title.y = element_text(size = 10, margin = margin(r = 20)),
    legend.position = "none",
    panel.grid = element_blank()
  ) 

p.confusion_matrix

# ---- Plot Accuracy Comparison ----
accuracy_colors <- c(
  actual = "#cfa3c2", 
  null   = "#883268"
)

p.accuracy <- ggplot(accuracy_df_final, aes(x = category, y = accuracy, fill = type, alpha = type)) +
  geom_col(position = "identity", width = 0.7) +
  geom_text(
    data = accuracy_df_final %>% mutate(
      y_pos = if_else(category == "Wound" & type == "actual", accuracy + 10, accuracy)
    ),
    aes(label = sprintf("%.1f%%", accuracy), y = y_pos), 
    hjust = -0.1, size = 4
  ) +
  scale_fill_manual(
    values = accuracy_colors,
    labels = c(actual = "Actual Accuracy", null = "Null Accuracy"),
    name = "Accuracy Type"
  ) +
  scale_alpha_manual(values = c(actual = 0.5, null = 1), guide = "none") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4))) + # increased right space
  labs(x = NULL, y = "Accuracy (%)") +
  theme_void(base_size = 18) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) 

p.accuracy

# ---- Combine Confusion Matrix and Accuracy Plots ----
p.conf_mat_combined <- (p.confusion_matrix | p.accuracy) + 
  plot_layout(widths = c(2, 1))

p.conf_mat_combined

# ---- Save Plot ----
ggsave(
  filename = data_here("00_data", "Figures", "Fig2_ConfusionMatrix_Final.png"),
  plot = p.conf_mat_combined,
  width = 9.7,
  height = 5,
  dpi = 300
)

# ---- Create Stats Table Figure ----
stats_table <- data.frame(
  Metric = c("Weighted F1", "Macro F1", "Accuracy", "Balanced Accuracy"),
  Value = c(
    sprintf("%.3f", stats_final$Micro.F1[1]),
    sprintf("%.3f", stats_final$Macro.F1[1]),
    sprintf("%.3f", stats_final$Accuracy[1]),
    sprintf("%.3f", stats_final$Balanced.Accuracy[1])
  )
)

p.stats_table <- ggplot(stats_table, aes(x = 1, y = 1)) +
  geom_blank() +
  annotate(
    "text",
    x = 0.5, y = 0.7,
    label = "Model Performance Metrics",
    size = 6, fontface = "bold",
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = 0.75, y = 0.6,
    label = paste(stats_table$Metric, stats_table$Value, sep = ": ", collapse = "\n"),
    size = 5,
    hjust = 1,
    vjust = 1
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void()

p.stats_table


# ---- Save Stats Table ----
ggsave(
  filename = data_here("00_data", "Figures", "Fig2_StatsTable.png"),
  plot = p.stats_table,
  width = 4,
  height = 3,
  dpi = 300
)








