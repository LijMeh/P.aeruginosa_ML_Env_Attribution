# ============================================================================
# Plot ANI depth vs. ENV diversity
# ============================================================================

# ---- Load Required Packages ----
library(install.load)
install_load("here")

# ---- Load configuration and setup ----
source(here("config.R"))
source(here("Iterative_Model_Approach", "00_setup.R"))
source(here("Iterative_Model_Approach", "00_color_palette.R"))
source(here("Iterative_Model_Approach", "00_functions.R"))

# ---- Load Processed Data ----
metadata_filtered <- readRDS(data_here("00_data", "Figure_Processed_Data", "Metadata_Filtered.rds"))

ani_data <- data_here("00_Universal", "ANI_Clusters/synthetic_group") %>%
  list.files(pattern = "\\.csv$", full.names = TRUE) %>%
  purrr::map(~{
    df <- read_csv(.x)
    cluster <- str_extract(basename(.x), "\\d+\\.\\d+")
    colnames(df)[2] <- cluster
    df
  }) %>%
  reduce(full_join, by = colnames(.[[1]])[1])

# ---- Merge with Metadata ----
metadata_filtered_prepped <- metadata_filtered %>%
  mutate(Genome_ID = str_replace_all(Genome_ID, "\\.", "_"))

key_col <- colnames(ani_data)[1]
merged_data <- left_join(
  metadata_filtered_prepped, 
  ani_data, 
  by = c("Genome_ID" = key_col)
)

# ---- Reshape to Long Format ----
ani_long <- merged_data %>%
  pivot_longer(
    cols = matches("^\\d+\\.\\d+$"),
    names_to = "ANI_Cluster",
    values_to = "ANI_Value"
  ) %>%
  filter(!is.na(ANI_Value)) %>%
  add_count(ANI_Cluster, ANI_Value, name = "value_count") %>%
  filter(value_count > 8) %>%  # Min 9 samples per group
  select(-value_count)

# ---- Calculate Niche Diversity per ANI Group ----
ani_summary <- ani_long %>%
  filter(!is.na(ANI_Value), !is.na(Niche)) %>%
  group_by(ANI_Cluster, ANI_Value) %>%
  summarise(
    num_unique_niches = pmin(n_distinct(Niche), 9),
    group_size = n(),
    .groups = "drop"
  )

# ---- Filter to Selected ANI Thresholds ----
ani_keep <- c(99.33, 99.5, 99.75, 99.96, 99.99)
ani_summary_filtered <- ani_summary %>% 
  filter(ANI_Cluster %in% ani_keep)

# ---- Calculate Weighted Mean per ANI Threshold ----
ani_mean_points <- ani_summary_filtered %>%
  group_by(ANI_Cluster) %>%
  summarise(
    mean_niches = weighted.mean(num_unique_niches, group_size),
    .groups = "drop"
  )

# ---- Plot ANI vs. ENV Diversity ----
p.ani_env <- ggplot(
  ani_summary_filtered, 
  aes(x = ANI_Cluster, y = num_unique_niches)
) +
  geom_jitter(aes(color = group_size), width = 0.2, alpha = 0.8) +
  scale_color_gradientn(
    colors = c("#F3E8F0", "#CFA3C2", "#9A4F7F", "#914174", "#883268"),
    trans = "log10",
    limits = c(10, 1000),
    oob = scales::squish,
    breaks = c(10, 25, 50, 100, 500),
    labels = c("10", "25", "50", "100", "500"),
    name = "Samples in group"
  ) +
  guides(color = guide_colorbar(barwidth = 0.8, barheight = 4)) +
  geom_point(
    data = ani_mean_points,
    aes(x = ANI_Cluster, y = mean_niches),
    color = "black", 
    shape = 8, 
    size = 2, 
    stroke = 1,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(breaks = 1:9, limits = c(1, 9)) +
  labs(x = "ANI Cluster", y = "Number of Environments") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 8),
    axis.title.x = element_text(margin = margin(t = 15))
  )

ggsave(
  filename = data_here("00_data", "Figures", "Fig1_ANIdepth.png"),
  plot = p.ani_env,
  width = 4,
  height = 3,
  dpi = 300
)

