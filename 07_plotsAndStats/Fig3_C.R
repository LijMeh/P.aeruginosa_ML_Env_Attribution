# ============================================================================
# Plot latent space PCA for Env XGBoost model
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
flat_df_annotated <- readRDS(
  data_here("00_data", "Figure_Processed_Data", "flat_df_annotated_SHAP_data.RDS")
)

# ---- Build Environment Network from Shared Gene Features ----

# Create all pairwise environment combinations
envs <- unique(flat_df_annotated$Env)
pair_indices <- combn(envs, 2, simplify = FALSE)

# Calculate edge weights based on shared genes between environment pairs
edge_stats <- lapply(pair_indices, function(pair) {
  # Get data for each environment
  df1 <- flat_df_annotated %>% filter(Env == pair[1])
  df2 <- flat_df_annotated %>% filter(Env == pair[2])
  
  # Merge on gene to find shared genes
  merged <- merge(
    df1[, c("gene", "Consensus_Pres_Summary", "mean_shap_pres", "mean_shap_abs")],
    df2[, c("gene", "Consensus_Pres_Summary", "mean_shap_pres", "mean_shap_abs")],
    by = "gene", suffixes = c("_1", "_2")
  )
  
  # Calculate presence expected weight (average SHAP for shared genes)
  pres_mask <- merged$Consensus_Pres_Summary_1 == "Presence_Expected" & 
               merged$Consensus_Pres_Summary_2 == "Presence_Expected"
  pres_weight <- sum(
    rowMeans(cbind(
      abs(merged$mean_shap_pres_1[pres_mask]),
      abs(merged$mean_shap_pres_2[pres_mask])
    ), na.rm = TRUE),
    na.rm = TRUE
  )
  
  # Calculate absence expected weight (average SHAP for shared genes)
  abs_mask <- merged$Consensus_Pres_Summary_1 == "Absence_Expected" & 
              merged$Consensus_Pres_Summary_2 == "Absence_Expected"
  abs_weight <- sum(
    rowMeans(cbind(
      abs(merged$mean_shap_abs_1[abs_mask]),
      abs(merged$mean_shap_abs_2[abs_mask])
    ), na.rm = TRUE),
    na.rm = TRUE
  )
  
  data.frame(
    from = pair[1],
    to = pair[2],
    pres_weight = pres_weight,
    abs_weight = abs_weight
  )
})

# Calculate edge weights based on NPV betweeen pairs
edge_stats <- lapply(pair_indices, function(pair) {
  # Get data for each environment
  df1 <- flat_df_annotated %>% filter(Env == pair[1])
  df2 <- flat_df_annotated %>% filter(Env == pair[2])
  
  # Merge on gene to find shared genes
  merged <- merge(
    df1[, c("gene", "Consensus_Pres_Summary", "Net_Predictive_Value")],
    df2[, c("gene", "Consensus_Pres_Summary", "Net_Predictive_Value")],
    by = "gene", suffixes = c("_1", "_2")
  )
  
  # Calculate NPV weight (average SHAP for shared genes with opposite presence/absence)
  npv_mask <- (merged$Consensus_Pres_Summary_1 == "Presence_Expected" & 
               merged$Consensus_Pres_Summary_2 == "Absence_Expected") |
              (merged$Consensus_Pres_Summary_1 == "Absence_Expected" & 
               merged$Consensus_Pres_Summary_2 == "Presence_Expected")
  npv_weight <- sum(
    rowMeans(cbind(
      abs(merged$Net_Predictive_Value_1[npv_mask]),
      abs(merged$Net_Predictive_Value_2[npv_mask])
    ), na.rm = TRUE),
    na.rm = TRUE
  )
  
  data.frame(
    from = pair[1],
    to = pair[2],
    weight = npv_weight
  )
})


# Combine edge statistics and calculate total weight
edge_df <- do.call(rbind, edge_stats) %>%
  mutate(weight = pres_weight + abs_weight) %>%
  filter(weight > 0)

# Combine edge statistics and calculate total weight (for NPV-based edges)
edge_df <- do.call(rbind, edge_stats) %>%
  mutate(weight = weight) %>%
  filter(weight > 0)

message("Environment pairs with shared genes: ", nrow(edge_df))
message("Edge weight range: ", round(min(edge_df$weight), 2), " to ", round(max(edge_df$weight), 2))

# Create network graph from edges
library(igraph)
g_tidy <- as_tbl_graph(edge_df[, c("from", "to", "weight")], directed = FALSE)

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

set.seed(44)

# ---- Create Network Plot ----
p.network <- g_tidy %>%
  ggraph(layout = "stress", weights = weight) +
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.7) +
  scale_edge_width_continuous(range = c(1, 6), guide = "none") +
  scale_edge_color_gradient(
    name = "Shared SHAP",
    low = "white",
    high = "#883268"
  ) +
  geom_node_point(aes(color = name), size = 15) +
  scale_color_manual(values = niche_colors, name = NULL, guide = "none") +
  geom_node_text(aes(label = str_replace_all(str_replace(name, "Bronchiectasis", "Bronch-\nectasis"), " ", "\n")), 
                 color = NA, size = 2.2) +
  geom_node_text(aes(label = str_replace_all(str_replace(name, "Bronchiectasis", "Bronch-\nectasis"), " ", "\n")),
                 data = function(x) filter(x, name %in% c("Adult CF", "Pediatric CF")),
                 color = "black", size = 2.2) +
  geom_node_text(aes(label = str_replace_all(str_replace(name, "Bronchiectasis", "Bronch-\nectasis"), " ", "\n")),
                 data = function(x) filter(x, !name %in% c("Adult CF", "Pediatric CF")),
                 color = "white", size = 2.2) +
  theme_graph(base_family = "sans") +
  theme(legend.position = "right") +
  scale_x_continuous(expand = expansion(mult = c(.15, .15))) +
  scale_y_continuous(expand = expansion(mult = c(.15, .15)))

# ---- Save Plot ----
ggsave(
  filename = data_here("00_data", "Figures", "Fig3_Env_Network.png"),
  plot = p.network,
  width = 6,
  height = 5,
  dpi = 300
)

