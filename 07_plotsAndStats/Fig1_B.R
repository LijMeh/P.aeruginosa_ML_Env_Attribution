# ============================================================================
# Plot Phylogenetic Tree with Niche Heatmap and ConsenTRAIT values
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
phylo_filtered_rooted_treedata <- readRDS(
  data_here("00_data", "Figure_Processed_Data", "phylo_filtered_rooted_treedata.RDS")
)

consenTrait_data <- readRDS(
  data_here("00_data", "Figure_Processed_Data", "consenTrait_data.RDS")
)

# ---- Define Color Palette ----
niche_colors <- c(
  "Hospital"  = "#727272",
  "Aquatic"               = "#9CBED2",
  "Adult CF"              = "#F3BA36",
  "Pediatric CF"          = "#F7D074",
  "Bronchiectasis"        = "#E36B2E",
  "Blood"                 = "#7F4348",
  "Urinary"               = "#A8441B",
  "Fecal"                 = "#B89C8D",
  "Wound"                 = "#C57066"
)

# ---- Define Niche Order ----
niche_order <- c(
  "Aquatic",
  "Hospital",
  "Adult CF", 
  "Pediatric CF", 
  "Bronchiectasis", 
  "Blood", 
  "Urinary", 
  "Fecal", 
  "Wound")

# ---- Log-transform Branch Lengths ----
phylo_filtered_rooted_treedata_log <- phylo_filtered_rooted_treedata
#phylo_filtered_rooted_treedata_log@phylo$edge.length <- log1p(
#  phylo_filtered_rooted_treedata_log@phylo$edge.length
#)

# NOT log transforming anymore

# ---- Remove Outgroup from Tree ----
outgroup_tips <- as_tibble(phylo_filtered_rooted_treedata_log) %>%
  filter(Niche == "Outgroup") %>%
  pull(label)

pruned_tree <- ape::drop.tip(phylo_filtered_rooted_treedata_log@phylo, outgroup_tips)

# ---- Filter to Tips Only (after removing outgroup) ----
tree_tips <- pruned_tree$tip.label

# ---- Get all niches excluding NA and Outgroup ----
all_niches <- as_tibble(phylo_filtered_rooted_treedata_log) %>%
  filter(!is.na(Niche), Niche != "Outgroup") %>%
  pull(Niche) %>%
  unique()

# ---- Create binary presence/absence matrix (excluding outgroup and NA) ----
tip_niche_combos <- as_tibble(phylo_filtered_rooted_treedata_log) %>%
  select(label, Niche) %>%
  filter(!is.na(Niche), Niche != "Outgroup")

heatmap_long <- expand.grid(
  label = tree_tips, 
  Niche = all_niches, 
  stringsAsFactors = FALSE
) %>%
  mutate(value = as.integer(paste(label, Niche) %in% 
                            paste(tip_niche_combos$label, tip_niche_combos$Niche)),
         Niche = factor(Niche, levels = niche_order)) %>%
  filter(value == 1)

# ---- Create Tree Plot ----
p.tree <- ggtree(pruned_tree, color = "darkgrey", size = .5) +
  geom_treescale(
    x = 0,
    y = -150,      # move lower (more negative = lower)
    width = 0.1,
    offset = 0.5,
    fontsize = 3
  ) +
  coord_cartesian(ylim = c(-25, length(tree_tips) + 0.5), clip = "off") +
  theme(plot.margin = margin(0, 5, 30, 5))

# Scale bar doesn't work with log transformed tree

# ---- Create Heatmap Plot ----
p.heatmap <- ggplot(
    heatmap_long,
    aes(
        x = Niche,
        y = factor(label, levels = rev(tree_tips)),
        fill = Niche,
        color = Niche
    )
) +
  geom_tile(width = 0.9, height = 1, size = 0.3) +
  scale_fill_manual(values = niche_colors, na.value = "white") +
  scale_color_manual(values = niche_colors) +
  scale_x_discrete(position = "top") +
  theme_void() +
  theme(
    legend.position = "none",
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0.1),
    axis.ticks.x.top = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0, 25, 35, 5)
  ) +
  labs(x = NULL, y = NULL) +
  # Add ConsenTRAIT Mean_depth with significance indicator below niches
  geom_text(
    data = consenTrait_data %>%
      filter(!is.na(Niche), Niche != "Outgroup") %>%
      mutate(
        Niche = factor(Niche, levels = niche_order),
        depth_label = paste0(format(Mean_depth, digits = 2), ifelse(P_value < 0.05, "*", "")),
        fontface = ifelse(P_value < 0.05, "bold", "plain")
      ),
    aes(
      x = Niche,
      y = -100,
      label = depth_label,
      fontface = fontface
    ),
    inherit.aes = FALSE,
    vjust = 1,
    hjust = .9,
    size = 3.5,
    angle = 45
  ) +
  coord_cartesian(clip = "off", ylim = c(0.5, length(tree_tips) + 0.5)) +
  theme(plot.margin = margin(0, 25, 35, 0)) +
  scale_y_discrete(limits = rev(tree_tips))

# ---- Combine Plots ----
p.30Ktree <- p.tree + p.heatmap + plot_layout(widths = c(2, 3))

# ---- Save Figure ----
ggsave(
  filename = data_here("00_data", "Figures", "Fig1_PhyloTree.png"),
  plot = p.30Ktree,
  width = 12,
  height = 7,
  dpi = 300
)















