# ============================================================================
# ---- Experimental Data Tree ----
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
filtered_experimental_tree <- readRDS(
  data_here("00_data", "Figure_Processed_Data", "filtered_experimental_tree.RDS")
)

filtered_experimental_metadata <- readRDS(
  data_here("00_data", "Figure_Processed_Data", "filtered_experimental_metadata.RDS")
)

# ---- Prepare Metadata for Tree Tips ----

tree_metadata <- filtered_experimental_metadata %>%
  filter(treeNames %in% filtered_experimental_tree$tip.label) %>%
  mutate(Env = recode(Env, "Non-CF Human Infection" = "Acute Human Infection")) %>%
  select(treeNames, Env, Clade) %>%
  as.data.frame()

rownames(tree_metadata) <- tree_metadata$treeNames

# ---- Prepare Bar Plot Data ----
env_colors <- c(
  "Cystic Fibrosis"                     = "#F3BA36",
  "Acute Human Infection"  = "#B6502B",
  "Environmental"            = "#9CA780"
)


cols <- c("Elastase Production", "Biofilm Fraction", "Growth Rate")

bar_data <- filtered_experimental_metadata %>%
  filter(treeNames %in% filtered_experimental_tree$tip.label) %>%
  select(treeNames, all_of(cols))

# Scale bar: 0.001 = 0.1% sequence divergence (1 substitution per 1000 sites)
scale_width <- 0.001

p.tree <- ggtree(filtered_experimental_tree, layout = "rectangular") %<+% tree_metadata +
  geom_tippoint(aes(color = Env), size = 2) +
  scale_color_manual(
    values = env_colors,
    breaks = c("Cystic Fibrosis", "Acute Human Infection", "Environmental")
  ) +
  scale_fill_manual(
    values = env_colors,
    breaks = c("Cystic Fibrosis", "Acute Human Infection", "Environmental")
  ) +
  geom_treescale(
    x = 0,
    y = -2,
    width = scale_width,
    linesize = 0.75,
    offset = 1,
    fontsize = 3
  ) +
  geom_fruit(
    data = bar_data,
    geom = geom_col,
    mapping = aes(y = treeNames, x = `Elastase Production`, fill = Env),
    orientation = "y",
    axis.params = list(axis = "x", title = "\nElastase\nProduction\nA₄₉₅/OD₆₀₀"),
    offset = 0.05
  ) +
  geom_fruit(
    data = bar_data,
    geom = geom_col,
    mapping = aes(y = treeNames, x = `Growth Rate`, fill = Env),
    orientation = "y",
    axis.params = list(axis = "x", title = "\nGrowth\nRate\nµ(h⁻¹)"),
    offset = 0.05
  ) +
  geom_fruit(
    data = bar_data,
    geom = geom_col,
    mapping = aes(y = treeNames, x = `Biofilm Fraction`, fill = Env),
    orientation = "y",
    axis.params = list(axis = "x", title = "\nBiofilm\n Fraction"),
    offset = 0.05
  )

# ---- Save experimental Tree Plot ----
ggsave(
  filename = data_here("00_data", "Figures", "Fig5_experimentalTree_v2.png"),
  plot = p.tree,
  width = 8,
  height = 5,
  dpi = 300
)


