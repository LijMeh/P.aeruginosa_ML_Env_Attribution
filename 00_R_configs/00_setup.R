# Set Posit Package Manager for faster installations on Ubuntu
options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/latest'))

# Functions
library(install.load)

install_load("here",
             "conflicted",
             "parallel",
             "tidyverse",
             "tidyr",
             "data.table",
             "dplyr",
             "patchwork",
             "ggplot2",
             "ggview",
             "ggnewscale",
             "ggimage",
             "ggforce",
             "ggrepel",
             "ggsankey",
             "ape",
             "tidytree",
             "treeio",
             "ggtree",
             "ggtreeExtra",
             "devtools",
             "knitr",
             "effsize",
             "reshape2",
             "ggraph",
             "tidygraph",
             "igraph",
             "scales",
             "AnnotationHub",
             "stringr",
             "ggtext",
             "ggpubr",
             "ggh4x",
             "scatterpie",
             "phytools",
             "phangorn", 
             "TreeTools",
             "phytools",
             "cowplot",
             "graphlayouts",
             "ggExtra",
             "ggside",
             "vegan",
             "ggpubr",
             "pheatmap",
             "gridExtra",
             "circlize",
             "RColorBrewer",
             "ComplexHeatmap",
             "kableExtra",
             "gghalves",
             "ComplexUpset",
             "ggrain",
             "ggVennDiagram",
             "flextable",
             "yardstick",
             "xgboost",
             "progressr",
             "future",
              "furrr",
              "data.table",
              "pROC",
              "ggpattern",
              "emmeans",
              "multcomp",
              "multcompView",
              "lme4"           
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("ggtree", quietly = TRUE)) {
  BiocManager::install("ggtree")
}

if (!requireNamespace("ggtreeExtra", quietly = TRUE)) {
  BiocManager::install("ggtreeExtra")
}

if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
  BiocManager::install("AnnotationHub")
}

if (!requireNamespace("ggsankey", quietly = TRUE)) {
  devtools::install_github("davidsjoberg/ggsankey")
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  devtools::install_github("raivokolde/pheatmap")
}

# Set Defaults

suppressMessages(conflict_prefer_all("dplyr"))
suppressMessages(conflict_prefer_all("tidyverse"))
conflict_prefer("make_long", "ggsankey")
conflict_prefer("merge", "base")
conflicts_prefer(purrr::map)

save_ggplots <- function(plots, prefix = "plot", out_dir = NULL) {
  # Ensure output directory exists
  out_dir <- here::here("Figures", out_dir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Helper to save a single plot
  save_one_plot <- function(plot, name) {
    ggsave(
      filename = file.path(out_dir, paste0(name, ".png")),
      plot = plot,
      width = 6,
      height = 4,
      dpi = 300
    )
  }
  
  # If input is a list, save each with its name
  if (is.list(plots)) {
    for (nm in names(plots)) {
      save_one_plot(plots[[nm]], nm)
    }
  } else {
    # Single plot
    save_one_plot(plots, prefix)
  }
}
