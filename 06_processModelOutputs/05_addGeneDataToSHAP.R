# ============================================================================
# Process extra data for SHAP genes
# ============================================================================
# For the genes that are going to be used in the rest of the analysis, 
# which are the gene set produced by 02_Proccess_SHAP_Breadth_Bias.R, I am 
# making sure that all extra data (ie. annotation etc.) is added to them
# so that I can make my final figures. 

# Note: I'm using the SHAP results from the blocked model, as I beleive its 
# better to leave the test data ONLY for validation, not for analysis. 
# ============================================================================

# ---- Load Required Packages ----
library(install.load)
install_load("here")

source(here("config.R"))
source(here("config.R"))
source(here("Iterative_Model_Approach", "00_setup.R"))
source(here("Iterative_Model_Approach", "00_functions.R"))

# ---- Load Data ----

results <- readRDS("00_data/model_output/niche_results.rds")

# Flatten all combined_data from results into one dataframe
results_flat <- lapply(names(results), function(nm) {
  results[[nm]]$combined_data %>%
    mutate(Env = nm, .after = gene)
}) %>% bind_rows()

gene_presence_absence_summary <- read.csv(
  data_here("00_Panaroo_100_UCF", "Model_Input", "gene_presence_absence_summary.csv")
) %>% as_tibble()

# ---- Join extra data to SHAP results ----
results_flat_annotated <- results_flat %>%
  left_join(gene_presence_absence_summary, by = c("gene" = "Gene"))

# ---- Process cluster info ----

# For each cluster, get Annotation from the gene with the longest protein length
annotation_inferred_df <- results_flat_annotated %>%
  filter(!is.na(aa_cluster)) %>%
  group_by(aa_cluster) %>%
  filter(aa_length == max(aa_length, na.rm = TRUE)) %>%
  slice(1) %>%  # in case of ties, take the first
  ungroup() %>%
  select(aa_cluster, Annotation) %>%
  rename(Annotation_Inferred = Annotation)

# Join annotation_inferred_df
results_flat_clustered <- results_flat_annotated %>%
  left_join(annotation_inferred_df, by = "aa_cluster")


# Remove columns with .x/.y suffixes if they exist
dup_cols <- names(results_flat_clustered)[stringr::str_ends(names(results_flat_clustered), "\\.x") | 
                                                  stringr::str_ends(names(results_flat_clustered), "\\.y")]
if (length(dup_cols) > 0) {
  base_names <- unique(sub("[.]x$|[.]y$", "", dup_cols))
  for (bn in base_names) {
    x_col <- paste0(bn, ".x")
    y_col <- paste0(bn, ".y")
    if (x_col %in% names(results_flat_clustered)) {
      results_flat_clustered[[bn]] <- results_flat_clustered[[x_col]]
    } else if (y_col %in% names(results_flat_clustered)) {
      results_flat_clustered[[bn]] <- results_flat_clustered[[y_col]]
    }
  }
  results_flat_clustered <- results_flat_clustered %>% select(-all_of(dup_cols))
}

# Infer gene name by cluster: for each cluster, if any gene has a non-NA Non.unique.Gene.name, assign to all in cluster
if ("Non.unique.Gene.name" %in% colnames(results_flat_clustered)) {
  gene_name_map <- results_flat_clustered %>%
    filter(!is.na(aa_cluster)) %>%
    group_by(aa_cluster) %>%
    summarise(gene_name_cluster_inferred = first(na.omit(Non.unique.Gene.name)), .groups = "drop")
  
  results_flat_clustered <- results_flat_clustered %>%
    left_join(gene_name_map, by = "aa_cluster")
}

results_flat_clustered %>% filter(Env == "Adult CF") %>% arrange(desc(mean_abs_shap)) %>% head(10)

results_flat_clustered %>% colnames()

# ---- Export ----

saveRDS(
  results_flat_clustered,
  data_here("00_data/model_output", "shap_genes_annotated_blocked_model.rds")
)

results_flat_clustered_ordered <- results_flat_clustered %>% arrange(desc(mean_abs_shap))

write.csv(
  results_flat_clustered_ordered,
  data_here("00_data/model_output", "shap_genes_annotated_blocked_model.csv"),
  row.names = FALSE
)

# ---- Get Genes of Intrest for Cluster Validation ----
results_flat_clustered_ordered %>%
    select(gene) %>%
    distinct() %>%
    write.csv(
      data_here("00_data/model_output", "shap_genes_of_intrest_blocked_model.csv"),
      row.names = FALSE
    )

# ---- Output Flat, Non-Redundant Gene List with average SHAP per gene ----

results_flat_clustered %>%
  group_by(gene) %>%
  summarise(
    mean_abs_shap_avg = mean(mean_abs_shap, na.rm = TRUE),
    aa_length = first(aa_length),
    longest_aa = first(longest_aa),
    aa_cluster = first(aa_cluster),
    Annotation_Inferred = first(Annotation_Inferred),
    gene_name_cluster_inferred = first(gene_name_cluster_inferred),
    Presence_Count_AllClasses = first(Presence_Count_AllClasses),
    Absence_Count_AllClasses = first(Absence_Count_AllClasses),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_shap_avg)) %>%
  write.csv(
    data_here("00_data/model_output", "shap_genes_avg_across_envs_blocked_model.csv"),
    row.names = FALSE
  )
