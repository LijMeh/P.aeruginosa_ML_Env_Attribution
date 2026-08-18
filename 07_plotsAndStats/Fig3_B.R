# ============================================================================
# Plot stats for Top genes mutation rates
# ============================================================================

# ---- Load Required Packages ----
library(install.load)
install_load("here")

# ---- Load configuration and setup ----
source(here("config.R"))
source(here("00_R_configs", "00_setup.R"))
source(here("00_R_configs", "00_color_palette.R"))
source(here("00_R_configs", "00_functions.R"))

# ---- Data ----

df_pao1 <- load_all_gene_files("00_data/Mutant-Calling-11K-CCS-Top27shap/results/ouputs_for_graphing") %>%
  mutate(reference = "PAO1", .after = 1)

flat_df_annotated_gene_names <- readRDS(
  file = data_here("00_data", "Figure_Processed_Data", "flat_df_annotated_gene_names_SHAP_data.RDS")
)

pa_gene_ids <- read.csv(
  data_here("00_data", "PAO1_blast_result_new.csv")) %>% as.tibble()


gene_info <- readRDS(
  file = data_here("00_data", "Figure_Processed_Data", "gene_variant_PA_Annotations.rds")
)

metadata_mapped <- readRDS(
  file = data_here("00_data", "Figure_Processed_Data", "metadata_mapped.rds")
)

metadata_mapped

stats_df_merged <- df_pao1 %>%
  left_join(metadata_mapped, by = c("Genome" = "Genome_ID")) %>%
  dplyr::filter(!is.na(Niche)) %>%
  filter(Clade == "CladeA") %>%
  left_join(pa_gene_ids, by = c("gene" = "Locus.Tag")) %>%
  mutate(
    Mut_Status = case_when(
      Mutation == "Contig split" ~ "Alternative function", # CALLING CONTIG SPLIT AS ALT FUNCTION FOR NOW.
      TRUE ~ as.character(Mut_Status)
    ),
    Locus.Tag = gene
  ) %>%
  group_by(Genome, gene) %>%
  slice(1) %>%
  ungroup()

# ---- DATA QUALITY CHECK ----
cat("\n")
cat(strrep("=", 80), "\n")
cat("DATA QUALITY CHECK: stats_df_merged (includes contig splits, deduplicated)\n")
cat(strrep("=", 80), "\n")
cat("Total rows: ", nrow(stats_df_merged), "\n", sep = "")
cat("Unique genomes: ", n_distinct(stats_df_merged$Genome), "\n", sep = "")
cat("Unique genes: ", n_distinct(stats_df_merged$gene), "\n", sep = "")
cat("Expected max n per gene: ", n_distinct(stats_df_merged$Genome), "\n", sep = "")
cat(strrep("=", 80), "\n\n")

stats_df_merged_noContigSplits <- df_pao1 %>%
  left_join(metadata_mapped, by = c("Genome" = "Genome_ID")) %>%
  dplyr::filter(!is.na(Niche)) %>%
  filter(Clade == "CladeA") %>%
  left_join(pa_gene_ids, by = c("gene" = "Locus.Tag")) %>%
  mutate(Locus.Tag = gene) %>%  # Preserve Locus.Tag column for duplicate checking
  filter(Mutation != "Contig split")

# ---- DATA QUALITY CHECK ----
cat("\n")
cat(strrep("=", 80), "\n")
cat("DATA QUALITY CHECK: stats_df_merged_noContigSplits\n")
cat(strrep("=", 80), "\n")
cat("Total rows: ", nrow(stats_df_merged_noContigSplits), "\n", sep = "")
cat("Unique genomes: ", n_distinct(stats_df_merged_noContigSplits$Genome), "\n", sep = "")
cat("Unique genes: ", n_distinct(stats_df_merged_noContigSplits$gene), "\n", sep = "")
cat("Expected max n per gene: ", n_distinct(stats_df_merged_noContigSplits$Genome), "\n", sep = "")
cat("\nObservations per gene (first 20):\n")
per_gene_check <- stats_df_merged_noContigSplits %>%
  group_by(gene) %>%
  summarise(
    n_rows = n(),
    n_unique_genomes = n_distinct(Genome),
    .groups = "drop"
  ) %>%
  arrange(desc(n_rows))
print(head(per_gene_check, 20), n = Inf)
cat("\nRows per genome per gene (checking for duplicates - first 10):\n")
per_genome_gene_check <- stats_df_merged_noContigSplits %>%
  group_by(Genome, gene) %>%
  summarise(
    n_rows = n(),
    .groups = "drop"
  ) %>%
  filter(n_rows > 1)
if (nrow(per_genome_gene_check) > 0) {
  cat("WARNING: Found ", nrow(per_genome_gene_check), " genome-gene combinations with >1 row (duplicates)\n", sep = "")
  print(head(per_genome_gene_check, 10))
} else {
  cat("✓ No duplicates found (each genome-gene pair appears once)\n")
}
cat(strrep("=", 80), "\n\n")

# ---- INSPECT DUPLICATE ROWS ----
cat("\n")
cat(strrep("=", 80), "\n")
cat("DEDUPLICATING: Keep only first row per Genome-Gene pair\n")
cat(strrep("=", 80), "\n\n")

# First, verify that all duplicates have identical values for key columns
cat("SAFETY CHECK: Verifying duplicate rows have identical Mutation and Mut_Status...\n\n")
dup_check_cols <- stats_df_merged_noContigSplits %>%
  group_by(Genome, gene) %>%
  filter(n() > 1) %>%
  summarise(
    n_rows = n(),
    n_unique_mutations = n_distinct(Mutation, na.rm = TRUE),
    n_unique_mut_status = n_distinct(Mut_Status, na.rm = TRUE),
    n_unique_niches = n_distinct(Niche, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_unique_mutations > 1 | n_unique_mut_status > 1 | n_unique_niches > 1)

if (nrow(dup_check_cols) > 0) {
  cat("⚠️  WARNING: Found duplicates with DIFFERENT Mutation/Mut_Status/Niche values!\n")
  cat("Cannot safely deduplicate without losing information:\n\n")
  print(dup_check_cols, n = Inf)
  stop("Duplicates have inconsistent data - cannot filter safely")
} else {
  cat("✓ PASSED: All duplicate rows have identical Mutation, Mut_Status, and Niche values\n")
  cat("   Safe to deduplicate by keeping only first row per Genome-Gene pair\n\n")
}

# Deduplicate: keep only first row per Genome-Gene pair
stats_df_merged_noContigSplits <- stats_df_merged_noContigSplits %>%
  group_by(Genome, gene) %>%
  slice(1) %>%
  ungroup()

# Verify deduplication worked
cat("DEDUPLICATION RESULTS:\n")
cat("Total rows: ", nrow(stats_df_merged_noContigSplits), "\n", sep = "")
cat("Unique genomes: ", n_distinct(stats_df_merged_noContigSplits$Genome), "\n", sep = "")
cat("Unique genes: ", n_distinct(stats_df_merged_noContigSplits$gene), "\n", sep = "")
per_genome_gene_check_after <- stats_df_merged_noContigSplits %>%
  group_by(Genome, gene) %>%
  summarise(n_rows = n(), .groups = "drop") %>%
  filter(n_rows > 1)
if (nrow(per_genome_gene_check_after) > 0) {
  cat("⚠️  ERROR: Still have ", nrow(per_genome_gene_check_after), " duplicate genome-gene pairs!\n", sep = "")
  stop("Deduplication failed")
} else {
  cat("✓ SUCCESS: Each genome-gene pair appears exactly once\n")
}
cat(strrep("=", 80), "\n\n")

mutation_counts_phzD1 <- stats_df_merged %>%
  filter(Gene.Name %in% c("phzD1")) %>%
  group_by(gene, Gene.Name, Mutation, Mut_Status, Niche) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  ) %>%
  arrange(gene, desc(Count)) %>%
  filter(Niche %in% "Bronchiectasis")

print(mutation_counts_phzD1, n = 20)
# A tibble: 10 × 6
#   gene   Gene.Name Mutation            Mut_Status           Niche         Count
#   <chr>  <chr>     <chr>               <chr>                <chr>         <int>
# 1 PA4213 phzD1     Contig split        Alternative function Bronchiectas…   647
# 2 PA4213 phzD1     WT                  Functional           Bronchiectas…   503
# 3 PA4213 phzD1     Full Deletion       No function          Bronchiectas…   363
# 4 PA4213 phzD1     R99C                Alternative function Bronchiectas…    57
# 5 PA4213 phzD1     A112S               Alternative function Bronchiectas…    32
# 6 PA4213 phzD1     V56L, P102A         Alternative function Bronchiectas…    14

mutation_counts_tonB2 <- stats_df_merged %>%
  filter(Gene.Name %in% c("tonB2")) %>%
  group_by(gene, Gene.Name, Mutation, Mut_Status, Niche) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  ) %>%
  arrange(gene, desc(Count)) %>%
  filter(Niche %in% "Pediatric CF")

print(mutation_counts_tonB2, n = 20)
#gene   Gene.Name Mutation                              Mut_Status Niche Count
#   <chr>  <chr>     <chr>                                 <chr>      <chr> <int>
# 1 PA0197 tonB2     Contig split                          Alternati… Pedi…   342
# 2 PA0197 tonB2     WT                                    Functional Pedi…    44
# 3 PA0197 tonB2     Deletion 170 to 175                   Alternati… Pedi…    10
# 4 PA0197 tonB2     V119I                                 Alternati… Pedi…    10
# 5 PA0197 tonB2     Q161P                                 Alternati… Pedi…     8


# ---- Run Analysis and Plot ----

top_genes <- c("PA0763", "PA2020", "PA0426", "PA1430", "PA0197", "PA4213")

glmer_results <- run_glmer_analysis(
  df_input = stats_df_merged,
  reference_name = "PAO1",
  random_effect_var = c("Genomovar_Cluster"),
  niche_ref = "Adult CF",
  mut_wt_values = c("Functional"),
  mut_case_values = NULL,
  genes_to_analyze = top_genes
)

# Generate statistical methods summary (for Methods section)
generate_statistical_methods_summary(
  glmer_results,
  stats_df_merged,
  exact_n_export_path = data_here(
    "00_data", "Figures",
    "Fig3_TopGenes_ExactNPerGroup.csv"
  )
)

# Format and print publication-ready statistics tables
format_publication_stats(
  glmer_results$coef_flat,
  export_path = data_here(
    "00_data", "Figures",
    "Fig3_TopGenes_PublicationStats.csv"
  )
)

format_pairwise_comparison_table(
  glmer_results$interpretation_summary,
  export_path = data_here(
    "00_data", "Figures",
    "Fig3_TopGenes_PairwiseComparisons.csv"
  )
)

# Define niche colors
niche_colors <- c(
    "Aquatic" = "#9CBED2",
  "Hospital" = "#727272",
  "Adult CF" = "#F3BA36",
  "Pediatric CF" = "#F7D074",
  "Bronchiectasis" = "#E36B2E",
  "Blood" = "#7F4348",
  "Urinary" = "#A8441B",
  "Fecal" = "#B89C8D",
  "Wound" = "#C57066"
)

cld_plot <- plot_mutation_probabilities_cld(
  glmer_results = glmer_results,
  genes = top_genes,
  niche_colors = niche_colors,
  ncol = 3
)

# ---- Save Plot ----

ggsave(
  filename = data_here("00_data", "Figures", "Fig3_TopGenes_MutationProbabilities_CLD.png"),
  plot = cld_plot$plot,
  width = 7,
  height = 4,
  dpi = 300
)

# ---- No Contig Splits ----

glmer_results_noContigSplits <- run_glmer_analysis(
  df_input = stats_df_merged_noContigSplits,
  reference_name = "PAO1",
  random_effect_var = c("Genomovar_Cluster"),
  niche_ref = "Adult CF",
  mut_wt_values = c("Functional"),
  mut_case_values = NULL,
  genes_to_analyze = top_genes
)

# Statistical summary for no contig splits analysis
generate_statistical_methods_summary(
  glmer_results_noContigSplits,
  stats_df_merged_noContigSplits,
  exact_n_export_path = data_here(
    "00_data", "Figures",
    "Fig3_TopGenes_NoContigSplits_ExactNPerGroup.csv"
  )
)

format_publication_stats(
  glmer_results_noContigSplits$coef_flat,
  export_path = data_here(
    "00_data", "Figures",
    "Fig3_TopGenes_NoContigSplits_PublicationStats.csv"
  )
)

format_pairwise_comparison_table(
  glmer_results_noContigSplits$interpretation_summary,
  export_path = data_here(
    "00_data", "Figures",
    "Fig3_TopGenes_NoContigSplits_PairwiseComparisons.csv"
  )
)

noContigSplits_cld_plot <- plot_mutation_probabilities_cld(
  glmer_results = glmer_results_noContigSplits,
  genes = top_genes,
  niche_colors = niche_colors,
  ncol = 3
)

ggsave(
  filename = data_here("00_data", "Figures", "Fig3_TopGenes_NoContigSplits_MutationProbabilities_CLD.png"),
  plot = noContigSplits_cld_plot$plot,
  width = 7,
  height = 4,
  dpi = 300
)

# ---- Loss-of-function-only GLMM + CLD ----

stats_df_merged_lof <- stats_df_merged %>%
  mutate(
    Mut_Status = case_when(
      Mut_Status == "No function" ~ "No function",
      Mut_Status != "No function" ~ "Functional",
      TRUE ~ NA_character_
    )
  )

# ---- DATA QUALITY CHECK ----
cat("\n")
cat(strrep("=", 80), "\n")
cat("DATA QUALITY CHECK: stats_df_merged_lof (Loss of function reclassified, deduplicated)\n")
cat(strrep("=", 80), "\n")
cat("Total rows: ", nrow(stats_df_merged_lof), "\n", sep = "")
cat("Unique genomes: ", n_distinct(stats_df_merged_lof$Genome), "\n", sep = "")
cat("Unique genes: ", n_distinct(stats_df_merged_lof$gene), "\n", sep = "")
cat("Expected max n per gene: ", n_distinct(stats_df_merged_lof$Genome), "\n", sep = "")
cat(strrep("=", 80), "\n\n")


lof_glmer_results <- run_glmer_analysis(
  df_input = stats_df_merged_lof,
  reference_name = "PAO1",
  random_effect_var = c("Genomovar_Cluster"),
  niche_ref = "Adult CF",
  mut_wt_values = c("Functional"),
  mut_case_values = c("No function"),
  genes_to_analyze = top_genes
)

# Statistical summary for loss-of-function analysis
generate_statistical_methods_summary(
  lof_glmer_results,
  stats_df_merged_lof,
  exact_n_export_path = data_here(
    "00_data", "Figures",
    "Fig3_TopGenes_LossOfFunction_ExactNPerGroup.csv"
  )
)

# Publication-ready tables for LOF analysis
format_publication_stats(
  lof_glmer_results$coef_flat,
  export_path = data_here(
    "00_data", "Figures",
    "Fig3_TopGenes_LossOfFunction_PublicationStats.csv"
  )
)

format_pairwise_comparison_table(
  lof_glmer_results$interpretation_summary,
  export_path = data_here(
    "00_data", "Figures",
    "Fig3_TopGenes_LossOfFunction_PairwiseComparisons.csv"
  )
)

lof_cld_plot <- plot_mutation_probabilities_cld(
  glmer_results = lof_glmer_results,
  genes = top_genes,
  niche_colors = niche_colors,
  ncol = 3
)

lof_zero_mut_counts <- purrr::imap_dfr(lof_glmer_results$results_list, function(res, gene_id) {
  if (is.null(res) || is.null(res$status) || res$status != "success" || is.null(res$niche_counts)) {
    return(NULL)
  }

  res$niche_counts %>%
    transmute(
      Gene = as.character(gene_id),
      Niche = as.character(Niche),
      n_mutated_niche = as.numeric(n_mutated_niche)
    )
})

lof_niche_levels <- levels(lof_cld_plot$data$Niche)

lof_cld_plot$data <- lof_cld_plot$data %>%
  left_join(lof_zero_mut_counts, by = c("Gene", "Niche")) %>%
  mutate(
    Lower = ifelse(!is.na(n_mutated_niche) & n_mutated_niche == 0, NA_real_, Lower),
    Upper = ifelse(!is.na(n_mutated_niche) & n_mutated_niche == 0, NA_real_, Upper),
    Lower = ifelse(is.finite(Lower), Lower, NA_real_),
    Upper = ifelse(is.finite(Upper), Upper, NA_real_)
  ) %>%
  select(-n_mutated_niche) %>%
  mutate(Niche = factor(Niche, levels = lof_niche_levels))
lof_cld_plot$plot <- lof_cld_plot$plot %+% lof_cld_plot$data

ggsave(
  filename = data_here("00_data", "Figures", "Fig3_TopGenes_LossOfFunction_MutationProbabilities_CLD.png"),
  plot = lof_cld_plot$plot,
  width = 7,
  height = 4,
  dpi = 300
)

lof_coef_summary <- lof_glmer_results$coef_flat %>%
  filter(Gene %in% top_genes, term != "(Intercept)") %>%
  arrange(Gene, term)

lof_interpretation_summary <- lof_glmer_results$interpretation_summary %>%
  filter(Gene %in% top_genes) %>%
  arrange(Gene, Comparison)

print(lof_coef_summary, n = Inf)
print(lof_interpretation_summary, n = Inf)

# ---- Print Percentage loss of function per gene per env ----

gene_mutation_summary <- stats_df_merged %>%
  filter(gene %in% top_genes) %>%
  group_by(gene, Niche, Gene.Name) %>%
  summarise(
    Total = n(),
    Mutated = sum(Mut_Status == "No function", na.rm = TRUE),
    Mutation_Rate = Mutated / Total,
    .groups = "drop"
  ) %>%
  mutate(
    gene = factor(gene, levels = top_genes),
    Gene_Label = ifelse(is.na(Gene.Name) | Gene.Name == "", as.character(gene), as.character(Gene.Name))
  ) %>%
  arrange(gene, Niche) %>%
  as_tibble()

max_genomes <- n_distinct(stats_df_merged$Genome)
cat("\nTotal max genomes: ", max_genomes, "\n", sep = "")

total_lines_per_gene <- gene_mutation_summary %>%
  group_by(gene, Gene_Label) %>%
  summarise(total_lines = sum(Total, na.rm = TRUE), .groups = "drop") %>%
  arrange(gene)

cat("Total lines per gene (sum across environments):\n")
print(total_lines_per_gene, n = Inf)

mutation_rate_plot <- ggplot(gene_mutation_summary, aes(x = Niche, y = Mutation_Rate, fill = Niche)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  facet_wrap(~ Gene_Label, ncol = 3) +
  scale_fill_manual(values = niche_colors, drop = FALSE) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = ggplot2::expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Environment",
    y = "Observed loss-of-function rate"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = data_here("00_data", "Figures", "Fig3_TopGenes_ObservedLossOfFunction.png"),
  plot = mutation_rate_plot,
  width = 7,
  height = 4,
  dpi = 300
)





