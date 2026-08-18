# ============================================================================
# Plot stats for Top genes mutation rates
# Includes stats for Sup Figs 7 & 9
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
# NOT PROCESSED YET, NEED TO MOVE TO `06_Process_Data_For_Figures.RMD`

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
    )
  )

stats_df_merged_noContigSplits <- df_pao1 %>%
  left_join(metadata_mapped, by = c("Genome" = "Genome_ID")) %>%
  dplyr::filter(!is.na(Niche)) %>%
  filter(Clade == "CladeA") %>%
  left_join(pa_gene_ids, by = c("gene" = "Locus.Tag")) %>%
  filter(Mutation != "Contig split")

# Summarize counts of individual mutations per gene
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

glmer_results <- run_glmer_analysis(
  df_input = stats_df_merged,
  reference_name = "PAO1",
  random_effect_var = c("Genomovar_Cluster"),
  niche_ref = "Adult CF",
  mut_wt_values = c("Functional"),
  mut_case_values = NULL,
  genes_to_analyze = NULL
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

top_genes <- c("PA0763", "PA2020", "PA0426", "PA1430", "PA0197", "PA4213")

cld_plot <- plot_mutation_probabilities_cld(
  glmer_results = glmer_results,
  genes = top_genes,
  niche_colors = niche_colors,
  ncol = 3
)

# ---- Save Plot ----

ggsave(
  filename = here("00_data", "Figures", "Fig3_TopGenes_MutationProbabilities_CLD.png"),
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

noContigSplits_cld_plot <- plot_mutation_probabilities_cld(
  glmer_results = glmer_results_noContigSplits,
  genes = top_genes,
  niche_colors = niche_colors,
  ncol = 3
)

ggsave(
  filename = here("00_data", "Figures", "Fig3_TopGenes_NoContigSplits_MutationProbabilities_CLD.png"),
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
    ))

lof_glmer_results <- run_glmer_analysis(
  df_input = stats_df_merged_lof,
  reference_name = "PAO1",
  random_effect_var = c("Genomovar_Cluster"),
  niche_ref = "Adult CF",
  mut_wt_values = c("Functional"),
  mut_case_values = c("No function"),
  genes_to_analyze = top_genes
)

lof_cld_plot <- plot_mutation_probabilities_cld(
  glmer_results = lof_glmer_results,
  genes = top_genes,
  niche_colors = niche_colors,
  ncol = 3
)

ggsave(
  filename = here("00_data", "Figures", "Fig3_TopGenes_LossOfFunction_MutationProbabilities_CLD.png"),
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
  filename = here("00_data", "Figures", "Fig3_TopGenes_ObservedLossOfFunction.png"),
  plot = mutation_rate_plot,
  width = 7,
  height = 4,
  dpi = 300
)

# ============================================================================
# EXPORT TABLES FOR SUPPLEMENTARY MATERIALS
# ============================================================================

# ---- Extract Tables for All Analyses ----

# Main analysis (all data)
model_details_all <- extract_model_details(glmer_results, stats_df_merged, "All_Data")
n_per_group_all <- extract_exact_n_per_group(stats_df_merged, "All_Data")
pairwise_stats_all <- extract_pairwise_stats(glmer_results, stats_df_merged, "All_Data")

# No contig splits analysis
model_details_noContigSplits <- extract_model_details(glmer_results_noContigSplits, stats_df_merged_noContigSplits, "No_Contig_Splits")
n_per_group_noContigSplits <- extract_exact_n_per_group(stats_df_merged_noContigSplits, "No_Contig_Splits")
pairwise_stats_noContigSplits <- extract_pairwise_stats(glmer_results_noContigSplits, stats_df_merged_noContigSplits, "No_Contig_Splits")

# Loss-of-function analysis
model_details_lof <- extract_model_details(lof_glmer_results, stats_df_merged_lof, "Loss_Of_Function")
n_per_group_lof <- extract_exact_n_per_group(stats_df_merged_lof, "Loss_Of_Function")
pairwise_stats_lof <- extract_pairwise_stats(lof_glmer_results, stats_df_merged_lof, "Loss_Of_Function")

# ---- Save export objects as RDS (separate files per analysis) ----

# All Data analysis
saveRDS(
  model_details_all,
  file = here("00_data", "Figure_Processed_Data", "model_details_All_Data.RDS")
)
saveRDS(
  n_per_group_all,
  file = here("00_data", "Figure_Processed_Data", "n_per_group_All_Data.RDS")
)
saveRDS(
  pairwise_stats_all,
  file = here("00_data", "Figure_Processed_Data", "pairwise_stats_All_Data.RDS")
)

# No Contig Splits analysis
saveRDS(
  model_details_noContigSplits,
  file = here("00_data", "Figure_Processed_Data", "model_details_No_Contig_Splits.RDS")
)
saveRDS(
  n_per_group_noContigSplits,
  file = here("00_data", "Figure_Processed_Data", "n_per_group_No_Contig_Splits.RDS")
)
saveRDS(
  pairwise_stats_noContigSplits,
  file = here("00_data", "Figure_Processed_Data", "pairwise_stats_No_Contig_Splits.RDS")
)

# Loss-of-Function analysis
saveRDS(
  model_details_lof,
  file = here("00_data", "Figure_Processed_Data", "model_details_Loss_Of_Function.RDS")
)
saveRDS(
  n_per_group_lof,
  file = here("00_data", "Figure_Processed_Data", "n_per_group_Loss_Of_Function.RDS")
)
saveRDS(
  pairwise_stats_lof,
  file = here("00_data", "Figure_Processed_Data", "pairwise_stats_Loss_Of_Function.RDS")
)



