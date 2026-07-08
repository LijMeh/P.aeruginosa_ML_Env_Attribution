# ============================================================================
# Process SHAP values, Breadth, and Study Bias for Blocked Models
# ============================================================================
# This script processes SHAP values, breadth, and study bias metrics to determine 
# which features to use for the final model.  
# ============================================================================

# ============================================================================
# Prep Data 
# ============================================================================

## ---- Load Required Packages ----
library(install.load)
install_load("here")

source(here("00_R_configs", "config.R"))
source(here("00_R_configs", "config.R"))
source(here("00_R_configs", "00_setup.R"))

## ---- Load Data Files ----

metadata <- read.csv(data_here("00_data", "blockedModelInputMetadata.csv"))

shap_dir <- data_here("00_data", "model_output", "SHAP_Values")
model_input_path <- data_here("00_data", "model_input", "ModelInput_NoDuplicates.csv")

## ---- Define Niche remapping ----
niche_map <- c(
  "Bronchiectasis" = "Bronchiectasis",
  "CF" = "Adult CF",
  "early.CF" = "Pediatric CF",
  "Early CF" = "Pediatric CF",
  "Urinary" = "Urinary",
  "Gastrointestinal" = "Gastrointestinal",
  "Wound" = "Wound",
  "Blood" = "Blood",
  "Aquatic" = "Aquatic",
  "Ocean" = "Aquatic",
  "Hospital Environment" = "Hospital",
  "Hospital.Environment" = "Hospital",
  "Hospital.Enviornment" = "Hospital",
  "Rectal/Feces" = "Fecal",
  "Rectal-Feces" = "Fecal"
)

remap_niches <- function(niche_vector) {
  recode(niche_vector, !!!niche_map, .default = NA_character_)
}

# ============================================================================
# Define Helper Functions For SHAP, Breadth, and Study Bias Calculation
# ============================================================================

## ---- Filter Model Input Data ----
load_model_input_for_genes <- function(model_input_path, genes, genomes_to_keep) {
  cat(sprintf("Loading model input for %d genes and %d genomes\n", 
              length(genes), length(genomes_to_keep)))
  
  model_data <- fread(
    model_input_path,
    select = c("Genome_ID", genes)
  ) %>%
    filter(Genome_ID %in% genomes_to_keep)
  
  cat(sprintf("Loaded model input: %d genomes x %d genes\n", 
              nrow(model_data), ncol(model_data) - 1))
  
  return(model_data)
}

## ---- Find SHAP Naming Pattern ----
# Parse filenames to extract niche identifiers
parse_shap_filenames <- function(filenames) {
  tibble(filename = filenames) %>%
    mutate(
      # Extract everything after "niche_" and before ".csv"
      niche_raw = str_extract(filename, "niche_(.+)\\.csv$") %>% 
        str_remove("niche_") %>% 
        str_remove("\\.csv$"),
      # Extract cluster if present (optional - might not always exist)
      cluster = str_extract(filename, "Cluster_(\\d+)") %>% 
        str_remove("Cluster_")
    )
}

## ---- Load all SHAP files for a specific niche ----
load_niche_shap <- function(shap_dir, niche_raw, genomes_to_keep, filter_mode = "common") {
  # filter_mode options:
  #   "common" = only genes in ALL clusters (no NAs)
  #   "union"  = all genes from ANY cluster (creates NAs where gene not evaluated)
  
  # Find all files for this niche
  all_files <- list.files(shap_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Filter to files matching this niche (check basename only)
  niche_pattern <- paste0("niche_", niche_raw, "\\.csv$")
  niche_files <- all_files[str_detect(basename(all_files), niche_pattern)]
  
  if (length(niche_files) == 0) {
    stop(sprintf("No SHAP files found for niche: %s", niche_raw))
  }
  
  cat(sprintf("Loading %d SHAP files for niche: %s\n", length(niche_files), niche_raw))
  
  # Load all files
  shap_list <- lapply(niche_files, function(f) {
    dt <- fread(f, colClasses = list(character = "Genome_ID"))
    dt <- dt[Genome_ID %in% genomes_to_keep]  # KEEP only these genomes
    return(dt)
  })
  
  if (filter_mode == "common") {
    # Find genes present in ALL clusters
    all_cols <- lapply(shap_list, colnames)
    common_genes <- Reduce(intersect, all_cols)
    
    # Calculate how many genes before filtering
    total_genes_before <- length(all_cols[[1]]) - 1  # -1 for Genome_ID
    genes_kept <- length(common_genes) - 1  # -1 for Genome_ID
    genes_removed <- total_genes_before - genes_kept
    
    cat(sprintf("Kept %d genes common to all %d clusters (removed %d cluster-specific genes)\n", 
                genes_kept, length(shap_list), genes_removed))
    
    # Subset each to common genes only
    shap_list <- lapply(shap_list, function(dt) dt[, ..common_genes])
  } else if (filter_mode == "union") {
    all_cols <- lapply(shap_list, colnames)
    all_unique_genes <- unique(unlist(all_cols))
    cat(sprintf("Keeping all %d unique genes across clusters (will create NAs for missing values)\n", 
                length(all_unique_genes) - 1))
  }
  
  # Combine all clusters for this niche
  shap_combined <- rbindlist(shap_list, fill = (filter_mode == "union"))
  
  cat(sprintf("Loaded %d total genomes for niche %s\n", nrow(shap_combined), niche_raw))
  
  return(shap_combined)
}

## ---- Map SHAP niche names to metadata niche names ----
shap_to_metadata_niche <- function(shap_niche) {
  mapping <- c(
    "Early CF" = "early.CF",
    "CF" = "CF",
    "Bronchiectasis" = "Bronchiectasis",
    "Urinary" = "Urinary",
    "Wound" = "Wound",
    "Blood" = "Blood",
    "Aquatic" = "Aquatic",
    "Ocean" = "Ocean",
    "Hospital Environment" = "Hospital.Environment",
    "Rectal-Feces" = "Rectal/Feces"
  )
  
  ifelse(shap_niche %in% names(mapping), mapping[shap_niche], shap_niche)
}

## ---- Prep Niche Metadata ----
prepare_niche_metadata <- function(metadata, niche_remapped, genomes_in_shap) {
  # Subset to genomes that are:
  # 1. Actually IN this niche (for breadth and study bias)
  # 2. Present in the SHAP data (were evaluated)
  
  meta_subset <- metadata %>%
    filter(Niche == niche_remapped, 
           Genome_ID %in% genomes_in_shap)
  
  cat(sprintf("Metadata subset for niche %s: %d genomes\n", 
              niche_remapped, nrow(meta_subset)))
  
  return(meta_subset)
}

## ---- Process SHAP Data ----
calculate_shap_stats <- function(shap_data, model_input) {
  # Get gene names (excluding Genome_ID)
  genes <- setdiff(colnames(shap_data), "Genome_ID")
  
  cat(sprintf("Calculating SHAP stats for %d genes across %d genomes\n", 
              length(genes), nrow(shap_data)))
  
  # Ensure same genome order
  model_input <- model_input[base::match(shap_data$Genome_ID, model_input$Genome_ID), ]
  
  # Calculate stats for each gene
  results <- lapply(genes, function(g) {
    shap_vals <- shap_data[[g]]
    presence <- model_input[[g]]
    
    # Mean SHAP when present (gene = 1)
    mean_shap_pres <- mean(shap_vals[presence == 1], na.rm = TRUE)
    
    # Mean SHAP when absent (gene = 0)
    mean_shap_abs <- mean(shap_vals[presence == 0], na.rm = TRUE)
    
    # Mean of absolute values
    mean_abs_shap <- (abs(mean_shap_pres) + abs(mean_shap_abs)) / 2
    
    data.table(
      gene = g,
      mean_shap_pres = ifelse(is.nan(mean_shap_pres), 0, mean_shap_pres),
      mean_shap_abs = ifelse(is.nan(mean_shap_abs), 0, mean_shap_abs),
      mean_abs_shap = mean_abs_shap
    )
  })
  
  results_dt <- rbindlist(results)
  
  cat(sprintf("Completed SHAP stats for %d genes\n", nrow(results_dt)))
  
  return(results_dt)
}

## ---- Calculate Genomovar Breadth ----
calculate_gene_breadth <- function(model_input, metadata_subset, genes_to_calc) {
  # model_input: gene presence/absence for genomes (already filtered to genomes_in_study)
  # metadata_subset: metadata for genomes ACTUALLY in this niche (from prepare_niche_metadata)
  # genes_to_calc: vector of genes to calculate breadth for (from SHAP data)
  
  # Filter model_input to only genomes in this niche
  model_niche <- model_input %>%
    filter(Genome_ID %in% metadata_subset$Genome_ID)
  
  # Add genomovar cluster info
  model_niche <- model_niche %>%
    left_join(metadata_subset %>% select(Genome_ID, Genomovar_Cluster), 
              by = "Genome_ID")
  
  # Get total number of unique genomovar clusters in this niche
  total_genomovars <- n_distinct(model_niche$Genomovar_Cluster)
  
  cat(sprintf("Calculating breadth for %d genes across %d genomovar clusters (%d genomes)\n",
              length(genes_to_calc), total_genomovars, nrow(model_niche)))
  
  # For each gene, calculate breadth
  results <- lapply(genes_to_calc, function(g) {
    gene_data <- model_niche[[g]]
    genomovar <- model_niche$Genomovar_Cluster
    
    # Group by genomovar and check if gene present/absent at least once
    genomovar_presence <- tapply(gene_data, genomovar, function(x) any(x == 1))
    genomovar_absence <- tapply(gene_data, genomovar, function(x) any(x == 0))
    
    # Breadth = proportion of genomovars with at least one occurrence
    breadth_present <- sum(genomovar_presence, na.rm = TRUE) / total_genomovars
    breadth_absent <- sum(genomovar_absence, na.rm = TRUE) / total_genomovars
    
    data.table(
      gene = g,
      breadth_present = breadth_present,
      breadth_absent = breadth_absent
    )
  })
  
  results_dt <- rbindlist(results)
  
  cat(sprintf("Completed breadth calculation for %d genes\n", nrow(results_dt)))
  
  return(results_dt)
}

## ---- Calculate Study Bias ----
calculate_study_bias <- function(model_input, metadata_subset, genes_to_calc,
                                  high_thresh = 0.5, min_diff = 0.5, min_samples = 10,
                                  potential_low_thresh = 0.1, n_cores = 30) {
  # Filter to genomes in this niche
  model_niche <- model_input %>%
    filter(Genome_ID %in% metadata_subset$Genome_ID) %>%
    left_join(metadata_subset %>% select(Genome_ID, Assembly_Submitter = `Assembly Submitter`),
              by = "Genome_ID")
  
  # Only keep genes we care about + metadata
  genes_present <- intersect(genes_to_calc, colnames(model_niche))
  model_niche <- model_niche %>% select(Genome_ID, Assembly_Submitter, all_of(genes_present))
  
  # Pre-filter to valid submitters
  model_niche <- model_niche %>% filter(!is.na(Assembly_Submitter))
  
  # Convert to data.table for speed
  setDT(model_niche)
  
  cat(sprintf("Calculating study bias for %d genes using %d cores\n", 
              length(genes_present), n_cores))
  
  # Parallelize - each worker only gets gene name, pulls column from shared model_niche
  results <- mclapply(genes_present, function(gene) {
    # Summarize by submitter
    gene_by_submitter <- model_niche[, .(
      n_present = sum(get(gene) == 1, na.rm = TRUE),
      n_total = .N,
      prop_present = sum(get(gene) == 1, na.rm = TRUE) / .N
    ), by = Assembly_Submitter][
      n_total >= min_samples
    ][order(-prop_present)]
    
    # Need at least 2 submitters
    if (nrow(gene_by_submitter) < 2) {
      return(data.table(
        gene = gene,
        submitter_high = NA_character_,
        prop_high = NA_real_,
        n_high = NA_integer_,
        submitter_second = NA_character_,
        prop_second = NA_real_,
        n_second = NA_integer_,
        unevenness = NA_real_,
        study_bias = "No"
      ))
    }
    
    # Calculate unevenness
    prop_high <- gene_by_submitter$prop_present[1]
    prop_second <- gene_by_submitter$prop_present[2]
    unevenness <- prop_high - prop_second
    
    # Determine bias category
    is_biased <- (prop_high >= high_thresh) & (unevenness >= min_diff)
    is_potential <- is_biased & (prop_second > potential_low_thresh)
    
    bias_category <- if (is_potential) {
      "Potential"
    } else if (is_biased) {
      "Biased"
    } else {
      "No"
    }
    
    data.table(
      gene = gene,
      submitter_high = gene_by_submitter$Assembly_Submitter[1],
      prop_high = prop_high,
      n_high = gene_by_submitter$n_total[1],
      submitter_second = gene_by_submitter$Assembly_Submitter[2],
      prop_second = prop_second,
      n_second = gene_by_submitter$n_total[2],
      unevenness = unevenness,
      study_bias = bias_category
    )
  }, mc.cores = n_cores)
  
  rbindlist(results)
}

## ---- Calculate Presence/Absence Counts ----
calc_presence_absence_counts <- function(model_input_subset,
                                        genes_to_calc,
                                        genomes_class = NULL) {
  # model_input_subset: data.table or data.frame with Genome_ID + gene columns (0/1)
  # genes_to_calc: genes to compute for
  # genomes_class: if provided, restrict to these Genome_IDs (class-specific)

  if (!is.null(genomes_class)) {
    model_dt <- model_input_subset[model_input_subset$Genome_ID %in% genomes_class, , drop = FALSE]
  } else {
    model_dt <- model_input_subset
  }

  genes_present <- intersect(genes_to_calc, colnames(model_dt))

  counts_list <- lapply(genes_present, function(g) {
    v <- model_dt[[g]]
    v <- v[v %in% c(0, 1)]  # ignore NA/other
    data.frame(
      gene = g,
      Presence_Count = sum(v == 1, na.rm = TRUE),
      Absence_Count  = sum(v == 0, na.rm = TRUE),
      Total_Count    = length(v)
    )
  })

  dplyr::bind_rows(counts_list)
}

## ---- Process A Niche ----
process_niche <- function(niche, shap_dir, model_input_path, metadata, genomes_in_study,
                          filter_mode = "common", n_cores = 30) {
  cat(sprintf("\n========== Processing niche: %s ==========\n", niche))

  # 1. Load SHAP data for this niche
  shap_data <- load_niche_shap(shap_dir, niche, genomes_in_study, filter_mode = filter_mode)
  genes_in_niche <- setdiff(colnames(shap_data), "Genome_ID")
  cat(sprintf("Loaded %d genes from SHAP data\n", length(genes_in_niche)))

  # 2. Load model input for these genes (only genomes_in_study)
  model_input <- load_model_input_for_genes(model_input_path, genes_in_niche, genomes_in_study)

  # 3. Prepare niche-specific metadata
  niche_remapped <- remap_niches(niche)
  niche_for_metadata <- shap_to_metadata_niche(niche)
  metadata_niche <- prepare_niche_metadata(metadata, niche_for_metadata, shap_data$Genome_ID)
  
  # CHECK: metadata_niche should not be empty
  if (nrow(metadata_niche) == 0) {
    stop(sprintf("ERROR in %s: metadata_niche has 0 rows. Check that '%s' exists in metadata$Niche", 
                 niche, niche_for_metadata))
  }
  
  # CHECK: metadata_niche must have required columns
  required_cols <- c("Genome_ID", "Genomovar_Cluster", "Assembly Submitter")
  missing_cols <- setdiff(required_cols, colnames(metadata_niche))
  if (length(missing_cols) > 0) {
    stop(sprintf("ERROR in %s: metadata_niche missing columns: %s", 
                 niche, paste(missing_cols, collapse = ", ")))
  }

  # 4. Calculate SHAP stats
  shap_stats <- calculate_shap_stats(shap_data, model_input)

  # Ensure join key is 'gene'
  if ("Gene" %in% colnames(shap_stats) && !"gene" %in% colnames(shap_stats)) {
    shap_stats <- dplyr::rename(shap_stats, gene = Gene)
  }

  # 5. Add Net Predictive Value + Consensus
  shap_stats <- shap_stats %>%
    dplyr::mutate(
      Net_Predictive_Value = mean_shap_pres - mean_shap_abs,
      Consensus_Pres_Summary = dplyr::case_when(
        Net_Predictive_Value > 0 ~ "Presence_Expected",
        Net_Predictive_Value < 0 ~ "Absence_Expected",
        TRUE ~ "neutral"
      )
    )

  # 6. Add presence/absence counts (class-specific + global within genomes_in_study)
  counts_class <- lapply(genes_in_niche, function(g) {
    if (!g %in% colnames(model_input)) return(NULL)

    idx <- model_input$Genome_ID %in% metadata_niche$Genome_ID
    v <- model_input[[g]][idx]
    v <- v[v %in% c(0, 1)]

    data.frame(
      gene = g,
      Presence_Count_ClassSpecific = sum(v == 1, na.rm = TRUE),
      Absence_Count_ClassSpecific  = sum(v == 0, na.rm = TRUE)
    )
  }) %>% dplyr::bind_rows()
  
  # CHECK: counts_class should have non-zero counts for at least some genes
  if (all(counts_class$Presence_Count_ClassSpecific == 0 & counts_class$Absence_Count_ClassSpecific == 0)) {
    stop(sprintf("ERROR in %s: All class-specific counts are 0. metadata_niche may have wrong genomes.", niche))
  }

  counts_global <- lapply(genes_in_niche, function(g) {
    if (!g %in% colnames(model_input)) return(NULL)

    v <- model_input[[g]]
    v <- v[v %in% c(0, 1)]

    data.frame(
      gene = g,
      Presence_Count_AllClasses = sum(v == 1, na.rm = TRUE),
      Absence_Count_AllClasses  = sum(v == 0, na.rm = TRUE)
    )
  }) %>% dplyr::bind_rows()

  # INTEGRATE counts into shap_stats
  shap_stats <- shap_stats %>%
    dplyr::left_join(counts_class, by = "gene") %>%
    dplyr::left_join(counts_global, by = "gene")

  # 7. Calculate gene breadth
  gene_breadth <- calculate_gene_breadth(model_input, metadata_niche, genes_in_niche)
  
  # CHECK: breadth should not be all NaN
  if (all(is.nan(gene_breadth$breadth_present)) && all(is.nan(gene_breadth$breadth_absent))) {
    stop(sprintf("ERROR in %s: All breadth values are NaN. Check metadata_niche Genomovar_Cluster.", niche))
  }

  # 8. Calculate study bias
  study_bias <- calculate_study_bias(model_input, metadata_niche, genes_in_niche, n_cores = n_cores)
  
  # CHECK: study_bias should have some non-NA submitter values
  n_with_submitters <- sum(!is.na(study_bias$submitter_high))
  if (n_with_submitters == 0) {
    warning(sprintf("WARNING in %s: No genes have Assembly Submitter data for bias calculation.", niche))
  }

  # 9. Combine all metrics
  combined <- shap_stats %>%
    dplyr::left_join(gene_breadth, by = "gene") %>%
    dplyr::left_join(study_bias, by = "gene")

  cat(sprintf("Completed %s: %d genes with all metrics\n", niche, nrow(combined)))

  return(list(
    niche = niche,
    combined_data = combined,
    n_genes = nrow(combined),
    n_biased = sum(combined$study_bias == "Biased", na.rm = TRUE)
  ))
}

# ============================================================================
# Process all niches
# ============================================================================

## ---- List SHAP files ----
shap_files_list <- list.files(
  path = shap_dir,
  pattern = "\\.csv$",     # all CSVs (tighten if your dir has other CSVs)
  full.names = TRUE,
  recursive = TRUE
)

file_info <- parse_shap_filenames(shap_files_list)

## Get unique raw niches (for loading files)
unique_niches_raw <- file_info %>%
  dplyr::pull(niche_raw) %>%
  unique() %>%
  sort()

niche_results <- setNames(
  lapply(unique_niches_raw, function(niche) {
    process_niche(
      niche = niche,
      shap_dir = shap_dir,                 # (make sure shap_dir == where the SHAP files actually are)
      model_input = model_input_path,      # your process_niche currently accepts a path here
      metadata = metadata,
      genomes_in_study = genomes_in_study,
      filter_mode = "common",
      n_cores = 30
    )
  }),
  unique_niches_raw
)

## ---- Sanity check ----

## First make sure the problem groups look good
niche_results$`Early CF`$combined_data %>% arrange(desc(mean_abs_shap)) %>% head(10) 
niche_results$`Bronchiectasis`$combined_data %>% arrange(desc(mean_abs_shap)) %>% head(10) 

na_summary <- lapply(names(niche_results), function(nm) {
  df <- niche_results[[nm]]$combined_data
  na_counts <- sapply(df, function(col) sum(is.na(col)))
  na_counts <- na_counts[na_counts > 0]
  
  if (length(na_counts) > 0) {
    data.frame(
      niche = nm,
      column = names(na_counts),
      n_na = as.integer(na_counts),
      total = nrow(df),
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
})

na_summary_df <- dplyr::bind_rows(na_summary)

niche_results$`Early CF`$combined_data %>% filter(if_any(everything(), is.na)) %>% select(gene, mean_shap_pres, mean_shap_abs, mean_abs_shap)

# There are 941 genes with NaN for mean_abs_shap, because their mean_shap_pres is 0 (and frequently their mean_shap_abs is also 0)
# Down to 243 after Cluster 1 is added

# ============================================================================
# Filter genes with zero predictive value (neutral in all environments)
# ============================================================================

# First, flatten niche_results to identify neutral genes across all envs
niche_results_flat_for_neutral_check <- lapply(names(niche_results), function(nm) {
  niche_results[[nm]]$combined_data %>%
    mutate(Env = nm, .after = gene)
}) %>% bind_rows()

# Identify genes that are neutral in ALL environments
genes_neutral_all_envs <- niche_results_flat_for_neutral_check %>%
  group_by(gene) %>%
  summarise(
    all_neutral = all(Consensus_Pres_Summary == "neutral"),
    n_envs = n(),
    .groups = "drop"
  ) %>%
  filter(all_neutral) %>%
  pull(gene)

cat(sprintf("\n=== Filtering Neutral Genes ===\n"))
cat(sprintf("Genes with zero predictive value in all environments: %d\n", length(genes_neutral_all_envs)))

# Remove these genes from niche_results
niche_results <- filter_niche_results_globally(niche_results, genes_neutral_all_envs)

# Print summary
for (nm in names(niche_results)) {
  cat(sprintf("%s: %d genes remaining\n", nm, nrow(niche_results[[nm]]$combined_data)))
}


niche_results$`Early CF`$combined_data %>% colnames()

if (nrow(na_summary_df) > 0) {
  print(na_summary_df)
} else {
  cat("No unexpected NAs found in any niche\n")
}

niche_results 

saveRDS(niche_results, data_here("00_data/model_output/", "niche_results_raw.rds"))

# USE THESE for the breadth calculation for comparing to the "full" model
niche_results_flat <- lapply(names(niche_results), function(nm) {
  niche_results[[nm]]$combined_data %>%
    mutate(Env = nm, .after = gene)
}) %>% bind_rows()

saveRDS(niche_results_flat, data_here("00_data/model_output/", "blocked_breadth.rds"))

# ============================================================================
# Helper functions to identify and remove globally biased genes
# ============================================================================

## ---- Identify Globally Biased Genes ----
get_globally_biased_genes <- function(niche_results) {
  # Returns a list with three data frames: biased, potential, and no
  
  biased_df <- dplyr::bind_rows(lapply(names(niche_results), function(nm) {
    x <- niche_results[[nm]]
    df <- x$combined_data
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df %>%
      dplyr::filter(study_bias == "Biased") %>%
      dplyr::transmute(
        gene = gene,
        Biased_In_Niche = nm
      )
  }))
  
  potential_df <- dplyr::bind_rows(lapply(names(niche_results), function(nm) {
    x <- niche_results[[nm]]
    df <- x$combined_data
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df %>%
      dplyr::filter(study_bias == "Potential") %>%
      dplyr::transmute(
        gene = gene,
        Potential_In_Niche = nm
      )
  }))
  
  no_df <- dplyr::bind_rows(lapply(names(niche_results), function(nm) {
    x <- niche_results[[nm]]
    df <- x$combined_data
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df %>%
      dplyr::filter(study_bias == "No") %>%
      dplyr::transmute(
        gene = gene,
        No_In_Niche = nm
      )
  }))
  
  list(
    biased = biased_df %>% dplyr::distinct(gene, Biased_In_Niche) %>% dplyr::arrange(gene, Biased_In_Niche),
    potential = potential_df %>% dplyr::distinct(gene, Potential_In_Niche) %>% dplyr::arrange(gene, Potential_In_Niche),
    no = no_df %>% dplyr::distinct(gene, No_In_Niche) %>% dplyr::arrange(gene, No_In_Niche)
  )
}

parse_cdhit_clstr_gene_map <- function(clstr_path) {
  lines <- readLines(clstr_path)
  
  current_cluster <- NA_integer_
  out_gene <- character(0)
  out_len  <- integer(0)
  out_cl   <- integer(0)
  
  for (line in lines) {
    if (startsWith(line, ">Cluster ")) {
      current_cluster <- as.integer(sub("^>Cluster\\s+", "", line))
      next
    }
    
    # Typical line: "0\t123aa, >gene_name... *"
    m <- regexec("([0-9]+)aa, >([A-Za-z0-9_-]+)", line)
    mm <- regmatches(line, m)[[1]]
    if (length(mm) == 3) {
      out_len  <- c(out_len,  as.integer(mm[2]))
      out_gene <- c(out_gene, mm[3])
      out_cl   <- c(out_cl,   current_cluster)
    }
  }
  
  # Calculate longest_aa per cluster
  df <- tibble(
    gene = out_gene,
    aa_length = out_len,
    aa_cluster = out_cl
  ) %>%
    filter(!is.na(aa_cluster), !is.na(gene)) %>%
    distinct(gene, aa_cluster, .keep_all = TRUE)
  
  # Add longest_aa column
  longest_aa_df <- df %>%
    group_by(aa_cluster) %>%
    summarise(longest_aa = max(aa_length, na.rm = TRUE), .groups = "drop")
  
  df <- df %>%
    left_join(longest_aa_df, by = "aa_cluster")
  
  return(df)
}

expand_biased_to_cluster_mates <- function(biased_genes, potential_genes, gene_cluster_df) {
  # gene_cluster_df must have columns: gene, aa_cluster
  
  # Get clusters containing biased genes
  biased_clusters <- gene_cluster_df %>%
    filter(gene %in% biased_genes) %>%
    pull(aa_cluster) %>%
    unique()
  
  # All genes in those clusters become biased
  cluster_mates <- gene_cluster_df %>%
    filter(aa_cluster %in% biased_clusters) %>%
    pull(gene) %>%
    unique()
  
  # Potential genes that cluster with biased genes are upgraded to biased
  potential_upgraded <- intersect(potential_genes, cluster_mates)
  
  # Final biased set
  final_biased <- sort(unique(c(biased_genes, cluster_mates)))
  
  # Remaining potential genes (not clustered with biased)
  final_potential <- setdiff(potential_genes, final_biased)
  
  list(
    biased = final_biased,
    potential = final_potential,
    potential_upgraded = potential_upgraded
  )
}

filter_niche_results_globally <- function(niche_results, genes_to_remove) {
  lapply(niche_results, function(x) {
    for (slot in c("combined_data", "shap_stats", "gene_breadth", "study_bias")) {
      if (!is.null(x[[slot]]) && is.data.frame(x[[slot]]) && "gene" %in% names(x[[slot]])) {
        x[[slot]] <- x[[slot]] %>% filter(!gene %in% genes_to_remove)
      }
    }
    x
  })
}

# ============================================================================
# Run Global Filtering of Biased Genes
# ============================================================================

niche_results_flat <- readRDS(data_here("00_data/model_output/", "blocked_breadth.rds"))

## ---- Run global filtering ----
clstr_path <- "00_data/feature_clusters/clusters.fasta.clstr"

gene_cluster_df <- parse_cdhit_clstr_gene_map(clstr_path)

# niche_results <- readRDS(data_here("00_data/model_output/", "niche_results_raw.rds"))

globally_biased_list <- get_globally_biased_genes(niche_results)

globally_biased_genes <- sort(unique(globally_biased_list$biased$gene))
globally_potential_genes <- sort(unique(globally_biased_list$potential$gene))

# Expand biased to cluster mates and upgrade potential genes clustered with biased
expanded_results <- expand_biased_to_cluster_mates(globally_biased_genes, globally_potential_genes, gene_cluster_df)

# Remove biased + cluster mates (including upgraded potential genes)
genes_to_remove <- expanded_results$biased

niche_results_filtered <- filter_niche_results_globally(niche_results, genes_to_remove)

# Quick sanity numbers
c(
  n_globally_biased = length(globally_biased_genes),
  n_globally_potential = length(globally_potential_genes),
  n_potential_upgraded_to_biased = length(expanded_results$potential_upgraded),
  n_remaining_potential = length(expanded_results$potential),
  n_removed_including_cluster_mates = length(genes_to_remove),
  n_clustered_genes_total = nrow(gene_cluster_df)
)

niche_results$`Bronchiectasis`$combined_data %>% arrange(desc(mean_abs_shap)) %>% select(gene, study_bias) %>% head(10) 

niche_results$`Early CF`$combined_data %>% arrange(desc(mean_abs_shap)) %>% select(gene, study_bias) %>% head(10) 

niche_results_filtered$`Early CF`$combined_data %>% arrange(desc(mean_abs_shap)) %>% select(gene, study_bias) %>% head(10)

niche_results_filtered$`Bronchiectasis`$combined_data %>% arrange(desc(mean_abs_shap)) %>% select(gene, study_bias) %>% head(10)

## ---- Rename Niches ----
old_names <- names(niche_results_filtered)
new_names <- remap_niches(old_names)

names(niche_results_filtered) <- new_names

niche_results_filtered$`Pediatric CF`$combined_data %>% arrange(desc(mean_abs_shap)) %>% head(10)

niche_results_filtered$`Fecal`$combined_data %>% arrange(desc(mean_abs_shap)) %>% head(10)

## ---- Add Cluster Info to Filtered Results ----
niche_results_filtered_clusters <- lapply(niche_results_filtered, function(nr) {
  nr$combined_data <- nr$combined_data %>%
    dplyr::left_join(updated_cluster_df %>% select(gene, aa_cluster, aa_length, longest_aa), by = "gene")
  nr
})

niche_results_filtered_clusters$Wound$combined_data %>% select(gene, aa_cluster, aa_length, longest_aa) %>% head(10)

# ============================================================================
# Check why certain genes are being removed by bias filtering
# ============================================================================

# Specifically looking at the CCON genes (as they were high importance before)

# Define your genes of interest
genes_of_interest <- c(
  "group_225623", "group_365760", "group_225624", "group_366794", 
  "group_220395", "group_369578", "group_778810", "group_661854", 
  "group_225647", "group_712352", "group_746570", "group_712365", 
  "group_661960", "group_464426", "group_752254", "group_662036", 
  "group_731039", "group_680779"
)

# Filter globally_biased_list by these genes
globally_biased_list_filtered <- list(
  biased = globally_biased_list$biased %>% filter(gene %in% genes_of_interest),
  potential = globally_biased_list$potential %>% filter(gene %in% genes_of_interest),
  no = globally_biased_list$no %>% filter(gene %in% genes_of_interest)
)

# Due to biased genes
#           gene      Biased_In_Niche
#         <char>               <char>
#1: group_225623       Bronchiectasis
#2: group_225624       Bronchiectasis
#3: group_225647                   CF
#4: group_369578 Hospital Environment

niche_results_flat %>%
  filter(gene %in% genes_of_interest) %>%
  arrange(gene, Env) %>%
  filter(study_bias != "No") %>%
  select(gene, Env, mean_abs_shap, study_bias, Consensus_Pres_Summary, breadth_present, breadth_absent, submitter_high, prop_high, submitter_second, prop_second, unevenness)

# ============================================================================
# Filter Genes that aren't "Broad" enough
# ============================================================================

# First, collect all genes across all niches
combined_data <- lapply(names(niche_results_filtered_clusters), function(nm) {
  niche_results_filtered_clusters[[nm]]$combined_data %>%
    mutate(niche = nm)
}) %>% bind_rows()

# Get total unique genes before filtering
total_genes_before <- combined_data %>% pull(gene) %>% unique() %>% length()

# Identify genes to remove (fail breadth in ANY niche, either direction)
genes_to_remove_breadth <- combined_data %>%
  filter(
    (Consensus_Pres_Summary == "Presence_Expected" & breadth_present < 0.25) |
    (Consensus_Pres_Summary == "Absence_Expected" & breadth_absent < 0.25)
  ) %>%
  pull(gene) %>%
  unique()

genes_to_keep <- combined_data %>%
  filter(!gene %in% genes_to_remove_breadth) %>%
  pull(gene) %>%
  unique()

cat(sprintf("\n=== Breadth Filtering Summary ===\n"))
cat(sprintf("Total genes before: %d\n", total_genes_before))
cat(sprintf("Genes removed: %d (%.1f%%)\n", 
            length(genes_to_remove_breadth), 
            100*length(genes_to_remove_breadth)/total_genes_before))
cat(sprintf("Genes kept: %d (%.1f%%)\n\n", 
            length(genes_to_keep),
            100*length(genes_to_keep)/total_genes_before))

# Apply filter to each niche
niche_results_filtered_breadth <- lapply(names(niche_results_filtered_clusters), function(nm) {
  nr <- niche_results_filtered_clusters[[nm]]
  n_before <- nrow(nr$combined_data)
  
  # Remove phylogenetically biased genes
  nr$combined_data <- nr$combined_data %>%
    filter(!gene %in% genes_to_remove_breadth)
  
  n_after <- nrow(nr$combined_data)
  n_removed <- n_before - n_after
  
  cat(sprintf("%s: removed %d genes (%.1f%%), %d remaining\n", 
              nm, n_removed, 100*n_removed/n_before, n_after))
  
  nr$n_genes <- n_after
  nr
})

names(niche_results_filtered_breadth) <- names(niche_results_filtered_clusters)

## ---- Sanity check after breadth filtering ----
niche_results_filtered_breadth$`Pediatric CF`$combined_data %>% arrange(desc(mean_abs_shap)) %>% head(10)

# ============================================================================
# Identify Which Genes Are Removed and Why (by niche)
# ============================================================================

flat_niche_results_filtered_clusters <- lapply(names(niche_results_filtered_clusters), function(nm) {
  niche_results_filtered_clusters[[nm]]$combined_data %>%
    mutate(Env = nm, .after = gene)
}) %>%
  bind_rows()

shap_importance_threshold <- 0.01  

genes_failing_breadth <- flat_niche_results_filtered_clusters %>%
  select(gene, Env, Consensus_Pres_Summary, mean_abs_shap, mean_shap_abs, mean_shap_pres, breadth_present, breadth_absent) %>%
  mutate(
    # Determine if fails WITHOUT threshold
    fails_without_threshold = (
      (Consensus_Pres_Summary == "Presence_Expected" & breadth_present < 0.25) |
      (Consensus_Pres_Summary == "Absence_Expected" & breadth_absent < 0.25)
    ),
    
    # Determine if fails WITH threshold
    fails_with_threshold = (
      mean_abs_shap >= shap_importance_threshold &  # Only if important
      (
        (Consensus_Pres_Summary == "Presence_Expected" & breadth_present < 0.25) |
        (Consensus_Pres_Summary == "Absence_Expected" & breadth_absent < 0.25)
      )
    ),
    
    # Add failure reason for those that fail
    failure_reason = case_when(
      Consensus_Pres_Summary == "Presence_Expected" & breadth_present < 0.25 ~ 
        sprintf("Presence_Expected but breadth_present=%.3f", breadth_present),
      Consensus_Pres_Summary == "Absence_Expected" & breadth_absent < 0.25 ~ 
        sprintf("Absence_Expected but breadth_absent=%.3f", breadth_absent),
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(gene, Env)

# ---- Summary statistics ----
genes_removed_no_threshold <- genes_failing_breadth %>% 
  filter(fails_without_threshold) %>% 
  pull(gene) %>% 
  unique()

genes_removed_with_threshold <- genes_failing_breadth %>% 
  filter(fails_with_threshold) %>% 
  pull(gene) %>% 
  unique()

genes_spared_by_threshold <- setdiff(genes_removed_no_threshold, genes_removed_with_threshold)

# ---- Look for some specific genes ----
# yecS, accC (group_254688), and mexA (group_349871) all showed up in our last analysis (without strict filtering), why are they lost now?
genes_to_check <- c("group_365760")

cat("\n=== Your Specific Genes of Interest ===\n")
cat("Which niches cause removal:\n")

specific_genes_comparison <- genes_failing_breadth %>%
  filter(gene %in% genes_to_check) %>%
  select(gene, Env, mean_shap_abs, Consensus_Pres_Summary, breadth_present, breadth_absent, 
         fails_without_threshold, fails_with_threshold, failure_reason) %>%
  arrange(gene, Env) %>%
  filter(fails_without_threshold == TRUE)

print(specific_genes_comparison)
print(specific_genes_removal)

if (nrow(specific_genes_removal) == 0) {
  cat("None of your genes of interest are being removed by breadth filtering!\n")
}

# ============================================================================
# Export ALL genes, with "removed for" categories
# ============================================================================

## ---- Track genes removed ----
bias_filtered_genes <- genes_to_remove_bias  # From expanded_results_updated$biased

breadth_failed_genes <- genes_to_remove_breadth

## ---- Create version with all genes + "Failed" column ----
niche_results_all_failed_passed <- lapply(names(niche_results), function(nm) {
  nr <- niche_results[[nm]]
  
  # Add "Failed" column based on pre-calculated gene lists
  nr$combined_data <- nr$combined_data %>%
    mutate(
      Failed = case_when(
        # If biased, mark as "Biased" (takes precedence)
        gene %in% bias_filtered_genes ~ "Biased",
        # If fails breadth
        gene %in% breadth_failed_genes ~ "Breadth",
        # Otherwise passes all filters
        TRUE ~ "No"
      )
    )
  
  nr
})

names(niche_results_all_failed_passed) <- names(niche_results)

niche_results_all_failed_passed_flat <- lapply(names(niche_results_all_failed_passed), function(nm) {
  niche_results_all_failed_passed[[nm]]$combined_data %>%
    mutate(Env = nm, .after = gene)
}) %>% bind_rows()

saveRDS(niche_results_all_failed_passed_flat, data_here("00_data/model_output/", "niche_results_all_failed_passed_flat.rds"))

# ============================================================================
# Helper Function to Save Niche Results Filtered
# ============================================================================

## ---- Save niche_results_filtered ----
save_niche_results <- function(niche_results_filtered,
                               out_dir,
                               rds_path,
                               slot = "combined_data",
                               file_prefix = "niche_",
                               compress_rds = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Save single RDS (entire object, for downstream reuse)
  saveRDS(niche_results_filtered, file = rds_path, compress = compress_rds)

  # Save per-niche CSVs (default: combined_data)
  niche_names <- names(niche_results_filtered)
  if (is.null(niche_names) || any(niche_names == "")) {
    niche_names <- as.character(seq_along(niche_results_filtered))
  }

  for (i in seq_along(niche_results_filtered)) {
    nm <- niche_names[i]
    x <- niche_results_filtered[[i]]

    df <- x[[slot]]
    if (is.null(df) || nrow(df) == 0) next
    
    # Arrange by mean_abs_shap descending
    if ("mean_abs_shap" %in% colnames(df)) {
      df <- df %>% dplyr::arrange(desc(mean_abs_shap))
    }

    # sanitize niche name for filename
    nm_file <- gsub("[^A-Za-z0-9._-]+", "_", nm)
    csv_path <- file.path(out_dir, paste0(file_prefix, nm_file, ".csv"))

    # fwrite is fast and handles large tables well
    data.table::fwrite(df, file = csv_path)
  }

  invisible(list(out_dir = out_dir, rds_path = rds_path))
}

# ============================================================================
# Save Processed Niche Results Filtered
# ============================================================================

out_dir <- "00_data/model_output/Filtered_Features"
rds_path <- "00_data/model_output/niche_results.rds"

save_niche_results(
  niche_results_filtered = niche_results_filtered_breadth,
  out_dir = out_dir,
  rds_path = rds_path,
  slot = "combined_data"
)
