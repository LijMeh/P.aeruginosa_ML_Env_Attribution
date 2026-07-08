# Longer Functions used to process data for figures

# ---- Functions for Blocked Model Evaluation ----

parse_blocked_accuracy <- function(df) {
  df %>%
    mutate(
      Cluster = str_extract(model, "^Cluster_\\d+"),
      Model.Type = case_when(
        str_detect(model, "_isolated$") ~ "Blocked",
        str_detect(model, "_baseline$") ~ "Full",
        TRUE ~ NA_character_
      ),
      Number.Isolates = str_count(`True.Labels`, "'") / 2
    ) %>%
    select(Cluster, Model.Type, `Micro.F1`, `Macro.F1`, Accuracy, `Balanced.Accuracy`, Number.Isolates)
}

calc_env_accuracy <- function(df, blocked = TRUE) {
  # Filter for isolated models
  if(blocked == TRUE){
    df <- df %>% filter(str_detect(model, "_isolated$"))
    } else {
        df <- df %>% filter(str_detect(model, "_baseline$"))
    }
  
  # For each row, split True and Predicted labels, then compare
  env_results <- lapply(seq_len(nrow(df)), function(i) {
    true_envs <- str_extract_all(df$True.Labels[i], "'([^']+)'")[[1]]
    pred_envs <- str_extract_all(df$Predicted.Labels[i], "'([^']+)'")[[1]]
    # Remove quotes
    true_envs <- gsub("'", "", true_envs)
    pred_envs <- gsub("'", "", pred_envs)
    data.frame(
      ENV = true_envs,
      Correct = true_envs == pred_envs,
      stringsAsFactors = FALSE
    )
  })

  env_results_df <- bind_rows(env_results)

  env_summary <- env_results_df %>%
    group_by(ENV) %>%
    summarise(
      n = n(),
      correct = sum(Correct),
      wrong = n - correct,
      accuracy = correct / n
    ) %>%
    arrange(desc(accuracy))
  
  env_summary
}

calc_normalized_accuracy <- function(env_df) {
  num_unique_envs <- env_df %>% distinct(ENV) %>% nrow()
  prob_random <- 100 / num_unique_envs
  
  env_df %>%
    mutate(
      normalized_accuracy = ((accuracy * 100 - prob_random) / (100 - prob_random))
    )
}

calc_overall_normalized_accuracy <- function(env_df) {
  overall_accuracy <- mean(env_df$accuracy) * 100
  num_unique_envs <- nrow(env_df)
  prob_random <- 100 / num_unique_envs
  ((overall_accuracy - prob_random) / (100 - prob_random))
}

# Alternative: return as wide format without pivot_wider
get_misclassifications <- function(df, blocked = TRUE) {
  if(blocked == TRUE){
    df <- df %>% filter(str_detect(model, "_isolated$"))
  } else {
    df <- df %>% filter(str_detect(model, "_baseline$"))
  }
  
  misclass_results <- lapply(seq_len(nrow(df)), function(i) {
    true_envs <- str_extract_all(df$True.Labels[i], "'([^']+)'")[[1]]
    pred_envs <- str_extract_all(df$Predicted.Labels[i], "'([^']+)'")[[1]]
    true_envs <- gsub("'", "", true_envs)
    pred_envs <- gsub("'", "", pred_envs)
    incorrect_mask <- true_envs != pred_envs
    data.frame(
      True_ENV = true_envs[incorrect_mask],
      Predicted_ENV = pred_envs[incorrect_mask],
      stringsAsFactors = FALSE
    )
  })
  
  misclass_df <- bind_rows(misclass_results)
  
  misclass_ranked <- misclass_df %>%
    group_by(True_ENV, Predicted_ENV) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(True_ENV) %>%
    mutate(
      total_misclass = sum(count),
      relative_frequency = count / total_misclass,
      rank = rank(-count)
    ) %>%
    filter(rank <= 3) %>%
    arrange(True_ENV, rank)
  
  return(misclass_ranked)
}

# Function to plot misclassifications
plot_misclassifications <- function(misclass_df, round_num) {
  misclass_df %>%
    ggplot(aes(x = reorder(True_ENV, -rank), y = relative_frequency, fill = Predicted_ENV)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = paste0(count)), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 3) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = paste0("Top Misclassifications by Environment (Round ", round_num, ")"),
      x = "True Environment",
      y = "Relative Frequency of Misclassification",
      fill = "Misclassified as",
      caption = "Labels show count; bars ordered by frequency (tallest = most common misclassification)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 13, face = "bold"),
      legend.position = "right"
    )
}

# ---- Calculate Genomovar Breadth for Each Gene-Niche Combination ----
calculate_genomovar_breadth_with_importance <- function(shap_data, metadata, model_input_filtered) {
  # Get unique niches from shap_data
  niches <- unique(shap_data$Niche)
  
  # Initialize output list
  breadth_results <- list()
  
  for (niche in niches) {
    message("Processing niche: ", niche)
    
    # Step 1: Get top genes for this niche
    genes_in_niche <- unique(shap_data$gene[shap_data$Niche == niche])
    
    # Step 2: Get genomes in this niche from metadata 
    genomes_in_niche <- metadata$Genome_ID[metadata$Niche == niche]
    
    # Step 3: Get genomovar clusters for these genomes
    metadata_subset <- metadata %>%
      filter(Niche == niche) %>%
      select(Genome_ID, Genomovar_Cluster)
    
    # Step 4: Get gene presence/absence for genomes in this niche
    model_data_subset <- model_input_filtered %>%
      filter(Genome_ID %in% genomes_in_niche) %>%
      left_join(metadata_subset, by = "Genome_ID")
    
    # Get all unique genomovar clusters for this niche
    unique_genomovars <- unique(model_data_subset$Genomovar_Cluster)
    
    # Step 5: For each top gene, calculate breadth
    for (gene in genes_in_niche) {
      if (gene %in% colnames(model_data_subset)) {
        # Get presence/absence data for this gene
        gene_data <- model_data_subset %>%
          select(Genome_ID, Genomovar_Cluster, all_of(gene)) %>%
          rename(Gene_Present = all_of(gene)) %>%
          mutate(Gene_Present = as.numeric(Gene_Present))
        
        # Count present/absent per genomovar cluster
        genomovar_stats <- gene_data %>%
          group_by(Genomovar_Cluster) %>%
          summarise(
            Present_Count = sum(Gene_Present, na.rm = TRUE),
            Absent_Count = sum(1 - Gene_Present, na.rm = TRUE),
            Total_Count = n(),
            .groups = "drop"
          )
        
        # Breadth: proportion of genomovars with at least one presence
        breadth <- sum(genomovar_stats$Present_Count > 0) / length(unique_genomovars)
        
        # Get SHAP importance from shap_data
        gene_info <- shap_data %>%
          filter(gene == !!gene & Niche == niche)
        
        result <- data.frame(
          Gene = gene,
          Niche = niche,
          Breadth = breadth,
          mean_shap_pres = gene_info$mean_shap_pres,
          mean_shap_abs = gene_info$mean_shap_abs,
          mean_abs_shap = gene_info$mean_abs_shap
        )
        
        breadth_results[[paste(niche, gene, sep = "_")]] <- result
      }
    }
  }
  
  # Combine results
  bind_rows(breadth_results) %>%
    select(Gene, Niche, Breadth, mean_shap_pres, mean_shap_abs, mean_abs_shap) %>%
    arrange(Niche, Gene)
}


# ---- Categorize Gene Annotations into Two-Tier Functional Categories ----
categorize_two_tier_per_part <- function(anno) {
  # Handle blank annotations
  if (is.na(anno) || trimws(anno) == "") {
    return(data.frame(Level1 = "Uncharacterized",
                     Level2 = "Unknown Functions",
                     Weight = 1))
  }
  
  parts <- strsplit(anno, ";")[[1]] %>% trimws()
  
  matched_list <- lapply(parts, function(p) {
    if (grepl("^hypothetical protein$", p, ignore.case = TRUE)) {
      return(data.frame(Level1 = "Uncharacterized",
                       Level2 = "Unknown Functions"))
    }
    for (i in seq_len(nrow(category_map_two_tier))) {
      if (grepl(category_map_two_tier$Pattern[i], p, ignore.case = TRUE)) {
        return(data.frame(Level1 = category_map_two_tier$Level1[i],
                         Level2 = category_map_two_tier$Level2[i]))
      }
    }
    return(NULL)
  })
  
  matched_list <- Filter(Negate(is.null), matched_list)
  
  if (length(matched_list) == 0) {
    return(data.frame(Level1 = "Other",
                     Level2 = "Other",
                     Weight = 1))
  }
  
  matched_df <- bind_rows(matched_list)
  real_rows <- matched_df$Level2 != "Unknown Functions" & !is.na(matched_df$Level2)
  
  if (any(real_rows)) {
    matched_df <- matched_df[real_rows, , drop = FALSE]
  }
  
  if (nrow(matched_df) == 0) {
    return(data.frame(Level1 = "Uncharacterized",
                     Level2 = "Unknown Functions",
                     Weight = 1))
  }
  
  matched_df <- distinct(matched_df)
  matched_df$Weight <- 1 / nrow(matched_df)
  
  matched_df
}

parse_blocked_accuracy <- function(df) {
  df %>%
    mutate(
      Cluster = str_extract(model, "^Cluster_\\d+"),
      Model.Type = case_when(
        str_detect(model, "_isolated$") ~ "Blocked",
        str_detect(model, "_baseline$") ~ "Full",
        TRUE ~ NA_character_
      ),
      Number.Isolates = str_count(`True.Labels`, "'") / 2
    ) %>%
    select(Cluster, Model.Type, `Micro.F1`, `Macro.F1`, Accuracy, `Balanced.Accuracy`, Number.Isolates)
}

calc_env_accuracy <- function(df, blocked = TRUE) {
  # Filter for isolated models
  if(blocked == TRUE){
    df <- df %>% filter(str_detect(model, "_isolated$"))
    } else {
        df <- df %>% filter(str_detect(model, "_baseline$"))
    }
  
  # For each row, split True and Predicted labels, then compare
  env_results <- lapply(seq_len(nrow(df)), function(i) {
    true_envs <- str_extract_all(df$True.Labels[i], "'([^']+)'")[[1]]
    pred_envs <- str_extract_all(df$Predicted.Labels[i], "'([^']+)'")[[1]]
    # Remove quotes
    true_envs <- gsub("'", "", true_envs)
    pred_envs <- gsub("'", "", pred_envs)
    data.frame(
      ENV = true_envs,
      Correct = true_envs == pred_envs,
      stringsAsFactors = FALSE
    )
  })

  env_results_df <- bind_rows(env_results)

  env_summary <- env_results_df %>%
    group_by(ENV) %>%
    summarise(
      n = n(),
      correct = sum(Correct),
      wrong = n - correct,
      accuracy = correct / n
    ) %>%
    arrange(desc(accuracy))
  
  env_summary
}

calc_normalized_accuracy <- function(env_df) {
  num_unique_envs <- env_df %>% distinct(ENV) %>% nrow()
  prob_random <- 100 / num_unique_envs
  
  env_df %>%
    mutate(
      normalized_accuracy = ((accuracy * 100 - prob_random) / (100 - prob_random))
    )
}

calc_overall_normalized_accuracy <- function(env_df) {
  overall_accuracy <- mean(env_df$accuracy) * 100
  num_unique_envs <- nrow(env_df)
  prob_random <- 100 / num_unique_envs
  ((overall_accuracy - prob_random) / (100 - prob_random))
}

# Alternative: return as wide format without pivot_wider
get_misclassifications <- function(df, blocked = TRUE) {
  if(blocked == TRUE){
    df <- df %>% filter(str_detect(model, "_isolated$"))
  } else {
    df <- df %>% filter(str_detect(model, "_baseline$"))
  }
  
  misclass_results <- lapply(seq_len(nrow(df)), function(i) {
    true_envs <- str_extract_all(df$True.Labels[i], "'([^']+)'")[[1]]
    pred_envs <- str_extract_all(df$Predicted.Labels[i], "'([^']+)'")[[1]]
    true_envs <- gsub("'", "", true_envs)
    pred_envs <- gsub("'", "", pred_envs)
    incorrect_mask <- true_envs != pred_envs
    data.frame(
      True_ENV = true_envs[incorrect_mask],
      Predicted_ENV = pred_envs[incorrect_mask],
      stringsAsFactors = FALSE
    )
  })
  
  misclass_df <- bind_rows(misclass_results)
  
  misclass_ranked <- misclass_df %>%
    group_by(True_ENV, Predicted_ENV) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(True_ENV) %>%
    mutate(
      total_misclass = sum(count),
      relative_frequency = count / total_misclass,
      rank = rank(-count)
    ) %>%
    filter(rank <= 3) %>%
    arrange(True_ENV, rank)
  
  return(misclass_ranked)
}

# Function to plot misclassifications
plot_misclassifications <- function(misclass_df, round_num) {
  misclass_df %>%
    ggplot(aes(x = reorder(True_ENV, -rank), y = relative_frequency, fill = Predicted_ENV)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = paste0(count)), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 3) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = paste0("Top Misclassifications by Environment (Round ", round_num, ")"),
      x = "True Environment",
      y = "Relative Frequency of Misclassification",
      fill = "Misclassified as",
      caption = "Labels show count; bars ordered by frequency (tallest = most common misclassification)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 13, face = "bold"),
      legend.position = "right"
    )
}

# ---- Process Confusion Matrix Data ----
process_confusion_matrix <- function(pairs_df) {
  # Core processing logic shared by both functions
  
  # Get actual counts FIRST (before setting factor levels)
  cm_long <- pairs_df %>%
    dplyr::count(True, Pred, name = "n")
  
  # Calculate row percentages from actual counts
  cm_row_pct <- cm_long %>%
    dplyr::group_by(True) %>%
    dplyr::mutate(
      row_total = sum(n),
      row_pct = n / row_total
    ) %>%
    dplyr::ungroup()
  
  # Set factor levels to include all possible combinations for complete matrix
  all_levels <- sort(unique(c(pairs_df$True, pairs_df$Pred)))
  
  # Complete the long format with zeros for missing combinations
  cm_long_complete <- cm_long %>%
    dplyr::mutate(
      True = factor(True, levels = all_levels),
      Pred = factor(Pred, levels = all_levels)
    ) %>%
    tidyr::complete(True, Pred, fill = list(n = 0))
  
  # Complete the row_pct with zeros for missing combinations
  cm_row_pct_complete <- cm_row_pct %>%
    dplyr::mutate(
      True = factor(True, levels = all_levels),
      Pred = factor(Pred, levels = all_levels)
    ) %>%
    tidyr::complete(True, Pred, fill = list(n = 0, row_total = 0, row_pct = 0))
  
  # Confusion matrix (wide format)
  cm_counts <- cm_long_complete %>%
    tidyr::pivot_wider(names_from = Pred, values_from = n, values_fill = 0)
  
  # Update pairs_df with factor levels
  pairs_df <- pairs_df %>%
    dplyr::mutate(
      True = factor(True, levels = all_levels),
      Pred = factor(Pred, levels = all_levels)
    )
  
  list(
    pairs = pairs_df,
    confusion_counts_wide = cm_counts,
    confusion_counts_long = cm_long_complete,
    confusion_row_pct = cm_row_pct_complete
  )
}


# ---- Get Niche assingment for specific clusters ----
get_niche_assignments_for_cluster <- function(df_blocked, cluster_name) {
  # Find the row for the isolated test of the given cluster
  pattern <- paste0("^", cluster_name, "_stats_isolated$")
  row_idx <- which(stringr::str_detect(df_blocked$model, pattern))
  if (length(row_idx) != 1) {
    stop("Could not uniquely identify the isolated row for this cluster.")
  }
  
  # Extract true and predicted labels
  true_envs <- stringr::str_extract_all(df_blocked$True.Labels[row_idx], "'([^']+)'")[[1]] %>% gsub("'", "", .)
  pred_envs <- stringr::str_extract_all(df_blocked$Predicted.Labels[row_idx], "'([^']+)'")[[1]] %>% gsub("'", "", .)
  
  # Safety check
  if (length(true_envs) != length(pred_envs)) {
    stop("Mismatch in number of true and predicted labels for this cluster.")
  }
  
  tibble::tibble(
    True_Niche = true_envs,
    Pred_Niche = pred_envs
  )
}

# ---- Calculate Accuracy For Confusion Matrix ----
calculate_accuracy_metrics <- function(cm_result, env_labels) {
  # Convert to wide matrix
  confusion_matrix <- cm_result$confusion_counts_long %>%
    filter(True %in% env_labels, Pred %in% env_labels) %>%
    mutate(
      True = factor(True, levels = env_labels),
      Pred = factor(Pred, levels = env_labels)
    ) %>%
    pivot_wider(names_from = Pred, values_from = n, values_fill = 0) %>%
    column_to_rownames("True") %>%
    as.matrix()
  
  # Calculate metrics
  actual_accuracy <- (diag(confusion_matrix) / rowSums(confusion_matrix) * 100)[env_labels]
  pred_totals <- colSums(confusion_matrix)
  null_accuracy <- (pred_totals / sum(pred_totals) * 100)[env_labels]
  
  # Return data frame
  data.frame(
    category = factor(env_labels, levels = env_labels),
    actual = actual_accuracy,
    null = null_accuracy
  ) %>%
    pivot_longer(
      cols = c(actual, null),
      names_to = "type",
      values_to = "accuracy"
    )
}

# ---- Calculate SHAP Variability Across ENV Groups ----
calculate_variability_plot_data <- function(
  df,
  niche_groups,
  min_niches_per_group = 2,
  require_consensus = "all"
) {
  
  # Helper function to assign group
  assign_group <- function(env, niche_groups) {
    for (group in names(niche_groups)) {
      if (env %in% niche_groups[[group]]) return(group)
    }
    return(NA_character_)
  }
  
  # Filter genes meeting criteria
  genes_filtered <- df %>%
    mutate(Group = sapply(Env, function(e) assign_group(e, niche_groups))) %>%
    filter(!is.na(Group)) %>%
    group_by(gene, Group, Annotation_Inferred) %>%
    summarise(
      niche_count = n_distinct(Env),
      n_consensus_states = n_distinct(Consensus_Pres_Summary, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(niche_count >= min_niches_per_group)
  
  # Apply consensus filtering
  if (require_consensus %in% c("min", "all")) {
    genes_filtered <- genes_filtered %>%
      filter(n_consensus_states == 1)
  }
  
  genes_filtered <- genes_filtered %>%
    select(gene, Group, Annotation_Inferred)
  
  # Calculate both metrics in one pipeline
  df %>%
    mutate(Group = sapply(Env, function(e) assign_group(e, niche_groups))) %>%
    filter(!is.na(Group)) %>%
    semi_join(genes_filtered, by = c("gene", "Group", "Annotation_Inferred")) %>%
    group_by(gene, Group, Annotation_Inferred) %>%
    summarise(
      variability = sd(Net_Predictive_Value, na.rm = TRUE),
      mean_SHAP_group = mean(Net_Predictive_Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Group,
      values_from = c(variability, mean_SHAP_group),
      names_sep = "_"
    ) %>%
    mutate(
      NPV_diff_Chronic = mean_SHAP_group_Chronic - mean_SHAP_group_Acute,
      NPV_diff_Acute = mean_SHAP_group_Acute - mean_SHAP_group_Chronic
    ) %>%
    pivot_longer(
      cols = starts_with("NPV_diff_"),
      names_to = "group_perspective",
      values_to = "NPV_difference"
    ) %>%
    mutate(
      group_perspective = str_remove(group_perspective, "NPV_diff_"),
      variability_focal = case_when(
        group_perspective == "Chronic" ~ variability_Chronic,
        group_perspective == "Acute" ~ variability_Acute
      )
    ) %>%
    select(gene, Annotation_Inferred, Group = group_perspective, variability_focal, NPV_difference)
}

# ---- Remap Niche Names Using Standard Mapping ----
remap_niches <- function(niche_vector) {
  niche_map <- c(
    "Bronchiectasis" = "Bronchiectasis",
    "CF" = "Adult CF",
    "early.CF" = "Pediatric CF",
    "Early CF" = "Pediatric CF",
    "Early.CF" = "Pediatric CF",
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
    "Rectal-Feces" = "Fecal",
    "Rectal.Feces" = "Fecal"
  )
  recode(niche_vector, !!!niche_map, .default = NA_character_)
}

# ---- Create Histogram Wave Visualization from NPV Data ----
make_hist_wave_npv_units <- function(x, y_npv, bins = 80, smooth_window = 5) {
  ok <- is.finite(x) & is.finite(y_npv)
  x <- pmax(x[ok], 0)
  y_npv <- y_npv[ok]

  if (length(x) < 2) {
    return(tibble::tibble(x = numeric(0), y = numeric(0), type = character(0)))
  }

  # Start explicitly at 0 for gradual buildup
  brks <- seq(0, max(x, na.rm = TRUE), length.out = bins + 1)
  mids <- (brks[-1] + brks[-length(brks)]) / 2
  bin_id <- cut(x, breaks = brks, include.lowest = TRUE, labels = FALSE)

  binned <- tibble::tibble(bin = bin_id, y_npv = y_npv) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      npv_pos = sum(pmax(y_npv, 0), na.rm = TRUE),
      npv_neg = sum(pmin(y_npv, 0), na.rm = TRUE),
      .groups = "drop"
    )

  all_bins <- tibble::tibble(bin = seq_len(bins), x = mids) %>%
    dplyr::left_join(binned, by = "bin") %>%
    dplyr::mutate(
      npv_pos = dplyr::coalesce(npv_pos, 0),
      npv_neg = dplyr::coalesce(npv_neg, 0)
    )

  # light smoothing (moving average)
  smooth_ma <- function(v, k = 5) {
    if (k <= 1) return(v)
    s <- as.numeric(stats::filter(v, rep(1 / k, k), sides = 2))
    s[is.na(s)] <- v[is.na(s)]  # keep edges stable
    s
  }

  pos_s <- smooth_ma(all_bins$npv_pos, smooth_window)
  neg_s <- -smooth_ma(abs(all_bins$npv_neg), smooth_window)

  dplyr::bind_rows(
    tibble::tibble(x = 0, y = 0, type = "npv_pos"),
    tibble::tibble(x = all_bins$x, y = pos_s, type = "npv_pos"),
    tibble::tibble(x = 0, y = 0, type = "npv_neg"),
    tibble::tibble(x = all_bins$x, y = neg_s, type = "npv_neg")
  )
}

# ---- Load All Gene CSV Files from Subdirectories ----
load_all_gene_files <- function(base_path) {
  # Get all subdirectories (gene folders)
  gene_dirs <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)
  
  # Load all CSV files matching the pattern
  map_df(gene_dirs, ~{
    files <- list.files(.x, pattern = "_funct_no_funct_alt\\.csv$", full.names = TRUE)
    if (length(files) > 0) {
      read.csv(files[1]) %>% 
        mutate(gene = basename(.x), .before = 1)
    }
  })
}

# ---- Fit GLMER Model for Single Gene ----
fit_gene_glmer <- function(gene_name, df_analysis, random_effect_vars = NULL) {
  
  gene_name <- as.character(gene_name)
  random_effect_vars <- random_effect_vars %>%
    unlist(use.names = FALSE) %>%
    as.character()
  random_effect_vars <- unique(random_effect_vars[!is.na(random_effect_vars) & nzchar(trimws(random_effect_vars))])
  use_random_effect <- length(random_effect_vars) > 0
  
  df_gene <- df_analysis %>%
    filter(gene == gene_name) %>%
    droplevels()
  
  cat(strrep("=", 80), "\n")
  cat("ANALYZING GENE:", gene_name, "\n")
  cat(strrep("=", 80), "\n")
  cat("Samples: n =", nrow(df_gene), "\n")
  cat("Mutation rate:", round(mean(df_gene$Mutated), 3), "\n\n")
    
  # Skip if too few samples
  if (nrow(df_gene) < 10) {
    cat("  ✗ Skipped: too few samples\n\n")
    return(list(
      gene = gene_name,
      status = "skipped",
      reason = "too few samples"
    ))
  }
  
  # Check if we have variation in both outcome and predictors
  if (length(unique(df_gene$Mutated)) < 2) {
    cat("  ✗ Skipped: no variation in mutation status\n\n")
    return(list(
      gene = gene_name,
      status = "skipped",
      reason = "no variation in outcome"
    ))
  }
  
  if (dplyr::n_distinct(stats::na.omit(df_gene$Niche)) < 2) {
    cat("  ✗ Skipped: only one environment represented\n\n")
    return(list(
      gene = gene_name,
      status = "skipped",
      reason = "insufficient environment variation"
    ))
  }

  if (use_random_effect && !all(random_effect_vars %in% names(df_gene))) {
    missing_vars <- setdiff(random_effect_vars, names(df_gene))
    cat("  ✗ Skipped: random effect column not found in analysis data\n\n")
    return(list(
      gene = gene_name,
      status = "skipped",
      reason = paste0("random effect column missing: ", paste(missing_vars, collapse = ", "))
    ))
  }

  if (use_random_effect) {
    low_var_cols <- random_effect_vars[vapply(
      random_effect_vars,
      function(v) dplyr::n_distinct(stats::na.omit(df_gene[[v]])) < 2,
      logical(1)
    )]
  } else {
    low_var_cols <- character(0)
  }

  if (length(low_var_cols) > 0) {
    cat("  ✗ Skipped: one or more random-effect groups have only one observed level\n\n")
    return(list(
      gene = gene_name,
      status = "skipped",
      reason = paste0("insufficient random-effect variation: ", paste(low_var_cols, collapse = ", "))
    ))
  }
  
  tryCatch({
    if (use_random_effect) {
      random_terms <- paste0("(1 | `", random_effect_vars, "`)", collapse = " + ")
      model_formula <- as.formula(paste0("Mutated ~ Niche + ", random_terms))

      model <- lme4::glmer(
        model_formula,
        data = df_gene,
        family = binomial(link = "logit"),
        control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
        nAGQ = 1
      )
      random_effects <- lme4::ranef(model)
      fixed_effects <- lme4::fixef(model)
    } else {
      model_formula <- Mutated ~ Niche
      model <- glm(
        model_formula,
        data = df_gene,
        family = binomial(link = "logit")
      )
      random_effects <- NULL
      fixed_effects <- coef(model)
    }
    
    cat("  ✓ Model fitted successfully\n\n")
    
    list(
      gene = gene_name,
      status = "success",
      model = model,
      coef_summary = coef(summary(model)),
      fixed_effects = fixed_effects,
      random_effects = random_effects,
      vcov = vcov(model),
      n_obs = nrow(df_gene),
      n_mutated = sum(df_gene$Mutated)
    )
    
  }, error = function(e) {
    cat("  ✗ Model fitting failed:", e$message, "\n\n")
    list(
      gene = gene_name,
      status = "error",
      error_message = e$message
    )
  })
}

# ---- Summarize Effects from Fitted GLMER Model ----
summarize_gene_effects <- function(result_obj, df_analysis, gene_name_lookup) {
  if (is.null(result_obj) || result_obj$status != "success") {
    return(NULL)
  }

  gene_label <- unname(gene_name_lookup[as.character(result_obj$gene)])
  if (length(gene_label) == 0 || is.na(gene_label) || gene_label == "") {
    gene_label <- as.character(result_obj$gene)
  }

  emm_obj <- tryCatch(
    emmeans::emmeans(result_obj$model, ~ Niche),
    error = function(e) NULL
  )
  if (is.null(emm_obj)) {
    return(NULL)
  }

  emm_probs <- tryCatch(
    as.data.frame(summary(emm_obj, type = "response")),
    error = function(e) NULL
  )
  if (is.null(emm_probs) || nrow(emm_probs) < 2) {
    return(NULL)
  }

  prob_col <- dplyr::case_when(
    "prob" %in% names(emm_probs) ~ "prob",
    "response" %in% names(emm_probs) ~ "response",
    "emmean" %in% names(emm_probs) ~ "emmean",
    TRUE ~ NA_character_
  )
  if (is.na(prob_col)) {
    return(NULL)
  }

  prob_tbl <- emm_probs %>%
    transmute(
      Niche = as.character(Niche),
      Probability = as.numeric(.data[[prob_col]])
    )

  pair_tbl <- tryCatch(
    as.data.frame(
      summary(
        emmeans::contrast(emm_obj, method = "pairwise", adjust = "tukey"),
        infer = c(TRUE, TRUE)
      )
    ),
    error = function(e) NULL
  )
  if (is.null(pair_tbl) || nrow(pair_tbl) == 0 || !"contrast" %in% names(pair_tbl)) {
    return(NULL)
  }

  pair_tbl <- pair_tbl %>%
    mutate(
      Comparison_Group = trimws(sub(" - .*", "", contrast)),
      Reference_Group = trimws(sub(".* - ", "", contrast))
    ) %>%
    left_join(
      prob_tbl %>% rename(Comparison_Group = Niche, Comparison_Probability = Probability),
      by = "Comparison_Group"
    ) %>%
    left_join(
      prob_tbl %>% rename(Reference_Group = Niche, Reference_Probability = Probability),
      by = "Reference_Group"
    )

  pair_tbl %>%
    mutate(
      Difference_Percentage_Points = (Comparison_Probability - Reference_Probability) * 100,
      P_Value = .data[["p.value"]],
      Significant = ifelse(!is.na(P_Value) & P_Value < 0.05, "Yes", "No"),
      Direction = ifelse(
        Difference_Percentage_Points > 0,
        "Higher",
        ifelse(Difference_Percentage_Points < 0, "Lower", "No change")
      ),
      Comparison = paste0(Comparison_Group, " vs ", Reference_Group),
      Interpretation = dplyr::case_when(
        Significant == "No" ~ paste0("No significant difference detected for ", Comparison_Group, " vs ", Reference_Group),
        Difference_Percentage_Points > 0 ~ paste0(Comparison_Group, " has a higher mutation probability than ", Reference_Group, " by ", round(abs(Difference_Percentage_Points), 1), " percentage points"),
        Difference_Percentage_Points < 0 ~ paste0(Comparison_Group, " has a lower mutation probability than ", Reference_Group, " by ", round(abs(Difference_Percentage_Points), 1), " percentage points"),
        TRUE ~ paste0(Comparison_Group, " and ", Reference_Group, " have similar mutation probability")
      ),
      Gene = result_obj$gene,
      Gene_Name = gene_label
    ) %>%
    transmute(
      Gene,
      Gene_Name,
      Comparison,
      Reference_Group,
      Comparison_Group,
      Reference_Probability = round(Reference_Probability, 3),
      Comparison_Probability = round(Comparison_Probability, 3),
      Difference_Percentage_Points = round(Difference_Percentage_Points, 1),
      Direction,
      Significant,
      P_Value,
      Interpretation
    )
}

# ---- Run Full GLMER Analysis Pipeline ----
run_glmer_analysis <- function(
  df_input,
  reference_name,
  niche_ref = NULL,
  random_effect_var = NULL,
  random_effect_vars = NULL,
  mut_wt_values = "Functional",
  mut_case_values = NULL,
  genes_to_analyze = NULL
) {
  random_effect_vars <- c(random_effect_vars, random_effect_var) %>%
    unlist(use.names = FALSE) %>%
    as.character()
  random_effect_vars <- unique(random_effect_vars[!is.na(random_effect_vars) & nzchar(trimws(random_effect_vars))])
  use_random_effect <- length(random_effect_vars) > 0
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("ANALYZING:", reference_name, "\n")
  cat(strrep("=", 80), "\n\n")
  
  # USER INPUT: Select genes to analyze (NULL = all genes)
  if (is.null(genes_to_analyze)) {
    genes_to_analyze <- df_input$gene %>% unique() %>% as.character()
  }

  # Validate requested reference level (optional)
  if (!is.null(niche_ref) && !(niche_ref %in% unique(df_input$Niche))) {
    stop(paste0("Requested niche_ref '", niche_ref, "' not found in df_input$Niche"))
  }

  # Validate requested random effect columns
  if (use_random_effect && !all(random_effect_vars %in% names(df_input))) {
    missing_vars <- setdiff(random_effect_vars, names(df_input))
    stop(paste0("Requested random_effect_vars not found in df_input: ", paste(missing_vars, collapse = ", ")))
  }

  # Validate mutation definitions
  if (!is.null(mut_case_values) && any(mut_case_values %in% mut_wt_values)) {
    stop("mut_case_values and mut_wt_values overlap. They must define disjoint groups.")
  }

  cat("Model settings:\n")
  cat("  Niche reference level:", ifelse(is.null(niche_ref), "none (existing factor order)", niche_ref), "\n")
  cat("  Random effect variable(s):", if (use_random_effect) paste(random_effect_vars, collapse = ", ") else "none", "\n")
  if (is.null(mut_case_values)) {
    cat("  Mutation definition: Mutated = Mut_Status NOT IN", paste(mut_wt_values, collapse = ", "), "\n\n")
  } else {
    cat("  Mutation definition: Mutated = Mut_Status IN", paste(mut_case_values, collapse = ", "), "\n")
    cat("                     Wild-Type = Mut_Status IN", paste(mut_wt_values, collapse = ", "), "\n")
    cat("                     Other statuses excluded from analysis\n\n")
  }
  
  niche_factor <- as.factor(df_input$Niche)
  if (!is.null(niche_ref)) {
    niche_factor <- relevel(niche_factor, ref = niche_ref)
  }

  # Prepare data with user-defined binary mutation status
  df_analysis <- df_input %>%
    mutate(
      Mut_Status_chr = as.character(Mut_Status),
      Mutated = case_when(
        Mut_Status_chr %in% mut_wt_values ~ 0,
        is.null(mut_case_values) & !is.na(Mut_Status_chr) ~ 1,
        !is.null(mut_case_values) & Mut_Status_chr %in% mut_case_values ~ 1,
        TRUE ~ NA_real_
      ),
      Niche = niche_factor,
      gene = as.factor(gene)
    )

  if (use_random_effect) {
    df_analysis <- df_analysis %>%
      mutate(across(all_of(random_effect_vars), as.factor)) %>%
      drop_na(Mutated, Niche, all_of(random_effect_vars))
  } else {
    df_analysis <- df_analysis %>%
      drop_na(Mutated, Niche)
  }

  df_analysis <- df_analysis %>%
    filter(gene %in% genes_to_analyze) %>%
    dplyr::select(-Mut_Status_chr)
  
  cat("Overall data summary:\n")
  cat("  Total observations:", nrow(df_analysis), "\n")
  cat("  Genes:", nlevels(df_analysis$gene), "\n")
  cat("  Environments:", nlevels(df_analysis$Niche), "\n")
  if (use_random_effect) {
    random_group_summary <- paste(
      paste0(random_effect_vars, "=", vapply(random_effect_vars, function(v) nlevels(df_analysis[[v]]), numeric(1))),
      collapse = ", "
    )
    cat("  Random effect groups:", random_group_summary, "\n")
  } else {
    cat("  Random effect groups: not used\n")
  }
  cat("  Overall mutation rate:", round(mean(df_analysis$Mutated), 3), "\n\n")
  
  # Fit models for selected genes
  results_list <- map(
    unique(df_analysis$gene),
    ~ fit_gene_glmer(.x, df_analysis, random_effect_vars = if (use_random_effect) random_effect_vars else NULL)
  )
  
  names(results_list) <- sapply(results_list, function(x) x$gene)
  
  # Print status summary
  status_summary <- table(sapply(results_list, function(x) x$status))
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("Model Fitting Summary:\n")
  cat(strrep("=", 80), "\n")
  print(status_summary)
  cat("\n")
  
  # Build a gene-name lookup from the analysis input
  gene_name_lookup <- df_input %>%
    mutate(
      gene = as.character(gene),
      Gene.Name = as.character(Gene.Name)
    ) %>%
    dplyr::select(gene, Gene.Name) %>%
    distinct() %>%
    mutate(Gene.Name = ifelse(is.na(Gene.Name) | Gene.Name == "", gene, Gene.Name)) %>%
    deframe()
  
  coef_list <- map(results_list, ~{
    if (.x$status == "success") {
      .x$coef_summary
    } else {
      NULL
    }
  })
  
  coef_list <- map(coef_list, ~{
    if (!is.null(.x)) {
      .x <- as.data.frame(.x)
      .x$Significant <- ifelse(.x$`Pr(>|z|)` < 0.05, "Yes", "No")
      .x$Odds_Ratio <- exp(.x$Estimate)
      .x$OR_CI_Lower <- exp(.x$Estimate - 1.96 * .x$`Std. Error`)
      .x$OR_CI_Upper <- exp(.x$Estimate + 1.96 * .x$`Std. Error`)
      .x$term <- rownames(.x)
      .x
    } else {
      NULL
    }
  })

  coef_list_non_null <- coef_list[!vapply(coef_list, is.null, logical(1))]

  if (length(coef_list_non_null) == 0) {
    cat("No models fit successfully. Returning empty result tables.\n\n")

    coef_flat <- tibble(
      Gene = character(),
      Estimate = numeric(),
      `Std. Error` = numeric(),
      `z value` = numeric(),
      `Pr(>|z|)` = numeric(),
      Significant = character(),
      Odds_Ratio = numeric(),
      OR_CI_Lower = numeric(),
      OR_CI_Upper = numeric(),
      term = character(),
      Gene_Name = character(),
      Reference = character(),
      Niche_Reference = character(),
      Random_Effect_Variable = character(),
      Mutation_WT_Definition = character(),
      Mutation_Case_Definition = character()
    )

    interpretation_summary <- tibble(
      Gene = character(),
      Gene_Name = character(),
      Comparison = character(),
      Reference_Group = character(),
      Comparison_Group = character(),
      Reference_Probability = numeric(),
      Comparison_Probability = numeric(),
      Difference_Percentage_Points = numeric(),
      Direction = character(),
      Significant = character(),
      P_Value = numeric(),
      Interpretation = character()
    )

    adapted_genes <- coef_flat

    return(list(
      coef_flat = coef_flat,
      interpretation_summary = interpretation_summary,
      results_list = results_list,
      adapted_genes = adapted_genes,
      settings = list(
        niche_ref = niche_ref,
        random_effect_var = if (length(random_effect_vars) == 1) random_effect_vars else NULL,
        random_effect_vars = random_effect_vars,
        mut_wt_values = mut_wt_values,
        mut_case_values = mut_case_values,
        genes_to_analyze = genes_to_analyze
      )
    ))
  }
  
  # Flatten coef_list and add Gene (PA tag) column
  coef_flat <- map_df(coef_list_non_null, identity, .id = "Gene")
  
  # Map gene identifier to gene name from the input data
  if (!"Gene" %in% names(coef_flat)) {
    coef_flat$Gene <- character(nrow(coef_flat))
  }

  coef_flat <- coef_flat %>%
    mutate(
      Gene_Name = unname(gene_name_lookup[as.character(Gene)]),
      Gene_Name = ifelse(is.na(Gene_Name) | Gene_Name == "", as.character(Gene), Gene_Name)
    )
  
  # Add reference information
  coef_flat <- coef_flat %>%
    mutate(Reference = reference_name, .after = Gene_Name)

  # Add analysis setting metadata
  coef_flat <- coef_flat %>%
    mutate(
      Niche_Reference = niche_ref,
      Random_Effect_Variable = ifelse(use_random_effect, paste(random_effect_vars, collapse = "|"), "NONE"),
      Mutation_WT_Definition = paste(mut_wt_values, collapse = "|"),
      Mutation_Case_Definition = ifelse(
        is.null(mut_case_values),
        "NOT_WT",
        paste(mut_case_values, collapse = "|")
      )
    )

  interpretation_summary <- map_dfr(
    results_list,
    ~ summarize_gene_effects(.x, df_analysis, gene_name_lookup)
  ) %>%
    mutate(
      Reference = reference_name,
      Niche_Reference = niche_ref,
      Random_Effect_Variable = ifelse(use_random_effect, paste(random_effect_vars, collapse = "|"), "NONE"),
      Mutation_WT_Definition = paste(mut_wt_values, collapse = "|"),
      Mutation_Case_Definition = ifelse(
        is.null(mut_case_values),
        "NOT_WT",
        paste(mut_case_values, collapse = "|")
      ),
      .after = Gene_Name
    )
  
  # Filter genes where ALL environment effects are significant AND negative
  adapted_genes <- coef_flat %>%
    filter(term != "(Intercept)") %>%
    group_by(Gene, Gene_Name) %>%
    filter(all(Estimate < 0) & all(Significant == "Yes")) %>%
    ungroup() %>%
    dplyr::select(Gene, Gene_Name, Reference, everything())
  
  cat("\nAdapted genes (all env. effects significant & negative):\n")
  print(adapted_genes %>% dplyr::select(Gene_Name) %>% distinct())
  cat("\n")
  
  return(list(
    coef_flat = coef_flat,
    interpretation_summary = interpretation_summary,
    results_list = results_list,
    adapted_genes = adapted_genes,
    settings = list(
      niche_ref = niche_ref,
      random_effect_var = if (length(random_effect_vars) == 1) random_effect_vars else NULL,
      random_effect_vars = if (use_random_effect) random_effect_vars else character(0),
      mut_wt_values = mut_wt_values,
      mut_case_values = mut_case_values,
      genes_to_analyze = genes_to_analyze
    )
  ))
}

# ---- Plot Mutation Probabilities with CLD Letters ----
plot_mutation_probabilities_cld <- function(
  glmer_results,
  genes = NULL,
  niche_colors,
  adjust_method = "tukey",
  ncol = 4
) {
  if (is.null(glmer_results$results_list) || length(glmer_results$results_list) == 0) {
    stop("glmer_results$results_list is empty.")
  }

  pick_first_col <- function(df, candidates) {
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) {
      return(NA_character_)
    }
    hit[1]
  }

  gene_labels <- glmer_results$interpretation_summary %>%
    dplyr::select(Gene, Gene_Name) %>%
    dplyr::distinct()

  plot_df <- purrr::imap_dfr(glmer_results$results_list, function(res, gene_id) {
    if (is.null(res) || is.null(res$status) || res$status != "success") {
      return(NULL)
    }

    emm_obj <- tryCatch(
      emmeans::emmeans(res$model, ~ Niche),
      error = function(e) NULL
    )
    if (is.null(emm_obj)) {
      return(NULL)
    }

    emm_df <- tryCatch(
      as.data.frame(summary(emm_obj, type = "response")),
      error = function(e) NULL
    )
    if (is.null(emm_df) || nrow(emm_df) == 0) {
      return(NULL)
    }

    prob_col <- pick_first_col(emm_df, c("prob", "response", "emmean"))
    lower_col <- pick_first_col(emm_df, c("asymp.LCL", "lower.CL", "LCL"))
    upper_col <- pick_first_col(emm_df, c("asymp.UCL", "upper.CL", "UCL"))
    if (is.na(prob_col)) {
      return(NULL)
    }

    out <- emm_df %>%
      dplyr::transmute(
        Gene = as.character(gene_id),
        Niche = as.character(Niche),
        Probability = as.numeric(.data[[prob_col]]),
        Lower = if (!is.na(lower_col)) as.numeric(.data[[lower_col]]) else NA_real_,
        Upper = if (!is.na(upper_col)) as.numeric(.data[[upper_col]]) else NA_real_
      ) %>%
      dplyr::mutate(
        Lower = ifelse(is.na(Lower), Probability, Lower),
        Upper = ifelse(is.na(Upper), Probability, Upper)
      )
    pair_tbl <- tryCatch(
      as.data.frame(
        summary(
          emmeans::contrast(emm_obj, method = "pairwise", adjust = adjust_method),
          infer = c(TRUE, TRUE)
        )
      ),
      error = function(e) NULL
    )

    if (!is.null(pair_tbl) && nrow(pair_tbl) > 0 && "contrast" %in% names(pair_tbl) && "p.value" %in% names(pair_tbl)) {
      pair_tbl <- pair_tbl %>%
        dplyr::mutate(
          Comparison_Group = trimws(sub(" - .*", "", contrast)),
          Reference_Group = trimws(sub(".* - ", "", contrast)),
          pair_name = paste0(Comparison_Group, "-", Reference_Group)
        )

      pvals <- pair_tbl$p.value
      names(pvals) <- pair_tbl$pair_name

      cld_letters <- tryCatch(
        multcompView::multcompLetters(pvals)$Letters,
        error = function(e) NULL
      )

      if (!is.null(cld_letters)) {
        letters_df <- tibble::tibble(
          Niche = names(cld_letters),
          CLD = as.character(cld_letters)
        )
        out <- out %>% dplyr::left_join(letters_df, by = "Niche")
      }
    }

    if (!"CLD" %in% names(out)) {
      out <- out %>% dplyr::mutate(CLD = "")
    }

    out
  }) %>%
    dplyr::left_join(gene_labels, by = "Gene") %>%
    dplyr::mutate(
      Gene_Name = ifelse(is.na(Gene_Name) | Gene_Name == "", Gene, Gene_Name),
      Gene_Label = paste0(Gene_Name)
    )

  if (!is.null(genes)) {
    genes_chr <- as.character(genes)
    plot_df <- plot_df %>%
      dplyr::filter(Gene %in% genes_chr | Gene_Name %in% genes_chr)
    
    # Reorder based on input gene order
    gene_order <- tibble(
      Gene = genes_chr,
      gene_order = seq_along(genes_chr)
    ) %>%
      bind_rows(
        tibble(
          Gene = plot_df %>% 
            dplyr::filter(!(Gene %in% genes_chr)) %>% 
            pull(Gene) %>% 
            unique(),
          gene_order = NA_real_
        )
      )
    
    plot_df <- plot_df %>%
      left_join(gene_order, by = "Gene") %>%
      arrange(gene_order, Gene) %>%
      mutate(Gene_Label = factor(Gene_Label, levels = unique(Gene_Label)))
  }

  if (nrow(plot_df) == 0) {
    stop("No data available to plot after filtering.")
  }

  niche_levels <- names(niche_colors)
  plot_df <- plot_df %>%
    dplyr::mutate(
      Niche = factor(Niche, levels = niche_levels),
      Probability = pmax(0, pmin(1, Probability)),
      Lower = pmax(0, pmin(1, Lower)),
      Upper = pmax(0, pmin(1, Upper))
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Niche, y = Probability, color = Niche)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Lower, ymax = Upper),
      width = 0.15,
      linewidth = 0.4,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::geom_text(
      ggplot2::aes(y = pmin(1, ifelse(is.na(Upper), Probability, Upper) + 0.05), label = CLD),
      color = "black",
      size = 3.2,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(~ Gene_Label, ncol = ncol) +
    ggplot2::scale_color_manual(values = niche_colors, drop = FALSE) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0.02, 0.14))
    ) +
    ggplot2::labs(
      x = "Environment",
      y = "Estimated mutation probability"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )

  list(plot = p, data = plot_df)
}

# ---- Extract Model Details ----
extract_model_details <- function(glmer_results, df_analysis, analysis_name) {
  "Extract model coefficients with additional statistics for supplementary tables"
  
  coef_flat <- glmer_results %>%
    select(coef_flat) %>%
    mutate(
      OR_95CI = paste0(
        "[",
        round(OR_CI_Lower, 3),
        ", ",
        round(OR_CI_Upper, 3),
        "]"
      ),
      Analysis = analysis_name,
      .after = Gene_Name
    ) %>%
    dplyr::select(
      Gene_Name,
      term,
      Estimate = Estimate,
      `Std. Error` = `Std. Error`,
      `z value` = `z value`,
      `Pr(>|z|)` = `Pr(>|z|)`,
      Odds_Ratio,
      OR_95CI,
      n_total = n_obs,
      Analysis
    ) %>%
    mutate(
      Estimate = round(Estimate, 4),
      `Std. Error` = round(`Std. Error`, 4),
      `z value` = round(`z value`, 4),
      `Pr(>|z|)` = round(`Pr(>|z|)`, 4),
      Odds_Ratio = round(Odds_Ratio, 4)
    )
  
  return(coef_flat)
}

# ---- Extract Exact N Per Group ----
extract_exact_n_per_group <- function(df_analysis, analysis_name) {
  "Extract exact sample counts per gene per niche per mutation status"
  
  n_per_group <- df_analysis %>%
    group_by(gene, Gene.Name, Niche) %>%
    summarise(
      n_niche = n(),
      n_mutated_niche = sum(Mutated == 1, na.rm = TRUE),
      n_wt_niche = sum(Mutated == 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(Gene = gene, Gene_Name = Gene.Name) %>%
    mutate(Analysis = analysis_name, .after = Gene_Name) %>%
    arrange(Gene, Niche)
  
  return(n_per_group)
}

# ---- Extract Pairwise Stats ----
extract_pairwise_stats <- function(glmer_results, df_analysis, analysis_name) {
  "Extract pairwise comparisons with sample sizes"
  
  # Get gene-level sample sizes per niche
  n_per_group <- df_analysis %>%
    group_by(gene, Gene.Name, Niche) %>%
    summarise(
      n_niche = n(),
      .groups = "drop"
    ) %>%
    rename(Gene = gene, Gene_Name = Gene.Name)
  
  # Get pairwise stats
  pairwise <- glmer_results$interpretation_summary %>%
    dplyr::select(
      Gene_Name,
      Comparison,
      Reference_Group,
      Comparison_Group,
      Reference_Probability,
      Comparison_Probability,
      Difference_Percentage_Points,
      P_Value,
      Significant
    ) %>%
    left_join(
      n_per_group %>% 
        dplyr::select(Gene_Name, Niche, n_niche) %>%
        rename(Reference_Group = Niche, n_Reference_Group = n_niche),
      by = c("Gene_Name" = "Gene_Name", "Reference_Group" = "Reference_Group")
    ) %>%
    left_join(
      n_per_group %>% 
        dplyr::select(Gene_Name, Niche, n_niche) %>%
        rename(Comparison_Group = Niche, n_Comparison_Group = n_niche),
      by = c("Gene_Name" = "Gene_Name", "Comparison_Group" = "Comparison_Group")
    ) %>%
    mutate(Analysis = analysis_name, .after = Gene_Name) %>%
    arrange(Gene_Name, Comparison)
  
  return(pairwise)
}

# ---- Generate Statistical Methods Summary ----
generate_statistical_methods_summary <- function(glmer_results, df_analysis, exact_n_export_path = NULL, print_exact_n_table = FALSE) {
  settings <- glmer_results$settings
  exact_niche_counts <- tibble::tibble()
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("STATISTICAL METHODS SUMMARY\n")
  cat(strrep("=", 80), "\n\n")
  
  cat("MODEL SPECIFICATION:\n")
  cat("  Family: Binomial with logit link\n")
  cat("  Optimizer: bobyqa\n")
  if (length(settings$random_effect_vars) > 0) {
    cat("  Formula: Mutated ~ Niche + (1 | ", paste(settings$random_effect_vars, collapse = ") + (1 | "), ")\n", sep = "")
  } else {
    cat("  Formula: Mutated ~ Niche\n")
  }
  
  cat("\nREPLICATE DEFINITION:\n")
  cat("  Replicates: Each genome (unique isolate) represents one independent replicate.\n")
  cat("  Total genomes analyzed (CladeA): ", nrow(df_analysis) / nlevels(df_analysis$gene), " (average per gene)\n", sep = "")
  
  cat("\nMUTATION OUTCOME DEFINITION:\n")
  cat("  Mutated: Mut_Status NOT IN ", paste(settings$mut_wt_values, collapse = ", "), "\n", sep = "")
  if (!is.null(settings$mut_case_values)) {
    cat("  (Specifically: Mut_Status IN ", paste(settings$mut_case_values, collapse = ", "), ")\n", sep = "")
  }
  
  cat("\nREFERENCE LEVEL:\n")
  cat("  Niche reference: ", ifelse(is.null(settings$niche_ref), "None (alphabetical order)", settings$niche_ref), "\n", sep = "")
  
  cat("\nRANDOM EFFECTS STRUCTURE:\n")
  if (length(settings$random_effect_vars) > 0) {
    cat("  Variables: ", paste(settings$random_effect_vars, collapse = ", "), "\n", sep = "")
  } else {
    cat("  No random effects (GLM used instead of GLMER)\n")
  }
  
  cat("\nFIXED EFFECTS:\n")
  cat("  Intercept: baseline log-odds of mutation at reference niche level\n")
  cat("  Niche coefficients: log-odds ratio relative to reference niche\n")
  cat("  All reported as z-scores with two-tailed p-values\n")
  
  cat("\nESTIMATED MARGINAL MEANS (emmeans):\n")
  cat("  Reported on: Response scale (probability of mutation, 0-1)\n")
  cat("  Confidence intervals: 95%\n")
  cat("  Pairwise comparisons: Tukey adjustment\n")
  cat("  Error bars in plots: 95% CI on probability scale\n")
  
  cat("\nSAMPLE SIZES:\n")
  coef_tbl <- glmer_results$coef_flat %>% 
    filter(term == "(Intercept)") %>%
    select(Gene, Gene_Name, n_total, n_mutated, n_wt, Converged)
  if (nrow(coef_tbl) > 0) {
    cat("  Genes analyzed: ", nrow(coef_tbl), "\n", sep = "")
    cat("  Range of n per gene: ", min(coef_tbl$n_total, na.rm = T), "-", max(coef_tbl$n_total, na.rm = T), "\n", sep = "")
    cat("  Convergence: ", sum(coef_tbl$Converged == TRUE, na.rm = T), " of ", sum(!is.na(coef_tbl$Converged)), " models converged\n", sep = "")

    gene_name_lookup <- coef_tbl %>%
      distinct(Gene, Gene_Name) %>%
      mutate(Gene_Name = ifelse(is.na(Gene_Name) | Gene_Name == "", Gene, Gene_Name))

    exact_niche_counts <- purrr::imap_dfr(glmer_results$results_list, function(res, gene_id) {
      if (is.null(res) || is.null(res$status) || res$status != "success" || is.null(res$niche_counts)) {
        return(NULL)
      }

      res$niche_counts %>%
        mutate(
          Gene = as.character(gene_id),
          n_wt_niche = n_niche - n_mutated_niche
        ) %>%
        select(Gene, Niche, n_niche, n_mutated_niche, n_wt_niche)
    }) %>%
      left_join(gene_name_lookup, by = "Gene") %>%
      mutate(Gene_Name = ifelse(is.na(Gene_Name) | Gene_Name == "", Gene, Gene_Name)) %>%
      select(Gene, Gene_Name, Niche, n_niche, n_mutated_niche, n_wt_niche) %>%
      arrange(Gene_Name, Gene, Niche)

    if (nrow(exact_niche_counts) > 0) {
      if (isTRUE(print_exact_n_table)) {
        cat("\n  Exact per-gene per-niche sample sizes (from results_list[[gene]]$niche_counts):\n")
        print(as.data.frame(exact_niche_counts), row.names = FALSE)
      }

      if (!is.null(exact_n_export_path) && nzchar(trimws(exact_n_export_path))) {
        write.csv(exact_niche_counts, exact_n_export_path, row.names = FALSE)
        cat("\n  Exported exact sample-size table to: ", exact_n_export_path, "\n", sep = "")
      }
    }
  }
  
  cat("\nP-VALUE REPORTING:\n")
  cat("  All p-values reported exactly (not rounded to p < 0.05 or p > 0.05)\n")
  cat("  Significance threshold: α = 0.05 (two-tailed)\n")
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("\n")

  invisible(exact_niche_counts)
}

# ---- Format Publication Stats ----
format_publication_stats <- function(coef_flat, export_path = NULL) {
  cat("\n")
  cat(strrep("=", 100), "\n")
  cat("PUBLICATION-READY FIXED EFFECTS TABLE (All genes, all terms)\n")
  cat(strrep("=", 100), "\n\n")
  
  pub_table <- coef_flat %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::select(
      Gene_Name, term, Estimate, `Std. Error`, `z value`, `Pr(>|z|)`,
      Odds_Ratio, OR_CI_Lower, OR_CI_Upper,
      n_total, n_mutated, n_wt, Converged
    ) %>%
    dplyr::mutate(
      Estimate = round(Estimate, 4),
      `Std. Error` = round(`Std. Error`, 4),
      `z value` = round(`z value`, 3),
      `Pr(>|z|)` = format(`Pr(>|z|)`, scientific = TRUE, digits = 3),
      Odds_Ratio = round(Odds_Ratio, 3),
      OR_CI_Lower = round(OR_CI_Lower, 3),
      OR_CI_Upper = round(OR_CI_Upper, 3),
      OR_95CI = paste0("[", OR_CI_Lower, ", ", OR_CI_Upper, "]")
    ) %>%
    dplyr::select(
      Gene_Name, term, Estimate, `Std. Error`, `z value`, `Pr(>|z|)`,
      Odds_Ratio, OR_95CI, n_total, n_mutated, n_wt, Converged
    )
  
  print(as.data.frame(pub_table), row.names = FALSE)
  cat("\n")

  if (!is.null(export_path) && nzchar(trimws(export_path))) {
    write.csv(pub_table, export_path, row.names = FALSE)
    cat("  Exported publication stats table to: ", export_path, "\n\n", sep = "")
  }

  invisible(pub_table)
}

# ---- Format Pairwise Comparison Table ----
format_pairwise_comparison_table <- function(interpretation_summary, export_path = NULL) {
  cat("\n")
  cat(strrep("=", 120), "\n")
  cat("PUBLICATION-READY PAIRWISE COMPARISONS TABLE (Tukey-adjusted)\n")
  cat(strrep("=", 120), "\n\n")
  
  comp_table <- interpretation_summary %>%
    dplyr::select(
      Gene_Name, Comparison, Reference_Group, n_Reference_Group,
      Comparison_Group, n_Comparison_Group,
      Reference_Probability, Comparison_Probability,
      Difference_Percentage_Points, `P_Value`, Significant
    ) %>%
    dplyr::mutate(
      Reference_Probability = paste0(round(Reference_Probability * 100, 1), "%"),
      Comparison_Probability = paste0(round(Comparison_Probability * 100, 1), "%"),
      Difference_Percentage_Points = paste0(round(Difference_Percentage_Points, 1), " pp"),
      `P_Value` = format(`P_Value`, scientific = TRUE, digits = 3)
    ) %>%
    dplyr::arrange(Gene_Name, `P_Value`)
  
  print(as.data.frame(comp_table), row.names = FALSE)
  cat("\n")

  if (!is.null(export_path) && nzchar(trimws(export_path))) {
    write.csv(comp_table, export_path, row.names = FALSE)
    cat("  Exported pairwise comparison table to: ", export_path, "\n\n", sep = "")
  }

  invisible(comp_table)
}

# ---- Calculate Constraint Depth ----
calculate_consentrait_depth <- function(niche_value, tree, metadata) {
  
  # Create binary vector: 1 = present, 0 = absent
  # CRITICAL: Must be in exact tree tip order by index
  tip_states <- as.numeric(metadata$Niche == niche_value)
  
  cat(sprintf('\n--- Constraint Depth: %s ---\n', niche_value))
  cat(sprintf('n_present: %d, n_absent: %d\n', sum(tip_states), length(tip_states) - sum(tip_states)))
  
  # Validate: need at least 1 present and 1 absent
  if(sum(tip_states) < 1 || sum(tip_states) == length(tip_states)) {
    warning(sprintf("No variation in trait for %s (all present or all absent)", niche_value))
    return(NULL)
  }
  
  tryCatch({
    # Call consentrait_depth with correct parameters
    result <- castor::consentrait_depth(
      tree = tree,
      tip_states = tip_states,
      min_fraction = 0.9,
      count_singletons = TRUE,
      singleton_resolution = 0,
      weighted = FALSE,
      Npermutations = 10000
    )
    
    cat(sprintf('Mean depth: %.4f\n', result$mean_depth))
    cat(sprintf('Var depth: %.4f\n', result$var_depth))
    cat(sprintf('N positives: %d\n', result$Npositives))
    if(!is.na(result$P)) {
      cat(sprintf('P-value: %.4f\n', result$P))
    }
    
    return(tibble(
      Niche = niche_value,
      Mean_depth = result$mean_depth,
      Var_depth = result$var_depth,
      Min_depth = result$min_depth,
      Max_depth = result$max_depth,
      Npositives = result$Npositives,
      P_value = ifelse(is.null(result$P), NA, result$P),
      Mean_random_depth = result$mean_random_depth
    ))
  }, error = function(e) {
    cat(sprintf('ERROR: %s\n', e$message))
    return(NULL)
  })
}