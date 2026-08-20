# ============================================================================
# Train and Evaluate Final (Non-Blocked) Model
# ============================================================================
# Trains a single XGBoost model on the full gene set that passed SHAP,
# breadth, and study bias filtering in 03_processShapBestModel.R (the "477
# genes"), using ALL training genomes rather than holding any genomovar out.
# The model is evaluated once on the held-out Test split, and the raw
# outputs needed downstream (performance metrics, ROC curves, per-genome
# class probabilities) are written for 07_processDataForFigs.RMD to process.
# ============================================================================

# ---- Load Required Packages ----
library(install.load)
install_load("here")

# ---- Load configuration and setup ----
source(here("00_R_configs", "config.R"))
source(here("00_R_configs", "00_setup.R"))
source(here("00_R_configs", "00_functions.R"))

# ============================================================================
# Load Filtered Feature Set ("477 genes")
# ============================================================================

# niche_results.rds is the output of 03_processShapBestModel.R:
# combined_data per niche already contains only genes that passed the
# bias and breadth filters, so no additional filtering is needed here.
niche_results_filtered <- readRDS(data_here("00_data/model_output", "niche_results.rds"))

top_genes <- unique(unlist(lapply(niche_results_filtered, function(nr) {
  nr$combined_data %>% pull(gene)
})))

write.csv(
  data.frame(gene = top_genes),
  "00_data/model_output/full_model_results/final_model_genes.csv",
  row.names = FALSE
)

cat(sprintf("Using %d genes for final model\n", length(top_genes)))

# ============================================================================
# Load Train/Test Split and Niche Labels
# ============================================================================

test_train_split <- read.csv(data_here("00_data", "testTrainMetadata.csv")) %>%
  filter(Clade == "CladeA") %>%
  rename(Genome_ID = Genome) %>%
  mutate(Genome_ID = str_replace_all(Genome_ID, "\\.", "_"))

train_genomes <- test_train_split %>%
  filter(Test_Train == "Train") %>%
  pull(Genome_ID)

test_genomes <- test_train_split %>%
  filter(Test_Train == "Test") %>%
  pull(Genome_ID)

# ---- Niche remapping (matches 03_processShapBestModel.R / 05_addGeneDataToSHAP.R) ----
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

metadata_remapped <- test_train_split %>%
  filter(Niche %in% names(niche_map)) %>%
  mutate(Niche = niche_map[Niche])

# ============================================================================
# Load and Filter Model Input
# ============================================================================

model_input <- fread(
  data_here("00_data", "model_input", "ModelInput_NoDuplicates.csv"),
  select = c("Genome_ID", top_genes)
) %>%
  filter(Genome_ID %in% c(train_genomes, test_genomes))

modeling_df <- model_input %>%
  left_join(metadata_remapped %>% select(Genome_ID, Niche), by = "Genome_ID") %>%
  filter(!is.na(Niche))

train_df <- modeling_df %>% filter(Genome_ID %in% train_genomes)
test_df <- modeling_df %>% filter(Genome_ID %in% test_genomes)

cat(sprintf("Train: %d genomes, Test: %d genomes\n", nrow(train_df), nrow(test_df)))

# ============================================================================
# Prepare Matrices and Labels
# ============================================================================

X_train <- train_df %>% select(all_of(top_genes)) %>% as.matrix()
X_test <- test_df %>% select(all_of(top_genes)) %>% as.matrix()

y_train <- as.factor(train_df$Niche)
y_test <- as.factor(test_df$Niche)

all_levels <- sort(unique(c(levels(y_train), levels(y_test))))
y_train <- factor(y_train, levels = all_levels)
y_test <- factor(y_test, levels = all_levels)

y_train_int <- as.integer(y_train) - 1
y_test_int <- as.integer(y_test) - 1

dtrain <- xgb.DMatrix(data = X_train, label = y_train_int)
dtest <- xgb.DMatrix(data = X_test, label = y_test_int)

# ============================================================================
# Train Final XGBoost Model
# ============================================================================
# NOTE: uses multi:softprob (not multi:softmax) since downstream ROC/raw
# probability outputs require per-class probabilities, not hard labels.
# Hyperparameters below match the defaults used elsewhere in this repo;
# swap in the winning round's tuned hyperparameters if those are recovered
# from the optuna study, rather than re-tuning here.

params <- list(
  objective = "multi:softprob",
  num_class = length(all_levels),
  eval_metric = "mlogloss",
  max_depth = 6,
  eta = 0.3,
  nthread = 30
)

xgb_fit <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 200,
  watchlist = list(train = dtrain, test = dtest),
  verbose = 1,
  early_stopping_rounds = 20
)

# ============================================================================
# Predict on Test Set
# ============================================================================

raw_preds <- predict(xgb_fit, X_test)
probs_test <- matrix(raw_preds, ncol = length(all_levels), byrow = TRUE)
colnames(probs_test) <- all_levels

preds_int <- max.col(probs_test) - 1
preds_niche <- all_levels[preds_int + 1]

# ============================================================================
# Performance Metrics
# ============================================================================
# Matches the schema written by trainBlocks.py's `eval()` (Micro F1/Macro
# F1/Accuracy/Balanced Accuracy columns, True/Predicted Labels as
# Python-style stringified lists) so the final model is directly comparable
# to the blocked-model holdout CSVs read in 07_processDataForFigs.RMD.
#
# "Micro F1" matches trainBlocks.py's micro_f1_score(), which despite the
# name is sklearn's `f1_score(average="weighted")`, not a true micro
# average - kept for consistency with the blocked model results.

metrics_df <- tibble(
  truth = y_test,
  prediction = factor(preds_niche, levels = all_levels)
)

micro_f1 <- yardstick::f_meas(metrics_df, truth, prediction, estimator = "macro_weighted")$.estimate
macro_f1 <- yardstick::f_meas(metrics_df, truth, prediction, estimator = "macro")$.estimate
accuracy <- yardstick::accuracy(metrics_df, truth, prediction)$.estimate

# Chance-adjusted balanced accuracy, matching sklearn's
# balanced_accuracy_score(adjusted=True): (macro recall - chance) / (1 - chance)
recall_per_class <- yardstick::recall(metrics_df, truth, prediction, estimator = "macro")$.estimate
chance <- 1 / length(all_levels)
balanced_accuracy <- (recall_per_class - chance) / (1 - chance)

cat(sprintf("\nFinal Model Test Set Performance:\n"))
cat(sprintf("Micro F1 (weighted): %.3f\n", micro_f1))
cat(sprintf("Macro F1: %.3f\n", macro_f1))
cat(sprintf("Accuracy: %.3f\n", accuracy))
cat(sprintf("Balanced Accuracy (adjusted): %.3f\n", balanced_accuracy))

py_list_str <- function(x) {
  paste0("[", paste0("'", x, "'", collapse = ", "), "]")
}

final_test_results <- data.frame(
  `Micro F1` = micro_f1,
  `Macro F1` = macro_f1,
  Accuracy = accuracy,
  `Balanced Accuracy` = balanced_accuracy,
  `True Labels` = py_list_str(as.character(y_test)),
  `Predicted Labels` = py_list_str(preds_niche),
  check.names = FALSE
)

dir.create("00_data/model_output/full_model_results", recursive = TRUE, showWarnings = FALSE)

write.csv(
  final_test_results,
  "00_data/model_output/full_model_results/477_filtered_test.csv",
  row.names = FALSE
)

# ============================================================================
# ROC Curves and Raw Probabilities (one-vs-rest per niche)
# ============================================================================

roc_list <- lapply(all_levels, function(niche) {
  response <- as.numeric(y_test == niche)
  roc_obj <- pROC::roc(response = response, predictor = probs_test[, niche], quiet = TRUE)
  list(
    fpr = rev(1 - roc_obj$specificities),
    tpr = rev(roc_obj$sensitivities)
  )
})
names(roc_list) <- all_levels

max_len <- max(vapply(roc_list, function(r) length(r$fpr), integer(1)))

pad_na <- function(v, len) c(v, rep(NA_real_, len - length(v)))

fpr_data <- as.data.frame(lapply(roc_list, function(r) pad_na(r$fpr, max_len)))
tpr_data <- as.data.frame(lapply(roc_list, function(r) pad_na(r$tpr, max_len)))
colnames(fpr_data) <- all_levels
colnames(tpr_data) <- all_levels

dir.create("00_data/model_output/roc", recursive = TRUE, showWarnings = FALSE)

write.csv(fpr_data, "00_data/model_output/roc/477model_fpr.csv")
write.csv(tpr_data, "00_data/model_output/roc/477model_tpr.csv")

raw_probs <- as.data.frame(probs_test)
raw_probs$Genome_IDs <- test_df$Genome_ID
raw_probs <- raw_probs %>% select(Genome_IDs, everything())

write.csv(raw_probs, "00_data/model_output/roc/477model_raw_probs.csv", row.names = FALSE)
