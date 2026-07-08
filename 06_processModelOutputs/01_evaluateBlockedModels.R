# ============================================================================
# Process Model Performance from Blocked Models
# ============================================================================
# This script processes model performance metrics to determine which environments
# to use for the final "blocked" model which gets preliminary features. 
# This version uses the verified no-validation set data
# ============================================================================

# ---- Load Required Packages ----
library(install.load)
install_load("here")

# ---- Load configuration and setup ----
source(here("00_R_configs", "config.R"))
source(here("00_R_configs", "00_setup.R"))
source(here("00_R_configs", "00_color_palette.R"))
source(here("00_R_configs", "00_functions.R"))

# ---- Load Blocked Data ----
localPath <- "00_data/Blocked_Cross_Results"

firstRound <- read.csv(data_here(localPath, "first_round_holdout_stats_newPanaroo.csv"))
secondRound <- read.csv(data_here(localPath, "second_round_holdout_stats_newPanaroo.csv"))
thirdRound <- read.csv(data_here(localPath, "third_round_holdout_stats_newPanaroo.csv"))
fourthRound <- read.csv(data_here(localPath, "fourth_round_holdout_stats_newPanaroo.csv"))
fifthRound <- read.csv(data_here(localPath, "fifth_round_holdout_stats_newPanaroo.csv"))
sixthRound <- read.csv(data_here(localPath, "sixth_round_holdout_stats_newPanaroo.csv"))


# ---- Parse Accuracies Per Round ----
parsed_firstRound_Envs <- firstRound %>% 
  calc_env_accuracy() %>% 
  calc_normalized_accuracy() %>%
  mutate(
    ENV = ifelse(ENV == "Gastrointestinal", "Rectal-Feces", ENV)
  )

parsed_firstRound_Envs %>% filter(normalized_accuracy > 0)

parsed_secondRound_Envs <- secondRound %>% 
  calc_env_accuracy() %>% 
  calc_normalized_accuracy()

parsed_secondRound_Envs %>% filter(normalized_accuracy > 0)

parsed_thirdRound_Envs <- thirdRound %>% 
  calc_env_accuracy() %>% 
  calc_normalized_accuracy()

parsed_thirdRound_Envs %>% filter(normalized_accuracy > 0)

parsed_fourthRound_Envs <- fourthRound %>% 
  calc_env_accuracy() %>% 
  calc_normalized_accuracy()

parsed_fourthRound_Envs %>% filter(normalized_accuracy > 0)

parsed_fifthRound_Envs <- fifthRound %>% 
  calc_env_accuracy() %>% 
  calc_normalized_accuracy()

parsed_fifthRound_Envs %>% filter(normalized_accuracy > 0)

parsed_sixthRound_Envs <- sixthRound %>% 
  calc_env_accuracy() %>% 
  calc_normalized_accuracy()

parsed_sixthRound_Envs %>% filter(normalized_accuracy > 0)

# ============================================================================
# Conclusions Summary
# ============================================================================
# Round 1 Conclusions: Clearly suboptimal, many categories with negative normalized accuracy
#
# Round 2-4 Conclusions:
#
# - Pneumonia does better combined into Other Lung than as its own category, or as "Other Clinical"
#   - Round 2: Pneumonia: -0.0275, Other Lung: 0.0815
#   - Round 3: Other Clinical (with Pneumonia): -0.000858, Other Lung: 0.116
#   - Round 4: Other Lung (with Pneumonia): 0.209
# - Aquatic doesn't really do better or worse when Ocean is added to it
#   - Round 2: Aquatic: 0.0741, Ocean: -0.0449
#   - Round 3: Aquatic: 0.0533, Other Environmental: -0.0692
#   - Round 4: Aquatic (with Ocean): 0.0470
#   - Probably will keep as it removes the Ocean category which did poorly alone
# - Hospital does worse when combined with Human Environment broadly
#   - Round 2: Hospital: 0.338, Built Environment: -0.0588, Waste Water: -0.0588
#   - Round 3: Hospital: 0.279, Human Environment (Built + Waste Water): -0.0769
#   - Round 4: Human Environment (with Hospital + Built + Waste Water): 0.207
#   - Keep as is, leave Human Environment from Round 2?, no DROP
# - Other Clinical does bad no matter what, DROP
# - Animal is barely at the cutoff, DROP
# - Other Environmental is bad, with Ocean or with Terrestrial alone (so move ocean to Aquatic), DROP
#
# Round 5-6 Conclusions:
# - Round 6 is best performing, Pneumonia had to be removed. All other niches (besides Aquatic? Improved as a result)
# ============================================================================
