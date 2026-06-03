#!/bin/bash
# ============================================================================
# Split Test Train
# ============================================================================
# Randomly selects 80% of genomes for training and 20% for testing. Adds
# a column to metadata to indicate train/test assignment for each clade.
# The split is performed separately for each niche to ensure representation 
# in both sets.
# ============================================================================

# ============================================================================
# Set Paths
# ============================================================================

# Path to metadata
metadata="00_data/aniMetadata.csv"

# Output file for Clade A with test/train labels
CladeA_OutputFile="00_data/cladeA_metadata_TestTrain.csv"

# Output file for Clade B with test/train labels
CladeB_OutputFile="00_data/cladeB_metadata_TestTrain.csv"

# Output metadata
output_metadata="00_data/finalMetadata.csv"

# ============================================================================
# Split Test Train
# ============================================================================

python <<EOF
import pandas as pd
from sklearn.model_selection import train_test_split

# Load the data
metadata = pd.read_csv("$metadata")

# Split metadata by Clade column
cladeA = metadata[metadata['Clade'] == 'CladeA']
cladeB = metadata[metadata['Clade'] == 'CladeB']

# Initialize an empty DataFrame for the final output
outputCladeA = pd.DataFrame()
outputCladeB = pd.DataFrame()

# Split the data by the 'Niche' column for Clade A
for Niche in cladeA['Niche'].unique():
    Niche_data = cladeA[cladeA['Niche'] == Niche]
    train, test = train_test_split(Niche_data, test_size=0.2, random_state=42)
    train['Test_Train'] = 'Train'
    test['Test_Train'] = 'Test'
    outputCladeA = pd.concat([outputCladeA, train, test])

# Split the data by the 'Niche' column for Clade B
for Niche in cladeB['Niche'].unique():
    Niche_data = cladeB[cladeB['Niche'] == Niche]
    train, test = train_test_split(Niche_data, test_size=0.2, random_state=42)
    train['Test_Train'] = 'Train'
    test['Test_Train'] = 'Test'
    outputCladeB = pd.concat([outputCladeB, train, test])

outputCladeA.to_csv("$CladeA_OutputFile", index=False)
outputCladeB.to_csv("$CladeB_OutputFile", index=False)

# Combine test/train assignments from both clades and merge with original metadata
outputCombined = pd.concat([outputCladeA[['Genome', 'Test_Train']], outputCladeB[['Genome', 'Test_Train']]])
metadata = metadata.merge(outputCombined, on='Genome', how='left')
metadata.to_csv("$output_metadata", index=False)

EOF
