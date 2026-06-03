#!/bin/bash
# ============================================================================
# Prep metadata for use with model
# ============================================================================
# This script assigns Env grouping for use in the model, and filters out
# the final set of genomes, such that, before running the model, the input is 
# filtered to ONLY genomes in this metadata file. 
# ============================================================================

# ============================================================================
# Set Paths
# ============================================================================

# Path to metadata with test/train labels
metadata="00_data/testTrainMetadata.csv"

# Env Block Groupings
blocked_metadata="00_data/envBlockedMetadata.csv"

# Output metadata for model input
model_input_metadata="00_data/blockedModelInputMetadata.csv"

# ============================================================================
# Add Environment Clusters
# ============================================================================

python 06_addEnvClusters.py "$metadata" "$blocked_metadata"

# ============================================================================
# Filter Metadata for Blocked Model Input
# ============================================================================

python <<EOF
import pandas as pd

# Load blocked metadata with environment clusters
metadata = pd.read_csv("$blocked_metadata")

# Filter for CladeA only
metadata = metadata[metadata['Clade'] == 'CladeA']

# Filter for Training set only
metadata = metadata[metadata['Test_Train'] == 'Train']

# Filter to remove NAs in Env_Cluster_1 (dropped environments)
metadata = metadata[metadata['Env_Cluster_1'].notna()]

# Export filtered metadata
metadata.to_csv("$model_input_metadata", index=False)

print(f"Filtered metadata exported: {len(metadata)} genomes")

EOF

