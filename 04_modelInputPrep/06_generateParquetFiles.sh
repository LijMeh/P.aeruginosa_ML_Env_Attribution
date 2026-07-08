#!/bin/bash
# ============================================================================
# Generate Parquet Files for Model Training
# ============================================================================
# Creates the three required parquet files from input data:
#   1. features.parquet - ML features with Genome_ID as index
#   2. splits.parquet   - Train/Val/Test assignments 
#   3. labels.parquet   - Niche and Cluster labels for each genome
# ============================================================================

set -e

echo "Generating parquet files for model training..."

python <<'EOF'
import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================================
# 1. Create features.parquet
# ============================================================================
print("\n[1/3] Creating features.parquet...")

# Load model input data
model_input = pd.read_csv("00_data/model_input/ModelInput_NoDuplicates.csv")

# First column is Genome_ID - set as index
model_input = model_input.set_index(model_input.columns[0])
model_input.index.name = "Genome_ID"

# All remaining columns are features
features_df = model_input

# Save as parquet
features_df.to_parquet("features.parquet")
print(f"  ✓ Saved {features_df.shape[0]} genomes × {features_df.shape[1]} features")

# ============================================================================
# 2. Create splits.parquet
# ============================================================================
print("\n[2/3] Creating splits.parquet...")

# Load the test/train metadata
metadata = pd.read_csv("00_data/cladeA_metadata_TestTrain.csv")

# Create splits dataframe
splits_df = pd.DataFrame(index=features_df.index)

# Map Test_Train to Train/Val/Test
# Using 80% Train, 10% Val, 10% Test stratification
# First, mark all "Test" samples from the csv as "Test"
# Then split the "Train" samples into Train/Val

def create_splits(metadata_df):
    """Create Train/Val/Test splits from metadata."""
    splits = {}
    
    for genome_id in metadata_df['Genome']:
        meta_row = metadata_df[metadata_df['Genome'] == genome_id]
        if meta_row.empty:
            continue
            
        test_train = meta_row['Test_Train'].values[0]
        
        if test_train == 'Test':
            splits[genome_id] = 'Test'
        else:  # Training samples
            # Further split into Train (80%) and Val (20%)
            # Using a deterministic hash-based split
            hash_val = hash(genome_id) % 100
            if hash_val < 80:
                splits[genome_id] = 'Train'
            else:
                splits[genome_id] = 'Val'
    
    return splits

split_dict = create_splits(metadata)

# Create split column
splits_df['split'] = splits_df.index.map(split_dict)

# Fill any missing values (shouldn't be any) with "Train" as default
splits_df['split'].fillna('Train', inplace=True)

# Save as parquet
splits_df.to_parquet("splits.parquet")

split_counts = splits_df['split'].value_counts()
print(f"  ✓ Split distribution:")
for split_type in ['Train', 'Val', 'Test']:
    count = split_counts.get(split_type, 0)
    pct = 100 * count / len(splits_df)
    print(f"    {split_type}: {count} ({pct:.1f}%)")

# ============================================================================
# 3. Create labels.parquet
# ============================================================================
print("\n[3/3] Creating labels.parquet...")

# Load metadata that contains Niche and Cluster information
metadata = pd.read_csv("00_data/blockedModelInputMetadata.csv")

# Create labels dataframe with Genome_ID as index
labels_df = pd.DataFrame(index=features_df.index)
labels_df.index.name = "Genome_ID"

# Map metadata columns to labels
# blockedModelInputMetadata has: Genome, Niche, Clade, Genomvar_cluster, etc.
# We need: Niche (target variable) and Cluster (for holdout validation)

labels_dict_niche = {}
labels_dict_cluster = {}

for _, row in metadata.iterrows():
    genome_id = row['Genome']
    niche = row.get('Niche', 'Unknown')
    # Use Genomvar_cluster as Cluster for genomovar holdout validation
    cluster = row.get('Genomvar_cluster', 'Unknown')
    
    labels_dict_niche[genome_id] = niche
    labels_dict_cluster[genome_id] = cluster

# Assign to dataframe
labels_df['Niche'] = labels_df.index.map(labels_dict_niche)
labels_df['Cluster'] = labels_df.index.map(labels_dict_cluster)

# Remove any rows with missing labels
labels_df = labels_df.dropna()

# Save as parquet
labels_df.to_parquet("labels.parquet")

unique_niches = labels_df['Niche'].nunique()
unique_clusters = labels_df['Cluster'].nunique()
print(f"  ✓ Saved {len(labels_df)} genomes")
print(f"    Unique Niches: {unique_niches}")
print(f"    Unique Clusters: {unique_clusters}")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*60)
print("✓ All parquet files created successfully!")
print("="*60)
print(f"  features.parquet: {features_df.shape}")
print(f"  splits.parquet:   {splits_df.shape}")
print(f"  labels.parquet:   {labels_df.shape}")
print("\nReady for training with trainBlocks.py")

EOF

echo "Done!"
