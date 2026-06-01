#!/bin/bash
# ============================================================================
# Assign Genomovars
# ============================================================================
# This script takes the fastANI output and identifies genomovars at a 99.5%
# ANI threshold (based on Conrad et al. 2024, Nat Coms).
# ============================================================================

# ============================================================================
# Set Paths
# ============================================================================

# Input fastANI output file
fastani_output_file="/00_data/fastANI_output/processed_merged_fastANI_output.csv"

# Output path for filtered fastANI output
filtered_fastani_output_path="/00_data/fastANI_output"

# ============================================================================
# Identify Genomovars
# ============================================================================

# ANI threshold to cluster
ani=99.5

python clusterANI.py -i "$fastani_output_file" -o "$filtered_fastani_output_path" -t "$ani"
