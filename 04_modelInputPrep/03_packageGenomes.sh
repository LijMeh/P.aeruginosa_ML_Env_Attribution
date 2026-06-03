#!/bin/bash
# ============================================================================
# Package Genomes
# ============================================================================
# Copies train/test set of genomes for Clades A and B and then zips them.
# ============================================================================

# ============================================================================
# Set Paths
# ============================================================================

# Path to Clade A metadata with test/train labels
CladeA="00_data/cladeA_metadata_TestTrain.csv"

# Path to Clade B metadata with test/train labels
CladeB="00_data/cladeB_metadata_TestTrain.csv"

# Source directory containing filtered genomes
src_dir="00_data/genomes/fna_files/filt_genomes"

# Destination directory for organized train/test splits
dest_dir="00_data/genomes/fna_files/Train_Test_Splits"

# ============================================================================
# Copy and Package Genomes
# ============================================================================

mkdir -p "$dest_dir/CladeA_Test" "$dest_dir/CladeA_Train" "$dest_dir/CladeB_Test" "$dest_dir/CladeB_Train"
awk -F',' 'NR>1 {print $1, $7}' "$CladeA" | while read -r genome_id test_train; do
    genome_file="${genome_id}.fna"
    if [[ $test_train == "Test" ]]; then
        cp "$src_dir/$genome_file" "$dest_dir/CladeA_Test/"
    else
        cp "$src_dir/$genome_file" "$dest_dir/CladeA_Train/"
    fi
done

awk -F',' 'NR>1 {print $1, $7}' "$CladeB" | while read -r genome_id test_train; do
    genome_file="${genome_id}.fna"
    if [[ $test_train == "Test" ]]; then
        cp "$src_dir/$genome_file" "$dest_dir/CladeB_Test/"
    else
        cp "$src_dir/$genome_file" "$dest_dir/CladeB_Train/"
    fi
done

find "$dest_dir/CladeA_Test" -type f ! -name '*.gz' -exec pigz {} +
find "$dest_dir/CladeA_Train" -type f ! -name '*.gz' -exec pigz {} +
find "$dest_dir/CladeB_Test" -type f ! -name '*.gz' -exec pigz {} +
find "$dest_dir/CladeB_Train" -type f ! -name '*.gz' -exec pigz {} +