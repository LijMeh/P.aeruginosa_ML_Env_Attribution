#!/bin/bash
# ============================================================================
# Script Runs fastANI
# ============================================================================
# This script takes the genomes and runs fastANI to compare them to each other.
# If some genomes have already been run, it will only add to the existing output.
# ============================================================================

# ============================================================================
# Set Paths
# ============================================================================

# Directory containing the new .fna files to use as query
source_directory="/00_data/genomes/fna_files/filt_genomes"

# Directory to store the fastANI output
output_directory="/00_data/fastANI_output"

# Directory with genomes already run through fastANI (optional)
# If provided, query genomes will be compared against existing + new genomes as reference.
# If left empty, pairwise fastANI will run on all genomes in source_directory only.
existing_genome_directory="/00_data/existing_genomes"

# ============================================================================
# Run fastANI
# ============================================================================

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Function to find all .fna files in a directory and save their paths to a .txt file
find_fna_files() {
    local directory=$1
    local output_file=$2
    find "$directory" -type f -name "*.fna" > "$output_file"
}

# Find all .fna files in the source directory and save their paths to a .txt file
fna_files_list="$output_directory/fastANI_input_genomes_list.txt"
find_fna_files "$source_directory" "$fna_files_list"

if [ -n "$existing_genome_directory" ]; then
    # Find all .fna files in the existing genome directory and save their paths to a .txt file
    existing_genomes_list="$output_directory/existing_fastANI_genomes_list.txt"
    find_fna_files "$existing_genome_directory" "$existing_genomes_list"

    # Append the new genome paths to the existing genomes list to use as reference
    cat "$fna_files_list" >> "$existing_genomes_list"
fi

# Set query and reference lists
query_list="$output_directory/fastANI_input_genomes_list.txt"

if [ -n "$existing_genome_directory" ]; then
    reference_list="$output_directory/existing_fastANI_genomes_list.txt"
else
    reference_list="$output_directory/fastANI_input_genomes_list.txt"
fi

# Run fastANI with the query list against the reference list
fastANI --ql "$query_list" --rl "$reference_list" -o "$output_directory/fastani_output.txt" -t 24