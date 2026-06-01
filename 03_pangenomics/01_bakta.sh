#!/bin/bash
# ============================================================================
# Script Runs Bakta
# ============================================================================
# This script runs Bakta on the genomes that passed the checkM quality filter. 
# ============================================================================

# ============================================================================
# Set Paths 
# ============================================================================

pathToGenomes="/00_data/genomes/fna_files/filt_genomes"
baktaDB="/00_data/bakta_db/db"
baktaOutput="/00_data/bakta_annotations"

# ============================================================================
# Run Bakta
# ============================================================================

cd $pathToGenomes # location of fasta files

# Specify the number of cores you want to use per Bakta run
coresPerRun=2

# Get the number of available cores and divide by the number of cores per run
numCores=$(($(nproc) / coresPerRun))

export baktaDB baktaOutput

# Use GNU Parallel to run the jobs in parallel
ls *.fna | parallel -j $numCores 'basename=$(basename {} .fna); echo "Processing: $basename"; bakta --db $baktaDB --threads '"$coresPerRun"' -o "$baktaOutput"/"${basename}" --force {}'