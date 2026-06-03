#!/bin/bash
# ============================================================================
# Run Panaroo
# ============================================================================
# This script runs a custom implementation of Panaroo on a set of GFF files.
# Specifically, it enforces a STRICT 100% protein identity clustering threshold
# by turning off the refind mode and manually setting strict clustering parameters
# throughout the pipeline, including in sections that are not usually exposed to
# the user. 
# ============================================================================

# ============================================================================
# Move GFF files to a single directory
# ============================================================================

mkdir -p /00_data/genomes/gff

find /00_data/bakta_annotations -name "*.gff3" -type f -exec mv {} /00_data/genomes/gff/ \;

# ============================================================================
# Run Panaroo
# ============================================================================

eval "$(micromamba shell hook --shell bash)"

micromamba activate Panaroo_ucf

panaroo \
    -i /00_data/genomes/gff/*.gff3 \
    -o /00_data/panarooOutput \
    --clean-mode strict \
    --threshold 1.0 \
    --family_threshold 1.0 \
    --len_dif_percent 1.0 \
    --family_len_dif_percent 1.0 \
    --refind-mode off \
    --threads 200 \
    --merge_paralogs \
    --remove-invalid-genes
