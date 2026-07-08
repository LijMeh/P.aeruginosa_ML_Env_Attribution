#!/bin/bash
# ============================================================================
# Download Pseudomonas aeruginosa Genome Metadata from NCBI
# ============================================================================
# Downloads assembly metadata for all P. aeruginosa genomes from NCBI using
# the datasets command-line tool. Filters for reference/representative genomes
# and extracts key metadata fields.
# ============================================================================

set -e

# Output directory
OUTPUT_DIR="00_data"
mkdir -p "$OUTPUT_DIR"

OUTPUT_FILE="$OUTPUT_DIR/aniMetadata.csv"

echo "Downloading P. aeruginosa genome metadata from NCBI..."
echo "This may take a few minutes..."

# Fetch metadata for all P. aeruginosa genomes
datasets summary genome taxon 'Pseudomonas aeruginosa' \
  --assembly-source genbank \
  --as-json-lines | \
dataformat tsv genome \
  --fields accession,organism-name,assminfo-biosample-isolation-source,assminfo-biosample-host-disease,assminfo-submitter,assminfo-sequencing-tech,assminfo-release-date,assmstats-number-of-contigs,assmstats-contig-l50,assmstats-contig-n50,assmstats-total-sequence-len,assmstats-total-ungapped-len,checkm-completeness,checkm-contamination > "$OUTPUT_FILE" || {
  echo "ERROR: Failed to download metadata"
  exit 1
}

# Count genomes
count=$(tail -n +2 "$OUTPUT_FILE" | wc -l)

echo ""
echo "✓ Successfully downloaded metadata for $count genomes"
echo "  Output: $OUTPUT_FILE"
echo ""


