#!/bin/bash
# ============================================================================
# Update Metadata
# ============================================================================
# Updates the metadata with Clade, Genomovar, and Clone details
# ============================================================================

# ============================================================================
# Set Paths
# ============================================================================

input_metadata="00_data/rawMetadata.csv"

clade_file="00_data/fastANI_output/cladeAssignment.csv"
clone_file="00_data/fastANI_output/ANI_cluster_99.99_output.csv"
genomovar_file="00_data/fastANI_output/ANI_cluster_99.5_output.csv"

output_metadata="00_data/aniMetadata.csv"

# ============================================================================
# Update Metadata
# ============================================================================
python <<EOF
import pandas as pd

# Load the raw metadata
metadata = pd.read_csv("$input_metadata")

# Load the clade, clone, and genomovar assignments
clade_assignments = pd.read_csv("$clade_file")
clone_assignments = pd.read_csv("$clone_file")
genomovar_assignments = pd.read_csv("$genomovar_file")

# Merge the clade, clone, and genomovar assignments into the metadata
metadata = metadata.merge(clade_assignments[['Genome', 'Clade']], on='Genome', how='left')
metadata = metadata.merge(clone_assignments[['Genome', 'Clone']], on='Genome', how='left')
metadata = metadata.merge(genomovar_assignments[['Genome', 'Genomovar']], on='Genome', how='left') 

# Save the updated metadata
metadata.to_csv("$output_metadata", index=False)

EOF