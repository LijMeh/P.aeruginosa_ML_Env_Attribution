#!/bin/bash
# ============================================================================
# Assign Clades
# ============================================================================
# This script splits genomes into clades based on the highest ANI value match
# to the three reference genomes (PA7, PA14, and PAO1).
# ============================================================================

# ============================================================================
# Set Paths
# ============================================================================

input_file="/00_data/fastANI_output/reference_Output.txt"
output_file="/00_data/fastANI_output/cladeAssignment.csv"

# ============================================================================
# Identify Clade Groups
# ============================================================================

python3 <<EOF

import pandas as pd

# Load the input file
df = pd.read_csv("$input_file", sep='\t', header=None, names=['ReferenceGenome', 'Genome2', 'ANI', 'Fragment_Matches', 'Total_Fragments'])

# Remove the paths from the ReferenceGenome and Genome2 columns
df['ReferenceGenome'] = df['ReferenceGenome'].apply(lambda x: x.split('/')[-1])
df['Genome2'] = df['Genome2'].apply(lambda x: x.split('/')[-1])

# Group by Genome2 and find the ReferenceGenome with the highest ANI
grouped = df.loc[df.groupby('Genome2')['ANI'].idxmax()]

# Add a Clade column based on the ReferenceGenome
def assign_clade(reference_genome):
    if reference_genome == 'PAO1_GCA_000006765.1_ASM676v1_genomic.fna':
        return 'CladeA'
    elif reference_genome == 'PA14_GCA_000014625.1_ASM1462v1_genomic.fna':
        return 'CladeB'
    elif reference_genome == 'PA7_GCA_000017205.1_ASM1720v1_genomic.fna':
        return 'CladeC'
    else:
        return 'Unknown'

grouped['Clade'] = grouped['ReferenceGenome'].apply(assign_clade)

# Rename Genome2 to Genome and select the required columns
grouped = grouped.rename(columns={'Genome2': 'Genome'})[['Genome', 'Clade', 'ANI']]

# Save the result to a CSV file
grouped.to_csv("$output_file", index=False)

print(grouped)

EOF