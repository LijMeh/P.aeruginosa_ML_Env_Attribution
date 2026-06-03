import pandas as pd
import numpy as np
from collections import defaultdict
import hashlib
import sys

if len(sys.argv) < 3:
    print("Usage: python 00_processPanarooOutput.py <input_path> <output_path> [rtab_path]")
    sys.exit(1)

input_path = sys.argv[1]
output_path = sys.argv[2]
# Set default rtab_path, allow override with optional argument
if len(sys.argv) > 3:
    rtab_path = sys.argv[3]
else:
    rtab_path = input_path + '/gene_presence_absence.Rtab'

# Function to read the CSV in chunks and transpose
def load_and_transpose_csv(file_path, chunk_size=100000):
    chunks = []
    for chunk in pd.read_csv(file_path, 
                             delimiter='\t', 
                             chunksize=chunk_size, 
                             index_col=0):
        chunks.append(chunk.T)
    return pd.concat(chunks, axis=1)

# Load and transpose the data
df = load_and_transpose_csv(rtab_path)

# Ensure the data contains only 0s and 1s
df = df.astype(bool).astype(int)

# Function to hash a column
def hash_column(column):
    return hashlib.md5(column.tobytes()).hexdigest()

# Create a dictionary to store gene groups
gene_groups = defaultdict(list)

print(df)

# Hash each column and group identical hashes
for gene in df.columns:
    column_hash = hash_column(df[gene].values)
    gene_groups[column_hash].append(gene)

# Collect columns for the new dataframe for gene groups
gene_group_columns = []

# Create a list to store gene group information
gene_group_info = []

# Counter for multi-gene groups
multi_gene_group_counter = 1

# Populate the gene group columns and info
for hash_value, genes in gene_groups.items():
    if len(genes) == 1:
        # For single-gene groups, keep the original gene name
        group_name = genes[0]
    else:
        # For multi-gene groups, use Gene_Group_X naming
        group_name = f"Gene_Group_{multi_gene_group_counter}"
        multi_gene_group_counter += 1
    
    gene_group_columns.append(pd.Series(df[genes[0]], name=group_name))  # Use the first gene's pattern for the group
    
    # Add information about this gene group
    gene_group_info.append({
        'Gene_Group': group_name,
        'Genes': ', '.join(genes),
        'Number_of_Genes': len(genes)
    })

# Create the gene group dataframe
gene_group_df = pd.concat(gene_group_columns, axis=1)

# Create a dataframe with gene group information
gene_group_info_df = pd.DataFrame(gene_group_info)

# Rename the index to Genome_ID and remove .fna from all values
gene_group_df.index = gene_group_df.index.str.replace('.fna', '', regex=False)
gene_group_df.index.name = 'Genome_ID'

# Save the final dataframe
gene_group_df.to_csv(output_path + '/grouped_gene_presence_absence.csv')

# Save the gene info dataframe
gene_group_info_df.to_csv(output_path + '/grouped_gene_info.csv')

# Calculate the percentage of 1s for each column
presence_percentage = (gene_group_df.sum() / len(gene_group_df)) * 100

# Filter out genes that are present in less than .1% or more than 99.9% of genomes
genes_to_keep = presence_percentage[(presence_percentage >= .1) & (presence_percentage <= 99.9)].index

gene_group_df_filtered = gene_group_df[genes_to_keep]

# Save the final filtered dataframe
gene_group_df_filtered.to_csv(output_path + '/ModelInput.csv')