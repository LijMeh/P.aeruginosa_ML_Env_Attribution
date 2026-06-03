import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import argparse
from joblib import Parallel, delayed

# Parse command line arguments
parser = argparse.ArgumentParser(description='Process fastANI output and cluster genomes.')
parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input fastANI file as a .csv')
parser.add_argument('-o', '--output_path', type=str, required=True, help='Path to save the processed output files')
parser.add_argument('-t', '--cluster_threshold', type=float, required=True, help='ANI threshold for clustering')
parser.add_argument('-p', '--num_threads', type=int, required=True, help='Number of threads to use for processing')
args = parser.parse_args()

# Read the input file into a DataFrame
df = pd.read_csv(args.input_file)

# Define the required columns
required_columns = ['Genome1', 'Genome2', 'ANI', 'Fragment_Matches', 'Total_Fragments'] 

# Check if the number of columns matches the required columns
if len(df.columns) != len(required_columns):
    raise ValueError(f"Input file must contain exactly {len(required_columns)} columns, but {len(df.columns)} columns were found.")

# Assign the required columns as the DataFrame's columns
df.columns = required_columns
print(f"Columns set to: {', '.join(df.columns)}")

# Function to extract the base name without extension
def get_base_name(path):
    return os.path.splitext(os.path.basename(path))[0]

# Apply the function to the Genome1 and Genome2 columns
df['Genome1'] = df['Genome1'].apply(get_base_name)
df['Genome2'] = df['Genome2'].apply(get_base_name)

# Print the processed DataFrame
print(df)

# Optionally, save the processed DataFrame to a new file
#output_file = os.path.join(args.output_path, 'processed_fastANI_output.csv')
#df.to_csv(output_file, index=False)

# Create a graph
G = nx.Graph()

# Function to determine edges to add
def determine_edges(row):
    if row['ANI'] >= args.cluster_threshold:
        return (row['Genome1'], row['Genome2'])
    return None

# Use parallel processing to determine edges
edges = Parallel(n_jobs=args.num_threads)(delayed(determine_edges)(row) for index, row in df.iterrows())

# Add edges to the graph
for edge in edges:
    if edge is not None:
        G.add_edge(*edge)

# Add nodes for genomes that do not meet the threshold with any other genome
for genome in set(df['Genome1']).union(set(df['Genome2'])):
    if genome not in G:
        G.add_node(genome)

# Find connected components (clusters)
clusters = list(nx.connected_components(G))

# Create a DataFrame to store the cluster information
cluster_data = []
for cluster_id, cluster in enumerate(clusters, start=1):
    for genome in cluster:
        cluster_data.append({'Genome': genome, 'Cluster': f'Cluster_{cluster_id}'})

cluster_df = pd.DataFrame(cluster_data)

# Print the number of genomes per cluster
genomes_per_cluster = cluster_df['Cluster'].value_counts()
# Print the whole table (ignore print limits)
with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(genomes_per_cluster)

# Define the output file path for the cluster DataFrame
cluster_output_file = os.path.join(args.output_path, f'ANI_cluster_{args.cluster_threshold}_output.csv')

# Export the cluster DataFrame to a CSV file
cluster_df.to_csv(cluster_output_file, index=False)