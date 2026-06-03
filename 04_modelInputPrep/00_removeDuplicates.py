import pandas as pd
import numpy as np
import sys
import os

# Usage: python 00_removeDuplicates.py <path>
if len(sys.argv) != 2:
	print("Usage: python 00_removeDuplicates.py <path>")
	sys.exit(1)

path = sys.argv[1]

# Construct file paths from input and output paths
input_file = os.path.join(path, "ModelInput.csv")
output_file_1 = os.path.join(path, "Unique_Feature_Duplicates.csv")
output_file_2 = os.path.join(path, "ModelInput_NoDuplicates.csv")

# Load the CSV file
df = pd.read_csv(input_file, low_memory=False)

# Identify the columns to compare (13 to end, assuming first column is 1)
columns_to_compare = df.columns[1:]  # Python uses 0-based indexing, so 12 corresponds to column 13

# Create a hash of the values in the specified columns for each row
df['hash'] = df[columns_to_compare].apply(lambda x: hash(tuple(x)), axis=1)

# Group by the hash and assign a Feature_Duplicate number
df['Feature_Duplicate'] = df.groupby('hash').ngroup()

# Create a new dataframe with Genome_ID, Feature_Duplicate, and additional columns
Feature_Duplicates_df = df[['Genome_ID', 'Feature_Duplicate']]

Feature_Duplicates_df.to_csv(output_file_1, index=False)
print(f"Clonal groups with additional information have been saved to {output_file_1}")

# Create a version of the input dataframe where all Feature_Duplicates (anything that shows up more than once) is removed
Feature_Duplicate_counts = df['Feature_Duplicate'].value_counts()
unique_Feature_Duplicates = Feature_Duplicate_counts[Feature_Duplicate_counts == 1].index
df_unique = df[df['Feature_Duplicate'].isin(unique_Feature_Duplicates)].drop(['hash', 'Feature_Duplicate'], axis=1)

df_unique.to_csv(output_file_2, index=False)
print(f"XGBoost Input Dataframe with unique features has been saved to {output_file_2}")

# Print some statistics
total_genomes = len(df)
unique_groups = df['Feature_Duplicate'].nunique()
unique_genomes = len(df_unique)
print(f"Total number of genomes: {total_genomes}")
print(f"Number of unique clonal groups: {unique_groups}")
print(f"Number of genomes after removing duplicates: {unique_genomes}")