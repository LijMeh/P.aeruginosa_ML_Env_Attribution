# ============================================================================
# Script to filter genome metadata for P. aeruginosa
# ============================================================================
# This script takes the metadata from NCBI datasets for P. aeruginosa and filters
# it based on quality metrics and environment of isolation metadata. 
# ============================================================================

# ============================================================================
# Import Libraries 
# ============================================================================

import pandas as pd
import matplotlib.pyplot as plt

# ============================================================================
# Load Data 
# ============================================================================

# Specify the path to the TSV file
dat_tsv = "/00_data/PA_30K_GenomeMetadata.tsv"

# Load the TSV file into a DataFrame with tab as the delimiter
dat = pd.read_csv(dat_tsv, sep='\t')

# Function to convert column names
def convert_column_names(df):
    df.columns = [col.replace(' ', '.').replace('-', '.') for col in df.columns]
    df.columns = [col[:-1] if col.endswith('.') else col for col in df.columns]
    return df

# Apply the function
dat = convert_column_names(dat)

# ============================================================================
# Initial Filter on Quality Metrics 
# ============================================================================

contigs = 300  # Set the maximum number of contigs allowed
n50 = 100000  # Set the minimum N50 value required


# Filter the DataFrame to retain rows where the number of contigs is less than or equal to the specified value
# and the N50 value is greater than or equal to the specified value
dat_filt = dat[(dat['Assembly.Stats.Number.of.Contigs'] <= contigs) & (dat['Assembly.Stats.Contig.N50'] >= n50)]

# ============================================================================
# Define Niche Distribution from Metadata 
# ============================================================================

dat_passed = dat_filt

# Define a dictionary to map keywords to categories
niche_mapping = {
    'early.CF': ['child cystic', 'cystic fibrosis infant'],
    'CF': ['cf', 'cystic', 'fibrosis'],
    'Sputum': ['sputum'],
    'Blood': ['bacteremia', 'blood', 'sepsis'],
    'Urinary': ['urin', 'catheter'],
    'Ear': ['ear'],
    'Throat': ['throat'],
    'Rectal/Feces': ['feces', 'rectal'],
    'Pneumonia': ['pneumonia'],
    'Lung': ['lung', 'bronchial', 'lower respiratory', 'bronch', 'trach', 'respiratory', 'airway', 'alveolar', 'copd'],
    'Burn': ['burn'],
    'Wound': ['wound', 'sore', 'abscess', 'tissue', 'ulcer', 'abcess', 'lesion'],
    'Skin': ['skin'],
    'Eye': ['eye', 'cornea', 'conjunctiva'],
    'Upper.Respiratory': ['nose', 'nasopharynx', 'sinus', 'nasal'],
    'Ocean': ['ocean', 'sea', 'submarine'],
    'Aquatic': ['water', 'pipe', 'river', 'spring', 'lilypad', 'puddle', 'saltern', 'pool', 'pond', 'stream'],
    'Waste.Water': ['runoff', 'pharma', 'waste', 'sewage', 'sludge'],
    'Hospital.Enviornment': ['hospital'],
    'Built.Enviornment': ['shower', 'kitchen', 'drain', 'sink', 'toilet', 'drinking', 'fountain', 'humidifier'],
    'Terrestrial': ['soil', 'sand', 'dirt'],
    'Animal': ['animal', 'dog', 'cat', 'fish', 'dolphin', 'chicken', 'poultry', 'canine', 'pork', 'pig', 'veterinary'],
    'Plant': ['plant', 'root', 'rhizo', 'tree', 'pepper', 'tomato', 'potato', 'onion', 'leaf', 'fruit', 'lettuce'],
}


# Function to determine the niche category
def determine_niche(isolation_source):
    if pd.isna(isolation_source):
        return None
    isolation_source_lower = isolation_source.lower()
    for niche, keywords in niche_mapping.items():
        if any(keyword in isolation_source_lower for keyword in keywords):
            return niche
    return None

# Ensure the 'Niche' column exists
if 'Niche_IsolationSource' not in dat_passed.columns:
    dat_passed['Niche_IsolationSource'] = None

# Apply the function to the DataFrame
dat_passed.loc[:, 'Niche_IsolationSource'] = dat_passed.apply(
    lambda row: determine_niche(row['Assembly.BioSample.Isolation.source']) 
                if pd.isna(row['Niche_IsolationSource']) 
                else row['Niche_IsolationSource'], axis=1)

# Copy 'Niche' to a new column 'Niche_HostDisease'
dat_passed.loc[:, 'Niche_HostDisease'] = dat_passed['Niche_IsolationSource']

# Overwrite the current value, but not if the new value will become NaN or 'None'
dat_passed.loc[:, 'Niche_HostDisease'] = dat_passed.apply(
    lambda row: (new_value := determine_niche(row['Assembly.BioSample.Host.disease'])) 
                if (new_value := determine_niche(row['Assembly.BioSample.Host.disease'])) not in [None, 'None']
                else row['Niche_HostDisease'], axis=1)

# Copy 'Niche_HostDisease' to a new column 'Niche', now that we know its a better representation of the niche
dat_passed.loc[:, 'Niche'] = dat_passed['Niche_HostDisease'] 

# ============================================================================
# Add Manually Curated Pediatric CF and Bronchiectasis Data
# ============================================================================

manualDat_csv = "/00_data/EarlyCFAndBronchiectasisMetadata.csv"

manualDat = pd.read_csv(manualDat_csv)

# Drop the redundant 'Nucc_Acc' column from the merged dataframe
merged_data.drop(columns=['Nucc_Acc'], inplace=True)

# Move the 'Niche' column to the second position
cols = list(merged_data.columns)
cols.insert(1, cols.pop(cols.index('Niche')))
merged_data = merged_data[cols]

# Function to update 'Niche' based on 'Acute/non-acute' values
def update_niche(df):
    conditions = [
        'Acute(16-20)', 'Acute(6-10)', 'Acute(11-15)', 'Acute(0-5)', 'Child', 'Children/young individuals (early)', 'Early'
    ]
    df.loc[df['Acute/non-acute'].isin(conditions), 'Niche'] = 'early.CF'
    return df

# Apply the function
merged_data = update_niche(merged_data)

# Filter out rows where the 'Niche' column is None or the value itself is "None"
final_data = merged_data[merged_data['Niche'].notna() & (merged_data['Niche'] != 'None')]

# ============================================================================
# Output Final List of Genomes to Download with Metadata Table
# ============================================================================

# Export the final_data dataframe as a .csv file
final_data.to_csv('rawMetadata.csv', index=False)

# Export the Assembly.Accession column as a .txt file with each accession on a new line
final_data['Assembly.Accession'].to_csv('assemblyAccessions.txt', index=False, header=False)