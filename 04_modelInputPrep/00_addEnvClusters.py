import pandas as pd
import sys

# Define environment clustering mappings
env_clusters = {
    'Env_Cluster_1': {
        'Lung': 'Other Lung',
        'Bronchiectasis': 'Bronchiectasis',
        'Pneumonia': 'Pneumonia',
        'CF': 'CF',
        'early.CF': 'Early CF',
        'Urinary': 'Urinary',
        'Gastrointestinal': 'Gastrointestinal',
        'Wound': 'Wound',
        'Burn': 'Burn',
        'Blood': 'Blood',
        'Eye': 'Eye',
        'Terrestrial': 'Terrestrial',
        'Aquatic': 'Aquatic',
        'Human_Environment': 'Human Environment',
        'Animal': 'Animal',
        'Ocean': 'Ocean',
        'Built.Environment': 'Built Environment',
        'Hospital.Environment': 'Hospital Environment',
        'Rectal/Feces': 'Rectal/Feces',
        'Waste.Water': 'Waste Water',
    },
    'Env_Cluster_2': {
        'Lung': 'Other Lung',
        'Bronchiectasis': 'Bronchiectasis',
        'Pneumonia': 'Other Clinical',
        'CF': 'CF',
        'early.CF': 'Early CF',
        'Urinary': 'Urinary',
        'Gastrointestinal': 'Gastrointestinal',
        'Wound': 'Wound',
        'Burn': 'Other Clinical',
        'Blood': 'Blood',
        'Eye': 'Other Clinical',
        'Terrestrial': 'Other Environmental',
        'Aquatic': 'Aquatic',
        'Human_Environment': 'Built Environment',
        'Animal': 'Animal',
        'Ocean': 'Other Environmental',
        'Built.Environment': 'Human Environment',
        'Hospital.Environment': 'Hospital Environment',
        'Rectal/Feces': 'Rectal/Feces',
        'Waste.Water': 'Human Environment',
    },
    'Env_Cluster_3': {
        'Lung': 'Other Lung',
        'Bronchiectasis': 'Bronchiectasis',
        'Pneumonia': 'Other Lung',
        'CF': 'CF',
        'early.CF': 'Early CF',
        'Urinary': 'Urinary',
        'Gastrointestinal': 'Gastrointestinal',
        'Wound': 'Wound',
        'Burn': 'Other Clinical',
        'Blood': 'Blood',
        'Eye': 'Other Clinical',
        'Terrestrial': 'Other Environmental',
        'Aquatic': 'Aquatic',
        'Human_Environment': 'Built Environment',
        'Animal': 'Animal',
        'Ocean': 'Aquatic',
        'Built.Environment': 'Human Environment',
        'Hospital.Environment': 'Human Environment',
        'Rectal/Feces': 'Rectal/Feces',
        'Waste.Water': 'Human Environment',
    },
    'Env_Cluster_4': {
        'Lung': 'Other Lung',
        'Bronchiectasis': 'Bronchiectasis',
        'Pneumonia': 'Other Lung',
        'CF': 'CF',
        'early.CF': 'Early CF',
        'Urinary': 'Urinary',
        'Gastrointestinal': 'Gastrointestinal',
        'Wound': 'Wound',
        'Blood': 'Blood',
        'Aquatic': 'Aquatic',
        'Ocean': 'Aquatic',
        'Hospital.Environment': 'Hospital Environment',
        'Rectal/Feces': 'Rectal/Feces',
    },
    'Env_Cluster_5': {
        'Bronchiectasis': 'Bronchiectasis',
        'Pneumonia': 'Pneumonia',
        'CF': 'CF',
        'early.CF': 'Early CF',
        'Urinary': 'Urinary',
        'Gastrointestinal': 'Gastrointestinal',
        'Wound': 'Wound',
        'Blood': 'Blood',
        'Aquatic': 'Aquatic',
        'Ocean': 'Aquatic',
        'Hospital.Environment': 'Hospital Environment',
        'Rectal/Feces': 'Rectal/Feces',
    },
    'Env_Cluster_6': {
        'Bronchiectasis': 'Bronchiectasis',
        'CF': 'CF',
        'early.CF': 'Early CF',
        'Urinary': 'Urinary',
        'Gastrointestinal': 'Gastrointestinal',
        'Wound': 'Wound',
        'Blood': 'Blood',
        'Aquatic': 'Aquatic',
        'Ocean': 'Aquatic',
        'Hospital.Environment': 'Hospital Environment',
        'Rectal/Feces': 'Rectal/Feces',
    },
}

input_file = sys.argv[1] if len(sys.argv) > 1 else "00_data/metadata.csv"
output_file = sys.argv[2] if len(sys.argv) > 2 else "00_data/envBlockedMetadata.csv"

metadata = pd.read_csv(input_file)

for cluster_name, mapping in env_clusters.items():
    metadata[cluster_name] = metadata['Niche'].map(mapping)

metadata.to_csv(output_file, index=False)

print(f"Environment clusters added successfully")
print(f"Columns added: {', '.join(env_clusters.keys())}")
