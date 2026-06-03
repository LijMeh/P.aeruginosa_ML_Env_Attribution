#!/bin/bash
# ============================================================================
# Script to run CheckM on downloaded genomes for P. aeruginosa
# ============================================================================
# This script runs CheckM on the downloaded genome sequences to assess their quality.
# This is not USUALLY strictly neccessary, as most genomes on NCBI currently have
# checkM quality and completeness data already available.
# ============================================================================

# ============================================================================
# Set Paths 
# ============================================================================

# Directory containing the folders with .fna files
source_directory="/00_data/genomes/ncbi_dataset/data"
# Directory to copy the .fna files to
destination_directory="/00_data/genomes/fna_files"

# ============================================================================
# Relocate Genomes
# ============================================================================

# Create the destination directory if it doesn't exist
mkdir -p "$destination_directory/unfilt_genomes"

# Find and copy all .fna files from the source directory to the destination directory
find "$source_directory" -type f -name "*.fna" -exec cp {} "$destination_directory/unfilt_genomes" \;

# Count the number of folders in the source directory
num_folders=$(find "$source_directory" -mindepth 1 -maxdepth 1 -type d | wc -l)

echo "Number of folders in the source directory: $num_folders"

# Count the number of .fna files copied to the destination directory
num_files=$(find "$destination_directory/unfilt_genomes" -type f -name "*.fna" | wc -l)

echo "Number of .fna files copied to the destination directory: $num_files"
echo "All .fna files have been copied to $destination_directory/unfilt_genomes"

# ============================================================================
# Rename Genomes to Consistent Standard
# ============================================================================

echo "Running run_rename_genomes"

# Directory containing the files
binpath="$destination_directory/unfilt_genomes"

# Check if any .fna files exist
if [ ! "$(ls -A "$binpath"/*.fna 2>/dev/null)" ]; then
    echo "No .fna files found in "$binpath""
fi

# Loop through all .fna files in the directory
for file in "$binpath"/*.fna; do
    if [ -f "$file" ]; then
        # Extract the base name of the file
        basename=$(basename "$file")
        
        # Extract the GCA part and version, replacing . with _
        gca_part=$(echo "$basename" | sed -E 's/(GCA_[0-9]+\.[0-9]+).*/\1/' | tr '.' '_')
        
        # Rename the file
        mv "$file" "$binpath/$gca_part.fna"
        
        echo "Renamed $basename to $gca_part.fna"
    fi
done

# ============================================================================
# Run checkM 
# ============================================================================
# This assumes there are a "reasonable" number of genomes to process, if you're 
# running 10K+ genomes you want to batch the files.  
# Make sure your checkM reference data is exported as a variable `export CHECKM_DATA_PATH=/00_data/checkM`

reference_file="/00_data/pAeruginosa.ms"

cd $destination_directory

export CHECKM_DATA_PATH="/00_resources/checkM"

binpath="$destination_directory/unfilt_genomes"
output="$destination_directory/CheckM_Output"

mkdir -p "$output"

# Get the total number of available cores
total_cores=$(nproc)

checkm analyze $reference_file -t "$total_cores" -x fna "$binpath" "$output"

checkm qa $reference_file -t "$total_cores" -o 2 --tab_table "$output" > "${output}/output.tsv"

echo "CheckM analysis complete"

# Count the number of folders in the CheckM_Output directory
num_output_folders=$(find "$output/bins" -mindepth 1 -maxdepth 1 -type d | wc -l)
echo "Number of folders in the CheckM_Output directory: $num_output_folders"