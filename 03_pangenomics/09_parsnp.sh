#!/bin/bash

# Path to the reference genome
referenceGenome="/00_data/referenceGenomes/PAO1_GCA_000006765.1_ASM676v1_genomic.fna"

# Path to the directory containing the query genomes
queryGenomesPath="/00_data/genomes/fna_files/filt_genomes"

# Path to the output directory
outputDirectory="00_data"

# Create a variable to store the list of files
genomeFiles=""

if [[ -n "$genomeListFile" ]]; then
    # Read the list of genome base names from the separate .txt file
    mapfile -t genomeList < "$genomeListFile"

    # Loop through each file in the queryGenomesPath directory
    for file in "$queryGenomesPath"/*; do
        # Get the base name of the file
        baseName=$(basename "$file")
        # Check if the base name is in the genome list
        if [[ " ${genomeList[@]} " =~ " ${baseName} " ]]; then
            genomeFiles+="$file "
        fi
    done
else
    # Use all genomes in the queryGenomesPath directory
    for file in "$queryGenomesPath"/*; do
        genomeFiles+="$file "
    done
fi

# Add the hardcoded outgroup genome to the list of genome files
outgroupGenome="/00_data/referenceGenomes/PA7_GCA_000017205.1_ASM1720v1_genomic.fna"
if [[ -f "$outgroupGenome" ]]; then
    genomeFiles+="$outgroupGenome "
else
    echo "Error: Outgroup genome file not found at $outgroupGenome"
    exit 1
fi

# Trim the trailing space
genomeFiles=$(echo $genomeFiles | sed 's/ *$//')

# Create a temporary directory for filtered genomes
tempDir=$(mktemp -d)

# Create symbolic links for the filtered genome files in the temporary directory
for file in $genomeFiles; do
    ln -s "$file" "$tempDir"
done

# Run Parsnp with the temporary directory
parsnp -r "$referenceGenome" -d "$tempDir" -o "$outputDirectory" -p 24

# Clean up the temporary directory
rm -r "$tempDir"