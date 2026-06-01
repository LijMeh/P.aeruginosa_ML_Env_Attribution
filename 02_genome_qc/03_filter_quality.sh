#!/bin/bash
# ============================================================================
# Script filters the downloaded genomes based on quality metrics from checkM
# ============================================================================
# This script takes the CheckM output and filters out genomes that do not 
# meet the quality criteria.
# ============================================================================

# ============================================================================
# Set Paths 
# ============================================================================

input_path="/00_data/genomes/fna_files/"
unfilt_directory="/00_data/genomes/fna_files/unfilt_genomes"
target_directory="/00_data/genomes/fna_files/filt_genomes"
tsv_file="/00_data/merged_checkM_output_Qualityfiltered.tsv"

# ============================================================================
# Merge and Filter CheckM Outputs
# ============================================================================
# This should work regardless of if you batched checkM or not.

# Initialize an empty file list
file_list=()

# Loop through each directory in the given folder
for dir in "$input_path"*/; do
    # Check if "output.tsv" exists in the directory
    if [ -f "${dir}output.tsv" ]; then
        # Add the path to the file list
        file_list+=("${dir}output.tsv")
    fi
done

cd $input_path

# Initialize variables
merged_file="merged_checkM_output.tsv"
filtered_file="${merged_file%.tsv}_Qualityfiltered.tsv"
first_file=true

# Create a temporary file to hold the merged output
> "$merged_file"

# Check if there are multiple files or a single file
if [ "${#file_list[@]}" -eq 1 ]; then
    # Only one file provided, use it directly
    input_file="${file_list[0]}"
    trimmed_data=$(sed '1,7d;$d' "$input_file")
    echo "$trimmed_data" > "$merged_file"
    echo "Processed single file $input_file"
else
    # Multiple files provided, merge them
    for input_file in "${file_list[@]}"; do
        # Trim first 7 lines and last line, then merge
        trimmed_data=$(sed '1,7d;$d' "$input_file")
        
        if $first_file; then
            # Include header for the first file
            echo "$trimmed_data" >> "$merged_file"
            first_file=false
        else
            # Skip header for subsequent files and append data
            echo "$trimmed_data" | tail -n +2 >> "$merged_file"
        fi
        
        echo "Processed $input_file"
    done
    echo "All files merged into $merged_file"
fi

# Extract header and determine column indices
header=$(head -n 1 "$merged_file")
IFS=$'\t' read -ra columns <<< "$header"

declare -A col_map=(["Completeness"]=-1 ["Contamination"]=-1 ["N50 (contigs)"]=-1 ["# contigs"]=-1)

for i in "${!columns[@]}"; do
    case "${columns[$i]}" in
        "Completeness") col_map["Completeness"]=$((i+1)) ;;
        "Contamination") col_map["Contamination"]=$((i+1)) ;;
        "N50 (contigs)") col_map["N50 (contigs)"]=$((i+1)) ;;
        "# contigs") col_map["# contigs"]=$((i+1)) ;;
    esac
done

# Validate column detection
for col in "${!col_map[@]}"; do
    if [[ ${col_map[$col]} -eq -1 ]]; then
        echo "Error: Missing required column '$col'"
        exit 1
    fi
done

echo "Column indices: Completeness=${col_map[Completeness]}, Contamination=${col_map[Contamination]}, N50=${col_map[N50 (contigs)]}, # contigs=${col_map[# contigs]}"

# Filter rows based on criteria
lines_kept=0
lines_removed=0

{
    read -r header_line # Read and write header to output file
    echo "$header_line" > "$filtered_file"

    while IFS=$'\t' read -r -a row; do
        completeness=${row[${col_map[Completeness]}-1]}
        contamination=${row[${col_map[Contamination]}-1]}
        n50=${row[${col_map[N50 (contigs)]}-1]}
        contigs=${row[${col_map[# contigs]}-1]}

        if (( $(echo "$completeness >= 98" | bc -l) )) &&
           (( $(echo "$contamination <= 2" | bc -l) )) &&
           (( $(echo "$n50 >= 100000" | bc -l) )) &&
           (( $(echo "$contigs <= 300" | bc -l) )); then
            # Row meets criteria, keep it
            echo "${row[*]}" | tr ' ' '\t' >> "$filtered_file"
            ((lines_kept++))
        else
            # Row does not meet criteria, log it to stderr
            echo "Filtered out: Completeness=$completeness, Contamination=$contamination, N50=$n50, # contigs=$contigs" >&2
            ((lines_removed++))
        fi
    done
} < <(tail -n +2 "$merged_file")

echo "Lines kept: $lines_kept"
echo "Lines removed: $lines_removed"
echo "Filtered file saved as $filtered_file"


# ============================================================================
# Relocate Genomes That Passed Filter
# ============================================================================

# Extract the first column from the provided TSV file and use it as the input list file
if [ ! -f "$tsv_file" ]; then
    echo "Error: TSV file '$tsv_file' not found."
    exit 1
fi

# Extract the first column and save it to the list file
list_file="$(pwd)/genome_list.txt"
cut -f1 "$tsv_file" > "$list_file"
echo "Extracted first column from $tsv_file and saved to $list_file"

# Check if the list file exists
if [ ! -f "$list_file" ]; then
    echo "Error: List file '$list_file' not found."
    exit 1
fi

# Create the target directory if it doesn't exist
mkdir -p "$target_directory"

# Read the list file and move the corresponding files
while IFS= read -r filename; do
    # Find files with the given name (without extension) and move them
    cd $unfilt_directory
    for file in "$filename".fna; do
        if [ -e "$file" ]; then
            cp "$file" "$target_directory"
            echo "Moved $file to $target_directory"
        else
            echo "File $file not found."
        fi
    done
done < "$list_file"

echo "All specified files have been moved."