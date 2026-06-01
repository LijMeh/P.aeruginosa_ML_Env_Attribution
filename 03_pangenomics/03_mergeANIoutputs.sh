#!/bin/bash
# ============================================================================
# Merge fastANI Outputs
# ============================================================================
# This script merges fastANI output files from a new run with previous output
# so that the complete set of genomes are compared to each other.
# ============================================================================

# ============================================================================
# Set Paths
# ============================================================================

# Path to the fastANI output directory
fastANI_output="/00_data/fastANI_output"

# ============================================================================
# Merge ANI Outputs
# ============================================================================

cd "$fastANI_output"

# Remove any existing merged output before regenerating
if [ -f merged_fastANI_output.txt ]; then
    rm merged_fastANI_output.txt
fi

# Concatenate all per-batch output files into a single merged file
cat output_*.txt > merged_fastANI_output.txt

python3 <<EOF
import os
from multiprocessing import Pool, cpu_count

# Function to extract the base name without extension
def get_base_name(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

# Process a chunk of lines
def process_chunk(lines):
    processed_lines = []
    for line in lines:
        genome1, genome2, ani, fragment_matches, total_fragments = line.strip().split()
        base_genome1 = get_base_name(genome1)
        base_genome2 = get_base_name(genome2)
        processed_lines.append(f"{base_genome1},{base_genome2},{ani},{fragment_matches},{total_fragments}\n")
    return processed_lines

# Process the merged_fastANI_output.txt to extract base names using multiple cores
def process_large_file_multicore(input_file, output_file, chunk_size=1024*1024, num_cores=cpu_count()):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        outfile.write('Genome1,Genome2,ANI,Fragment_Matches,Total_Fragments\n')
        pool = Pool(num_cores)
        while True:
            lines = infile.readlines(chunk_size)
            if not lines:
                break
            processed_chunks = pool.map(process_chunk, [lines])
            for chunk in processed_chunks:
                outfile.writelines(chunk)
        pool.close()
        pool.join()

num_cores = 24

process_large_file_multicore('merged_fastANI_output.txt', 'processed_merged_fastANI_output.csv', num_cores=num_cores)
EOF

# ============================================================================
# Quick Sanity Check
# ============================================================================

# Count rows in the processed output as a completeness check
processed_merged_ANI_output="$fastANI_output/processed_merged_fastANI_output.csv"

line_count=$(wc -l < "$processed_merged_ANI_output")
echo "Total number of rows: $line_count"