#!/bin/bash
# ============================================================================
# ANI Clustering
# ============================================================================
# This script takes the fastANI output and identifies clusters at multiple
# ANI thresholds ranging from 99.5% to 99.99%.
# ============================================================================

# ============================================================================
# Set Paths
# ============================================================================

# Input fastANI output file
fastani_output_file="/00_data/fastANI_output/processed_merged_fastANI_output.csv"

# Output path for ANI clustering results
filtered_fastani_output_path="/00_data/fastANI_output/ANI_clusters"

# ============================================================================
# Identify Clusters Across ANI Thresholds
# ============================================================================

mkdir -p "$filtered_fastani_output_path"

# Number of parallel jobs
num_jobs=15
if [ "$num_jobs" -lt 1 ]; then num_jobs=1; fi

# Function to run clusterANI.py for a given threshold
run_cluster() {
    ani="$1"
    ani_rounded=$(printf "%.2f" "$ani")
    out_dir="${filtered_fastani_output_path}/ANI_${ani_rounded}"
    mkdir -p "$out_dir"
    echo "Clustering at ANI threshold: $ani_rounded%"
    python 00_clusterANI.py -i "$fastani_output_file" -o "$out_dir" -t "$ani_rounded" -p 10
}

export -f run_cluster
export fastani_output_file
export filtered_fastani_output_path

# Specify range and number of steps
ani_start=99.5
ani_end=99.99
ani_steps=15

# Generate the threshold sequence using awk
seq_list=$(awk -v start="$ani_start" -v end="$ani_end" -v steps="$ani_steps" \
    'BEGIN{for(i=0;i<steps;i++) printf "%.8f\n", start+i*(end-start)/(steps-1)}')

# Run clustering in parallel across all thresholds
echo "$seq_list" | xargs -P "$num_jobs" -I {} bash -c 'run_cluster "$@"' _ {}