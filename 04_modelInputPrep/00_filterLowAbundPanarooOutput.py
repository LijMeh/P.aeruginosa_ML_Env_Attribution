import pandas as pd
import numpy as np
import sys
from concurrent.futures import ThreadPoolExecutor

def process_chunk(chunk, cutoff, total_cols):
    print(f"Processing chunk with {len(chunk)} rows...")
    # Exclude first column (gene name) for sum
    sums = chunk.iloc[:, 1:].sum(axis=1)
    percentages = (sums / total_cols)
    keep_mask = percentages >= cutoff
    filtered = chunk[keep_mask]
    removed = chunk[~keep_mask]
    removed_info = pd.DataFrame({
        'Gene': removed.iloc[:, 0],
        'Percentage': (percentages[~keep_mask] * 100)
    })
    return filtered, removed_info

def main(input_file, output_path, cutoff, threads=4, chunksize=10000):
    # Get total number of columns (excluding gene name)
    print(f"Loading header from {input_file} to determine total columns...")
    df_head = pd.read_csv(input_file, sep='\t', nrows=1)
    total_cols = len(df_head.columns) - 1
    print(f"Total columns (excluding gene name): {total_cols}")

    filtered_chunks = []
    removed_chunks = []

    print(f"Starting chunked processing with {threads} threads and chunk size {chunksize}...")
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        for i, chunk in enumerate(pd.read_csv(input_file, sep='\t', chunksize=chunksize)):
            print(f"Submitting chunk {i+1} to thread pool...")
            futures.append(executor.submit(process_chunk, chunk, cutoff, total_cols))
        for i, future in enumerate(futures):
            filtered, removed_info = future.result()
            print(f"Chunk {i+1} processed: {len(filtered)} kept, {len(removed_info)} removed.")
            filtered_chunks.append(filtered)
            removed_chunks.append(removed_info)

    print("Concatenating filtered and removed chunks...")
    filtered_df = pd.concat(filtered_chunks)
    removed_df = pd.concat(removed_chunks)

    print(f"Total rows kept: {len(filtered_df)}")
    print(f"Total rows removed: {len(removed_df)}")

    # Save outputs
    percent_cutoff = cutoff * 100
    cutoff_str = str(percent_cutoff).rstrip('0').rstrip('.') if '.' in str(percent_cutoff) else str(int(percent_cutoff))
    output_file = os.path.join(output_path, f"gene_presence_absence_cutoff_{cutoff_str}%.Rtab")
    removed_file = os.path.join(output_path, f"gene_presence_absence_cutoff_{cutoff_str}%_removed_genes.csv")
    print(f"Saving filtered output to {output_file}...")
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print(f"Saving removed genes info to {removed_file}...")
    removed_df.to_csv(removed_file, sep='\t', index=False)

if __name__ == "__main__":
    import os
    print("Starting filter_rtab.py...")
    if len(sys.argv) != 5:
        print("Usage: python filter_rtab.py <input_rtab> <cutoff_fraction> <threads> <output_path>")
        print("cutoff_fraction should be between 0 and 1 (e.g., 0.1 for 10%)")
        sys.exit(1)
    input_file = sys.argv[1]
    cutoff = float(sys.argv[2])
    threads = int(sys.argv[3])
    output_path = sys.argv[4]
    print(f"Input file: {input_file}")
    print(f"Cutoff fraction: {cutoff}")
    print(f"Threads: {threads}")
    print(f"Output path: {output_path}")
    if not os.path.exists(output_path):
        print(f"Output path {output_path} does not exist. Creating it...")
        os.makedirs(output_path, exist_ok=True)
    main(input_file, output_path, cutoff, threads=threads)