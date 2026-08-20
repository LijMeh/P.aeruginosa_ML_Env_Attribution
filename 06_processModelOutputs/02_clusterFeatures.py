"""
Cluster protein sequences for genes evaluated by the blocked models.

Gathers the full set of genes that trainBlocks.py computed SHAP values for
(one wide-format CSV per genomovar holdout x niche in
00_data/model_output/SHAP_Values, genes as columns), pulls their protein
sequences from a reference FASTA, and clusters them with CD-HIT so that
sequence-similar/paralogous genes can be treated as one signal downstream.

Writes 00_data/feature_clusters/clusters.fasta and clusters.fasta.clstr,
which 03_processShapBestModel.R loads via parse_cdhit_clstr_gene_map() to
propagate study-bias flags across cluster mates.
"""

import argparse
import glob
import os
import subprocess

from Bio import SeqIO


def gather_genes_of_interest(shap_dir):
    """Union of gene names (column headers) across all per-niche SHAP CSVs."""
    shap_files = glob.glob(os.path.join(shap_dir, "*.csv"))
    if not shap_files:
        raise FileNotFoundError(f"No SHAP CSVs found in {shap_dir}")

    genes = set()
    for f in shap_files:
        with open(f) as fh:
            header = fh.readline().strip().split(",")
        # First column is the Genome_ID index, the rest are gene names
        genes.update(header[1:])

    print(f"Found {len(shap_files)} SHAP files, {len(genes)} unique genes")
    return sorted(genes)


def extract_gene_sequences(reference_fasta, genes, output_fasta):
    """Pull the protein sequences for `genes` out of `reference_fasta`."""
    genes = set(genes)
    written = 0

    with open(output_fasta, "w") as out_fh:
        for record in SeqIO.parse(reference_fasta, "fasta"):
            if record.id in genes:
                SeqIO.write(record, out_fh, "fasta")
                written += 1

    print(f"Wrote {written}/{len(genes)} requested gene sequences to {output_fasta}")
    if written < len(genes):
        print(f"WARNING: {len(genes) - written} genes had no match in {reference_fasta}")

    return written


def run_cdhit(input_fasta, output_prefix, identity, threads, memory, word_size):
    """Run CD-HIT to cluster protein sequences at ortholog-level identity."""
    output_fasta = f"{output_prefix}.fasta"
    cluster_file = f"{output_prefix}.fasta.clstr"

    cmd = [
        "cd-hit",
        "-i", input_fasta,
        "-o", output_fasta,
        "-c", str(identity),
        "-n", str(word_size),
        "-T", str(threads),
        "-M", str(memory),
        "-d", "0",     # full length description
        "-g", "1",     # slower but more accurate clustering
        "-aL", "0.0",  # no alignment length requirement
        "-aS", "0.0",  # no alignment coverage requirement for the shorter sequence
    ]

    print(f"Running CD-HIT at identity threshold {identity}: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    return output_fasta, cluster_file


def summarize_clusters(cluster_file):
    """Print basic cluster count/size stats for a sanity check."""
    clusters = {}
    current_cluster = None

    with open(cluster_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_cluster = int(line.split()[1])
                clusters[current_cluster] = []
            elif line:
                seq_id = line.split(">")[1].split("...")[0]
                clusters[current_cluster].append(seq_id)

    sizes = [len(v) for v in clusters.values()]
    print(f"Total clusters: {len(clusters)}")
    print(f"Average cluster size: {sum(sizes) / len(clusters):.2f}")
    print(f"Largest cluster size: {max(sizes)}")


def word_size_for_identity(identity):
    # CD-HIT protein word-size recommendations: 5 (>=0.7), 4 (>=0.6), 3 (>=0.5), 2 (<0.5)
    if identity >= 0.7:
        return 5
    if identity >= 0.6:
        return 4
    if identity >= 0.5:
        return 3
    return 2


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Cluster the protein sequences of SHAP-evaluated genes with CD-HIT."
    )
    parser.add_argument(
        "--reference-fasta", required=True,
        help="Protein FASTA with one sequence per gene, headers matching the model's gene names"
    )
    parser.add_argument(
        "--shap-dir", default="00_data/model_output/SHAP_Values",
        help="Directory of per-niche SHAP CSVs from trainBlocks.py (default: %(default)s)"
    )
    parser.add_argument(
        "--output-prefix", default="00_data/feature_clusters/clusters",
        help="Output prefix; writes <prefix>.fasta and <prefix>.fasta.clstr (default: %(default)s)"
    )
    parser.add_argument("--identity", type=float, default=0.70, help="CD-HIT sequence identity threshold (default: %(default)s)")
    parser.add_argument("--threads", type=int, default=50, help="CPU threads for CD-HIT (default: %(default)s)")
    parser.add_argument("--memory", type=int, default=100000, help="Memory limit in MB for CD-HIT (default: %(default)s)")
    parser.add_argument("--word-size", type=int, default=None, help="CD-HIT word size (default: auto from --identity)")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_prefix), exist_ok=True)

    genes = gather_genes_of_interest(args.shap_dir)

    genes_fasta = os.path.join(os.path.dirname(args.output_prefix), "genes_of_interest.fasta")
    extract_gene_sequences(args.reference_fasta, genes, genes_fasta)

    word_size = args.word_size if args.word_size is not None else word_size_for_identity(args.identity)

    output_fasta, cluster_file = run_cdhit(
        genes_fasta, args.output_prefix, args.identity, args.threads, args.memory, word_size
    )

    summarize_clusters(cluster_file)
    print(f"Representative sequences: {output_fasta}")
    print(f"Cluster file: {cluster_file}")
