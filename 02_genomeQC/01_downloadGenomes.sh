#!/bin/bash
# ============================================================================
# Script to download genomes for P. aeruginosa
# ============================================================================
# This script uses the assembly accessions from the filtered metadata to download
# the corresponding genome sequences from NCBI.
# ============================================================================

# ============================================================================
# Set Paths 
# ============================================================================

input_file="00_data/assemblyAccessions.txt"
download_folder="00_data/genomes"
genomesToDownload="00_data/genomesToDownload.txt"

# ============================================================================
#  Download Genomes
# ============================================================================

mkdir -p $download_folder

datasetsInput=$(<$genomesToDownload)

datasets download genome accession $datasetsInput --dehydrated --filename $download_folder/new_PA.zip

unzip $download_folder/new_PA.zip -d $download_folder/

datasets rehydrate --directory $download_folder/