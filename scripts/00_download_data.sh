#!/usr/bin/env bash
# 00_download_data.sh
# Download example RNA-seq FASTQ files into data/raw/
# Edit the SAMPLE_IDS and URLs/SRA IDs for the dataset you choose.

set -euo pipefail

# Directory where FASTQ files will be stored
RAW_DIR="data/raw"
mkdir -p "$RAW_DIR"

echo "Downloading RNA-seq data into: ${RAW_DIR}"
echo

###############################################################################
# OPTION 1: Download FASTQ files directly with wget/curl (e.g., from ENA/GEO)
###############################################################################
# Uncomment and edit this section if you have direct URLs to FASTQ files.
#
# Example:
#   SAMPLE_IDS=("sample1" "sample2" "sample3" "sample4")
#   FASTQ_URLS_R1=("URL_for_sample1_R1" "URL_for_sample2_R1" ...)
#   FASTQ_URLS_R2=("URL_for_sample1_R2" "URL_for_sample2_R2" ...)
#
# SAMPLE_IDS=("sample1" "sample2")
# FASTQ_URLS_R1=(
#   "https://example.org/path/to/sample1_R1.fastq.gz"
#   "https://example.org/path/to/sample2_R1.fastq.gz"
# )
# FASTQ_URLS_R2=(
#   "https://example.org/path/to/sample1_R2.fastq.gz"
#   "https://example.org/path/to/sample2_R2.fastq.gz"
# )
#
# for i in "${!SAMPLE_IDS[@]}"; do
#   SAMPLE="${SAMPLE_IDS[$i]}"
#   URL_R1="${FASTQ_URLS_R1[$i]}"
#   URL_R2="${FASTQ_URLS_R2[$i]}"
#
#   echo "Downloading ${SAMPLE}..."
#   wget -O "${RAW_DIR}/${SAMPLE}_R1.fastq.gz" "${URL_R1}"
#   wget -O "${RAW_DIR}/${SAMPLE}_R2.fastq.gz" "${URL_R2}"
# done
#
# echo "Download finished (Option 1)."
# exit 0

###############################################################################
# OPTION 2: Download from SRA using the SRA Toolkit (prefetch + fasterq-dump)
###############################################################################
# Uncomment and edit this section if you want to pull data by SRA accessions.
# Requires: prefetch, fasterq-dump (from SRA Toolkit) in your PATH.
###############################################################################

# List of SRA run accessions (replace with real ones)
SRA_IDS=(
  "SRRXXXXXXX"  # replace with real SRA ID
  "SRRYYYYYYY"  # replace with real SRA ID
)

echo "Using SRA Toolkit to download and convert SRA files..."
echo "SRA IDs: ${SRA_IDS[*]}"
echo

for SRA_ID in "${SRA_IDS[@]}"; do
  echo "=== Processing ${SRA_ID} ==="
  # Download .sra file to the default SRA directory
  prefetch "${SRA_ID}"

  # Convert .sra to FASTQ; adjust options if your data is single-end
  fasterq-dump "${SRA_ID}" -O "${RAW_DIR}" --split-files

  # This will produce files like:
  #   data/raw/SRRXXXXXXX_1.fastq
  #   data/raw/SRRXXXXXXX_2.fastq
  #
  # Optionally compress them:
  gzip -f "${RAW_DIR}/${SRA_ID}"_1.fastq
  gzip -f "${RAW_DIR}/${SRA_ID}"_2.fastq

  echo "Finished ${SRA_ID}"
  echo
done

echo "All downloads complete. FASTQ files are in: ${RAW_DIR}"

