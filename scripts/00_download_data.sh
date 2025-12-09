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

