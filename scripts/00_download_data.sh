#!/usr/bin/env bash
# 00_download_data.sh
# Download example RNA-seq FASTQ files into data/raw/ using SRA Toolkit

set -euo pipefail

# Directory where FASTQ files will be stored
RAW_DIR="data/raw"
mkdir -p "$RAW_DIR"

echo "Downloading RNA-seq data into: ${RAW_DIR}"
echo

# List of SRA run accessions (replace with your real IDs)
SRA_IDS=(
  "SRRWTWTWTWT"
  "SRRWTWTWTWT"
  "SRRWTWTWTWT"
  "SRRWTWTWTWT"
  "SRRKOKOKOKO"
  "SRRKOKOKOKO"
  "SRRKOKOKOKO"
  "SRRKOKOKOKO"
)

echo "Using SRA Toolkit to download and convert SRA files..."
echo "SRA IDs: ${SRA_IDS[*]}"
echo

for SRA_ID in "${SRA_IDS[@]}"; do
  echo "=== Processing ${SRA_ID} ==="

  # Download .sra file to the default SRA directory
  prefetch "${SRA_ID}"

  # Convert .sra to FASTQ; --split-files for paired-end
  fasterq-dump "${SRA_ID}" -O "${RAW_DIR}" --split-files

  # Compress FASTQ files
  if [[ -f "${RAW_DIR}/${SRA_ID}_1.fastq" ]]; then
    gzip -f "${RAW_DIR}/${SRA_ID}_1.fastq"
  fi
  if [[ -f "${RAW_DIR}/${SRA_ID}_2.fastq" ]]; then
    gzip -f "${RAW_DIR}/${SRA_ID}_2.fastq"
  fi

  echo "Finished ${SRA_ID}"
  echo
done

echo "All downloads complete. FASTQ files are in: ${RAW_DIR}"
