#!/usr/bin/env bash
# 01_qc.sh
# Run FastQC on all FASTQ files in data/raw and save reports to results/qc

set -euo pipefail

# Directories (can be overridden by passing arguments)
RAW_DIR="${1:-data/raw}"
OUT_DIR="${2:-results/qc}"

mkdir -p "${OUT_DIR}"

echo "Running FastQC..."
echo "Input directory : ${RAW_DIR}"
echo "Output directory: ${OUT_DIR}"
echo

# Run FastQC on all FASTQ/FASTQ.GZ files
for fq in "${RAW_DIR}"/*.fastq "${RAW_DIR}"/*.fastq.gz 2>/dev/null; do
  # Skip if no files match
  [ -e "$fq" ] || continue
  echo "Processing: $(basename "$fq")"
  fastqc "$fq" -o "${OUT_DIR}"
done

echo
echo "FastQC complete. Reports are in: ${OUT_DIR}"

