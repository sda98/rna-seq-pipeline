#!/usr/bin/env bash
# 02_trim.sh
# Run Trimmomatic on paired-end FASTQ files in data/raw and write trimmed reads to data/trimmed
# Assumes filenames like: sample_X_R1.fastq.gz and sample_X_R2.fastq.gz

set -euo pipefail

#########################
# Configuration
#########################

# Where Trimmomatic .jar is located (override by exporting TRIMMOMATIC_JAR)
TRIMMOMATIC_JAR="${TRIMMOMATIC_JAR:-/usr/share/java/trimmomatic-0.39.jar}"

# Input / output directories (can override with positional args)
RAW_DIR="${1:-data/raw}"
OUT_DIR="${2:-data/trimmed}"

# Number of threads for Trimmomatic
THREADS="${THREADS:-4}"

#########################
# Setup
#########################

if [[ ! -f "${TRIMMOMATIC_JAR}" ]]; then
  echo "ERROR: Trimmomatic jar not found at: ${TRIMMOMATIC_JAR}" >&2
  echo "Set TRIMMOMATIC_JAR to the correct path before running this script." >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

echo "Running Trimmomatic on paired-end reads..."
echo "Input directory : ${RAW_DIR}"
echo "Output directory: ${OUT_DIR}"
echo "Trimmomatic JAR : ${TRIMMOMATIC_JAR}"
echo

#########################
# Main loop
#########################

# Loop over R1 files (supports both .fastq and .fastq.gz)
shopt -s nullglob

for R1 in "${RAW_DIR}"/*_R1.fastq.gz "${RAW_DIR}"/*_R1.fastq; do
  [[ -e "$R1" ]] || continue  # skip if no match

  # Derive R2 path by replacing _R1 with _R2
  R2="${R1/_R1/_R2}"

  if [[ ! -f "$R2" ]]; then
    echo "WARNING: Could not find matching R2 for ${R1}, skipping..." >&2
    continue
  fi

  fname="$(basename "$R1")"

  # Determine sample name (strip _R1.fastq[.gz])
  if [[ "$fname" == *_R1.fastq.gz ]]; then
    sample="${fname%_R1.fastq.gz}"
    ext=".fastq.gz"
  elif [[ "$fname" == *_R1.fastq ]]; then
    sample="${fname%_R1.fastq}"
    ext=".fastq"
  else
    echo "WARNING: Unexpected file name pattern: ${fname}, skipping..." >&2
    continue
  fi

  echo "=== Trimming sample: ${sample} ==="
  echo "R1: $R1"
  echo "R2: $R2"
