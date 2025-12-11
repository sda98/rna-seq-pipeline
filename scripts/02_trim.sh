#!/usr/bin/env bash
# 02_trim.sh
# Run Trimmomatic on FASTQ files in data/raw and write trimmed reads to data/trimmed
# Handles both paired-end and single-end reads.
# Paired naming patterns supported:
#   sample_R1.fastq(.gz) <-> sample_R2.fastq(.gz)
#   sample_1.fastq(.gz)  <-> sample_2.fastq(.gz)

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

echo "Running Trimmomatic..."
echo "Input directory : ${RAW_DIR}"
echo "Output directory: ${OUT_DIR}"
echo "Trimmomatic JAR : ${TRIMMOMATIC_JAR}"
echo

shopt -s nullglob

# Track which files we've already processed (for single-end pass)
declare -A PROCESSED

#########################
# Pass 1: Paired-end reads
#########################

echo "=== Pass 1: Paired-end trimming ==="

# Candidate R1 patterns
R1_PATTERNS=(
  "*_R1.fastq.gz"
  "*_R1.fastq"
  "*_1.fastq.gz"
  "*_1.fastq"
)

for pattern in "${R1_PATTERNS[@]}"; do
  for R1 in "${RAW_DIR}"/${pattern}; do
    [[ -e "$R1" ]] || continue

    # Skip if we already processed this file
    if [[ -n "${PROCESSED[$R1]:-}" ]]; then
      continue
    fi

    fname="$(basename "$R1")"
    ext=""
    sample=""
    R2=""

    # Determine pattern and mate
    if [[ "$fname" == *_R1.fastq.gz ]]; then
      sample="${fname%_R1.fastq.gz}"
      ext=".fastq.gz"
      R2="${RAW_DIR}/${sample}_R2${ext}"
    elif [[ "$fname" == *_R1.fastq ]]; then
      sample="${fname%_R1.fastq}"
      ext=".fastq"
      R2="${RAW_DIR}/${sample}_R2${ext}"
    elif [[ "$fname" == *_1.fastq.gz ]]; then
      sample="${fname%_1.fastq.gz}"
      ext=".fastq.gz"
      R2="${RAW_DIR}/${sample}_2${ext}"
    elif [[ "$fname" == *_1.fastq ]]; then
      sample="${fname%_1.fastq}"
      ext=".fastq"
      R2="${RAW_DIR}/${sample}_2${ext}"
    else
      # Not a recognized R1 pattern
      continue
    fi

    if [[ ! -f "$R2" ]]; then
      echo "WARNING: Could not find matching R2 for ${R1}, will treat as single-end later." >&2
      continue
    fi

    echo "=== Trimming paired sample: ${sample} ==="
    echo "R1: $R1"
    echo "R2: $R2"

    # Output filenames
    out_R1_paired="${OUT_DIR}/${sample}_R1.trimmed.paired${ext}"
    out_R1_unpaired="${OUT_DIR}/${sample}_R1.trimmed.unpaired${ext}"
    out_R2_paired="${OUT_DIR}/${sample}_R2.trimmed.paired${ext}"
    out_R2_unpaired="${OUT_DIR}/${sample}_R2.trimmed.unpaired${ext}"

    # Run Trimmomatic in paired-end mode
    java -jar "${TRIMMOMATIC_JAR}" PE \
      -threads "${THREADS}" \
      "$R1" "$R2" \
      "$out_R1_paired" "$out_R1_unpaired" \
      "$out_R2_paired" "$out_R2_unpaired" \
      LEADING:26 TRAILING:26 SLIDINGWINDOW:4:26 MINLEN:140

    PROCESSED["$R1"]=1
    PROCESSED["$R2"]=1

    echo
  done
done

#########################
# Pass 2: Single-end reads
#########################

echo "=== Pass 2: Single-end trimming ==="

for fq in "${RAW_DIR}"/*.fastq "${RAW_DIR}"/*.fastq.gz; do
  [[ -e "$fq" ]] || continue

  # Skip if already processed as part of a pair
  if [[ -n "${PROCESSED[$fq]:-}" ]]; then
    continue
  fi

  fname="$(basename "$fq")"
  ext=""
  sample=""

  if [[ "$fname" == *.fastq.gz ]]; then
    sample="${fname%.fastq.gz}"
    ext=".fastq.gz"
  elif [[ "$fname" == *.fastq ]]; then
    sample="${fname%.fastq}"
    ext=".fastq"
  else
    echo "WARNING: Unexpected file name pattern for single-end: ${fname}, skipping..." >&2
    continue
  fi

  echo "=== Trimming single-end sample: ${sample} ==="
  echo "Input: $fq"

  out_SE="${OUT_DIR}/${sample}.trimmed${ext}"

  java -jar "${TRIMMOMATIC_JAR}" SE \
    -threads "${THREADS}" \
    "$fq" \
    "$out_SE" \
    LEADING:26 TRAILING:26 SLIDINGWINDOW:4:26 MINLEN:140

  PROCESSED["$fq"]=1

  echo
done

echo "Trimmomatic trimming complete. Trimmed files are in: ${OUT_DIR}"
