#!/usr/bin/env bash
# 03_align.sh
# Align paired-end reads with HISAT2 and produce sorted BAM files
# Default input: data/trimmed (from 02_trim.sh)
# Default output: results/alignments
#
# Assumes filenames like: sample_X_R1.fastq.gz and sample_X_R2.fastq.gz
# Requires:
#   - hisat2
#   - samtools
#   - a pre-built HISAT2 index (set via HISAT2_INDEX_PREFIX or argument)

set -euo pipefail

#########################
# Configuration
#########################

# HISAT2 index prefix (without .1.ht2 etc.)
# You can set this via environment variable or as the third argument.
HISAT2_INDEX_PREFIX="${HISAT2_INDEX_PREFIX:-${3:-/path/to/hisat2/index/prefix}}"

# Input / output directories (can override with positional args)
IN_DIR="${1:-data/trimmed}"
OUT_DIR="${2:-results/alignments}"

# Number of threads
THREADS="${THREADS:-4}"

#########################
# Checks
#########################

if ! command -v hisat2 >/dev/null 2>&1; then
  echo "ERROR: hisat2 not found in PATH." >&2
  exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools not found in PATH." >&2
  exit 1
fi

if [[ "${HISAT2_INDEX_PREFIX}" == "/path/to/hisat2/index/prefix" ]]; then
  echo "ERROR: HISAT2 index prefix is not set." >&2
  echo "Set HISAT2_INDEX_PREFIX env var or pass as 3rd argument, e.g.:" >&2
  echo "  bash scripts/03_align.sh data/trimmed results/alignments /path/to/index/prefix" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

echo "Running HISAT2 alignment..."
echo "Input directory     : ${IN_DIR}"
echo "Output directory    : ${OUT_DIR}"
echo "HISAT2 index prefix : ${HISAT2_INDEX_PREFIX}"
echo "Threads             : ${THREADS}"
echo

#########################
# Main loop
#########################

shopt -s nullglob

# Look for paired-end files with _R1
for R1 in "${IN_DIR}"/*_R1.fastq.gz "${IN_DIR}"/*_R1.fastq; do
  [[ -e "$R1" ]] || continue  # skip if no match

  # Derive R2 by replacing _R1 with _R2
  R2="${R1/_R1/_R2}"

  if [[ ! -f "$R2" ]]; then
    echo "WARNING: No matching R2 for ${R1}, skipping..." >&2
    continue
  fi

  fname="$(basename "$R1")"

  # Determine sample name (strip _R1.fastq[.gz])
  if [[ "$fname" == *_R1.fastq.gz ]]; then
    sample="${fname%_R1.fastq.gz}"
  elif [[ "$fname" == *_R1.fastq ]]; then
    sample="${fname%_R1.fastq}"
  else
    echo "WARNING: Unexpected filename pattern: ${fname}, skipping..." >&2
    continue
  fi

  OUT_BAM="${OUT_DIR}/${sample}.hisat2.sorted.bam"

  echo "=== Aligning sample: ${sample} ==="
  echo "R1: $R1"
  echo "R2: $R2"
  echo "Output BAM: $OUT_BAM"

  # HISAT2 alignment piped directly into samtools sort (coordinate-sorted BAM)
  hisat2 \
    -x "${HISAT2_INDEX_PREFIX}" \
    -1 "${R1}" \
    -2 "${R2}" \
    -p "${THREADS}" \
  | samtools sort -@ "${THREADS}" -o "${OUT_BAM}"

  # Index BAM
  samtools index "${OUT_BAM}"

  echo "Finished sample: ${sample}"
  echo
done

echo "All alignments complete. Sorted BAM + index files are in: ${OUT_DIR}"

