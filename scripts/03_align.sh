#!/usr/bin/env bash
# 03_align.sh
# Align trimmed reads with HISAT2 and produce coordinate-sorted BAM files.
#
# Default input : data/trimmed (outputs from 02_trim.sh)
# Default output: results/alignments
#
# Handles:
#   - Paired-end: *_R1.trimmed.paired.fastq[.gz] + *_R2.trimmed.paired.fastq[.gz]
#   - Single-end: *.trimmed.fastq[.gz], *_R1.trimmed.unpaired.fastq[.gz], *_R2.trimmed.unpaired.fastq[.gz]
#
# Reference:
#   - By default, uses mouse GRCm39 (mm39) primary assembly from GENCODE (GRCm39.primary_assembly.genome.fa)
#   - Builds a HISAT2 index under ref/mm39/hisat2/mm39 on first run.
#   - You can override the index path with HISAT2_INDEX_PREFIX or a 3rd argument.
#
# Requirements:
#   - hisat2
#   - hisat2-build
#   - samtools
#   - wget (for auto-download of mm39)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

IN_DIR="${1:-${ROOT_DIR}/data/trimmed}"
OUT_DIR="${2:-${ROOT_DIR}/results/alignments}"
TMP_DIR="${TMP_DIR:-/tmp/rnaseq_samtools}"
mkdir -p "${OUT_DIR}" "${TMP_DIR}"

#########################
# Configuration
#########################

# Input / output directories (can override with positional args)
IN_DIR="${1:-data/trimmed}"
OUT_DIR="${2:-results/alignments}"

# Number of threads
THREADS="${THREADS:-4}"

# Default mm39 reference layout
DEFAULT_REF_DIR="ref/mm39"
DEFAULT_INDEX_PREFIX="${DEFAULT_REF_DIR}/hisat2/mm39"

# HISAT2 index prefix (can override via env or 3rd arg)
HISAT2_INDEX_PREFIX="${HISAT2_INDEX_PREFIX:-${3:-$DEFAULT_INDEX_PREFIX}}"

# GRCm39 genome FASTA URL (Gencode, primary assembly)
MM39_FASTA_URL="${MM39_FASTA_URL:-https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M29/GRCm39.primary_assembly.genome.fa.gz}"

# Where we'll store the genome FASTA if using default mm39 layout
GENOME_FASTA="${DEFAULT_REF_DIR}/GRCm39.primary_assembly.genome.fa"
GENOME_FASTA_GZ="${GENOME_FASTA}.gz"

#########################
# Checks
#########################

if ! command -v hisat2 >/dev/null 2>&1; then
  echo "ERROR: hisat2 not found in PATH." >&2
  exit 1
fi

if ! command -v hisat2-build >/dev/null 2>&1; then
  echo "ERROR: hisat2-build not found in PATH." >&2
  exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools not found in PATH." >&2
  exit 1
fi

# Only need wget if we have to auto-download mm39
if [[ "${HISAT2_INDEX_PREFIX}" == "${DEFAULT_INDEX_PREFIX}" ]]; then
  if ! command -v wget >/dev/null 2>&1; then
    echo "ERROR: wget not found in PATH (needed to auto-download GRCm39)." >&2
    echo "Install it with: sudo apt install wget" >&2
    exit 1
  fi
fi

mkdir -p "${OUT_DIR}"

echo "Running HISAT2 alignment..."
echo "Input directory     : ${IN_DIR}"
echo "Output directory    : ${OUT_DIR}"
echo "HISAT2 index prefix : ${HISAT2_INDEX_PREFIX}"
echo "Threads             : ${THREADS}"
echo

#########################
# Reference setup
#########################

ensure_hisat2_index() {
  local prefix="$1"

  # If index already exists, do nothing
  if [[ -f "${prefix}.1.ht2" || -f "${prefix}.1.ht2l" ]]; then
    echo "Found existing HISAT2 index: ${prefix}"
    echo
    return
  fi

  # If using the default mm39 layout, auto-download and build
  if [[ "$prefix" == "$DEFAULT_INDEX_PREFIX" ]]; then
    echo "HISAT2 index not found for prefix: ${prefix}"
    echo "Preparing GRCm39 (mm39) reference and building index..."
    echo

    mkdir -p "${DEFAULT_REF_DIR}"
    mkdir -p "$(dirname "$prefix")"

    # Download genome FASTA if needed
    if [[ ! -f "${GENOME_FASTA}" ]]; then
      if [[ ! -f "${GENOME_FASTA_GZ}" ]]; then
        echo "Downloading GRCm39 FASTA from:"
        echo "  ${MM39_FASTA_URL}"
        wget -O "${GENOME_FASTA_GZ}" "${MM39_FASTA_URL}"
      fi

      echo "Decompressing genome FASTA..."
      gunzip -f "${GENOME_FASTA_GZ}"
    fi

    echo "Building HISAT2 index (this may take a while)..."
    hisat2-build -p "${THREADS}" "${GENOME_FASTA}" "${prefix}"
    echo "HISAT2 index built at: ${prefix}"
    echo
  else
    echo "ERROR: HISAT2 index files not found for prefix: ${prefix}" >&2
    echo "Please build them manually, e.g.:" >&2
    echo "  hisat2-build -p ${THREADS} genome.fa ${prefix}" >&2
    exit 1
  fi
}

ensure_hisat2_index "${HISAT2_INDEX_PREFIX}"

#########################
# Alignment
#########################

shopt -s nullglob

# Track which FASTQ files have already been used (so we can handle single-end separately)
declare -A ALIGNED

echo "=== Pass 1: Paired-end alignment (trimmed.paired) ==="

# Paired-end R1 patterns from 02_trim.sh
for R1 in "${IN_DIR}"/*_R1.trimmed.paired.fastq.gz "${IN_DIR}"/*_R1.trimmed.paired.fastq; do
  [[ -e "$R1" ]] || continue

  # Find matching R2 by replacing _R1.trimmed.paired with _R2.trimmed.paired
  R2="${R1/_R1.trimmed.paired/_R2.trimmed.paired}"

  if [[ ! -f "$R2" ]]; then
    echo "WARNING: No matching R2 for ${R1}, skipping paired alignment for this file..." >&2
    continue
  fi

  fname="$(basename "$R1")"

  # Derive sample name for BAM prefix
  if [[ "$fname" == *_R1.trimmed.paired.fastq.gz ]]; then
    sample="${fname%_R1.trimmed.paired.fastq.gz}"
  elif [[ "$fname" == *_R1.trimmed.paired.fastq ]]; then
    sample="${fname%_R1.trimmed.paired.fastq}"
  else
    echo "WARNING: Unexpected R1 filename pattern: ${fname}, skipping..." >&2
    continue
  fi

  OUT_BAM="${OUT_DIR}/${sample}.hisat2.sorted.bam"

  echo "=== Aligning paired sample: ${sample} ==="
  echo "R1: $R1"
  echo "R2: $R2"
  echo "Output BAM: $OUT_BAM"

  hisat2 \
    -x "${HISAT2_INDEX_PREFIX}" \
    -1 "${R1}" \
    -2 "${R2}" \
    -p "${THREADS}" \
  | samtools sort -@ "${THREADS}" -o "${OUT_BAM}"

  samtools index "${OUT_BAM}"

  ALIGNED["$R1"]=1
  ALIGNED["$R2"]=1

  echo "Finished paired sample: ${sample}"
  echo
done

echo "=== Pass 2: Single-end alignment (trimmed/unpaired) ==="

# Single-end candidates:
#   - *.trimmed.fastq[.gz]       (true SE input)
#   - *_R1.trimmed.unpaired.*    (lost mate)
#   - *_R2.trimmed.unpaired.*    (lost mate)
for fq in \
  "${IN_DIR}"/*.trimmed.fastq.gz \
  "${IN_DIR}"/*.trimmed.fastq \
  "${IN_DIR}"/*trimmed.unpaired.fastq.gz \
  "${IN_DIR}"/*trimmed.unpaired.fastq; do

  [[ -e "$fq" ]] || continue

  # If this file was already used in a paired alignment, skip it
  if [[ -n "${ALIGNED[$fq]:-}" ]]; then
    continue
  fi

  base="$(basename "$fq")"
  sample="${base%.fastq.gz}"
  sample="${sample%.fastq}"

  OUT_BAM="${OUT_DIR}/${sample}.hisat2.sorted.bam"

  echo "=== Aligning single-end: ${sample} ==="
  echo "Input FASTQ: $fq"
  echo "Output BAM : $OUT_BAM"

  hisat2 \
    -x "${HISAT2_INDEX_PREFIX}" \
    -U "${fq}" \
    -p "${THREADS}" \
  | samtools sort -@ "${THREADS}" -o "${OUT_BAM}"

  samtools index "${OUT_BAM}"

  ALIGNED["$fq"]=1

  echo "Finished single-end: ${sample}"
  echo
done

echo "All alignments complete. Sorted BAM + index (*.bai) are in: ${OUT_DIR}"
