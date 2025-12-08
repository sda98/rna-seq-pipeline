#!/usr/bin/env bash
# 04_count.sh
# Count reads per gene using htseq-count on coordinate-sorted BAM files
#
# Default input : results/alignments  (from 03_align.sh)
# Default output: results/counts
#
# Requires:
#   - htseq-count (HTSeq)
#   - a GTF annotation file (set via GTF_ANNOTATION env var or 3rd argument)
#
# Assumes:
#   - BAM files are coordinate-sorted (samtools sort)
#   - Data is unstranded (change -s no if needed)

set -euo pipefail

#########################
# Configuration
#########################

# Directories (can override with positional args)
BAM_DIR="${1:-results/alignments}"
OUT_DIR="${2:-results/counts}"

# GTF annotation file: set via env var or 3rd argument
GTF_ANNOTATION="${GTF_ANNOTATION:-${3:-/path/to/annotation.gtf}}"

# htseq-count options (edit if needed)
STRANDED="${STRANDED:-no}"   # "yes", "no", or "reverse"
FEATURE_TYPE="${FEATURE_TYPE:-exon}"
ID_ATTRIBUTE="${ID_ATTRIBUTE:-gene_id}"

#########################
# Checks
#########################

if ! command -v htseq-count >/dev/null 2>&1; then
  echo "ERROR: htseq-count not found in PATH." >&2
  echo "Install HTSeq (e.g., pip install HTSeq) before running this script." >&2
  exit 1
fi

if [[ "${GTF_ANNOTATION}" == "/path/to/annotation.gtf" ]]; then
  echo "ERROR: GTF annotation file is not set." >&2
  echo "Set GTF_ANNOTATION env var or pass as 3rd argument, e.g.:" >&2
  echo "  bash scripts/04_count.sh results/alignments results/counts /path/to/annotation.gtf" >&2
  exit 1
fi

if [[ ! -f "${GTF_ANNOTATION}" ]]; then
  echo "ERROR: GTF annotation file not found at: ${GTF_ANNOTATION}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

echo "Running htseq-count..."
echo "BAM directory      : ${BAM_DIR}"
echo "Output directory   : ${OUT_DIR}"
echo "GTF annotation     : ${GTF_ANNOTATION}"
echo "Strandedness (-s)  : ${STRANDED}"
echo "Feature type (-t)  : ${FEATURE_TYPE}"
echo "ID attribute (-i)  : ${ID_ATTRIBUTE}"
echo

#########################
# Main loop
#########################

shopt -s nullglob

for BAM in "${BAM_DIR}"/*.bam; do
  [[ -e "$BAM" ]] || continue  # skip if no BAMs

  # Skip index files just in case
  [[ "$BAM" == *.bai ]] && continue

  fname="$(basename "$BAM")"
  sample="${fname%.bam}"

  OUT_FILE="${OUT_DIR}/${sample}.counts.txt"

  echo "=== Counting for sample: ${sample} ==="
  echo "Input BAM : $BAM"
  echo "Output    : $OUT_FILE"

  htseq-count \
    -f bam \
    -r pos \
    -s "${STRANDED}" \
    -t "${FEATURE_TYPE}" \
    -i "${ID_ATTRIBUTE}" \
    "${BAM}" \
    "${GTF_ANNOTATION}" \
    > "${OUT_FILE}"

  echo "Finished sample: ${sample}"
  echo
done

echo "All counting complete. Per-sample count files are in: ${OUT_DIR}"
