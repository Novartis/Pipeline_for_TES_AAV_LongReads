#!/bin/bash
"""
Copyright 2025 Novartis Institutes for BioMedical Research Inc.
 
Licensed under the MIT License (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
 
https://www.mit.edu/~amini/LICENSE.md
 
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
# run_enrichment.sh
# Runs permutation-based enrichment tests for AAV integration sites against
# liver-expressed gene features. Outputs a TSV file with results.
#
# Usage:
#   ./run_enrichment.sh <input.bed> [output_results.tsv]
#
# Arguments:
#   input.bed            BED file of AAV integration sites (see README for format)
#   output_results.tsv   Output path (default: output/enrichment_results.tsv)

set -e

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input.bed> [output_results.tsv]" >&2
    echo "" >&2
    echo "  input.bed   BED6 file of AAV integration sites" >&2
    echo "  See README.md for input format requirements." >&2
    exit 1
fi

INPUT_FILE="$1"

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input file not found: $INPUT_FILE" >&2
    exit 1
fi

# Resolve script directory
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

AAV_BED="$INPUT_FILE"

GENOME_SIZES="$BASE_DIR/reference/hg38.chrom.sizes"
BLACKLIST="$BASE_DIR/reference/blacklist.v2.bed"
ITERATIONS=1000

RESULTS_TSV="${2:-$BASE_DIR/output/enrichment_results.tsv}"
mkdir -p "$(dirname "$RESULTS_TSV")"

# Validate reference files
for REF_FILE in "$GENOME_SIZES" "$BLACKLIST"; do
    if [[ ! -f "$REF_FILE" ]]; then
        echo "ERROR: Reference file not found: $REF_FILE" >&2
        exit 1
    fi
done

# --- Restrict genome sizes to standard chromosomes (chr1-22, X, Y) ---
GENOME_SIZES_STD=$(mktemp)
grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)\b' "$GENOME_SIZES" > "$GENOME_SIZES_STD"
trap 'rm -f "$GENOME_SIZES_STD"' EXIT

# --- Feature definitions: "Label\tBED file\tCategory" ---
declare -A FEATURE_BED
declare -A FEATURE_CAT
FEATURE_ORDER=()

add_feature() {
    local label="$1" bed="$2" cat="$3"
    FEATURE_ORDER+=("$label")
    FEATURE_BED["$label"]="$bed"
    FEATURE_CAT["$label"]="$cat"
}

# Liver-expressed gene features (user must provide these BED files in features/):
add_feature "Liver Expressed Genes (high)"                "$BASE_DIR/features/liver_expressed_genes_high.bed"              "positive"
add_feature "Liver Expressed Genes (high TSS +/-10kb)"    "$BASE_DIR/features/liver_expressed_genes_high_tss_10kb.bed"     "positive"
add_feature "Liver Expressed Genes (high TSS +/-50kb)"    "$BASE_DIR/features/liver_expressed_genes_high_tss_50kb.bed"     "positive"
add_feature "Liver Expressed Genes (high TSS +/-100kb)"   "$BASE_DIR/features/liver_expressed_genes_high_tss_100kb.bed"    "positive"
add_feature "Liver Expressed Genes (mid+high)"            "$BASE_DIR/features/liver_expressed_genes_mid_high.bed"          "positive"
add_feature "Liver Expressed Genes (mid+high TSS +/-10kb)"  "$BASE_DIR/features/liver_expressed_genes_mid_high_tss_10kb.bed" "positive"
add_feature "Liver Expressed Genes (mid+high TSS +/-50kb)"  "$BASE_DIR/features/liver_expressed_genes_mid_high_tss_50kb.bed" "positive"
add_feature "Liver Expressed Genes (mid+high TSS +/-100kb)" "$BASE_DIR/features/liver_expressed_genes_mid_high_tss_100kb.bed" "positive"

# -------------------------------------------------------
# Write TSV header
# -------------------------------------------------------
AAV_COUNT=$(grep -c '' "$AAV_BED" 2>/dev/null || echo 0)
echo "# aav_sites=$AAV_COUNT iterations=$ITERATIONS generated=$(date '+%Y-%m-%dT%H:%M:%S')" > "$RESULTS_TSV"
printf 'label\tcategory\tstatus\tobserved\tmean_random\tgreater_equal\titerations\n' >> "$RESULTS_TSV"

# -------------------------------------------------------
# Run each feature
# -------------------------------------------------------
for label in "${FEATURE_ORDER[@]}"; do
    bed="${FEATURE_BED[$label]}"
    cat="${FEATURE_CAT[$label]}"

    echo "Running [$label]..."

    if [[ ! -f "$bed" ]]; then
        echo "  WARNING: Feature BED not found: $bed — skipping"
        printf '%s\t%s\tMISSING\tNA\tNA\tNA\tNA\n' "$label" "$cat" >> "$RESULTS_TSV"
        continue
    fi

    # Observed overlap count
    OBSERVED=$(bedtools intersect -u -a "$AAV_BED" -b "$bed" | wc -l)

    # Permutation test
    RANDOM_COUNTS_FILE=$(mktemp)
    for i in $(seq 1 $ITERATIONS); do
        RANDOM_OVERLAP=$(bedtools shuffle -i "$AAV_BED" -g "$GENOME_SIZES_STD" -excl "$BLACKLIST" -noOverlapping \
            | bedtools intersect -u -a stdin -b "$bed" \
            | wc -l)
        echo "$RANDOM_OVERLAP" >> "$RANDOM_COUNTS_FILE"
    done

    # Statistics
    MEAN_RANDOM=$(awk '{ sum += $1 } END { printf "%.2f", sum/NR }' "$RANDOM_COUNTS_FILE")
    GREATER_EQUAL=$(awk -v obs="$OBSERVED" '$1 >= obs { count++ } END { print count+0 }' "$RANDOM_COUNTS_FILE")

    printf '%s\t%s\tok\t%d\t%s\t%d\t%d\n' "$label" "$cat" "$OBSERVED" "$MEAN_RANDOM" "$GREATER_EQUAL" "$ITERATIONS" >> "$RESULTS_TSV"

    echo "  Observed=$OBSERVED, Mean_random=$MEAN_RANDOM, >= observed: $GREATER_EQUAL/$ITERATIONS"
    rm -f "$RANDOM_COUNTS_FILE"
done

echo ""
echo "Results written to: $RESULTS_TSV"
