# Copyright 2025 Novartis Institutes for BioMedical Research Inc.
 
# Licensed under the MIT License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
 
# https://www.mit.edu/~amini/LICENSE.md
 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Script to perform structural variant inference #
# For each detected insertion, the script does the following:
# 1: Extracts the reads that have been marked as containing an insertion.
# 2: Map these reads to the custom reference containing both genome and plasmid sequence.
# 3: Run "sniffles" to detect structural variants.


# Read input:
INSERTIONS=$1
FASTQ_FILE=$2
REFERENCE=$3
AAV_NAMES=$4
OUTDIR=$5
SCRIPT_DIR=$6
mkdir -p $OUTDIR
THREADS=$7
MINIMAP_PARAMS=$8
SNIFFLES_PARAMS=$9
TOUCH=${10}

# echo $INSERTIONS
# echo $FASTQ_FILE
# echo $REFERENCE
# echo $AAV_NAMES
# echo $OUTDIR
# echo $SCRIPT_DIR
# echo $THREADS
# echo $MINIMAP_PARAMS
# echo $SNIFFLES_PARAMS
# echo $TOUCH

# For each insertion:
while IFS=$'\t' read -r chrom start end count mapq strand read_ids breakpoint; do
  # Make insertion prefix for filenames:
  INS_PREFIX="${OUTDIR}/region_${chrom}_${start}_${end}"
  BAM_FILE="${INS_PREFIX}.bam"                # temporary file
  VCF_FILE="${INS_PREFIX}.vcf"                # output file
  COV_FILE="${INS_PREFIX}.coverage.txt"       # output file
  PATTERN_FILE="${INS_PREFIX}_patterns.txt"   # temporary file

  # Extract the read IDs that have been marked as containing an insertion:
  echo "$read_ids" | tr ',' '\n' > "$PATTERN_FILE"
  # Then find and extract the reads corresponding to the read IDs:
  seqkit grep --pattern-file "$PATTERN_FILE" "$FASTQ_FILE" | \

  # Align extracted reads:
  eval "minimap2 -t $THREADS -a $MINIMAP_PARAMS $REFERENCE -" | \
  # Keep at most four top scoring (AS:i tag) alignment segments,
  # two overlapping the AAV, the other two, the highest scoring non-AAV segment:
  # samtools sort -n -O sam | \
  # python "$SCRIPT_DIR/filter_sam.py" --keep_Nseg_chr_chr 2 AS dsc "${AAV_NAMES%% *}" | \
  samtools sort -o "$BAM_FILE"
  samtools index "$BAM_FILE"
  samtools depth -g SECONDARY "$BAM_FILE" > "$COV_FILE"

  # Check to make sure that all read IDs were found, no more no less:
  # First get line counts:
  read_id_lines=$(wc -l < "$PATTERN_FILE")
  bam_lines=$(samtools view "$BAM_FILE" | cut -f1 | sort | uniq | wc -l)

  # Then, compare line counts:
  if [ "$read_id_lines" -ne "$bam_lines" ]; then
      echo "Error: attempted to extract read IDs but got inconsistent results."
      echo "$PATTERN_FILE: $read_id_lines read IDs"
      echo "$BAM_FILE: $bam_lines read IDs"
      exit 1
  fi

  # Get breakpoints and structrual variants:
  eval "sniffles --input $BAM_FILE --vcf $VCF_FILE $SNIFFLES_PARAMS"

  # Clean-up:
  rm "$BAM_FILE"
  rm "$BAM_FILE.bai"
  rm "$PATTERN_FILE"
  echo "Processed: $INS_PREFIX (SVs saved to $VCF_FILE)"

done < "$INSERTIONS"

# Touch file to indicate success:
touch "$TOUCH"



