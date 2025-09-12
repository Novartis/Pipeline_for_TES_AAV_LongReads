# Script to perform structural variant inference #
# For each detected insertion, the script does the following:
# 1: Extracts the reads that have been marked as containing an insertion.
# 2: Map these reads to the custom reference containing both genome and plasmid sequence.
# 3: Run "sniffles" to detect structural variants.


# Read input:
INSERTIONS=$1
FASTQ_FILE=$2
REFERENCE=$3
OUTDIR=$4
mkdir -p $OUTDIR
THREADS=$5
SNIFFLES_PARAMS=$6
TOUCH=$7


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
  minimap2 -t "$THREADS" -a -x map-hifi -s25 -m25 "$REFERENCE" - | \
  samtools sort -o "$BAM_FILE"
  samtools index "$BAM_FILE"
  samtools depth "$BAM_FILE" > "$COV_FILE"

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
  # sniffles --input "$BAM_FILE" --vcf "$VCF_FILE" --mosaic --mosaic-include-germline --minsvlen 20 --minsupport 1 --bnd-min-split-length 100 --qc-output-all
  eval "sniffles --input $BAM_FILE --vcf $VCF_FILE $SNIFFLES_PARAMS"

  # Clean-up:
  rm "$BAM_FILE"
  rm "$PATTERN_FILE"
  echo "Processed: $INS_PREFIX (SVs saved to $VCF_FILE)"

done < "$INSERTIONS"

# Touch file to indicate success:
touch "$TOUCH"



