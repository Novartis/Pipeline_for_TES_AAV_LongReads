
# Extract the sequence surrounding the GOT1/2 genes:
GENOME="/Users/davidkr2/genomes/GRCh38.primary_assembly.genome.fa"
OUTPUT="test_genome.fa"


# Define coordinates for GOT1:
CHR="chr10"
START=99176870
END=99650624
# Extract sequence using samtools:
samtools faidx "$GENOME" "${CHR}:${START}-${END}" > "$OUTPUT"

# Define coordinates for GOT2:
CHR="chr16"
START=58487131
END=58954316
# Extract sequence using samtools:
samtools faidx "$GENOME" "${CHR}:${START}-${END}" >> "$OUTPUT"



# Generate minimap index:
minimap2 -k15 -w10 -d test_genome_k15_w10.mmi test_genome.fa



# Then simulate test data:
mv ../simulation/
snakemake all --configfile config/simulate_test_data.yaml
mv ../test/



