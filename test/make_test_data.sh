
# Extract the sequence surrounding the GOT1/2 genes:
GENOME="/Users/davidkr2/genomes/GRCh38.primary_assembly.genome.fa"
OUTPUT="test_genome.fa"


# Define coordinates for GOT1:
CHR="chr10"
START=99176870
END=99650624
# Extract sequence using samtools:
samtools faidx "$GENOME" "${CHR}:${START}-${END}" | sed 's/^>chr10.*$/>chr10.GOT1/' > "$OUTPUT"

# Define coordinates for GOT2:
CHR="chr16"
START=58487131
END=58954316
# Extract sequence using samtools:
samtools faidx "$GENOME" "${CHR}:${START}-${END}" | sed 's/^>chr16.*$/>chr16.GOT2/' >> "$OUTPUT"



# Generate minimap index:
minimap2 -x map-hifi -d test_genome_k19_w19.mmi test_genome.fa
minimap2 -x splice:hq -d test_genome_k15_w5.mmi test_genome.fa


# Then simulate test data:
# cd ../simulation/
snakemake all --snakefile ../simulation/workflow/Snakefile.smk --configfile ../simulation/config/simulate_test_data.yaml
# cd ../test/



