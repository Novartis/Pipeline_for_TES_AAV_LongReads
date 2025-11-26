# HESI AAV IS Simulation

Setup for simulating AAV integration site sequencing data for the HESI CGT-TRACS ISA working group.


## Overview

This simulation is setup to test four variables affecting AAV integration site detection:
1. **Human Genome Backgrounds:** Two individual-specific genomes (HG001, HG002)
2. **AAV Constructs:** Three AAV constructs (eGFP, microdystrophin, FXN)
3. **Integration Site Bias:** Two integration site biases (Uniform across the genome, Liver transcriptionally active sites)
4. **Clonality Levels:** Four clonality patterns (90%, 50%, 10% and no clonality).

This results in 48 unique simulation conditions (2 genomes × 3 AAVs × 2 site biases × 4 clonality levels) with 3 technical replicates each, totaling 144 simulation runs.





## External Software

The following command-line tools must be installed and available in your PATH:

| Tool | Purpose | Used In | Repository | Version (SHA) |
|------|---------|---------|------------|---------------|
| `seqtk` | Sequence toolkit for FASTA/FASTQ manipulation | Genome preparation | [lh3/seqtk](https://github.com/lh3/seqtk) | `94e7070` |
| `samtools` | SAM/BAM file processing and genome indexing | Genome preparation | [samtools/samtools](https://github.com/samtools/samtools) | `906a3b5` |
| `bcftools` | VCF/BCF file processing and variant manipulation | Genome preparation, haplotype generation | [samtools/bcftools](https://github.com/samtools/bcftools) | `02ee548` |
| `htslib` | C library for high-throughput sequencing data | Dependency for samtools/bcftools | [samtools/htslib](https://github.com/samtools/htslib) | `0cadce2` |







## Setup (One-time)

### 1. Generate reference genomes

This step downloads GRCh38 reference and generates individual-specific pseudo-haplotypes from Genome in a Bottle (GIAB) samples (HG001, HG002, HG005). SNPs are applied to create realistic human genome variants.

```bash
cd genomes
./prepare_genomes.sh > prepare_genomes.log 2>&1
cd ..
```

**What this does:**
- Downloads GRCh38 reference genome (no-alt analysis set)
- Downloads GIAB high-confidence variant calls
- Filters for passing SNPs only
- Applies genotype 1 variants to reference using bcftools consensus
- Generates genome-specific FASTA files with indices

**Output files:**
- `GRCh38.fa` + `.fai` (reference genome)
- `HG001_GRCh38.hap1.fa` + `.fai` (individual genome)
- `HG002_GRCh38.hap1.fa` + `.fai` (individual genome)
- `HG005_GRCh38.hap1.fa` + `.fai` (individual genome)
- VCF files with variant calls


### 2. Generate integration site bias

This step creates BED files with integration site sampling probabilities based on tissue-specific gene expression from GTEx. Higher expression correlates with higher integration probability.
These files are pre-generated and should not be regenerated, but the script is included for reference.

```bash
### Only for reference, no need to rerun!!!
cd integration_site_bias
./prepare_integration_sites.sh
cd ..
```

**What this does:**
- Downloads GTEx v10 median TPM data
- Extracts liver and muscle skeletal tissue expression
- Downloads GENCODE v49 gene annotations
- Maps gene expression to transcript coordinates
- Generates probability BED files scaled proportionally to expression level

**Output files:**
- `liver_integration_bias.bed` (liver tissue bias)
- `muscle_integration_bias.bed` (muscle tissue bias)
- `GTEx_TPM_liver_muscle.tsv.gz` (expression data)
- `transcript_ranges.bed.gz` (gene coordinates)


## AAV Constructs

The pipeline includes three AAV constructs:

### CMV-EGFP-WPRE
CMV driven eGFP reporter construct  
**Source:** Novartis manuscript (Coratella et al. 2025, unpublished)  

### rAAVrh74.MHCK7.microdystrophin
MHCK7 driven micro-dystrophin  
**Source:** US patent US20220364117A1 (SEQ ID NO: 3)
**Link:** https://patents.google.com/patent/US20220364117A1  

### pFZ-CBI-CAG-FXNco3
CAG driven codon-optimized frataxin  
**Source:** Pfizer paper (Sheehan et al. 2024)  

Each AAV folder (`AAVs/{construct_name}/`) contains:
- `{construct}.fa` - Full plasmid sequence
- `{construct}_breakend_prob.bed` - AAV breakpoint probabilities
- `{construct}_ISlen.txt` - Empirical integration fragment lengths
- `{construct}_tes_panel.fa` - Target enrichment panel sequences


## Parameterization

Parameters for the simulation are derived from empirical AAV integration data (Sheehan et al. 2024).

### Estimate parameters from real data

The Jupyter notebook analyzes real AAV integration datasets to estimate:
- Translocation frequency
- Genomic deletion size distributions  
- Junction insertion size distributions
- AAV fragment size distributions
- AAV breakpoint preferences

Furthermore, it generates a tiles target-enrichment panel from the three plasmid sequences. 
These files are already generated but can be recreated by running the notebook:

```bash
jupyter notebook parameterization/estimate_params.ipynb
```


## Running Simulations

### Quick start

Run all 48 simulations with 3 replicates each (144 total runs):

```bash
./run_all_simulations.sh
```

**Output:** `simulation_output/HESI_ISA_simulation_*`  
**Logs:** `logs/HESI_ISA_simulation_*.log`  
**Run Table** `all_simulations_runtable.tsv`

### Configuration

Edit variables in `run_all_simulations.sh`:

```bash
# Number of simulations to run simultaneously
N_CONCURRENT_RUNS=2

# Number of CPU cores per simulation
N_SMK_CORES=8

# Dry-run mode (test without executing)
DRY_RUN=false  # or true, 1, 2
```

### Dry-run mode

Test the pipeline without executing:

```bash
# Print commands without running snakemake
DRY_RUN=true ./run_all_simulations.sh

# Run snakemake in dry-run mode
DRY_RUN=2 ./run_all_simulations.sh
```

### What the script does

1. Runs `setup_simulations.py` to generate 48 config files
2. Creates `all_simulations_runtable.tsv` with all run parameters
3. Launches simulations in parallel (N_CONCURRENT_RUNS at a time)
4. Logs each simulation to `logs/HESI_ISA_simulation_*.log`


### Job Scheduler

The above `./run_all_simulations.sh` can be pushed to the UGE cluster using
```bash
qsub qsub_run_all_simulations.sh
```


Alternatively, each simulation can be submitted as a separate job using the provided `run_all_simulations_qsub.sh` script.
This allows running multiple simulations in parallel.
Run the script:
```bash
./run_all_simulations_qsub.sh
```


## Simulation Design

### Parameter Space

**Full factorial design:** 2 × 3 × 2 × 4 = 48 conditions

| Variable | Levels | Values |
|----------|--------|--------|
| Genome | 2 | HG001, HG002 |
| AAV Construct | 3 | CMV-EGFP-WPRE, microdystrophin, FXNco3 |
| Site Bias | 2 | Uniform, Liver |
| Clonality | 4 | 900/100, 500/500, 100/900, 1000 |

### Clonality Levels

- **900/100:** 101 unique events, one event expanded 900× (high clonality, 90%)
- **500/500:** 501 unique events, one event expanded 500× (medium clonality, 50%)
- **100/900:** 901 unique events, one event expanded 100× (low clonality, 10%)
- **1000:** 1000 unique events, no expansion (polyclonal)

### Replicates

Each condition is run 3 times with different random seeds:
- Simulation 1: seeds 0, 1, 2
- Simulation 2: seeds 100, 101, 102
- ...
- Simulation 48: seeds 4700, 4701, 4702



## Post-Simulation Processing

### Calculate MD5 checksums

It is recommended to calculate the MD5 checksums for all simulated FASTQ files and add them to the run table.
This can be repeated, inserting a new column name with md5sum using the `--md5_colname` option with the following command:

```bash
python add_md5sum.py \
  --run_table_in all_simulations_runtable.tsv \
  --md5_colname fastq_md5sum \
  --verbose
```

**What this does:**
- Reads `all_simulations_runtable.tsv`
- For each FASTQ file, calculates MD5 checksum (decompressing .gz files)
- Adds checksum as new column to the table
- Leaves cells empty if FASTQ file doesn't exist

The updated run table can be used to verify data integrity after transfer.
The `all_simulations_runtable_NVS.tsv` contains the simulation run table and MD5 checksums for the full simulation run on Novartis' servers.

