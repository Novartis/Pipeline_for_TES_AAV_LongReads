# Code for analysis of TES data

This code is for AAV integration site analysis using Target Enrichment Sequencing (TES).
The pipeline is for data generated with the [Twist Long-Read Library Preparation and Standard Hyb v2 Enrichment protocol](https://www.twistbioscience.com/resources/protocol/long-read-library-preparation-and-standard-hyb-v2-enrichment) and sequenced with PacBio.
This code was optimized with data simulating random AAV integration and real clonal data from samples each containing data from a single expanded clone.



## Repo content
```
project-root/
├── workflow/
│   ├── Snakemake.smk (main snakemake file)
│   ├── rules/
│   │   └── *.smk (secondary snakemake files)
│   └── scripts/
│       └── Scripts invoked by pipeline
│
├── config/
│   ├── default_config.yaml (default config file)
│   ├── *.yaml (project specific config files)
│   └── *.tsv (project specific sample sheets)
│
├── simulation/
│   └── Snakemake simulation pipeline
│
├── qsub_IO/
│   └── qsub scripts for job submission
│
└── test/
    ├── Input/output for pipeline test
    └── simulated_data/
        └── Small test dataset
```


## Pipeline Workflow

The workflow follows these main steps:

1. **Adapter Trimming** (`clean_reads`)
   - Removes sequencing adapters using cutadapt

2. **AAV Read Filtering** (`extract_AAV_mapped_reads`)
   - Maps reads to AAV sequence using minimap2
   - Retains only reads containing AAV sequences

3. **Genome Alignment** (`align_to_genome`)
   - Aligns filtered reads to the reference genome
   - Filters alignments for quality and chimeric reads

4. **Insertion Identification** (`identify_insertions`)
   - Converts alignments to genomic intervals
   - Identifies potential insertion sites based on read coverage

5. **Structural Variant Detection** (`identify_SVs`)
   - Creates a custom reference combining AAV and genome
   - Uses Sniffles to identify structural variants at AAV-genome junctions


## Key Features

- Handles long reads from PacBio sequencing
- Optimized for detecting AAV integration sites
- Includes parameters for controlling analysis sensitivity
- Comes with a simulation module for generating test data
- Provides a test dataset for validating the pipeline



## Setup
The default pipeline configurations are stored in: `config/default_config.yaml`
Most of these configurations should remain unchanged but a few should change based on the input data:
* `ref_genome`: Path to reference genome fasta file
* `ref_mmi`: Path to minimap2 index of reference genome
* `outer_adapter_5` and `outer_adapter_3`: Outer adapter sequences


### Reference genome
It is not irrelevant which reference genome file is used.
See [Heng Li's blog post](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) on this. 
We use his recommended human reference genome.

The reference genome must be provided as an uncompressed fasta file along with its minimap2 index.
The index can be generated like so:
`minimap2 -x map-hifi -d ref_k19_w19.mmi ref.fa`


### Requirements
The following programms need to be installed and available in the environment:
  + [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
  + [cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)
  + [minimap2](https://github.com/lh3/minimap2)
  + [bedtools](https://github.com/arq5x/bedtools2)
  + [samtools](https://www.htslib.org/download/)
  + [sniffles](https://github.com/fritzsedlazeck/Sniffles)
  + [seqkit](https://github.com/shenwei356/seqkit)

Additionally, the following Python packages needs to be installed:
  + pandas


### Test run
The `test` folder contains a tiny test example to test if the pipeline can run to completion.
If the above mentioned requirements are correctly installed, this test should run successfully in a few seconds without requiring changes to any config files.
Use this test as a template for setting up a run on a set of real samples.



## Running samples
To run a set of samples, it is recommended to generate a yaml config file with the following minimum keys:
```
# Example sample config file  #
samplesheet: "mydir/mysample_samplesheet.tsv"
working_dir: "mydir/TES_pipeline_output"
out_folder: "mysamples"
```
This file can then be extended with other keys e.g. defining the reference genome (`ref_genome`) and these definitions will supersede the definitions in the default pipeline configurations.

To start the pipeline, the following command can be used:
```
snakemake all --printshellcmds --jobs 4 --snakefile ../workflow/Snakefile.smk --configfile ../config/<mysamples_configfile>.yaml
```
Change the `--jobs` parameter to reflect the number of CPU cores to run on.
Also change the `--snakefile` and `--configfile` parameters to reflect the paths to the main snakemake file and the sample set config file, respectively.
The pipeline will then generate a folder, specified by the `out_folder` key in the config file, and populate with processed data.



