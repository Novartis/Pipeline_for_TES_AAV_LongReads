"""
Author: M. Lemmens
Affiliation: NIBR
Aim: TES workflow
Date: August, 2024
Run: snakemake   -s Snakefile
"""
import sys, os, glob, copy
from pathlib import Path
from os.path import join
import pandas as pd


# Todo:
# Add some default plotting (see R script for inspiration)




#############
# Configure #
#############

### The pipeline repo has the follow folder structure:
# pipeline repo
#   - workflow
#     - Snakefile
#     - rules/*.smk
#     - scripts/*
#   - config/*.yaml
#   - test/*

# Extract paths based on the above structure:
snakefile_path = Path(workflow.snakefile).resolve()
pipeline_dir = snakefile_path.parent.parent

configfile: pipeline_dir / "config/default_config.yaml"
workdir: config['working_dir']
include: pipeline_dir / "workflow/rules/common.smk"
config['pipeline_dir'] = pipeline_dir
config['scripts_dir'] = pipeline_dir / "workflow/scripts"

# Parse the config file to validate, add information and return a sample sheet:
global ss_df
ss_df = parse_config()

shell.executable("/bin/bash")
shell.prefix("set -o pipefail; ")




##############################
# Prepare reads for analysis #
##############################

# Remove adapter sequences:
rule clean_reads:
    input:
        get_ss_fastq
    output:
        str(config["out_folder"] / "01_cleaned/{sample}.cleaned.fastq.gz"),
    threads: 1
    resources:
        mem_gb=120,
    params:
        # Adapter configuration:
        # <outer_adapter_5> <sample_adapter_5> <sequence> <sample_adapter_3> <outer_adapter_3>
        outer_adapter_5=config["outer_adapter_5"],
        sample_adapters=get_ss_adapters,
        outer_adapter_3=config["outer_adapter_3"],
        cutadapt_flags=' '.join(config["cutadapt_flags"].values()),
    log:
        str(config["out_folder"] / "01_cleaned/{sample}.clean_reads.log"),
    shell:
        "cutadapt {params.cutadapt_flags} "
        "-a ^{params.sample_adapters[0]}{params.outer_adapter_5}...{params.outer_adapter_3}{params.sample_adapters[1]} "
        "-o {output} {input} 2>&1 | tee {log}"




################
# Filter reads #
################

# Map reads to the AAV genome and extract mapped reads:
rule extract_AAV_mapped_reads:
    input:
        str(config["out_folder"] / "01_cleaned/{sample}.cleaned.fastq.gz"),
    output:
        str(config["out_folder"] / "02_filtered/{sample}.filtered.fastq.gz"),
    threads:
        config["threads"]["minimap2"]
    resources:
        mem_gb=120,   
    params:
        plasmid_seq=get_ss_plasmid,
    log:
        str(config["out_folder"] / "02_filtered/{sample}.extract_AAV_mapped_reads.log"),
    shell:
        "(minimap2 -t {threads} -a -x map-hifi -s25 -m25  {params.plasmid_seq} {input} | "
        "samtools fastq -F 4 | gzip) > {output} 2> >(tee {log} >&2)"
        # Align using map-hifi preset for pacbio hifi reads.
        # Additionally, the -m and -s parameters were decreased to 25
        # to allow small AAV integrations and therefore short alignments.
        # -s score of 25 means that the secondary alignment must have a score
        # that is at least 25% of the primary alignment score.
        # -m score describes the chaining score: Chaining score equals
        # the approximate number of matching bases minus a concave gap penalty.
        # See minimap2 for reference: https://lh3.github.io/minimap2/minimap2.html

        # samtools "-F 4" excludes unmapped reads.




#############################
# Align reads to the genome #
#############################

rule align_to_genome:
    input:
        str(config["out_folder"] / "02_filtered/{sample}.filtered.fastq.gz"),
    output:
        str(config["out_folder"] / "03_aligned/{sample}.bam"),
    params:
        genome_index=config["ref_mmi"],
        scripts_dir=str(config["scripts_dir"]),
    threads:
        config["threads"]["minimap2"]
    log:
        str(config["out_folder"] / "03_aligned/{sample}.align_to_genome.log"),
    shell:
        "(minimap2 -t {threads} -a -x map-hifi --cs {params.genome_index} {input} | "
        "python {params.scripts_dir}/filter_sam.py | "
        "samtools view -b -F 4) > {output} 2> >(tee {log} >&2)"
        # Align, to genome, stream sam output to filtering script,
        # then convert to mapped reads BAM format and write to file.




#######################
# Identify insertions #
#######################

# Identify reads that appear to cover an AAV insertion:
rule identify_insertions:
    input:
        str(config["out_folder"] / "03_aligned/{sample}.bam"),
    output:
        str(config["out_folder"] / "04_insertions/{sample}_identified_insertions.txt"),
    params:
        min_MAPQ=config["identify_insertions"]["min_MAPQ"],
        min_len=config["identify_insertions"]["min_len"],
    threads: 1
    log:
        str(config["out_folder"] / "04_insertions/{sample}.identify_insertions.log"),
    shell:
        """
        # Sort BAM:
        (samtools sort {input} | \

        # Convert the BAM file to BED format:
        bedtools bamtobed | \

        # Filter for MAPQ >= 50 and length >= 30:
        awk '{{ if($5 >= {params.min_MAPQ} && ($3 - $2) >= {params.min_len}) {{print $0}} }}' | \

        # Sort the BED file by chromosome and start position:
        sort -k1,1 -k2,2n | \

        # Merge close intervals, calculate counts, mean quality,
        # and distinct strand information (d -1, only mege overlapping features),
        # keep row IDs:
        bedtools merge -d 1 -c 4,5,6,4 -o count_distinct,mean,distinct,distinct | \

        # Sort the merged intervals by the count of distinct intervals in descending order:
        sort -k4,4nr) > {output} 2> >(tee {log} >&2)
        """




################################################
# Identify breakpoints and structrual variants #
################################################

# First make custom reference composed of the genome with
# the AAV sequence added as an extra chromosome.
rule make_custom_index:
    output:
        str(config["out_folder"] / "custom_index/{plasmid}_custom_genome_index.mmi"),
    params:
        plasmid_seq=get_plasmid_seq,
        ref_genome=config["ref_genome"],
        tmp_fa=str(config["out_folder"] / "custom_index/{plasmid}_custom_genome_index.fa"),
    threads:
        config["threads"]["minimap2"]
    log:
        str(config["out_folder"] / "custom_index/{plasmid}.make_custom_index.log"),
    shell:
        "(cat {params.plasmid_seq} {params.ref_genome} > {params.tmp_fa} && "
        "minimap2 -x map-hifi -d {output} {params.tmp_fa} && "
        "rm {params.tmp_fa}) 2>&1 | tee {log}"
        # Concatenate the plasmid and genome, then make index
        # and finally remove the concatenated fasta file.
        # Notice, that these commands are chained using parenthesis
        # encapsulation to allow all stderr and stdout to be redirected.


# Then, identify structural variants using the custom reference:
def custom_index_name(wildcards):
    return(str(config["out_folder"] / ("custom_index/" + get_plasmid_name(wildcards) + "_custom_genome_index.mmi")))

rule identify_SVs:
    input:
        insertions=str(config["out_folder"] / "04_insertions/{sample}_identified_insertions.txt"),
        clean_fastq=str(config["out_folder"] / "01_cleaned/{sample}.cleaned.fastq.gz"),
        genome_index=custom_index_name,
    output:
        str(config["out_folder"] / "05_SVs/{sample}/SV_inference_done"),
    params:
        outdir=str(config["out_folder"] / "05_SVs/{sample}"),
        sniffles_flags=' '.join(config["sniffles_flags"].values()),
        scripts_dir=str(config["scripts_dir"]),
    threads:
        config["threads"]["sniffles"]
    log:
        str(config["out_folder"] / "05_SVs/{sample}/identify_SVs.log"),
    shell:
        "bash {params.scripts_dir}/SV_inference.sh {input.insertions} {input.clean_fastq} "
        "{input.genome_index} {params.outdir} {threads} \"{params.sniffles_flags}\" {output} 2>&1 | tee {log}"


# Remove custom index after use to save space:
all_SV_done = list()
for sample, *_ in get_ss_samples():
    all_SV_done.append(str(config["out_folder"] / f"05_SVs/{sample}/SV_inference_done"))
rule cleanup_temp_folder:
    input:
        all_SV_done
    output:
        temp(str(config["out_folder"] / "custom_index/.cleanup_marker"))
    params:
        folder_to_rm=str(config["out_folder"] / "custom_index")
    shell:
        """
        rm -rf {params.folder_to_rm}
        mkdir -p {params.folder_to_rm}
        touch {output}
        """









###############################
# Extract statistics and plot #
###############################

# Add stuff here










################
# Run all rule #
################

# Make a list of files to generate if the "all" rule is invoked:
all_rule_input = list()
for sample, *_ in get_ss_samples():
    all_rule_input.append(str(config["out_folder"] / f"01_cleaned/{sample}.cleaned.fastq.gz"))
    all_rule_input.append(str(config["out_folder"] / f"02_filtered/{sample}.filtered.fastq.gz"))
    all_rule_input.append(str(config["out_folder"] / f"03_aligned/{sample}.bam"))
    all_rule_input.append(str(config["out_folder"] / f"04_insertions/{sample}_identified_insertions.txt"))
    all_rule_input.append(str(config["out_folder"] / f"05_SVs/{sample}/SV_inference_done"))
    all_rule_input.append(str(config["out_folder"] / "custom_index/.cleanup_marker"))

rule all:
    input:
        all_rule_input


# Create rule DAG:
# snakemake --dag all --configfile config.yaml | dot -Tsvg > dag.svg


