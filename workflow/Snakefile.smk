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
import sys, os, glob, copy
from pathlib import Path
from os.path import join
import pandas as pd




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

# Map reads to the AAV sequence and the genome and extract all mapped reads:
rule extract_mapped_reads:
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
        genome_index=config["ref_mmi"],
        m=config["minimap_settings"]["m"],
        s=config["minimap_settings"]["s"],
        minimap_flags=config["minimap_settings"]["flags_filter_aln"],
        filter_flags=config["minimap_settings"]["flags_filter_aln_filt"],
        scripts_dir=str(config["scripts_dir"]),
    log:
        str(config["out_folder"] / "02_filtered/{sample}.extract_mapped_reads.log"),
    shell:
        "(minimap2 -t {threads} -a {params.minimap_flags} -w2 -m{params.m} -s{params.s} {params.plasmid_seq} {input} | "
        "python {params.scripts_dir}/filter_sam.py {params.filter_flags} | "
        "samtools fastq -F 4 | "
        "minimap2 -t {threads} -a {params.minimap_flags} -m{params.m} -s{params.s} {params.genome_index} - | "
        "python {params.scripts_dir}/filter_sam.py {params.filter_flags} | "
        "samtools fastq -F 4 | "
        "gzip) > {output} 2> >(tee {log} >&2)"
        # Align reads using minimap2.
        # The minimizer window size is adjusted down to -w2 to allow -m/-s to control shortest segment size
        # for AAV alignments.
        # The -m and -s parameters were decreased to 50 (or specified in config)
        # to allow small (50 bp) AAV integrations and therefore short alignments.
        # See minimap2 for reference: https://lh3.github.io/minimap2/minimap2.html
        # This alignment is solely to filter reads that do not map both to
        # the AAV sequence and the genome.
        # samtools "-F 4" excludes unmapped reads.


# Extract the discarded reads:
rule extract_discarded_reads:
    input:
        orig_fastq=str(config["out_folder"] / "01_cleaned/{sample}.cleaned.fastq.gz"),
        filt_fastq=str(config["out_folder"] / "02_filtered/{sample}.filtered.fastq.gz"),
    output:
        filt_fastq_ids=str(config["out_folder"] / "02_filtered/{sample}.filtered.fastq_ids.txt"),
        disc_fastq=str(config["out_folder"] / "02_filtered/{sample}.discarded.fastq.gz"),
    threads: 1
    params:
        outdir=str(config["out_folder"] / "02_filtered")
    log:
        str(config["out_folder"] / "02_filtered/{sample}.extract_discarded_reads.log"),
    shell:
        """
        # Extract read IDs that passed filters:
        (zgrep '^@' {input.filt_fastq} | cut -f1 | sed 's/^.'// > {output.filt_fastq_ids} &&

        # Extract reads that did not pass filters:
        seqkit grep --invert-match --pattern-file {output.filt_fastq_ids} {input.orig_fastq} | 
        gzip) > {output.disc_fastq} 2> >(tee {log} >&2)
        """




####################################
# Align reads to the custom genome #
####################################

# First make custom reference composed of the genome with
# the AAV sequence added as an extra chromosome.
rule make_custom_index:
    output:
        custom_fa=str(config["out_folder"] / "custom_index/{plasmid}_custom_genome_index.fa"),
        custom_mmi=str(config["out_folder"] / "custom_index/{plasmid}_custom_genome_index.mmi"),
    params:
        plasmid_seq=get_plasmid_seq,
        ref_genome=config["ref_genome"],
        minimap_flags=config["minimap_settings"]["flags_custom_genome_aln"],
        tmp_fa=str(config["out_folder"] / "custom_index/{plasmid}_custom_genome_index.fa"),
    threads:
        config["threads"]["minimap2"]
    log:
        str(config["out_folder"] / "custom_index/{plasmid}.make_custom_index.log"),
    shell:
        "(cat {params.plasmid_seq} {params.ref_genome} > {output.custom_fa} && "
        "minimap2 {params.minimap_flags} -d {output.custom_mmi} {output.custom_fa} && "
        "samtools faidx {output.custom_fa}) 2>&1 | tee {log}"
        # Concatenate the plasmid and genome, then make index
        # and finally make fasta index.
        # Notice, that these commands are chained using parenthesis
        # encapsulation to allow all stderr and stdout to be redirected.


# Then, align to the custom reference:
def custom_index_name(wildcards):
    return(str(config["out_folder"] / ("custom_index/" + get_plasmid_name(wildcards) + "_custom_genome_index.mmi")))

rule align_to_genome:
    input:
        genome_index=custom_index_name,
        filt_fastq=str(config["out_folder"] / "02_filtered/{sample}.filtered.fastq.gz"),
    output:
        str(config["out_folder"] / "03_aligned/both/{sample}.bam"),
    params:
        flags=config["minimap_settings"]["flags_custom_genome_aln"],
        scripts_dir=str(config["scripts_dir"]),
    threads:
        config["threads"]["minimap2"]
    log:
        str(config["out_folder"] / "03_aligned/{sample}.align_to_genome.log"),
    shell:
        "(minimap2 -t {threads} -a {params.flags} {input.genome_index} {input.filt_fastq} | "
        "samtools sort | "
        "samtools view -b -F 4) > {output} 2> >(tee {log} >&2)"


# Split the alignments into separate files for AAV and genome alignments:
rule split_alignments:
    input:
        str(config["out_folder"] / "03_aligned/both/{sample}.bam"),
    output:
        aav_aln=str(config["out_folder"] / "03_aligned/AAV/{sample}_AAV_aln.bam"),
        genome_aln=str(config["out_folder"] / "03_aligned/genome/{sample}_genome_aln.bam"),
    params:
        aav_chr_names=get_plasmid_chr_name,
        scripts_dir=str(config["scripts_dir"]),
    threads: 1
    log:
        str(config["out_folder"] / "03_aligned/{sample}.split_alignments.log"),
    shell:
        """
        # Index the BAM file first to allow extraction of alignments:
        (samtools index {input} &&

        # Extract the alignments to the AAV sequence:
        samtools view -h -b {input} {params.aav_chr_names} | samtools sort > {output.aav_aln} &&
        samtools index {output.aav_aln} &&

        # Extract the alignments to the genome (excluding AAV):
        # samtools view -h {input} | awk '$3 != "{params.aav_chr_names}"' | samtools sort > {output.genome_aln} &&
        samtools view -h {input} | python {params.scripts_dir}/filter_sam.py --ref_excl {params.aav_chr_names} | samtools sort > {output.genome_aln} &&
        samtools index {output.genome_aln}) 2> >(tee {log} >&2)
        """




#######################
# Identify insertions #
#######################

# Identify reads that appear to cover an AAV insertion:
rule identify_insertions:
    input:
        genome_aln=str(config["out_folder"] / "03_aligned/genome/{sample}_genome_aln.bam"),
    output:
        str(config["out_folder"] / "04_insertions/{sample}_identified_insertions.txt"),
    params:
        min_MAPQ=config["identify_insertions"]["min_MAPQ"],
        min_len=config["identify_insertions"]["min_len"],
        max_merge_distance=config["identify_insertions"]["max_merge_distance"],
    threads: 1
    log:
        str(config["out_folder"] / "04_insertions/{sample}.identify_insertions.log"),
    shell:
        """
        # Sort BAM:
        (samtools sort {input.genome_aln} | \

        # Convert the BAM file to BED format:
        bedtools bamtobed | \

        # Filter for MAPQ >= 50 and length >= 30:
        awk '{{ if($5 >= {params.min_MAPQ} && ($3 - $2) >= {params.min_len}) {{print $0}} }}' | \

        # Sort the BED file by chromosome and start position:
        sort -k1,1 -k2,2n | \

        # Merge close intervals, calculate counts, mean quality,
        # and distinct strand information (d -1, only mege overlapping features),
        # keep row IDs:
        bedtools merge -d {params.max_merge_distance} -c 4,5,6,4 -o count_distinct,mean,distinct,distinct | \

        # Sort the merged intervals by the count of distinct intervals in descending order:
        sort -k4,4nr) > {output} 2> >(tee {log} >&2)
        """




################################################
# Identify breakpoints and structrual variants #
################################################

rule identify_SVs:
    input:
        insertions=str(config["out_folder"] / "04_insertions/{sample}_identified_insertions.txt"),
        filt_fastq=str(config["out_folder"] / "02_filtered/{sample}.filtered.fastq.gz"),
        genome_index=custom_index_name,
    output:
        str(config["out_folder"] / "05_SVs/{sample}/SV_inference_done"),
    params:
        outdir=str(config["out_folder"] / "05_SVs/{sample}"),
        minimap_flags=config["minimap_settings"]["flags_custom_genome_aln"],
        sniffles_flags=' '.join(config["sniffles_flags"].values()),
        scripts_dir=str(config["scripts_dir"]),
        aav_chr_names=get_plasmid_chr_name,
    threads:
        config["threads"]["sniffles"]
    log:
        str(config["out_folder"] / "05_SVs/{sample}/identify_SVs.log"),
    shell:
        "bash {params.scripts_dir}/SV_inference.sh {input.insertions} {input.filt_fastq} "
        "{input.genome_index} \"{params.aav_chr_names}\" {params.outdir} {params.scripts_dir} {threads} "
        "\"{params.minimap_flags}\" \"{params.sniffles_flags}\" {output} 2>&1 | tee {log}"


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




################
# Run all rule #
################

# Make a list of files to generate if the "all" rule is invoked:
all_rule_input = list()
for sample, *_ in get_ss_samples():
    # all_rule_input.append(str(config["out_folder"] / f"01_cleaned/{sample}.cleaned.fastq.gz"))
    # all_rule_input.append(str(config["out_folder"] / f"02_filtered/{sample}.filtered.fastq.gz"))
    # all_rule_input.append(str(config["out_folder"] / f"02_filtered/{sample}.discarded.fastq.gz"))
    # all_rule_input.append(str(config["out_folder"] / f"03_aligned/both/{sample}.bam"))
    all_rule_input.append(str(config["out_folder"] / f"04_insertions/{sample}_identified_insertions.txt"))
    all_rule_input.append(str(config["out_folder"] / f"05_SVs/{sample}/SV_inference_done"))
    all_rule_input.append(str(config["out_folder"] / "custom_index/.cleanup_marker"))


rule all:
    input:
        all_rule_input


# Create rule DAG:
# snakemake --dag all --configfile config.yaml | dot -Tsvg > dag.svg


