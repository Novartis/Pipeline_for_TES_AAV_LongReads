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


# vim: set syntax=python:

import sys, os, glob, copy
from pathlib import Path
import pandas as pd
from snakemake.utils import min_version
min_version("8.0")

# Set shell prefix for reproducibility
shell.prefix("set -euo pipefail; export PYTHONHASHSEED=0; ")


### The pipeline repo has the following folder structure:
# pipeline repo
#   - workflow
#     - Snakefile
#     - rules/*.smk
#     - scripts/*
#   - config/*.yaml
#   - test_run/*

# Extract paths based on the above structure:
snakefile_path = Path(workflow.snakefile).resolve()
pipeline_dir = snakefile_path.parent.parent

# Allow user to specify alternative default config file via --config default_configfile=path/to/config.yaml
# Otherwise use the standard default_config.yaml
default_config_path = config.get('default_configfile', pipeline_dir / "config/default_config.yaml")
if not Path(default_config_path).is_file():
    sys.exit(f"Error: Default config file {default_config_path} does not exist.")
configfile: default_config_path

# Change directories:
original_workdir = os.getcwd()
config['workdir'] = str(Path(original_workdir) / config['workdir'])
workdir: config['workdir']
include: pipeline_dir / "workflow/rules/common.smk"
config['pipeline_dir'] = pipeline_dir
config['scripts_dir'] = pipeline_dir / "workflow/scripts"
parse_config()



###############
# Input check #
###############

# Check to see that the reference genome conforms with
# the requirements of the downstream programs.
# This should produce an empty file, if not the genome
# has chromosomes too short.
# Also outputs chromosome lengths for breakpoint sampling.
empty_file_sha256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
rule check_reference:
    output:
        check=ensure(config['run_folder'] + "/reference_check", sha256=empty_file_sha256),
        chrom_lengths=config['run_folder'] + "/chromosome_lengths.tsv"
    params:
        max_len=11 * config["read_sim"]["max_len"] + 2 * config["integration"]["deletion"]["max_size"],
        genome_fnam=config["genome_fnam"],
        scrip_dir=config['pipeline_dir'] / "workflow/scripts"
    threads: 1
    shell:
        "python {params.scrip_dir}/ref_length.py --max_len {params.max_len} "
        "--fa {params.genome_fnam} --chrom_lengths_out {output.chrom_lengths} > {output.check}"




###################################
# Simulate AAV integration events #
###################################

# Generate the loci for integration:
rule draw_integration_sites:
    input:
        genome_check=config['run_folder'] + "/reference_check",
        chrom_lengths=config['run_folder'] + "/chromosome_lengths.tsv"
    output:
        integration_sites=config['run_folder'] + "/{repeat}_integration_sites.fa.gz"
    params:
        events=2 * config["integration"]["events"],
        max_len=10 * config["read_sim"]["max_len"] + 2 * config["integration"]["deletion"]["max_size"],
        genome_fnam=config["genome_fnam"],
        scrip_dir=config['pipeline_dir'] / "workflow/scripts",
        prob_bed=lambda wildcards: f"--prob_bed {config['integration']['sampling_prob_bed']}" if config.get('integration', {}).get('sampling_prob_bed', '') else "",
        seed_arg=get_repeat_seed,
        vsearch_seed=lambda wildcards: get_repeat_seed(wildcards, "vsearch")
    threads: 1
    # 1. Use sample_breakpoints.py to sample random genomic breakpoints,
    #    2x to allow for consuming two fragments for translocations
    # 2. Expand each breakpoint to a fragment of length 10x max_len + max deletion size on each side
    # 3. Use bedtools getfasta to extract sequences from the genome
    # 4. Shuffle the fragments with vsearch and return one fragment per AAV integration event
    # 5. gzip output fasta containing the integration site fragments
    # Note: If sampling_prob_bed is specified in config, use weighted sampling
    shell:
        "python {params.scrip_dir}/sample_breakpoints.py "
        "--chrom_lengths {input.chrom_lengths} --n_events {params.events} "
        "--fragment_len {params.max_len} {params.prob_bed} "
        "{params.seed_arg} | "
        "bedtools getfasta -fi {params.genome_fnam} -bed - | "
        "vsearch --fasta_width=0 --maxseqlengt={params.max_len} "
        "{params.vsearch_seed} --output /dev/stdout --shuffle - | "
        "gzip > {output.integration_sites}"


# Add AAV to integration sites and simulate the reads:
rule sim_AAV_loci:
    input:
        integration_sites=config['run_folder'] + "/{repeat}_integration_sites.fa.gz"
    output:
        events_tsv=config['run_folder'] + "/{repeat}_integration_events.tsv",
        sim_reads=config['run_folder'] + "/{repeat}_integration_loci.fa.gz",
    params:
        # For "add_integration.py":
        events=config["integration"]["events"],
        mean_frag_len=config["read_sim"]["mean_frag_len"],
        stdev_frag_len=config["read_sim"]["stdev_frag_len"],
        aav_method=config["integration"]["aav"]["method"],
        aav_prob_file = f"--aav_prob_file {config['integration']['aav']['prob_file']}" if 'prob_file' in config['integration']['aav'] and config['integration']['aav']['prob_file'] else "",
        aav_empirical_sizes = f"--aav_empirical_sizes {config['integration']['aav']['empirical_sizes']}" if 'empirical_sizes' in config['integration']['aav'] and config['integration']['aav']['empirical_sizes'] else "",
        aav_kde_plot = lambda wildcards: f"--aav_kde_plot {config['run_folder']}/aav_kde_plot.png" if 'empirical_sizes' in config['integration']['aav'] and config['integration']['aav']['empirical_sizes'] else "",
        expansion=lambda wildcards: ' '.join(map(str, config["integration"]["expansion"])) if isinstance(config["integration"]["expansion"], list) else str(config["integration"]["expansion"]),
        AAV_genome_fnam=config["AAV_genome_fnam"],
        del_max_len=config["integration"]["deletion"]["max_size"],
        del_flag=config["integration"]["deletion"]["size_dist"],
        ins_flag=config["integration"]["insertion"],
        p_trans=config["integration"]["p_translocation"],

        # For minimap/target enrichment:
        k=config["minimap"]["k"],
        enrich_min_len=config["enrichment"]["min_len"],
        tes_panel_fnam=config["tes_panel_fnam"],
        enrich_model=config["enrichment"]["enrich_model"],
        amplification=config["enrichment"]["amplification"],
        leakage=config["enrichment"]["leakage"],

        scrip_dir=config['pipeline_dir'] / "workflow/scripts",
        seed_arg=get_repeat_seed,
        minimap2_seed=lambda wildcards: get_repeat_seed(wildcards, "minimap2")
    threads: 1
    # 1. Pipe integration site fragments into add_integration.py
    # 2. Sample AAV fragments and integrate them into the 
    #    integration site fragments, then fragment again to
    #    the right size and print these in fastq format to stdout.
    #    Meanwhile, write integration events information to tsv file.
    #    This is done by add_integration.py
    # 3. Align these fastq fragment to the target enrichment panel with minimap2
    # 4. Simulate target enrichment and amplification of the fragments
    #    using target_enrichment.py
    # 5. gzip output fasta containing the target enriched loci
    shell:
        "gunzip -dc {input.integration_sites} | "
        "python {params.scrip_dir}/add_integration.py --fq_out --events {params.events} "
        "--mfl {params.mean_frag_len} --sfl {params.stdev_frag_len} "
        "--aav_method \"{params.aav_method}\" {params.aav_prob_file} {params.aav_empirical_sizes} {params.aav_kde_plot} --exp {params.expansion} "
        "--aav {params.AAV_genome_fnam} --int_out {output.events_tsv} "
        "--del_max {params.del_max_len} --del_flag \"{params.del_flag}\" "
        "--ins_flag \"{params.ins_flag}\" --p_trans {params.p_trans} "
        "{params.seed_arg} | "
        "minimap2 {params.minimap2_seed} --eqx -t1 -k{params.k} -w2 --sr --frag=yes "
        "-A1 -B120 -O20,40 -E40,25 -b0 -r100 -f1000,5000 "
        "-n2 -m{params.enrich_min_len} -s{params.enrich_min_len} "
        "--secondary=no -a {params.tes_panel_fnam} - | "
        "python {params.scrip_dir}/target_enrichment.py --ep_model {params.enrich_model} "
        "--amp {params.amplification} --leak {params.leakage} "
        "{params.seed_arg} | "
        "gzip > {output.sim_reads}"




######################################
# Simulate fragments from the genome #
######################################

# Draw loci from the genome and perform target enrichment:
rule sim_genome_loci:
    input:
        genome_check=config['run_folder'] + "/reference_check"
    output:
        sim_reads=config['run_folder'] + "/{repeat}_genome_loci.fa.gz"
    params:
        # For wgsim fragment simulation:
        N=config["read_sim"]["genome_reads"],
        max_len=config["read_sim"]["max_len"],
        genome_fnam=config["genome_fnam"],

        # For fragment_resizing.py fragment resizing:
        mean_frag_len=config["read_sim"]["mean_frag_len"],
        stdev_frag_len=config["read_sim"]["stdev_frag_len"],
        
        # For minimap/target enrichment:
        k=config["minimap"]["k"],
        enrich_min_len=config["enrichment"]["min_len"],
        tes_panel_fnam=config["tes_panel_fnam"],
        enrich_model=config["enrichment"]["enrich_model"],
        leakage=config["enrichment"]["leakage"],

        scrip_dir=config['pipeline_dir'] / "workflow/scripts",
        seed_arg=get_repeat_seed,
        wgsim_seed=lambda wildcards: get_repeat_seed(wildcards, "wgsim"),
        minimap2_seed=lambda wildcards: get_repeat_seed(wildcards, "minimap2")
    threads: 3
    # 1. Use wgsim to sample same size fragments of the genome
    # 2. Use fragment_resizing.py to randomly resize these fragments
    #    such that the resulting fragment size follow a gamma distriubtion
    # 3. Align these fastq fragment to the target enrichment panel with minimap2
    # 4. Simulate target enrichment and amplification of the fragments
    #    using target_enrichment.py
    # 5. gzip output fasta containing the target enriched loci
    shell:
        "wgsim -e0 -d0 -s0 -r0 -R0 -X0 -N {params.N} -1 {params.max_len} "
        "{params.wgsim_seed} {params.genome_fnam} /dev/stdout /dev/null | "
        "python {params.scrip_dir}/fragment_resizing.py "
        "--mfl {params.mean_frag_len} --sfl {params.stdev_frag_len} "
        "{params.seed_arg} | "
        "minimap2 {params.minimap2_seed} -t4 -k{params.k} -w2 --sr --frag=yes "
        "-A1 -B120 -O20,40 -E40,25 -b0 -r100 -f1000,5000 "
        "-n2 -m{params.enrich_min_len} -s{params.enrich_min_len} "
        "--secondary=no -a {params.tes_panel_fnam} - | "
        "python {params.scrip_dir}/target_enrichment.py --ep_model {params.enrich_model} "
        "--amp 1 --leak {params.leakage} "
        "{params.seed_arg} | "
        "gzip > {output.sim_reads}"




#########################################
# Simulate episomal AAV DNA fragments   #
#########################################

# Generate episomal AAV fragments and perform target enrichment:
rule sim_episomes:
    output:
        sim_reads=config['run_folder'] + "/{repeat}_episome_loci.fa.gz"
    params:
        # For episomal AAV generation:
        N=config["read_sim"]["episomal_reads"],
        AAV_genome_fnam=config["AAV_genome_fnam"],
        aav_method=config["integration"]["aav"]["method"],
        aav_prob_file = f"--aav_prob_file {config['integration']['aav']['prob_file']}" if 'prob_file' in config['integration']['aav'] and config['integration']['aav']['prob_file'] else "",
        aav_empirical_sizes = f"--aav_empirical_sizes {config['integration']['aav']['empirical_sizes']}" if 'empirical_sizes' in config['integration']['aav'] and config['integration']['aav']['empirical_sizes'] else "",
        
        # For minimap/target enrichment:
        k=config["minimap"]["k"],
        enrich_min_len=config["enrichment"]["min_len"],
        tes_panel_fnam=config["tes_panel_fnam"],
        enrich_model=config["enrichment"]["enrich_model"],
        leakage=config["enrichment"]["leakage"],

        scrip_dir=config['pipeline_dir'] / "workflow/scripts",
        seed_arg=get_repeat_seed,
        minimap2_seed=lambda wildcards: get_repeat_seed(wildcards, "minimap2")
    threads: 1
    # 1. Use generate_episomes.py to sample AAV episomal fragments
    # 2. Align these fastq fragments to the target enrichment panel with minimap2
    # 3. Simulate target enrichment and amplification of the fragments
    #    using target_enrichment.py
    # 4. gzip output fasta containing the target enriched episomal loci
    shell:
        "python {params.scrip_dir}/generate_episomes.py "
        "--aav {params.AAV_genome_fnam} --aav_method \"{params.aav_method}\" "
        "{params.aav_prob_file} {params.aav_empirical_sizes} --events {params.N} "
        "{params.seed_arg} | "
        "minimap2 {params.minimap2_seed} --eqx -t1 -k{params.k} -w2 --sr --frag=yes "
        "-A1 -B120 -O20,40 -E40,25 -b0 -r100 -f1000,5000 "
        "-n2 -m{params.enrich_min_len} -s{params.enrich_min_len} "
        "--secondary=no -a {params.tes_panel_fnam} - | "
        "python {params.scrip_dir}/target_enrichment.py --ep_model {params.enrich_model} "
        "--amp 1 --leak {params.leakage} "
        "{params.seed_arg} | "
        "gzip > {output.sim_reads}"




########################################
# Combine fragments and simulate reads #
########################################

# Combine and shuffle loci,
# then simulate reads using badreads,
# then clean up:
rule sim_reads:
    input:
        integration_sites=config['run_folder'] + "/{repeat}_integration_sites.fa.gz",
        integration_loci=config['run_folder'] + "/{repeat}_integration_loci.fa.gz",
        genome_loci=config['run_folder'] + "/{repeat}_genome_loci.fa.gz",
        episome_loci=config['run_folder'] + "/{repeat}_episome_loci.fa.gz"
    output:
        sim_reads=config['run_folder'] + "/{repeat}_simulated_reads.fastq.gz",
        sim_loci=config['run_folder'] + "/{repeat}_all_loci.fa.gz"
    params:
        max_len=config["read_sim"]["max_len"],

        # Badread parameters:
        coverage=config["read_sim"]["coverage"],
        identity=config["badread"]["identity"],
        error_model=config["badread"]["error_model"],
        qscore_model=config["badread"]["qscore_model"],
        start_adapter=config["badread"]["start_adapter"],
        end_adapter=config["badread"]["end_adapter"],
        start_adapter_seq=config["badread"]["start_adapter_seq"],
        end_adapter_seq=config["badread"]["end_adapter_seq"],
        chimeras=config["badread"]["chimeras"],
        glitches=config["badread"]["glitches"],
        junk_reads=config["badread"]["junk_reads"],
        random_reads=config["badread"]["random_reads"],
        vsearch_seed=lambda wildcards: get_repeat_seed(wildcards, "vsearch"),
        badread_seed=lambda wildcards: get_repeat_seed(wildcards, "badread")
    threads: 1
    # 1. Merge the target enriched loci w/wo AAV integration and episomal AAV
    # 2. Suffle loci with VSEARCH
    # 3. gzip merged loci
    # 4. Run badreads to sample reads from these loci
    # 5. [optional] Remove all intermediate files
    shell:
        "gunzip -dc {input.integration_loci} {input.genome_loci} {input.episome_loci} | "
        "vsearch --fasta_width=0 --maxseqlengt={params.max_len} "
        "{params.vsearch_seed} --output /dev/stdout --shuffle - | "
        "gzip > {output.sim_loci} && "
        "badread simulate --reference {output.sim_loci} --quantity {params.coverage}x "
        "--length {params.max_len},0 --identity {params.identity} "
        "--error_model {params.error_model} --qscore_model {params.qscore_model} "
        "--start_adapter {params.start_adapter} --end_adapter {params.end_adapter} "
        "--start_adapter_seq {params.start_adapter_seq} "
        "--end_adapter_seq {params.end_adapter_seq} --chimeras {params.chimeras} "
        "--glitches {params.glitches} --junk_reads {params.junk_reads} "
        "--random_reads {params.random_reads} {params.badread_seed} | gzip > {output.sim_reads}"


# A downsampling rule can be added here if needed.


#########
# Clean #
#########

# Remove intermediate files:
if config["keep_intermediate"]:
    rule clean:
        input:
            sim_reads=config['run_folder'] + "/{repeat}_simulated_reads.fastq.gz",
            sim_loci=config['run_folder'] + "/{repeat}_all_loci.fa.gz",
            integration_sites=config['run_folder'] + "/{repeat}_integration_sites.fa.gz",
            integration_loci=config['run_folder'] + "/{repeat}_integration_loci.fa.gz",
            genome_loci=config['run_folder'] + "/{repeat}_genome_loci.fa.gz",
            episome_loci=config['run_folder'] + "/{repeat}_episome_loci.fa.gz"
        output:
            config['run_folder'] + "/{repeat}_simulation_done"
        shell:
            "touch {output}"
else:
    rule clean:
        input:
            sim_reads=config['run_folder'] + "/{repeat}_simulated_reads.fastq.gz",
            sim_loci=config['run_folder'] + "/{repeat}_all_loci.fa.gz",
            integration_sites=config['run_folder'] + "/{repeat}_integration_sites.fa.gz",
            integration_loci=config['run_folder'] + "/{repeat}_integration_loci.fa.gz",
            genome_loci=config['run_folder'] + "/{repeat}_genome_loci.fa.gz",
            episome_loci=config['run_folder'] + "/{repeat}_episome_loci.fa.gz"
        output:
            config['run_folder'] + "/{repeat}_simulation_done"
        params:
            kde_plot=config['run_folder'] + "/aav_kde_plot.png"
        shell:
            "rm {input.sim_loci} {input.integration_loci} {input.genome_loci} {input.episome_loci} {input.integration_sites} && "
            "rm -f {params.kde_plot} && "
            "touch {output}"






all_lst = list()
for i in range(config["repeats"]):
    repeat = f"rep{i+1}"
    simulation_done = config['run_folder'] + f"/{repeat}_simulation_done"
    all_lst.append(simulation_done)


# Make rule that traverse the whole DAG:
rule all:
    input:
        all_lst


# Create rule DAG:
# snakemake all --snakefile workflow/Snakefile.smk --printshellcmds --verbose --dag all --configfile config/test_config.yaml | dot -Tsvg > dag_test.svg

# Run all rules:
# snakemake all --configfile config/test_config.yaml --snakefile workflow/Snakefile.smk --printshellcmds --verbose

# Run all rules, specifying run name:
# snakemake all --configfile config/test_config.yaml --snakefile workflow/Snakefile.smk --printshellcmds --verbose --config run_name="testrun"





