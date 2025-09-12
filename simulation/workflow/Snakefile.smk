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


### The pipeline repo has the follow folder structure:
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

configfile: pipeline_dir / "config/default_config.yaml"
workdir: config['workdir']
config['workdir'] = os.getcwd()
include: pipeline_dir / "workflow/rules/common.smk"
config['pipeline_dir'] = pipeline_dir
config['scripts_dir'] = pipeline_dir / "workflow/scripts"
parse_config()






# The minimap2 alignment to the AAV reference should
# probably be forward only. 





###############
# Input check #
###############

# Check to see that the reference genome conforms with
# the requirements of the downstream programs.
# This should produce an empty file, if not the genome
# has chromosomes too short:
empty_file_sha256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
rule check_reference:
    output:
        ensure(config['run_folder'] + "/reference_check", sha256=empty_file_sha256)
    params:
        max_len=11 * config["read_sim"]["max_len"],
        genome_fnam=config["genome_fnam"],
        scrip_dir=config['pipeline_dir'] / "workflow/scripts"
    threads: 1
    shell:
        "python {params.scrip_dir}/ref_length.py --max_len {params.max_len} "
        "--fa {params.genome_fnam} > {output}"




###################################
# Simulate AAV integration events #
###################################

# Generate the loci for integration:
rule draw_integration_sites:
    input:
        genome_check=config['run_folder'] + "/reference_check"
    output:
        integration_sites=config['run_folder'] + "/{repeat}_integration_sites.fa.gz"
    params:
        N=config["integration"]["site_draws"],
        events=config["integration"]["events"],
        max_len=10 * config["read_sim"]["max_len"],
        genome_fnam=config["genome_fnam"]
    threads: 1
    # 1. Use wgsim to sample long (10x max read length) fragments of the genome
    # 2. Convert to fragments to fasta witj seqkit
    # 3. Shuffle the order of the fragments with vsearch and return one fragment
    #    per AAV integration event
    # 4. gzip output fasta containing the integration site fragments
    shell:
        "wgsim -e0 -d0 -s0 -r0 -R0 -X0 -N {params.N} -1 {params.max_len} "
        "{params.genome_fnam} /dev/stdout /dev/null | "
        "seqkit fq2fa | vsearch --fasta_width=0 --maxseqlengt={params.max_len} "
        "--topn={params.events} --output /dev/stdout --shuffle - | "
        "gzip > {output.integration_sites}"


# Add AAV to integration sites and simulate the reads:
rule sim_AAV_loci:
    input:
        integration_sites=config['run_folder'] + "/{repeat}_integration_sites.fa.gz"
    output:
        events_tsv=config['run_folder'] + "/{repeat}_integration_events.tsv",
        sim_reads=config['run_folder'] + "/{repeat}_integration_loci.fa.gz"
    params:
        # For "add_integration.py":
        mean_frag_len=config["read_sim"]["mean_frag_len"],
        stdev_frag_len=config["read_sim"]["stdev_frag_len"],
        mean_aav_len=config["integration"]["mean_aav_len"],
        stdev_aav_len=config["integration"]["stdev_aav_len"],
        min_aav_len=config["integration"]["min_aav_len"],
        expansion=config["integration"]["expansion"],
        AAV_genome_fnam=config["AAV_genome_fnam"],

        # For minimap/target enrichment:
        k=config["minimap"]["k"],
        enrich_min_len=config["enrichment"]["min_len"],
        tes_panel_fnam=config["tes_panel_fnam"],
        enrich_prob=config["enrichment"]["enrich_prob"],
        amplification=config["enrichment"]["amplification"],
        leakage=config["enrichment"]["leakage"],
        additive=config["enrichment"]["additive_flag"],

        scrip_dir=config['pipeline_dir'] / "workflow/scripts"
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
        "python {params.scrip_dir}/add_integration.py --fq_out "
        "--mfl {params.mean_frag_len} --sfl {params.stdev_frag_len} "
        "--mal {params.mean_aav_len} --sal {params.stdev_aav_len} "
        "--min_aav {params.min_aav_len} --exp {params.expansion} "
        "--aav {params.AAV_genome_fnam} --int_out {output.events_tsv} | "
        "minimap2 -t1 -k{params.k} -w2 --sr --frag=yes "
        "-A1 -B120 -O20,40 -E40,25 -b0 -r100 -p0.5 -N100 -f1000,5000 "
        "-n2 -m{params.enrich_min_len} -s{params.enrich_min_len} "
        "--secondary=no -a {params.tes_panel_fnam} - | "
        "python {params.scrip_dir}/target_enrichment.py --ep {params.enrich_prob} "
        "--amp {params.amplification} --leak {params.leakage} {params.additive} | "
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
        enrich_prob=config["enrichment"]["enrich_prob"],
        amplification=config["enrichment"]["amplification"],
        leakage=config["enrichment"]["leakage"],
        additive=config["enrichment"]["additive_flag"],

        scrip_dir=config['pipeline_dir'] / "workflow/scripts"
    threads: 5
    # 1. Use wgsim to sample same size fragments of the genome
    # 2. Use fragment_resizing.py to randomly resize these fragment
    #    such that the resulting fragment size follow a gamma distriubtion
    # 3. Align these fastq fragment to the target enrichment panel with minimap2
    # 4. Simulate target enrichment and amplification of the fragments
    #    using target_enrichment.py
    # 5. gzip output fasta containing the target enriched loci
    shell:
        "wgsim -e0 -d0 -s0 -r0 -R0 -X0 -N {params.N} -1 {params.max_len} "
        "{params.genome_fnam} /dev/stdout /dev/null | "
        "python {params.scrip_dir}/fragment_resizing.py "
        "--mfl {params.mean_frag_len} --sfl {params.stdev_frag_len} | "
        "minimap2 -t4 -k{params.k} -w2 --sr --frag=yes "
        "-A1 -B120 -O20,40 -E40,25 -b0 -r100 -p0.5 -N100 -f1000,5000 "
        "-n2 -m{params.enrich_min_len} -s{params.enrich_min_len} "
        "--secondary=no -a {params.tes_panel_fnam} - | "
        "python {params.scrip_dir}/target_enrichment.py --ep {params.enrich_prob} "
        "--amp {params.amplification} --leak {params.leakage} {params.additive} | "
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
        genome_loci=config['run_folder'] + "/{repeat}_genome_loci.fa.gz"
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
        random_reads=config["badread"]["random_reads"]
    threads: 1
    # 1. Merge the target enriched loci w/wo AAV integration
    # 2. Suffle loci with VSEARCH
    # 3. gzip merged loci
    # 4. Run badreads to sample reads from these loci
    # 5. [optional] Remove all intermediate files
    shell:
        "gunzip -dc {input.integration_loci} {input.genome_loci} | "
        "vsearch --fasta_width=0 --maxseqlengt={params.max_len} "
        "--output /dev/stdout --shuffle - | "
        "gzip > {output.sim_loci} && "
        "badread simulate --reference {output.sim_loci} --quantity {params.coverage}x "
        "--length {params.max_len},0 --identity {params.identity} "
        "--error_model {params.error_model} --qscore_model {params.qscore_model} "
        "--start_adapter {params.start_adapter} --end_adapter {params.end_adapter} "
        "--start_adapter_seq {params.start_adapter_seq} "
        "--end_adapter_seq {params.end_adapter_seq} --chimeras {params.chimeras} "
        "--glitches {params.glitches} --junk_reads {params.junk_reads} "
        "--random_reads {params.random_reads} | gzip > {output.sim_reads}"




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
            genome_loci=config['run_folder'] + "/{repeat}_genome_loci.fa.gz"
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
            genome_loci=config['run_folder'] + "/{repeat}_genome_loci.fa.gz"
        output:
            config['run_folder'] + "/{repeat}_simulation_done"
        shell:
            "rm {input.sim_loci} {input.integration_loci} {input.genome_loci} {input.integration_sites} && "
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
# snakemake --dag all --configfile config/test_config.yaml | dot -Tsvg > dag_test.svg

# Run all rules:
# snakemake all --configfile config/test_config.yaml

# Run all rules, specifying run name:
# snakemake all --configfile config/test_config.yaml --config run_name="testrun"





