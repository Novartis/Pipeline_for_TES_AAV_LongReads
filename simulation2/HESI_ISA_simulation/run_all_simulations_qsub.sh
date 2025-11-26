#!/bin/bash
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

### Bash script to start all simulations for HESI ISA read simulation ###

# Glob pattern for config files:
CONFIG_GLOB="config_files/HESI_ISA_simulation_*"

# UGE job parameters:
N_SMK_CORES=5
MEMORY="10G"
RUNTIME=36000

DEFAULT_CONFIG="config_files/default_config.yaml"

# Dry-run toggle: set to true to run snakemake in dry-run mode:
# false = normal run mode
# true/1 = prints commands, no snakemake execution
# 2 = dry-run snakemake execution
DRY_RUN=false

# Prepare extra snakemake args (no environment variables)
if [ "$DRY_RUN" = "2" ]; then
  SNAKE_EXTRA=(--dryrun)
else
  SNAKE_EXTRA=()
fi



# Setup all the simulation defined in the all_simulations_params.yaml.
# Thus script will generate all config files and a run table for all simulations.
python setup_simulations.py \
  --params_yaml all_simulations_params.yaml \
  --template_yaml config_files/all_sim_template.yaml \
  --seed_space 100 \
  --config_out_dir config_files \
  --workdir simulation_output \
  --run_table_out all_simulations_runtable.tsv



# Prepare to start simulation runs:
mkdir -p logs
mkdir -p simulation_output
shopt -s nullglob
configs=($CONFIG_GLOB)
if [ ${#configs[@]} -eq 0 ]; then
  echo "No config files found matching $CONFIG_GLOB" >&2
  exit 1
fi


for config_file in "${configs[@]}"; do
  # Derive logfile name from config file basename without any extensions:
  base="$(basename "$config_file")"
  base_no_ext="${base%%.*}"
  logfile="logs/${base_no_ext}.log"

  # Prepare snakemake command:
  cmd="snakemake all"
  [ ${#SNAKE_EXTRA[@]} -gt 0 ] && cmd="$cmd ${SNAKE_EXTRA[*]}"
  cmd="$cmd --config default_configfile=${DEFAULT_CONFIG}"
  cmd="$cmd --configfile $config_file"
  cmd="$cmd --snakefile ../workflow/Snakefile.smk"
  cmd="$cmd --cores $((2 * N_SMK_CORES))"
  cmd="$cmd --printshellcmds --verbose"

  # If dry-run mode true/1, just print the command that would be run:
  if [ "$DRY_RUN" = true ] || [ "$DRY_RUN" = "1" ]; then
    echo "Dry-run mode: would submit job for $config_file with command:"
    echo "qsub -pe smp $N_SMK_CORES -l m_mem_free=$MEMORY,h_rt=$RUNTIME -S /bin/bash -o $logfile -j y -cwd -V -b y $cmd"
    continue
  fi

  # Submit the simulation as a UGE job:
  echo "Submitting job: $config_file -> $logfile"
  qsub -pe smp "$N_SMK_CORES" -l "m_mem_free=$MEMORY,h_rt=$RUNTIME" -S /bin/bash -o "$logfile" -j y -cwd -V -b y $cmd
done




