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

# Change to run more or less simulations in parallel:
N_CONCURRENT_RUNS=10
# Snakemake cores per simulation:
N_SMK_CORES=20

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
  # wait while there are already N_CONCURRENT_RUNS background jobs
  while [ "$(jobs -rp | wc -l)" -ge "$N_CONCURRENT_RUNS" ]; do
    sleep 10
  done

  # Derive logfile name from config file basename without any extensions:
  base="$(basename "$config_file")"
  base_no_ext="${base%%.*}"
  logfile="logs/${base_no_ext}.log"

  # Prepare snakemake command:
  cmd=(
    snakemake
    all
    "${SNAKE_EXTRA[@]}"
    --config
    "default_configfile=${DEFAULT_CONFIG}"
    --configfile
    "$config_file"
    --snakefile
    ../workflow/Snakefile.smk
    --cores
    "$N_SMK_CORES"
    --printshellcmds
    --verbose
  )

  # If dry-run mode true/1, just print the command that would be run:
  if [ "$DRY_RUN" = true ] || [ "$DRY_RUN" = "1" ]; then
    echo "Dry-run mode: would start simulation with command:"
    echo "${cmd[*]} > $logfile 2>&1 &"
    continue
  fi

  # Start the simulation in the background, redirecting output to logfile:
  echo "Starting: $config_file -> $logfile"
  echo "Running: ${cmd[*]} > $logfile 2>&1 &"
  echo "Running: ${cmd[*]} > $logfile 2>&1 &" > "$logfile"
  "${cmd[@]}" >> "$logfile" 2>&1 &
done
# wait for remaining background jobs to finish
wait




