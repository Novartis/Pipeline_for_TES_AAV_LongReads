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

import os, sys, glob, uuid, json
import pandas as pd




def parse_config():
    # Convert all paths to absolute paths:
    make_abs_path()
    # Check if all files exist:
    check_file_exists()

    # Parse config file to read simulation length specification:
    config["read_sim"]["mean_frag_len"], \
      config["read_sim"]["stdev_frag_len"], \
      config["read_sim"]["max_len"] = \
      list(map(int, config["read_sim"]["length"].split(',')))

    # Similarly for AAV length specification:
    aav_len_params = list(map(int, config["integration"]["size"].split(',')))
    config["integration"]["mean_aav_len"] = aav_len_params[0]
    config["integration"]["stdev_aav_len"] = aav_len_params[1]
    if len(aav_len_params) == 3:
        config["integration"]["min_aav_len"] = aav_len_params[2]
    else:
        config["integration"]["min_aav_len"] = 25


    if config["enrichment"]["additive"]:
        config["enrichment"]["additive_flag"] = "--add"
    else:
        config["enrichment"]["additive_flag"] = ''


    # Extract the minimum length of enrichment:
    config["enrichment"]["min_len"] = int(config["enrichment"]["enrich_prob"].split(',')[1])
    if config["enrichment"]["min_len"] < 14:
        raise Exception("The minimum length of complementarity for target enrichment must be 14 or more.") 


    # minimap kmer length:
    config["minimap"] = dict()
    if config["enrichment"]["min_len"] > 21:
        config["minimap"]["k"] = 21
    else:
        config["minimap"]["k"] = config["enrichment"]["min_len"]


    # The number of sites drawn from the genome.
    # This number needs to be larger than the
    # number of chromosomes in the reference:
    if config["integration"]["events"] < 1000:
        config["integration"]["site_draws"] = 1000
    else:
        config["integration"]["site_draws"] = config["integration"]["events"]

    # Create random name for a run folder:
    if not "run_name" in config:
        config["run_name"] = str(uuid.uuid4())
    config["run_folder"] = f"{config['workdir']}/{config['run_name']}"

    # Write config file to JSON:
    config_json = f"{config['run_folder']}/config.json"
    config_json_dir = '/'.join(config_json.split('/')[:-1])
    config['config_json'] = config_json
    os.system(f"mkdir -p {config_json_dir}")
    with open(config_json, 'w') as fh_out:
        json.dump(make_json_serializable(config), fh_out, indent=4)


# Convert all values to strings if they are not JSON serializable:
def make_json_serializable(obj):
    if isinstance(obj, dict):
        return({k: make_json_serializable(v) for k, v in obj.items()})
    elif isinstance(obj, list):
        return([make_json_serializable(i) for i in obj])
    elif isinstance(obj, (Path,)):
        return(str(obj))
    else:
        return(obj)


def make_abs_path():
    file_keys = ["genome_fnam", "AAV_genome_fnam", "tes_panel_fnam"]
    for k in file_keys:
        config[k] = os.path.abspath(config[k])


def check_file_exists():
    file_keys = ["genome_fnam", "AAV_genome_fnam", "tes_panel_fnam"]
    for k in file_keys:
        if not os.path.exists(config[k]):
            raise Exception(f"Filepath ({config[k]}) not found for {k} key in config file.")




