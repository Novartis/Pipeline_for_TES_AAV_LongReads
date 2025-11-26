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
from pathlib import Path
import pandas as pd

# Add scripts directory to path to import bed_validator
scripts_dir = os.path.join(os.path.dirname(workflow.snakefile), '..', 'scripts')
sys.path.insert(0, os.path.abspath(scripts_dir))
from bed_validator import (
    validate_integration_probability_file,
    validate_aav_breakpoint_probability_file,
    BEDValidationError
)




def parse_config():
    # Convert all paths to absolute paths:
    make_abs_path()
    # Check if all files exist:
    check_file_exists()
    # Validate BED probability files:
    validate_probability_files()

    # Validate and normalize expansion parameter
    expansion = config["integration"]["expansion"]
    if isinstance(expansion, int):
        # Single integer: convert to list for consistency
        config["integration"]["expansion"] = [expansion]
    elif isinstance(expansion, list):
        # Validate list length and contents
        if len(expansion) == 0:
            raise ValueError("expansion must have at least 1 value")
        if len(expansion) > 2:
            raise ValueError(f"expansion can have at most 2 values, got {len(expansion)}: {expansion}")
        if not all(isinstance(x, int) for x in expansion):
            raise ValueError(f"All expansion values must be integers, got: {expansion}")
    else:
        raise ValueError(f"expansion must be an integer or list of integers, got {type(expansion)}: {expansion}")

    # Parse config file to read simulation length specification:
    config["read_sim"]["mean_frag_len"], \
      config["read_sim"]["stdev_frag_len"], \
      config["read_sim"]["max_len"] = \
      list(map(int, config["read_sim"]["length"].split(',')))

    # Check if genomic integration fragment length exceeds 1e6 bp
    genomic_frag_len = 10 * config["read_sim"]["max_len"] + 2 * config["integration"]["deletion"]["max_size"]
    if genomic_frag_len > 1e6:
        print(f"WARNING: Genomic integration fragment length ({genomic_frag_len} bp) exceeds 1e6 bp.", file=sys.stderr)
        print(f"         This results no integrations in large parts of the beginning and end of each chromosome.", file=sys.stderr)
        print(f"         Fragment length = 10 × max_read_length ({config['read_sim']['max_len']}) + 2 × max_deletion_size ({config['integration']['deletion']['max_size']})", file=sys.stderr)

    # Extract the minimum length of enrichment:
    if config["enrichment"]["min_len"] < 14:
        raise Exception("The minimum length of complementarity for target enrichment must be 14 or more.") 


    # minimap kmer length:
    config["minimap"] = dict()
    if config["enrichment"]["min_len"] > 21:
        config["minimap"]["k"] = 21
    else:
        config["minimap"]["k"] = config["enrichment"]["min_len"]

    # Set up seed-related parameters for reproducibility:
    if config.get("seed") is not None:
        config["seed_arg"] = f"--seed {config['seed']}"
        config["wgsim_seed"] = f"-S {config['seed']}"
        config["vsearch_seed"] = f"--randseed {config['seed']}"
        config["minimap2_seed"] = f"--seed {config['seed']}"
        config["badread_seed"] = f"--seed {config['seed']}"
        config["base_seed"] = config["seed"]
    else:
        config["seed_arg"] = ""
        config["wgsim_seed"] = ""
        config["vsearch_seed"] = ""
        config["minimap2_seed"] = ""
        config["badread_seed"] = ""
        config["base_seed"] = None

    # Create random name for a run folder:
    if not "run_name" in config:
        config["run_name"] = "sim_run-" + str(uuid.uuid4())
    config["run_folder"] = f"{config['workdir']}/{config['run_name']}"

    # Write config file to JSON:
    config_json = f"{config['run_folder']}/config.json"
    config_json_dir = '/'.join(config_json.split('/')[:-1])
    config['config_json'] = config_json
    os.system(f"mkdir -p {config_json_dir}")
    with open(config_json, 'w') as fh_out:
        json.dump(make_json_serializable(config), fh_out, indent=4)


def get_repeat_seed(wildcards, tool=""):
    """Generate repeat-specific seed by extracting repeat number and adding to base seed."""
    if config.get("base_seed") is None:
        return ""
    
    # Extract repeat number from wildcard (format: "rep1", "rep2", etc.)
    repeat_num = int(wildcards.repeat.replace("rep", ""))
    repeat_seed = config["base_seed"] + repeat_num - 1
    
    # Return tool-specific seed format
    if tool == "wgsim":
        return f"-S {repeat_seed}"
    elif tool == "vsearch":
        return f"--randseed {repeat_seed}"
    elif tool == "minimap2" or tool == "badread" or tool == "":
        return f"--seed {repeat_seed}"
    else:
        return f"--seed {repeat_seed}"


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
    """Recursively walk the config object and convert any path-like entries
    (strings or Path objects) to absolute paths, expanding user (~) and
    resolving relative segments. We keep original non-path scalars unchanged.

    A value is considered path-like if:
      - It is a Path instance, OR
      - It is a string containing a path separator, starts with '.', '~', or '/',
        or points to an existing file/dir.
    """

    def _normalize(v):
        # Path objects
        if isinstance(v, Path):
            return v.expanduser().resolve()
        # Strings that might be paths
        if isinstance(v, str):
            looks_like_path = (
                os.sep in v or
                v.startswith('.') or
                v.startswith('~') or
                v.startswith('/') or
                os.path.exists(os.path.expanduser(v))
            )
            if looks_like_path:
                # expand ~ first
                expanded = os.path.expanduser(v)
                # If empty after expansion, return as-is
                if expanded == '':
                    return v
                return Path(expanded).resolve()
        # Lists / tuples
        if isinstance(v, list):
            return [ _normalize(i) for i in v ]
        if isinstance(v, tuple):
            return tuple(_normalize(i) for i in v)
        # Dicts
        if isinstance(v, dict):
            return { kk: _normalize(vv) for kk, vv in v.items() }
        return v

    # Don't normalize certain config keys that are not file paths
    keys_to_skip = {'run_name', 'workdir'}
    
    for k in list(config.keys()):
        if k in keys_to_skip:
            continue
        config[k] = _normalize(config[k])

    # Maintain backwards compatibility: ensure primary expected file keys are absolute strings
    file_keys = ["genome_fnam", "AAV_genome_fnam", "tes_panel_fnam"]
    for k in file_keys:
        if k in config:
            # Convert Path back to str for downstream tools expecting string
            if isinstance(config[k], Path):
                config[k] = str(config[k])


def check_file_exists():
    file_keys = ["genome_fnam", "AAV_genome_fnam", "tes_panel_fnam"]
    for k in file_keys:
        if not os.path.exists(config[k]):
            raise Exception(f"Filepath ({config[k]}) not found for {k} key in config file.")


def validate_probability_files():
    """
    Validate BED-like probability files for integration sites and AAV breakpoints.
    Checks file format, ensures start < end, and detects overlapping regions.
    """
    try:
        # Validate integration site probability file if specified
        if "integration" in config and "sampling_prob_bed" in config["integration"]:
            prob_file = config["integration"]["sampling_prob_bed"]
            if prob_file and prob_file != "":
                validate_integration_probability_file(prob_file)
        
        # Validate AAV breakpoint probability file if specified
        if "integration" in config and "aav" in config["integration"] and "prob_file" in config["integration"]["aav"]:
            prob_file = config["integration"]["aav"]["prob_file"]
            if prob_file and prob_file != "":
                validate_aav_breakpoint_probability_file(prob_file)
    
    except BEDValidationError as e:
        raise Exception(f"BED probability file validation failed: {str(e)}")





