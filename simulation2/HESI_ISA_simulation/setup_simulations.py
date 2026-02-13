#!/usr/bin/env python3
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
"""
Generate simulation configuration files from template and parameter combinations.

This script reads a template YAML file and a parameters YAML file, then generates
all combinations of the variable parameters to create individual simulation configs.
"""

import argparse
import yaml
import itertools
from pathlib import Path
import pandas as pd
import sys


# Custom list class to force flow style in YAML
class FlowStyleList(list):
    """Custom list class to force flow style in YAML."""
    pass


def represent_flow_style_list(dumper, data):
    """YAML representer for flow-style lists."""
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)


# Register the representer globally
yaml.add_representer(FlowStyleList, represent_flow_style_list)


def str_representer(dumper, data):
    """YAML representer for strings to always use double quotes."""
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='"')


# Register string representer to use double quotes
yaml.add_representer(str, str_representer)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate simulation config files from template and parameters"
    )
    parser.add_argument(
        "--template_yaml",
        required=True,
        help="Template YAML file for simulation config"
    )
    parser.add_argument(
        "--params_yaml",
        required=True,
        help="YAML file with all variable parameters for simulations"
    )
    parser.add_argument(
        "--config_out_dir",
        required=True,
        help="Output directory to store generated simulation config files"
    )
    parser.add_argument(
        "--workdir",
        required=True,
        help="Snakemake workdir."
    )
    parser.add_argument(
        "--run_table_out",
        required=True,
        help="Output table (TSV) with all simulation runs and their parameters"
    )
    parser.add_argument(
        "--seed_space",
        type=int,
        help="Spacing between seeds for each simulation run. If specified, run 1 gets seed=seed_space*1, run 2 gets seed=seed_space*2, etc."
    )
    return parser.parse_args()


def read_yaml(filepath):
    """Read and parse a YAML file."""
    with open(filepath, 'r') as f:
        return yaml.safe_load(f)


def write_yaml(data, filepath):
    """Write data to a YAML file with proper formatting."""
    with open(filepath, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def fill_template(template_str, params):
    """Fill template string with parameter values."""
    filled = template_str
    for key, value in params.items():
        placeholder = f"<{key}>"
        # Handle list values by converting to appropriate format
        if isinstance(value, list):
            value_str = str(value)
        else:
            value_str = str(value)
        filled = filled.replace(placeholder, value_str)
    return filled


def generate_combinations(params):
    """
    Generate all combinations of variable parameters.
    
    Returns:
        list of dicts: Each dict contains one combination of parameters
    """
    # Extract the variable parameters and their values
    genomes = params.get('genomes', [])
    aavs = params.get('AAVs', [])
    is_site_bias = params.get('IS_site_bias', [])
    clonality = params.get('clonality', [])
    
    # Generate all combinations
    combinations = []
    for genome, aav, bias, clon in itertools.product(genomes, aavs, is_site_bias, clonality):
        combo = {
            'genome': genome,
            'AAV': aav,
            'IS_site_bias': bias,
            'clonality': clon
        }
        combinations.append(combo)
    
    return combinations


def check_file_exists(filepath, description, config_dir):
    """
    Check if a file exists. Relative paths are resolved from config_dir.
    
    Args:
        filepath: Path to check (can be relative or absolute)
        description: Description of the file for error messages
        config_dir: Base directory for resolving relative paths
    
    Returns:
        Absolute path if file exists, raises error otherwise
    """
    if not filepath or filepath == "":
        # Empty string is allowed for optional files
        return filepath
    
    path = Path(filepath)
    if not path.is_absolute():
        # Resolve relative to config directory
        path = (config_dir / path).resolve()
    
    if not path.exists():
        raise FileNotFoundError(f"{description} not found: {filepath} (resolved to {path})")
    
    return str(path)


def fill_config_values(template_data, combo, params, config_dir, run_name):
    """
    Fill template configuration with specific parameter values.
    
    Args:
        template_data: The template configuration dict
        combo: Dictionary with the current parameter combination
        params: Full parameters dictionary with mappings
        config_dir: Base directory for resolving relative file paths
        run_name: Unique name for this simulation run
    
    Returns:
        Filled configuration dictionary
    """
    import copy
    config = copy.deepcopy(template_data)
    
    # Fill in the values based on the combination
    genome = combo['genome']
    aav = combo['AAV']
    bias = combo['IS_site_bias']
    clon = combo['clonality']
    
    # Map to actual file paths and values with file existence checks
    genome_file = params['genome_fnam'][genome]
    check_file_exists(genome_file, f"Genome file for {genome}", config_dir)
    config['genome_fnam'] = genome_file
    
    aav_genome_file = params['AAV_genome_fnam'][aav]
    check_file_exists(aav_genome_file, f"AAV genome file for {aav}", config_dir)
    config['AAV_genome_fnam'] = aav_genome_file
    
    tes_panel_file = params['tes_panel_fnam'][aav]
    check_file_exists(tes_panel_file, f"TES panel file for {aav}", config_dir)
    config['tes_panel_fnam'] = tes_panel_file
    
    config['integration']['events'] = params['integration']['events'][clon]
    
    # Convert expansion to FlowStyleList for proper YAML formatting
    expansion_value = params['integration']['expansion'][clon]
    if isinstance(expansion_value, list):
        config['integration']['expansion'] = FlowStyleList(expansion_value)
    else:
        config['integration']['expansion'] = expansion_value
    
    prob_file = params['integration']['aav']['prob_file'][aav]
    check_file_exists(prob_file, f"AAV breakpoint probability file for {aav}", config_dir)
    config['integration']['aav']['prob_file'] = prob_file
    
    empirical_sizes_file = params['integration']['aav']['empirical_sizes'][aav]
    check_file_exists(empirical_sizes_file, f"AAV empirical sizes file for {aav}", config_dir)
    config['integration']['aav']['empirical_sizes'] = empirical_sizes_file
    
    sampling_prob_file = params['integration']['sampling_prob_bed'][bias]
    check_file_exists(sampling_prob_file, f"Integration site sampling probability file for {bias}", config_dir)
    config['integration']['sampling_prob_bed'] = sampling_prob_file
    
    return config


def main():
    """Main execution function."""
    args = parse_args()
    
    # Read input files
    print(f"Reading template from: {args.template_yaml}")
    template_data = read_yaml(args.template_yaml)
    
    print(f"Reading parameters from: {args.params_yaml}")
    params = read_yaml(args.params_yaml)
    
    # Create output directory if it doesn't exist
    config_out_dir = Path(args.config_out_dir)
    config_out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {config_out_dir}")
    
    # Use the config output directory for resolving relative paths
    config_out_dir_resolved = config_out_dir.resolve()
    
    # Get run basename
    run_basename = params.get('run_basename', 'simulation')
    
    # Get number of repeats (use from template or default to 3)
    repeats = template_data.get('repeats', 3)
    
    # Generate all combinations
    combinations = generate_combinations(params)
    print(f"Generated {len(combinations)} simulation combinations")
    print(f"Each simulation will have {repeats} replicates")
    
    # Prepare data for run table
    run_table_data = []
    
    # Process each combination
    for idx, combo in enumerate(combinations, start=1):
        # Create unique run name
        run_name = f"{run_basename}_{idx:03d}"
        
        # Fill configuration with file validation
        try:
            config = fill_config_values(template_data, combo, params, config_out_dir_resolved, 
                                       run_name)
            config['run_name'] = run_name
            config['repeats'] = repeats
            
            # Add seed if seed_space is specified
            if args.seed_space is not None:
                config['seed'] = idx * args.seed_space
        except FileNotFoundError as e:
            print(f"\nError in combination {idx}: {e}")
            print(f"Skipping: {combo}")
            continue
        
        # Define output config file path
        config_filename = f"{run_name}.yaml"
        config_filepath = config_out_dir / config_filename
        
        # Write config file
        write_yaml(config, config_filepath)
        
        # Create one row per replicate in the run table
        for rep_num in range(1, repeats + 1):
            rep_name = f"rep{rep_num}"
            fastq_filename = f"{rep_name}_simulated_reads.fastq.gz"
            
            # Determine workdir: prefer per-config, fall back to template default, then to output dir
            # workdir_val = config.get('workdir') or template_data.get('workdir')
            # print(workdir_val)
            result_dir = str(Path(args.workdir) / run_name)

            # Adjust seed for each replicate if seed is specified (as done in the Snakemake pipeline):
            seed = config.get('seed', '')
            if seed != '':
                seed += rep_num

            run_record = {
                'run_index': idx,
                'run_name': run_name,
                'replicate': rep_name,
                'config_file': str(config_filepath),
                'result_dir': result_dir,
                'fastq_file': fastq_filename,
                'genome': combo['genome'],
                'AAV': combo['AAV'],
                'IS_site_bias': combo['IS_site_bias'],
                'clonality': combo['clonality'],
                'events': config['integration']['events'],
                'expansion': str(config['integration']['expansion']),
                'genome_fnam': config['genome_fnam'],
                'AAV_genome_fnam': config['AAV_genome_fnam'],
                'tes_panel_fnam': config['tes_panel_fnam'],
                'prob_file': config['integration']['aav']['prob_file'],
                'empirical_sizes': config['integration']['aav']['empirical_sizes'],
                'sampling_prob_bed': config['integration']['sampling_prob_bed'],
                'seed': seed
            }
            run_table_data.append(run_record)
        
        print(f"  [{idx}/{len(combinations)}] Created: {config_filename} (with {repeats} replicates)")
    
    # Create and save run table
    run_table = pd.DataFrame(run_table_data)
    # Reorder columns for better readability:
    column_order = [
        'run_index', 'run_name', 'replicate', 'config_file', 'result_dir', 'fastq_file',
        'genome', 'AAV', 'IS_site_bias', 'clonality',
        'genome_fnam', 'AAV_genome_fnam', 'tes_panel_fnam', 
        'prob_file', 'empirical_sizes', 'sampling_prob_bed',
        'expansion', 'events', 'seed'
    ]
    run_table = run_table.loc[:, column_order].copy()

    run_table.to_csv(args.run_table_out, sep='\t', index=False)
    print(f"\nRun table saved to: {args.run_table_out}")
    print(f"Total simulation configurations created: {len(combinations)}")
    print(f"Total simulation runs (including replicates): {len(run_table_data)}")
    print(f"Expected output: {len(run_table_data)} FASTQ files")


if __name__ == "__main__":
    main()
