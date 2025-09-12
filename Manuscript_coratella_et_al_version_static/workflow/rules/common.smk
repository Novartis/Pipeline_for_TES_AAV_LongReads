import os, sys, glob, uuid, json
from pathlib import Path
import pandas as pd




def parse_config():
    # Parse config file #
    # Notice, config is a global and thus is available for all functions.

    # Parse sample sheet:
    ss_df = parse_sample_sheet()

    # Check plasmid sequence basenames:
    check_plasmid_names(ss_df)

    # Enforce certains config keys to be defined:
    required_config_keys()

    # Check file existence:
    check_file_exists()

    # Parse flags settings:
    parse_flags()


    # Check if out_folder is absolute or relative.
    config['working_dir'] = Path(config['working_dir']).resolve()
    out_folder_path = Path(config['out_folder'])
    # If relative add to working_dir to make absolute:
    if not out_folder_path.is_absolute():
        config['out_folder'] = config['working_dir'] / out_folder_path
    else:
        # If absolute keep as is:
        config['out_folder'] = out_folder_path

    # Write config file to JSON:
    os.system(f"mkdir -p {config['out_folder']}")
    config_json = config['out_folder'] / "config.json"
    config['config_json'] = config_json
    with open(config_json, 'w') as fh_out:
        json.dump(make_json_serializable(config), fh_out, indent=4)

    return(ss_df)


def check_plasmid_names(ss_df):
    # Check the plasmid sequence basenames listed in ss_df
    # map to the same path, otherwise "get_plasmid_seq" function
    # will fail.
    plasmid_seq_path = list(ss_df["plasmid_seq"].values)
    plasmid_seq_names = [str(Path(p).stem) for p in plasmid_seq_path]
    if not len(set(plasmid_seq_path)) == len(set(plasmid_seq_names)):
        raise Exception("Plasmid sequence basenames are not uniquely mapping back to their paths. This must be fixed.")

def parse_flags():
    # If flags definition is missing, add an empty one:
    if not "cutadapt_flags" in config:
        config["cutadapt_flags"] = dict()
    if not "sniffles_flags" in config:
        config["sniffles_flags"] = dict()
    
    # The flag string for sniffles cannot contain double quote:
    sniffles_flags_str = ' '.join(config["sniffles_flags"].values())
    if '"' in sniffles_flags_str:
        raise Exception(f"Double quotes not allowed in sniffles flags: {sniffles_flags_str}")

def required_config_keys():
    required_keys = ["ref_genome", "ref_mmi"]
    for k in required_keys:
        if not k in config:
            raise Exception(f"Key missing in config file: {k}")

def check_file_exists():
    file_keys = ["ref_genome", "ref_mmi"]
    for k in file_keys:
        if not os.path.exists(config[k]):
            raise Exception(f"Filepath ({config[k]}) not found for {k} key in config file.")

def parse_sample_sheet():
    # Check that samplesheet is defined in config file:
    if not "samplesheet" in config:
        raise Exception("The \"samplesheet\" key must be defined in the config file. This should contain a tab separated sample sheet.")
    
    # Define required columns and validity tests:
    required_ss_cols = {"sample_name": {},
                        "fastq_path":  {"test": [[lambda x: os.path.exists(x), "File missing"], ]},
                        "plasmid_seq": {"test": [[lambda x: os.path.exists(x), "File missing"], ]},
                        "adapter5":{},
                        "adapter3":{}}

    # Read and validate sample sheet:
    ss_df = pd.read_csv(config["samplesheet"], sep="\t")
    for col in required_ss_cols.keys():
        if not col in ss_df:
            raise Exception(f"Column name \"{col}\" missing in sample sheet.")
        if "test" in required_ss_cols[col]:
            for ci, cell in enumerate(ss_df[col]):
                for test in required_ss_cols[col]["test"]:
                    if not test[0](cell):
                        raise Exception(f"Failed to validate sample sheet value in column \"{col}\" and row {ci+0}: {test[0]}")

    # Enforce unique sample names:
    if not ss_df["sample_name"].is_unique:
        raise Exception(f"Failed to validate sample sheet. Values in the \"sample_name\" column must be unique.")

    # No adapter is empty string:
    ss_df["adapter5"] = ss_df["adapter5"].fillna('')
    ss_df["adapter3"] = ss_df["adapter3"].fillna('')

    return(ss_df)



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


def get_ss_samples():
    return(ss_df.loc[:, ["sample_name", "fastq_path", "plasmid_seq", "adapter5", "adapter3"]].values)

def get_ss_fastq(wildcards):
    mask = (ss_df["sample_name"] == wildcards.sample)
    return(ss_df[mask]["fastq_path"].values[0])

def get_ss_plasmid(wildcards):
    mask = (ss_df["sample_name"] == wildcards.sample)
    return(ss_df[mask]["plasmid_seq"].values[0])

def get_ss_adapters(wildcards):
    mask = (ss_df["sample_name"] == wildcards.sample)
    return([ss_df[mask]["adapter5"].values[0], ss_df[mask]["adapter3"].values[0]])

def get_plasmid_name(wildcards):
    # Extract name of plasmid:
    mask = (ss_df["sample_name"] == wildcards.sample)
    plasmid_seq_path = Path(ss_df[mask]["plasmid_seq"].values[0])
    plasmid_seq_name = plasmid_seq_path.stem
    return(plasmid_seq_name)

def get_plasmid_seq(wildcards):
    # Given the plasmid name, get the filepath to its sequence:
    mask = (ss_df["plasmid_seq"].apply(lambda x: Path(x).stem == wildcards.plasmid))
    plasmid_seq_paths = ss_df[mask]["plasmid_seq"].values
    assert(len(set(plasmid_seq_paths)) == 1)
    return(plasmid_seq_paths[0])




