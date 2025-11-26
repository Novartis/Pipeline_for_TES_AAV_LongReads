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
Add MD5 checksums for FASTQ files to the simulation run table.

This script reads a run table TSV file, calculates MD5 checksums for each
FASTQ file (decompressing gzipped files on the fly), and adds the checksums
as a new column.
"""

import argparse
import pandas as pd
import hashlib
import gzip
from pathlib import Path
import sys


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Add MD5 checksums for FASTQ files to run table"
    )
    parser.add_argument(
        "--run_table_in",
        required=True,
        help="Input TSV file with simulation runs (e.g., all_simulations_runtable.tsv)"
    )
    parser.add_argument(
        "--md5_colname",
        default="fastq_md5",
        help="Name for the MD5 checksum column (default: fastq_md5)"
    )
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=8192,
        help="Chunk size for reading files in bytes (default: 8192)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print progress information"
    )
    return parser.parse_args()


def calculate_md5_gzip(filepath, chunk_size=8192):
    """
    Calculate MD5 checksum of a gzipped file after decompression.
    
    Args:
        filepath: Path to the gzipped file
        chunk_size: Size of chunks to read at a time (bytes)
    
    Returns:
        MD5 hex digest string, or None if file doesn't exist
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        return None
    
    md5_hash = hashlib.md5()
    
    try:
        with gzip.open(filepath, 'rb') as f:
            # Read file in chunks to handle large files efficiently
            while True:
                chunk = f.read(chunk_size)
                if not chunk:
                    break
                md5_hash.update(chunk)
        
        return md5_hash.hexdigest()
    
    except Exception as e:
        print(f"Error calculating MD5 for {filepath}: {e}", file=sys.stderr)
        return None


def add_md5_checksums(run_table_path, md5_colname, chunk_size=8192, verbose=False):
    """
    Add MD5 checksums to the run table.
    
    Args:
        run_table_path: Path to the input TSV file
        md5_colname: Name for the MD5 column
        chunk_size: Size of chunks for reading files
        verbose: Whether to print progress
    
    Returns:
        Updated DataFrame
    """
    # Read the run table
    if verbose:
        print(f"Reading run table from: {run_table_path}")
    
    df = pd.read_csv(run_table_path, sep='\t')
    
    if verbose:
        print(f"Found {len(df)} rows in run table")
    
    # Check if required columns exist
    if 'result_dir' not in df.columns:
        raise ValueError("Column 'result_dir' not found in run table")
    if 'fastq_file' not in df.columns:
        raise ValueError("Column 'fastq_file' not found in run table")
    
    # Initialize MD5 column
    md5_values = []
    
    # Calculate MD5 for each FASTQ file
    for idx, row in df.iterrows():
        result_dir = row['result_dir']
        fastq_file = row['fastq_file']
        
        # Construct full path to FASTQ file
        fastq_path = Path(result_dir) / fastq_file
        
        if verbose:
            print(f"[{idx + 1}/{len(df)}] Processing: {fastq_path}")
        
        # Calculate MD5
        md5sum = calculate_md5_gzip(fastq_path, chunk_size)
        
        if md5sum is None:
            if verbose:
                print(f"  → File not found or error, leaving empty")
            md5_values.append('')
        else:
            if verbose:
                print(f"  → MD5: {md5sum}")
            md5_values.append(md5sum)
    
    # Add MD5 column to dataframe
    df[md5_colname] = md5_values
    
    return df


def main():
    """Main execution function."""
    args = parse_args()
    
    # Check if input file exists
    if not Path(args.run_table_in).exists():
        print(f"Error: Input file not found: {args.run_table_in}", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Add MD5 checksums
        df = add_md5_checksums(
            args.run_table_in,
            args.md5_colname,
            args.chunk_size,
            args.verbose
        )
        
        # Save updated table (overwrite original)
        if args.verbose:
            print(f"\nSaving updated run table to: {args.run_table_in}")
        
        df.to_csv(args.run_table_in, sep='\t', index=False)
        
        # Print summary
        total = len(df)
        filled = (df[args.md5_colname] != '').sum()
        missing = total - filled
        
        print(f"\nSummary:")
        print(f"  Total rows: {total}")
        print(f"  MD5 calculated: {filled}")
        print(f"  Files not found: {missing}")
        print(f"  Column added: '{args.md5_colname}'")
        print(f"\nRun table updated successfully!")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
