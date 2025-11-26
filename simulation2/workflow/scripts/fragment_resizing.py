#!/usr/bin/env python
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

import sys
import os
import argparse
import random
import numpy as np
from badread_dependencies import FragmentLengths


def process_fastq_entry(header, seq, plus, qual, frag_lengths):
    """
    Process a single FASTQ entry and trim to gamma-distributed fragment length.
    
    Args:
        header: FASTQ header line (without @)
        seq: Sequence line
        plus: Plus line (should be '+')
        qual: Quality score line
        frag_lengths: FragmentLengths object for sampling fragment lengths
    
    Returns:
        str: Formatted FASTQ entry (4 lines joined with newlines), or None to skip
    """
    fragment_length = frag_lengths.get_fragment_length()
    
    # If fragment is longer than input, keep as-is
    if len(seq) <= fragment_length:
        return f"{header}\n{seq}\n{plus}\n{qual}"
    
    # Calculate trim length and randomly choose left or right trim
    trim_len = len(seq) - fragment_length
    header_parts = header.split('_')
    
    # Keep only first 3 parts (chrom_start_end)
    header_parts = header_parts[0:3]
    
    # Randomly trim from left or right (Bernoulli sampling)
    if np.random.binomial(1, 0.5):
        # Trim from left
        header_parts[1] = str(int(header_parts[1]) + trim_len)
        new_header = '_'.join(header_parts)
        return f"{new_header}\n{seq[trim_len:]}\n{plus}\n{qual[trim_len:]}"
    else:
        # Trim from right
        header_parts[2] = str(int(header_parts[2]) - trim_len)
        new_header = '_'.join(header_parts)
        return f"{new_header}\n{seq[:-trim_len]}\n{plus}\n{qual[:-trim_len]}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Takes a stream of FASTQ lines from wgsim and trims them so the
                       fragment lengths become gamma distributed."""
    )

    # Input options:
    parser.add_argument("--mfl", default=9000, type=int,
        help="Mean fragment length.")

    parser.add_argument("--sfl", default=3000, type=int,
        help="Standard deviation of fragment length.")

    parser.add_argument("--verbose", action='store_true')

    parser.add_argument("--seed", default=None, type=int,
        help="Random seed for reproducibility.")

    # Parse the arguments:
    args = parser.parse_args()
    
    # Set random seed for reproducibility
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    if args.verbose:
        screen_output = sys.stderr
        print(args, file=screen_output)
    else:
        screen_output = open(os.devnull, 'w')

    # Initiate the fragment length objective:
    frag_lengths = FragmentLengths(args.mfl, args.sfl, screen_output)

    # Read FASTQ entries from stdin until EOF:
    while True:
        header = sys.stdin.readline().strip()
        if not header:
            break
        seq = sys.stdin.readline().strip()
        plus = sys.stdin.readline().strip()
        qual = sys.stdin.readline().strip()

        # Process and print the trimmed entry
        result = process_fastq_entry(header, seq, plus, qual, frag_lengths)
        if result:
            print(result)



