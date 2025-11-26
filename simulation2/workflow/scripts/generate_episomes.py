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
import argparse
import random
import numpy as np
import os
from Bio import SeqIO
from Bio.Seq import Seq
from aav_fragment_generator import AAVFragmentGenerator


def generate_episome_fragment(aav_fragment_gen, event_num):
    """
    Generate a single episomal AAV fragment.
    
    Args:
        aav_fragment_gen: AAVFragmentGenerator instance
        event_num: Event number for naming (1-based)
    
    Returns:
        str: FASTQ entry (4 lines joined with newlines)
    """
    # Sample AAV fragment using the configured generation method
    aav_frag, aav_start, aav_end, aav_len = aav_fragment_gen.generate()
    
    # Randomly decide if AAV fragment should be reverse complemented
    if random.choice([True, False]):
        aav_frag = str(Seq(aav_frag).reverse_complement())
    
    # Output in FASTQ format (required for minimap2)
    header = f"episome_{event_num}"
    frag_qual = '~' * len(aav_frag)
    return f"@{header}\n{aav_frag}\n+\n{frag_qual}"


def main():
    """Generate episomal AAV DNA fragments for simulation."""
    parser = argparse.ArgumentParser(
        description="Generate episomal AAV DNA fragments for simulation."
    )
    
    parser.add_argument("--aav", required=True, type=str,
        help="AAV reference genome file (FASTA).")
    
    parser.add_argument("--aav_method", required=True, type=str,
        help="AAV fragment generation method string.")
    
    parser.add_argument("--aav_prob_file", default="", type=str,
        help="Optional BED file with AAV breakpoint probabilities.")
    
    parser.add_argument("--aav_empirical_sizes", default="", type=str,
        help="Optional file with empirical fragment sizes for empirical_kde method (one integer per line).")
    
    parser.add_argument("--events", required=True, type=int,
        help="Number of episomal AAV fragments to generate.")
    
    parser.add_argument("--verbose", action='store_true',
        help="Enable verbose output.")
    
    parser.add_argument("--seed", default=None, type=int,
        help="Random seed for reproducibility.")
    
    args = parser.parse_args()
    
    # Set random seed for reproducibility
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)
    
    # Set up output streams
    if args.verbose:
        screen_output = sys.stderr
        print(args, file=screen_output)
    else:
        screen_output = open(os.devnull, 'w')
    
    # Read in the AAV sequence
    aav_seq = [str(record.seq).upper() for record in SeqIO.parse(args.aav, "fasta")]
    if len(aav_seq) > 1:
        print(f"Warning: Multiple sequences found for the AAV: {args.aav}\nOnly the first sequence will be used.",
              file=screen_output)
    aav_seq = aav_seq[0]
    
    # Initialize the AAV fragment generator
    print("AAV fragment generation method:", args.aav_method, file=screen_output)
    if args.aav_prob_file:
        print(f"AAV breakpoint probabilities file: {args.aav_prob_file}", file=screen_output)
    if args.aav_empirical_sizes:
        print(f"AAV empirical sizes file: {args.aav_empirical_sizes}", file=screen_output)
    
    aav_fragment_gen = AAVFragmentGenerator(
        aav_seq=aav_seq,
        method_str=args.aav_method,
        verbose_output=screen_output,
        prob_file=args.aav_prob_file if args.aav_prob_file else None,
        empirical_sizes=args.aav_empirical_sizes if args.aav_empirical_sizes else None,
        seed=args.seed
    )
    
    # Generate episomal AAV fragments
    for i in range(args.events):
        fastq_entry = generate_episome_fragment(aav_fragment_gen, i + 1)
        print(fastq_entry)
    
    if not args.verbose:
        screen_output.close()


if __name__ == "__main__":
    main()
