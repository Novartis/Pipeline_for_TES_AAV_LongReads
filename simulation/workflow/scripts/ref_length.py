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

import sys, os, argparse


def calculate_chromosome_lengths(fasta_file):
    chromosome_lengths = {}
    current_chromosome = None
    current_length = 0

    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_chromosome is not None:
                    chromosome_lengths[current_chromosome] = current_length
                current_chromosome = line[1:]
                current_length = 0
            else:
                current_length += len(line)
        
        if current_chromosome is not None:
            chromosome_lengths[current_chromosome] = current_length
    return(chromosome_lengths)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Print the length of each entry in a fasta file.""")

    # Input options:
    parser.add_argument("--fa", required=True, type=str,
        help="Input fasta file.")
    parser.add_argument("--min_len", default=0, type=int,
        help="Only print if over minimum length.")
    parser.add_argument("--max_len", default=int(1e9), type=int,
        help="Only print if less than maximum length.")
    parser.add_argument("--total", action='store_true', help="Print total length.")
    args = parser.parse_args()

    # Print:
    chromosome_lengths = calculate_chromosome_lengths(args.fa)
    total_length = 0
    for chromosome, length in chromosome_lengths.items():
        if args.min_len < length and args.max_len > length:
            print(f'{chromosome}: {length} bp')
        total_length += length
    
    if args.total:
        print(f"Total length: {total_length}")




