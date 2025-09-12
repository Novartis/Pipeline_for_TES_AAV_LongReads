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

import sys, os, argparse, re


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
        description="""Trim fasta file to remove short sequences or entries
                       with specific naming.""")

    # Input options:
    parser.add_argument("--fa", required=True, type=str,
        help="Input fasta file.")
    parser.add_argument("", default=0, type=int,
        help="Only print entry if over minimum length.")
    parser.add_argument("--max_len", default=int(1e9), type=int,
        help="Only print entry if less than maximum length.")
    parser.add_argument("--header_regex", default='', type=str,
        help="Only print if header is found by regex.")
    args = parser.parse_args()


    # Filter and print:
    chromosome_lengths = calculate_chromosome_lengths(args.fa)
    header_pass = False
    with open(args.fa, 'r') as file:
        for line in file:
            if line.startswith('>'):
                chr_len = chromosome_lengths[line[1:].strip()]
                regex_found = bool(re.search(args.header_regex, line[1:].strip()))
                if regex_found and chr_len > args.min_len and chr_len < args.max_len:
                    header_pass = True
                else:
                    header_pass = False

            if header_pass:
                print(line, end='')


