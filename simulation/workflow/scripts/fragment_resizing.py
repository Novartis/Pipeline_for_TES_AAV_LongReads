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


import sys, os, argparse, random
import numpy as np
from badread_dependencies import FragmentLengths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This takes a stream of fastq lines from wgsim and trims them so the
                       fragment lengths become gamma distributed.""")

    # Input options:
    parser.add_argument("--mfl", default=9000, type=int,
        help="Mean fragment length.")

    parser.add_argument("--sfl", default=3000, type=int,
        help="Standard deviation of fragment length.")

    parser.add_argument("--verbose", action='store_true')

    # Parse the arguments:
    args = parser.parse_args()

    if args.verbose:
        screen_output = sys.stderr
        print(args, file=screen_output)
    else:
        screen_output = open(os.devnull, 'w')

    # Initiate the fragment length objective:
    frag_lengths = FragmentLengths(args.mfl, args.sfl, screen_output)

    # Read fastq entries from stdin until EOF:
    while True:
        header = sys.stdin.readline().strip()
        if not header:
            break
        seq = sys.stdin.readline().strip()
        plus = sys.stdin.readline().strip()
        qual = sys.stdin.readline().strip()

        fragment_length = frag_lengths.get_fragment_length()
        # If fragment is longer than input, skip:
        if len(seq) <= fragment_length:
            print(f"{header}\n{seq}\n{plus}\n{qual}")
            continue

        # Write trimmed reads:
        header_split = header.split('_')
        # Trim header junk:
        header_split = header_split[0:3]
        trim_len = len(seq) - fragment_length
        # Bernoulli sampling if trimming from left or right:
        if np.random.binomial(1, 0.5):
            # Trim from left:
            header_split[1] = str(int(header_split[1]) + trim_len)
            header = '_'.join(header_split)
            print(f"{header}\n{seq[trim_len:]}\n{plus}\n{qual[trim_len:]}")
        else:
            # Trim from right:
            header_split[2] = str(int(header_split[2]) - trim_len)
            header = '_'.join(header_split)
            print(f"{header}\n{seq[:-trim_len]}\n{plus}\n{qual[:-trim_len]}")



