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

import sys, argparse, random
from pulldown_probability import PulldownProbability



if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This takes a stream of SAM lines and applies a target enrichment
                       on the primary alignments to simulate to AAV integration site target enrichment.""")

    # Input options:
    parser.add_argument("--amp", default=10, type=int,
        help="Fold amplification, pre-enrichment.")

    parser.add_argument("--leak", default=1, type=float,
        help="Percent of reads that pass through without enrichment.")

    parser.add_argument("--ep_model", default="Twist", type=str,
        help="""Which empirical model to use for enrichment probability.""")

    parser.add_argument("--verbose", action='store_true')

    parser.add_argument("--seed", default=None, type=int,
        help="Random seed for reproducibility.")


    # Parse the arguments.
    args = parser.parse_args()
    ep_model = PulldownProbability(args.ep_model)
    
    # Set random seed for reproducibility
    if args.seed is not None:
        random.seed(args.seed)

    if args.verbose:
        print(args, file=sys.stderr)


    # Read each SAM formatted line streamed from stdin:
    for sam_line in sys.stdin:
        # Skip header lines:
        sam_line = sam_line.strip()
        if sam_line.startswith('@'):
            continue

        # The alignment must be the primary (bitwise flag = 0)
        # or the primary reverse complementary (bitwise flag = 16)
        # or the read must be unmapped (bitwise flag = 4):
        sam_cols = sam_line.split('\t')
        if sam_cols[1] != '0' and sam_cols[1] != '16' and sam_cols[1] != '4':
            continue

        # Convert the CIGAR string(s) to a probability of pulldown.
        # Unmapped reads have no CIGAR string, assign zero probability:
        if sam_cols[1] == '4':
            prob = 0.0
        # Else apply pulldown model:
        else:
            prob = ep_model.probability_function(sam_cols[5])

        # Output reads if they are pulled down or leaked:
        Npl = 0
        for amp_i in range(args.amp): # pre-pulldown PCR amplification
            # If the sequence got pulled down:
            if prob >= random.random():
                Npl += 1
                print(f">{sam_cols[0]}_amp{Npl}\n{sam_cols[9]}")

            # If the sequence did not get pulled down, apply leakage:
            elif args.leak >= (100 * random.random()):
                Npl += 1
                print(f">{sam_cols[0]}_leak{Npl}\n{sam_cols[9]}")



