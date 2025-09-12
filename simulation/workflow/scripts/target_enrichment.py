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

import sys, os, argparse, re, random, math


def enrichment_probability_all(q, m, sub_str_len):
    p_exl = math.log(1)
    for x in sub_str_len:
        one_ep = enrichment_probability(q, m, x)
        if one_ep == 1:
            return(1)
        p_exl += math.log(1 - one_ep)

    p = 1 - math.exp(p_exl)
    return(p)


def enrichment_probability(q, m, x):
    p = 1 - 1 / (math.exp(q*(x-m)))
    if p > 0:
        return(p)
    else:
        return(0)


def cs_lengths(sam_cols):
    # Find CIGAR string of primary and all supplementary alignments:
    cigar_lst = [sam_cols[5]]
    if len(sam_cols) > 11:
        for tag in sam_cols[11:]:
            if tag.startswith("SA:Z:"):
                cigar_lst.extend([sa.split(',')[3] for sa in tag[5:].split(';') if sa])
                # Break after seeing the SA tag, nothing more to do:
                break

    # Convert CIGAR strings to match lengths:
    lengths = list()
    for cigar in cigar_lst:
        lengths.extend(cigar_match_lengths(cigar))
    return(lengths)


def cigar_match_lengths(cigar):
    # Find the match length(s) from the alignment CIGAR:
    match_len = list()
    if cigar == '*':
        return(match_len)

    cigar_split = re.split(r'([MIDNSH])', cigar)[:-1]
    while len(cigar_split):
        op = cigar_split.pop()
        num = int(cigar_split.pop())
        if op == 'M':
            match_len.append(num)
    return(match_len)




### This does not quite work b/c alignments are joined
# even when the alignment score is low.
# Modify by testign AS:i: >= min alignment length
# ms:i: likely contains the max score of a continuous alignment

def cigar_lcs_length(cigar):
    # Find the length of the longest common substring
    # from the alignment CIGAR:
    if cigar == '*':
        return(0)

    cigar_split = re.split(r'([MIDNSH])', cigar)[:-1]
    max_match = 0
    while len(cigar_split):
        op = cigar_split.pop()
        num = int(cigar_split.pop())
        if op == 'M' and num > max_match:
            max_match = num
    return(max_match)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This takes a stream of SAM lines and applies a target enrichment
                       to simulate to AAV integration site target enrichment.""")

    # Input options:
    parser.add_argument("--amp", default=10, type=int,
        help="Fold amplification, pre-enrichment.")

    parser.add_argument("--leak", default=1, type=int,
        help="Percent of reads that pass through without enrichment.")

    parser.add_argument("--ep", default="0.05,15", type=str,
        help="""Enrichment probability string . This defined the probability of pulling
                down a sequence, which depeneds on a minimum length of complementarity (m)
                and a scaling parameter (q). The probability is defined by this function:
                f(x) = 1 - 1 / (exp(q*(x-m))) if x>m else 0
             """)

    parser.add_argument("--add", action='store_true',
        help="""Let the enrichment probability be additive with respect to multiple
                binding sequences, found as supplementary alignments.""")

    parser.add_argument("--verbose", action='store_true')


    # Parse the arguments.
    # Convert the ep string to numbers:
    args = parser.parse_args()
    args.ep_q, args.ep_m = args.ep.split(',')
    args.ep_q = float(args.ep_q)
    args.ep_m = int(args.ep_m)

    if args.verbose:
        print(args, file=sys.stderr)


    # Read each SAM formatted line streamed from stdin:
    for sam_line in sys.stdin:
        # Skip header lines:
        sam_line = sam_line.strip()
        if sam_line.startswith('@'):
            continue

        # The alignment must be the primary (bitwise flag = 0)
        # or the primary reverse complementary (bitwise flag = 16):
        sam_cols = sam_line.split('\t')
        if sam_cols[1] != '0' and sam_cols[1] != '16':
            continue

        # Find lengths of common substrings:
        if args.add:
            raise Exception("Not implemented")
            sub_str_len = cs_lengths(sam_cols)
        else:
            # Only use the longest common substring:
            lcs_length = cigar_lcs_length(sam_cols[5])
            sub_str_len = [lcs_length]


        # The lengths of common substrings is proportional
        # to the binding potential. Convert this to a
        # probability and sample if passing enrichment or not:
        prob = enrichment_probability_all(args.ep_q, args.ep_m, sub_str_len)
        Npass = 0
        if prob > 0:
            for amp_i in range(args.amp):
                if prob >= random.random():
                    Npass += 1
                    print(f">{sam_cols[0]}_enrich{Npass}\n{sam_cols[9]}")

        # If the sequence did not pass, apply leakage:
        if args.leak >= (100 * random.random()):
            Npass += 1
            print(f">{sam_cols[0]}_leak\n{sam_cols[9]}")



