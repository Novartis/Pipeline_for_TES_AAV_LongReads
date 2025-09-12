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
from Bio import SeqIO
from Bio.Seq import Seq
from badread_dependencies import FragmentLengths



# Randomly break a sequence into pieces.
# Return all the pieces when/if their length=k:
def break_seq(seq, k):
    # Base case: if the sequence is too short, return an empty list
    if len(seq) < k:
        return([])

    # Randomly select a breaking point:
    break_point = random.randint(1, len(seq) - 1)

    # Define the two pieces:
    piece1 = seq[0:break_point]
    piece2 = seq[break_point:]

    # Initialize a list to store fragments of length k:
    fragments = list()

    # Check if any piece has length k and add to fragments list,
    # else start recursive process of further breaking:
    if len(piece1) == k:
        fragments.append(piece1)
    else:
        fragments.extend(break_seq(piece1, k))

    if len(piece2) == k:
        fragments.append(piece2)
    else:
        fragments.extend(break_seq(piece2, k))

    return(fragments)


# Generate a single fragment given a requested size:
def generate_fragment(seq, size):
    if size >= len(seq):
        return(seq)

    while True:
        broken_seqs = break_seq(seq, size)
        if broken_seqs:
            return(random.choice(broken_seqs))


def fragment_integration(intgr_seq, intgr_start, intgr_end, frag_lengths):
    debug = False

    # Generate fragments until the AAV insertion has been covered:
    frag_seq_lst = list()
    frag_side_lst = list()
    flank_left = 0
    read_split = 0
    flank_right = 0
    frag_idx_start = 0
    frag_idx_end = 0
    while frag_idx_end < len(intgr_seq):
        frag_len = frag_lengths.get_fragment_length()
        frag_idx_end = frag_idx_start + frag_len
        if frag_idx_end > len(intgr_seq):
            frag_idx_end = len(intgr_seq)

        if debug:
            print(frag_idx_start, frag_idx_end)

        # Overlap between fragment and AAV insertion:
        if frag_idx_end > intgr_start and frag_idx_start < intgr_end:
            frag_seq = intgr_seq[frag_idx_start:frag_idx_end]
            frag_seq_lst.append(frag_seq)

            # Fragment is in the middle of the insertion:
            if frag_idx_start > intgr_start and frag_idx_end < intgr_end:
                if debug:
                    print("Fragment is in the middle of the insertion")
                frag_seq_lst.pop()  # Discard fragment
            # Fragment __is__ the insertion:
            elif frag_idx_start == intgr_start and frag_idx_end == intgr_end:
                if debug:
                    print("Fragment __is__ the insertion")
                return()
            # Fragment spans the whole insertion:
            elif frag_idx_start < intgr_start and frag_idx_end > intgr_end:
                if debug:
                    print("Fragment spans the whole insertion")
                flank_left = intgr_start - frag_idx_start
                read_split = 0
                flank_right = frag_idx_end - intgr_end
                frag_side_lst.append('M')
                return(frag_seq_lst, frag_side_lst, flank_left, read_split, flank_right)
            # Fragment is the left fragment with no flank:
            elif frag_idx_start == intgr_start and frag_idx_end < intgr_end:
                if debug:
                    print("Fragment is the left fragment with no flank")
                frag_seq_lst.pop()  # Discard fragment
            # Fragment is the left fragment:
            elif frag_idx_start < intgr_start and frag_idx_end <= intgr_end:
                if debug:
                    print("Fragment is the left fragment")
                flank_left = intgr_start - frag_idx_start
                read_split = frag_idx_end - intgr_start
                frag_side_lst.append('L')
                if frag_idx_end == intgr_end:
                    return(frag_seq_lst, frag_side_lst, flank_left, read_split, flank_right)
            # Fragment is the right fragment with no flank:
            elif frag_idx_start > intgr_start and frag_idx_end == intgr_end:
                if debug:
                    print("Fragment is the right fragment with no flank")
                frag_seq_lst.pop()  # Discard fragment
                return(frag_seq_lst, frag_side_lst, flank_left, read_split, flank_right)
            # Fragment is the right fragment:
            elif frag_idx_start >= intgr_start and frag_idx_end > intgr_end:
                if debug:
                    print("Fragment is the right fragment")
                read_split_tmp = frag_idx_start - intgr_start
                # Left fragment with no flanks:
                if not read_split:
                    read_split = read_split_tmp
                # If the insertion has a middle fragment:
                elif read_split != read_split_tmp:
                    read_split = (read_split, read_split_tmp)
                flank_right = frag_idx_end - intgr_end
                frag_side_lst.append('R')
                assert(len(frag_seq_lst) <= 2)
                return(frag_seq_lst, frag_side_lst, flank_left, read_split, flank_right)
            else:
                # Shouldn't get here:
                raise Exception("Adding integration failed. Somehow the code is broken.")

        frag_idx_start = frag_idx_end





if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This takes a stream of fastq lines from wgsim and trims them so the
                       fragment lengths become gamma distributed.""")

    # Input/output options:
    parser.add_argument("--mfl", default=9000, type=int,
        help="Mean fragment length.")

    parser.add_argument("--sfl", default=3000, type=int,
        help="Standard deviation of fragment length.")

    parser.add_argument("--aav", required=True, type=str,
        help="Fasta file of the AAV genome.")

    parser.add_argument("--mal", required=True, type=int,
        help="Mean integrated AAV length.")

    parser.add_argument("--sal", required=True, type=int,
        help="Standard deviation of integrated AAV length.")

    parser.add_argument("--min_aav", default=25, type=int,
        help="Standard deviation of integrated AAV length.")

    parser.add_argument("--exp", default=1, type=int,
        help="Fold expansion of each event.")

    parser.add_argument("--int_out", required=True, type=str,
        help="Output file for information about the AAV integrations.")

    parser.add_argument("--fq_out", action='store_true',
        help="Print output in fastq format.")

    parser.add_argument("--verbose", action='store_true')

    debug = False

    # Parse the arguments:
    args = parser.parse_args()

    devnull = open(os.devnull, 'w')
    if args.verbose:
        screen_output = sys.stderr
        print(args, file=screen_output)
    else:
        screen_output = devnull


    # Initiate the integrated AAV length objective:
    print("AAV integration length distribution", file=screen_output)
    aav_lengths = FragmentLengths(args.mal, args.sal, screen_output)

    # Initiate the fragment length objective (only for verbose):
    print("Read length distribution", file=screen_output)
    frag_lengths = FragmentLengths(args.mfl, args.sfl, screen_output)

    # Read in the AAV sequence:
    aav_seq = [str(record.seq).upper() for record in SeqIO.parse(args.aav, "fasta")]
    # Warn if multiple AAV sequences are found:
    if len(aav_seq) > 1:
        print(f"Warning: Multiple sequences found for the AAV: {args.aav}\nOnly the first sequence will be used.")
        # raise Exception(f"Multiple sequences found for the AAV: {args.aav}\nOnly one sequence allowed.")
    aav_seq = aav_seq[0]


    # Write integration information to tab separated file:
    with open(args.int_out, 'w') as fh_out:
        int_out_header = ["chrom", "start", "event", "expansion",
                          "aav_start", "aav_end", "flank_left",
                          "read_split", "flank_right", "frag_sides", "frag_seqs"]
        # Here is how left_flank, read_split and right_flank works.
        # Left marked by "<", AAV insertion by "-" and right by ">":
        # <<<<<<<<<<<------>>>>>>>>>>>  flank_left=11, read_split=0, flank_right=11
        # <<<<<<---   ---->>>>>>>>>>    flank_left=6, read_split=3, flank_right=11
        # ------->>>>>>>>>>    flank_left=0, read_split=0, flank_right=11
        # <<<<<<---   --  -->>>>>>>>>>    flank_left=6, read_split=(3,5), flank_right=11
        int_out_header_str = '\t'.join(int_out_header)
        print(int_out_header_str, file=fh_out)


        # Read the from stdin until EOF:
        event = 0
        while True:
            event += 1
            header = sys.stdin.readline().strip()
            if not header:
                break
            seq = sys.stdin.readline().strip()

            # Split header:
            header_split = header.split('_')
            chrom, start, end = header_split[0:3]
            chrom = chrom[1:]
            start = int(start)
            end = int(end)


            # Sample AAV integration:
            while True:
                # Enforce minimum AAV length:
                aav_len = aav_lengths.get_fragment_length()
                if args.min_aav <= aav_len:
                    break
            aav_frag = generate_fragment(aav_seq, aav_len)
            aav_start = aav_seq.index(aav_frag) + 1
            aav_end = aav_start + len(aav_frag) - 1

            if debug:
                print("###### AAV fragment #######")
                print(aav_frag)
                print("###########################")

            # Randomly offset the start of the insertion:
            event_offset = random.randint(0, args.mfl-1)
            seq_left_len = len(seq) // 2 - args.mfl // 2 + event_offset
            event_start = start + seq_left_len
            intgr_seq = seq[:seq_left_len] + aav_frag + seq[seq_left_len:]
            intgr_start = seq_left_len
            intgr_end = intgr_start + len(aav_frag)

            # For each expansion, fragment the AAV integrated sequence:
            for exp_i in range(args.exp):
                frag_seq_lst, frag_side_lst, flank_left, read_split, flank_right = fragment_integration(intgr_seq, intgr_start, intgr_end, frag_lengths)

                # Print fragments:
                for frag_seq, lr in zip(frag_seq_lst, frag_side_lst):
                    header = f"{chrom}_{event_start}_AAVint_{event}_{exp_i+1}_{lr}"
                    if args.fq_out:
                        frag_qual = '~'*len(frag_seq)
                        print(f"@{header}\n{frag_seq}\n+\n{frag_qual}")
                    else:
                        print(f">{header}\n{frag_seq}")

                # Write integration table:
                int_out_lst = [chrom, event_start, event, exp_i+1, aav_start, aav_end,
                               flank_left, read_split, flank_right, frag_side_lst, frag_seq_lst]
                int_out_str = '\t'.join(map(str, int_out_lst))
                print(int_out_str, file=fh_out)



