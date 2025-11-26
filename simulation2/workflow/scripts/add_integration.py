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
from indel_generator import DeletionGenerator, InsertionGenerator
from aav_fragment_generator import AAVFragmentGenerator



def fragment_integration(intgr_seq, intgr_start, intgr_end, frag_lengths):
    """
    Fragment a sequence containing an AAV insertion and identify which fragments
    overlap the insertion boundaries.
    
    Returns fragments that span the insertion boundaries (have both genomic and AAV content).
    
    Args:
        intgr_seq: Sequence with AAV inserted (left_flank + AAV + right_flank)
        intgr_start: Start position of AAV insertion in intgr_seq
        intgr_end: End position of AAV insertion in intgr_seq
        frag_lengths: Fragment length generator
    
    Returns:
        tuple: (frag_seq_lst, frag_side_lst, flank_left, frag_split, flank_right)
        - frag_seq_lst: List of fragment sequences that cross insertion boundaries
        - frag_side_lst: List of 'L', 'R', or 'M' indicating which side(s) each fragment spans
        - flank_left: Length of genomic sequence on left side of insertion
        - frag_split: Position(s) where AAV is split (0 if fragment spans entire insertion)
        - flank_right: Length of genomic sequence on right side of insertion
    """
    debug = False
    
    # Output collections
    frag_seq_lst = []
    frag_side_lst = []
    flank_left = 0
    frag_split = 0
    flank_right = 0
    
    # Fragment the entire sequence
    frag_idx_start = 0
    while frag_idx_start < len(intgr_seq):
        # Get next fragment length and calculate end position
        frag_len = frag_lengths.get_fragment_length()
        frag_idx_end = min(frag_idx_start + frag_len, len(intgr_seq))
        
        if debug:
            print(f"Fragment: {frag_idx_start}-{frag_idx_end}, Insertion: {intgr_start}-{intgr_end}")
        
        # Check if this fragment overlaps the insertion
        overlaps_insertion = (frag_idx_end > intgr_start and frag_idx_start < intgr_end)
        
        if overlaps_insertion:
            # Calculate overlap characteristics
            has_left_flank = frag_idx_start < intgr_start
            has_right_flank = frag_idx_end > intgr_end
            entirely_within_insertion = frag_idx_start >= intgr_start and frag_idx_end <= intgr_end
            spans_entire_insertion = frag_idx_start < intgr_start and frag_idx_end > intgr_end
            
            # Case 1: Fragment spans entire insertion (has both flanks)
            if spans_entire_insertion:
                if debug:
                    print("  → Fragment spans the whole insertion")
                frag_seq = intgr_seq[frag_idx_start:frag_idx_end]
                frag_seq_lst.append(frag_seq)
                frag_side_lst.append('M')  # Middle - has both sides
                flank_left = intgr_start - frag_idx_start
                frag_split = 0  # No split - entire insertion in one fragment
                flank_right = frag_idx_end - intgr_end
                return (frag_seq_lst, frag_side_lst, flank_left, frag_split, flank_right)
            
            # Case 2: Fragment is entirely within insertion (no genomic flanks)
            elif entirely_within_insertion:
                if debug:
                    print("  → Fragment entirely within insertion - discarding")
                # Don't keep fragments that are purely AAV with no genomic sequence
                pass
            
            # Case 3: Fragment has left flank (crosses left boundary)
            elif has_left_flank and not has_right_flank:
                if debug:
                    print("  → Fragment has left flank")
                frag_seq = intgr_seq[frag_idx_start:frag_idx_end]
                frag_seq_lst.append(frag_seq)
                frag_side_lst.append('L')
                flank_left = intgr_start - frag_idx_start
                frag_split = frag_idx_end - intgr_start  # How much AAV in this fragment
                
                # If this fragment ends exactly at insertion end, we're done
                if frag_idx_end == intgr_end:
                    return (frag_seq_lst, frag_side_lst, flank_left, frag_split, flank_right)
            
            # Case 4: Fragment has right flank (crosses right boundary)
            elif has_right_flank and not has_left_flank:
                if debug:
                    print("  → Fragment has right flank")
                frag_seq = intgr_seq[frag_idx_start:frag_idx_end]
                frag_seq_lst.append(frag_seq)
                frag_side_lst.append('R')
                
                # Calculate where this right fragment starts within the insertion
                aav_portion_start = frag_idx_start - intgr_start
                
                # Handle frag_split for cases with or without a left fragment
                if frag_split == 0:
                    # No left fragment - this is the only AAV-spanning fragment
                    frag_split = aav_portion_start
                elif frag_split != aav_portion_start:
                    # There was a left fragment AND middle fragments
                    frag_split = (frag_split, aav_portion_start)
                # else: frag_split already equals aav_portion_start (no middle fragments)
                
                flank_right = frag_idx_end - intgr_end
                
                # Validate we have at most left and right fragments
                if len(frag_seq_lst) > 2:
                    raise RuntimeError(f"Expected at most 2 fragments (left and right), got {len(frag_seq_lst)}")
                
                return (frag_seq_lst, frag_side_lst, flank_left, frag_split, flank_right)
        
        # Move to next fragment
        frag_idx_start = frag_idx_end
    
    # If we get here, something unexpected happened
    # This could happen if insertion is exactly the size of the sequence
    return (frag_seq_lst, frag_side_lst, flank_left, frag_split, flank_right)



def parse_header(header):
    """
    Parse FASTA header to extract chromosome and coordinates.
    
    Handles two formats:
    - bedtools getfasta format: >chrom:start-end
    - wgsim format: >chrom_start_end_...
    
    Args:
        header: FASTA header string starting with '>'
    
    Returns:
        tuple: (chromosome, start_position, end_position)
    
    Raises:
        ValueError: If header format is invalid
    """
    if not header or not header.startswith('>'):
        raise ValueError(f"Invalid FASTA header format: {header}")
    
    # Remove '>' prefix
    header_no_prefix = header[1:]
    
    # bedtools getfasta format: chrom:start-end
    if ':' in header_no_prefix and '-' in header_no_prefix:
        try:
            chrom, coords = header_no_prefix.split(':', 1)
            start_str, end_str = coords.split('-', 1)
            return chrom, int(start_str), int(end_str)
        except ValueError as e:
            raise ValueError(f"Invalid bedtools header format '{header}': {e}")
    
    # wgsim format: chrom_start_end_...
    else:
        header_parts = header_no_prefix.split('_')
        if len(header_parts) < 3:
            raise ValueError(f"Invalid wgsim header format '{header}': expected at least 3 underscore-separated fields")
        try:
            chrom, start_str, end_str = header_parts[0:3]
            return chrom, int(start_str), int(end_str)
        except ValueError as e:
            raise ValueError(f"Invalid wgsim header coordinates in '{header}': {e}")





if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This takes a stream of fasta lines from bedtools or wgsim adds the AAV integration
                       and trims the sequences so the fragment lengths become gamma distributed.""")

    # Input/output options:
    parser.add_argument("--events", required=True, type=int,
        help="Number of integration events.")

    parser.add_argument("--mfl", default=9000, type=int,
        help="Mean fragment length.")

    parser.add_argument("--sfl", default=3000, type=int,
        help="Standard deviation of fragment length.")

    parser.add_argument("--aav", required=True, type=str,
        help="Fasta file of the AAV genome.")

    parser.add_argument("--aav_method", required=True, type=str,
        help="AAV fragment generation method with parameters (e.g., 'gamma_stick,1000,500,25').")

    parser.add_argument("--aav_prob_file", default="", type=str,
        help="Optional BED file with breakpoint probabilities for gamma_biased_stick method.")

    parser.add_argument("--aav_empirical_sizes", default="", type=str,
        help="Optional file with empirical fragment sizes for empirical_kde/empirical_kde_biased methods (one integer per line).")

    parser.add_argument("--aav_kde_plot", default="", type=str,
        help="Optional output file for KDE plot when using empirical_kde/empirical_kde_biased methods (e.g., 'kde_plot.png').")

    parser.add_argument("--del_max", default=0, type=int,
        help="Max deletion length.")

    parser.add_argument("--del_flag", default="", type=str,
        help="Deletion distribution flag.")

    parser.add_argument("--ins_flag", default="", type=str,
        help="Insertion distribution flag.")

    parser.add_argument("--p_trans", default=0.0, type=float,
        help="Probability of translocation.")

    parser.add_argument("--exp", nargs='+', type=int, default=[1],
        help="Fold expansion of each event. Can be a single value (applied to all events), "
             "or a list of 2 values where the first applies to event 1 and the second to all remaining events.")

    parser.add_argument("--int_out", required=True, type=str,
        help="Output file for information about the AAV integrations.")

    parser.add_argument("--fq_out", action='store_true',
        help="Print output in fastq format.")

    parser.add_argument("--verbose", action='store_true')

    parser.add_argument("--seed", default=None, type=int,
        help="Random seed for reproducibility.")

    debug = False

    # Parse the arguments:
    args = parser.parse_args()
    
    # Validate arguments
    if not (0.0 <= args.p_trans <= 1.0):
        raise ValueError(f"Translocation probability must be between 0 and 1, got {args.p_trans}")
    
    # Validate and process expansion argument
    if len(args.exp) > 2:
        raise ValueError(f"--exp can have at most 2 values, got {len(args.exp)}: {args.exp}")
    if len(args.exp) == 0:
        raise ValueError("--exp must have at least 1 value")
    
    # Convert expansion list to a function that returns expansion for each event
    if len(args.exp) == 1:
        # Single value: use for all events
        get_expansion = lambda event_num: args.exp[0]
    else:
        # Two values: first for event 1, second for all others
        get_expansion = lambda event_num: args.exp[0] if event_num == 1 else args.exp[1]
    
    # Set random seeds for reproducibility
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    devnull = open(os.devnull, 'w')
    if args.verbose:
        screen_output = sys.stderr
        print(args, file=screen_output)
    else:
        screen_output = devnull


    # Initiate the fragment length objective (for verbose output):
    print("Read length distribution", file=screen_output)
    frag_lengths = FragmentLengths(args.mfl, args.sfl, screen_output)

    # Read in the AAV sequence:
    aav_seq = [str(record.seq).upper() for record in SeqIO.parse(args.aav, "fasta")]
    # Warn if multiple AAV sequences are found:
    if len(aav_seq) > 1:
        print(f"Warning: Multiple sequences found for the AAV: {args.aav}\nOnly the first sequence will be used.")
        # raise Exception(f"Multiple sequences found for the AAV: {args.aav}\nOnly one sequence allowed.")
    aav_seq = aav_seq[0]

    # Initialize the AAV fragment generator:
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
    
    # Generate KDE plot if requested and using empirical_kde or empirical_kde_biased method
    if args.aav_kde_plot and args.aav_empirical_sizes:
        try:
            aav_fragment_gen.plot_empirical_kde(
                output_file=args.aav_kde_plot,
                title="AAV Fragment Size Distribution (Empirical KDE)"
            )
            print(f"KDE plot saved to: {args.aav_kde_plot}", file=screen_output)
        except Exception as e:
            print(f"Warning: Could not generate KDE plot: {e}", file=screen_output)

    # Initiate deletion and insertion generators:
    gen_deletion = DeletionGenerator(max_size=int(args.del_max), size_dist_str=args.del_flag)
    gen_insertion = InsertionGenerator(size_dist_str=args.ins_flag)

    # Write integration information to tab separated file:
    with open(args.int_out, 'w') as fh_out:
        int_out_header = ["start", "end", "event", "expansion",
                          "aav_start", "aav_end", "aav_dir", "flank_left", "ins_left",
                          "frag_split", "ins_right", "flank_right", "frag_sides", "frag_seqs"]
        # Here is how left_flank, frag_split and right_flank works.
        # Left marked by "<", AAV insertion by "-" and right by ">":
        # <<<<<<<<<<<------>>>>>>>>>>>  flank_left=11, frag_split=0, flank_right=11
        # <<<<<<---   ---->>>>>>>>>>    flank_left=6, frag_split=3, flank_right=11
        # ------->>>>>>>>>>    flank_left=0, frag_split=0, flank_right=11
        # <<<<<<---   --  -->>>>>>>>>>    flank_left=6, frag_split=(3,5), flank_right=11
        int_out_header_str = '\t'.join(int_out_header)
        print(int_out_header_str, file=fh_out)


        # Read from stdin until EOF:
        event = 0
        while event < args.events:
            event += 1

            # Fetch genome fragment to integrate into:
            header = sys.stdin.readline().strip()
            seq = sys.stdin.readline().strip()
            chrom_left, start_left, end_left = parse_header(header)
            dir_left = "fwd"  # Default direction for left fragment


            # Sample AAV fragment using the configured generation method:
            aav_frag, aav_start, aav_end, aav_len = aav_fragment_gen.generate()
            
            # Decide if AAV fragment is inverted or not:
            aav_dir = "fwd"  # Default direction for AAV fragment
            if random.choice([True, False]):
                # Invert AAV fragment:
                aav_dir = "rev"
                aav_frag = str(Seq(aav_frag).reverse_complement())

            # Add insertions at both junctions:
            ins_left, ins_right = gen_insertion.draw()
            aav_frag = ins_left + aav_frag + ins_right

            if debug:
                print("###### AAV fragment #######")
                print(aav_frag)
                print("###########################")


            # Insert AAV into genome fragment, either as simple insertion or translocation:
            if random.random() < args.p_trans:
                # Translocation event
                # Fetch second genome fragment to integrate into:
                header2 = sys.stdin.readline().strip()
                seq2 = sys.stdin.readline().strip()
                chrom_right, start_right, end_right = parse_header(header2)
                dir_right = "fwd"

                # Decide if second fragment is inverted:
                invert_second = random.choice([True, False])
                if invert_second:
                    dir_right = "rev"
                    start_right, end_right = end_right, start_right
                    seq2 = str(Seq(seq2).reverse_complement())
                
                # Create integrated sequence from two fragments (i.e. translocation).
                # Offset the start of the first fragment to avoid biased integration sites:
                event_offset = random.randint(0, args.mfl-1)
                intgr_seq = seq[event_offset:] + aav_frag + seq2
                intgr_start = len(seq) - event_offset
                intgr_end = intgr_start + len(aav_frag)

            else:
                # No translocation event:
                chrom_right = chrom_left
                end_right = end_left
                dir_right = dir_left

                # Insert AAV into the middle of the sequence,
                # with a randomly offset start of the insertion
                # to avoid biased integration sites:
                event_offset = random.randint(0, args.mfl-1)
                seq_left_len = len(seq) // 2 - args.mfl // 2 + event_offset
                end_left = start_left + seq_left_len
                deletion_length = gen_deletion.draw()
                start_right = end_left + deletion_length + 1
                intgr_seq = seq[:seq_left_len] + aav_frag + seq[(seq_left_len + deletion_length):]
                intgr_start = seq_left_len
                intgr_end = intgr_start + len(aav_frag)


            # For each expansion, fragment the AAV integrated sequence:
            expansion_count = get_expansion(event)
            for exp_i in range(expansion_count):
                frag_seq_lst, frag_side_lst, flank_left, frag_split, flank_right = fragment_integration(intgr_seq, intgr_start, intgr_end, frag_lengths)
                
                # Genome start position and direction for each side of the integration:
                left_pos_str = f"{chrom_left}-{end_left}-{dir_left}"
                right_pos_str = f"{chrom_right}-{start_right}-{dir_right}"

                # Print fragments:
                for frag_seq, lr in zip(frag_seq_lst, frag_side_lst):
                    header = f"{left_pos_str}_{right_pos_str}_AAVint_{event}_{exp_i+1}_{lr}"
                    if args.fq_out:
                        frag_qual = '~'*len(frag_seq)
                        print(f"@{header}\n{frag_seq}\n+\n{frag_qual}")
                    else:
                        print(f">{header}\n{frag_seq}")

                # Write integration table:
                int_out_lst = [left_pos_str, right_pos_str, event, exp_i+1, aav_start,
                               aav_end, aav_dir, flank_left, ins_left, frag_split, ins_right,
                               flank_right, frag_side_lst, frag_seq_lst]
                int_out_str = '\t'.join(map(str, int_out_lst))
                print(int_out_str, file=fh_out)


        # Read the remaining lines from stdin to EOF:
        while sys.stdin.readline():
            pass




