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
from collections import defaultdict


def read_chrom_lengths(chrom_file):
    """Read chromosome lengths from a TSV file (chrom\tlength format)."""
    chrom_lengths = {}
    with open(chrom_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                chrom = parts[0]
                length = int(parts[1])
                chrom_lengths[chrom] = length
    return chrom_lengths


def read_probability_bed(bed_file):
    """
    Read a BED file with sampling probabilities.
    Format: chrom\tstart\tend\tprobability
    Returns a dict mapping chromosome to list of (start, end, prob) tuples.
    """
    prob_regions = defaultdict(list)
    with open(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 4:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                prob = float(parts[3])
                prob_regions[chrom].append((start, end, prob))
    
    # Sort regions by start position for efficient lookup
    for chrom in prob_regions:
        prob_regions[chrom].sort(key=lambda x: x[0])
    
    return prob_regions


def calculate_genome_weights(chrom_lengths, prob_regions, baseline_prob=1.0):
    """
    Calculate cumulative weights for each chromosome position.
    This allows for efficient weighted sampling.
    
    Returns:
        - total_weight: sum of all weights across genome
        - chrom_weights: dict mapping chromosome to (cumulative_weight, length)
    """
    total_weight = 0.0
    chrom_weights = {}
    
    for chrom in sorted(chrom_lengths.keys()):
        length = chrom_lengths[chrom]
        
        # Calculate weight for this chromosome
        if chrom in prob_regions:
            # Sum weights for regions with custom probabilities
            chrom_weight = 0.0
            last_end = 0
            
            for start, end, prob in prob_regions[chrom]:
                # Add baseline weight for gap before this region
                if start > last_end:
                    chrom_weight += (start - last_end) * baseline_prob
                
                # Add custom weight for this region
                chrom_weight += (end - start) * prob
                last_end = end
            
            # Add baseline weight for remaining part of chromosome
            if last_end < length:
                chrom_weight += (length - last_end) * baseline_prob
        else:
            # Entire chromosome has baseline probability
            chrom_weight = length * baseline_prob
        
        chrom_weights[chrom] = (total_weight + chrom_weight, length)
        total_weight += chrom_weight
    
    return total_weight, chrom_weights


def sample_breakpoints(chrom_lengths, n_events, fragment_len, prob_regions=None, 
                       baseline_prob=1.0, seed=None):
    """
    Sample random breakpoints from the genome and expand them to fragment_len.
    
    Args:
        chrom_lengths: dict mapping chromosome names to their lengths
        n_events: number of breakpoints to sample
        fragment_len: desired fragment length (will be centered on breakpoint if possible)
        prob_regions: dict mapping chrom to list of (start, end, prob) for weighted sampling
        baseline_prob: probability weight for regions not in prob_regions
        seed: random seed for reproducibility
    
    Yields:
        tuples of (chrom, start, end) in 0-based coordinates
    """
    if seed is not None:
        random.seed(seed)
    
    # If no probability regions provided, use uniform sampling
    if prob_regions is None:
        prob_regions = {}
    
    # Calculate chromosome weights for weighted sampling
    total_weight, chrom_weights = calculate_genome_weights(
        chrom_lengths, prob_regions, baseline_prob
    )
    
    # Sample n_events breakpoints
    for _ in range(n_events):
        # Sample a random weighted position in the genome
        target_weight = random.uniform(0, total_weight)
        
        # Find which chromosome this weight falls on
        selected_chrom = None
        chrom_start_weight = 0.0
        
        for chrom in sorted(chrom_lengths.keys()):
            cumulative_weight, length = chrom_weights[chrom]
            if target_weight <= cumulative_weight:
                selected_chrom = chrom
                break
            chrom_start_weight = cumulative_weight
        
        if selected_chrom is None:
            # Shouldn't happen, but fallback to last chromosome
            selected_chrom = sorted(chrom_lengths.keys())[-1]
        
        # Now sample a position within this chromosome using weighted sampling
        chrom_length = chrom_lengths[selected_chrom]
        chrom_target_weight = target_weight - chrom_start_weight
        
        # Find position within chromosome that matches this weight
        chrom_pos = sample_position_in_chrom(
            selected_chrom, chrom_length, chrom_target_weight,
            prob_regions, baseline_prob
        )
        
        # The breakpoint is a 2bp region
        breakpoint_start = chrom_pos
        breakpoint_end = min(chrom_pos + 2, chrom_length)
        
        # Expand to fragment_len, centered on breakpoint if possible
        half_len = fragment_len // 2
        
        # Try to center the fragment on the breakpoint
        frag_start = breakpoint_start - half_len
        frag_end = frag_start + fragment_len
        
        # Handle edge cases at chromosome boundaries
        if frag_start < 0:
            # Extend further to the right
            frag_start = 0
            frag_end = min(fragment_len, chrom_length)
        elif frag_end > chrom_length:
            # Extend further to the left
            frag_end = chrom_length
            frag_start = max(0, chrom_length - fragment_len)
        
        # Ensure we still have the desired fragment length (or as close as possible)
        # In case chromosome is shorter than fragment_len
        if frag_end - frag_start < fragment_len and frag_start == 0:
            frag_end = min(fragment_len, chrom_length)
        elif frag_end - frag_start < fragment_len and frag_end == chrom_length:
            frag_start = max(0, chrom_length - fragment_len)
        
        yield (selected_chrom, frag_start, frag_end)


def sample_position_in_chrom(chrom, chrom_length, target_weight, 
                             prob_regions, baseline_prob):
    """
    Sample a position within a chromosome given a target cumulative weight.
    
    Args:
        chrom: Chromosome name
        chrom_length: Total length of chromosome
        target_weight: Target cumulative weight to find position for
        prob_regions: Dictionary mapping chromosome to list of (start, end, prob) tuples
        baseline_prob: Baseline probability for regions not in prob_regions
    
    Returns:
        int: Sampled position (0-based) within the chromosome
    """
    cumulative_weight = 0.0
    
    if chrom in prob_regions:
        last_end = 0
        
        for start, end, prob in prob_regions[chrom]:
            # Check baseline region before this custom region
            if start > last_end:
                gap_weight = (start - last_end) * baseline_prob
                if cumulative_weight + gap_weight >= target_weight:
                    # Position is in this baseline gap
                    offset = (target_weight - cumulative_weight) / baseline_prob
                    return int(last_end + offset)
                cumulative_weight += gap_weight
            
            # Check custom probability region
            region_weight = (end - start) * prob
            if cumulative_weight + region_weight >= target_weight:
                # Position is in this custom region
                offset = (target_weight - cumulative_weight) / prob
                return int(start + offset)
            cumulative_weight += region_weight
            last_end = end
        
        # Check remaining baseline region
        if last_end < chrom_length:
            gap_weight = (chrom_length - last_end) * baseline_prob
            if cumulative_weight + gap_weight >= target_weight:
                offset = (target_weight - cumulative_weight) / baseline_prob
                return int(last_end + offset)
    else:
        # Entire chromosome is baseline
        offset = target_weight / baseline_prob
        return int(offset)
    
    # Fallback (shouldn't happen with correct weights)
    return int(chrom_length * random.random())


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Sample random genomic breakpoints and expand them to fragments."
    )
    
    parser.add_argument(
        "--chrom_lengths",
        required=True,
        type=str,
        help="File containing chromosome lengths (TSV format: chrom<tab>length)"
    )
    parser.add_argument(
        "--n_events",
        required=True,
        type=int,
        help="Number of breakpoints/fragments to sample"
    )
    parser.add_argument(
        "--fragment_len",
        required=True,
        type=int,
        help="Length of fragment to extract around each breakpoint"
    )
    parser.add_argument(
        "--prob_bed",
        type=str,
        default=None,
        help="Optional BED file with sampling probabilities (format: chrom<tab>start<tab>end<tab>prob). "
             "Regions not in this file will have baseline probability of 1.0"
    )
    parser.add_argument(
        "--baseline_prob",
        type=float,
        default=1.0,
        help="Baseline sampling probability for regions not in prob_bed file"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output BED file (default: stdout)"
    )
    
    args = parser.parse_args()
    
    # Read chromosome lengths
    chrom_lengths = read_chrom_lengths(args.chrom_lengths)
    
    if not chrom_lengths:
        sys.stderr.write("Error: No chromosome lengths found in input file\n")
        sys.exit(1)
    
    # Read probability regions if provided
    prob_regions = None
    if args.prob_bed:
        prob_regions = read_probability_bed(args.prob_bed)
        sys.stderr.write(f"Loaded probability regions for {len(prob_regions)} chromosomes\n")
    
    # Open output file or use stdout
    if args.output:
        outf = open(args.output, 'w')
    else:
        outf = sys.stdout
    
    try:
        # Sample breakpoints and write BED format
        for chrom, start, end in sample_breakpoints(
            chrom_lengths, args.n_events, args.fragment_len, 
            prob_regions, args.baseline_prob, args.seed
        ):
            outf.write(f"{chrom}\t{start}\t{end}\n")
    finally:
        if args.output:
            outf.close()


if __name__ == "__main__":
    main()
