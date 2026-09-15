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
"""
Convert *_integration_sites.tsv -> BED6 for the enrichment analysis.

Produces a standard BED6 file (chrom, start, end, name, score, strand) where:
  - each site is a 1-bp point interval at the midpoint of the two breakpoints
    (or at the single resolved breakpoint when only one is available)
  - positions are integers
  - score = supporting read count
  - strand = orientation column (fwd -> +, rev -> -, mixed/other -> .)
  - translocations, rows with missing positions, non-standard chromosomes,
    and sites with unreasonably large breakpoint distances are excluded

Usage:
    python preprocess_integration_site.py <integration_sites.tsv> [output.bed]

If output.bed is omitted, writes to ../input/aav_sites_filtered.bed
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Convert integration_sites.tsv to a valid BED6 file."
    )
    parser.add_argument(
        "input_tsv",
        help="Path to *_integration_sites.tsv",
    )
    parser.add_argument(
        "output_bed",
        nargs="?",
        default=None,
        help="Output BED path (default: ../input/aav_sites_filtered.bed)",
    )
    parser.add_argument(
        "--min-reads",
        type=int,
        default=1,
        help="Minimum supporting read count to include a site (default: 1)",
    )
    parser.add_argument(
        "--max-distance",
        type=int,
        default=10000,
        help="Maximum allowed distance (bp) between bp1 and bp2. Sites exceeding "
             "this are likely merge artifacts and are excluded (default: 10000)",
    )
    parser.add_argument(
        "--keep-non-standard-chroms",
        action="store_true",
        default=False,
        help="Keep sites on non-standard chromosomes (chrM, random, alt, Un). "
             "By default only chr1-22, chrX, chrY are retained.",
    )
    parser.add_argument(
        "--exclude-confidence",
        nargs="*",
        default=None,
        help="Exclude sites with these confidence values (e.g., --exclude-confidence low). "
             "Requires a 'confidence' column in the input TSV.",
    )
    args = parser.parse_args()

    if args.output_bed is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        args.output_bed = os.path.join(script_dir, "..", "input", "aav_sites_filtered.bed")

    # -- Load --
    df = pd.read_csv(args.input_tsv, sep="\t")
    print(f"Loaded {len(df)} sites from {args.input_tsv}")

    # -- Filter --
    # Exclude translocations (chr_bp_upstream != chr_bp_downstream)
    df = df[df["is_translocation"] == False].copy()
    print(f"  After removing translocations: {len(df)}")

    # Drop rows where NEITHER breakpoint is fully resolved
    has_bp1 = df["chr_bp_upstream"].notna() & df["pos_bp_upstream"].notna()
    has_bp2 = df["chr_bp_downstream"].notna() & df["pos_bp_downstream"].notna()
    df = df[has_bp1 | has_bp2].copy()
    print(f"  After removing rows with no resolved breakpoint: {len(df)}")

    # Confidence filter
    if args.exclude_confidence and "confidence" in df.columns:
        exclude_set = set(args.exclude_confidence)
        n_before = len(df)
        df = df[~df["confidence"].isin(exclude_set)].copy()
        print(f"  After excluding confidence={list(exclude_set)}: {len(df)} "
              f"(dropped {n_before - len(df)})")

    # Minimum read count filter
    df = df[df["count"] >= args.min_reads].copy()
    print(f"  After min-reads filter (>={args.min_reads}): {len(df)}")

    # -- Build BED6 --
    has_bp1 = df["chr_bp_upstream"].notna() & df["pos_bp_upstream"].notna()
    has_bp2 = df["chr_bp_downstream"].notna() & df["pos_bp_downstream"].notna()

    # Chromosome: prefer bp1; fall back to bp2
    chrom = df["chr_bp_upstream"].where(has_bp1, df["chr_bp_downstream"])

    # Filter to standard chromosomes (chr1-22, chrX, chrY)
    if not args.keep_non_standard_chroms:
        standard = chrom.str.match(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y)$")
        n_before = len(df)
        df = df[standard.values].copy()
        has_bp1 = has_bp1[standard.values]
        has_bp2 = has_bp2[standard.values]
        chrom = chrom[standard.values]
        print(f"  After removing non-standard chromosomes: {len(df)} "
              f"(dropped {n_before - len(df)})")

    # Positions - compute midpoint for point intervals
    p1 = pd.to_numeric(df["pos_bp_upstream"], errors="coerce")
    p2 = pd.to_numeric(df["pos_bp_downstream"], errors="coerce")

    # For sites with both breakpoints, filter by max distance
    both = has_bp1 & has_bp2
    bp_distance = np.abs(p1 - p2)
    too_far = both & (bp_distance > args.max_distance)
    if too_far.any():
        n_dropped = too_far.sum()
        keep = ~too_far
        df = df[keep.values].copy()
        has_bp1 = has_bp1[keep.values]
        has_bp2 = has_bp2[keep.values]
        chrom = chrom[keep.values]
        p1 = p1[keep.values]
        p2 = p2[keep.values]
        both = has_bp1 & has_bp2
        print(f"  After max-distance filter (<={args.max_distance} bp): {len(df)} "
              f"(dropped {n_dropped})")

    # Point interval: midpoint when two breakpoints, single breakpoint otherwise
    p1_filled = p1.fillna(0)
    p2_filled = p2.fillna(0)
    midpoint = np.where(
        both,
        ((p1_filled + p2_filled) / 2).astype(int),
        np.where(has_bp1, p1_filled, p2_filled),
    ).astype(int)

    bed = pd.DataFrame()
    bed["chrom"] = chrom.values
    bed["start"] = midpoint
    bed["end"]   = midpoint + 1

    bed["name"] = "AAV_site_" + df["site_id"].astype(str).values
    bed["score"] = df["count"].astype(int).values

    # Map orientation to BED strand; use "." for missing/NaN/mixed values
    strand_map = {"fwd": "+", "rev": "-"}
    bed["strand"] = df["orientation"].map(strand_map).where(
        df["orientation"].isin(strand_map), "."
    ).values

    # Sort by chrom, start
    bed = bed.sort_values(["chrom", "start"]).reset_index(drop=True)

    # -- Write --
    os.makedirs(os.path.dirname(os.path.abspath(args.output_bed)), exist_ok=True)
    bed.to_csv(args.output_bed, sep="\t", header=False, index=False)

    print(f"\nWrote {len(bed)} sites to {args.output_bed}")


if __name__ == "__main__":
    main()
