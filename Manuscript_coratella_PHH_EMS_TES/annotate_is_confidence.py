#!/usr/bin/env python3
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
annotate_is_confidence.py

Post-processing script to annotate integration sites with confidence flags.
Run this AFTER the main pipeline (resolve_integration_sites.py) to flag
potential artifacts caused by vector-genome sequence homology.

This is not part of the Snakemake pipeline. It is run as a standalone
post-processing step because the repeat regions are vector-specific and
must be defined by the user.

Requirements:
    - Python 3
    - pandas
    - pysam (pip install pysam)
    - mappy (pip install mappy) — Python bindings for minimap2

Annotates integration sites with confidence flags based on two checks:

1. Direct homology check (always run):
   Extracts genomic sequence at each breakpoint and aligns it back to the
   plasmid. If the breakpoint region has significant alignment to the
   vector, the site is likely an artifact caused by sequence homology
   (e.g., hGH polyA → GH cluster, HBB intron → HBB locus, Alu in
   vector → genomic Alu elements).

2. Repeat region check (optional):
   Checks if the AAV breakpoint falls within a user-specified repeat
   region on the plasmid (e.g., an Alu element). If yes, the artifact
   may be caused by PCR template switching at a repetitive sequence.

Output columns added:
    homology_flag       - "vector_homology" if genomic breakpoint aligns to plasmid
    repeat_flag         - "repeat_mediated" if AAV breakpoint is in a plasmid repeat region
    confidence          - "high", "caution", or "low":
                            "low" = vector_homology (high-identity, e.g. HBB/GH cluster)
                            "caution" = repeat_mediated + divergent Alu homology (Tier 2)
                            "high" = everything else (including repeat_mediated without homology)
    homology_detail     - alignment details (plasmid position, length, identity)

Usage:
    # Basic usage (homology check only):
    python annotate_is_confidence.py \
        --sites integration_sites.tsv \
        --plasmid plasmid.fa \
        --ref-genome ref.fa \
        --output integration_sites_annotated.tsv

    # With repeat region annotation (e.g., Alu element at positions 2900-3139):
    python annotate_is_confidence.py \
        --sites integration_sites.tsv \
        --plasmid plasmid.fa \
        --ref-genome ref.fa \
        --output integration_sites_annotated.tsv \
        --repeat-regions 2900-3139 \
        --window 300 \
        --min-homology 25
"""

import argparse
import os
import sys
from collections import defaultdict

import pandas as pd
import pysam


# ---------------------------------------------------------------------------
# Repeat region check (user-specified coordinate ranges on the plasmid)
# ---------------------------------------------------------------------------

def parse_repeat_regions(region_strings):
    """
    Parse repeat region strings like '2900-3139' into a list of
    (start, end) tuples (0-based half-open).
    """
    regions = []
    if not region_strings:
        return regions
    for rs in region_strings:
        parts = rs.split("-")
        if len(parts) != 2:
            print(
                f"[annotate_confidence] WARNING: invalid repeat region '{rs}', "
                f"expected format 'start-end'. Skipping.",
                file=sys.stderr,
            )
            continue
        try:
            start = int(parts[0])
            end = int(parts[1])
            regions.append((start, end))
        except ValueError:
            print(
                f"[annotate_confidence] WARNING: invalid repeat region '{rs}'. Skipping.",
                file=sys.stderr,
            )
    return regions


def position_in_repeat_regions(pos, regions):
    """
    Check if a 0-based position falls within any of the repeat regions.
    Returns True/False.
    """
    if pos is None:
        return False
    for start, end in regions:
        if start <= pos < end:
            return True
    return False


# ---------------------------------------------------------------------------
# Direct homology check using minimap2
# ---------------------------------------------------------------------------

def _check_has_homology(row, ref_fasta, plasmid_index, args):
    """Quick check: does either breakpoint have homology with this index?"""
    for col_chr, col_pos in [
        ("chr_bp_upstream", "pos_bp_upstream"),
        ("chr_bp_downstream", "pos_bp_downstream"),
    ]:
        hit = check_breakpoint_homology(
            row.get(col_chr), row.get(col_pos),
            ref_fasta, plasmid_index,
            window=args.window,
            min_homology=args.min_homology,
            bp_tolerance=args.bp_tolerance,
        )
        if hit and hit["has_homology"] and hit["breakpoint_covered"]:
            return True
    return False


def check_breakpoint_homology(
    chrom, pos, ref_fasta, plasmid_index, window=150, min_homology=25,
    bp_tolerance=0,
):
    """
    Extract genomic sequence at [pos-window, pos+window] and align to the
    plasmid. Returns alignment details if the breakpoint region has
    significant homology to the vector.

    Uses mappy (minimap2 Python API) for fast alignment.

    The alignment must come within bp_tolerance bp of the breakpoint
    position to count as covering it.

    Returns:
        dict with keys: has_homology, plasmid_chr, plasmid_start, plasmid_end,
                        align_len, identity, breakpoint_covered
        or None if no significant homology found.
    """
    if chrom is None or pos is None:
        return None

    try:
        # Extract genomic sequence around the breakpoint:
        chrom_len = ref_fasta.get_reference_length(chrom)
        start = max(0, int(pos) - window)
        end = min(chrom_len, int(pos) + window)
        seq = ref_fasta.fetch(chrom, start, end)
    except (KeyError, ValueError):
        return None

    if not seq or len(seq) < min_homology:
        return None

    # Align the extracted sequence to the plasmid:
    best_hit = None
    best_len = 0

    for hit in plasmid_index.map(seq):
        # hit attributes: ctg, r_st, r_en, q_st, q_en, mapq, NM, strand, etc.
        align_len = hit.r_en - hit.r_st
        if align_len < min_homology:
            continue

        # Check if the alignment is within bp_tolerance of the breakpoint
        # The breakpoint is at position (pos - start) in the extracted sequence
        bp_in_query = int(pos) - start
        query_covers_bp = (
            (hit.q_st - bp_tolerance) <= bp_in_query <= (hit.q_en + bp_tolerance)
        )

        if align_len > best_len:
            best_len = align_len
            mlen = hit.mlen  # number of matching bases
            blen = hit.blen  # alignment block length
            identity = mlen / blen if blen > 0 else 0
            best_hit = {
                "has_homology": True,
                "plasmid_chr": hit.ctg,
                "plasmid_start": hit.r_st,
                "plasmid_end": hit.r_en,
                "align_len": align_len,
                "identity": round(identity, 3),
                "breakpoint_covered": query_covers_bp,
            }

    return best_hit


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Annotate integration sites with confidence flags based on "
            "vector-genome homology and repeat content."
        )
    )
    ap.add_argument(
        "--sites", required=True,
        help="Integration sites TSV (from resolve_integration_sites.py)",
    )
    ap.add_argument(
        "--plasmid", required=True,
        help="Plasmid/vector FASTA file",
    )
    ap.add_argument(
        "--ref-genome", required=True,
        help="Reference genome FASTA (must be indexed with .fai)",
    )
    ap.add_argument(
        "--output", required=True,
        help="Output annotated TSV path",
    )
    ap.add_argument(
        "--repeat-regions", nargs="*", default=None,
        help=(
            "Optional: 0-based coordinate ranges on the plasmid that contain "
            "repeat elements (e.g., Alu). Format: start-end. "
            "Example: --repeat-regions 2900-3139 500-700"
        ),
    )
    ap.add_argument(
        "--window", type=int, default=300,
        help="Window size (bp) around breakpoint for homology check (default: 300)",
    )
    ap.add_argument(
        "--min-homology", type=int, default=25,
        help="Minimum alignment length (bp) to flag as homology (default: 25)",
    )
    ap.add_argument(
        "--bp-tolerance", type=int, default=50,
        help=(
            "Tolerance (bp) for breakpoint coverage check. The alignment to "
            "the plasmid must come within this distance of the breakpoint "
            "to count as covering it (default: 50)"
        ),
    )
    args = ap.parse_args()

    # -- Import mappy (minimap2 Python API) --------------------------------
    try:
        import mappy
    except ImportError:
        print(
            "[annotate_confidence] ERROR: 'mappy' package not found. "
            "Install with: pip install mappy",
            file=sys.stderr,
        )
        sys.exit(1)

    # -- Load integration sites --------------------------------------------
    sites_df = pd.read_csv(args.sites, sep="\t")
    n_total = len(sites_df)
    print(f"[annotate_confidence] {n_total} integration sites loaded", file=sys.stderr)

    if n_total == 0:
        sites_df["homology_flag"] = pd.Series(dtype=str)
        sites_df["repeat_flag"] = pd.Series(dtype=str)
        sites_df["confidence"] = pd.Series(dtype=str)
        sites_df["homology_detail"] = pd.Series(dtype=str)
        sites_df.to_csv(args.output, sep="\t", index=False)
        print("[annotate_confidence] No sites to annotate. Done.", file=sys.stderr)
        return

    # -- Build minimap2 indices for the plasmid -----------------------------
    # Tier 1: preset="sr" (k=15) for all sites — catches high-identity homology
    # Tier 2: k=11 (sensitive) only for repeat_mediated sites — confirms
    #         divergent Alu homology that sr cannot seed
    print("[annotate_confidence] Building plasmid indices ...", file=sys.stderr)
    plasmid_index_sr = mappy.Aligner(args.plasmid, preset="sr", best_n=5)
    plasmid_index_sensitive = mappy.Aligner(args.plasmid, k=11, w=5, best_n=5)
    if not plasmid_index_sr or not plasmid_index_sensitive:
        print(
            f"[annotate_confidence] ERROR: could not build index for {args.plasmid}",
            file=sys.stderr,
        )
        sys.exit(1)

    # -- Open reference genome FASTA ---------------------------------------
    ref_fasta = pysam.FastaFile(args.ref_genome)

    # -- Parse repeat regions (optional) ------------------------------------
    repeat_regions = parse_repeat_regions(args.repeat_regions)
    if repeat_regions:
        print(
            f"[annotate_confidence] {len(repeat_regions)} repeat region(s) on plasmid: "
            + ", ".join(f"{s}-{e}" for s, e in repeat_regions),
            file=sys.stderr,
        )

    # -- Annotate each integration site ------------------------------------
    homology_flags = []
    repeat_flags = []
    confidence_values = []
    homology_details = []

    for idx, row in sites_df.iterrows():
        # --- Repeat region check (done FIRST to guide tier selection) ---
        repeat_flag = ""
        if repeat_regions:
            aav_bp_up = row.get("aav_bp_upstream")
            aav_bp_down = row.get("aav_bp_downstream")

            up_in_repeat = pd.notna(aav_bp_up) and position_in_repeat_regions(
                int(aav_bp_up), repeat_regions
            )
            down_in_repeat = pd.notna(aav_bp_down) and position_in_repeat_regions(
                int(aav_bp_down), repeat_regions
            )

            if up_in_repeat or down_in_repeat:
                repeat_flag = "repeat_mediated"

        # --- Direct homology check (two-tier) ---
        # Tier 1: preset="sr" for all sites
        # Tier 2: k=11 for repeat_mediated sites (confirms divergent Alu homology)
        plasmid_index = plasmid_index_sr
        if repeat_flag and not _check_has_homology(
            row, ref_fasta, plasmid_index_sr, args
        ):
            # Tier 1 missed it; try sensitive index for repeat_mediated sites
            plasmid_index = plasmid_index_sensitive

        up_hit = check_breakpoint_homology(
            row.get("chr_bp_upstream"),
            row.get("pos_bp_upstream"),
            ref_fasta, plasmid_index,
            window=args.window,
            min_homology=args.min_homology,
            bp_tolerance=args.bp_tolerance,
        )
        down_hit = check_breakpoint_homology(
            row.get("chr_bp_downstream"),
            row.get("pos_bp_downstream"),
            ref_fasta, plasmid_index,
            window=args.window,
            min_homology=args.min_homology,
            bp_tolerance=args.bp_tolerance,
        )

        # Determine homology flag:
        has_homology = False
        detail_parts = []

        for label, hit in [("upstream", up_hit), ("downstream", down_hit)]:
            if hit and hit["has_homology"] and hit["breakpoint_covered"]:
                has_homology = True
                detail_parts.append(
                    f"{label}:{hit['plasmid_chr']}:{hit['plasmid_start']}-"
                    f"{hit['plasmid_end']}({hit['align_len']}bp,"
                    f"{hit['identity']*100:.0f}%id)"
                )

        homology_flag = "vector_homology" if has_homology else ""
        homology_detail = ";".join(detail_parts) if detail_parts else ""

        # --- Add repeat region info to detail ---
        if repeat_flag:
            aav_bp_up = row.get("aav_bp_upstream")
            aav_bp_down = row.get("aav_bp_downstream")
            up_in_repeat = pd.notna(aav_bp_up) and position_in_repeat_regions(
                int(aav_bp_up), repeat_regions
            )
            bp_val = int(aav_bp_up) if up_in_repeat else int(row.get("aav_bp_downstream"))
            repeat_info = ""
            for s, e in repeat_regions:
                if s <= bp_val < e:
                    repeat_info = f"plasmid_repeat:{s}-{e}"
                    break
            if homology_detail:
                homology_detail += f";{repeat_info}"
            else:
                homology_detail = repeat_info

        # --- Confidence assignment ---
        # "low": vector_homology from Tier 1 (high-identity, e.g. HBB, GH cluster)
        # "caution": repeat_mediated confirmed by Tier 2 (divergent Alu homology)
        # "high": everything else (including repeat_mediated without homology)
        if has_homology and not repeat_flag:
            # Pure vector homology (HBB, GH cluster) — clear artifact
            confidence = "low"
        elif has_homology and repeat_flag and plasmid_index == plasmid_index_sr:
            # High-identity Alu caught by Tier 1 — clear artifact
            confidence = "low"
        elif has_homology and repeat_flag and plasmid_index == plasmid_index_sensitive:
            # Divergent Alu confirmed by Tier 2 — likely artifact
            confidence = "caution"
        else:
            confidence = "high"

        homology_flags.append(homology_flag)
        repeat_flags.append(repeat_flag)
        confidence_values.append(confidence)
        homology_details.append(homology_detail)

    # -- Add columns to dataframe ------------------------------------------
    sites_df["confidence"] = confidence_values
    sites_df["homology_flag"] = homology_flags
    sites_df["repeat_flag"] = repeat_flags
    sites_df["homology_detail"] = homology_details

    # -- Summary -----------------------------------------------------------
    n_low = sum(1 for c in confidence_values if c == "low")
    n_caution = sum(1 for c in confidence_values if c == "caution")
    n_high = sum(1 for c in confidence_values if c == "high")
    n_homology = sum(1 for f in homology_flags if f)
    n_repeat = sum(1 for f in repeat_flags if f)

    print(
        f"[annotate_confidence] Results: {n_high} high, "
        f"{n_caution} caution, {n_low} low confidence sites",
        file=sys.stderr,
    )
    print(
        f"[annotate_confidence]   - {n_homology} flagged with vector homology",
        file=sys.stderr,
    )
    if repeat_regions:
        print(
            f"[annotate_confidence]   - {n_repeat} flagged as repeat-mediated",
            file=sys.stderr,
        )

    # -- Write output ------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    sites_df.to_csv(args.output, sep="\t", index=False)
    print(f"[annotate_confidence] Written to {args.output}", file=sys.stderr)

    ref_fasta.close()


if __name__ == "__main__":
    main()
