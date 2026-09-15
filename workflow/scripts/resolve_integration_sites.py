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
resolve_integration_sites.py

Takes the output of step 4 (identified_insertions.txt) and the full
alignment BAM (03_aligned/both/{sample}.bam) and resolves each insertion
cluster into precise chromosomal and AAV breakpoints using split-read
evidence.

The key insight: reads that span an AAV integration site appear as
chimeric alignments (primary + supplementary records) in the BAM.  For
each read the segments are ordered along the read's 5'→3' axis; a
genome↔AAV transition marks a breakpoint.

Before BAM interrogation, step-4 clusters are transitively merged using
union-find: any two clusters that share even one read ID are joined into
a single site.  This mirrors the igraph-based merging in the original R
implementation and correctly handles reads that span large integrations
and therefore land in different bedtools-merge windows.

Breakpoints are classified as upstream (lower genomic coordinate) and
downstream (higher genomic coordinate) regardless of read strand.

Output TSV columns
------------------
site_id              serial number (1-based)
chr_bp_upstream      chromosome of upstream genomic breakpoint
pos_bp_upstream      0-based position of upstream genomic breakpoint
chr_bp_downstream    chromosome of downstream genomic breakpoint
pos_bp_downstream    0-based position of downstream genomic breakpoint
aav_bp_upstream      AAV coordinate at upstream junction
aav_bp_downstream    AAV coordinate at downstream junction
count                supporting read count (UMI-deduped when RX tag present)
aav_insertion_len    sum of AAV segment lengths per read (median across reads);
                     NA when only one junction is detected
is_translocation     True when breakpoints are on different chromosomes
n_clusters_merged    number of step-4 clusters collapsed into this site
cluster_chr          semicolon-joined chromosome(s) of source clusters
cluster_start        min start across source clusters
cluster_end          max end across source clusters
bp_evidence          breakpoint evidence quality: 'both_junctions',
                     'upstream_only', 'downstream_only', or 'no_junction'
bp_upstream_spread   IQR (bp) of upstream breakpoint positions across reads
bp_downstream_spread IQR (bp) of downstream breakpoint positions across reads
orientation          AAV insertion direction relative to genome:
                     'fwd', 'rev', or 'mixed'
"""

import argparse
import os
import statistics
import sys
from collections import Counter, defaultdict

import pandas as pd
import pysam


# ---------------------------------------------------------------------------
# Transitive cluster merging (union-find)
# ---------------------------------------------------------------------------

class _UnionFind:
    """Simple union-find with path compression."""
    def __init__(self):
        self._parent: dict = {}

    def find(self, x):
        if x not in self._parent:
            self._parent[x] = x
        while self._parent[x] != x:
            self._parent[x] = self._parent[self._parent[x]]  # path halving
            x = self._parent[x]
        return x

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx != ry:
            self._parent[ry] = rx


def merge_clusters_by_shared_reads(
    ins_df: pd.DataFrame,
    aav_chrs: set,
    min_cross_chr_reads: int = 2,
    max_same_chr_distance: int = 100_000,
) -> list:
    """
    Transitively merge step-4 clusters that share read IDs, using a
    two-phase strategy that distinguishes confident merges from ambiguous
    ones.

    Phase 1 — nearby same-chromosome merging (aggressive):
    Clusters on the same chromosome whose intervals are within
    max_same_chr_distance bp of each other are merged if they share
    >=1 read ID.

    Phase 2 — distant / cross-chromosome merging (conservative):
    All remaining component pairs are merged only if they share
    >= min_cross_chr_reads reads.

    Rows whose chromosome is an AAV reference are excluded from clustering.

    Returns a list of dicts:
        read_ids          – set of all read IDs in the merged site
        source_intervals  – list of (chrom, start, end) for BAM fetching
        cluster_chr       – semicolon-joined unique chromosomes
        cluster_start     – min start across source intervals
        cluster_end       – max end across source intervals
        is_translocation  – True when >1 distinct genomic chromosome
        n_clusters_merged – number of step-4 rows collapsed
    """
    genomic_mask = ~ins_df["chrom"].isin(aav_chrs)

    row_reads: dict = {}
    row_chrom: dict = {}
    row_start: dict = {}
    row_end: dict = {}
    for i, row in ins_df[genomic_mask].iterrows():
        reads = [r.strip() for r in str(row["read_ids"]).split(",") if r.strip()]
        row_reads[i] = reads
        row_chrom[i] = str(row["chrom"])
        row_start[i] = int(row["start"])
        row_end[i] = int(row["end"])

    # Phase 1: same-chromosome, nearby clusters
    uf = _UnionFind()
    read_to_rows: dict = defaultdict(list)

    for i, reads in row_reads.items():
        for read_id in reads:
            read_to_rows[read_id].append(i)

    for read_id, row_indices in read_to_rows.items():
        if len(row_indices) < 2:
            continue
        for a_idx in range(len(row_indices)):
            for b_idx in range(a_idx + 1, len(row_indices)):
                a, b = row_indices[a_idx], row_indices[b_idx]
                if row_chrom[a] == row_chrom[b]:
                    gap = max(0, max(row_start[a], row_start[b])
                              - min(row_end[a], row_end[b]))
                    if gap <= max_same_chr_distance:
                        uf.union(a, b)

    # Build Phase-1 components
    phase1_groups: dict = defaultdict(list)
    for i in row_reads:
        phase1_groups[uf.find(i)].append(i)

    comp_reads: dict = {}
    for root, indices in phase1_groups.items():
        reads = set()
        for i in indices:
            reads.update(row_reads[i])
        comp_reads[root] = reads

    # Phase 2: cross-chromosome / distant merging
    if min_cross_chr_reads <= 1:
        for read_id, row_indices in read_to_rows.items():
            for idx in row_indices[1:]:
                uf.union(row_indices[0], idx)
    else:
        roots = list(phase1_groups.keys())
        changed = True
        while changed:
            changed = False
            for a_idx in range(len(roots)):
                for b_idx in range(a_idx + 1, len(roots)):
                    ra_root = uf.find(roots[a_idx])
                    rb_root = uf.find(roots[b_idx])
                    if ra_root == rb_root:
                        continue
                    reads_a = comp_reads.get(ra_root, set())
                    reads_b = comp_reads.get(rb_root, set())
                    shared = len(reads_a & reads_b)
                    if shared >= min_cross_chr_reads:
                        uf.union(ra_root, rb_root)
                        new_root = uf.find(ra_root)
                        pooled = reads_a | reads_b
                        comp_reads.pop(ra_root, None)
                        comp_reads.pop(rb_root, None)
                        comp_reads[new_root] = pooled
                        changed = True

    # Build final groups
    final_groups: dict = defaultdict(list)
    for i in row_reads:
        final_groups[uf.find(i)].append(i)

    merged = []
    for row_indices in final_groups.values():
        all_reads: set = set()
        chroms: set = set()
        starts: list = []
        ends: list = []
        source_intervals: list = []

        for i in row_indices:
            row = ins_df.loc[i]
            all_reads.update(row_reads[i])
            chrom = str(row["chrom"])
            if chrom not in aav_chrs:
                chroms.add(chrom)
            starts.append(int(row["start"]))
            ends.append(int(row["end"]))
            source_intervals.append((chrom, int(row["start"]), int(row["end"])))

        merged.append({
            "read_ids":          all_reads,
            "source_intervals":  source_intervals,
            "cluster_chr":       ";".join(sorted(chroms)),
            "cluster_start":     min(starts),
            "cluster_end":       max(ends),
            "is_translocation":  len(chroms) > 1,
            "n_clusters_merged": len(row_indices),
        })

    return merged

    # ---------------------------------------------------------------------------
# Query-coordinate utilities
# ---------------------------------------------------------------------------

def query_coords_5prime(rec: pysam.AlignedSegment):
    """
    Return (q_start, q_end) of an alignment in the 5'→3' coordinates of
    the *original* full-length read (re-adding hard-clipped bases).

    Works correctly for primary and hard-clipped supplementary records on
    either strand.
    """
    cig = rec.cigartuples
    if not cig:
        return None, None

    leading_h  = cig[0][1]  if cig[0][0]  == 5 else 0
    trailing_h = cig[-1][1] if cig[-1][0] == 5 else 0

    total = rec.infer_read_length()
    if total is None:
        return None, None

    stored_len = total - leading_h - trailing_h
    if stored_len <= 0:
        return None, None

    q_s = rec.query_alignment_start
    q_e = rec.query_alignment_end

    if rec.is_reverse:
        return trailing_h + (stored_len - q_e), trailing_h + (stored_len - q_s)
    else:
        return leading_h + q_s, leading_h + q_e


# ---------------------------------------------------------------------------
# Non-overlapping segment chain selection (NEW FIX)
# ---------------------------------------------------------------------------

def _select_nonoverlapping_chain(segs, max_query_overlap_frac=0.5):
    """
    From all alignment segments of a read, select a non-overlapping set
    of segments in query (read) coordinates.

    Strategy: greedy, longest-first.  A candidate segment is rejected if
    it overlaps an already-selected segment by more than
    max_query_overlap_frac of the SHORTER segment's query span.

    This removes spurious supplementary alignments where the aligner
    re-maps a portion of the read that is already covered by a longer,
    more reliable alignment.

    Parameters
    ----------
    segs : list of dicts
        Each dict must have 'qs' and 'qe' (query start/end in 5'->3'
        coordinates), plus all other fields needed downstream.
    max_query_overlap_frac : float
        Maximum allowed overlap as a fraction of the shorter segment.
        Default 0.5 means a segment is rejected if >50% of its query
        span is already covered by a selected segment.

    Returns
    -------
    list of dicts - the selected non-overlapping segments, sorted by qs.
    """
    if not segs:
        return []

    # Sort by query span descending (longest = most reliable first),
    # then by MAPQ descending as tiebreaker
    segs_sorted = sorted(
        segs,
        key=lambda s: (s["qe"] - s["qs"], s.get("mapq", 0)),
        reverse=True,
    )

    selected = []
    for candidate in segs_sorted:
        c_qs, c_qe = candidate["qs"], candidate["qe"]
        c_len = c_qe - c_qs
        if c_len <= 0:
            continue

        dominated = False
        for kept in selected:
            k_qs, k_qe = kept["qs"], kept["qe"]
            k_len = k_qe - k_qs

            # Compute overlap in query coordinates
            overlap_start = max(c_qs, k_qs)
            overlap_end = min(c_qe, k_qe)
            overlap = max(0, overlap_end - overlap_start)

            if overlap <= 0:
                continue

            # Fraction of the SHORTER segment that is overlapped
            shorter_len = min(c_len, k_len)
            if shorter_len <= 0:
                continue
            overlap_frac = overlap / shorter_len

            if overlap_frac > max_query_overlap_frac:
                dominated = True
                break

        if not dominated:
            selected.append(candidate)

    # Return sorted by query position for junction detection
    selected.sort(key=lambda s: s["qs"])
    return selected


# ---------------------------------------------------------------------------
# Junction detection (FIXED: chain selection + correct edge logic)
# ---------------------------------------------------------------------------

def find_junctions(
    records: list,
    aav_chrs: set,
) -> tuple:
    """
    Given all aligned records for one read (primary + supplementary),
    detect genome<->AAV boundary transitions and return:

        junctions       - list of dicts, each with:
                            genome_chr:  chromosome of genomic breakpoint
                            genome_pos:  0-based genomic breakpoint position
                            aav_pos:     AAV coordinate at junction
                            ori:         'fwd' or 'rev' (AAV orientation
                                         relative to genome)
        total_aav_len   - sum of reference bases covered by every AAV
                          segment in this read.

    Key insight for breakpoint coordinates:

    Segments are ordered along the read 5'->3'. At a junction between
    segment A (left in read) and segment B (right in read), the breakpoint
    is at:
      - Segment A: the edge corresponding to its RIGHT end in read coords.
          Forward-mapped (read and ref go same direction):
              right end in read = high ref coord = reference_end
          Reverse-mapped (read and ref go opposite directions):
              right end in read = low ref coord = reference_start
      - Segment B: the edge corresponding to its LEFT end in read coords.
          Forward-mapped: left end in read = low ref coord = reference_start
          Reverse-mapped: left end in read = high ref coord = reference_end
    """
    segs = []
    for rec in records:
        if rec.is_unmapped:
            continue
        qs, qe = query_coords_5prime(rec)
        if qs is None:
            continue
        segs.append({
            "qs":     qs,
            "qe":     qe,
            "chr":    rec.reference_name,
            "rs":     rec.reference_start,
            "re":     rec.reference_end,
            "is_aav": rec.reference_name in aav_chrs,
            "is_rev": rec.is_reverse,
            "mapq":   rec.mapping_quality,
        })

    if not segs:
        return [], 0

    # ── KEY FIX: remove overlapping/redundant segments ────────────────
    segs = _select_nonoverlapping_chain(segs)

    if not segs:
        return [], 0

    # Sum the reference length of every AAV segment
    total_aav_len = sum(s["re"] - s["rs"] for s in segs if s["is_aav"])

    junctions = []

    for i in range(len(segs) - 1):
        a = segs[i]      # left segment in read (5' side of junction)
        b = segs[i + 1]  # right segment in read (3' side of junction)

        # We only care about genome<->AAV transitions
        if a["is_aav"] == b["is_aav"]:
            continue

        # ── Identify which is genome and which is AAV ─────────────────
        if not a["is_aav"] and b["is_aav"]:
            # genome (a) -> AAV (b): read EXITS genome, ENTERS AAV
            genome_seg = a
            aav_seg = b
            genome_is_left = True   # genome is on LEFT side of junction in read

        elif a["is_aav"] and not b["is_aav"]:
            # AAV (a) -> genome (b): read EXITS AAV, ENTERS genome
            genome_seg = b
            aav_seg = a
            genome_is_left = False  # genome is on RIGHT side of junction in read

        else:
            continue

        # ── Genomic breakpoint ────────────────────────────────────────
        # The breakpoint is at the edge of the genomic alignment that is
        # closest to the junction in READ coordinates.
        if genome_is_left:
            # Genome is LEFT of junction -> breakpoint at its RIGHT edge in read
            # Forward: right edge in read = reference_end
            # Reverse: right edge in read = reference_start
            if genome_seg["is_rev"]:
                genome_pos = genome_seg["rs"]
            else:
                genome_pos = genome_seg["re"]
        else:
            # Genome is RIGHT of junction -> breakpoint at its LEFT edge in read
            # Forward: left edge in read = reference_start
            # Reverse: left edge in read = reference_end
            if genome_seg["is_rev"]:
                genome_pos = genome_seg["re"]
            else:
                genome_pos = genome_seg["rs"]

        # ── AAV breakpoint ────────────────────────────────────────────
        # Same logic applied to the AAV segment.
        if genome_is_left:
            # AAV is RIGHT of junction -> breakpoint at its LEFT edge in read
            # Forward: left edge in read = reference_start
            # Reverse: left edge in read = reference_end
            if aav_seg["is_rev"]:
                aav_pos = aav_seg["re"]
            else:
                aav_pos = aav_seg["rs"]
        else:
            # AAV is LEFT of junction -> breakpoint at its RIGHT edge in read
            # Forward: right edge in read = reference_end
            # Reverse: right edge in read = reference_start
            if aav_seg["is_rev"]:
                aav_pos = aav_seg["rs"]
            else:
                aav_pos = aav_seg["re"]

        # ── Orientation ───────────────────────────────────────────────
        # "fwd" = AAV inserted in same direction as genome at this locus
        #         (both segments have same is_rev value)
        # "rev" = AAV inserted in reverse complement
        #         (segments have opposite is_rev values)
        ori = "fwd" if genome_seg["is_rev"] == aav_seg["is_rev"] else "rev"

        # ── Junction type ─────────────────────────────────────────────
        # Upstream boundary: genome breakpoint is at the HIGH end
        #   (reference_end) of the genomic segment — it's the rightmost
        #   base of the left-flanking genomic sequence.
        # Downstream boundary: genome breakpoint is at the LOW end
        #   (reference_start) of the genomic segment — it's the leftmost
        #   base of the right-flanking genomic sequence.
        #
        # The XOR of genome_is_left and genome_seg["is_rev"] correctly
        # identifies this regardless of read strand:
        #   fwd read + genome_is_left  → pos=re (HIGH) → upstream
        #   rev read + genome_is_right → pos=re (HIGH) → upstream
        #   fwd read + genome_is_right → pos=rs (LOW)  → downstream
        #   rev read + genome_is_left  → pos=rs (LOW)  → downstream
        jtype = "upstream" if (genome_is_left != genome_seg["is_rev"]) else "downstream"

        junctions.append({
            "genome_chr":    genome_seg["chr"],
            "genome_pos":    genome_pos,
            "aav_pos":       aav_pos,
            "ori":           ori,
            "jtype":         jtype,
            "genome_is_left": genome_is_left,
        })

    return junctions, total_aav_len

# ---------------------------------------------------------------------------
# Aggregation helpers
# ---------------------------------------------------------------------------

def _median_int(vals: list):
    vals = [v for v in vals if v is not None]
    return int(statistics.median(vals)) if vals else None


def _mode_str(vals: list):
    vals = [v for v in vals if v]
    return Counter(vals).most_common(1)[0][0] if vals else None


def _iqr(vals: list):
    """Interquartile range; returns None when fewer than 2 values."""
    vals = sorted(v for v in vals if v is not None)
    n = len(vals)
    if n < 2:
        return None
    q1 = vals[n // 4]
    q3 = vals[(3 * n) // 4]
    return q3 - q1


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Resolve precise AAV integration-site breakpoints from "
            "split-read evidence in the 03_aligned/both BAM."
        )
    )
    ap.add_argument(
        "--insertions", required=True,
        help="Step-4 output file (identified_insertions.txt)",
    )
    ap.add_argument(
        "--bam", required=True,
        help="Path to 03_aligned/both/{sample}.bam (must be indexed with .bai)",
    )
    ap.add_argument(
        "--aav-chr", required=True, nargs="+",
        help="AAV chromosome name(s) in the custom reference",
    )
    ap.add_argument(
        "--output", required=True,
        help="Output TSV path",
    )
    ap.add_argument(
        "--min-reads", type=int, default=1,
        help="Minimum supporting reads to report a site (default: 1)",
    )
    ap.add_argument(
        "--min-cross-chr-reads", type=int, default=2,
        help=(
            "Minimum shared reads required to merge clusters on different "
            "chromosomes (default: 2)."
        ),
    )
    ap.add_argument(
        "--max-same-chr-distance", type=int, default=100_000,
        help=(
            "Maximum gap (bp) between two same-chromosome clusters for "
            "Phase-1 (single-read) merging (default: 100000)."
        ),
    )
    args = ap.parse_args()

    aav_chrs = set(args.aav_chr)

    # -- Load step-4 clusters ---------------------------------------------
    ins_cols = ["chrom", "start", "end", "count", "mapq", "strand", "read_ids"]
    ins_df = pd.read_csv(
        args.insertions, sep="\t", header=None, names=ins_cols,
    )
    print(
        f"[resolve] {len(ins_df)} insertion clusters loaded from {args.insertions}",
        file=sys.stderr,
    )

    # -- Transitively merge clusters that share read IDs -------------------
    merged_sites = merge_clusters_by_shared_reads(
        ins_df, aav_chrs,
        min_cross_chr_reads=args.min_cross_chr_reads,
        max_same_chr_distance=args.max_same_chr_distance,
    )
    n_merged = len(ins_df) - len(merged_sites)
    print(
        f"[resolve] {len(merged_sites)} sites after merging "
        f"({n_merged} cluster(s) collapsed into larger sites)",
        file=sys.stderr,
    )
    n_trans = sum(1 for s in merged_sites if s["is_translocation"])
    if n_trans:
        print(f"[resolve] {n_trans} putative translocation site(s) detected "
              f"(clusters spanning >1 chromosome)", file=sys.stderr)

    # -- Filter chimera-artifact singletons --------------------------------
    _CHIMERA_MAX_READS = 2
    _CHIMERA_SIZE_RATIO = 10

    read_to_site_sizes: dict = defaultdict(list)
    for idx, site in enumerate(merged_sites):
        n = len(site["read_ids"])
        for rid in site["read_ids"]:
            read_to_site_sizes[rid].append((idx, n))

    chimera_artifact_indices: set = set()
    for idx, site in enumerate(merged_sites):
        n_site_reads = len(site["read_ids"])
        if n_site_reads > _CHIMERA_MAX_READS:
            continue
        all_subsumed = True
        for rid in site["read_ids"]:
            max_other = max(
                (sz for j, sz in read_to_site_sizes[rid] if j != idx),
                default=0,
            )
            if max_other < n_site_reads * _CHIMERA_SIZE_RATIO:
                all_subsumed = False
                break
        if all_subsumed:
            chimera_artifact_indices.add(idx)

    if chimera_artifact_indices:
        print(
            f"[resolve] Skipping {len(chimera_artifact_indices)} likely "
            f"chimera-artifact site(s) (<={_CHIMERA_MAX_READS} reads, all "
            f"subsumed by >={_CHIMERA_SIZE_RATIO}x larger site)",
            file=sys.stderr,
        )

    # Collect ALL read IDs for pre-scanning AAV chromosomes.
    all_read_ids: set = set()
    for site in merged_sites:
        all_read_ids.update(site["read_ids"])

    bam = pysam.AlignmentFile(args.bam, "rb")

    # -- Pre-scan AAV chromosome(s) once -----------------------------------
    print("[resolve] Pre-scanning AAV chromosome(s) for split-read partners ...",
          file=sys.stderr)
    aav_recs_by_id: dict = defaultdict(list)
    for aav_chr in aav_chrs:
        try:
            for rec in bam.fetch(aav_chr):
                if rec.query_name in all_read_ids:
                    aav_recs_by_id[rec.query_name].append(rec)
        except (ValueError, KeyError):
            print(
                f"[resolve] WARNING: could not fetch '{aav_chr}' from BAM",
                file=sys.stderr,
            )

    # -- Process each merged site ------------------------------------------
    rows = []
    site_id = 0

    for site_idx, site in enumerate(merged_sites):
        if site_idx in chimera_artifact_indices:
            continue
        read_ids = site["read_ids"]

        # Fetch genomic-side records from every source interval.
        genomic_recs: dict = defaultdict(list)
        for chrom, start, end in site["source_intervals"]:
            try:
                for rec in bam.fetch(chrom, start, end):
                    if rec.query_name in read_ids:
                        genomic_recs[rec.query_name].append(rec)
            except (ValueError, KeyError):
                pass

        # Detect UMI tag (RX:Z:)
        use_umi = any(
            rec.has_tag("RX")
            for recs in genomic_recs.values()
            for rec in recs
        )

        # -- Resolve junctions per read ------------------------------------
        all_junctions: list = []    # list of junction dicts across all reads
        all_aav_lens: list = []     # per-read total AAV segment lengths
        contributing_reads: set = set()
        seen_umis: set = set()

        for read_id in read_ids:
            # Combine genomic and AAV records, deduplicating
            recs = list(genomic_recs.get(read_id, []))
            seen_keys = {
                (r.reference_name, r.reference_start, r.reference_end,
                 r.is_reverse, r.is_supplementary)
                for r in recs
            }
            for r in aav_recs_by_id.get(read_id, []):
                key = (r.reference_name, r.reference_start,
                       r.reference_end, r.is_reverse, r.is_supplementary)
                if key not in seen_keys:
                    recs.append(r)
                    seen_keys.add(key)

            # Drop secondary alignments (keep primary + supplementary)
            recs = [r for r in recs if not r.is_secondary]
            if not recs:
                continue

            # UMI deduplication
            if use_umi:
                umi = next(
                    (r.get_tag("RX") for r in recs if r.has_tag("RX")),
                    None,
                )
                if umi in seen_umis:
                    continue
                if umi:
                    seen_umis.add(umi)

            junctions, aav_len = find_junctions(recs, aav_chrs)

            if junctions:
                contributing_reads.add(read_id)
                all_junctions.extend(junctions)
                if aav_len > 0:
                    all_aav_lens.append(aav_len)

        # Count: UMI-deduped if tags present, else unique read IDs
        n_reads = len(seen_umis) if use_umi else len(contributing_reads)
        if n_reads < args.min_reads:
            continue

        # -- Assign upstream/downstream by genomic position ----------------
        if not all_junctions:
            # No junction evidence at all
            site_id += 1
            rows.append({
                "site_id":               site_id,
                "chr_bp_upstream":       None,
                "pos_bp_upstream":       None,
                "chr_bp_downstream":     None,
                "pos_bp_downstream":     None,
                "aav_bp_upstream":       None,
                "aav_bp_downstream":     None,
                "count":                 n_reads,
                "aav_insertion_len":     None,
                "is_translocation":      site["is_translocation"],
                "n_clusters_merged":     site["n_clusters_merged"],
                "cluster_chr":           site["cluster_chr"],
                "cluster_start":         site["cluster_start"],
                "cluster_end":           site["cluster_end"],
                "bp_evidence":           "no_junction",
                "bp_upstream_spread":    None,
                "bp_downstream_spread":  None,
                "orientation":           None,
            })
            continue

        # Cluster junctions into upstream and downstream groups using
        # junction type, which accounts for read strand via the XOR of
        # genome_is_left and genome_seg["is_rev"].
        upstream_jns = [j for j in all_junctions if j["jtype"] == "upstream"]
        downstream_jns = [j for j in all_junctions if j["jtype"] == "downstream"]

        # Handle translocations: the jtype formula can assign the same
        # type to both boundaries when one flank is forward-strand and
        # the other is reverse-strand (both end up at reference_end).
        # Detect this by checking if one group spans multiple chromosomes
        # and split by genome_is_left (entry vs exit side of AAV) instead.
        if upstream_jns and not downstream_jns:
            chroms = set(j["genome_chr"] for j in upstream_jns)
            if len(chroms) > 1:
                entry_jns = [j for j in upstream_jns if j["genome_is_left"]]
                exit_jns = [j for j in upstream_jns if not j["genome_is_left"]]
                if entry_jns and exit_jns:
                    upstream_jns = entry_jns
                    downstream_jns = exit_jns
        elif downstream_jns and not upstream_jns:
            chroms = set(j["genome_chr"] for j in downstream_jns)
            if len(chroms) > 1:
                entry_jns = [j for j in downstream_jns if j["genome_is_left"]]
                exit_jns = [j for j in downstream_jns if not j["genome_is_left"]]
                if entry_jns and exit_jns:
                    upstream_jns = entry_jns
                    downstream_jns = exit_jns

        # Compute breakpoints
        chr_bp_up = _mode_str([j["genome_chr"] for j in upstream_jns]) if upstream_jns else None
        pos_bp_up = _median_int([j["genome_pos"] for j in upstream_jns]) if upstream_jns else None
        aav_bp_up = _median_int([j["aav_pos"] for j in upstream_jns]) if upstream_jns else None

        chr_bp_down = _mode_str([j["genome_chr"] for j in downstream_jns]) if downstream_jns else None
        pos_bp_down = _median_int([j["genome_pos"] for j in downstream_jns]) if downstream_jns else None
        aav_bp_down = _median_int([j["aav_pos"] for j in downstream_jns]) if downstream_jns else None

        # -- Determine bp_evidence ----------------------------------------
        has_upstream = pos_bp_up is not None
        has_downstream = pos_bp_down is not None

        if has_upstream and has_downstream:
            bp_evidence = "both_junctions"
        elif has_upstream:
            bp_evidence = "upstream_only"
        elif has_downstream:
            bp_evidence = "downstream_only"
        else:
            bp_evidence = "no_junction"

        # -- aav_insertion_len: only when both junctions detected ----------
        if has_upstream and has_downstream:
            aav_ins_len = _median_int(all_aav_lens) if all_aav_lens else None
        else:
            aav_ins_len = None

        # -- Breakpoint spread (IQR) --------------------------------------
        bp_up_spread = _iqr([j["genome_pos"] for j in upstream_jns]) if upstream_jns else None
        bp_down_spread = _iqr([j["genome_pos"] for j in downstream_jns]) if downstream_jns else None

        # -- Orientation: majority vote across all junctions ---------------
        all_oris = [j["ori"] for j in all_junctions]
        if all_oris:
            ori_counts = Counter(all_oris)
            top_ori, top_n = ori_counts.most_common(1)[0]
            if top_n == sum(ori_counts.values()):
                orientation = top_ori
            else:
                orientation = "both"
        else:
            orientation = None

        # -- Translocation: resolved from breakpoint chromosomes -----------
        if chr_bp_up is not None and chr_bp_down is not None:
            resolved_is_translocation = (chr_bp_up != chr_bp_down)
        else:
            resolved_is_translocation = site["is_translocation"]

        site_id += 1
        rows.append({
            "site_id":               site_id,
            "chr_bp_upstream":       chr_bp_up,
            "pos_bp_upstream":       pos_bp_up,
            "chr_bp_downstream":     chr_bp_down,
            "pos_bp_downstream":     pos_bp_down,
            "aav_bp_upstream":       aav_bp_up,
            "aav_bp_downstream":     aav_bp_down,
            "count":                 n_reads,
            "aav_insertion_len":     aav_ins_len,
            "is_translocation":      resolved_is_translocation,
            "n_clusters_merged":     site["n_clusters_merged"],
            "cluster_chr":           site["cluster_chr"],
            "cluster_start":         site["cluster_start"],
            "cluster_end":           site["cluster_end"],
            "bp_evidence":           bp_evidence,
            "bp_upstream_spread":    bp_up_spread,
            "bp_downstream_spread":  bp_down_spread,
            "orientation":           orientation,
        })

    bam.close()

    # -- Write output ------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    out_df = pd.DataFrame(rows, columns=[
        "site_id", "chr_bp_upstream", "pos_bp_upstream",
        "chr_bp_downstream", "pos_bp_downstream",
        "aav_bp_upstream", "aav_bp_downstream",
        "count", "aav_insertion_len",
        "is_translocation", "n_clusters_merged",
        "cluster_chr", "cluster_start", "cluster_end",
        "bp_evidence",
        "bp_upstream_spread", "bp_downstream_spread",
        "orientation",
    ])
    # Use nullable integer dtype so integer columns with NaN
    # render without trailing '.0'
    for col in ["pos_bp_upstream", "pos_bp_downstream",
                "aav_bp_upstream", "aav_bp_downstream",
                "aav_insertion_len",
                "bp_upstream_spread", "bp_downstream_spread"]:
        out_df[col] = out_df[col].astype(pd.Int64Dtype())
    out_df.to_csv(args.output, sep="\t", index=False)
    print(
        f"[resolve] {len(out_df)} integration sites written to {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
