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

import sys
import gzip

def die(msg, code=1):
    sys.stderr.write(msg + "\n")
    sys.exit(code)

def parse_args():
    if len(sys.argv) != 6:
        die(
            "Usage: python merge_bed_tpm.py <BED.gz> <TPM.tsv.gz> <X:int> <scale:float>\n"
            "  - BED.gz: gzipped BED file; transcript ID in column 4\n"
            "  - TPM.tsv.gz: gzipped TSV; transcript ID in column 1\n"
            "  - X: 1-based column index in TPM file to extract TPM\n"
            "  - scale: float scaling factor to divide TPM values\n"
            "  - baseline: float baseline probability"
        )
    bed_path = sys.argv[1]
    tpm_path = sys.argv[2]
    try:
        tpm_col = int(sys.argv[3])
        if tpm_col < 1:
            die("X must be a positive 1-based column index.")
    except ValueError:
        die("X must be an integer (1-based).")
    try:
        scale = float(sys.argv[4])
        if scale == 0.0:
            die("scale must be non-zero.")
    except ValueError:
        die("scale must be a float.")
    try:
        baseline = float(sys.argv[5])
    except ValueError:
        die("baseline must be a float.")
    return bed_path, tpm_path, tpm_col, scale, baseline

def load_tpm_map(tpm_path, tpm_col_1based, scale, baseline):
    """
    Load transcript TPMs from gzipped TSV.
    Returns dict: transcript_id -> scaled_tpm (float).
    Skips lines where the selected column is missing or non-numeric.
    """
    col_idx = tpm_col_1based - 1  # convert to 0-based
    tpm_map = {}
    with gzip.open(tpm_path, "rt") as fh:
        for line_num, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("Name"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 0:
                continue
            tid = parts[0]
            # Ensure the column exists
            if col_idx >= len(parts):
                # silently skip malformed line
                continue
            val_str = parts[col_idx]
            try:
                val = float(val_str)
            except ValueError:
                # If it's a header or non-numeric, skip
                continue
            scaled = val / scale + baseline
            # If duplicates exist, keep the first occurrence; warn softly if needed.
            if tid not in tpm_map:
                tpm_map[tid] = scaled
    return tpm_map

def load_intervals_with_scaled_tpm(bed_path, tpm_map):
    """
    Read gzipped BED; return dict chrom -> list of (start, end, scaled_tpm).
    Skips BED lines whose transcript ID is not found in tpm_map or invalid coords.
    """
    per_chr = {}
    with gzip.open(bed_path, "rt") as fh:
        for line_num, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("track") or s.startswith("browser"):
                continue
            parts = s.split("\t")
            if len(parts) < 4:
                # Need at least chr, start, end, id
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                # malformed coordinates
                continue
            if end <= start:
                # zero/negative length; skip
                continue
            tid = parts[3]
            scaled_tpm = tpm_map.get(tid)
            if scaled_tpm is None:
                # No TPM for this transcript ID; skip
                continue
            per_chr.setdefault(chrom, []).append((start, end, scaled_tpm))
    # Sort intervals per chromosome by start (then end)
    for chrom in per_chr:
        per_chr[chrom].sort(key=lambda x: (x[0], x[1]))
    return per_chr

def merge_weighted(per_chr):
    """
    Merge overlapping intervals per chromosome.
    For each merged cluster, compute:
        weighted_tpm = sum(len_i * tpm_i) / sum(len_i)
    Return list of (chrom, merged_start, merged_end, weighted_tpm).
    """
    merged = []
    for chrom, intervals in per_chr.items():
        if not intervals:
            continue
        # Initialize cluster with first interval
        cur_start, cur_end, _ = intervals[0]
        # Accumulators for weighted average over the cluster
        # Note: length weighting uses each interval's full length; overlaps count for each contributing interval.
        num = (cur_end - cur_start) * intervals[0][2]
        den = (cur_end - cur_start)

        for (start, end, tpm) in intervals[1:]:
            if start <= cur_end:  # overlap or touch
                # Extend cluster
                cur_end = max(cur_end, end)
                # Accumulate weighted sums for this interval
                length = end - start
                num += length * tpm
                den += length
            else:
                # Close previous cluster
                weighted = num / den if den > 0 else 0.0
                merged.append((chrom, cur_start, cur_end, weighted))
                # Start new cluster
                cur_start, cur_end = start, end
                num = (end - start) * tpm
                den = (end - start)

        # Emit last cluster
        weighted = num / den if den > 0 else 0.0
        merged.append((chrom, cur_start, cur_end, weighted))
    return merged

def main():
    bed_path, tpm_path, tpm_col, scale, baseline = parse_args()
    # print(bed_path, tpm_path, tpm_col, scale)
    tpm_map = load_tpm_map(tpm_path, tpm_col, scale, baseline)
    # print(len(tpm_map))
    # print(tpm_map["ENSG00000237613"])
    if not tpm_map:
        die("No TPM entries were loaded—check your TSV file and column index (X).")
    per_chr = load_intervals_with_scaled_tpm(bed_path, tpm_map)
    merged = merge_weighted(per_chr)
    # Output as BED-like: chrom \t start \t end \t weighted_scaled_TPM
    out = sys.stdout

    valid_chr = {f"chr{i}" for i in range(1, 23)}
    valid_chr.add("chrY")
    valid_chr.add("chrX")
    for chrom, start, end, w in merged:
        if chrom in valid_chr:
            out.write(f"{chrom}\t{start}\t{end}\t{w:.6f}\n")

if __name__ == "__main__":
    main()


