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

import sys, gzip
import pandas as pd
from io import StringIO

if len(sys.argv) < 2:
    sys.exit("Usage: python script.py <gtf_file>")

gtf_file = sys.argv[1]

# Read gene IDs from stdin
gene_ids = [line.strip().split('.')[0] for line in sys.stdin if line.strip()]

# Load GTF file
gtf_cols = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
with gzip.open(gtf_file, 'rt') as fh_in:
    gtf = pd.read_csv(fh_in, sep="\t", comment="#", names=gtf_cols)


# Extract gene_id
def get_gene_id(attr):
    for item in attr.split(";"):
        if "gene_id" in item:
            return item.split('"')[1].split('.')[0]
    return None

gtf["gene_id"] = gtf["attributes"].apply(get_gene_id)

# Build dictionary: gene_id -> (chrom, strand, min_start, max_end)
gene_ranges = {}
for _, row in gtf.iterrows():
    gid = row["gene_id"]
    if gid not in gene_ranges:
        gene_ranges[gid] = [row["chrom"], row["strand"], row["start"], row["end"]]
    else:
        gene_ranges[gid][2] = min(gene_ranges[gid][2], row["start"])
        gene_ranges[gid][3] = max(gene_ranges[gid][3], row["end"])

# Output BED for requested IDs
for gid in gene_ids:
    if gid in gene_ranges:
        chrom, strand, start, end = gene_ranges[gid]
        print(f"{chrom}\t{start-1}\t{end}\t{gid}\t{strand}")




