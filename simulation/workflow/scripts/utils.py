#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

import sys, os, gzip, time, re, json, copy
from io import StringIO
from collections import Counter, defaultdict
import numpy as np
import pandas as pd



def weighted_hamming_dist(seq1, seq2, weights, strand='+'):
    if weights is None:
        weights = list(np.ones(len(seq1)))
    if strand == '-':
        weights = weights[::-1]
    return(sum(w for s1, s2, w in zip(seq1, seq2, weights) if s1!=s2))


def find_target(seq, target_seq, target_dist_w=None):
    assert(len(seq) >= len(target_seq))
    if not target_dist_w is None:
        assert(len(target_dist_w) == len(target_seq))

    tlen = len(target_seq)
    hdist = list()
    for i in range(0, len(seq)+1-tlen):
        d = weighted_hamming_dist(seq[i:(i+tlen)], target_seq, target_dist_w)
        dr = weighted_hamming_dist(seq[i:(i+tlen)], str(Seq(target_seq).reverse_complement()), target_dist_w, strand='-')
        hdist.append((i, '+', d))
        hdist.append((i, '-', dr))
    target_index, target_strand, target_dist = sorted(hdist, key=lambda x: x[2])[0]
    return(target_index, target_strand, target_dist)


def cigar2tlen(cigar):
    cigar_op = re.split(r"([MDI])", cigar)
    tlen = 0
    i = 0
    while i < (len(cigar_op)-1):
        if cigar_op[i+1] in ['M', 'D']:
            tlen += int(cigar_op[i])
        i += 1
    return(tlen)


def cigar2indel_len(cigar):
    cigar_op = re.split(r"([MDI])", cigar)
    tlen = 0
    i = 0
    while i < (len(cigar_op)-1):
        if cigar_op[i+1] in ['I', 'D']:
            tlen += int(cigar_op[i])
        i += 1
    return(tlen)


def most_freq_len(lst):
    lc = defaultdict(int)
    for s in lst:
        lc[len(s)]+=1
    lc = dict(lc)
    fl = sorted(lc.items(), key=lambda x: x[1], reverse=True)[0][0]
    return(fl)


def find_consensus(sequences, min_pos_frac):
    # Transpose the list of sequences to get columns of nucleotides
    transposed = list(zip(*sequences))

    consensus = list()
    for column in transposed:
        # Count the occurrences of each nucleotide in the column:
        counts = Counter(column)
        # Find the most common nucleotide:
        consensus_nucleotide = counts.most_common(1)[0][0]
        consensus_count = counts.most_common(1)[0][1]
        if (consensus_count / len(column)) >= min_pos_frac:
            consensus.append(consensus_nucleotide)
        else:
            consensus.append('N')

    return(''.join(consensus))










