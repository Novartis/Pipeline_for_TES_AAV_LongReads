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
import pysam

aav_ids_file = sys.argv[1]
genome_ids_file = sys.argv[2]
input_bam = sys.argv[3]
output_bam = sys.argv[4]

# Load read IDs
with open(aav_ids_file) as f:
    aav_ids = set(line.strip() for line in f)

with open(genome_ids_file) as f:
    genome_ids = set(line.strip() for line in f)

shared_ids = aav_ids & genome_ids

# Filter BAM
with pysam.AlignmentFile(input_bam, "rb") as bam_in, \
     pysam.AlignmentFile(output_bam, "wb", template=bam_in) as bam_out:
    for read in bam_in.fetch(until_eof=True):
        if read.query_name in shared_ids:
            bam_out.write(read)

# Index
pysam.index(output_bam)
