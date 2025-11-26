#!/bin/bash
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

# Download median TPM from GTEx:
wget https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz

# Isolate the Liver and Muscle data:
gunzip -dc GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz | python clean_GTEx_data.py | gzip > GTEx_TPM_liver_muscle.tsv.gz


# Download transcript exon annotations from GENCODE (GTF format):
 wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.basic.annotation.gtf.gz

# Isolate the start and end for each Ensembl ids:
gunzip -dc GTEx_TPM_liver_muscle.tsv.gz | cut -f1 | grep "^ENSG" | python isolate_transcript_ranges.py gencode.v49.primary_assembly.basic.annotation.gtf.gz | gzip > transcript_ranges.bed.gz


# Prepare BED file with integration site sampling probabilities based on TPM values:
python make_bias_bed.py transcript_ranges.bed.gz GTEx_TPM_liver_muscle.tsv.gz 2 10 1 > liver_integration_bias.bed
python make_bias_bed.py transcript_ranges.bed.gz GTEx_TPM_liver_muscle.tsv.gz 3 10 1 > muscle_integration_bias.bed



