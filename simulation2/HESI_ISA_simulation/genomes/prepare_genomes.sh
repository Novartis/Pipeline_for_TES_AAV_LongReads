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

#: <<'END_COMMENT'

# Download and prepare human reference genome (GRCh38):
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz GRCh38.fa.gz
gunzip GRCh38.fa.gz
(for i in {1..22}; do echo chr$i; done; echo chrX; echo chrY) | sort | uniq | seqtk subseq -l 60 GRCh38.fa - > GRCh38.fa_clean
mv GRCh38.fa_clean GRCh38.fa
samtools faidx GRCh38.fa



# Download VCF files for genome in a bottle individual HG001, HG002 and HG005:
wget -O HG001_GRCh38.vcf.gz https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -O HG002_GRCh38.vcf.gz https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -O HG005_GRCh38.vcf.gz https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz


#END_COMMENT


# SAMPLE="HG001"
for SAMPLE in HG001 HG002 HG005; do

  # Isolate SNPs (-v snps), that are either homozygous or heterozygous (-m1 -M2), passing filter (-f PASS):
  bcftools view -v snps -m1 -M2 -f PASS -Oz -o "${SAMPLE}_GRCh38.snps.pass.vcf.gz" "${SAMPLE}_GRCh38.vcf.gz"
  bcftools index "${SAMPLE}_GRCh38.snps.pass.vcf.gz"

  # Apply SNPs from allele index 1 onto the reference genome:
  bcftools consensus -f GRCh38.fa -H 1 -o "${SAMPLE}_GRCh38.hap1.fa" "${SAMPLE}_GRCh38.snps.pass.vcf.gz"
  samtools faidx "${SAMPLE}_GRCh38.hap1.fa"

  if [ -z "$(diff GRCh38.fa.fai ${SAMPLE}_GRCh38.hap1.fa.fai)" ]; then
      echo "Made reference for ${SAMPLE}"
  else
      echo "Differences detected for ${SAMPLE}! Aborting."
      exit 1
  fi

done




