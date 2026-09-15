# Enrichment Analysis for AAV Integration Sites

This module tests whether AAV integration sites are enriched near
liver-expressed genes using a permutation-based approach with `bedtools`.

## Folder Structure

```
enrichment/
├── scripts/
│   ├── run_enrichment.sh              # Main enrichment pipeline
│   ├── preprocess_integration_site.py # Convert integration_sites.tsv → BED
├── reference/
│   ├── hg38.chrom.sizes              # Chromosome sizes for hg38
│   └── blacklist.v2.bed              # ENCODE blacklist v2
├── features/                          # Feature BED files (user-provided)
│   └── README.md                      # Describes required feature files
└── output/                            # Results (created by pipeline)
```

## Quick Start

```bash
# 1. Place feature BED files in features/ (see "Preparing Feature BED Files" below)

# 2. Convert pipeline output to BED format
python scripts/preprocess_integration_site.py sample_integration_sites.tsv input_sites.bed

# 3. Run enrichment analysis
bash scripts/run_enrichment.sh input_sites.bed output/enrichment_results.tsv


## Input Requirements

### Integration Sites TSV

The input to `preprocess_integration_site.py` is the `*_integration_sites.tsv`
file produced by the main TES pipeline (`resolve_integration_sites` rule).
Required columns:

| Column | Description |
|--------|-------------|
| `site_id` | Unique identifier for the integration site |
| `chr_bp_upstream` | Chromosome of upstream breakpoint |
| `pos_bp_upstream` | Position of upstream breakpoint |
| `chr_bp_downstream` | Chromosome of downstream breakpoint |
| `pos_bp_downstream` | Position of downstream breakpoint |
| `count` | Number of supporting reads |
| `orientation` | Orientation (fwd/rev) |
| `is_translocation` | Boolean: whether site spans different chromosomes |

Optional column: `confidence` (used with `--exclude-confidence` flag).

### Integration Sites BED (input to run_enrichment.sh)

If you already have a BED file of integration sites (not from this pipeline),
it must be **BED6 format** with standard chromosomes:

```
chr1    12345    12346    site_1    5    +
chr2    67890    67891    site_2    3    -
```

| Column | Description |
|--------|-------------|
| 1. chrom | Chromosome (chr1-22, chrX, chrY) |
| 2. start | 0-based start position |
| 3. end | End position (start + 1 for point intervals) |
| 4. name | Site identifier |
| 5. score | Supporting read count |
| 6. strand | Strand (+, -, or .) |

## Preparing Feature BED Files

The enrichment script expects BED files of genomic features in `features/`.
These represent liver-expressed gene regions to test for overlap with
integration sites.

### Required feature BED files

Place the following files in `features/`:

| File | Description |
|------|-------------|
| `liver_expressed_genes_high.bed` | Gene bodies of highly expressed liver genes |
| `liver_expressed_genes_high_tss_10kb.bed` | TSS ±10kb windows around highly expressed genes |
| `liver_expressed_genes_high_tss_50kb.bed` | TSS ±50kb windows|
| `liver_expressed_genes_high_tss_100kb.bed` | TSS ±100kb windows |
| `liver_expressed_genes_mid_high.bed` | Gene bodies of mid+high expressed liver genes |
| `liver_expressed_genes_mid_high_tss_10kb.bed` | TSS ±10kb windows around mid+high expressed genes |
| `liver_expressed_genes_mid_high_tss_50kb.bed` | TSS ±50kb windows  |
| `liver_expressed_genes_mid_high_tss_100kb.bed` | TSS ±100kb windows |

### How to generate these files

1. **Identify liver-expressed genes** from an expression dataset such as:
   - [GTEx v8 median TPM](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression)
   - Your own RNA-seq data from the target tissue (e.g., hepatocytes)

2. **Classify genes by expression level** into tiers:
   - **high**: highly expressed
   - **mid+high**: moderately to highly expressed

3. **Map gene symbols to genomic coordinates** using a GENCODE annotation:
   - Download [gencode.v49.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz)
   - Extract gene-level entries (feature type = "gene")
   - Convert GTF 1-based coordinates to BED 0-based

4. **Create gene body BED files** — one row per gene:
   ```
   chr1    11869    14409    DDX11L1
   chr1    14404    29570    WASH7P
   ```

5. **Create TSS window BED files** — for each gene, take the transcription
   start site (start coordinate) and extend ±10kb/±50kb/±100kb.
   Clamp to chromosome boundaries using `hg38.chrom.sizes`, then merge
   overlapping intervals (e.g., with `bedtools merge`).

### Feature BED format

All files must be **tab-separated BED** with minimum 3 columns:

```
chrom    start    end    [name]
```

| Column | Description |
|--------|-------------|
| 1. chrom | Chromosome (e.g., chr1) |
| 2. start | 0-based start position |
| 3. end | End position |
| 4. name | (optional) Gene name or feature ID |

All coordinates must be **hg38 (GRCh38)**.

## Methodology

- **Permutation test:** Integration sites are shuffled randomly across the
  genome 1000 times using `bedtools shuffle`
- **Blacklist exclusion:** ENCODE blacklist regions are excluded from shuffling
- **No stacking:** `-noOverlapping` prevents random sites from piling up
- **Overlap counting:** `bedtools intersect -u` counts each site at most once
- **P-value:** Empirical p = (M + 1) / (N + 1) where M = random iterations
  with overlaps >= observed, N = total iterations
- **Fold enrichment:** observed overlaps / mean random overlaps

## Requirements

- Python 3 with: `pandas`, `numpy`, `matplotlib`
- [bedtools]