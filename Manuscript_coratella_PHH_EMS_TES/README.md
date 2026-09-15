# Study-Specific Analysis — Coratella et al.

This folder contains the downstream analysis scripts used in Coratella et al. (Investigation of DNA damage response and viral integration profile in rAAV-transduced primary human hepatocytes)
for AAV integration site characterization in EMS-treated primary human
hepatocytes.

## Contents

| File | Description |
|------|-------------|
| `annotate_is_confidence.py` | Post-processing script to annotate integration sites with confidence flags. |
| `CIS_analysis.Rmd` | Common Integration Site (CIS) detection using a sliding-window Poisson test with data-driven minimum order thresholds. Compares AAV-only vs EMS-treated groups. |
| `../enrichment/` | Permutation-based enrichment analysis for genomic features. See [`enrichment/README.md`](../enrichment/github_upload/README.md) for details. |

## Input

Enrichment and CIS analyses take as input the `*_integration_sites.tsv` files produced by
the TES-AAV pipeline (step 5: `resolve_integration_sites.py`), optionally
post-processed with `annotate_is_confidence.py` to add confidence flags.

## Requirements

### CIS Analysis (`CIS_analysis.Rmd`)

R packages:

- `tidyverse`
- `GenomicRanges` (Bioconductor)
- `kableExtra`
- `gt23`
- `lisat`

### Enrichment Analysis

- `bedtools` (≥ 2.30)
- Python 3 with `matplotlib` and `numpy`
- Genomic feature BED files (see `enrichment/README.md` for preparation instructions)

## Usage

### CIS Analysis

The `is_dir` parameter should point to a directory containing
`*_integration_sites.tsv` files (one per sample). The sample metadata
parsing block in the Rmd should be adapted to match your file naming
convention.

### Enrichment Analysis

```bash
cd enrichment/
bash scripts/run_enrichment.sh
```

See `enrichment/README.md` for full instructions on preparing feature BED
files and configuring the analysis.
