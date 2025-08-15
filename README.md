# QTL Mapping Scripts

This repository contains utility scripts for quantitative trait loci (QTL) analysis using [R/qtl](https://rqtl.org/).

## mqtl_pipeline.R

`mqtl_pipeline.R` performs genome scans and fits multiple QTL models for every
phenotype present in a cross file. The script is designed to be reused from the
command line with different datasets.

### Requirements
- R with the **qtl** and **optparse** packages

### Usage
```bash
Rscript mqtl_pipeline.R -i my_cross.csv -o results --permutations 1000
```

**Options**
- `-i`, `--input` : path to the input cross file (CSV rotated format)
- `-o`, `--outdir` : directory for output files (default `results`)
- `-p`, `--permutations` : number of permutations for LOD threshold (default 1000)
- `--bc-gen` : number of backcross generations (default 1)
- `--f-gen` : number of inbreeding generations (default 6)
- `--step` : step size for genotype probability calculation (default 0)
- `--na` : string used for missing values in the input file (default `NA`)

The script produces CSV and text files summarizing single- and multiple-QTL
analyses as well as diagnostic plots in the specified output directory.

## Other scripts
- `RQTL_script.r`: original long-form R/qtl workflow.
- `mqm_ghxu.r`, `mqm_with_gene_exp_cov_ghxu.r`: alternative mapping routines.
- `vcf_to_csv.pl`: Perl helper to convert VCF files into cross-compatible CSVs.

These scripts are kept for reference and may require manual editing for specific
use cases.
