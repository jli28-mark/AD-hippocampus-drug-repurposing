# GSE5281: Differential Expression Analysis

## Overview

This directory contains the R scripts used for preprocessing, quality control and differential expression analysis of the GSE5281 hippocampal transcriptomic dataset.

The analysis includes 23 hippocampal neuronal samples, comprising 13 controls and 10 Alzheimer's disease (AD) samples.

## Scripts

**01_check_GSE5281_input.R**

- Checks sample identifiers and clinical metadata.
- Examines the expression data scale and distribution.
- Reviews available GEO preprocessing information.
- Performs quality checks without modifying the original expression data.

**02_GSE5281_DEG_age_sex.R**

- Performs log2 transformation and quantile normalization.
- Maps microarray probes to gene symbols using `hgu133plus2.db`.
- Averages retained probes corresponding to the same gene.
- Performs differential expression analysis using `limma`, adjusting for age and sex.
- Generates differential expression results, volcano plots, PCA plots and sample-level quality control outputs.

## Data Source

- Dataset: GSE5281
- Database: NCBI Gene Expression Omnibus (GEO)
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281
- Comparison: AD versus Control
- Tissue: Hippocampus
- Samples: 23 (AD = 10; Control = 13)

## Statistical Analysis

Differential expression was assessed using the limma package with the following design:

`~ Group + Age_c + Sex`

Differentially expressed genes were identified using:

- Benjamini–Hochberg-adjusted P < 0.05
- Absolute log2 fold change ≥ 0.58

Positive log2 fold changes indicate higher expression in AD samples.

## Input Requirements

The scripts require:

- `GSE5281_series_matrix.txt`: GEO probe-level expression matrix.
- `sample.xlsx` or `sample(1).xlsx`: sample metadata containing sample identifiers, diagnostic groups, age, sex and brain region.

The input-checking script additionally requires sample titles and cell-type annotations.

## Execution

The scripts contain default local file paths from the original analysis environment. Users should specify their own data directory using the `base_dir` argument when calling the analysis functions.

Run the input-checking script before the differential expression script.

## Outputs

The analysis generates gene-level differential expression tables, sample characteristics, normalization summaries, quality control plots, analysis settings and R session information.

Outputs are saved in timestamped subdirectories under `input_checks/` and `adjusted_results/`.
