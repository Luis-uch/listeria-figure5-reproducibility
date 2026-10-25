# Reproducibility code for Figure 5: *Listeria monocytogenes* growth and ovine gut microbiota

## Overview

This repository contains R scripts supporting the statistical analyses presented in Figure 5 of a manuscript examining associations between changes in ex vivo cultured ovine gut microbiota and the growth of *Listeria monocytogenes*.

The repository covers the genus-level association analyses for Figure 5 only. It does not reproduce primary 16S rRNA sequence processing or the other analyses in the manuscript.

## Analyses

- **Figure 5a:** Spearman rank correlations between genus-level compositional changes (ΔCLR = CLR at 72 h − CLR at 0 h) and *L. monocytogenes* growth, summarized by area under the growth curve (AUC). Correlations are calculated separately for each treatment and sampling time.
- **Figure 5b:** Per-treatment linear mixed-effects models of standardized AUC against standardized ΔCLR, with sampling time as a fixed effect and animal identity as a random intercept.

P-values are adjusted using the Benjamini–Hochberg false discovery rate (FDR) procedure within each treatment–sampling-time combination (Figure 5a) or treatment (Figure 5b). Heatmaps display genera with at least one nominal association (p < 0.05); significance symbols distinguish nominal p-values from FDR-adjusted results.

## Repository structure

```text
scripts/
  00_helpers.R                  Shared functions for importing, transforming, and pairing data
  01_prepare_data.R             Prepare microbiota abundance data
  02_growth_auc.R               Calculate growth AUC and summary tables
  03_spearman_correlations.R    Figure 5a correlations and heatmap
  04_mixed_effects_models.R     Figure 5b mixed-effects models and heatmap
data/                            Input files (availability pending)
results/                         Generated outputs (not required in the public repository)
```

The workflow expects input microbiota abundance tables, a sample mapping file, and growth measurements in `data/`. File names and required columns are defined in the scripts. The analysis uses paired microbiota observations from the beginning (0 h) and end (72 h) of ex vivo incubation.

## Software requirements

Analyses were developed using **R 4.5.2**. R package dependencies are loaded by the scripts and include packages from the tidyverse, `lme4`, `lmerTest`, `broom.mixed`, `ggplot2`, and `svglite`. Install any missing dependencies before running the workflow. Exact package versions are not pinned in this repository.

## Running the analyses

From the repository root, after placing the required input files in `data/`, run:

```bash
Rscript scripts/01_prepare_data.R
Rscript scripts/02_growth_auc.R
Rscript scripts/03_spearman_correlations.R
Rscript scripts/04_mixed_effects_models.R
```

`00_helpers.R` is sourced by the analysis scripts and does not need to be run separately. Generated tables and SVG figures are written to `results/figure_5a/` and `results/figure_5b/`.

**Note:** The analysis scripts cannot be executed end-to-end without the corresponding input data.

## Statistical considerations

Mixed-effects models retain a random intercept for animal identity across treatments. Models for the penicillin treatment yielded singular fits, with the estimated animal-level random-intercept variance at or near zero. Singularity is recorded in the model results and should be considered when interpreting these estimates. These models were retained to preserve a consistent analytical specification.

The analyses identify statistical associations, not causal effects of individual bacterial genera on *L. monocytogenes* growth.

## Code availability

The R scripts used for the statistical analyses presented in Figure 5 are provided in this repository.

## Data availability

The raw sequencing data generated in this study are available from [repository name] under accession number [accession number]. The processed genus-level abundance tables, sample metadata, and *Listeria monocytogenes* growth measurements required to reproduce the Figure 5 analyses are available from the corresponding author upon reasonable request.

<!-- TODO: Confirm whether the underlying data will be deposited in a public repository or made available from the corresponding author upon reasonable request. Update this section before publication. -->
