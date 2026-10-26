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
data/                            Input files included with the repository
results/                         Generated outputs (not required in the public repository)
```

## Input data

The input files required to reproduce the Figure 5 analyses are included in `data/`:

- `mapping_file.tsv`: tab-separated sample metadata. It links each microbiota `SampleID` to treatment, in vivo sampling time, ex vivo incubation time (`t0` or `t72`), and experimental animal (`Sheep`).
- `growth.csv`: semicolon-separated *L. monocytogenes* growth measurements. It contains the sample identifier, treatment, sampling time, and logCFU values measured after 0, 24, 48, and 72 hours of ex vivo incubation.
- `Microbial Metagenomics - t0_abundance.csv`: taxon-by-sample count table at the beginning of ex vivo incubation. Taxonomic annotations are stored in the first column and sample identifiers in the column headers.
- `Microbial Metagenomics - t72_abundance.csv`: taxon-by-sample count table after 72 hours of ex vivo incubation.
- `Microbial Metagenomics - t0_rel_abundance.csv`: relative-abundance counterpart of the `t0` count table. Values are expressed as fractions rather than percentages.
- `Microbial Metagenomics - t72_rel_abundance.csv`: relative-abundance counterpart of the `t72` count table. Values are expressed as fractions rather than percentages.

The workflow aggregates the taxonomic tables at genus level and pairs observations from the beginning (`t0`) and end (`t72`) of ex vivo incubation. The files contain experimental sample identifiers only.

## Software requirements

Analyses were developed using **R 4.5.2**. R package dependencies are loaded by the scripts and include packages from the tidyverse, `lme4`, `lmerTest`, `broom.mixed`, `ggplot2`, and `svglite`. Install any missing dependencies before running the workflow. Exact package versions are not pinned in this repository.

## Running the analyses

From the repository root, run:

```bash
Rscript scripts/01_prepare_data.R
Rscript scripts/02_growth_auc.R
Rscript scripts/03_spearman_correlations.R
Rscript scripts/04_mixed_effects_models.R
```

`00_helpers.R` is sourced by the analysis scripts and does not need to be run separately. Generated tables and SVG figures are written to `results/figure_5a/` and `results/figure_5b/`.

The required input files are provided in `data/`. Generated outputs are written to `results/`, `results/figure_5a/`, and `results/figure_5b/`.

## Statistical considerations

Mixed-effects models retain a random intercept for animal identity across treatments. Models for the penicillin treatment yielded singular fits, with the estimated animal-level random-intercept variance at or near zero. Singularity is recorded in the model results and should be considered when interpreting these estimates. These models were retained to preserve a consistent analytical specification.

The analyses identify statistical associations, not causal effects of individual bacterial genera on *L. monocytogenes* growth.

## Code availability

The R scripts used for the statistical analyses presented in Figure 5 are provided in this repository.

## Data availability

The raw sequencing data generated in this study are available from [repository name] under accession number [accession number]. The processed abundance tables, sample metadata, and *Listeria monocytogenes* growth measurements required to reproduce the Figure 5 analyses are provided in the `data/` directory of this repository.

<!-- TODO: Replace the repository name and accession number before publication. -->
