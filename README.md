# MCL Thesis Analysis

Tumor microenvironment characterization of mantle cell lymphoma (MCL) using multi-omics approaches.

## Overview

This repository contains the analysis code for my doctoral dissertation, investigating MCL tumor microenvironment subtypes using TMT proteomics and RNA-seq data with BayesDeBulk deconvolution and cola consensus clustering.

## Repository Structure
```
MCL_thesis_repo/
├── MCL_thesis_analysis.Rmd    # Main analysis script
├── MCL_thesis_analysis.md     # Rendered output (viewable on GitHub)
├── data/
│   └── processed_data/        # Processed data files
├── figure-gfm/                # Generated plots
└── docs/                      # Interactive reports
```

## Interactive Reports

The cola report could not be displayed fully in my docs/ folder. Here is the link to the full report:
- [Cola Consensus Clustering Report](file:///Users/heona/git-repos/MCL_thesis_repo/docs/cola_bayesdb_log/cola_report.html)

## Methods

- **Proteomics normalization**: Sample loading normalization, internal reference scaling, HarmonizR batch correction
- **Deconvolution**: BayesDeBulk for cell type estimation
- **Clustering**: Cola consensus clustering with multiple methods (ATC, SD, MAD × hclust, kmeans, skmeans)
- **Visualization**: PCA, UMAP, heatmaps
- **GSEA (Gene set enrichment analysis)**: heatmaps, volcano plots
- **GSVA (Gene set variation analysis)**: heatmaps

## Data Availability

Raw data is not included due to size and privacy constraints.

## Author

Heona Lee
