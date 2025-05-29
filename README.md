# TFM_MarioMasero

# Membrane-Associated Proteins Induced by Irradiation – Proteomic Analysis Pipeline

This repository contains the full analytical pipeline used for the identification and characterization of membrane-associated proteins induced by irradiation in two human lung cancer cell lines (A549 and H2009), based on label-free quantitative proteomics (LFQ).

## Repository Structure

├── GeneCards/               # Gene annotation outputs (e.g., subcellular location, drugs)  
├── MetaboAnalyst/           # Results from alternative statistical analysis with MetaboAnalyst  
├── Scripts/                 # All R scripts used in the workflow  
└── README.md                # Project description and usage instructions

## Required Software and Packages

All scripts are written in R. The following packages are required (install from CRAN or Bioconductor as needed):

- `readr`, `readxl`, `dplyr`, `tidyr`, `stringr`
- `ggplot2`, `ggrepel`, `pheatmap`, `VennDiagram`, `grid`
- `limma`, `preprocessCore`
- `clusterProfiler`, `enrichplot`, `org.Hs.eg.db`, `DOSE`
- `HPAanalyze`, `BiocStyle`

## Scripts Overview

### `Scripts/Limma.R`
Performs data preprocessing, quantile normalization, and differential expression analysis using the **limma** package. Outputs include:

- DEG tables 
- Volcano plots
- Venn diagrams comparing DEPs
- Heatmaps of selected DEPs

### `Scripts/DataVisualization.R`
Performs exploratory data analysis of raw and normalized LFQ intensity distributions (histograms and abundance summaries).

### `Scripts/EnrichmentAnalysis.R`
Performs **Gene Ontology enrichment (GO)** and **Gene Set Enrichment Analysis (GSEA)** using DEPs from both cell lines. 

- `enrichGO` 
- `gseGO` analyses
- `compareCluster` comparisons for A549 vs H2009

### `Scripts/SubCellularLocation.R`
Retrieves and integrates subcellular localization data from:

- Human Protein Atlas (HPA)
- GeneCards: UniProt, Compartments, Cellular Components
- Filters for membrane-associated proteins and visualizes overlap using Venn diagrams

### `Scripts/Targeteable.R`
Identifies **targetable proteins** (approved or investigational drugs) among the membrane-associated DEPs. Outputs:

## Notes

- Some scripts use custom fixes (e.g., `hpaDownload_fixed`) due to outdated URLs in the HPAanalyze package.
- All plots are saved in the `Images/` folder (create manually if missing).
