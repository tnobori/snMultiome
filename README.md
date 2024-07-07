# snMultiome
This repository contains key scripts used for single-nucleus multiomics (RNA and ATAC) in [Nobori et al., 2023, bioRxiv](https://www.biorxiv.org/content/10.1101/2023.04.10.536170v1).

## **Scripts**
**[0_config_multiome.R](scripts/0_config_multiome.R)**\
This script installs necessary libraries and functions

**[1_qc_data_integration.R](scripts/1_qc_data_integration.R)**\
QC based on snRNA-seq and snATAC-seq data

**[2_clustering.R](scripts/2_clustering.R)**\
Clustering analysis based on snRNA-seq, snATAC-seq, or joint data

**[3_linkage_analysis.R](scripts/3_linkage_analysis.R)**\
Correlation analysis between mRNA expression and chromatin accessibility

**[4_motif_analysis.R](scripts/4_motif_analysis.R)**\
Motif enrichment analysis at the single-cell resolution

**[5_subclustering_pseudobulking.R](scripts/5_subclustering_pseudobulking.R)**\
Sub-clustering analysis of individual major clusters

