# snMultiome
This repository contains key scripts used for single-nucleus multiomics (RNA and ATAC) in [Nobori et al., 2023, bioRxiv](https://www.biorxiv.org/content/10.1101/2023.04.10.536170v1).

## **Scripts**
**[snMultiome_preprocessing_cellranger.sh](scripts/snMultiome_preprocessing_cellranger.sh)**\
cellranger-arc for processing snMultiome raw data

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

**[6_gt3a_bulkRNAseq_fig6.R](scripts/6_gt3a_bulkRNAseq_fig6.R)**\
Bulk RNA-seq analysis of a GT-3A overexpression line

**[7_comparison_with_bulk_omics_figS1_figS3.R](scripts/7_comparison_with_bulk_omics_figS1_figS3.R)**\
Comparisons between snRNA-seq/snATAC-seq and bulk RNA-seq/ATACseq

**[8_snRNA-seq_of_gt3aKO_figS8.R](scripts/8_snRNA-seq_of_gt3aKO_figS8.R)**\
snRNA-seq analysis of a GT-3A knowckout mutant 