#!/bin/bash


# Reference Files:
# 
# Genome and annotation files:
#    - Download the genome and annotation files used in this study from the following URL:
#      http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/arabidopsis_genome/
#
# Cell Ranger reference file:
#    - Download the Cell Ranger reference file from the following URL:
#      http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/cellranger_reference/


# Instructions for Creating a Sample Sheet
# 
# 1. Download FASTQ files:
#    - Use a command like `wget` in the terminal to download the FASTQ files from the following URL:
#      http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/fastq_files/
#
# 2. Create a sample sheet:
#    - Follow the guideline provided by 10x Genomics for specifying input FASTQ files:
#      https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/inputs/specifying-input-fastq-count
#
# 3. The sample sheet should be a CSV file with the following format:
#
#     fastq,sample,library_type
#     path_to_RNA_fastq,sample_name,Gene Expression
#     path_to_ATAC_fastq,sample_name,Chromatin Accessibility
#
# Replace `path_to_RNA_fastq` and `path_to_ATAC_fastq` with the actual paths to your RNA and ATAC FASTQ files, respectively.
# Replace `sample_name` with the name of your sample.


Path="/path/to/your/directory"

# This is your output directory
cd ${Path}/count

export PATH=/path/to/cellranger-arc/cellranger-arc-2.0.0:$PATH

cellranger-arc count --id=sample_name \
                 --reference=/path/to/reference_files/TAIR10_fixed_211122_CHremoved \
                 --libraries=samplesheet.csv #sample sheet specifying the paths to RNA and ATAC fastqs
