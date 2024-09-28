#!/bin/bash
Path="/path/to/your/directory"

# This is your output directory
cd ${Path}/count

export PATH=/path/to/cellranger-arc/cellranger-arc-2.0.0:$PATH

cellranger-arc count --id=sample_name \
                 --reference=/path/to/reference_files/TAIR10_fixed_211122_CHremoved \
                 --libraries=samplesheet.csv #sample sheet specifying the paths to RNA and ATAC fastqs
