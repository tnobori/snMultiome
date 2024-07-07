##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## This script performs comparisons between single-cell RNA/ATAC-seq data and publicly available bulk RNA/ATAC-seq data
## Figures produced with this script: 

## A fully processed Seurat object can be downloaded at http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/processed_seurat_object/combined_filtered.rds

## Contact: Tatsuya Nobori (tatsuyanobori@gmail.com)
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

rm(list = ls())
source("http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/scripts/_config_multiome.R")

## Assuming you are in the directory where this and other R scripts are stored.
mainDir <- file.path(getwd(), "data_out")
subDir <- "7_bulk_omics_comparison"

# Create the nested directory structure in one step
dir.create(file.path(mainDir, subDir), recursive = TRUE, showWarnings = FALSE)

# Set the working directory to the deepest sub-directory
setwd(file.path(mainDir, subDir))

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Functions
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

optns <- theme(
  axis.text.x = element_text(margin = margin(c(20, 0, 0, 0)), size = 30, face = "bold", family = "Helvetica"),
  axis.text.y = element_text(margin = margin(c(0, 20, 0, 0)), size = 30, face = "bold", family = "Helvetica"), 
  axis.title = element_text(size = 40, face = "bold", family = "Helvetica"),
  axis.title.x = element_text(margin = margin(c(30, 0, 0, 0))),
  axis.title.y = element_text(margin = margin(c(0, 30, 0, 0))),
  axis.ticks = element_line(size = 3),
  axis.ticks.length = unit(.5, "cm"),
  axis.line  = element_line(size = 2),
  panel.background = element_blank(),
  plot.background = element_blank(),
  legend.background = element_blank(),
  legend.text = element_text(size = 15, face = "bold", family = "Helvetica", margin = margin(c(20, 0, 20, 0))),
  legend.title = element_text(size = 15, face = "bold", family = "Helvetica"),
  legend.key = element_rect(fill = NA, size = 5)
)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## snATAC vs bulk ATAC
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
data <- readRDS(paste(mainDir, "_seurat_object","combined_filtered.rds", sep = "/")) # This object is created in 1_qc_data_integration_figS1.R

library(rtracklayer)
library(GenomicRanges)


url <- "https://zenodo.org/records/10972575/files/DataSet2.Consensus_peaks.zip?download=1"
destfile <- "consensus_peaks.zip"
localPath <- file.path(mainDir, subDir, destfile)
download.file(url, destfile, mode = "wb")
bulk_ATAC_path <- file.path(mainDir, subDir)
unzip(localPath, exdir =  bulk_ATAC_path)
file.remove(destfile)


loadBED <- function(fname){
  bed_data <- read.table(file.path(bulk_ATAC_path, "Consensus_peaks",fname), 
                         header = FALSE, 
                         stringsAsFactors = FALSE, 
                         sep = "\t")
  
  bed_data$V1 <- gsub("Chr", "", bed_data$V1)
  
  # Convert to GRuanges
  bulk_data <- GRanges(seqnames = bed_data$V1,
                       ranges = IRanges(start = as.integer(bed_data$V2), end = as.integer(bed_data$V3)),
                       strand = "*")  # Assuming strand information is not necessary
  
  bulk_data
}

bedfiles <- list.files(file.path(bulk_ATAC_path, "Consensus_peaks"))

bed_data <- loadBED(bedfiles[1])
bulk_peaks_comb <- bed_data

for (i in c(2:length(bedfiles))){
  bed_data <- loadBED(bedfiles[i])
  bulk_peaks_comb <- c(bulk_peaks_comb, bed_data)
  
}

bulk_peaks_comb <- GenomicRanges::reduce(bulk_peaks_comb)
bulk_peaks_comb <-  GenomicRanges::sort(bulk_peaks_comb)

#single-cell peaks 
multiome_peaks <- data[['peaks']]@ranges
multiome_peaks <- GenomicRanges::reduce(multiome_peaks)
multiome_peaks <-  GenomicRanges::sort(multiome_peaks)

# overlap analysis 
overlapping_indices <- findOverlaps(multiome_peaks, bulk_peaks_comb)
common_peaks <- multiome_peaks[queryHits(overlapping_indices)]
unique_multiome <- GenomicRanges::setdiff(multiome_peaks, bulk_peaks_comb)


# # Count how many times each peak in multiome_peaks overlaps with any peak in bulk_peaks_comb
# overlap_counts <- countOverlaps(multiome_peaks, bulk_peaks_comb)
# 
# # See if any peak is counted more than once
# sum(overlap_counts > 1)

# Find all overlaps
all_overlaps <- findOverlaps(multiome_peaks, bulk_peaks_comb, type = "any")

# Extract overlapping peaks
common_peaks_corrected <- unique(multiome_peaks[queryHits(all_overlaps)])

# Extract unique peaks
unique_multiome_corrected <- multiome_peaks[setdiff(seq_along(multiome_peaks), queryHits(all_overlaps))]

# Check new counts
print(length(common_peaks_corrected))
print(length(unique_multiome_corrected))
print(length(common_peaks_corrected) + length(unique_multiome_corrected))


#======#======#======#======#======#======#======
#== same for bulk data======
# # Count how many times each peak in multiome_peaks overlaps with any peak in bulk_peaks_comb
# overlap_counts2 <- countOverlaps(bulk_peaks_comb, multiome_peaks)
# 
# # See if any peak is counted more than once
# sum(overlap_counts2 > 1)

# Find all overlaps
all_overlaps2 <- findOverlaps(bulk_peaks_comb, multiome_peaks, type = "any")

# Extract overlapping peaks
common_peaks_corrected2 <- unique(bulk_peaks_comb[queryHits(all_overlaps2)])

# Extract unique peaks
unique_bulk_corrected <- bulk_peaks_comb[setdiff(seq_along(bulk_peaks_comb), queryHits(all_overlaps2))]

# Check new counts
print(length(common_peaks_corrected2))
print(length(unique_bulk_corrected))
print(length(common_peaks_corrected2) + length(unique_bulk_corrected))

#======#======#======#======#======#======#======

#summerizing overlaps 
out_overlap <- data.frame(bulk=rep(0, 2), singlecell=rep(0, 2))
rownames(out_overlap) <- c("unique", "shared")

out_overlap[, 1] <- c(length(unique_bulk_corrected), length(common_peaks_corrected2))
out_overlap[, 2] <- c(length(unique_multiome_corrected), length(common_peaks_corrected))

# plot
# Transform the data frame to long format
long_data <- pivot_longer(out_overlap, cols = c(bulk, singlecell), names_to = "category", values_to = "count")

# Adding a row identification for stacking
long_data$type <- rep(c("unique", "shared"), each = 2)

# Print the transformed data
print(long_data)

# Create the stacked bar plot
q <- ggplot(long_data, aes(x = category, y = count, fill = type )) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Pastel1")  + optns

ggsave("figS1e.pdf", q)


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## snRNA vs bulk RNA
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==








