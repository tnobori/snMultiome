##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## This script combines cellranger-arc outputs from individual samples and performs QC and cell filtering
## Figures produced with this script: FigS1

## A fully processed Seurat object can be downloaded at http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/processed_seurat_object/combined_filtered.rds

## Contact: Tatsuya Nobori (tatsuyanobori@gmail.com)
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

rm(list = ls())
source("http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/scripts/_config_multiome.R")

options(timeout = 6000) # Set timeout to 100 minutes


## Assuming you are in the directory where this and other R scripts are stored.
mainDir <- file.path(getwd(), "data_out")
dir.create(file.path(mainDir), showWarnings = FALSE)

subDir <- "1_qc_data_integration"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

setwd(file.path(mainDir, subDir))
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Functions
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
create.seurat <- function(countpath, fragpath, annotation){

  counts <- Read10X_h5(countpath)

  # create a Seurat object containing the RNA adata
  data <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA"
  )

  # create ATAC assay and add it to the object
  data[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )

  data

}

# You need to install MACS2. This can be done using pip or conda, or by building the package from source. 
# Specify the path to macs2 in the following function.
loadBEDandGenomeData.mod <- function(bed, ann, sizes, attribute="Parent", verbose=T, is.fragment=T, genomesize=0.8e8,loadbed=T, macs=F,downsample=0.01,
                                    shift= -50,
                                    extsize=100,
                                    output="bulk_peaks",
                                    tempdir="./macs2_temp",
                                    fdr=0.05,
                                    macs2.path = "/opt/anaconda3/bin/macs2"){
  
  mac2temp <- tempdir
 
  if(loadbed){
    if(verbose){message(" - loading data (this may take a while for big BED files) ...")}
    if(grepl(".gz$", bed)){
      bedData <- read.table(gzfile(as.character(bed)))
    }else{
      bedData <- read.table(as.character(bed))
    }
  }
  # toc()
  #loading: 195.606 sec elapsed
  if(grepl(".gtf", ann)){
    anntype <- "gtf"
  }else{
    anntype <- "gff3"
  }

  # tic("converting")
  if(verbose){message(" - converting fragment file to single-bp resolution Tn5 insertions sites ...")}
  start.coordinates <- data.frame(V1=bedData$V1, V2=(bedData$V2), V3=(bedData$V2+1), V4=bedData$V4, V5="+")
  start.coordinates <- start.coordinates[!duplicated(start.coordinates),]
  end.coordinates <- data.frame(V1=bedData$V1, V2=bedData$V2, V3=(bedData$V3-1), V4=bedData$V4, V5="-")
  end.coordinates <- end.coordinates[!duplicated(end.coordinates),]
  all.coordinates <- rbind(start.coordinates, end.coordinates)
  all.coordinates2 <- all.coordinates[sample(1:nrow(all.coordinates), nrow(all.coordinates)*downsample) ,]

    # load Gff
  gff <- suppressWarnings(suppressMessages(makeTxDbFromGFF(as.character(ann), format=anntype, dbxrefTag=attribute)))
  chrom <- read.delim(sizes, header = F)

  if(macs){
    # verbose
    if(verbose){message(" - running MACS2 on bulk BED file ...")}
    #create temp dif
    mac2temp <- tempdir
    if(file.exists(mac2temp)){
      unlink(mac2temp)
      dir.create(mac2temp)
    }else{
      dir.create(mac2temp)
    }
    
    # build command
    cmdline <- paste0(macs2.path, " callpeak -t ", bed, " -f BED -g ", genomesize, " --keep-dup all -n ", output,
                      " --nomodel --shift ",shift, " --extsize ", extsize, " --outdir ",tempdir, " --qvalue ", fdr)
    # run macs2
    suppressMessages(system(cmdline))
  }
 
  # load peaks
  peaks <- read.table(paste0(mac2temp,"/",output,"_peaks.narrowPeak"))

  #' buildMetaData_OG
  skip.PtMt = TRUE
  
  bed2 <- all.coordinates2
  acr <- peaks
  
  #convert fragment bed file to Granges
  bed.gr <- GRanges(seqnames=as.character(bed2$V1),
                    ranges=IRanges(start=as.numeric(bed2$V2),
                                   end=as.numeric(bed2$V3)),
                    strand=as.character(bed2$V5),
                    names=as.character(bed2$V4))

  # convert ACRs to Granges
  acr.gr <- GRanges(seqnames=as.character(acr$V1),
                    ranges=IRanges(start=as.numeric(acr$V2),
                                   end=as.numeric(acr$V3)))

  # get up and down stream of TSS
  tss <- promoters(gff, upstream=2000, downstream=1000)
  tss <- tss[!duplicated(tss),]

  # get reads that overlap tss
  if(verbose){message(" - counting Tn5 sites at TSSs per barcode ...")}
  tss.reads <- suppressWarnings(subsetByOverlaps(bed.gr, tss, ignore.strand=T))
  tss.reads <- as.data.frame(tss.reads)
  tss.counts <- table(tss.reads$names)

  # get reads overlapping ACRs
  if(verbose){message(" - counting Tn5 sites within ACRs per barcode ...")}
  acr.reads <- suppressWarnings(subsetByOverlaps(bed.gr, acr.gr, ignore.strand=T))
  acr.reads <- as.data.frame(acr.reads)
  acr.counts <- table(acr.reads$names)

  if(verbose){message(" - finalizing meta data creation ...")}

  # merge
  ids <- unique(as.character(bed2$V4))
  counts <- table(bed2$V4)
  counts <- counts[ids]
  tss.counts <- tss.counts[ids]
  acr.counts <- acr.counts[ids]
  if(!skip.PtMt){
    org.counts <- org_val[ids]
  }
  counts[is.na(counts)] <- 0
  tss.counts[is.na(tss.counts)] <- 0
  acr.counts[is.na(acr.counts)] <- 0
  if(!skip.PtMt){
    org.counts[is.na(org.counts)] <- 0
  }else{
    org.counts <- rep(NA, length(ids))
    names(org.counts) <- ids
  }


  df <- data.frame(cellID=ids,
                   total=as.numeric(counts),
                   tss=as.numeric(tss.counts),
                   acrs=as.numeric(acr.counts),
                   row.names=ids)

  return(df)
}


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Preparing an annotation file; THIS NEEDS TO BE RUN ONLY ONCE
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Download GTF file from server
gffUrl <- "http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/arabidopsis_genome/Arabidopsis_thaliana.TAIR10.45_ChRemoved.gtf"
gffPath <- "./Arabidopsis_thaliana.TAIR10.45_ChRemoved.gtf"
if (!file.exists(gffPath)) {
  download.file(gffUrl, destfile = gffPath)
}

# In the GTF file, there are cases where a row has two gene_name columns. 
# An example is LRR XI-23 and several others. 
# Actual gene names seem to be used for cellranger output, but gene IDs are loaded as annotation.
# The following code fixes this issue.
gtf_lines <- readLines(gffPath)

# First, AT4G04890 needs to be treated separately because there are two PDF2, while there are also PDF2.1-PDF2.6, which
# made one of PDF2s PDF2.7. 
modify_gene_name <- function(line) {
  if (grepl('gene_id "AT4G04890";', line)) {
    # Pattern to find the gene_name attribute
    gene_name_pattern <- 'gene_name "[^"]+"'
    # Replacement string with the new gene_name
    replacement <- 'gene_name "PDF2.7"'
    # Replace the gene_name in the line
    line <- sub(gene_name_pattern, replacement, line)
  }
  return(line)
}
gtf_lines <- sapply(gtf_lines, modify_gene_name, USE.NAMES = FALSE)


correct_line <- function(line) {
  # Find all matches for gene_name occurrences
  matches <- gregexpr('gene_name "[^"]+"', line, perl = TRUE)[[1]]
  match_lengths <- attr(matches, "match.length")
  
  # Check if there are multiple gene_name occurrences
  if(length(matches) > 1 && matches[1] != -1) {
    # Extract the first gene_name, correctly handling the full match
    first_gene_name_match <- substr(line, matches[1], matches[1] + match_lengths[1] - 1)
    
    # Replace all subsequent gene_name occurrences with the first one's full match
    for (i in 2:length(matches)) {
      line <- substr_replace(line, first_gene_name_match, 
                             start = matches[i], stop = matches[i] + match_lengths[i] - 1)
    }
  }
  return(line)
}

# Helper function to replace substrings within a string
substr_replace <- function(string, replacement, start, stop) {
  prefix <- substr(string, 1, start - 1)
  suffix <- substr(string, stop + 1, nchar(string))
  return(paste0(prefix, replacement, suffix))
}

# Apply the correction to each line
corrected_gtf_lines <- sapply(gtf_lines, correct_line, USE.NAMES = FALSE)

# Path for the corrected GTF file
corrected_gtf_path <- "Corrected_Arabidopsis_thaliana.TAIR10.45_ChRemoved.gtf"

# Write the corrected lines to a new GTF file
writeLines(corrected_gtf_lines, corrected_gtf_path)

# Inport corrected annotation with rtracklayer
annotation <- rtracklayer::import(corrected_gtf_path)

# colnames(annotation@elementMetadata)[colnames(annotation@elementMetadata)=="Name"] <- "gene_name" #this change is needed for AnnotationPlot to work
# colnames(annotation@elementMetadata)[colnames(annotation@elementMetadata)=="biotype"] <- "gene_biotype" #this change is needed for AnnotationPlot to work
levels(annotation$type)[levels(annotation$type)=="gene"]="body"
annotation$gene_name[grepl("exon", annotation$gene_name)] <- gsub("\\.[0-9].exon[0-9]{1,}", "", annotation$gene_name[grepl("exon", annotation$gene_name)]) #remove exon number
annotation$tx_id <- annotation$transcript_id # tx_id column is needed later for Signac analysis 
annotation$gene_name <- gsub("_", "-", annotation$gene_name) # Seurat does not accept "_" in gene names, and it replaces them with "-". So, the annotation should match these changes.

## The original GTF contains duplicated gene_name. Seurat deals with this by adding ".1" ".2" to the duplicated gene names.
## The following modifications are to match annotation gene_name to cellranger outputs.
## Modify duplicated gene names by adding ".1", ".2" ... at the end in an order same as gene IDs. 
## The first gene will remain unchanged.
## Function to correct a single line

# Define the function to make character vector names unique
makeUniqueNames <- function(charVector) {
  uniqueNames <- vector("character", length(charVector)) # Initialize an empty character vector to store unique names
  countedValues <- table(charVector) # Count occurrences of each value
  
  for (value in names(countedValues)) {
    occurrences <- countedValues[[value]] # Number of occurrences of the current value
    if (occurrences > 1) {
      # Generate unique names for all occurrences except the first one
      uniqueNames[charVector == value] <- c(value, paste0(value, ".", 1:(occurrences - 1)))
    } else {
      uniqueNames[charVector == value] <- value
    }
  }
  
  return(uniqueNames)
}

genenames <-annotation$gene_name %>% unique()
genebody <- annotation$gene_name[annotation$type=="body"]
dupgenes <- table(genebody)[table(genebody)>1] %>% names() # These are duplicated gene names to be fixed


# Step 1: Filter 'body' type annotations for only those genes in 'dupgenes'
body_annotations <- annotation[annotation$type == "body" & annotation$gene_name %in% dupgenes, ]

# Apply makeUniqueNames to these 'body' type annotations

total_genes <- length(dupgenes)
pb <- txtProgressBar(min = 0, max = total_genes, style = 3)

for (i in seq_along(dupgenes)) {
  genes <- dupgenes[i]
  body_annotations$gene_name[body_annotations$gene_name %in% genes ] <- 
    makeUniqueNames(body_annotations$gene_name[body_annotations$gene_name %in% genes ])
  
  # Update progress bar
  setTxtProgressBar(pb, i)
}
close(pb)

# Create a mapping from gene_id to the new unique gene_name for 'body' type annotations
gene_id_to_unique_name <- setNames( body_annotations$gene_name, body_annotations$gene_id)

# Step 2: Update gene_name for all annotations using the mapping, but only for genes in 'dupgenes'
# Filter the main annotation dataframe to include only entries with gene_id that need updates
annotations_to_update <- annotation[annotation$gene_id %in% names(gene_id_to_unique_name), ]

# Iterate through these filtered annotations
total_genes <- length(annotations_to_update)
pb <- txtProgressBar(min = 0, max = total_genes, style = 3)
for (i in c(1:length(annotations_to_update))) {
  gene_id <- annotations_to_update$gene_id[i]
  # Update the gene_name in the original 'annotation' DataFrame
  annotation$gene_name[annotation$gene_id == gene_id] <- gene_id_to_unique_name[[gene_id]]
  setTxtProgressBar(pb, i)
}
close(pb)

# Here is another exceptional case that needs to be addressed separetely
annotation$gene_name[annotation$gene_name %in% "LHCB1.1"  &  annotation$type=="body"] <- makeUniqueNames(annotation$gene_name[annotation$gene_name %in% "LHCB1.1"  &  annotation$type=="body"])

export(annotation, con = "Corrected_Arabidopsis_thaliana.TAIR10.45_ChRemoved.gtf", format = "gtf")

# load corrected annotation
annotation <- rtracklayer::import("Corrected_Arabidopsis_thaliana.TAIR10.45_ChRemoved.gtf")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Download required files to your local device.
## Alternatively, you may download these files with wget()
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

gffUrl <- "http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/arabidopsis_genome/Arabidopsis_thaliana.TAIR10.45_ChRemoved.gtf"
gffPath <- "./Arabidopsis_thaliana.TAIR10.45_ChRemoved.gtf"
if (!file.exists(gffPath)) {
  download.file(gffUrl, destfile = gffPath)
}

chromUrl <- "http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/arabidopsis_genome/at_chromosom_file.txt"
chromPath <- "./at_chromosom_file.txt"
if (!file.exists(chromPath)) {
  download.file(chromUrl, destfile = chromPath)
}

sampleListUrl <- "http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/misc/sample_list.txt"
sampleListPath <- "./sample_list.txt"
if (!file.exists(sampleListPath)) {
  download.file(sampleListUrl, destfile = sampleListPath)
}

samples <- read.delim(sampleListPath, header = F)

# Downloading raw 10x outputs
for (i in c(1:nrow(samples))) {
  sampleDir <- file.path(mainDir, samples[i, 1])
  dir.create(sampleDir, showWarnings = FALSE)
  url_list <- paste0("http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/10x_outputs/", samples[i, 1], "/", 
                     c("atac_fragments.tsv.gz", "atac_fragments.tsv.gz.tbi", "filtered_feature_bc_matrix.h5"))
  
  for (url in url_list) {
    localPath <- file.path(sampleDir, basename(url))
    if (!file.exists(localPath)) {
      download.file(url, destfile = localPath)
    }
  }
}


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## creating Seurat object for each sample
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
# create a list of seurat objects
list_seurat <- list()
for (i in c(1:nrow(samples))){
  try(list_seurat[[i]] <- create.seurat(countpath = file.path(samples[i, 1], "filtered_feature_bc_matrix.h5"), 
                                        fragpath = file.path(samples[i, 1], "atac_fragments.tsv.gz"), 
                                        annotation = annotation ))
}
names(list_seurat) <- samples[, 1]

# rename cells using object names as prefix
for (i in names(list_seurat)) {
  try(
    list_seurat[[i]] <- RenameCells(list_seurat[[i]],
                                         add.cell.id = i)
  )
}

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
# doublet removal; this should be done before data integration
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# function for running doubletFinder
doubletFinder_run <- function(data){
  data <- SCTransform(data, verbose = FALSE)
  data <- RunPCA(data, verbose = FALSE)
  data <- FindNeighbors(data, dims = 1:20, verbose = FALSE)
  data <- FindClusters(data, resolution = 1.5, verbose = FALSE, algorithm = 3)
  data <- RunUMAP(data, dims = 1:20, verbose = FALSE)

  ## pK Identification (no ground-truth)
  sweep.res.data <- paramSweep_v3(data, PCs = 1:20, sct = TRUE)
  sweep.data <- summarizeSweep(sweep.res.data, GT = FALSE)
  bcmvn_data <- find.pK(sweep.data)
  pK <- bcmvn_data$pK %>% levels() %>% as.numeric() %>% mean()
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  clusters <- data@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(clusters)       
  nExp_poi <- round(0.075*nrow(data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  data <- doubletFinder_v3(data, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  data
}

for (i in c(1:nrow(samples))){
  tic(msg = paste0("sample ", i, " of ", nrow(samples)))
  list_seurat[[i]] <- doubletFinder_run(list_seurat[[i]])
  toc()
}

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
#combining seurat objects
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

data <- Merge_Seurat_List(
  list_seurat,
  add.cell.ids = NULL,
  merge.data = TRUE,
  project = "PathogenMultiome"
)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
#adding other metadata
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
#adding other metadata
list.sample <- samples[,1] %>% as.list()
samplenames <- lapply(list.sample, function(x){rep(x, times = length(grep(x, colnames(data))))})
data$sample <- unlist(samplenames)
data$sample <- gsub("col_", "", data$sample)

#time points
sample.time <- gsub("_rep*.", "", data$sample)
sample.time <- gsub(".*_", "",sample.time)
data$time <- sample.time

#MT score
mt.gene.idx <- c(grep("ORF153A", rownames(data)):grep("ORF204", rownames(data))) #these are MT genes
mt.gene.idx <- c(mt.gene.idx, grep("ATMG", rownames(data)))
mt.count <- data@assays$RNA@counts[mt.gene.idx, ] %>% colSums()
percent.mt <- (mt.count / data$nCount_RNA) * 100
data$percent.mt <- percent.mt

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
#save data
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
# Adding annotation file to the seurat object
# In the future, you can load annotation by: annotation <- Annotation(data[['ATAC']])
Annotation(data@assays$ATAC) <- annotation 

out.dir.seurat <- file.path(mainDir, "_seurat_object")
dir.create(out.dir.seurat, showWarnings = FALSE)
saveRDS(data, paste0(out.dir.seurat, "/_data_combined_unfiltered.rds"))

#loading unfiltered data
data <- readRDS(paste(mainDir, "_seurat_object","_data_combined_unfiltered.rds", sep = "/"))

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## ATAC QC analysis; this takes time to finish.
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
out.dir <- file.path(mainDir, subDir, "ATAC_QC_files")
dir.create(out.dir, showWarnings = FALSE)

for (i in c(1:nrow(samples))){
  
  out <- loadBEDandGenomeData.mod(bed = file.path(samples[i, 1], "atac_fragments.tsv.gz"), ann =gffpath, sizes = chrompath, macs = T)
  write.table(out, file = paste0(out.dir ,"/acr_tss_QC_", samples[i, 1],".txt"), row.names=T, col.names=NA, sep="\t", quote=F)
  
}

##plotting QC results

optns <- theme(
  axis.text.x = element_text(margin = margin(c(5, 0, 0, 0)), size = 20, face = "bold", family = "Helvetica",),
  axis.text.y = element_text(margin = margin(c(0, 5, 0, 0)) ,size = 20, face = "bold", family = "Helvetica"), 
  axis.title = element_text(size = 20, face = "bold", family = "Helvetica"),
  axis.ticks.y = element_line(size = 1),
  axis.ticks.x = element_line(size = 1),
  axis.ticks.length = unit(.3, "cm"),
  axis.line  = element_line(size = 1),
  panel.background = element_blank(),
  legend.text = element_text(size = 15, face = "bold", family = "Helvetica", margin = margin(c(20, 0, 20, 0))),
  legend.title = element_text(size = 15, face = "bold", family = "Helvetica"),
  legend.key = element_rect(fill = NA , size= 2.5)
) 

#figS1cd
for (i in c(1:nrow(samples))){
  qc <- read.delim( file = paste0(out.dir ,"/acr_tss_QC_", samples[i,1],".txt"), header = T, row.names = 1)
  qc$tss.frac <- qc$tss / qc$total
  qc$total.log <- log10(qc$total)
  qc$acr.frac <- qc$acrs / qc$total
  cells <- colnames(data)[grep(samples[i,1], colnames(data))]
  cells <- gsub(paste0(samples[i,1], "_"), "", cells)
  qc.filtered <- qc[cells ,] 
  
  q <- ggplot(qc, aes(x = total.log, y = tss.frac)) + geom_bin2d(bins = 70) + scale_fill_continuous(type = "viridis") + xlim(1.5, 4) + optns
  ggsave(paste0(out.dir ,"/tss_percent_", samples[i,1],".pdf"), height = 5, width = 7 , q)
  q <- ggplot(qc, aes(x = total.log, y = acr.frac)) + geom_bin2d(bins = 70) + scale_fill_continuous(type = "viridis") + xlim(1.5, 4) + optns
  ggsave(paste0(out.dir ,"/FRiP_", samples[i,1],".pdf"), height = 5, width = 7 , q)
}


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##Cell filtering 
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Doublet filtering 
x <- data@meta.data[, grep("DF.classifications", colnames(data@meta.data))]
data$doublet <- x[!is.na(x)]
singlets <- colnames(data)[data$doublet=="Singlet"]
data <- subset(x = data, cells = singlets)

# Additional filtering
a <- rownames(data@meta.data)[data@meta.data[, "nCount_RNA"] < 7000 & 
                                data@meta.data[, "nCount_RNA"] > 200 &
                                data@meta.data[, "nFeature_RNA"] > 180 &
                              data@meta.data[, "nCount_ATAC"] > 200 &
                              data@meta.data[, "nCount_ATAC"] < 20000 &
                                data@meta.data[, "nFeature_ATAC"] > 150 &
                                data@meta.data[, "nCount_ATAC"] > 200 &
                                data@meta.data[, "percent.mt"] < 10 ]

data.filter <- subset(x = data, cells = a)

saveRDS(data.filter, paste0(out.dir.seurat, "/combined_filtered.rds"))

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## post filtering and plots; FigS1
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
data <- readRDS(paste(mainDir, "_seurat_object","combined_filtered.rds", sep = "/"))

col.sample2 <- c("Mock_rep1" = "#005b96",
                "Mock_rep2" = "#005b96",
                "00_Mock" = "#005b96",
                "DC3000_04h" = "#d9c99e",
                "DC3000_06h" = "#bda155",
                "DC3000_09h" = "#a2790d",
                "DC3000_24h"= "#513c06",
                "AvrRpt2_04h"= "#66b2b2",
                "AvrRpt2_06h"= "#008080",
                "AvrRpt2_09h"= "#006666",
                "AvrRpt2_24h"= "#004c4c",
                "AvrRpt2_04h_rep2"= "#66b2b2",
                "AvrRpt2_06h_rep2"= "#008080",
                "AvrRpt2_09h_rep2"= "#006666",
                "AvrRpt2_24h_rep2"= "#004c4c",
                "AvrRpm1_04h"= "#ca6666",
                "AvrRpm1_06h"= "#af1919",
                "AvrRpm1_09h"= "#740000",
                "AvrRpm1_24h"= "#420000"
                
)


data2 <- data
data2$sample <- gsub("00_Mock", "Mock", data2$sample)
Idents(data2) <- data2$sample2
q <- VlnPlot(
  object = data2,
  features = c("nCount_RNA",   "nCount_ATAC", "nCount_peaks", "nFeature_RNA","nFeature_ATAC", "nFeature_peaks"),
  ncol = 3,
  pt.size = 0
 
  
)
q <-q & scale_fill_manual(values = col.sample2)
ggsave("figS1a.pdf", height = 5, width = 12 ,q)


# Number of cells 
# Set theme options for ggplot2
optns2 <- theme(
  axis.text.x = element_text(
    margin = margin(c(5, 0, 0, 0)), 
    size = 10, 
    face = "bold",
    angle = 45, 
    hjust = 1
  ),
  axis.text.y = element_text(
    margin = margin(c(0, 5, 0, 0)),
    size = 10, 
    face = "bold"
  ),
  axis.title = element_blank(),
  axis.ticks = element_line(size = .5),
  axis.ticks.length = unit(.1, "cm"),
  axis.line = element_line(size = .5),
  panel.background = element_blank(),
  legend.text = element_text(
    size = 15, 
    face = "bold", 
    family = "Helvetica", 
    margin = margin(c(20, 0, 20, 0))
  ),
  legend.title = element_text(
    size = 15, 
    face = "bold", 
    family = "Helvetica"
  ),
  legend.key = element_rect(
    fill = NA, 
    size = 2.5
  )
) 

# Reordering samples and calculating the number of cells
idx <- c(1, 2, 12:15, 7:11, 3:6) # Sample order
ncells <- table(data$sample)[idx] %>% as.data.frame()

# Adjusting frequencies
ncells[2, 2] <- ncells[2, 2] + ncells[1, 2]
ncells[10, 2] <- ncells[10, 2] + ncells[9, 2]
ncells <- ncells[-c(1, 9), ]

# Cleaning sample names
ncells$Var1 <- gsub("_rep*.", "", ncells$Var1)
levels(ncells$Var1)[1] <- "Mock"

# Plotting
q <- ggplot(ncells, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  optns2 +
  NoLegend() +
  scale_fill_manual(values = col.sample2)

# Saving the plot
ggsave("figS1b.pdf", plot = q, height = 2.5, width = 3.3)


