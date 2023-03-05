source(".../_config_multiome.R")

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
      a <- read.table(gzfile(as.character(bed)))
    }else{
      a <- read.table(as.character(bed))
    }
  }
  
  if(grepl(".gtf", ann)){
    anntype <- "gtf"
  }else{
    anntype <- "gff3"
  }
  
  if(verbose){message(" - converting fragment file to single-bp resolution Tn5 insertions sites ...")}
  start.coordinates <- data.frame(V1=a$V1, V2=(a$V2), V3=(a$V2+1), V4=a$V4, V5="+")
  start.coordinates <- start.coordinates[!duplicated(start.coordinates),]
  
  end.coordinates <- data.frame(V1=a$V1, V2=a$V2, V3=(a$V3-1), V4=a$V4, V5="-")
  end.coordinates <- end.coordinates[!duplicated(end.coordinates),]
  all.coordinates <- rbind(start.coordinates, end.coordinates)
  all.coordinates2 <- all.coordinates[sample(1:nrow(all.coordinates), nrow(all.coordinates)*downsample) ,]
  aa <- all.coordinates2
  
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
  
  
  # building table"
  skip.PtMt = TRUE
  bed2 <- aa 
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
                   #ptmt=as.numeric(org.counts),
                   row.names=ids)

  return(df)
}



##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Color scheme
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

col.sample <- c("Mock" = "#005b96",
                "DC3000_04h" = "#d9c99e",
                "DC3000_06h" = "#bda155",
                "DC3000_09h" = "#a2790d",
                "DC3000_24h"= "#513c06",
                "AvrRpt2_04h"= "#66b2b2",
                "AvrRpt2_06h"= "#008080",
                "AvrRpt2_09h"= "#006666",
                "AvrRpt2_24h"= "#004c4c",
                "AvrRpm1_04h"= "#ca6666",
                "AvrRpm1_06h"= "#af1919",
                "AvrRpm1_09h"= "#740000",
                "AvrRpm1_24h"= "#420000")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## load and reformat annotation
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

annotation <- rtracklayer::import( ".../At_reference/Arabidopsis_thaliana.TAIR10.52.gff3")
colnames(annotation@elementMetadata)[colnames(annotation@elementMetadata)=="Name"] <- "gene_name" 
colnames(annotation@elementMetadata)[colnames(annotation@elementMetadata)=="biotype"] <- "gene_biotype" 
levels(annotation$type)[levels(annotation$type)=="gene"]="body"
annotation$gene_name[grepl("exon", annotation$gene_name)] <- gsub("\\.[0-9].exon[0-9]{1,}", "", annotation$gene_name[grepl("exon", annotation$gene_name)]) #remove exon number

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## creating Seurat object for each sample
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# load the RNA and ATAC data

#a list of paths to count and fragment files
path.list <- read.delim(file = ".../raw_path_list.txt", header = T, row.names = 1)

# create a list of seurat objects
list_seurat <- list()

for (i in c(1:nrow(path.list))){
  list_seurat[[i]] <- create.seurat(countpath = path.list$counts[i], fragpath = path.list$fragments[i], annotation = annotation )
}

names(list_seurat) <- rownames(path.list)

# rename cells using object names as prefix
for (i in names(list_seurat)) {
  list_seurat[[i]] <- RenameCells(list_seurat[[i]],
                                         add.cell.id = i)
}

#combining seurat objects
data <- Merge_Seurat_List(
  list_seurat,
  add.cell.ids = NULL,
  merge.data = TRUE,
  project = "PathogenMultiome"
)

#adding other metadata
list.sample <- rownames(path.list) %>% as.list()
samplenames <- lapply(list.sample, function(x){rep(x, times = length(grep(x, colnames(data))))})
data$sample <- unlist(samplenames)
sample.time <- gsub(".*_", "", data$sample)
data$time <- sample.time

#MT score
mt.gene.idx <- c(grep("ORF153A", rownames(data)):grep("ORF204", rownames(data))) #these are MT genes
mt.count <- data@assays$RNA@counts[mt.gene.idx, ] %>% colSums()
percent.mt <- (mt.count / data$nCount_RNA) * 100
data$percent.mt <- percent.mt

#saving
saveRDS(data, "_data_combined_unfiltered.rds")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## ATAC QC analysis
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

gffpath <- ".../fixed_Arabidopsis_thaliana.TAIR10.45_modified_211122.gtf"
chrompath <- ".../at_chromosom_file.txt" #chromosome size file
out.dir <- "path/to/output/directory"

dir.create(out.dir, showWarnings = FALSE)

# run ATAC QC
for (i in c()){
  out <- loadBEDandGenomeData.mod(bed = path.list$fragments[i], ann =gffpath, sizes = chrompath)
  write.table(out, file = paste0(out.dir ,"acr_tss_QC_", rownames(path.list)[i],".txt"), row.names=T, col.names=NA, sep="\t", quote=F)
}

##plotting

optns <- theme(
  axis.text.x = element_text(margin = margin(c(5, 0, 0, 0)), size = 20, face = "bold", family = "Helvetica"),
  axis.text.y = element_text(margin = margin(c(0, 5, 0, 0)) ,size = 20, face = "bold", family = "Helvetica"), 
  axis.title = element_text(size = 20, face = "bold", family = "Helvetica"),
  axis.ticks.y = element_line(size = 1),
  axis.ticks.x = element_line(size = 1),
  axis.ticks.length = unit(.3, "cm"),
  axis.line  = element_line(size = 1),
  panel.background = element_blank(),
  legend.text = element_text(size = 15, face = "bold", family = "Helvetica", margin = margin(c(20, 0, 20, 0))),
  legend.title = element_text(size = 15, face = "bold", family = "Helvetica"),
  legend.key = element_rect(fill = NA , size= 2.5),
  ) 


for (i in c(11:nrow(path.list))){
  qc <- read.delim( file = paste0(out.dir ,"/acr_tss_QC_", rownames(path.list)[i],".txt"), header = T, row.names = 1)
  qc$tss.frac <- qc$tss / qc$total
  qc$total.log <- log10(qc$total)
  qc$acr.frac <- qc$acrs / qc$total
  cells <- colnames(data)[grep(rownames(path.list)[i], colnames(data))]
  cells <- gsub(paste0(rownames(path.list)[i], "_"), "", cells)
  qc.filtered <- qc[cells ,] 
  
  q <- ggplot(qc, aes(x = total.log, y = tss.frac)) + geom_bin2d(bins = 70) + scale_fill_continuous(type = "viridis") + xlim(1.5, 4) + optns
  ggsave(paste0(out.dir ,"/tss_percent_", rownames(path.list)[i],".pdf"), height = 5, width = 7 , q)
  q <- ggplot(qc, aes(x = total.log, y = acr.frac)) + geom_bin2d(bins = 70) + scale_fill_continuous(type = "viridis") + xlim(1.5, 4) + optns
  ggsave(paste0(out.dir ,"/FRiP_", rownames(path.list)[i],".pdf"), height = 5, width = 7 , q)
}

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##Cell filtering 
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

data <- readRDS("_data_combined_unfiltered.rds")

a <- rownames(data@meta.data)[data@meta.data[, "nCount_RNA"] < 7000 & 
                                data@meta.data[, "nCount_RNA"] > 200 &
                                data@meta.data[, "nFeature_RNA"] > 180 &
                              data@meta.data[, "nCount_ATAC"] > 200 &
                              data@meta.data[, "nCount_ATAC"] < 20000 &
                                data@meta.data[, "nFeature_ATAC"] > 150 &
                                data@meta.data[, "nCount_ATAC"] > 200 &
                                data@meta.data[, "percent.mt"] < 5 ]

data.filter <- subset(x = data, cells = a)
saveRDS(data.filter,  "combined_filtered.rds")
