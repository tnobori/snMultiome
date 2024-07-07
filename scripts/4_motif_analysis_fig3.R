##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## This script performs motif analyses
## Figures produced with this script: Fig3, FigS5

## A fully processed Seurat object can be downloaded at http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/processed_seurat_object/combined_filtered.rds

## Contact: Tatsuya Nobori (tatsuyanobori@gmail.com)
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Clear all objects from R environment
rm(list = ls())

# Load configuration script for multiome analysis
source("http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/scripts/_config_multiome.R")

## Assuming you are in the directory where this and other R scripts are stored.
mainDir <- file.path(getwd(), "data_out")
subDir <- "4_motif"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

setwd(file.path(mainDir, subDir))

# Theme options for plots
optns <- theme(
  axis.text.x = element_text(margin = margin(c(20, 0, 0, 0)), size = 20, face = "bold", family = "Helvetica"),
  axis.text.y = element_text(margin = margin(c(0, 20, 0, 0)), size = 20, face = "bold", family = "Helvetica"),
  axis.title = element_text(size = 20, face = "bold", family = "Helvetica"),
  axis.ticks.length = unit(.3, "cm"),
  axis.line = element_line(size = 1),
  panel.background = element_blank(),
  legend.text = element_text(size = 15, face = "bold", family = "Helvetica", margin = margin(c(20, 0, 20, 0))),
  legend.title = element_text(size = 15, face = "bold", family = "Helvetica"),
  legend.key = element_rect(fill = NA, size = 2.5)
)
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## data loading
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
# Load dataset
data <- readRDS(paste(mainDir, "_seurat_object","combined_filtered.rds", sep = "/"))

# Load annotation
annotation <- Annotation(data@assays$ATAC)

# Directory path for linkage data
link_dir <- file.path(mainDir, "3_linkage/")

# Load linkage data for mock and treatments
links_mock_all <- read.delim(file = paste0(link_dir ,"linkage_mock_ALL.txt"), head = TRUE, row.names = 1)
links_kt56_all <- read.delim(file = paste0(link_dir ,"linkage_kt56_ALL.txt"), head = TRUE, row.names = 1)
links_kt57_all <- read.delim(file = paste0(link_dir ,"linkage_kt57_ALL.txt"), head = TRUE, row.names = 1)
links_kt58_all <- read.delim(file = paste0(link_dir ,"linkage_kt58_ALL.txt"), head = TRUE, row.names = 1)

# Combine linkage data and filter for significance
links.comb <- rbind(links_mock_all, links_kt56_all, links_kt57_all, links_kt58_all)
links.comb.significant <- links.comb[links.comb$score > 0.1 & links.comb$width < 5000 ,]


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## pre processing
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

DefaultAssay(data) <- "peaks"
seqnames(BSgenome.Athaliana.TAIR.TAIR9) <- gsub("^Chr", "", seqnames(BSgenome.Athaliana.TAIR.TAIR9)) ##need this change to avoid an error 

# Load JASPAR2020 PWMs for Arabidopsis thaliana
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 3702, all_versions = FALSE)
)

# Add motif information to the dataset
data <- AddMotifs(data, genome = BSgenome.Athaliana.TAIR.TAIR9, pfm = pwm)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Function for motif analysis between two clusters
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

motif.analysis <- function(clst1, clst2) {
  tic()
  # Setting the active identity and default assay for the analysis
  Idents(data) <- data$SCT_snn_res.1
  DefaultAssay(data) <- "peaks"
  
  # Prepare output directory
  out_dir <- paste0("clst", clst1, "vs", clst2, "/")
  dir.create(out_dir, showWarnings = TRUE)
  
  print("Motif identification in progress...")
  
  # Identifying differentially accessible peaks
  da_peaks <- FindMarkers(
    object = data,
    ident.1 = clst1,
    ident.2 = clst2,
    only.pos = TRUE,
    test.use = 'LR',
    min.pct = 0.05,
    latent.vars = 'nCount_peaks'
  )
  
  print("Motif enrichment analysis (part 1/2)...")
  
  # Test enrichment of motifs in differentially accessible peaks
  enriched.motifs <- FindMotifs(
    object = data,
    features = rownames(da_peaks)
  )
  write.table(enriched.motifs, file = paste0(out_dir, "enriched_motif.txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  
  # Plot enriched motifs
  q <- MotifPlot(
    object = data,
    motifs = head(rownames(enriched.motifs), n = 5) # Assuming you want to plot top 5 motifs
  )
  ggsave(paste0(out_dir, "enriched_motif.pdf"), plot = q, height = 5, width = 10)
  
  print("Motif enrichment analysis (part 2/2)...")
  
  # Overlap marker peaks with CREs and test for enrichment
  links.comb.marker <- da_peaks[rownames(da_peaks) %in% links.comb.significant$peak,]
  
  enriched.motifs.CRES <- FindMotifs(
    object = data,
    features = rownames(links.comb.marker)
  )
  write.table(enriched.motifs.CRES, file = paste0(out_dir, "enriched_motif_CREs.txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  
  q <- MotifPlot(
    object = data,
    motifs = head(rownames(enriched.motifs.CRES), n = 5) # Assuming you want to plot top 5 motifs
  )
  ggsave(paste0(out_dir, "enriched_motif_CREs.pdf"), plot = q, height = 5, width = 10)
  
  dm.table <- read.delim(file = paste0(out_dir, "enriched_motif_CREs.txt"), header = TRUE, row.names = 1)
  dm.table$name.show <- NA
  gene_to_show <- dm.table$motif.name[1:10]
  dm.table$name.show[match(gene_to_show, dm.table$motif.name)] <- gene_to_show
  
  q <- ggplot(dm.table, aes(x = fold.enrichment, y = -log10(pvalue), label = name.show)) + 
    geom_point(color = as.vector(cols[clst1 + 1])) + 
    geom_text_repel(force = 50, size = 3, hjust = 0.5) +
    optns
  ggsave(paste0(out_dir, "dif_motif_plot.pdf"), plot = q)
  
  q <- MotifPlot(
    object = data,
    motifs = rownames(dm.table)[match(gene_to_show[1:4], dm.table$motif.name)],
    assay = "peaks"
  )
  ggsave(paste0(out_dir, "dif_motif_plot_topLogo.pdf"), plot = q)
  
  toc()
}

# Fig3a
motif.analysis(3, 1)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## building CROMVAR assay
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
DefaultAssay(data) <- "peaks"

# Preparing data matrices and adjusting chromosome names for CROMVAR analysis
motif.matrix <- GetMotifData(object = data, slot = "data")
peak.matrix <- GetAssayData(object = data, slot = "counts")
idx.keep <- rowSums(x = peak.matrix) > 0
peak.matrix <- peak.matrix[idx.keep, , drop = FALSE]
motif.matrix <- motif.matrix[idx.keep, , drop = FALSE]
peak.ranges <- granges(x = data)
peak.ranges <- peak.ranges[idx.keep]

# Creating a SummarizedExperiment object for chromVAR analysis
chromvar.obj <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = peak.matrix),
  rowRanges = peak.ranges)
# Adding GC bias to the SummarizedExperiment object
chromvar.obj <- chromVAR::addGCBias(
  object = chromvar.obj,
  genome = BSgenome.Athaliana.TAIR.TAIR9)

# Fixing NA values which can cause errors in background calculation
row.data <- data.frame(rowData(chromvar.obj))
row.data[is.na(row.data)] <- 0
rowData(chromvar.obj) <- row.data

# Calculating background peaks and deviations for chromVAR analysis
bg <- chromVAR::getBackgroundPeaks(
  object = chromvar.obj)

dev <- chromVAR::computeDeviations(
  object = chromvar.obj,
  annotations = motif.matrix,
  background_peaks = bg
)

# Extracting chromVAR z-scores and updating the dataset with chromvar assay
chromvar.z <- SummarizedExperiment::assays(dev)[[2]]
rownames(x = chromvar.z) <- colnames(x = motif.matrix)
obj <- CreateAssayObject(data = chromvar.z)#counts = chromvar.z was not in the original script
data[["chromvar"]] <- obj

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##  CROMVAR analysis
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

DefaultAssay(data) <- 'chromvar'

data@assays$chromvar@data[is.na(data@assays$chromvar@data)] <- 0 #replacing NA to 0 as the FeaturePlot function does not take NAs

# Generate combined motif-gene name identifiers
motifs <- data[["peaks"]]@motifs@motif.names
motif.id <- names(motifs)
motif.comb <- paste(motif.id, motifs, sep = "_")
rownames(data@assays$chromvar@data)[match(motif.id, rownames(data@assays$chromvar@data))] <- motif.comb #updating motif ids with combined names

out.dir.seurat <- file.path(mainDir, "/_seurat_object")
saveRDS(data, paste0(out.dir.seurat, "/combined_filtered.rds"))

###making an object for each strain
cell_kt56 <- rownames(data@meta.data)[ grepl("DC3000", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] #4100
cell_kt57 <- rownames(data@meta.data)[ grepl("AvrRpt2", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] #4100
cell_kt58 <- rownames(data@meta.data)[ grepl("AvrRpm1", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] #4100

data_kt56 <- subset(x = data, cells = cell_kt56)
data_kt57 <- subset(x = data, cells = cell_kt57)
data_kt58 <- subset(x = data, cells = cell_kt58)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##  plot motif activity and TF expression together 
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

motif.tf.exp <- function(gene_list, fname){
  
  input_motif <- motif.comb[match(gene_list, motifs)]
  ##################################
  DefaultAssay(data_kt56) <- 'chromvar'  
  
  tryCatch({
    q1 <- FeaturePlot(
      object = data_kt56,
      features = input_motif,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = .3,
      split.by = "sample2",
      order = TRUE,
      cols = c("lightgrey","blue"),
      reduction = "umap"
    )
    q1 <- q1 & scale_colour_viridis_c()& NoAxes()
    
    DefaultAssay(data_kt56) <- 'SCT'
    q2 <- FeaturePlot(data_kt56, features = gene_list, split.by = "sample2", pt.size = .3, order=TRUE,
                      min.cutoff = 'q1',
                      max.cutoff = 'q99',
                      cols = c("lightgrey","darkmagenta"),
                      reduction = "umap"
    )
    q2 <- q2 & scale_colour_gradientn(
      colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
      na.value = "lightgray"
    ) & NoAxes()
    p1 <- CombinePlots(list(q1, q2), nrow  =2)
    ggsave(filename = paste("motif_tf_plot/","_", fname,"_feture_KT56.pdf", sep = ""), width = 16, height = 8, p1) 
  }, error=function(e){})
  
  ################################## ################################## ################################## ##################################
  DefaultAssay(data_kt57) <- 'chromvar'  
  tryCatch({
    q1 <- FeaturePlot(
      object = data_kt57,
      features = input_motif,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = .3,
      split.by = "sample2",
      order = TRUE,
      cols = c("lightgrey","blue"),
      reduction = "umap"
    )
    q1 <- q1 & scale_colour_viridis_c()& NoAxes()
    
    DefaultAssay(data_kt57) <- 'SCT'
    q2 <- FeaturePlot(data_kt57, features = gene_list, split.by = "sample2", pt.size = .3, order=TRUE,
                      min.cutoff = 'q1',
                      max.cutoff = 'q99', cols = c("lightgrey","darkmagenta"),
                      reduction = "umap"
    )
    q2 <- q2 & scale_colour_gradientn(
      colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
      na.value = "lightgray"
    ) & NoAxes()
    p1 <- CombinePlots(list(q1, q2), nrow  =2)
    ggsave(filename = paste("motif_tf_plot/","_", fname,"_feture_KT57.pdf", sep = ""), width = 16, height = 8, p1)
  }, error=function(e){})
  ################################## ################################## ##################################
  DefaultAssay(data_kt58) <- 'chromvar'  
  tryCatch({
    q1 <- FeaturePlot(
      object = data_kt58,
      features = input_motif,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = .3,
      split.by = "sample2",
      order = TRUE,
      cols = c("lightgrey","blue"),
      reduction = "umap"
    )
    q1 <- q1 & scale_colour_viridis_c()& NoAxes()
    
    DefaultAssay(data_kt58) <- 'SCT'
    q2 <- FeaturePlot(data_kt58, features = gene_list, split.by = "sample2", pt.size = .3, order=TRUE,
                      min.cutoff = 'q1', cols = c("lightgrey","darkmagenta"),
                      max.cutoff = 'q90',
                      reduction = "umap"
    )
    q2 <- q2 & scale_colour_gradientn(
      colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
      na.value = "lightgray"
    ) & NoAxes()
    p1 <- CombinePlots(list(q1, q2), nrow  =2)
    ggsave(filename = paste("motif_tf_plot/","_", fname,"_feture_KT58.pdf", sep = ""), width = 16, height = 8, p1)
  }, error=function(e){})
  
  
  
  ################################## ################################## ################################## ##################################
  ##plots for legend
  DefaultAssay(data_kt57) <- 'chromvar'  
  tryCatch({
    q1 <- FeaturePlot(
      object = data_kt57,
      features = input_motif,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = 1,
      # split.by = "sample",
      order = TRUE,
      cols = c("lightgrey","blue"),
      reduction = "umap"
    )
    q1 <- q1 & scale_colour_viridis_c()& NoAxes()
    
    DefaultAssay(data_kt57) <- 'SCT'
    q2 <- FeaturePlot(data_kt57, features = gene_list, 
                      # split.by = "sample", 
                      pt.size = 1, order=TRUE,
                      min.cutoff = 'q1',
                      max.cutoff = 'q99', cols = c("lightgrey","darkmagenta"),
                      reduction = "umap"
    )
    q2 <- q2 & scale_colour_gradientn(
      colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
      na.value = "lightgray"
    ) & NoAxes()
    p1 <- CombinePlots(list(q1, q2), nrow  =2)
    ggsave(filename = paste("motif_tf_plot/","_", fname,"_feture_KT57_legend.pdf", sep = ""), width = 5, height = 8, p1)
  }, error=function(e){})
  ################################## ################################## ##################################
}

#fig4e
motif.tf.exp("WRKY46", "WRKY46")


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##  cluster specific motif analysis 
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
dir.create("marker_motif_plot", showWarnings = FALSE)
Idents(data) <- data$SCT_snn_res.1
marker.motif <- FindAllMarkers(data, assay = "chromvar", only.pos = TRUE )
write.table(marker.motif, file = "marker_motif_plot/markers_RNA.txt", row.names=T, col.names=NA, sep="\t", quote=F)

marker.motif.top <- marker.motif %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

motif.activity <- function(gene_list, fname, subdir){
  
  input_motif <- gene_list
  ##################################
  DefaultAssay(data_kt56) <- 'chromvar'  
  
  tryCatch({
    q1 <- FeaturePlot(
      object = data_kt56,
      features = input_motif,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = 1,
      split.by = "sample2",
      order = TRUE,
      cols = c("lightgrey","blue"),
      reduction = "umap"
    )
    q1 <- q1 & scale_colour_viridis_c()& NoAxes()
    
    ggsave(filename = paste("marker_motif_plot/", subdir, "/_", fname,"_feture_KT56.png", sep = ""), width = 16, height = 4, q1) 
  }, error=function(e){})
  
  ################################## ################################## ################################## ##################################
  DefaultAssay(data_kt57) <- 'chromvar'  
  tryCatch({
    q1 <- FeaturePlot(
      object = data_kt57,
      features = input_motif,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = 1,
      split.by = "sample2",
      order = TRUE,
      cols = c("lightgrey","blue"),
      reduction = "umap"
    )
    q1 <- q1 & scale_colour_viridis_c()& NoAxes()
    
    ggsave(filename = paste("marker_motif_plot/", subdir, "/_", fname,"_feture_KT57.png", sep = ""), width = 16, height = 4, q1)
  }, error=function(e){})
  ################################## ################################## ##################################
  DefaultAssay(data_kt58) <- 'chromvar'  
  tryCatch({
    q1 <- FeaturePlot(
      object = data_kt58,
      features = input_motif,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = 1,
      split.by = "sample2",
      order = TRUE,
      cols = c("lightgrey","blue"),
      reduction = "umap"
    )
    q1 <- q1 & scale_colour_viridis_c()& NoAxes()
    
    ggsave(filename = paste("marker_motif_plot/", subdir, "/_", fname,"_feture_KT58.png", sep = ""), width = 16, height = 4, q1)
  }, error=function(e){})
}

for (i in c(1:nrow(marker.motif.top))){
  motif.activity(gene_list = marker.motif.top$gene[i], 
               fname = marker.motif.top$gene[i],
               subdir =marker.motif.top$cluster[i] )
}




optns <- theme(
  axis.text.x = element_text(margin = margin(c(5, 0, 0, 0)), size = 10, face = "bold", family = "Helvetica"),
  axis.text.y = element_text(margin = margin(c(0, 5, 0, 0)) ,size = 10, face = "bold", family = "Helvetica"), 
  axis.title = element_text(size = 5, face = "bold", family = "Helvetica"),
  axis.ticks.y = element_line(size = 1),
  axis.ticks.x = element_line(size = 1),
  axis.ticks.length = unit(.3, "cm"),
  axis.line  = element_line(size = 1),
  panel.background = element_blank(),
  legend.text = element_text(size = 15, face = "bold", family = "Helvetica", margin = margin(c(20, 0, 20, 0))),
  legend.title = element_text(size = 15, face = "bold", family = "Helvetica"),
  legend.key = element_rect(fill = NA , size= 2.5)
) 
marker.motif.top1 <- marker.motif %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
q <- DotPlot(data,assay = "chromvar", features = unique(marker.motif.top1$gene), cols = c("#ffcfff", "#221100"), dot.scale = 10) + coord_flip() + optns 
ggsave(paste0( "marker_motif_plot/fig3b.pdf"), height = 5, width = 11 , q)


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##  CROMVAR correlation between chromVAR motif activity and TF expression...to further filter interesting TFs
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
dir.create("motif_tf_correlation", showWarnings = FALSE)
DefaultAssay(data)<- "SCT"

# JASPAR motif names had to be fixed to match gene names used in our analysis. 
fixed <- read.delim("http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/misc/motif_name_fix_added.txt")
rownames(data@assays$chromvar@data) <- paste(fixed$X, fixed$fixed, sep = "_")
saveRDS(data, paste0(out.dir.seurat, "/combined_filtered.rds"))

motif.name <- data@assays$peaks@motifs@motif.names %>% unlist() %>% as.data.frame()
motif.name$comb <- paste(rownames(motif.name), motif.name[, 1], sep = "_")

data$celltype <- "Mesophyll"
data$celltype[data$SCT_snn_res.1 %in% c(0,12,19,21,29)] <- "Epidermis"
data$celltype[data$SCT_snn_res.1 %in% c(6,9,10,14)] <- "Vasculature"
data$celltype[data$SCT_snn_res.1 %in% c(16)] <- "Unknown"

## tissue wise
DefaultAssay(data) <- "RNA"
data2 <- data
data2 <- NormalizeData(data)

DefaultAssay(data2) <- "RNA"
tissue <- data2$celltype %>% unique() %>% sort()


for (i in c(1:length(tissue))){
  data.clst <- data2[, data2$celltype==tissue[i]] 
  data.clst.chromvar <- data.clst@assays$chromvar@data
  data.clst.rna <- data.clst@assays$RNA@data
  
  a <- data.clst.chromvar[motif.name[, 1] %in% rownames(data.clst.rna), ]
  b <- data.clst.rna[motif.name[motif.name[, 1] %in% rownames(data.clst.rna) ,1] ,]
  
  out <- c()
  for (j in c(1:nrow(a))){
    out[j] <- cor(b[j,], a[j,])
  }
  names(out) <- rownames(b)
  out <- sort(out, decreasing = T)
  write.table(out, file=paste0("motif_tf_correlation/correlation_motif_TFexp_celltype_", tissue[i], ".txt"), row.names=T, col.names=NA, sep="\t", quote=F)
}

##combining correlation in each cell type 
clst.summary <- matrix(data = NA, nrow = nrow(b), ncol = length(tissue))
rownames(clst.summary) <- rownames(b)
colnames(clst.summary) <- tissue

for (i in c(1:length(tissue))){
  data.clst <- data2[, data2$celltype==tissue[i]] 
  data.clst.chromvar <- data.clst@assays$chromvar@data
  data.clst.rna <- data.clst@assays$RNA@data
  
  a <- data.clst.chromvar[motif.name[, 1] %in% rownames(data.clst.rna), ]
  b <- data.clst.rna[motif.name[motif.name[, 1] %in% rownames(data.clst.rna) ,1] ,]
  
  out <- c()
  for (j in c(1:nrow(a))){
    out[j] <- cor(b[j,], a[j,])
  }
  names(out) <- rownames(b)
  
  clst.summary[, i] <- out
}

clst.summary2 <- clst.summary
clst.summary2[is.na(clst.summary2)] <- 0

#plotting 

pdf("motif_tf_correlation/fig3d.pdf")
pheatmap(clst.summary2[apply(clst.summary2, 1, max)>0.12,],
         #scale = "row",
         cluster_cols=TRUE,
         cluster_rows = TRUE,
         cellwidth = 10,
         #cellheight = .07,
         # color = colorRampPalette(c("magenta", "black",  "green"))(100),
         color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),
         #color = hcl.colors(50, "Inferno"),
         breaks= seq(-0.1, 0.2, length.out = 100),
         #annotation_row = an_deg,
         #annotation_colors = an_color,
         fontsize_row = 10,
         fontsize_col = 10
)
dev.off()


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##  motif enrichment score heatmap
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
DefaultAssay(data) <- "chromvar"

motif.mat <- data@assays$chromvar@data
mean <- apply(motif.mat, 1, mean)
sd <- apply(motif.mat, 1, sd)
motif.mat.z <- (motif.mat - mean)/sd


library(viridis)


ann_col <- as.data.frame(colnames(motif.mat.z)) 
colnames(ann_col) <- "rowname"
ann_col$clst <- data$SCT_snn_res.1
ann_col <- ann_col[, -1] %>% as.data.frame()
rownames(ann_col) <- colnames(data)


list_color <- list(. = cols)

png("motif_heatmap/figS5a.png")
pheatmap(motif.mat.z[, order(data$SCT_snn_res.1)],
         #scale = "row",
         cluster_cols=FALSE,
         cluster_rows = TRUE,
         #cellwidth = 70,
         #cellheight = .07,
         # color = colorRampPalette(c("magenta", "black",  "green"))(100),
         # color = magma(11),
         color = colorRampPalette(brewer.pal(n = 9, name ="YlGnBu"))(100),
         #color = hcl.colors(50, "Inferno"),
         breaks= seq(-2, 2, length.out = 100),
         #annotation_row = an_deg,
         annotation_col = ann_col,
         annotation_colors = list_color,
         fontsize_row = .00001,
         fontsize_col = .00001,
         # gaps_col = c(3,5)
         # gaps_col = as.vector(table(data$SCT_snn_res.1)),
         legend = F
         
)
dev.off()

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## TF heatmap
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
tfs <- read.delim("http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/misc/Thale_cress-transcription%20factor.txt", header = F)
tfs <- tfs[grep("\\.1", tfs[,1]) ,]
tfs[,1] <- gsub("\\.[1-9]", "", tfs[,1])
tfs <- tfs[, c(1:2)]
rownames(tfs) <- tfs[,1]
tfs_ann <- as.data.frame(tfs[,2])
rownames(tfs_ann) <- rownames(tfs)


x = rownames(tfs_ann)
xx <- c()
for (i in c(1:length(x))){
  if (x[i] %in% rownames(data)){
    xx[i] <- x[i]
    
  } else {
    xx[i] <- annotation$gene_name[grep(x[i], annotation$ID)][1]
  }
}
xx <- xx %>% na.omit()

DefaultAssay(data) <- "SCT"

q <- DoHeatmap(data, features = xx, raster = T, group.by = "SCT_snn_res.1", disp.max= 1, disp.min = -1,  group.colors = cols, size = 3,  draw.lines = F )
ggsave(filename = "TF_heatmap/figS5b.pdf", q)


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## genes linked with motifs; fig3gh related
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
dir.create("motif_linked_genes", showWarnings = FALSE)

# Function to link genes with specific motifs
motif.link.genes <- function(target.motif){
  DefaultAssay(data) <- "peaks"
  motifID <- rownames(motif.name)[motif.name$.==target.motif]
  # Extract peaks associated with the motif
  motif.peak.table <- Motifs(data)@data
  target.peaks <- motif.peak.table[, motifID]
  target.peaks <- target.peaks[target.peaks==1]
 
   # Find links between target peaks and genes
  target.links <- links.comb[links.comb$peak %in% names(target.peaks) ,]
  target.links.filter <- target.links[target.links$score > 0.1 & target.links$width < 5000 ,]
  write.table(target.links.filter, file=paste0("motif_linked_genes/", target.motif,".txt"), row.names=T, col.names=NA, sep="\t", quote=F)
  target.links.filter$gene %>% table() %>% sort(decreasing = T)
}

motif.link.genes("WRKY46") #Fig3g is generated based on this output using Cytoscape 

# Perform analysis for top motifs and all motifs
top.motifs <- rownames(clst.summary2[apply(clst.summary2, 1, max)>0.12,])
all.motifs <- motif.name$.

for (motif in top.motifs){
  motif.link.genes(motif)
}

for (motif in all.motifs){
  motif.link.genes(motif)
}

##GO analysis of motif linked genes

# Function for GO enrichment analysis
go.motif.link <- function(target.motif){
  gene <- read.delim(file = paste0("motif_linked_genes/", target.motif,".txt"), header = T, row.names = 1)
  gene <- gene$gene %>% unique()
  
  in.gene <- c()
  
  print("gene names to ID...")
  #gene names to gene ID
  for (j in c(1:length(gene))){
    in.gene<- c(in.gene,  genes_to_id(gene[j]))
    
  }
  
  print("GO enrichment...")
  ego <- enrichGO(gene          = in.gene,
                  OrgDb         = org.Athaliana.eg.db,
                  ont           = "BP",
                  keyType = 'GID',
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  pdf(paste0("motif_linked_genes/GO_analysis/", target.motif,".pdf"), 
      # width = 1000, height = 1000
  )
  try(print(barplot(ego, showCategory=10)))
  dev.off()
  
  write.table(ego@result,file = paste0("motif_linked_genes/GO_analysis/",target.motif, ".txt") ,row.names=T, col.names=NA, sep="\t", quote=F)
}

for (motif in top.motifs){
  try(go.motif.link(motif))
}


#motif GO summary plot
#go data viz; take top 2 go terms per cluster
top.motifs <- rownames(clst.summary2[apply(clst.summary2, 1, max)>0.12,])
tfs.sel <- top.motifs %>% sort()
out <- matrix(data = NA, nrow = length(tfs.sel)*2, ncol = 3) %>% as.data.frame()
colnames(out) <- c("Description", "Size", "adjusted pvalue")

optns3 <- theme(
  axis.text.x = element_text(margin = margin(c(10, 0, 0, 0)), size = 15, face = "bold", family = "Helvetica", angle = 90, hjust = 1),
  axis.text.y = element_text(margin = margin(c(0,10, 0, 0)) ,size = 15, face = "bold", family = "Helvetica"), 
  axis.title = element_text(size = 10, face = "bold", family = "Helvetica"),
  axis.title.x = element_text(margin = margin(c(10, 0, 0, 0))),
  axis.title.y = element_text(margin = margin(c(0, 10, 0, 0))),
  axis.ticks.y = element_line(size = 1),
  axis.ticks.x = element_line(size = 1),
  axis.ticks.length = unit(.25, "cm"),
  axis.line  = element_line(size = 1),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(color = 'gray', size = .3, linetype = "dotted"),
  legend.box.background = element_blank(), # get rid of legend panel bg,
  legend.background = element_blank() # get rid of legend bg
) 

for (i in c(1:length(tfs.sel))){
  tryCatch({go.data <- read.delim(file = paste0("motif_linked_genes/GO_analysis/", tfs.sel[i],".txt"), header = T, row.names = 1)
  
  s <- 2*(i-1)+1
  
  out[s, ] <- c(go.data$Description[1],as.numeric(go.data$Count[1]), as.numeric(go.data$p.adjust[1]))
  out[s+1, ] <- c(go.data$Description[2],as.numeric(go.data$Count[2]), as.numeric(go.data$p.adjust[2]))
  }, error=function(e){})
}

out$clst <- as.character(rep(tfs.sel, each = 2))
#Then turn it back into a factor with the levels in the correct order
out$clst <- factor(out$clst, levels=unique(out$clst))

out$log.p <- -log10(as.numeric(out$`adjusted pvalue`))
out$Size <- as.numeric(out$Size)
out <- out %>% na.omit()

go.sel <- c("regulation of immune response", "systemic acquired resistance", "regulation of jasmonic acid mediated signaling pathway", "response to salicylic acid",
            "aging")

q <- ggplot(out, aes(x = clst, y = Description, size = Size , color = log.p)) + geom_point()+ optns3 + 
  scale_color_gradient(low = "grey",high = "#418557",limits = c(0,max(out$log.p)))+
  scale_size(range = c(0,as.numeric(quantile(out$Size)[3])))
ggsave("motif_linked_genes/GO_analysis/figS5d.pdf", height = 6, width = 16 ,q)

#genes commonly targeted by top TFs...
motif.link.genes2 <- function(target.motif){
  DefaultAssay(data) <- "peaks"
  motifID <- rownames(motif.name)[motif.name$.==target.motif]
  
  motif.peak.table <- Motifs(data)@data
  target.peaks <- motif.peak.table[, motifID]
  target.peaks <- target.peaks[target.peaks==1]
  
  target.links <- links.comb[links.comb$peak %in% names(target.peaks) ,]
  target.links.filter <- target.links[target.links$score > 0.15 & target.links$width < 5000 ,]
  # write.table(target.links.filter, file=paste0("motif_linked_genes/", target.motif,".txt"), row.names=T, col.names=NA, sep="\t", quote=F)
  table(target.links.filter$gene)
}

x <- clst.summary2[apply(clst.summary2, 1, max)>0.12,] #topmotifs
motif.sel<- rownames(x)

l <- list()
for (motif in motif.sel){
  l[[motif]] <- motif.link.genes2(motif)
}

lnames <- names(l)

table.list <- data.frame(TF = rep(names(l[1]), times = length(l[[1]])), target = names(l[[1]]), nlinks = as.vector(l[[1]]))
for (i in c(2:length(l))){
  ntable <- data.frame(TF = rep(names(l[i]), times = length(l[[i]])), target = names(l[[i]]), nlinks = as.vector(l[[i]]))
  table.list <- rbind(table.list, ntable)
}
table.list$interaction <- "tf_gene"
table.list$directed <- "true"


write.table(table.list, file=paste0("motif_linked_genes/table_for_network_cor0-15.txt"), row.names=T, col.names=NA, sep="\t", quote=F)
# write.table(table.list, file=paste0("motif_linked_genes/table_for_network",motif.sel,".txt"), row.names=T, col.names=NA, sep="\t", quote=F)

table.list <- read.delim(file =paste0("motif_linked_genes/table_for_network_cor0-15.txt"), header = T, row.names = 1 )

out <- matrix(data=0, nrow = length(unique(table.list$target)), ncol = length(unique(table.list$TF))) %>% as.data.frame()
colnames(out) <- unique(table.list$TF)
rownames(out) <- unique(table.list$target)

tfs.in <- colnames(out)

for (i in c(1:ncol(out))){
  target.in <- table.list$target[table.list$TF==tfs.in[i]]
  out[rownames(out) %in% target.in, i] <- 1
} 


out.table <- out

#selecting TFs
tf.sel <- c("IDD7", "IDD4", "WRKY46", "WRKY33", "WRKY8", "TGA7", "NAC029", "MYR2", "NAC055", "ERF7"
            , "ABF3", "MYC2", "WRKY75", "IDD5", "WRKY6", "GT-3A","CAMTA3")
out.table <- out.table[, tf.sel]
out.table <- out.table[rowSums(out.table)>0, ]

labs_row.rna <- rep("", times = length(labs))
labs_row.rna[na.omit(match(labs.sel, rownames(motif)))] <- labs.sel

#choose labels to show
d=dist(out.table)
h <- hclust(d)
labs <- h$labels[h$order]

labs.sel <- labs[seq(1, length(labs), by = 4)]

labs_row.rna <- rep("", times = length(labs))
labs_row.rna[na.omit(match(labs.sel, rownames(out.table)))] <- labs.sel

pdf("motif_tf_correlation/fig3h.pdf", height = 8,width = 15 )
pheatmap(out.table,
         cluster_cols=TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(brewer.pal(n = 9, name ="YlGnBu"))(100),
         breaks= seq(0, 1, length.out = 100),
         fontsize_row = 10,
         fontsize_col = 10,
         legend = F,
         treeheight_row = 0,
         treeheight_col = 0,
         labels_row = labs_row.rna
         
)
dev.off()
