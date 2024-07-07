##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## This script performs motif analyses
## Figures produced with this script: FigS8

## A fully processed Seurat object can be downloaded at http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/processed_seurat_object/combined_filtered.rds

## Contact: Tatsuya Nobori (tatsuyanobori@gmail.com)
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Clear all objects from R environment
rm(list = ls())

# Load configuration script for multiome analysis
source("http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/scripts/_config_multiome.R")

## Assuming you are in the directory where this and other R scripts are stored.
mainDir <- file.path(getwd(), "data_out")
subDir <- "8_gt3aKO_snRNAseq"
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

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##### Data Loading
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Load original dataset
data_orig <- readRDS(paste(mainDir, "_seurat_object", "combined_filtered.rds", sep = "/"))## A fully processed Seurat object can be downloaded at http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/processed_seurat_object/combined_filtered.rds. Or you can create one from scratch by running 1_qc_data_integration_figS1.R

# Load annotation
annotation <- Annotation(data_orig@assays$ATAC)

# Load GT3a data
data <- readRDS(paste(mainDir, "_seurat_object", "gt3ako.rds", sep = "/")) ##download this dataset at http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/processed_seurat_object/gt3ako.rds

# Extract AvrRpt2 data
cell_select <- colnames(data_orig)[c(grep("Mock", data_orig$sample), grep("AvrRpt2", data_orig$sample))]
data_orig_select <- subset(x = data_orig, cells = cell_select)
DefaultAssay(data_orig_select) <- "RNA"
DefaultAssay(data) <- "RNA"
data[['RNA']] <- JoinLayers(data[['RNA']])

# Add missing metadata columns
for (col in setdiff(colnames(data_orig_select@meta.data), colnames(data@meta.data))) {
  data[[col]] <- NA
}
for (col in setdiff(colnames(data@meta.data), colnames(data_orig_select@meta.data))) {
  data_orig_select[[col]] <- NA
}

# Remove unnecessary assays
data_orig_select[['ATAC']] <- NULL
data_orig_select[['peaks']] <- NULL
data_orig_select[['chromvar']] <- NULL
data_orig_select[['atacRNA_400bp']] <- NULL
data_orig_select[['SCT']] <- NULL
data[['SCT']] <- NULL

# Merge datasets
data <- merge(data, data_orig_select)

# Clean up metadata
data@meta.data <- data@meta.data[, -grep("DF.class", colnames(data@meta.data))]
data@meta.data <- data@meta.data[, -grep("pANN_", colnames(data@meta.data))]
data@meta.data <- data@meta.data[, -grep("peaks", colnames(data@meta.data))]

data$sample[grep("gt3a", data$orig.ident)] <- data$orig.ident[grep("gt3a", data$orig.ident)]

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##### Clustering
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# SCTransform and PCA
DefaultAssay(data) <- "RNA"
data <- SCTransform(data, vars.to.regress = "percent.mt")
data <- RunPCA(data)

# Check PCs
DepthCor(data, reduction = "pca")
ElbowPlot(data, reduction = "pca")

# Find neighbors and clusters
data <- FindNeighbors(data, dims = 1:20)
data <- FindClusters(data, resolution = 1, verbose = TRUE)
data <- RunUMAP(data, dims = 1:20, n.neighbors = 20, min.dist = 0.01)

# Plot UMAP
q <- DimPlot(data, pt.size = 2) + optns + NoLegend() + NoAxes() + ggtitle("")
ggsave(filename = "umap_RNA.png", width = 10, height = 10, q)

# Harmony integration
DefaultAssay(data) <- "SCT"
p1 <- DimPlot(object = data, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = data, features = "PC_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1, p2)

data <- data %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, assay.use = "SCT", reduction = "pca", project.dim = FALSE, reduction.save = "harmony.rna", dims.use = 1:20)

harmony_embeddings <- Embeddings(data, 'harmony.rna')

p1 <- DimPlot(object = data, reduction = "harmony.rna", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = data, features = "harmonyrna_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1, p2)

# Check Harmony PCs
DepthCor(data, reduction = "harmony.rna")
ElbowPlot(data, reduction = "harmony.rna")

# Run UMAP with Harmony
data <- data %>% 
  RunUMAP(reduction = "harmony.rna", dims = 1:20, n.neighbors = 30L, min.dist = 0.01) %>% 
  FindNeighbors(reduction = "harmony.rna", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

q <- DimPlot(data, pt.size = .3, reduction = "umap", label = TRUE, cols = cols, group.by = "SCT_snn_res.1", label.size = 2, label.color = "#201c1d", repel = TRUE) + optns + NoLegend() + NoAxes() + ggtitle("")
ggsave(filename = "umap_harmony_RNA_noLegend2.pdf", width = 3, height = 3, q)

# Save Harmony-integrated data
saveRDS(data, file = paste(mainDir, "_seurat_object", "data_harmony_origIntegration.rds", sep = "/"))

# Find all markers
markers.rna <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers.rna, file = "markers_RNA.txt", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##### Subclustering
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Select immune active mesophyll cells
data_sel <- data[, data$SCT_snn_res.1 %in% c(4, 12, 21, 23)]
dir.create("4_12_21_23", showWarnings = FALSE)

DefaultAssay(data_sel) <- "RNA"
data_sel <- SCTransform(data_sel, verbose = FALSE)
data_sel <- RunPCA(data_sel, verbose = FALSE)

data_sel <- FindNeighbors(data_sel, dims = 1:15, verbose = FALSE)
data_sel <- FindClusters(data_sel, resolution = 1, verbose = FALSE)
data_sel <- RunUMAP(data_sel, dims = 1:15, verbose = FALSE)

# Plot UMAP for subclusters
q <- DimPlot(data_sel, pt.size = 3, cols = cols, label = TRUE, repel = TRUE) + optns + NoLegend() + NoAxes() + ggtitle("")
ggsave(filename = paste0("4_12_21_23", "/umap_RNA.png"), width = 10, height = 10, q)

q <- DimPlot(data_sel, reduction = "umap", pt.size = 4, split.by = 'sample', cols = cols, label = TRUE, repel = TRUE) + optns + NoLegend() + NoAxes() + ggtitle("")
ggsave(filename = paste0("4_12_21_23", "/umap_RNA2.png"), width = 40, height = 20, q)

q <- DimPlot(data_sel, reduction = "umap", group.by = "sample", pt.size = 4, label = TRUE, repel = TRUE) + optns + NoLegend() + NoAxes() + ggtitle("")
ggsave(filename = paste0("4_12_21_23", "/umap_RNA3.png"), width = 20, height = 20, q)

# Harmony integration for subclusters
DefaultAssay(data_sel) <- "SCT"
p1 <- DimPlot(object = data_sel, reduction = "pca", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = data_sel, features = "PC_1", group.by = "sample", pt.size = .1)
plot_grid(p1, p2)

data_sel <- data_sel %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT", verbose = FALSE)

harmony_embeddings <- Embeddings(data_sel, 'harmony')

p1 <- DimPlot(object = data_sel, reduction = "harmony", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = data_sel, features = "harmony_1", group.by = "sample", pt.size = .1)
plot_grid(p1, p2)

data_sel <- data_sel %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = FALSE, n.neighbors = 30L, min.dist = 0.01) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20, verbose = FALSE) %>% 
  FindClusters(resolution = 1, verbose = FALSE) %>% 
  identity()

# Save UMAP plot for Harmony-integrated subclusters
q <- DimPlot(data_sel, pt.size = .5, reduction = "umap", cols = cols, label = TRUE, repel = TRUE) + optns + NoLegend() + NoAxes() + ggtitle("")
ggsave(filename = paste0("4_12_21_23", "/umap_harmony_RNA_noLegend.pdf"), width = 5, height = 5, q)

# Save Harmony-integrated subcluster data
saveRDS(data_sel, file = paste0("4_12_21_23", "/harmony.rds"))

# Marker gene analysis
print("Marker gene analysis...")
Idents(data_sel) <- data_sel$seurat_clusters
data_sel.markers <- FindAllMarkers(data_sel, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
write.table(data_sel.markers, file = paste0("4_12_21_23", "/markers_ch_removed.txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

# Subsetting data for visualization
data_sel$sample2 <- gsub("_rep*.", "", data_sel$sample)
desired_samples <- c("gt3a_kt57_9h", "gt3a_kt57_24h", "AvrRpt2_09h", "AvrRpt2_24h")
data_sel_sub6_forFIG <- subset(data_sel, sample2 %in% desired_samples)

q <- VlnPlot(data_sel_sub6_forFIG, features = "PUB36", group.by = "sample2", ncol = 1)
ggsave(filename = paste("4_12_21_23", "/_", paste("PUB36", collapse = "_"), "_clst6_violin_FIG.pdf", sep = ""), width = 5, height = 4, q)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##### Differential Expression Analysis
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Prepare for differential expression analysis
data_sel$clst.cond <- paste(Idents(data_sel), data_sel$sample2, sep = "_")
Idents(data_sel) <- "clst.cond"
data_sel <- PrepSCTFindMarkers(object = data_sel)

# Define functions for DE analysis
DE_Cluster_Genes_up <- function(a, b) {
  response_cluster <- FindMarkers(data_sel, ident.1 = a, ident.2 = b, verbose = FALSE)
  response_cluster$specificity <- response_cluster$pct.1 / response_cluster$pct.2
  response_significant <- subset(response_cluster, p_val_adj < 0.05 & avg_log2FC > 0.25)
  response_significant$gene <- rownames(response_significant)
  response_significant
}

DE_Cluster_Genes_down <- function(a, b) {
  response_cluster <- FindMarkers(data_sel, ident.1 = a, ident.2 = b, verbose = FALSE)
  response_cluster$specificity <- response_cluster$pct.1 / response_cluster$pct.2
  response_significant <- subset(response_cluster, p_val_adj < 0.05 & avg_log2FC < -0.25)
  response_significant$gene <- rownames(response_significant)
  response_significant
}

# Perform DE analysis for 9h
DE_Pairwise_up <- list(NA)
DE_Pairwise_down <- list(NA)
clusters <- unique(data_sel$SCT_snn_res.1) %>% sort()

for (i in seq_along(clusters)) {
  print(i)
  tryCatch({
    DE_Pairwise_up[[i]] <- DE_Cluster_Genes_up(paste(clusters[i], 'gt3a_kt57_9h', sep = '_'), paste(clusters[i], 'AvrRpt2_09h', sep = '_'))
    DE_Pairwise_down[[i]] <- DE_Cluster_Genes_down(paste(clusters[i], 'gt3a_kt57_9h', sep = '_'), paste(clusters[i], 'AvrRpt2_09h', sep = '_'))
    DE_Pairwise_up[[i]]$cluster <- clusters[i]
    DE_Pairwise_down[[i]]$cluster <- clusters[i]
  }, error = function(e) {})
}

write.table(do.call(rbind.data.frame, DE_Pairwise_up), file = "DEup_gt3a_kt57_9h_vs_AvrRpt2_09h.txt", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
write.table(do.call(rbind.data.frame, DE_Pairwise_down), file = "DEdown_gt3a_kt57_9h_vs_AvrRpt2_09h.txt", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

# Perform DE analysis for 24h
DE_Pairwise_up <- list(NA)
DE_Pairwise_down <- list(NA)

for (i in seq_along(clusters)) {
  print(i)
  tryCatch({
    DE_Pairwise_up[[i]] <- DE_Cluster_Genes_up(paste(clusters[i], 'gt3a_kt57_24h', sep = '_'), paste(clusters[i], 'AvrRpt2_24h', sep = '_'))
    DE_Pairwise_down[[i]] <- DE_Cluster_Genes_down(paste(clusters[i], 'gt3a_kt57_24h', sep = '_'), paste(clusters[i], 'AvrRpt2_24h', sep = '_'))
    DE_Pairwise_up[[i]]$cluster <- clusters[i]
    DE_Pairwise_down[[i]]$cluster <- clusters[i]
  }, error = function(e) {})
}

write.table(do.call(rbind.data.frame, DE_Pairwise_up), file = "DEup_gt3a_kt57_24h_vs_AvrRpt2_24h.txt", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
write.table(do.call(rbind.data.frame, DE_Pairwise_down), file = "DEdown_gt3a_kt57_24h_vs_AvrRpt2_24h.txt", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##### Analysis of PRIMER Cell DEGs
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Load DEGs for 9h and plot subclusters
up_9h <- read.delim("4_12_21_23/DEup_gt3a_kt57_9h_vs_AvrRpt2_09h.txt")
targets <- up_9h$gene[up_9h$cluster == 6]
targets <- up_9h$gene[up_9h$cluster == 12]


down_9h <- read.delim("4_12_21_23/DEdown_gt3a_kt57_9h_vs_AvrRpt2_09h.txt")
targets <- down_9h$gene[down_9h$cluster == 6] %>% na.omit() %>% sort()
targets <- down_9h$gene[down_9h$cluster == 12]

gene2 <- sapply(targets, genes_to_id)
ego <- enrichGO(
  gene = gene2,
  OrgDb = org.Athaliana.eg.db,
  ont = "BP",
  keyType = 'GID',
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE
)

pdf(paste("GO_mutantDOWN_9h_clst12", ".pdf", sep = ""))
try(print(barplot(ego, showCategory = 3)))
dev.off()
write.table(ego@result, file = paste0("GO_mutantDOWN_9h_clst12", ".txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

# Plot violin plot for ALD1
q <- VlnPlot(data_sel_sub12, features = "ALD1", group.by = "sample2", ncol = 1)
ggsave(filename = paste("4_12_21_23", "/_", paste("ALD1", collapse = "_"), "_clst12_violin_FIG.pdf", sep = ""), width = 5, height = 4, q)
