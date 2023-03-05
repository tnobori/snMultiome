source(".../_config_multiome.R")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Load annotation
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

annotation <- rtracklayer::import( ".../At_reference/Arabidopsis_thaliana.TAIR10.52.gff3")
colnames(annotation@elementMetadata)[colnames(annotation@elementMetadata)=="Name"] <- "gene_name" 
colnames(annotation@elementMetadata)[colnames(annotation@elementMetadata)=="biotype"] <- "gene_biotype" 
levels(annotation$type)[levels(annotation$type)=="gene"]="body"
annotation$gene_name[grepl("exon", annotation$gene_name)] <- gsub("\\.[0-9].exon[0-9]{1,}", "", annotation$gene_name[grepl("exon", annotation$gene_name)]) #remove exon number

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Load single-cell data
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

data <- readRDS("combined_filtered.rds")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## analysis
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

###each strain separately
cell_dc3000 <- rownames(data@meta.data)[ grepl("DC3000", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] #4100
cell_avrrpt2 <- rownames(data@meta.data)[ grepl("AvrRpt2", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] #4100
cell_avrrpm1 <- rownames(data@meta.data)[ grepl("AvrRpm1", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] #4100

data_dc3000 <- subset(x = data, cells = cell_dc3000)
data_avrrpt2 <- subset(x = data, cells = cell_avrrpt2)
data_avrrpm1 <- subset(x = data, cells = cell_avrrpm1)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Functions
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

##pseudobulk function
pb <- function(data_sel){
  DefaultAssay(data_sel) <- "RNA"
  clusters <- data_sel$SCT_snn_res.1 %>% unique() %>% sort(decreasing = F)
  out <- matrix(data = NA, nrow = nrow(data@assays$RNA), ncol = length(clusters))
  rownames(out) <- rownames(data_sel)
  colnames(out) <- clusters
  
  for (i in c(1:length(clusters))){
    data_sel.clst <- data_sel[, data_sel$SCT_snn_res.1==clusters[i]]
    pb <- data_sel.clst@assays$RNA@counts %>% rowSums()
    tpm <- pb/sum(pb)*1000000 
    tpm <- log2(tpm)
    tpm[is.infinite(tpm)] <-0
    out[, i] <- tpm
  }
  return(out)
}





kclust <- function(input, n){
  #input <- na.omit(input)
  km<- kmeans(input,n)
  m.kmeans<- cbind(input, km$cluster)
  colnames(m.kmeans)[ncol(m.kmeans)] <- "kmean"
  m.kmeans <- as.data.frame(m.kmeans)
  m.kmeans <- tibble::rownames_to_column(m.kmeans)
  o <-
    m.kmeans %>%
    arrange(kmean)
  
  #data_km <- o[, -ncol(o)]
  data_km <- column_to_rownames(o)
  data_km
}

coverage.plot2 <- function(data, target, size ,out.dir){
  data <- LinkPeaks(
    object = data,
    peak.assay = "peaks",
    expression.assay = "SCT",
    genes.use =target ##can be a single gene
  )
  q <- CoveragePlot(
    object = data,
    region = target,
    features = target,
    expression.assay = "SCT",
    extend.upstream = size,
    extend.downstream = size, 
    group.by = "sample.order",
    split.assays = TRUE, 
    heights = c(2,2,2,.5)
    
    
  )
  q <-q & scale_fill_manual(values = col.sample) #change colors 
  ggsave(filename = paste0(out.dir,"/coveragePlot_marker_", target, "_", size ,"bp.pdf"), width = 20, height = 10, q)
  
  q <- CoveragePlot(
    object = data,
    region = target,
    features = target,
    expression.assay = "SCT",
    extend.upstream = size,
    extend.downstream =size, 
    group.by = "SCT_snn_res.1",
    heights = c(2,2,2,.5)
  )
  q <-q & scale_fill_manual(values = cols) #change colors 
  ggsave(filename = paste0(out.dir,"/coveragePlot_marker_", target, "_", size ,"bp_RNAcluster.pdf"), width = 20, height = 10, q)
 
}

optns <- theme(
  axis.text.x = element_text(margin = margin(c(20, 0, 0, 0)), size = 30, face = "bold", family = "Helvetica"),
  axis.text.y = element_text(margin = margin(c(0, 20, 0, 0)) ,size = 30, face = "bold", family = "Helvetica"), 
  axis.title = element_text(size = 40, face = "bold", family = "Helvetica"),
  axis.title.x = element_text(margin = margin(c(30, 0, 0, 0))),
  axis.title.y = element_text(margin = margin(c(0, 30, 0, 0))),
  axis.ticks.y = element_line(size = 3),
  axis.ticks.x = element_line(size = 3),
  axis.ticks.length = unit(.5, "cm"),
  axis.line  = element_line(size = 2),
  legend.text = element_text(size = 15, face = "bold", family = "Helvetica", margin = margin(c(20, 0, 20, 0))),
  legend.title = element_text(size = 15, face = "bold", family = "Helvetica"),
  legend.key = element_rect(fill = NA , size= 2.5)) 

subclst.pb <- function(clst){
  print("subsetting data...")
  data_sel <- data[, data$SCT_snn_res.1 == clst ]
  
  print("SCT and PCA...")
  DefaultAssay(data_sel) <- "RNA"
  data_sel <- SCTransform(data_sel, verbose = F)
  data_sel <- RunPCA(data_sel, verbose = F)
  
  data_sel <- FindNeighbors(data_sel, dims = 1:15, verbose = F)
  data_sel <- FindClusters(data_sel, resolution = 1, verbose = FALSE)
  data_sel <- RunUMAP(data_sel, dims = 1:15, verbose = F)
  
  print("saving plot...")
  q <- DimPlot(data_sel, pt.size = 3, cols = cols, label = T, repel = T)+optns +NoLegend()+NoAxes()+ggtitle("")
  ggsave(filename = paste0(clst, "/umap_RNA.png"), width = 10, height = 10,q)
  
  q <- DimPlot(data_sel, reduction = "umap", pt.size = 4, split.by = 'sample', cols = cols, label = T, repel = T)+optns +NoLegend()+NoAxes()+ggtitle("")
  ggsave(filename = paste0(clst, "/umap_RNA2.png"), width = 40, height = 20, q)
  
  q <- DimPlot(data_sel, reduction = "umap", group.by = "sample", pt.size = 4, label = T, repel = T)+optns +NoLegend()+NoAxes()+ggtitle("")
  ggsave(filename = paste0(clst, "/umap_RNA3.png"), width = 20, height = 20, q)
  
  print("run harmony...")
  ################
  DefaultAssay(data_sel) <- "SCT"
  
  #visualize difference between samples 
  options(repr.plot.height = 5, repr.plot.width = 12)
  p1 <- DimPlot(object = data_sel, reduction = "pca", pt.size = .1, group.by = "sample")
  p2 <- VlnPlot(object = data_sel, features = "PC_1", group.by = "sample",  pt.size = .1)
  plot_grid(p1,p2)
  
  #Run Harmony 
  options(repr.plot.height = 2.5, repr.plot.width = 6)
  data_sel <- data_sel %>% 
    RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT", verbose = F)
  
  harmony_embeddings <- Embeddings(data_sel, 'harmony')

  options(repr.plot.height = 5, repr.plot.width = 12)
  p1 <- DimPlot(object = data_sel, reduction = "harmony", pt.size = .1, group.by = "sample") #removed " do.return = TRUE"
  p2 <- VlnPlot(object = data_sel, features = "harmony_1", group.by = "sample", pt.size = .1)
  plot_grid(p1,p2)
  
  
  data_sel <- data_sel %>% 
    RunUMAP(reduction = "harmony", dims = 1:20, verbose = F) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20, verbose = F) %>% 
    FindClusters(resolution = 1, verbose = F) %>% 
    identity()
  
  print("saving data...")
  saveRDS(data_sel, file = paste0(clst, "/harmony.rds"))
  
  print("marker gene analysis...")
  data_sel.markers <- FindAllMarkers(data_sel, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = F)
  out_sel_markers <- data_sel.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
  write.table(data_sel.markers, file = paste0(clst, "/markers.txt"), row.names=T, col.names=NA, sep="\t", quote=F)
  
  print("pseudobulking...")
  ##pseudobulk
  DefaultAssay(data_sel) <- "RNA"
  clusters <- data_sel$seurat_clusters %>% unique() %>% sort(decreasing = F)
  out <- matrix(data = NA, nrow = nrow(data@assays$RNA), ncol = length(clusters))
  rownames(out) <- rownames(data_sel)
  colnames(out) <- clusters
  
  for (i in c(1:length(clusters))){
    data_sel.clst <- data_sel[, data_sel$seurat_clusters==clusters[i]]
    pb <- data_sel.clst@assays$RNA@counts %>% rowSums()
    tpm <- pb/sum(pb)*1000000 
    tpm <- log2(tpm)
    tpm[is.infinite(tpm)] <-0
    out[, i] <- tpm
  }
  write.table(out, file = paste0(clst, "/pseudobulk"), row.names=T, col.names=NA, sep="\t", quote=F)
}

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## pseudobulking each cluster
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

clusters <- unique(data$SCT_snn_res.1) %>% sort()

#run subclst.pb function for each cluster 
for (i in c(1:length(clusters))){
  tryCatch({
    subclst.pb(clusters[i])
  }, error=function(e){})
}

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## sub cluster annotation
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

c <- data$SCT_snn_res.1 %>% unique() %>% sort()

data$sub_clst_rna <- NA
for (i in c(1:length(c))){
  tryCatch({
    data.sub <- readRDS(file = paste0( c[i], "/harmony.rds"))
    clst.comb <- paste0( c[i],"_" ,data.sub$seurat_clusters)
    data$sub_clst_rna[colnames(data.sub)] <- clst.comb
  }, error=function(e){})
}
