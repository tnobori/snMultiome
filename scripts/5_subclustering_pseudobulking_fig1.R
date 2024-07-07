##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## This script performs subclustering of each major cluster.
## Figures produced with this script: Fig1, Fig6, FigS3, Table S2

## A fully processed Seurat object can be downloaded at http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/processed_seurat_object/combined_filtered.rds

## Contact: Tatsuya Nobori (tatsuyanobori@gmail.com)
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Clearing the workspace
rm(list = ls())

# Loading project configuration
source("http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/scripts/_config_multiome.R")

## Assuming you are in the directory where this and other R scripts are stored.
mainDir <- file.path(getwd(), "data_out")
subDir <- "5_subclustering"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

setwd(file.path(mainDir, subDir))

# Load dataset
data <- readRDS(paste(mainDir, "_seurat_object","combined_filtered.rds", sep = "/"))

# Load annotation
annotation <- Annotation(data@assays$ATAC)

# Making an object for each strain
cell_kt56 <- rownames(data@meta.data)[ grepl("DC3000", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] 
cell_kt57 <- rownames(data@meta.data)[ grepl("AvrRpt2", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ]
cell_kt58 <- rownames(data@meta.data)[ grepl("AvrRpm1", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ]
cell_mock <- rownames(data@meta.data)[grepl("Mock", data@meta.data[, "sample"]) ]

data_kt56 <- subset(x = data, cells = cell_kt56)
data_kt57 <- subset(x = data, cells = cell_kt57)
data_kt58 <- subset(x = data, cells = cell_kt58)
data_mock <- subset(x = data, cells = cell_mock)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Functions
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Pseudobulking function
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
    #idents = "04h",
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
    #idents = idents.plot,
    extend.upstream = size,
    extend.downstream =size, 
    group.by = "SCT_snn_res.1",
    #split.assays = TRUE,
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
  legend.key = element_rect(fill = NA , size= 2.5)
) 

# Subclustering and pseudbulking 
subclst.pb <- function(clst){
 
  
  tic(paste0("cluster ", clst))
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
  
  dir.create(file.path(mainDir, subDir, clst), showWarnings = FALSE)
  
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
  # harmony_embeddings[1:5, 1:5]
  
  options(repr.plot.height = 5, repr.plot.width = 12)
  p1 <- DimPlot(object = data_sel, reduction = "harmony", pt.size = .1, group.by = "sample") #removed " do.return = TRUE"
  p2 <- VlnPlot(object = data_sel, features = "harmony_1", group.by = "sample", pt.size = .1)
  plot_grid(p1,p2)
  
  
  data_sel <- data_sel %>% 
    RunUMAP(reduction = "harmony", dims = 1:20, verbose = F, n.neighbors = 30L, min.dist = 0.01) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20, verbose = F) %>% 
    FindClusters(resolution = 1, verbose = F) %>% 
    identity()
  
  print("saving data...")
  saveRDS(data_sel, file = paste0(clst, "/harmony.rds"))
  
  
  q <- DimPlot(data_sel, pt.size = .5, reduction = "umap", cols = cols, label = T, repel = T)+optns +NoLegend()+NoAxes()+ggtitle("")
  ggsave(filename = paste0(clst, "/umap_harmony_RNA_noLegend.png"), width = 5, height = 5,q)
  
  print("marker gene analysis...")
  data_sel.markers <- FindAllMarkers(data_sel, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = F)
  out_sel_markers <- data_sel.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)

  write.table(data_sel.markers, file = paste0(clst, "/markers_ch_removed.txt"), row.names=T, col.names=NA, sep="\t", quote=F)

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
  toc()
}

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## pseudobulking each cluster; related to fig1b
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

clusters <- unique(data$SCT_snn_res.1) %>% sort()

#run subclst.pb function for each cluster 
for (i in c(1:length(clusters))){
  tryCatch({
    subclst.pb(clusters[i])
  }, error=function(e){})
}

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## sub cluster annotation; combine major and sub cluster labels to create new labels.
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

c <- data$SCT_snn_res.1 %>% unique() %>% sort()

data$sub_clst_rna <- NA
for (i in c(1:length(c))){
  tic(paste0("###cluster_", c[i]) )
  tryCatch({
    data.sub <- readRDS(file = paste0( c[i], "/harmony.rds"))
    clst.comb <- paste0( c[i],"_" ,data.sub$seurat_clusters)
    data$sub_clst_rna[colnames(data.sub)] <- clst.comb
  }, error=function(e){})
  toc()
}

out.dir.seurat <- file.path(mainDir, "/_seurat_object")
saveRDS(data, paste0(out.dir.seurat, "/combined_filtered.rds"))


##make subcluster pseudobulk table

sub.comb <- read.delim(file = paste0("0", "/pseudobulk"),header = T, row.names = 1)
x <- ncol(sub.comb) - 1
colnames(sub.comb) <- paste("0", c(0:x), sep = "_")

clusters <- data$seurat_clusters %>% unique() %>% sort(decreasing = F)
for (i in c(2:length(clusters))){
  tryCatch({pb <- read.delim(file = paste0(clusters[i], "/pseudobulk"), header = T, row.names = 1)
  x <- ncol(pb) - 1
  colnames(pb) <- paste(clusters[i], c(0:x), sep = "_")
  sub.comb <- cbind(sub.comb, pb)}, error=function(e){})
}

data.comb <-sub.comb

cor.mat <- cor(sub.comb)

annotation_col_sub <- colnames(cor.mat) %>% as.data.frame()
rownames(annotation_col_sub) <- annotation_col_sub$.
annotation_col_sub$. <- gsub("_.*", "",annotation_col_sub$. )

list_color <- list(. = cols)

pdf("figS3a.pdf")
pheatmap(cor.mat,
         #scale = "row",
         cluster_cols=T,
         cluster_rows = T,
         #cellwidth = 70,
         #cellheight = .07,
         #color = colorRampPalette(c("magenta", "black",  "green"))(100),
         #color = colorRampPalette(brewer.pal(n = 11, name ="RdBu"))(100),
         color = hcl.colors(50, "Cividis"),
         breaks= seq(0.75, 1, length.out = 50),
          annotation_row = annotation_col_sub,
          annotation_col = annotation_col_sub,
         annotation_colors = list_color,
         fontsize_row = .00001,
         fontsize_col = .00001,
         
)
dev.off()

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
###heatmap analysis using marker genes
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

mainDir <-  "/data_out"
subDir <- "5_subclustering/heatmap_using_MajSub_markers"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

setwd(file.path(mainDir, subDir))

#####assemble sub cluster markers 
out <- c()
clusters <- data$SCT_snn_res.1 %>% unique() %>% sort()
for (i in c(1:length(clusters))){
  marker.sub <- read.delim(file = paste0(mainDir, "/5_subclustering/", clusters[i], "/markers_ch_removed.txt"), header = T, row.names = 1)
  marker.sub <- marker.sub[marker.sub$avg_log2FC > 0.5 ,]
  out <- c(out, marker.sub$gene)
}
out <- out %>% unique()

write.table(out, file= "subclst_markers.txt", row.names=T, col.names=NA, sep="\t", quote=F)

major.marker <- read.delim("/data_out/2_clustering/RNA/markers_RNA.txt")
major.marker <- major.marker$gene %>% unique()

markers_all <- c(out$x, major.marker) %>% unique()

input <- sub.comb[markers_all, ] %>% na.omit()

input <- input[!rowSums(input==0)==ncol(input) ,]


##kmeans clustering 
message("kmean clustering")

pdf("kplot_genes.pdf")
kplot(input)
dev.off()

n.cluster.gene <- 12 #how many clusters?

input_clst <- kclust(input, n.cluster.gene)
an_row_km <- input_clst[, ncol(input_clst)]
names(an_row_km) <- rownames(input_clst)
an_row_km <- as.data.frame(an_row_km)
an_row_km[, 1] <- as.character(an_row_km[, 1])

input_clst2 <- input_clst[, -ncol(input_clst)]

write.table(input_clst, file=paste( "kmean_clustering_sub.txt", sep = ""), row.names=T, col.names=NA, sep="\t", quote=F)

input_clst <- read.delim("kmean_clustering_sub.txt", header = T, row.names = 1)
an_row_km <- input_clst$kmean %>% as.character() %>% as.data.frame()
rownames(an_row_km) <- rownames(input_clst)
colnames(an_row_km) <- "kmean"

input_clst2 <- input_clst2[!rowSums(input_clst2==0)==ncol(input_clst2), ]
colnames(input_clst2) <- gsub("X", "", colnames(input_clst2))
pdf(paste("heatmap_subclst_kmean2", ".pdf", sep = ""), width = 10, height = 4)

pheatmap(input_clst2,
         scale = "row",
         cluster_cols=TRUE,
         cluster_rows = FALSE,
         color = rev(colorRampPalette(brewer.pal(n = 11, name ="RdYlBu"))(100)),
         breaks= seq(-2, 2, length.out = 101),
         annotation_row = an_row_km,
         annotation_col = ann_col,
         annotation_colors = list_color,
         fontsize_row = .00001,
         fontsize_col = .00001,
         legend = T,
         treeheight_col = 0
)
dev.off()



genes_to_id <- function(genes){
  if (length(grep(genes, annotation$ID))==0 ){
    tryCatch({id <- annotation$ID[annotation$gene_name==genes] %>% na.omit
    id <- gsub("gene:", "", id)
    id[1] }, error=function(e){})
  }else{
    genes
  }
}
an_row_km$geneID <- sapply(rownames(an_row_km), genes_to_id) #gene names to ID

x <- c()
for (i in c(1:nrow(an_row_km))){
  x <- c(x, genes_to_id(rownames(an_row_km)[i]))
  if (mod(i, 10)==0){
    print(i)
  }
}

an_row_km$geneID <- x
gene_master_list <- an_row_km #use this for future gene ID conversion for this analysis

for (i in c(1:n.cluster.gene)){
  gene <- an_row_km$geneID[an_row_km[ ,1]==i]
  
  ego <- enrichGO(gene          = gene,
                  #universe      = names(geneList),
                  OrgDb         = org.Athaliana.eg.db,
                  ont           = "BP",
                  keyType = 'GID',
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  png(paste("GO_clst_", i,".png", sep = ""), 
      # width = 1000, height = 1000
  )
  try(print(barplot(ego, showCategory=10)))
  dev.off()
  
  write.table(ego@result, file= paste0("go_clst_", i, ".txt"), row.names=T, col.names=NA, sep="\t", quote=F)
  
}

#GO enrichment summary 

out <- matrix(data = NA, nrow =n.cluster.gene, ncol = 3) %>% as.data.frame()
colnames(out) <- c("Description", "Size", "adjusted pvalue")


optns3 <- theme(
  axis.text.x = element_text(margin = margin(c(10, 0, 0, 0)), size = 15, family = "Helvetica", angle = 0),
  axis.text.y = element_text(margin = margin(c(0,10, 0, 0)) ,size = 15, family = "Helvetica"), 
  axis.title = element_text(size = 10,  family = "Helvetica"),
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

for (i in c(1:n.cluster.gene)){
  go.data <- read.delim(file = paste0("GO_clst_",i,".txt"), header = T, row.names = 1)
  
  
  
  out[i, ] <- c(go.data$Description[1],as.numeric(go.data$Count[1]), as.numeric(go.data$p.adjust[1]))
  
}

out$clst <- as.character(c(1:n.cluster.gene))
#Then turn it back into a factor with the levels in the correct order
out$clst <- factor(out$clst, levels=unique(out$clst))

out$log.p <- -log10(as.numeric(out$`adjusted pvalue`))
out$Size <- as.numeric(out$Size)

q <- ggplot(out, aes(x = clst, y = Description, size = Size , color = log.p)) + geom_point()+ optns3 + 
  scale_color_gradient(low = "grey",high = "#418557",limits = c(0,51))+
  scale_size(range = c(0,10)) 
ggsave("go_summary_plot_for_paper.pdf", height = 6, width = 8 ,q)


# find immune cluster based on the GO analysis and use these genes as immune genes
immune.genes <- rownames(an_row_km)[an_row_km$kmean==1|an_row_km$kmean==2|an_row_km$kmean==10]

#z score
input <- input[immune.genes, ]
mean <- apply(input, 1, mean)
sd <- apply(input, 1, sd)
z <- (input - mean)/sd


ann_col <- gsub("_[0-9]+","", colnames(input_clst2)) %>% as.data.frame()
rownames(ann_col) <- colnames(input_clst2)

tissues <- c("#418557","#af2b3d","#ff9f22","#2c2cff" )
names(tissues) <- c("Mesophyll", "Epidermis", "Vasculature", "Unknown")
library(rstatix)
x <- data@meta.data[c("SCT_snn_res.1", "celltype")]
x$color <- NA
x$color[x$celltype=="Mesophyll"] <- "#418557"
x$color[x$celltype=="Epidermis"] <- "#af2b3d"
x$color[x$celltype=="Vasculature"] <- "#ff9f22"
x$color[x$celltype=="Unknown"] <- "#2c2cff"
y <- paste(x[,1], x[,2], x[,3], sep = "_") %>% unique()
DF <- data.frame(do.call(rbind, strsplit(y, "_", fixed=TRUE)))
rownames(DF) <- DF$X1

ann_col$celltype <- NA
ann_col$celltype <- DF[ann_col$. , 2]

celltype_color <- list("celltype" = c("Mesophyll" = "#418557", 
                                      "Epidermis" = "#af2b3d", 
                                      "Vasculature" = "#ff9f22", 
                                      "Unknown" = "#2c2cff"))
list_color <- c(list_color, celltype_color)

pdf(paste("figS3c", ".pdf", sep = ""), width = 20, height = 8)
pheatmap(z,
         cluster_cols=TRUE,
         cluster_rows = TRUE,
         color = rev(colorRampPalette(brewer.pal(n = 11, name ="RdYlBu"))(100)),
         breaks= seq(-2, 2, length.out = 101),
         annotation_col = ann_col,
         annotation_colors = list_color,
         fontsize_row = .00001,
         fontsize_col = .00001,
         legend = F,
         annotation_legend = F,
         treeheight_row = 0,
         treeheight_col = 0
)
dev.off()

#select genes for fig1d

gene_select <- c(
  "BGLU46", "EDS1", "ALD1", "CYP81F1", "PR1",
  "FMO1", "CYP79B2", "RBOHF", "UGT74B1",
  "EDS16", "RBOHD", "FRK1", "PBS3", "MPK3", "MPK6",
  "EIN2", "MYC2", "VSP1"
)


#z score
input <- sub.comb[gene_select, ] %>% na.omit()
mean <- apply(input, 1, mean)
sd <- apply(input, 1, sd)
z <- (input - mean)/sd

pdf(paste("fig1d", ".pdf", sep = ""), width = 8, height = 5)
pheatmap(z,
         cluster_cols=TRUE,
         cluster_rows = TRUE,
         color = rev(colorRampPalette(brewer.pal(n = 11, name ="RdYlBu"))(100)),
         breaks= seq(-2, 2, length.out = 101),
         annotation_col = ann_col,
         annotation_colors = list_color,
         fontsize_row = 8,
         fontsize_col = 0.0001,
         legend = F,
         annotation_legend = F,
         treeheight_row = 0,
         treeheight_col = 0
)
dev.off()

# Function for generating feature plots for each strain
plotting <- function(gene_list, fname){
  q <- FeaturePlot(data_kt56, features = gene_list, split.by = "sample", pt.size = 2, order=TRUE, max.cutoff = "q99", min.cutoff = "q1", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste("RNA_featurePlot/","_", fname,"_feture_KT56.pdf", sep = ""), width = 20, height = 5, q)
  
  
  q <- FeaturePlot(data_kt57, features = gene_list, split.by = "sample", pt.size = 2, order=TRUE, max.cutoff = "q99", min.cutoff = "q1", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste( "RNA_featurePlot/","_", fname,"_feture_KT57.pdf", sep = ""), width = 20, height = 5, q)
  
  
  q <- FeaturePlot(data_kt58, features = gene_list, split.by = "sample", pt.size = 2, order=TRUE, max.cutoff = "q99", min.cutoff = "q1", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste( "RNA_featurePlot/","_", fname,"_feture_KT58.pdf", sep = ""), width = 20, height = 5, q)
  
  q <- FeaturePlot(data, features = gene_list, order=TRUE, pt.size = 2, max.cutoff = "q99", min.cutoff = "q1", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste("RNA_featurePlot/","_", fname,".pdf", sep = ""), width = 10, height = 10, q)
  
}



##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
#####sub clustering immune active clusters (mesophyll): clusters 3, 7, 11
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

data_sel <- data[, data$SCT_snn_res.1 %in% c(3,7,11)  ]
dir.create("3_7_11", showWarnings = FALSE)

print("SCT and PCA...")
DefaultAssay(data_sel) <- "RNA"
data_sel <- SCTransform(data_sel, verbose = F)
data_sel <- RunPCA(data_sel, verbose = F)

data_sel <- FindNeighbors(data_sel, dims = 1:15, verbose = F)
data_sel <- FindClusters(data_sel, resolution = 1, verbose = FALSE)
data_sel <- RunUMAP(data_sel, dims = 1:15, verbose = F)

print("saving plot...")
q <- DimPlot(data_sel, pt.size = 3, cols = cols, label = T, repel = T)+optns +NoLegend()+NoAxes()+ggtitle("")
ggsave(filename = paste0("3_7_11", "/umap_RNA.png"), width = 10, height = 10,q)

q <- DimPlot(data_sel, reduction = "umap", pt.size = 4, split.by = 'sample', cols = cols, label = T, repel = T)+optns +NoLegend()+NoAxes()+ggtitle("")
ggsave(filename = paste0("3_7_11", "/umap_RNA2.png"), width = 40, height = 20, q)

q <- DimPlot(data_sel, reduction = "umap", group.by = "sample", pt.size = 4, label = T, repel = T)+optns +NoLegend()+NoAxes()+ggtitle("")
ggsave(filename = paste0("3_7_11", "/umap_RNA3.png"), width = 20, height = 20, q)

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
# harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = data_sel, reduction = "harmony", pt.size = .1, group.by = "sample") #removed " do.return = TRUE"
p2 <- VlnPlot(object = data_sel, features = "harmony_1", group.by = "sample", pt.size = .1)
plot_grid(p1,p2)


data_sel <- data_sel %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = F, n.neighbors = 30L, min.dist = 0.01) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20, verbose = F) %>% 
  FindClusters(resolution = 1, verbose = F) %>% 
  identity()

print("saving data...")
saveRDS(data_sel, file = paste0("3_7_11", "/harmony.rds"))
data_sel <- readRDS(paste0("3_7_11", "/harmony.rds"))
#data_sel <- readRDS(file = "RNA_featurePlot_clst34101113/harmony.rds")

q <- DimPlot(data_sel, pt.size = .5, reduction = "umap", cols = cols, label = T, repel = T)+optns +NoLegend()+NoAxes()+ggtitle("")
ggsave(filename = paste0("3_7_11", "/umap_harmony_RNA_noLegend.pdf"), width = 5, height = 5,q)

q <- DimPlot(data_sel, pt.size = .5, reduction = "umap", cols = tol22rainbow, label = F, repel = T)+optns +NoLegend()+NoAxes()+ggtitle("")
ggsave(filename = paste0("3_7_11", "/umap_harmony_RNA_noLegend_noLabel.pdf"), width = 5, height = 5,q)


q <- DimPlot(data_sel, reduction = "umap", pt.size = 1, split.by = 'sample', cols = cols, label = T, repel = T)+optns +NoLegend()+NoAxes()+ggtitle("")+optns
ggsave(filename = paste0("3_7_11", "/umap_harmony_RNA2.pdf"), width = 40, height = 6, q)

q <- DimPlot(data_sel, reduction = "umap", group.by = "sample", pt.size = 4,  label = T, repel = T)+optns +NoLegend()+NoAxes()+ggtitle("")+optns
ggsave(filename =paste0("3_7_11", "/umap_harmony_RNA3.pdf"), width = 20, height = 20, q)

print("marker gene analysis...")
data_sel.markers <- FindAllMarkers(data_sel, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = F)
out_sel_markers <- data_sel.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
# 
write.table(data_sel.markers, file = paste0("3_7_11", "/markers_ch_removed.txt"), row.names=T, col.names=NA, sep="\t", quote=F)
out_sel_markers <- read.delim(file = paste0("3_7_11", "/markers_ch_removed.txt"))

plot_subclst <- function(targets){
  q <- FeaturePlot(data_sel, features = targets, order=TRUE, pt.size = 1, max.cutoff = 'q99', min.cutoff = 'q1')
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste("3_7_11","/_", paste(targets, collapse  = "_"),".pdf", sep = ""), width = 10, height = 10, q)
}

targets <- c("CYP71A13", "CYP79B2", "CYP79B3", "IGMT2", "CYP81F2", "BGLU26", "ABCG36", "CYP71B15" )
plot_subclst(targets) #fig1g
targets <- c("ALD1","FMO1","BON3", "GT-3A")
plot_subclst(targets) #fig6ef

targets <- c("WRKY8","LSD1")
plot_subclst(targets) #fig8b



## motif analysis 
dir.create("3_7_11/marker_motif_plot", showWarnings = FALSE)

motif.activity <- function(gene_list, fname){
  
  input_motif <- gene_list
  ##################################
  DefaultAssay(data_sel) <- 'chromvar'  
  
  tryCatch({
    q1 <- FeaturePlot(
      object = data_sel,
      features = input_motif,
      min.cutoff = 'q10',
      max.cutoff = 'q99',
      pt.size = 1,
      
      order = TRUE,
      cols = c("lightgrey","blue"),
      reduction = "umap"
    )
    q1 <- q1 & scale_colour_viridis_c()& NoAxes()
    
    ggsave(filename = paste("3_7_11/marker_motif_plot/", "_", fname,"_feture.pdf", sep = ""), width = 5, height = 5, q1) 
  }, error=function(e){})
 
}

# fig6h
motif.activity(gene_list = "MA1207.1_GT-A3", 
               fname = "MA1207.1_GT-A3")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##PCC (Phloem Companion Cells) analysis; fig1e related
# NOTE: In our analysis, cluster 6 was annotated as PCC.
# This may be different if you re-run clustering. 
# Carefully inspect PCC marker gene expression before moving forward.
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
pcc_cluster <- 6

data_sel <- readRDS(paste0(pcc_cluster, "/harmony.rds"))

#motif enrichment 
Idents(data_sel) <- data_sel$seurat_clusters
DefaultAssay(data_sel) <- "peaks"
da_peaks <- FindAllMarkers(
  object = data_sel,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

clusters <- da_peaks$cluster %>% unique()

pcc.marker <- read.delim(file.path(pcc_cluster, "markers_ch_removed.txt"))
pcc.marker_sel <- pcc.marker[pcc.marker$cluster==8, ] #sub cluster 8 is enriched with SAR genes 
pcc.marker_sel <- pcc.marker_sel[pcc.marker_sel$avg_log2FC > 0.5, ]

major.markers <- read.delim(file.path(mainDir, "2_clustering/RNA/markers_RNA.txt"))
major.markers.immune <- major.markers[major.markers$cluster %in% c(3,4, 7,11,12,24,29) ,]
major.markers.immune <- major.markers.immune[major.markers.immune$avg_log2FC > 0.8, ]
major.markers.gene <- major.markers.immune$gene %>% unique()

# PCC subcluster 8 markers that are not markers for other immune active major clusters
pcc.unique <- pcc.marker_sel$gene[!(pcc.marker_sel$gene %in% major.markers.gene)]

plot_subclst <- function(targets){
  q <- FeaturePlot(data_sel, features = targets, order=TRUE, pt.size = 1, max.cutoff = 'q99', min.cutoff = 'q1')
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste(pcc_cluster,"/_", paste(targets, collapse  = "_"),".pdf", sep = ""), width = 10, height = 10, q)
}

DefaultAssay(data_sel) <- "SCT"
targets <- c("EDS16","SARD1","ALD1", "PAD4", "FMO1", "ILL6")
plot_subclst(targets) #fig1e


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
#####sub-cluster marker genes; Table S2
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

out <- read.delim(paste0( "0", "/markers_ch_removed.txt"))
out$cluster <- paste("0", out$cluster, sep = "_")

for (clst in c(1:29)){
  marker_table <- read.delim(paste0( clst, "/markers_ch_removed.txt"))
  marker_table$cluster <- paste(clst, marker_table$cluster, sep = "_")
  
  out <- rbind(out, marker_table)
}

out <- out[, -1]
write.table(out, file=paste0("subcluster_marker_combined.txt"), row.names=T, col.names=NA, sep="\t", quote=F)




##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##### Comparisons between PRIMER cells and bystander cells, Figure 6i
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Load all peaks from the data
allpeaks <- data@assays$peaks@ranges

# Filter genes by annotation type "body"
genes <- annotation[annotation$type == "body", ]

# Calculate distances between each peak and its nearest gene
distances_2 <- distanceToNearest(allpeaks, genes)

# Set threshold for distance
threshold <- 2000
withinThreshold_2 <- distances_2[distances_2@elementMetadata$distance <= threshold]

# Get peak-gene pairs within the threshold
genesWithinThreshold_2 <- genes[withinThreshold_2@to]
genesWithinThreshold_2_comb <- genesWithinThreshold_2$gene_name
genesWithinThreshold_2_comb[is.na(genesWithinThreshold_2$gene_name)] <- genesWithinThreshold_2$gene_id[is.na(genesWithinThreshold_2$gene_name)]
peaksWithinThreshold_2 <- allpeaks[withinThreshold_2@from]

# Function to perform motif enrichment analysis
motif_enrich_from_genes <- function(input_genes, fname){
  peak_match <- peaksWithinThreshold_2[genesWithinThreshold_2_comb %in% input_genes]
  
  # Convert peak match to string format
  peak_match_string <- mapply(function(seq, start, end) {
    paste(seq, start, end, sep = "-")
  }, as.character(seqnames(peak_match)), start(peak_match), end(peak_match))
  
  peak_match_string <- as.character(peak_match_string)
  
  # Perform motif enrichment analysis
  enriched.motifs <- FindMotifs(
    object = data,
    features = peak_match_string
  )
  
  # Save results
  dir.create("./motif_enrich_from_gene_set", showWarnings = FALSE)
  write.table(enriched.motifs, file = paste0("motif_enrich_from_gene_set/", fname, ".txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  
  # Plot the top enriched motifs
  q <- MotifPlot(
    object = data,
    motifs = head(rownames(enriched.motifs), 24)
  )
  ggsave(paste0("motif_enrich_from_gene_set/", fname, ".pdf"), height = 10, width = 10, q)
}

# Load the dataset
data_sel <- readRDS(paste0("3_7_11", "/harmony.rds"))

# Define cell states based on clusters
data_sel$cellstate <- NA
data_sel$cellstate[data_sel$seurat_clusters == 4] <- "PRIMER" # PRIMER cell cluster
data_sel$cellstate[data_sel$seurat_clusters == 14 | data_sel$seurat_clusters == 18 | data_sel$seurat_clusters == 2] <- "bystander" # Bystander cell clusters

# Set cell identities and find marker genes
data_sel2 <- data_sel
Idents(data_sel2) <- data_sel2$cellstate
bystander_primer_comp_marker <- FindMarkers(data_sel2, ident.1 = "PRIMER", ident.2 = "bystander")

# Define upregulated genes
primer_up <- rownames(bystander_primer_comp_marker)[bystander_primer_comp_marker$avg_log2FC > 0.5]
bystander_up <- rownames(bystander_primer_comp_marker)[bystander_primer_comp_marker$avg_log2FC < -0.5]

# Perform motif enrichment analysis
motif_enrich_from_genes(primer_up, "PRIMER_up_vs_bystander") # Figure 6i
motif_enrich_from_genes(bystander_up, "bystander_up_vs_PRIMER") # Figure 6i




##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##### Generalized Function to Find Genes with Nearby Motifs
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

motifGenes <- function(input_motif, peak_distance = 2000, full = FALSE) {
  # Create and set working directory for output
  dir.create(file.path(mainDir, subDir, paste(input_motif, peak_distance, "bp", sep = "_")), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir, paste(input_motif, peak_distance, "bp", sep = "_")))
  
  # Select motifs and convert to GRanges
  motif.sel <- data@assays$peaks@motifs@positions[grep(input_motif, data@assays$peaks@motifs@motif.names)] %>% as.data.frame()
  motif.sel.ranges <- paste(motif.sel$seqnames, motif.sel$start, motif.sel$end, sep = "-")
  ranges.show <- StringToGRanges(motif.sel.ranges)
  
  # Get all peaks and find overlaps with selected motifs
  allpeaks <- data@assays$peaks@ranges
  overlaps <- subsetByOverlaps(allpeaks, ranges.show)
  
  # Filter genes by annotation type "body"
  genes <- annotation[annotation$type == "body", ]
  
  # Calculate distances and filter by threshold
  distances <- distanceToNearest(overlaps, genes)
  threshold <- peak_distance
  withinThreshold <- distances[distances@elementMetadata$distance <= threshold]
  
  # Get genes within the threshold distance
  genesWithinThreshold <- genes[withinThreshold@to]
  genesWithinThreshold_genes <- c(genesWithinThreshold$gene_name[!is.na(genesWithinThreshold$gene_name)], genesWithinThreshold$gene_id[is.na(genesWithinThreshold$gene_name)])
  
  # Save genes with the selected motif within the threshold distance
  write.table(genesWithinThreshold, file = paste0("genes_with_", input_motif, "_motif_within_", threshold, "bp.txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  write.table(sort(table(genesWithinThreshold_genes), decreasing = TRUE), file = paste0(input_motif, "_motif_containing_genes_with_N_motifs.txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  
}

# Run the motifGenes function for specified motifs
motifGenes("CMTA3", peak_distance = 5000, full = TRUE)

                                                      
                                                      

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##### Analysis of bulk RNA-seq data of CAMTA3 mutant; Figure 6j
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Load bulk RNA-seq data
bulk <- read.delim(file = "/path/to/GSE92702_CountTable_At_2278_raw.txt", row.names = 1) # Download the raw count file from GEO (Jacob et al., PMID: 29226970)

# Normalize and log-transform the data
bulk_norm <- sweep(bulk, 2, colSums(bulk), "/") * 1000000
bulk_norm_log <- log2(bulk_norm)
bulk_norm_log[bulk_norm_log == -Inf] <- 0

# Average replicates
bulk_norm_log_ave <- cbind(
  ave(bulk_norm_log[, 1], bulk_norm_log[, 2]),
  ave(bulk_norm_log[, 3], bulk_norm_log[, 4]),
  ave(bulk_norm_log[, 5], bulk_norm_log[, 6]),
  ave(bulk_norm_log[, 7], bulk_norm_log[, 8]),
  ave(bulk_norm_log[, 9], bulk_norm_log[, 10]),
  ave(bulk_norm_log[, 11], bulk_norm_log[, 12])
)
rownames(bulk_norm_log_ave) <- rownames(bulk_norm_log)
colnames(bulk_norm_log_ave) <- gsub("_[0-9]", "", colnames(bulk_norm_log)) %>% unique()

# Calculate fold changes
fc_wt_vs_camta3D <- cbind(
  bulk_norm_log_ave[, 2] - bulk_norm_log_ave[, 1],
  bulk_norm_log_ave[, 6] - bulk_norm_log_ave[, 3],
  bulk_norm_log_ave[, 5] - bulk_norm_log_ave[, 4]
)
colnames(fc_wt_vs_camta3D) <- c("Rpm1", "Rps4", "Flg22")

# Identify genes suppressed in ETI
suppressed_eti <- rownames(fc_wt_vs_camta3D)[fc_wt_vs_camta3D[, 1] < -1 & fc_wt_vs_camta3D[, 2] < -1]

# Save suppressed ETI genes
write.table(suppressed_eti, file = "camta3d_suppressed_in_ETI.txt", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

# Function to convert gene IDs to gene names
id_to_gene <- function(genes) {
  if (length(na.omit(annotation$gene_name[annotation$gene_id == genes])) != 0) {
    tryCatch({id <- na.omit(annotation$gene_name[annotation$gene_id == genes])
    id[1]}, error = function(e) {})
  } else {
    genes
  }
}

# Convert suppressed ETI gene IDs to gene names
suppressed_eti <- lapply(suppressed_eti, id_to_gene)
suppressed_eti <- unlist(suppressed_eti)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##### Overlap with Cell State Marker Genes
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Load subcluster marker genes
submarker <- read.delim("3_7_11/markers_ch_removed.txt")
primer <- submarker$gene[submarker$cluster == 4] # PRIMER cell genes
bystander <- unique(submarker$gene[submarker$cluster %in% c(18, 14, 2)]) # Bystander cell genes

# Define specific marker genes
bystander_specific <- sort(bystander[!bystander %in% primer])
primer_specific <- primer[!primer %in% bystander]

# Identify overlap with suppressed ETI genes
suppressed_eti_primer <- suppressed_eti[suppressed_eti %in% primer]
suppressed_eti_primer_specific <- suppressed_eti[suppressed_eti %in% primer_specific]
suppressed_eti_bystander <- suppressed_eti[suppressed_eti %in% bystander]
suppressed_eti_bystander_specific <- suppressed_eti[suppressed_eti %in% bystander_specific]

# Load CAMTA motif genes
camta_motif_genes <- read.delim("CMTA3_5000_bp/CMTA3_motif_containing_genes_with_N_motifs.txt", row.names = 1) #out put from motifGenes("CMTA3", peak_distance = 5000, full = TRUE)
camta_motif_genes <- camta_motif_genes$genesWithinThreshold_genes

# Find overlap with CAMTA motif genes
suppressed_eti_bystander_specific[suppressed_eti_bystander_specific %in% camta_motif_genes]

# Save results
write.table(suppressed_eti_primer, file = "genes_suppressed_by_CAMTA3D_PRIMER.txt", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
write.table(suppressed_eti_bystander_specific, file = "genes_suppressed_by_CAMTA3D_bystander_specific.txt", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

## Using DEGs between PRIMER and bystander
primer_up <- rownames(bystander_primer_comp_marker)[bystander_primer_comp_marker$avg_log2FC > 0.5]
bystander_up <- rownames(bystander_primer_comp_marker)[bystander_primer_comp_marker$avg_log2FC < -0.5]

suppressed_eti_primer <- suppressed_eti[suppressed_eti %in% primer_up]
suppressed_eti_bystander <- suppressed_eti[suppressed_eti %in% bystander_up]

## Venn Diagram
library(VennDiagram)

# Define the lengths for Venn diagram
primer_up <- 4180
bystander_up <- 2636
suppressed_eti_primer <- 32
suppressed_eti_bystander <- 123
suppressed_eti <- 393

# Create and save the Venn diagram
venn.plot <- draw.triple.venn(
  area1 = primer_up,
  area2 = suppressed_eti,
  area3 = bystander_up,
  n12 = suppressed_eti_primer,
  n23 = suppressed_eti_bystander,
  n13 = 0,
  n123 = 0,
  category = c("PRIMER", "CAMTA3 suppressed", "Bystander"),
  fill = c("skyblue", "pink1", "mediumorchid"),
  cex = 2,
  cat.cex = 2,
  cat.col = c("skyblue", "pink1", "mediumorchid")
)

png("venn_diagram.png")
grid.draw(venn.plot)
dev.off()

## Euler Diagram
install.packages("eulerr")
library(eulerr)

# Define the lengths and overlaps for Euler diagram
fit <- euler(c(
  "PRIMER" = 4180,
  "Bystander" = 2636,
  "Suppressed by CAMTA3" = 393,
  "PRIMER&Bystander" = 0,
  "PRIMER&Suppressed by CAMTA3" = 32,
  "Bystander&Suppressed by CAMTA3" = 123,
  "PRIMER&Bystander&Suppressed by CAMTA3" = 0
))

# Plot and save the Euler diagram
pdf("venn_diagram.pdf", width = 3, height = 3)
plot(fit, 
     fills = list(fill = c("lightskyblue", "lightcoral", "thistle")),
     labels = list(font = 2, cex = 0.3),
     edges = list(lwd = 0.5),
     quantities = list(font = 2, cex = 0.3)
)
dev.off()

## Hypergeometric Test
q <- 123  # number of white balls drawn
m <- 393  # number of white balls in the box
n <- 35560  # number of black balls in the box
k <- 2636  # number of balls drawn

# Calculate p-value and fold enrichment
pval <- phyper(q, m, n, k, log.p = FALSE, lower.tail = FALSE)
fold_enrich <- (q / k) / (m / (m + n))
fdr <- p.adjust(pval, method = "fdr")

q <- 32  # number of white balls drawn
m <- 393  # number of white balls in the box
n <- 35560  # number of black balls in the box
k <- 4180  # number of balls drawn

# Calculate p-value and fold enrichment
pval <- phyper(q, m, n, k, log.p = FALSE, lower.tail = FALSE)
fold_enrich <- (q / k) / (m / (m + n))
fdr <- p.adjust(pval, method = "fdr")
