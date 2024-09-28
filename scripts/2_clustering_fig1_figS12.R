##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## This script performs 1st round clustering (major clustering) based on snRNA-seq, snATAC-seq, or both.
## Figures produced with this script: Fig1, FigS1, FigS2, Table S1

## A fully processed Seurat object can be downloaded at http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/processed_seurat_object/combined_filtered.rds

## Contact: Tatsuya Nobori (tatsuyanobori@gmail.com)
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

rm(list = ls())
source("http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/scripts/_config_multiome.R")

## Assuming you are in the directory where this and other R scripts are stored.
mainDir <- file.path(getwd(), "data_out")
subDir <- "2_clustering/RNA"

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
## RNA clustering
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
data <- readRDS(paste(mainDir, "_seurat_object","combined_filtered.rds", sep = "/")) ## A fully processed Seurat object can be downloaded at http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/processed_seurat_object/combined_filtered.rds. Or you can create one from scratch by running 1_qc_data_integration_figS1.R

# Load annotation
annotation <- Annotation(data@assays$ATAC)

#Gene expression data processing
DefaultAssay(data) <- "RNA"
data <- SCTransform(data, vars.to.regress = "percent.mt")
data <- RunPCA(data)

# checking PCs
DepthCor(data, reduction = "pca" )
ElbowPlot(data, reduction = "pca")

# setting parameters
nn = 20
mindist = 0.01

data <- FindNeighbors(data, dims = 1:20)
data <- FindClusters(data, resolution = 1, verbose = TRUE)
data <- RunUMAP(data, dims = 1:20,
                n.neighbors = nn, 
                min.dist = mindist
                )

q <- DimPlot(data, pt.size = 2, label = T)+optns +NoLegend()+NoAxes()+ggtitle("")
ggsave(filename = "umap_RNA_noHarmony.png", width = 10, height = 10,q)

# Harmony integration
DefaultAssay(data) <- "SCT"

#Run Harmony 
data <- data %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT",
             reduction = "pca", project.dim = FALSE, 
             reduction.save = "harmony.rna", 
             dims.use = 1:20)

harmony_embeddings <- Embeddings(data, 'harmony.rna')

# checking PCs
DepthCor(data, reduction = "harmony.rna" )
ElbowPlot(data, reduction = "harmony.rna")

nn = 30L
mindist = 0.01

data <- data %>% 
  RunUMAP(reduction = "harmony.rna", dims = 1:20 , n.neighbors = nn,  min.dist = mindist) %>% 
  FindNeighbors(reduction = "harmony.rna", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

q <- DimPlot(data, pt.size = .3, reduction = "umap", label = TRUE, cols = cols, group.by = "SCT_snn_res.1", label.size = 2, label.color = "#201c1d", repel = T, raster = FALSE )+optns +NoLegend()+NoAxes()+ggtitle("")
ggsave(filename = "fig1b.pdf", width = 10, height = 10,q)

# Adding cell type labels
# Tthis requires manual annotation of clusters
data$celltype <- ""
data$celltype[data$SCT_snn_res.1 %in% c(0, 4, 3,6, 23, 2, 5, 10, 15, 19, 24, 11, 21)]<- "mesophyll"
data$celltype[data$SCT_snn_res.1 %in% c(9, 12, 13,8, 18)]<- "vasculature"
data$celltype[data$SCT_snn_res.1 %in% c(1,14,16,22,7)]<- "epidermis"
data$celltype[data$SCT_snn_res.1 %in% c(17)]<- "undifferentiated"

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Cluster composition analysis
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

optns2 <- theme(
  axis.text.x = element_text(margin = margin(c(5, 0, 0, 0)), size = 10, 
                             face = "bold", hjust = 1),
  axis.text.y = element_text(margin = margin(c(0, 5, 0, 0)), size = 10, 
                             face = "bold"),
  axis.title = element_blank(),
  axis.ticks.y = element_line(size = .5),
  axis.ticks.x = element_line(size = .5),
  axis.ticks.length = unit(.1, "cm"),
  axis.line = element_line(size = .5),
  panel.background = element_blank(),
  legend.text = element_text(size = 15, face = "bold", family = "Helvetica", 
                             margin = margin(c(20, 0, 20, 0))),
  legend.title = element_text(size = 15, face = "bold", family = "Helvetica"),
  legend.key = element_rect(fill = NA, size = 1)
)

col.sample <- c("00_Mock_rep1" = "#005b96",
                "00_Mock_rep2" = "#005b96",
                "DC3000_04h" = "#d9c99e",
                "DC3000_06h" = "#bda155",
                "DC3000_09h" = "#a2790d",
                "DC3000_24h"= "#513c06",
                "AvrRpt2_04h_rep1"= "#66b2b2",
                "AvrRpt2_06h_rep1"= "#008080",
                "AvrRpt2_09h_rep1"= "#006666",
                "AvrRpt2_24h_rep1"= "#004c4c",
                "AvrRpt2_04h_rep2"= "#66b2b2",
                "AvrRpt2_06h_rep2"= "#008080",
                "AvrRpt2_09h_rep2"= "#006666",
                "AvrRpt2_24h_rep2"= "#004c4c",
                "AvrRpm1_04h"= "#ca6666",
                "AvrRpm1_06h"= "#af1919",
                "AvrRpm1_09h"= "#740000",
                "AvrRpm1_24h"= "#420000"
                
)

samples <- data$sample %>% unique()

l <- list()

for (i in c(1:length(samples))){
  clst_table <- data$SCT_snn_res.1[data$sample %in% samples[i]] %>% table()
  clst_table_norm <- clst_table/sum(clst_table)
  l[[i]]<-clst_table_norm 
}

clst_table <- do.call(rbind.data.frame, l) %>% t()
colnames(clst_table) <- samples
rownames(clst_table) <- names(l[[1]])

clst_table.m <- clst_table %>% melt()
clst_table.m$X1 <- as.factor(clst_table.m$X1)

q <- ggplot(clst_table.m, aes(x = X1, y = value, fill = X2)) + geom_bar(position = 'fill', stat = "identity") + optns2 
q <-q & scale_fill_manual(values = col.sample) #change colors 
ggsave("figS1i.pdf", width = 30, height = 15 ,q) 

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## ATAC Peak calling for each cluster
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
# call peaks using MACS2
DefaultAssay(data) <- "ATAC"
peaks <- CallPeaks(data, 
                   group.by = "seurat_clusters",
                   macs2.path = "/opt/anaconda3/bin/macs2",
                   effective.genome.size = 1.35e8,
                   combine.peaks = TRUE,
                   extsize = 150,
                   shift = -75
)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(data),
  features = peaks,
  cells = colnames(data)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
data[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  annotation = annotation
)

sample.list <- data$sample %>% unique()
frac.comb <- c()
for (i in c(1:length(fragpath.list))){
  
  cell.names <- colnames(data[, data$sample==sample.list[i]])
  cell.names <- gsub(".*_", "", cell.names) #cell names should match the original cell barcode, not concatenated names.
  names(cell.names) <- colnames(data[, data$sample==sample.list[i]]) #names should match the concatenated names.
  
  fragments <- CreateFragmentObject(
    path = fragpath.list[i],
    cells = cell.names,
    validate.fragments = FALSE
  )
  frac.comb <- c(frac.comb, fragments)
}

Fragments(data) <- NULL
Fragments(data) <- frac.comb


out.dir.seurat <- file.path(mainDir, "/_seurat_object")
dir.create(out.dir.seurat, showWarnings = FALSE)
saveRDS(data, paste0(out.dir.seurat, "/combined_filtered.rds"))

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## marker gene/peak calling; on RNA clustering
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
# Set the identity of your cells to the desired column
Idents(data) <- data$SCT_snn_res.1

DefaultAssay(data) <- "SCT"
markers.rna <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers.rna, file = "markers_RNA.txt", row.names=T, col.names=NA, sep="\t", quote=F) #Table S1

#marker gene plot
optns2 <- theme(
  axis.text.x = element_text(margin = margin(c(5, 0, 0, 0)), size = 15, face = "bold", family = "Helvetica"),
  axis.text.y = element_text(margin = margin(c(0, 5, 0, 0)) ,size = 15, face = "bold", family = "Helvetica"), 
  axis.title = element_text(size = 10, face = "bold", family = "Helvetica"),
  axis.title.x = element_text(margin = margin(c(10, 0, 0, 0))),
  axis.title.y = element_text(margin = margin(c(0, 10, 0, 0))),
  axis.ticks.y = element_line(size = 1),
  axis.ticks.x = element_line(size = 1),
  axis.ticks.length = unit(.25, "cm"),
  axis.line  = element_line(size = 1)
) 

idt <- Idents(data)
levels(idt) <- sort(as.numeric(levels(idt)))
Idents(data) <- idt
top.markers <- markers.rna %>% group_by(cluster) %>% top_n(n = 1, wt =avg_log2FC )
top.marker.gene <- top.markers$gene
top.markers2 <- markers.rna %>% group_by(cluster) %>% top_n(n = 2, wt =avg_log2FC )
q <- DotPlot(data, features = unique(top.marker.gene), cols = c("#ffcfff", "#221100")) + RotatedAxis() + optns2 + NoLegend()
ggsave(filename = "figS2a.pdf", width = 8, height = 8 ,q)

#marker peaks
DefaultAssay(data) <- "peaks"
marker.peaks <- FindAllMarkers(
  object = data,
  min.pct = 0.05,
  test.use = 'wilcox',
  only.pos = TRUE,
  logfc.threshold = 0.2
)

write.table(marker.peaks, file =  "markers_peaks.txt", row.names=T, col.names=NA, sep="\t", quote=F)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## GO analysis
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
dir.create("GO_analysis", showWarnings = FALSE)

genes_to_id <- function(genes){
  if (length(grep(genes, annotation$ID))==0 ){
    id <- annotation$ID[annotation$gene_name==genes] %>% na.omit
    id <- gsub("gene:", "", id)
    id[1]
  }else{
    genes
  }
}

clst <- unique(markers.rna$cluster)
for (i in c(1:length(clst))){
  gene <- markers.rna$gene[markers.rna$cluster==clst[i]]
  
  in.gene <- c()
  
  for (j in c(1:length(gene))){
    in.gene<- c(in.gene,  genes_to_id(gene[j]))
    
  }
  
  
  ego <- enrichGO(gene          = in.gene,
                  OrgDb         = org.Athaliana.eg.db,
                  ont           = "BP",
                  keyType = 'GID',
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  png(paste("GO_analysis/GO_clst_", clst[i],".png", sep = "")) # width = 1000, height = 1000
  try(print(barplot(ego, showCategory=10)))
  dev.off()
  write.table(ego@result,file = paste0("GO_analysis/go_clst", clst[i], ".txt") ,row.names=T, col.names=NA, sep="\t", quote=F)
}

#For figure; take top 2 go terms per cluster
clst <- data$SCT_snn_res.1 %>% unique()
clst.sel <- clst %>% sort(decreasing = FALSE)
n.terms <- 10
out <- matrix(data = NA, nrow = length(clst.sel)*n.terms, ncol = 3) %>% as.data.frame()
colnames(out) <- c("Description", "Size", "adjusted pvalue")


optns3 <- theme(
  axis.text.x = element_text(margin = margin(c(10, 0, 0, 0)), size = 15,  family = "Helvetica"),
  axis.text.y = element_text(margin = margin(c(0,10, 0, 0)) ,size = 15,  family = "Helvetica"), 
  axis.title = element_text(size = 10, family = "Helvetica"),
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

for (i in c(1:length(clst.sel))){
  go.data <- read.delim(file = paste0("GO_analysis/GO_clst", clst.sel[i],".txt"), header = T, row.names = 1)
  
  s <- n.terms*(i-1)+1
  
  for (j in c(1:n.terms)){
    out[s+j-1, ] <- c(go.data$Description[j],as.numeric(go.data$Count[j]), as.numeric(go.data$p.adjust[j]))
  }
}

out$clst <- as.character(rep(clst.sel, each = n.terms))
out$clst <- factor(out$clst, levels=unique(out$clst))
out$log.p <- -log10(as.numeric(out$`adjusted pvalue`))
out$Size <- as.numeric(out$Size)
out$Size <- log2(out$Size)

out <- out[out$log.p > 1 ,]

out$Description[grep("photosynthesis", out$Description)] <- "photosynthesis"
out$Description[grep("jasmonic acid", out$Description)] <- "jasmonic acid pathway"
out$Description[grep("cellular response to decreased oxygen levels", out$Description)] <- "response to hypoxia"
out$Description[grep("hypoxia", out$Description)] <- "response to hypoxia"
out$Description[grep("defense response to oomycetes", out$Description)] <- "defense response to oomycetes"
out$Description[grep( "ion homeostasis", out$Description)] <- "ion homeostasis"

out <- out[-grep("regulation of response to biotic stimulus|oomycetes|metabolic|reticulum|circulatory|acclimation|xenobiotic", out$Description) ,]

#sort GO terms
idx_immune <- grep("defense|immun|salicylic|jasmonic|resistance", out$Description)
idx_transport <- grep("transport", out$Description)
idx_metab <- grep("metabo", out$Description)
for_others <- c(idx_immune,idx_transport,idx_metab ) 

out2 <- rbind(out[for_others, ], out[-for_others, ])
out2 <- out
out2$Description <- factor(out2$Description, levels=unique(out2$Description))

q <- ggplot(out2, aes(x = clst, y = Description, size = Size , color = log.p)) + geom_point()+ optns3 + 
  scale_color_gradient(low = "grey",high = "#418557",limits = c(0,30))+
  scale_size(range = c(0,10)) 
ggsave("GO_analysis/go_summary_plot_allClst2.pdf", height = 5, width = 12 ,q)

#for main Fig. showing only immune related...
out3 <- out2
out3 <- out2[grep("salicylic acid|immune|jasmonic" ,out2$Description) ,]
y <- out2[-grep("salicylic acid|immune|jasmonic" ,out2$Description) ,]
y$Description <- "response to salicylic acid"
y$Size <- NA
y$log.p <- NA

out4 <- rbind(out3, y)

q <- ggplot(out4, aes(x = Description, y = clst, size = Size , color = log.p)) + geom_point()+ optns3 + 
  scale_color_gradient(low = "grey",high = "#418557",limits = c(0,30))+
  scale_size(range = c(0,12)) 
ggsave("GO_analysis/fig1c.pdf", height = 10, width = 6 ,q)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## RNA clustering - gene expression analysis
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

#combining reps
sample2 <- data$sample
sample2[grep("00_Mock_rep1|00_Mock_rep2", sample2)] <- "00_Mock"
sample2[grep("AvrRpt2_04h_rep1|AvrRpt2_04h_rep2", sample2)] <- "AvrRpt2_04h"
sample2[grep("AvrRpt2_06h_rep1|AvrRpt2_06h_rep2", sample2)] <- "AvrRpt2_06h"
sample2[grep("AvrRpt2_09h_rep1|AvrRpt2_09h_rep2", sample2)] <- "AvrRpt2_09h"
sample2[grep("AvrRpt2_24h_rep1|AvrRpt2_24h_rep2", sample2)] <- "AvrRpt2_24h"

data$sample2 <- sample2

# each strain separately
cell_DC3000 <- rownames(data@meta.data)[ grepl("DC3000", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] 
cell_AvrRpt2 <- rownames(data@meta.data)[ grepl("AvrRpt2", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ]
cell_AvrRpm1 <- rownames(data@meta.data)[ grepl("AvrRpm1", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ]
cell_mock <- rownames(data@meta.data)[ grepl("Mock", data@meta.data[, "sample"]) ]

data_kt56 <- subset(x = data, cells = cell_DC3000)
data_kt57 <- subset(x = data, cells = cell_AvrRpt2)
data_kt58 <- subset(x = data, cells = cell_AvrRpm1)
data_mock <- subset(x = data, cells = cell_mock)

DefaultAssay(data) <- "SCT"
DefaultAssay(data_kt56) <- "SCT"
DefaultAssay(data_kt57) <- "SCT"
DefaultAssay(data_kt58) <- "SCT"
DefaultAssay(data_mock) <- "SCT"

dir.create("RNA_featurePlot", showWarnings = FALSE)

#plot each replicate separately. This is a function to make FeaturePlots, such as Fig2d
plotting <- function(gene_list, fname){
  q <- FeaturePlot(data_kt56, features = gene_list, split.by = "sample", pt.size = 0.5, order=TRUE, max.cutoff = "q99", min.cutoff = "q1", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste("RNA_featurePlot/","_", fname,"_feture_KT56.pdf", sep = ""), width = 20, height = 20*length(gene_list)/length(levels(factor(data_kt56@meta.data[,"sample"]))), q)
  
  
  q <- FeaturePlot(data_kt57, features = gene_list, split.by = "sample", pt.size = 0.5, order=TRUE, max.cutoff = "q99", min.cutoff = "q1", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste( "RNA_featurePlot/","_", fname,"_feture_KT57.pdf", sep = ""), width = 20, height = 20*length(gene_list)/length(levels(factor(data_kt57@meta.data[,"sample"]))), q)
  
  
  q <- FeaturePlot(data_kt58, features = gene_list, split.by = "sample", pt.size = 0.5, order=TRUE, max.cutoff = "q99", min.cutoff = "q1", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste( "RNA_featurePlot/","_", fname,"_feture_KT58.pdf", sep = ""), width = 20, height = 20*length(gene_list)/length(levels(factor(data_kt58@meta.data[,"sample"]))), q)
  
  q <- FeaturePlot(data, features = gene_list, order=TRUE, pt.size = .5, max.cutoff = "q99", min.cutoff = "q1", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste( "RNA_featurePlot/","_", fname,".pdf", sep = ""), width = 10, height = 10, q)
  
}

# examples
plotting(c( "CRK7.1" ), "CRK7.1")
plotting(c( "MAM1" ), "MAM1")
plotting(c( "UGT85A1" ), "UGT85A1")
plotting(c( "VSP2" ), "VSP2")
plotting(c( "CBP60G" ), "CBP60G")
plotting(c( "WRKY46" ), "WRKY46")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## ATAC feature plots
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
DefaultAssay(data) <- "peaks"
DefaultAssay(data_kt56) <- "peaks"
DefaultAssay(data_kt57) <- "peaks"
DefaultAssay(data_kt58) <- "peaks"

dir.create("ATAC_featurePlot", showWarnings = FALSE)

atac.feature <- function(gene_list, fname){
  q <- FeaturePlot(data, features = gene_list, order=TRUE, pt.size = .3, max.cutoff = "q99", min.cutoff = "q20", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray",  "#FCC5C0" ,"#FA9FB5" ,"#F768A1", "#DD3497" ,"#AE017E" ,"#7A0177" ,"#49006A"),
    na.value = "lightgray"
  )& NoAxes()
  ggsave(filename = paste("ATAC_featurePlot/","_", fname,".png", sep = ""), width = 5, height = 5, q)
}

region_name <- c("2-11172821-11173529")

atac.feature(region_name, "FDH_figS2c")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## coverage plots; FigS2b
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
region_name <- "2-11172834-11173492"
Annotation(data) <- annotation
q <- CoveragePlot(data, region = region_name,
                  assay = "peaks",
                  annotation = TRUE,
                  extend.upstream = 5000,
                  extend.downstream =5000, 
                  heights = c(5,1,1),
                  group.by = "SCT_snn_res.1"
                  
)
q <-q & scale_fill_manual(values = cols) #change colors 
ggsave(filename = paste0("coveragePlot/", region_name, ".pdf") ,width = 15, height = 10,q) #figS2b

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## ATAC clustering
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

out.dir.atac <- file.path(mainDir, subDir ,"ATAC")
dir.create(out.dir.atac, showWarnings = FALSE)
setwd(out.dir.atac)

DefaultAssay(data) <- "peaks"
data <- FindTopFeatures(data, min.cutoff = "q5")
data <- RunTFIDF(data, scale.factor = 100000)
data <- RunSVD(data)

DepthCor(data, reduction = "lsi" )
ElbowPlot(data, reduction = "lsi")

data <- RunUMAP(data, dims = 2:7, reduction = "lsi",
                n.neighbors = 20, #Alex paper 
                min.dist = 0.01 #Alex paper
)
data <- FindNeighbors(data, dims = 2:7, reduction = "lsi")
data <- FindClusters(data, resolution = 1, verbose = FALSE, algorithm = 1)

#harmony
DefaultAssay(data) <- "peaks"

#Run Harmony 
options(repr.plot.height = 2.5, repr.plot.width = 6)
data <- data %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "peaks", 
             reduction = "lsi", project.dim = FALSE, 
             reduction.save = "harmony.lsi", 
             dims.use = 2:10 
  ) #project.dim = FALSE is needed to avoid error 

harmony_embeddings <- Embeddings(data, 'harmony.lsi')
harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = data, reduction = "harmony.lsi", pt.size = .1, group.by = "sample") #removed " do.return = TRUE"
p2 <- VlnPlot(object = data, features = "harmony_1", group.by = "sample", pt.size = .1)
plot_grid(p1,p2)


DepthCor(data, reduction = "harmony.lsi" )
ElbowPlot(data, reduction = "harmony.lsi")


data <- data %>% 
  RunUMAP(reduction = "harmony.lsi", dims = 2:20, n.neighbors = 30L,  min.dist = 0.01 ) %>% 
  FindNeighbors(reduction = "harmony.lsi", dims = 2:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()



q <- DimPlot(data, pt.size = .3, reduction = "umap", label = TRUE, group.by = "peaks_snn_res.0.8", cols = cols,label.size = 7, label.color = "#505051", raster = FALSE)+optns +NoLegend()+NoAxes()+ggtitle("")
ggsave(filename = "figS1i",  width = 10, height = 10, q)

saveRDS(data, paste0(out.dir.seurat, "/combined_filtered.rds"))


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## ATAC gene activity score 
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

#update annotation file by adding gene_ids whereever gene_name is absent
annotation2 <- annotation
id <- (annotation$gene_biotype %in% "protein_coding") * is.na(annotation$gene_name)
id2 <- c(1:length(annotation))[!!id]
annotation2$gene_name[id2] <- annotation$gene_id[id2]
Annotation(data) <- annotation2

DefaultAssay(data) <- "peaks"
gene.activities <- GeneActivity(data,
                                extend.upstream = 400, #default = 2000,
                                biotypes = "protein_coding"
)

data[['atacRNA_400bp']] <- CreateAssayObject(counts = gene.activities)
data <- NormalizeData(
  object = data,
  assay = 'atacRNA_400bp',
  normalization.method = 'LogNormalize',
  scale.factor = median(data$nCount_RNA)
)

saveRDS(data, paste0(out.dir.seurat, "/combined_filtered.rds"))

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## rep1 vs rep2 correlation
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

optns <- theme(
  axis.text.x = element_text(margin = margin(c(10, 0, 0, 0)), size = 10, family = "Helvetica"),
  axis.text.y = element_text(margin = margin(c(0, 10, 0, 0)) ,size = 10,  family = "Helvetica"), 
  axis.title = element_text(size = 20, family = "Helvetica"),
  axis.title.x = element_text(margin = margin(c(0, 0, 0, 0))),
  axis.title.y = element_text(margin = margin(c(0, 0, 0, 0))),
  axis.ticks.y = element_line(size = 1),
  axis.ticks.x = element_line(size = 1),
  axis.ticks.length = unit(.1, "cm"),
  axis.line  = element_line(size = 1),
  panel.background = element_blank(),
  plot.background =element_blank(),
  legend.background = element_blank(),
  legend.text = element_text(size = 15, face = "bold", family = "Helvetica", margin = margin(c(20, 0, 20, 0))),
  legend.title = element_text(size = 15, face = "bold", family = "Helvetica"),
  legend.key = element_rect(fill = NA , size= 5)
) 

#making pseudobulk cout for each sample
DefaultAssay(data) <- "RNA"
samples <- data$sample %>% unique() %>% sort(decreasing = F)
out <- matrix(data = NA, nrow = nrow(data@assays$RNA), ncol = length(samples))
rownames(out) <- rownames(data)
colnames(out) <- samples

for (i in c(1:length(samples))){
  data.sample <- data[, data$sample==samples[i]]
  pb <- data.sample@assays$RNA@counts %>% rowSums()
  tpm <- pb/sum(pb)*1000000 
  tpm <- log2(tpm)
  tpm[is.infinite(tpm)] <-0
  out[, i] <- tpm
}

out2 <- out %>% as.data.frame()

pdf("mock_rep1_rep2_scatter.pdf")
plot(out2[,1], out2[,2])
dev.off()
cor(out2[,1], out2[,2]) 

pdf("KT57_rep1_rep2_scatter.pdf")
plot(out2[,9], out2[,10])
dev.off()
cor(out2[,9], out2[,10])

q1 <- qplot(out2[,1], out2[,2], data = out2)+optns
q2 <- qplot(out2[,9], out2[,10], data = out2)+optns
q3 <- qplot(out2[,1], out2[,9], data = out2)+optns
q4 <- qplot(out2[,1], out2[,10], data = out2)+optns
q5 <- qplot(out2[,2], out2[,9], data = out2)+optns
q6 <- qplot(out2[,2], out2[,10], data = out2)+optns
q <- arrangeGrob(q1, q2,q3, q4,q5, q6, nrow = 2)
ggsave("figS1h.pdf",  height = 10, width = 15, q)

out.filter <- out[rowSums(out)>0, ]
pca <- prcomp(t(out.filter), scale. = T)
pca2 <- pca$x
summary(pca)
pca3 <- as.data.frame(pca2)
rownames(pca3) <- gsub("00_M", "M", rownames(pca3))

col.sample.pca <- col.sample[rownames(pca3)]

q <- ggplot(pca3, aes(x=PC1, y = PC2, label = rownames(pca3) )) 
q <- q + optns
q <- q + xlab("PC1 (36.1 %)") + ylab("PC2 (14.1 %)") #summary(pca)
q <- q + geom_point(size = 10, color = col.sample.pca)+ scale_size(guide = 'none', range = c(2, 20))
q <- q + guides(colour = guide_legend(override.aes = list(size= 5, alpha = 1, shape = 15, keywidth = .1, keyheight = .1, default.unit = "inch"), order = 10)) + 
  guides(shape = guide_legend(override.aes = list(size= 5, keywidth = .1, keyheight = .1, default.unit = "inch"), reverse = F)) 
q <- q + geom_text_repel(force = 30, size = 8, hjust=0.5) 
ggsave("figS1g.pdf",  height = 15, width = 15, q)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## clustering for KT57 alone (for MERFISH integration)
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

DefaultAssay(data_kt57) <- "RNA"

data_kt57 <- SCTransform(data_kt57, vars.to.regress = "percent.mt")
data_kt57 <- RunPCA(data_kt57, reduction.name = "pca_kt57")

#Harmony
DefaultAssay(data_kt57) <- "SCT"

#Run Harmony 
options(repr.plot.height = 2.5, repr.plot.width = 6)
data_kt57 <- data_kt57 %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT",
             reduction = "pca_kt57", project.dim = FALSE, 
             reduction.save = "harmony.rna.kt57", 
             dims.use = 1:20)

harmony_embeddings <- Embeddings(data_kt57, 'harmony.rna.kt57')

data_kt57 <- data_kt57 %>% 
  RunUMAP(reduction = "harmony.rna.kt57", dims = 1:15 , n.neighbors = 30L,  min.dist = 0.01, 
          # reduction.name = "umap.harmony_mindist_0.0013_nn_10"
  ) %>% 
  FindNeighbors(reduction = "harmony.rna.kt57", dims = 1:15) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

tol24rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", "grey50", "grey25","black")
names(tol24rainbow) <- c(0:23)

q <- DimPlot(data_kt57, pt.size = .5, reduction = "umap", label = TRUE, cols = tol24rainbow, group.by = "SCT_snn_res.1", label.size = 7, label.color = "#201c1d", repel = T )+optns +NoLegend()+NoAxes()+ggtitle("")
ggsave(filename = "umap_harmony_RNA_noLegend2.pdf", width = 10, height = 10,q)

saveRDS(data_kt57, file = "AvrRpt2_alone2.rds")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## joint embedding; figS1k
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
out.dir.joint <- file.path(mainDir, subDir ,"joint")
dir.create(out.dir.joint, showWarnings = FALSE)
setwd(out.dir.joint)

DefaultAssay(data) <- "SCT"

# build a joint neighbor graph using both assays
data <- FindMultiModalNeighbors(
  object = data,
  reduction.list = list("harmony.rna", "harmony.lsi"),
  dims.list = list(1:20, 2:20),
  modality.weight.name = "RNA.weight",
  verbose = TRUE,
  return.intermediate = TRUE
)

modality.weight.ratio <- data@misc$modality.weight@modality.weight.list$harmony.rna/data@misc$modality.weight@modality.weight.list$harmony.lsi  #RNA/ATAC

# build a joint UMAP visualization
data <- RunUMAP(
  object = data,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE,
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  n.neighbors = 30L,
  min.dist = 0.1
)

data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)

cols <- c("0" = "#CCCCCC",
          "1" = "#721c28",
          "2" = "#666666",
          "3" = "#18392b",
          "4" = "#5d009f",
          "5" = "#eacd01",
          "6" = "#B2DF8A",
          "7" = "#33A02C",
          "8" = "#E78AC3",
          "9" = "#2c2cff",
          "10" = "#A65628",
          "11" = "#FB8072",
          "12" = "#D95F02",
          "13" = "#004242",
          "14" = "#F4CAE4",
          "15" = "#984EA3",
          "16" = "#418557",
          "17" = "#009088",
          "18" = "#08ff00",
          "19" = "#FB9A99",
          "20" = "#e36700",
          "21" = "#008fff",
          "22" = "#A6CEE3",
          "23" = "#1B9E77",
          "24" = "#819d26",
          "25" = "#ff3377",
          "26" = "#f70058",
          "27" = "#fe5757",
          "28" = "#E41A1C",
          "29" = "#6A3D9A",
          "30" = "#333333",
          "31" = "#1F78B4",
          "32" = "#7570B3",
          "33" = "#B3DE69",
          "34" = "#8DD3C7",
          "35" = "#3a451c",
          "36" = "#666666",
          "37" = "#A6761D"
          
)

q <- DimPlot(data, reduction = "wnn.umap", pt.size = .3,  label = T, group.by = "wsnn_res.0.5", cols = cols, label.size = 7, label.color = "#201c1d", repel = T, raster = FALSE)+NoLegend()+NoAxes()+ggtitle("")
ggsave("figS1k.pdf", width = 10, height = 10, q)


# ATAC marker peaks 
DefaultAssay(data) <- "peaks"
marker.peaks <- FindAllMarkers(
  object = data,
  min.pct = 0.05,
  test.use = 'wilcox',
  only.pos = TRUE,
  logfc.threshold = 0.2
)
write.table(marker.peaks, file =  "markers_peaks.txt", row.names=T, col.names=NA, sep="\t", quote=F)
