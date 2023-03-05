source(".../_config_multiome.R")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Color schemes
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

cols <- c("0" = "#91a3b0",
          "1" = "#721c28",
          "2" = "#2e5073",
          "3" = "#94ac78",
          "4" = "#3f5246",
          "5" = "#ff4c4c",
          "6" = "#418557",
          "7" = "#af2b3d",
          "8" = "#2c2cff",
          "9" = "#eacd01",
          "10" = "#9f59f8",
          "11" = "#e06666",
          "12" = "#7c5914",
          "13" = "#ff9f22",
          "14" = "#d1bac8",
          "15" = "#ad3074",
          "16" = "#af5b2b",
          "17" = "#009088",
          "18" = "#93d0fc",
          "19" = "#752f9a",
          "20" = "#009580",
          "21" = "#6dc831",
          "22" = "#d58282",
          "23" = "#2aaae1",
          "24" = "#fc0ed5")


col.sample <- c("00_Mock" = "#005b96",
                "01_DC3000_04h" = "#d9c99e",
                "01_DC3000_06h" = "#bda155",
                "01_DC3000_09h" = "#a2790d",
                "01_DC3000_24h"= "#513c06",
                "02_AvrRpt2_04h"= "#66b2b2",
                "02_AvrRpt2_06h"= "#008080",
                "02_AvrRpt2_09h"= "#006666",
                "02_AvrRpt2_24h"= "#004c4c",
                "03_AvrRpm1_04h"= "#ca6666",
                "03_AvrRpm1_06h"= "#af1919",
                "03_AvrRpm1_09h"= "#740000",
                "03_AvrRpm1_24h"= "#420000")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Functions
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
coverage.plot <- function(data, target, out.dir){
  tryCatch({data <- LinkPeaks(
    object = data,
    peak.assay = "peaks",
    expression.assay = "SCT",
    genes.use =target
  )
  q <- CoveragePlot(
    object = data,
    region = target,
    features = target,
    expression.assay = "SCT",
    extend.upstream = 3000,
    extend.downstream = 3000,
    group.by = "sample.order",
    split.assays = TRUE,
    heights = c(2,2,2,.5)



  )
  q <-q & scale_fill_manual(values = col.sample) #change colors
  ggsave(filename = paste0(out.dir,"/coveragePlot_marker_", target ,"_RNAcluster_sample.png"), width = 20, height = 10, q)
  
  q <- CoveragePlot(
    object = data,
    region = target,
    features = target,
    expression.assay = "SCT",
    #idents = idents.plot,
    extend.upstream = 3000,
    extend.downstream =3000, 
    group.by = "SCT_snn_res.1",
    #split.assays = TRUE,
    heights = c(2,2,2,.5)
    
  )
  q <-q & scale_fill_manual(values = cols) #change colors 
  ggsave(filename = paste0(out.dir,"/coveragePlot_marker_", target ,"_RNAcluster.png"), width = 20, height = 10, q)
  }, error=function(e){})
}

feature.plot.atac <- function(peak, fname){
  q <- FeaturePlot(data_dc3000, 
                   features = peak, 
                   split.by = "sample",
                   pt.size = 2, 
                   order=TRUE,
                   min.cutoff = 'q1',
                   max.cutoff = 'q99',
                   reduction = "wnn.umap",
                   cols = c("lightgrey","darkmagenta"),
  )
  ggsave(filename = paste("ATAC_featurePlot/","_", fname,"_feture_dc3000.png", sep = ""), width = 20, height = 5, q)
  q <- FeaturePlot(data_avrrpt2, 
                   features = peak, 
                   split.by = "sample",
                   pt.size = 2, 
                   order=TRUE,
                   min.cutoff = 'q1',
                   max.cutoff = 'q99',
                   reduction = "wnn.umap",
                   cols = c("lightgrey","darkmagenta"),
  )
  ggsave(filename = paste("ATAC_featurePlot/","_", fname,"_feture_avrrpt2.png", sep = ""), width = 20, height = 5, q)
  q <- FeaturePlot(data_avrrpm1, 
                   features = peak, 
                   split.by = "sample",
                   pt.size = 2, 
                   order=TRUE,
                   min.cutoff = 'q1',
                   max.cutoff = 'q99',
                   reduction = "wnn.umap",
                   cols = c("lightgrey","darkmagenta"),
  )
  ggsave(filename = paste("ATAC_featurePlot/","_", fname,"_feture_avrrpm1.png", sep = ""), width = 20, height = 5, q)
  
  q <- FeaturePlot(data, 
                   features = peak, 
                   # split.by = "sample", 
                   pt.size = 2, 
                   order=TRUE,
                   min.cutoff = 'q1',
                   max.cutoff = 'q99',
                   reduction = "wnn.umap",
                   cols = c("lightgrey","darkmagenta"),
  )
  ggsave(filename = paste("ATAC_featurePlot/","_", fname,".png", sep = ""), width = 10, height = 10, q)
  
}

sct.atac.coplot <- function(gene_list, fname){
  
  tryCatch({
    DefaultAssay(data_dc3000) <- "SCT"
    q1 <- FeaturePlot(
      object = data_dc3000,
      features = gene_list,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = 1,
      split.by = "sample",
      order = TRUE,
      reduction = "wnn.umap",
      cols = c("lightgrey","darkgreen")
    )
    
    DefaultAssay(data_dc3000) <- "atacRNA_400bp"
    q2 <- FeaturePlot(data_dc3000, 
                      features = gene_list, 
                      split.by = "sample", 
                      pt.size = 1, 
                      order=TRUE,
                      min.cutoff = 'q1',
                      max.cutoff = 'q99',
                      reduction = "wnn.umap",
                      cols = c("lightgrey","darkmagenta")
    )
    p1 <- ggarrange(q1, q2, nrow  =2)
    ggsave(filename = paste("RNA_ATACactivity_plot/" ,fname,"_feture_dc3000.png", sep = ""), width = 16, height = 8, p1) 
  }, error=function(e){})
  
  tryCatch({
    ##########
    DefaultAssay(data_avrrpt2) <- "SCT"
    
    q1 <- FeaturePlot(
      object = data_avrrpt2,
      features = gene_list,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = 1,
      split.by = "sample",
      order = TRUE,
      reduction = "wnn.umap",
      cols = c("lightgrey","darkgreen")
    )
    
    DefaultAssay(data_avrrpt2) <- "atacRNA_400bp"
    q2 <- FeaturePlot(data_avrrpt2, 
                      features = gene_list, 
                      split.by = "sample", 
                      pt.size = 1, 
                      order=TRUE,
                      min.cutoff = 'q1',
                      max.cutoff = 'q99',
                      reduction = "wnn.umap",
                      cols = c("lightgrey","darkmagenta")
    )
    p1 <- ggarrange(q1, q2, nrow  =2)
    ggsave(filename = paste("RNA_ATACactivity_plot/" , fname,"_feture_avrrpt2.png", sep = ""), width = 16, height = 8, p1) 
  }, error=function(e){})
  
  ##########
  tryCatch({
    DefaultAssay(data_avrrpm1) <- "SCT"
    q1 <- FeaturePlot(
      object = data_avrrpm1,
      features = gene_list,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = 1,
      split.by = "sample",
      order = TRUE,
      reduction = "wnn.umap",
      cols = c("lightgrey","darkgreen")
    )
    
    DefaultAssay(data_avrrpm1) <- "atacRNA_400bp"
    q2 <- FeaturePlot(data_avrrpm1, 
                      features = gene_list, 
                      split.by = "sample", 
                      pt.size = 1, 
                      order=TRUE,
                      min.cutoff = 'q1',
                      max.cutoff = 'q99',
                      reduction = "wnn.umap",
                      cols = c("lightgrey","darkmagenta")
    )
    p1 <- ggarrange(q1, q2, nrow  =2)
    ggsave(filename = paste("RNA_ATACactivity_plot/" , fname,"_feture_avrrpm1.png", sep = ""), width = 16, height = 8, p1) 
  }, error=function(e){})
  
  ##########
  tryCatch({
    DefaultAssay(data) <- "SCT"
    q1 <- FeaturePlot(
      object = data,
      features = gene_list,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = 1,
      # split.by = "sample",
      order = TRUE,
      reduction = "wnn.umap",
      cols = c("lightgrey","darkgreen")
    )
    
    DefaultAssay(data) <- "atacRNA_400bp"
    q2 <- FeaturePlot(data, 
                      features = gene_list, 
                      # split.by = "sample", 
                      pt.size = 1, 
                      order=TRUE,
                      min.cutoff = 'q1',
                      max.cutoff = 'q99',
                      reduction = "wnn.umap",
                      cols = c("lightgrey","darkmagenta")
    )
    p1 <- ggarrange(q1, q2, nrow  =2)
    ggsave(filename = paste("RNA_ATACactivity_plot/" , fname,"_feture.png", sep = ""), width = 4, height = 8, p1) 
  }, error=function(e){})
  
}
  
kclust <- function(input, n){
  km<- kmeans(input,n)
  m.kmeans<- cbind(input, km$cluster)
  colnames(m.kmeans)[ncol(m.kmeans)] <- "kmean"
  m.kmeans <- as.data.frame(m.kmeans)
  m.kmeans <- tibble::rownames_to_column(m.kmeans)
  o <-
    m.kmeans %>%
    arrange(kmean)
  
  data_km <- column_to_rownames(o)
  data_km
}


optns <- theme(
  axis.text.x = element_text(margin = margin(c(20, 0, 0, 0)), size = 20, face = "bold", family = "Helvetica"),
  axis.text.y = element_text(margin = margin(c(0, 20, 0, 0)) ,size = 20, face = "bold", family = "Helvetica"), 
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

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Load a combined seurat object
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

data <- readRDS(".../combined_filtered.rds")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Link detection
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
seqnames(BSgenome.Athaliana.TAIR.TAIR9) <- gsub("^Chr", "", seqnames(BSgenome.Athaliana.TAIR.TAIR9))
data <- RegionStats(data, genome = BSgenome.Athaliana.TAIR.TAIR9)

###each strain separately
cell_dc3000 <- rownames(data@meta.data)[ grepl("DC3000", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] #4100
cell_avrrpt2 <- rownames(data@meta.data)[ grepl("AvrRpt2", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] #4100
cell_avrrpm1 <- rownames(data@meta.data)[ grepl("AvrRpm1", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] #4100
cell_mock <- rownames(data@meta.data)[grepl("Mock", data@meta.data[, "sample"]) ] #4100

data_dc3000 <- subset(x = data, cells = cell_dc3000)
data_avrrpt2 <- subset(x = data, cells = cell_avrrpt2)
data_avrrpm1 <- subset(x = data, cells = cell_avrrpm1)
data_mock <- subset(x = data, cells = cell_mock)

DefaultAssay(data) <- "peaks"
DefaultAssay(data_dc3000) <- "peaks"
DefaultAssay(data_avrrpt2) <- "peaks"
DefaultAssay(data_avrrpm1) <- "peaks"
DefaultAssay(data_mock) <- "peaks"

#genome wide analysis of linkages in each condition
data_mock <- LinkPeaks(
  object = data_mock,
  peak.assay = "peaks",
  expression.assay = "SCT",
)

links_mock_all <- data_mock@assays$peaks@links
write.table(links_mock_all, file =  "linkage_mock_ALL.txt", row.names=T, col.names=NA, sep="\t", quote=F)

data_dc3000 <- LinkPeaks(
  object = data_dc3000,
  peak.assay = "peaks",
  expression.assay = "SCT")

links_dc3000_all <- data_dc3000@assays$peaks@links
write.table(links_dc3000_all, file =  "linkage_dc3000_ALL.txt", row.names=T, col.names=NA, sep="\t", quote=F)

data_avrrpt2 <- LinkPeaks(
  object = data_avrrpt2,
  peak.assay = "peaks",
  expression.assay = "SCT")

links_avrrpt2_all <- data_avrrpt2@assays$peaks@links
write.table(links_avrrpt2_all, file =  "linkage_avrrpt2_ALL.txt", row.names=T, col.names=NA, sep="\t", quote=F)

data_avrrpm1 <- LinkPeaks(
  object = data_avrrpm1,
  peak.assay = "peaks",
  expression.assay = "SCT")

links_avrrpm1_all <- data_avrrpm1@assays$peaks@links
write.table(links_avrrpm1_all, file =  "linkage_avrrpm1_ALL.txt", row.names=T, col.names=NA, sep="\t", quote=F)


data <- LinkPeaks(
  object = data,
  peak.assay = "peaks",
  expression.assay = "SCT")

links_all <- data@assays$peaks@links
write.table(links_all, file =  "linkage_combined_ALL.txt", row.names=T, col.names=NA, sep="\t", quote=F)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Linkage meta analysis
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
links_mock_all$sample <- "Mock"
links_dc3000_all$sample <- "DC3000"
links_avrrpt2_all$sample <- "AvrRpt2"
links_avrrpm1_all$sample <- "AvrRpm1"
links_all$sample <- "All"


score.threshod <- 0.05
top.links_mock_all <- links_mock_all[links_mock_all$score > score.threshod ,]
top.links_dc3000_all <- links_dc3000_all[links_dc3000_all$score > score.threshod ,]
top.links_avrrpt2_all <- links_avrrpt2_all[links_avrrpt2_all$score > score.threshod ,]
top.links_avrrpm1_all <- links_avrrpm1_all[links_avrrpm1_all$score > score.threshod ,]
top.links_all <- links_all[links_all$score > score.threshod ,]

##==##==####==##==####==##==####==##==##
indata <- rbind(as.data.frame(top.links_dc3000_all),
                as.data.frame(top.links_avrrpt2_all),
                as.data.frame(top.links_avrrpm1_all),
                as.data.frame(top.links_mock_all),
                as.data.frame(top.links_all)
)
q <- ggplot(indata, aes(x=width, y=score, color=sample)) + geom_point(size=.5) + xlim(0, 5000) + optns
ggsave("score_width_plot_short.png", height = 5, width = 8, q)
q <- ggplot(indata, aes(x=width, y=score, color=sample)) + geom_point(size=.5)  + optns
ggsave("score_width_plot.png", height = 5, width = 8, q)


#density plot of scores
q <- ggplot(indata, aes(x=width, color=sample)) + geom_density()+ xlim(0, 5000) + optns
ggsave("linkWidth_density_short.png", height = 5, width = 8, q)
q <- ggplot(indata, aes(x=width, color=sample)) + geom_density() + optns
ggsave("linkWidth_density_all.png", height = 5, width = 8, q)
