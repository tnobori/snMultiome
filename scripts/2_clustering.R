
source(".../_config_multiome.R")

out_dir <- "/path/to/output/directory"

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Load annotation
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

annotation <- rtracklayer::import( ".../At_reference/Arabidopsis_thaliana.TAIR10.52.gff3")
colnames(annotation@elementMetadata)[colnames(annotation@elementMetadata)=="Name"] <- "gene_name" 
colnames(annotation@elementMetadata)[colnames(annotation@elementMetadata)=="biotype"] <- "gene_biotype" 
levels(annotation$type)[levels(annotation$type)=="gene"]="body"
annotation$gene_name[grepl("exon", annotation$gene_name)] <- gsub("\\.[0-9].exon[0-9]{1,}", "", annotation$gene_name[grepl("exon", annotation$gene_name)]) #remove exon number

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Plot parameters
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

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
  legend.key = element_rect(fill = NA , size= 5),
 ) 



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
          "24" = "#fc0ed5"
          
)


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

cols.celltype <- c("epidermis" = "#721c28",
                   "mesophyll" = "#94ac78",
                   "undifferentiated" = "#009088",
                   "vasculature" = "#ff9f22"
)


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## RNA clustering
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
odata <- readRDS("combined_filtered.rds")

#Gene expression data processing
DefaultAssay(data) <- "RNA"
data <- SCTransform(data, vars.to.regress = "percent.mt")
data <- RunPCA(data)

# checking PCs
DepthCor(data, reduction = "pca" )
ElbowPlot(data, reduction = "pca")

data <- FindNeighbors(data, dims = 1:20)
data <- FindClusters(data, resolution = 1, verbose = TRUE)
data <- RunUMAP(data, dims = 1:20,
                n.neighbors = 20, 
                min.dist = 0.01)

#Harmony
DefaultAssay(data) <- "SCT"

#Run Harmony 
options(repr.plot.height = 2.5, repr.plot.width = 6)
data <- data %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT",
             reduction = "pca", project.dim = FALSE, 
             reduction.save = "harmony.rna", 
             dims.use = 1:20)

harmony_embeddings <- Embeddings(data, 'harmony.rna')

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = data, reduction = "harmony.rna", pt.size = .1, group.by = "sample") #removed " do.return = TRUE"
p2 <- VlnPlot(object = data, features = "harmonyrna_1", group.by = "sample", pt.size = .1)
plot_grid(p1,p2)

# checking PCs
DepthCor(data, reduction = "harmony.rna" )
ElbowPlot(data, reduction = "harmony.rna")

data <- data %>% 
  RunUMAP(reduction = "harmony.rna", dims = 1:20 , n.neighbors = 30L,  min.dist = 0.01) %>% 
  FindNeighbors(reduction = "harmony.rna", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## RNA clustering - gene expression
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

###each strain separately
cell_DC3000 <- rownames(data@meta.data)[ grepl("DC3000", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] #4100
cell_AvrRpt2 <- rownames(data@meta.data)[ grepl("AvrRpt2", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] #4100
cell_AvrRpm1 <- rownames(data@meta.data)[ grepl("AvrRpm1", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] #4100
cell_mock <- rownames(data@meta.data)[  grepl("Mock", data@meta.data[, "sample"]) ] #4100

data_dc3000 <- subset(x = data, cells = cell_DC3000)
data_avrrpt2 <- subset(x = data, cells = cell_AvrRpt2)
data_avrrpm1 <- subset(x = data, cells = cell_AvrRpm1)
data_mock <- subset(x = data, cells = cell_mock)

DefaultAssay(data) <- "SCT"
DefaultAssay(data_dc3000) <- "SCT"
DefaultAssay(data_avrrpt2) <- "SCT"
DefaultAssay(data_avrrpm1) <- "SCT"
DefaultAssay(data_mock) <- "SCT"

# function to generate RNA feature plots
plotting <- function(gene_list, fname){
  q <- FeaturePlot(data_dc3000, features = gene_list, split.by = "sample", pt.size = 2, order=TRUE, max.cutoff = "q99", min.cutoff = "q1", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste(out.dir.rna, "/RNA_featurePlot/","_", fname,"_feture_dc3000.png", sep = ""), width = 20, height = 20*length(gene_list)/length(levels(factor(data_dc3000@meta.data[,"sample"]))), q)
  
  
  q <- FeaturePlot(data_avrrpt2, features = gene_list, split.by = "sample", pt.size = 2, order=TRUE, max.cutoff = "q99", min.cutoff = "q1", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste(out.dir.rna, "/RNA_featurePlot/","_", fname,"_feture_avrrpt2.png", sep = ""), width = 20, height = 20*length(gene_list)/length(levels(factor(data_avrrpt2@meta.data[,"sample"]))), q)
  
  
  q <- FeaturePlot(data_avrrpm1, features = gene_list, split.by = "sample", pt.size = 2, order=TRUE, max.cutoff = "q99", min.cutoff = "q1", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste(out.dir.rna, "/RNA_featurePlot/","_", fname,"_feture_avrrpm1.png", sep = ""), width = 20, height = 20*length(gene_list)/length(levels(factor(data_avrrpm1@meta.data[,"sample"]))), q)
  
  q <- FeaturePlot(data, features = gene_list, order=TRUE, pt.size = 2, max.cutoff = "q99", min.cutoff = "q1", reduction = "umap")
  q <- q & scale_colour_gradientn(
    colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
    na.value = "lightgray"
  ) & NoAxes()
  ggsave(filename = paste(out.dir.rna, "/RNA_featurePlot/","_", fname,".png", sep = ""), width = 10, height = 10, q)
  
}
plotting(c( "genename" ), "genename")




##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## ATAC
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Peak calling for each cluster
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

fragpath.list <- c(".../atac_fragments1.tsv.gz",
                   ".../atac_fragments2.tsv.gz")

# create a new assay using the MACS2 peak set and add it to the Seurat object
data[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  annotation = annotation
)

sample.list <- data$sample %>% unique()
frac.comb <- c()
for (i in c(1:length(fragpath.list))){
  
  cell.names <- colnames(data[, data$sample==sample.list[i]])
  cell.names <- gsub(".*_", "", cell.names) 
  names(cell.names) <- colnames(data[, data$sample==sample.list[i]]) 
  
  fragments <- CreateFragmentObject(
    path = fragpath.list[i],
    cells = cell.names,
    validate.fragments = FALSE
  )
  frac.comb <- c(frac.comb, fragments)
}

Fragments(data) <- NULL
Fragments(data) <- frac.comb

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## ATAC clustering
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

DefaultAssay(data) <- "peaks"
data <- FindTopFeatures(data, min.cutoff = "q5")
data <- RunTFIDF(data, scale.factor = 100000)
data <- RunSVD(data)

DepthCor(data, reduction = "lsi" ) 
ElbowPlot(data, reduction = "lsi")

data <- RunUMAP(data, dims = 2:10, reduction = "lsi",
                n.neighbors = 20, 
                min.dist = 0.01)
data <- FindNeighbors(data, dims = 2:10, reduction = "lsi")
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
  ) 

harmony_embeddings <- Embeddings(data, 'harmony.lsi')
harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = data, reduction = "harmony.lsi", pt.size = .1, group.by = "sample") 
p2 <- VlnPlot(object = data, features = "harmony_1", group.by = "sample", pt.size = .1)
plot_grid(p1,p2)


DepthCor(data, reduction = "harmony.lsi" ) #see if some LSIs correlate with sequence depth...
ElbowPlot(data, reduction = "harmony.lsi")

data <- data %>% 
  RunUMAP(reduction = "harmony.lsi", dims = 2:20, n.neighbors = 30L,  min.dist = 0.01 ) %>% 
  FindNeighbors(reduction = "harmony.lsi", dims = 2:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## ATAC gene activity score 
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
DefaultAssay(data) <- "peaks"
gene.activities <- GeneActivity(data,
                                extend.upstream = 400
)

data[['atacRNA_400bp']] <- CreateAssayObject(counts = gene.activities)
data <- NormalizeData(
  object = data,
  assay = 'atacRNA_400bp',
  normalization.method = 'LogNormalize',
  scale.factor = median(data$nCount_RNA)
)

saveRDS(data, paste0(out.dir.seurat, "combined_filtered.rds"))

#gene expression vs gene activity 
DefaultAssay(data) <- "atacRNA_400bp"

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
      reduction = "umap",
      cols = c("lightgrey","darkgreen")
    )
    q1 <- q1 & scale_colour_gradientn(
      colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
      na.value = "lightgray"
    ) & NoAxes()
    
    DefaultAssay(data_dc3000) <- "atacRNA_400bp"
    q2 <- FeaturePlot(data_dc3000, 
                      features = gene_list, 
                      split.by = "sample", 
                      pt.size = 1, 
                      order=TRUE,
                      min.cutoff = 'q1',
                      max.cutoff = 'q99',
                      reduction = "umap",
                      cols = c("lightgrey","darkmagenta")
    )
    q2 <- q2 & scale_colour_gradientn(
      colours =c("lightgray",  "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D" ,"#238443" ,"#006837", "#004529"),
      na.value = "lightgray"
    )& NoAxes()
    p1 <- ggarrange(q1, q2, nrow  =2)
    ggsave(filename = paste("RNA_ATACactivity_plot/" ,fname,"_feture_dc3000.pdf", sep = ""), width = 16, height = 8, p1) 
  }, error=function(e){})
  
  tryCatch({
    
    DefaultAssay(data_avrrpt2) <- "SCT"
    
    q1 <- FeaturePlot(
      object = data_avrrpt2,
      features = gene_list,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = 1,
      split.by = "sample",
      order = TRUE,
      reduction = "umap",
      cols = c("lightgrey","darkgreen")
    )
    q1 <- q1 & scale_colour_gradientn(
      colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
      na.value = "lightgray"
    ) & NoAxes()
    
    DefaultAssay(data_avrrpt2) <- "atacRNA_400bp"
    q2 <- FeaturePlot(data_avrrpt2, 
                      features = gene_list, 
                      split.by = "sample", 
                      pt.size = 1, 
                      order=TRUE,
                      min.cutoff = 'q1',
                      max.cutoff = 'q99',
                      reduction = "umap",
                      cols = c("lightgrey","darkmagenta")
    )
    q2 <- q2 & scale_colour_gradientn(
      colours =c("lightgray",  "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D" ,"#238443" ,"#006837", "#004529"),
      na.value = "lightgray"
    )& NoAxes()
    p1 <- ggarrange(q1, q2, nrow  =2)
    ggsave(filename = paste("RNA_ATACactivity_plot/" , fname,"_feture_avrrpt2.pdf", sep = ""), width = 16, height = 8, p1) 
  }, error=function(e){})
  
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
      reduction = "umap",
      cols = c("lightgrey","darkgreen")
    )
    q1 <- q1 & scale_colour_gradientn(
      colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
      na.value = "lightgray"
    ) & NoAxes()
    
    DefaultAssay(data_avrrpm1) <- "atacRNA_400bp"
    q2 <- FeaturePlot(data_avrrpm1, 
                      features = gene_list, 
                      split.by = "sample", 
                      pt.size = 1, 
                      order=TRUE,
                      min.cutoff = 'q1',
                      max.cutoff = 'q99',
                      reduction = "umap",
                      cols = c("lightgrey","darkmagenta")
    )
    q2 <- q2 & scale_colour_gradientn(
      colours =c("lightgray",  "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D" ,"#238443" ,"#006837", "#004529"),
      na.value = "lightgray"
    )& NoAxes()
    p1 <- ggarrange(q1, q2, nrow  =2)
    ggsave(filename = paste("RNA_ATACactivity_plot/" , fname,"_feture_avrrpm1.pdf", sep = ""), width = 16, height = 8, p1) 
  }, error=function(e){})
  
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
      reduction = "umap",
      cols = c("lightgrey","darkgreen")
    )
    q1 <- q1 & scale_colour_gradientn(
      colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
      na.value = "lightgray"
    ) & NoAxes()
    
    DefaultAssay(data) <- "atacRNA_400bp"
    q2 <- FeaturePlot(data, 
                      features = gene_list, 
                      # split.by = "sample", 
                      pt.size = 1, 
                      order=TRUE,
                      # min.cutoff = 'q1',
                      max.cutoff = 'q99',
                      min.cutoff = 0.3,
                      reduction = "umap",
                      cols = c("lightgrey","darkmagenta")
    )
    q2 <- q2 & scale_colour_gradientn(
      colours =c("lightgray",  "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D" ,"#238443" ,"#006837", "#004529"),
      na.value = "lightgray"
    )& NoAxes()
    p1 <- ggarrange(q1, q2, nrow  =2)
    ggsave(filename = paste("RNA_ATACactivity_plot/" , fname,"_feture.pdf", sep = ""), width = 4, height = 8, p1) 
  }, error=function(e){})
  
}

sct.atac.coplot("SWEET12", "SWEET12")

#marker genes based on ATAC activity 
marker.atac.activity <- FindAllMarkers(
  object = data,
  min.pct = 0.05,
  test.use = 'wilcox',
  only.pos = TRUE,
  logfc.threshold = 0.2
)
#giving gene_id for gene_names
gene_id <- annotation$gene_id[match(marker.atac.activity$gene, annotation$gene_name)]
gene_id[is.na(gene_id)] <-  marker.atac.activity$gene[is.na(gene_id)]
marker.atac.activity$gene_id <- gene_id
write.table(marker.atac.activity, file =  "markers_peaks.txt", row.names=T, col.names=NA, sep="\t", quote=F)


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## coverage plots
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
region_name <- "2-11172834-11173492"

coverage.peak <- function(region_name){
  q <- CoveragePlot(data2, region = region_name,
                    #features = data@assays$peaks@links$gene[i],
                    assay = "peaks",
                    annotation = TRUE,
                    extend.upstream = 4000,
                    extend.downstream =4000, 
                    heights = c(5,.5,.5),
                    group.by = "SCT_snn_res.1"
  )
  q <-q & scale_fill_manual(values = cols) #change colors 
  ggsave(filename = paste0("coveragePlot/", region_name, ".pdf") ,width = 5, height = 5,q)
}

coverage.peak(region_name)


##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## joint embedding
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

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

DepthCor(data, reduction = "wnn.umap" ) #see if some LSIs correlate with sequence depth...
ElbowPlot(data, reduction = "wnn.umap")

data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)
