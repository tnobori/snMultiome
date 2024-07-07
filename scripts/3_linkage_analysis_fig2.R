##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## This script integrates snRNA-seq and snATAC-seq to analyze ACR-gene links
## Figures produced with this script: Fig2, FigS4

## A fully processed Seurat object can be downloaded at http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/processed_seurat_object/combined_filtered.rds

## Contact: Tatsuya Nobori (tatsuyanobori@gmail.com)
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# Clearing the workspace
rm(list = ls())

# Loading project configuration
source("http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/scripts/_config_multiome.R")

## Assuming you are in the directory where this and other R scripts are stored.
mainDir <- file.path(getwd(), "data_out")
subDir <- "3_linkage"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Functions
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

coverage.plot2 <- function(data, target, size, size2=size ,out.dir){
  tryCatch({data <- LinkPeaks(
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
    extend.downstream = size2, 
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
    #idents = "04h",
    extend.upstream = size,
    extend.downstream = size2, 
    group.by = "sample.order2",
    split.assays = TRUE, 
    heights = c(2,2,2,.5)
    
    
  )
  q <-q & scale_fill_manual(values = col.sample) #change colors 
  ggsave(filename = paste0(out.dir,"/coveragePlot_marker_", target , "_", size, "bp_noRep.pdf"), width = 20, height = 10, q)
  
  q <- CoveragePlot(
    object = data,
    region = target,
    features = target,
    expression.assay = "SCT",
    #idents = idents.plot,
    extend.upstream = size,
    extend.downstream =size2, 
    group.by = "SCT_snn_res.1",
    #split.assays = TRUE,
    heights = c(2,2,2,.5)
  )
  q <-q & scale_fill_manual(values = cols) #change colors 
  ggsave(filename = paste0(out.dir,"/coveragePlot_marker_", target, "_", size ,"bp_RNAcluster.pdf"), width = 20, height = 10, q)
  }, error=function(e){})
}


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
  legend.key = element_rect(fill = NA, size = 1)
)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
# Loading seurat object
data <- readRDS(paste(mainDir, "_seurat_object","combined_filtered.rds", sep = "/"))

# Load annotation
annotation <- Annotation(data@assays$ATAC)
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Link calling
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
DefaultAssay(data) <- "peaks"
seqnames(BSgenome.Athaliana.TAIR.TAIR9) <- gsub("^Chr", "", seqnames(BSgenome.Athaliana.TAIR.TAIR9))

# Computing GC content for each peak
data <- RegionStats(data, genome = BSgenome.Athaliana.TAIR.TAIR9)

### each strain separately
cell_kt56 <- rownames(data@meta.data)[ grepl("DC3000", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ] 
cell_kt57 <- rownames(data@meta.data)[ grepl("AvrRpt2", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ]
cell_kt58 <- rownames(data@meta.data)[ grepl("AvrRpm1", data@meta.data[, "sample"]) | grepl("Mock", data@meta.data[, "sample"]) ]
cell_mock <- rownames(data@meta.data)[grepl("Mock", data@meta.data[, "sample"]) ]

data_kt56 <- subset(x = data, cells = cell_kt56)
data_kt57 <- subset(x = data, cells = cell_kt57)
data_kt58 <- subset(x = data, cells = cell_kt58)
data_mock <- subset(x = data, cells = cell_mock)

# Setting the default assay for each subset
for (data_subset in list(data, data_kt56, data_kt57, data_kt58, data_mock)) {
  DefaultAssay(data_subset) <- "peaks"
}

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
# Genome-wide analysis of linkages for each condition; this takes time to finish
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
linkPeaks <- function(data_object, file_suffix) {
  data_object <- LinkPeaks(
    object = data_object,
    peak.assay = "peaks",
    expression.assay = "SCT"
  )
  links_all <- data_object@assays$peaks@links
  write.table(links_all, file = paste0("linkage_", file_suffix, "_ALL.txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
}

# These might take time to complete
linkPeaks(data_mock, "mock")
linkPeaks(data_kt56, "kt56")
linkPeaks(data_kt57, "kt57")
linkPeaks(data_kt58, "kt58")
linkPeaks(data, "combined")

# Save the updated Seurat object
out.dir.seurat <- file.path(mainDir, "/_seurat_object")
saveRDS(data, paste0(out.dir.seurat, "/combined_filtered.rds"))

## Loading linkage data
links_mock_all <- read.delim(file = "linkage_mock_ALL.txt", head = TRUE, row.names = 1)
links_kt56_all <- read.delim(file = "linkage_kt56_ALL.txt", head = TRUE, row.names = 1)
links_kt57_all <- read.delim(file = "linkage_kt57_ALL.txt", head = TRUE, row.names = 1)
links_kt58_all <- read.delim(file = "linkage_kt58_ALL.txt", head = TRUE, row.names = 1)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Linkage meta analysis
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
dir.create("stats_basics", showWarnings = FALSE)

# Add sample labels to linkage data
links_mock_all$sample <- "Mock"
links_kt56_all$sample <- "DC3000"
links_kt57_all$sample <- "AvrRpt2"
links_kt58_all$sample <- "AvrRpm1"

#adding peak size
peaksize <- function(input){
  link_split <- str_split(input$peak, "-")
  peaksize <- lapply(link_split, function(x){as.numeric(x[3])-as.numeric(x[2])})
  peaksize %>% unlist()
}
links_mock_all$peaksize <- peaksize(links_mock_all)
links_kt56_all$peaksize <- peaksize(links_kt56_all)
links_kt57_all$peaksize <- peaksize(links_kt57_all)
links_kt58_all$peaksize <- peaksize(links_kt58_all)


# Filter links based on score threshold
score.threshold <- 0.05
filterLinks <- function(link_data) {
  link_data[link_data$score > score.threshold, ]
}

top.links_mock_all <- filterLinks(links_mock_all)
top.links_kt56_all <- filterLinks(links_kt56_all)
top.links_kt57_all <- filterLinks(links_kt57_all)
top.links_kt58_all <- filterLinks(links_kt58_all)

# Combining filtered links
indata <- rbind(top.links_kt56_all, top.links_kt57_all, top.links_kt58_all, top.links_mock_all)

# Plotting
q <- ggplot(indata, aes(x=width, y=score, color=sample)) + geom_point(size=.2) + xlim(0, 5000) + optns
q <- q + scale_color_manual(values= c( "#af1919","#006666",  "#a2790d", "#005b96"))
ggsave("stats_basics/figS4a.pdf", height = 5, width = 8, q)

q <- ggplot(indata, aes(x=width, color=sample)) + geom_density(size=1.5)+ xlim(0, 10000) + optns
q <- q + scale_color_manual(values= c( "#af1919","#006666",  "#a2790d", "#005b96"))
ggsave("stats_basics/figS4c.pdf", height = 5, width = 8, q)

q <- ggplot(indata, aes(x=width, color=sample)) + geom_density(size=1.5)+ xlim(0, 1000) + optns
q <- q + scale_color_manual(values= c( "#af1919","#006666",  "#a2790d", "#005b96"))
ggsave("stats_basics/fig2b.pdf", height = 5, width = 8, q)

##analyze peak size for all peaks; compare with non-linked peaks
linked_peaksize <- indata[, c("peak", "peaksize")] %>% unique()
linked_peaksize$link <- "linked"
colnames(linked_peaksize)[1] <- "peak"

allpeaks <- paste(data@assays$peaks@ranges@seqnames, data@assays$peaks@ranges@ranges, sep = "-")
nonlinked_peaks <- allpeaks[!allpeaks %in% linked_peaksize$peak]

link_split <- str_split(nonlinked_peaks, "-")
peaksize <- lapply(link_split, function(x){as.numeric(x[3])-as.numeric(x[2])})
nonlinked_peaksize <- peaksize %>% unlist() %>% as.data.frame()
nonlinked_peaksize <- cbind(nonlinked_peaks, nonlinked_peaksize)
colnames(nonlinked_peaksize)[2] <- "peaksize"
nonlinked_peaksize$link <- "nonlinked"
colnames(nonlinked_peaksize)[1] <- "peak"

peaksize_comb <- rbind(linked_peaksize, nonlinked_peaksize)


q <- ggplot(peaksize_comb, aes(x=peaksize, color=link)) + geom_density(size=1.5) + optns
q <- q + scale_color_manual(values= c( "#af1919","#006666",  "#a2790d", "#005b96"))
ggsave("stats_basics/figS4d.pdf", height = 5, width = 8, q)


## distance from TSS for marker genes?
markers.rna <- read.delim(file.path(mainDir, "2_clustering/RNA" ,"markers_RNA.txt"), header = T) #Table S1
markers.rna <- markers.rna$gene %>% unique()

indata$marker <- ifelse(indata$gene %in% markers.rna, "Marker", "Non-marker")

# To specifically analyze distances within a certain range, for example, within 5kb from TSS
indata.sel <- indata %>% filter(width <= 5000)

q <- ggplot(indata.sel, aes(x=width, color=marker)) + geom_density(size=1.5)+  optns
q <- q + scale_color_manual(values= c( "#af1919","#006666",  "#a2790d", "#005b96"))
ggsave(paste0("stats_basics/figS2c.pdf"), height = 5, width = 8, q)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## MaxPearson analysis for all samples 
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
# Thresholds for filtering
score.threshod <- 0.1
width.threshod <- 500000

links_kt56_all <- read.delim(file = "linkage_kt56_ALL.txt", head = TRUE, row.names = 1)
links_kt57_all <- read.delim(file = "linkage_kt57_ALL.txt", head = TRUE, row.names = 1)
links_kt58_all <- read.delim(file = "linkage_kt58_ALL.txt", head = TRUE, row.names = 1)
links_mock_all <- read.delim(file = "linkage_mock_ALL.txt", head = TRUE, row.names = 1)

top.links_mock_all <- links_mock_all[links_mock_all$score > score.threshod ,]
top.links_kt56_all <- links_kt56_all[links_kt56_all$score > score.threshod ,]
top.links_kt57_all <- links_kt57_all[links_kt57_all$score > score.threshod ,]
top.links_kt58_all <- links_kt58_all[links_kt58_all$score > score.threshod ,]

top.links_mock_all <- top.links_mock_all[top.links_mock_all$width < width.threshod ,]
top.links_kt56_all <- top.links_kt56_all[top.links_kt56_all$width < width.threshod ,]
top.links_kt57_all <- top.links_kt57_all[top.links_kt57_all$width < width.threshod ,]
top.links_kt58_all <- top.links_kt58_all[top.links_kt58_all$width < width.threshod ,]

make.link.table <- function(indata){
nlinks <- table(indata$gene) %>% as.data.frame()
maxPearson <- indata %>% group_by(gene) %>% summarise_at(vars(score), list(name = max))%>% as.data.frame()
nlinks$maxPearson <- maxPearson$name
nlinks
}

link_table_mock <- make.link.table(top.links_mock_all)
link_table_kt56 <- make.link.table(top.links_kt56_all)
link_table_kt57 <- make.link.table(top.links_kt57_all)
link_table_kt58 <- make.link.table(top.links_kt58_all)

gene.list <- c(as.character(link_table_kt56$Var1),as.character(link_table_kt57$Var1),as.character(link_table_kt58$Var1),as.character(link_table_mock$Var1) ) %>% unique
link.table.comb <- matrix(data = 0, ncol = 4, nrow = length(gene.list))
colnames(link.table.comb) <- c("Mock", "DC3000", "AvrRpt2", "AvrRpm1")
rownames(link.table.comb) <- gene.list

link.table.comb[as.character(link_table_mock$Var1) ,1] <- link_table_mock$maxPearson
link.table.comb[as.character(link_table_kt56$Var1) ,2] <- link_table_kt56$maxPearson
link.table.comb[as.character(link_table_kt57$Var1) ,3] <- link_table_kt57$maxPearson
link.table.comb[as.character(link_table_kt58$Var1) ,4] <- link_table_kt58$maxPearson
link.table.comb <- link.table.comb %>% as.data.frame()
write.table(link.table.comb, file = paste0("link_table_maxPearson_score", score.threshod, "_width", width.threshod,".txt"), row.names=T, col.names=NA, sep="\t", quote=F)


link.table.comb.freq <- matrix(data = 0, ncol = 4, nrow = length(gene.list))
colnames(link.table.comb.freq) <- c("Mock", "DC3000", "AvrRpt2", "AvrRpm1")
rownames(link.table.comb.freq) <- gene.list
link.table.comb.freq[as.character(link_table_mock$Var1) ,1] <- link_table_mock$Freq
link.table.comb.freq[as.character(link_table_kt56$Var1) ,2] <- link_table_kt56$Freq
link.table.comb.freq[as.character(link_table_kt57$Var1) ,3] <- link_table_kt57$Freq
link.table.comb.freq[as.character(link_table_kt58$Var1) ,4] <- link_table_kt58$Freq
link.table.comb.freq <- link.table.comb.freq %>% as.data.frame()
write.table(link.table.comb.freq, file = paste0("link_table_Nlinks_", score.threshod, "_width", width.threshod,".txt"), row.names=T, col.names=NA, sep="\t", quote=F)

library(viridis)
pdf(paste0("fig2c_noKmean.pdf"))
pheatmap(link.table.comb,
         cluster_cols=FALSE,
         cluster_rows = TRUE,
         cellwidth = 20,
         color = viridis(100),
         breaks= seq(0, .3, length.out = 101),
         fontsize_row = 6,
         fontsize_col = 10
)
dev.off()

#Kmean clustering

input <- link.table.comb

# determining an optimal K value 
png(filename = "kplot_link.table.comb.png", width = 2000, height = 2000)
kplot(input)
dev.off()

n.cluster_gene <- 8 #number of clusters

input_clst <- kclust(input, n.cluster_gene)
an_row_km <- input_clst[, ncol(input_clst)]
names(an_row_km) <- rownames(input_clst)
an_row_km <- as.data.frame(an_row_km)
an_row_km[, 1] <- as.character(an_row_km[, 1])

input_clst2 <- input_clst[, -ncol(input_clst)]
labels_row <- rownames(input_clst2)
labs_row <- rep("", times = length(labels_row))
gene.sel <- c("ALD1",  "BCA2", "FDH", "MAM1", "WRKY75", "WRKY46", "VSP1", "LOX2", "AIG1", "NAC019", "CYP71A13", "CBP60G", "SWEET12",  "BSMT1", "VSP2", "ACD6", "STP13", "UGT85A1")
labs_row[match(gene.sel, labels_row)] <- gene.sel

an_color <- cbind(c(1:n.cluster_gene),  brewer.pal(n=n.cluster_gene, name="Dark2")) %>% as.data.frame()
colnames(an_color)[2] <- "an_row_km"
an_color <- list(an_row_km = c("1" = "#1B9E77",
                               "2" = "#D95F02",
                               "3" = "#7570B3",
                               "4" = "#E7298A",
                               "5" = "#66A61E",
                               "6" = "#E6AB02",
                               "7" = "#A6761D",
                               "8" = "#011b4b"
                               ))

pdf(paste0("fig2c.pdf"))

pheatmap(input_clst2,
         cluster_cols=TRUE,
         cluster_rows = FALSE,
         cellwidth = 20,
         color = viridis(50),
         breaks= seq(0, .3, length.out = 50),
         annotation_row = an_row_km,
         annotation_colors = an_color,
         fontsize_row = 6,
         fontsize_col = 10,
         labels_row = labs_row
)
dev.off()

##GO, with new database
# annotation <- Annotation(data)
input_clst$geneID <- sapply(rownames(input_clst), genes_to_id) #gene names to ID
for (i in c(1:n.cluster_gene)){
  gene <- input_clst$geneID[input_clst$kmean==i]
  
  ego <- enrichGO(gene          = gene,
                  #universe      = names(geneList),
                  OrgDb         = org.Athaliana.eg.db,
                  ont           = "BP",
                  keyType = 'GID',
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  pdf(paste("GO_linkage_maxPearson/clst_", i,".pdf", sep = ""), 
      # width = 1000, height = 1000
  )
  try(print(barplot(ego, showCategory=10)))
  dev.off()
  write.table(ego@result, file =  paste0("GO_linkage_maxPearson/clst_", i,".txt"), row.names=T, col.names=NA, sep="\t", quote=F)
  
}

#GO analysis summary plot
#Taking top 2 GO terms per cluster for visualization
clst.sel <- c(1:n.cluster_gene ) #clusters to show
out <- matrix(data = NA, nrow = length(clst.sel)*2, ncol = 3) %>% as.data.frame()
colnames(out) <- c("Description", "Size", "adjusted pvalue")

optns3 <- theme(
  axis.text.x = element_text(margin = margin(c(10, 0, 0, 0)), size = 15, face = "bold", family = "Helvetica", angle = 0),
  axis.text.y = element_text(margin = margin(c(0,10, 0, 0)) ,size = 15, face = "bold", family = "Helvetica"), 
  axis.title = element_text(size = 10, face = "bold", family = "Helvetica"),
  axis.title.x = element_text(margin = margin(c(10, 0, 0, 0))),
  axis.title.y = element_text(margin = margin(c(0, 10, 0, 0))),
  axis.ticks.y = element_line(size = 1),
  axis.ticks.x = element_line(size = 1),
  axis.ticks.length = unit(.25, "cm"),
  #axis.ticks.margin = unit(1, "pt"),
  axis.line  = element_line(size = 1),
  # panel.background = element_blank(),
  panel.background = element_rect(fill = "white"),
  #panel.grid.minor = element_line(color = 'gray', size = .2, linetype = "dashed"),
  panel.grid.major = element_line(color = 'gray', size = .3, linetype = "dotted"),
  # legend.text = element_text(size = 10, face = "bold", family = "Helvetica", margin = margin(c(20, 0, 20, 0))),
  # legend.title = element_text(size = 10, face = "bold", family = "Helvetica"),
  # # # legend.key = element_rect(fill = NA , size= 3),
  # legend.key.size = unit(1, 'cm'),
  # legend.key.height  = unit(1, 'cm'),
  # legend.key.width  = unit(1, 'cm'),
  # #panel.background = element_rect(fill = "black"), # bg of the panel
  #plot.background = element_rect(fill = "black", color = NA), # bg of the plot
  #legend.background = element_rect(fill = "black"), # get rid of legend bg
  #legend.box.background = element_rect(fill = "black") # get rid of legend panel bg
  legend.box.background = element_blank(), # get rid of legend panel bg,
  legend.background = element_blank() # get rid of legend bg
) 

for (i in c(1:length(clst.sel))){
  go.data <- read.delim(file = paste0("GO_linkage_maxPearson/clst_", clst.sel[i],".txt"), header = T, row.names = 1)
  
  s <- 2*(i-1)+1
  
  out[s, ] <- c(go.data$Description[1],as.numeric(go.data$Count[1]), as.numeric(go.data$p.adjust[1]))
  out[s+1, ] <- c(go.data$Description[2],as.numeric(go.data$Count[2]), as.numeric(go.data$p.adjust[2]))
  
}

out$clst <- as.character(rep(clst.sel, each = 2))
#Then turn it back into a factor with the levels in the correct order
out$clst <- factor(out$clst, levels=unique(out$clst))

out$log.p <- -log10(as.numeric(out$`adjusted pvalue`))
out$Size <- as.numeric(out$Size)

q <- ggplot(out, aes(x = clst, y = Description, size = Size , color = log.p)) + geom_point()+ optns3 + 
  scale_color_gradient(low = "grey",high = "#418557",limits = c(0,15))+
  scale_size(range = c(0,10)) 
ggsave("GO_linkage_maxPearson/figS4b.pdf", height = 6, width = 10 ,q)





##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Coverage plots (such as Fig2e)
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==



col.sample <- c("00_00_Mock_rep1" = "#005b96",
                "00_00_Mock_rep2" = "#005b96",
                "00_00_Mock" = "#005b96",
                "01_DC3000_04h" = "#d9c99e",
                "01_DC3000_06h" = "#bda155",
                "01_DC3000_09h" = "#a2790d",
                "01_DC3000_24h"= "#513c06",
                "02_AvrRpt2_04h"= "#66b2b2",
                "02_AvrRpt2_06h"= "#008080",
                "02_AvrRpt2_09h_rep1"= "#006666",
                "02_AvrRpt2_24h"= "#004c4c",
                "02_AvrRpt2_09h_rep2"= "#006666",
                "02_AvrRpt2_09h"= "#006666",
                "03_AvrRpm1_04h"= "#ca6666",
                "03_AvrRpm1_06h"= "#af1919",
                "03_AvrRpm1_09h"= "#740000",
                "03_AvrRpm1_24h"= "#420000"
                
)


coverage.plot2(data, "CBP60G", 3000,"coveragePlot")# Fig2e
coverage.plot2(data, "EDS16", 3000,"coveragePlot")# FigS1e
coverage.plot2(data, "ACT2", 3000,"coveragePlot")# FigS1e

########################################################################################################################
##find motifs from Linked or non-linked genes
########################################################################################################################
markers.rna <- read.delim("markers_RNA.txt", header = T) # Table S1
defense.markers.rna <- markers.rna[markers.rna$cluster==3|markers.rna$cluster==4|markers.rna$cluster==7|markers.rna$cluster==11|markers.rna$cluster==29|markers.rna$cluster==24|markers.rna$cluster==12,] #selecting immune-active clusters
defense.markers.rna2 <- defense.markers.rna$gene %>% unique()

linked_genes <- indata$gene[indata$width < 3000] %>% unique()
defense.markers.rna2 <- defense.markers.rna2 %>% as.data.frame() 
colnames(defense.markers.rna2) <- "rowname"
defense.markers.rna2$link <- "non-Linked"
defense.markers.rna2$link[defense.markers.rna2$rowname %in% linked_genes ] <- "Linked <3kb"

DefaultAssay(data) <- "peaks"
nolink.sel <- rownames(defense.markers.rna2)[defense.markers.rna2=="non-Linked"]
link.sel <- rownames(defense.markers.rna2)[defense.markers.rna2=="Linked <3kb"]

dir.create("fig2GH", showWarnings = FALSE)
write.table(nolink.sel, file = paste0("fig2GH/noLink_defense_genes.txt"), row.names=T, col.names=NA, sep="\t", quote=F)
write.table(link.sel, file = paste0("fig2GH/Link_defense_genes.txt"), row.names=T, col.names=NA, sep="\t", quote=F)

#analyze 3kb upstream of linked or nonlinked defense genes
#find peaks associated with genes of interests...
ranges.nolink <- promoters(annotation, upstream = 3000)[annotation$gene_name %in% nolink.sel, ] 
ranges.link <- promoters(annotation, upstream = 3000)[annotation$gene_name %in% link.sel, ] 

peaks.data <- StringToGRanges(rownames(data))

#overlap selected promoter regions with peaks in my data
peak.nolink <- findOverlaps(query = ranges.nolink, 
                            subject = peaks.data,
                            type = 'any',
                            select = 'all') %>% as.data.frame()

peak.nolink.ranges <- rownames(data)[peak.nolink$subjectHits]

peak.link <- findOverlaps(query = ranges.link, 
                          subject = peaks.data,
                          type = 'any',
                          select = 'all') %>% as.data.frame()

peak.link.ranges <- rownames(data)[peak.link$subjectHits]

#motif enrichment 
enriched.motifs <- FindMotifs(
  object = data,
  features = peak.nolink.ranges,
)

q <- MotifPlot(
  object = data,
  motifs = head(rownames(enriched.motifs))
)
ggsave(filename = paste0("fig2GH/fig2h.pdf"), width = 10, height = 5, q)


#motif enrichment 
enriched.motifs <- FindMotifs(
  object = data,
  features = peak.link.ranges,
)

q <- MotifPlot(
  object = data,
  motifs = head(rownames(enriched.motifs)
  )
)
ggsave(filename = paste0("nonlinked_fig4related/fig2g.pdf"), width = 10, height = 5, q)













#delete later 
annotation$gene_name <- gsub("_", "-", annotation$gene_name) # Seurat does not accept "_" in gene names, and it replaces them with "-". So, the annotation should match these changes.

# In the GTF file, there are cases where a row has two gene_name columns. An example is LRR XI-23 and several others. Actual names seem to be used for cellranger output, but IDs are loaded as annotation.

an_comb <- c(annotation$gene_id, annotation$gene_name) %>% na.omit() %>% unique()
x = rownames(data)[!(rownames(data) %in% an_comb)]
xx <- gsub("\\.[0-9]" ,"" ,x) %>% unique()# remove .1 to match gene names..

annotation[annotation$gene_name %in% xx ,]


annotation$gene_name_orig <- annotation$gene_name


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


for (dupgene in xx){
  annotation$gene_name[annotation$gene_name %in% dupgene &  annotation$type=="body"] <- makeUniqueNames(annotation$gene_name[annotation$gene_name %in% dupgene  &  annotation$type=="body"])
  # print(makeUniqueNames(annotation$gene_name[annotation$gene_name %in% dupgene &  annotation$type=="body"]))
}

annotation$gene_name[annotation$gene_name %in% "LHCB1.1"  &  annotation$type=="body"] <- makeUniqueNames(annotation$gene_name[annotation$gene_name %in% "LHCB1.1"  &  annotation$type=="body"])
