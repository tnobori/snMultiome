##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## This script analyzes bulk RNA-seq data
## Figures produced with this script: Fig6
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

rm(list = ls())
source("http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/scripts/_config_multiome.R")

mainDir <- "/data_out"
subDir <- "6_bulkRNAseq"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

setwd(file.path(mainDir, subDir))

# Load dataset
data <- readRDS(paste(mainDir, "_seurat_object","combined_filtered.rds", sep = "/"))## A fully processed Seurat object can be downloaded at http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/processed_seurat_object/combined_filtered.rds. Or you can create one from scratch by running 1_qc_data_integration_figS1.R

id_to_gene <- function(genes){
  if (length(annotation$gene_name[annotation$gene_id==genes] %>% na.omit)!=0 ){
    tryCatch({id <- annotation$gene_name[annotation$gene_id==genes] %>% na.omit
    id[1] }, error=function(e){})
  }else{
    genes
  }
}

gt3aox <- read.delim("http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/misc/fitted_mean_Col_GT3a_DC_q0.05_2fold.txt", header = T, row.names = 1)
down <- rownames(gt3aox)[((gt3aox[, 4] - gt3aox[, 3]) - (gt3aox[, 2] - gt3aox[, 1])) < -1]

up <- rownames(gt3aox)[(gt3aox[, 3] - gt3aox[, 1]) >1]

marker_sel <- out_sel_markers$gene[out_sel_markers$cluster==4]

down_name <- c()
for (gene in down){
  down_name <- c(down_name, id_to_gene(gene))
}


up_name <- c()
for (gene in up){
  up_name <- c(up_name, id_to_gene(gene))
}

marker_sel %in% down_name
marker_sel %in% up_name


# heatmap
input <- gt3aox
png(filename = "6_bulkRNAseq/kplot.png", width = 2000, height = 2000)
kplot(input) #K = 10
dev.off()

# Perform k-means clustering on genes
set.seed(123) # Set seed for reproducibility

# Scale the rows before computing the hierarchical clustering
input <- t(scale(t(gt3aox)))

# Perform hierarchical clustering
hc <- hclust(dist(input), method = "complete")

# Cut the tree into a specified number of clusters
# Replace `k` with your desired number of clusters
k <- 10
clusters <- cutree(hc, k)

# Create a data frame with cluster assignments
cluster_assignments <- data.frame(Cluster = as.factor(clusters))
rownames(cluster_assignments) <- rownames(input)

# Extract gene names for each cluster
gene_clusters <- by(rownames(input), clusters, function(x) x)

#renaming clusters to order correctly 
x <- cluster_assignments[hc$order ,1] #this is current order. wnat to make it 1-k
f <- factor(x, levels = unique(x))
y <- as.integer(f)
cluster_assignments[hc$order ,1]  <- y

# Create a heatmap with the hierarchical clustering and clusters annotation
pdf(paste("6_bulkRNAseq/gt3aOX_mock_DC3000_heatmap.pdf", sep = ""))

pheatmap(input,
         annotation_row = cluster_assignments,
         show_rownames = TRUE,
         show_colnames = TRUE,
         # clustering_method = "complete",
         color = rev(colorRampPalette(brewer.pal(n = 11, name ="RdYlBu"))(100)),
         breaks= seq(-1.5, 1.5, length.out = 101),
         scale = 'row',
         fontsize_row = 0.0001,
         cluster_rows = hc,
         cluster_cols = FALSE
)

dev.off()



#go enrich on akiras rna-seq data

for (i in c(1:k)){
  gene <- rownames(cluster_assignments)[cluster_assignments$Cluster==i]
  
  ego <- enrichGO(gene          = gene,
                  #universe      = names(geneList),
                  OrgDb         = org.Athaliana.eg.db,
                  ont           = "BP",
                  keyType = 'GID',
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  png(paste("6_bulkRNAseq/GO_clst_", i,".png", sep = ""), 
      # width = 1000, height = 1000
  )
  try(print(barplot(ego, showCategory=10)))
  dev.off()
  
  write.table(ego@result, file= paste0("6_bulkRNAseq/go_clst_", i, ".txt"), row.names=T, col.names=NA, sep="\t", quote=F)
  
}

#GO enrichment summary 

out <- matrix(data = NA, nrow =k*3, ncol = 3) %>% as.data.frame()
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

for (i in c(1:k)){
  go.data <- read.delim(file = paste0("6_bulkRNAseq/GO_clst_",i,".txt"), header = T, row.names = 1)
  
  s <- (i-1)*3+1
  
  out[s, ] <- c(go.data$Description[1],as.numeric(go.data$Count[1]), as.numeric(go.data$p.adjust[1]))
  out[s+1, ] <- c(go.data$Description[2],as.numeric(go.data$Count[2]), as.numeric(go.data$p.adjust[2]))
  out[s+2, ] <- c(go.data$Description[3],as.numeric(go.data$Count[3]), as.numeric(go.data$p.adjust[3]))
  
}

out$clst <- as.character(rep(c(1:k), each = 3))
#Then turn it back into a factor with the levels in the correct order
out$clst <- factor(out$clst, levels=unique(out$clst))

out$log.p <- -log10(as.numeric(out$`adjusted pvalue`))
out$Size <- as.numeric(out$Size)

q <- ggplot(out, aes(x = clst, y = Description, size = Size , color = log.p)) + geom_point()+ optns3 + 
  scale_color_gradient(low = "grey",high = "#418557",limits = c(0,51))+
  scale_size(range = c(0,10)) 
ggsave("6_bulkRNAseq/fig6i.pdf", height = 6, width = 8 ,q)



