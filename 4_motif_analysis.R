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
## Plotting parameters
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
          "24" = "#fc0ed5"
          
)

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
## data loading
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

data <- readRDS(".../combined_filtered.rds")


link_dir <- "/Users/tatsuyanobori/Dropbox/SALK_clowd/Projects/SA_PTI_ETI_single_cell/SA_039_85_multiome_publication/out/3_linkage/"
links_mock_all <- read.delim(file = paste0(link_dir ,"linkage_mock_ALL.txt"), head = TRUE, row.names = 1)
links_dc3000_all <- read.delim(file = paste0(link_dir ,"linkage_dc3000_ALL.txt"), head = TRUE, row.names = 1)
links_avrrpt2_all <- read.delim(file = paste0(link_dir ,"linkage_avrrpt2_ALL.txt"), head = TRUE, row.names = 1)
links_avrrpm1_all <- read.delim(file = paste0(link_dir ,"linkage_avrrpm1_ALL.txt"), head = TRUE, row.names = 1)

links.comb <- rbind(links_mock_all, links_dc3000_all, links_avrrpt2_all, links_avrrpm1_all)
links.comb.significant <- links.comb[links.comb$score > 0.1 & links.comb$width < 5000 ,]

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## pre processing
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

DefaultAssay(data) <- "peaks"
seqnames(BSgenome.Athaliana.TAIR.TAIR9) <- gsub("^Chr", "", seqnames(BSgenome.Athaliana.TAIR.TAIR9))

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 3702, all_versions = FALSE)
)

# add motif information
data <- AddMotifs(data, genome = BSgenome.Athaliana.TAIR.TAIR9, pfm = pwm)

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## motif enrichment analysis
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

# motif enrichment analysis between cluster 1 and cluster 2
motif.analysis <- function(clst1, clst2){
  Idents(data) <- data$SCT_snn_res.1
  DefaultAssay(data) <- "peaks"
  out_dir <- paste0(".../clst", clst1 , "vs", clst2, "/")
  dir.create(out_dir, showWarnings = FALSE)
  
  print("motif identification...")
  
  da_peaks <- FindMarkers(
    object = data,
    ident.1 = clst1,
    ident.2 = clst2,
    only.pos = TRUE,
    test.use = 'LR',
    min.pct = 0.05,
    latent.vars = 'nCount_peaks'
  )
  
  print("motif enrichment 1/2")
  # test enrichment
  enriched.motifs <- FindMotifs(
    object = data,
    features = rownames(da_peaks)
  )
  write.table(enriched.motifs, file = paste0(out_dir ,"enriched_motif.txt"), row.names=T, col.names=NA, sep="\t", quote=F)
  
  
  q <- MotifPlot(
    object = data,
    motifs = head(rownames(enriched.motifs))
  )
  ggsave(paste0(out_dir, "enriched_motif.pdf"), height = 5, width = 10 , q)
  
  ##CREs
  #overlap marker peaks with CREs
  links.comb.marker <- da_peaks[rownames(da_peaks) %in%  links.comb.significant$peak  ,]
  
  print("motif enrichment 2/2")
  
  # test enrichment
  enriched.motifs.CRES <- FindMotifs(
    object = data,
    features = rownames(links.comb.marker)
  )
  write.table(enriched.motifs, file = paste0(out_dir ,"enriched_motif_CREs.txt"), row.names=T, col.names=NA, sep="\t", quote=F)
  
  q <- MotifPlot(
    object = data,
    motifs = head(rownames(enriched.motifs.CRES))
  )
  ggsave(paste0(out_dir, "enriched_motif_CREs.pdf"), height = 5, width = 10 , q)
  toc()
}

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## building chromVAR assay
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

motif.matrix <- GetMotifData(object = data, slot = "data")
peak.matrix <- GetAssayData(object = data, slot = "counts")
idx.keep <- rowSums(x = peak.matrix) > 0
peak.matrix <- peak.matrix[idx.keep, , drop = FALSE]
motif.matrix <- motif.matrix[idx.keep, , drop = FALSE]
peak.ranges <- granges(x = data)
peak.ranges <- peak.ranges[idx.keep]
chromvar.obj <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = peak.matrix),
  rowRanges = peak.ranges)
chromvar.obj <- chromVAR::addGCBias(
  object = chromvar.obj,
  genome = BSgenome.Athaliana.TAIR.TAIR9)

row.data <- data.frame(rowData(chromvar.obj))
row.data[is.na(row.data)] <- 0
rowData(chromvar.obj) <- row.data

bg <- chromVAR::getBackgroundPeaks(
  object = chromvar.obj)

dev <- chromVAR::computeDeviations(
  object = chromvar.obj,
  annotations = motif.matrix,
  background_peaks = bg
)

chromvar.z <- SummarizedExperiment::assays(dev)[[2]]
rownames(x = chromvar.z) <- colnames(x = motif.matrix)
obj <- CreateAssayObject(data = chromvar.z)
data[["chromvar"]] <- obj

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##  CROMVAR analysis
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

DefaultAssay(data) <- 'chromvar'

data@assays$chromvar@data[is.na(data@assays$chromvar@data)] <- 0

enriched.motifs <- read.delim(file = ".../clst2vs0/enriched_motif_CREs.txt", header = T, row.names = 1)

#make motif ID_gene name combined
enriched.motifs$combname <- paste(enriched.motifs$motif, enriched.motifs$motif.name, sep = "_")
rownames(data@assays$chromvar@data)[match(enriched.motifs$motif, rownames(data@assays$chromvar@data))] <- enriched.motifs$combname

saveRDS(data, "combined_filtered.rds")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
##  plot motif enrichment and TF expression together 
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

motif.tf.exp <- function(gene_list, fname){
  
  input_motif <- enriched.motifs$combname[match(gene_list, enriched.motifs$motif.name)]
  ##################################
  DefaultAssay(data_dc3000) <- 'chromvar'  
  
  tryCatch({
    q1 <- FeaturePlot(
      object = data_dc3000,
      features = input_motif,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = 1,
      split.by = "sample",
      order = TRUE,
      cols = c("lightgrey","blue"),
      reduction = "umap"
    )
    q1 <- q1 & scale_colour_viridis_c()& NoAxes()
    
    DefaultAssay(data_dc3000) <- 'SCT'
    q2 <- FeaturePlot(data_dc3000, features = gene_list, split.by = "sample", pt.size = 1, order=TRUE,
                      min.cutoff = 'q1',
                      max.cutoff = 'q99',
                      cols = c("lightgrey","darkmagenta"),
                      reduction = "umap"
    )
    q2 <- q2 & scale_colour_gradientn(
      colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
      na.value = "lightgray"
    ) & NoAxes()
    p1 <- ggarrange(q1, q2, nrow  =2)
    ggsave(filename = paste("motif_tf_plot/","_", fname,"_feture_dc3000.png", sep = ""), width = 16, height = 8, p1) 
  }, error=function(e){})
  
  ################################## ################################## ################################## ##################################
  DefaultAssay(data_avrrpt2) <- 'chromvar'  
  tryCatch({
    q1 <- FeaturePlot(
      object = data_avrrpt2,
      features = input_motif,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = 1,
      split.by = "sample",
      order = TRUE,
      cols = c("lightgrey","blue"),
      reduction = "umap"
    )
    q1 <- q1 & scale_colour_viridis_c()& NoAxes()
    
    DefaultAssay(data_avrrpt2) <- 'SCT'
    q2 <- FeaturePlot(data_avrrpt2, features = gene_list, split.by = "sample", pt.size = 1, order=TRUE,
                      min.cutoff = 'q1',
                      max.cutoff = 'q99', cols = c("lightgrey","darkmagenta"),
                      reduction = "umap"
    )
    q2 <- q2 & scale_colour_gradientn(
      colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
      na.value = "lightgray"
    ) & NoAxes()
    p1 <- ggarrange(q1, q2, nrow  =2)
    ggsave(filename = paste("motif_tf_plot/","_", fname,"_feture_avrrpt2.png", sep = ""), width = 16, height = 8, p1)
  }, error=function(e){})
  ################################## ################################## ##################################
  DefaultAssay(data_avrrpm1) <- 'chromvar'  
  tryCatch({
    q1 <- FeaturePlot(
      object = data_avrrpm1,
      features = input_motif,
      min.cutoff = 'q1',
      max.cutoff = 'q99',
      pt.size = 1,
      split.by = "sample",
      order = TRUE,
      cols = c("lightgrey","blue"),
      reduction = "umap"
    )
    q1 <- q1 & scale_colour_viridis_c()& NoAxes()
    
    DefaultAssay(data_avrrpm1) <- 'SCT'
    q2 <- FeaturePlot(data_avrrpm1, features = gene_list, split.by = "sample", pt.size = 1, order=TRUE,
                      min.cutoff = 'q1', cols = c("lightgrey","darkmagenta"),
                      max.cutoff = 'q90',
                      reduction = "umap"
    )
    q2 <- q2 & scale_colour_gradientn(
      colours =c("lightgray", "#BFD3E6" ,"#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D" ,"#810F7C" , "#4D004B"),
      na.value = "lightgray"
    ) & NoAxes()
    p1 <- ggarrange(q1, q2, nrow  =2)
    ggsave(filename = paste("motif_tf_plot/","_", fname,"_feture_avrrpm1.png", sep = ""), width = 16, height = 8, p1)
  }, error=function(e){})
  
  
  
  ################################## ################################## ################################## ##################################
  ##plots for legend
  DefaultAssay(data_avrrpt2) <- 'chromvar'  
  tryCatch({
    q1 <- FeaturePlot(
      object = data_avrrpt2,
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
    
    DefaultAssay(data_avrrpt2) <- 'SCT'
    q2 <- FeaturePlot(data_avrrpt2, features = gene_list, 
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
    p1 <- ggarrange(q1, q2, nrow  =2)
    ggsave(filename = paste("motif_tf_plot/","_", fname,"_feture_avrrpt2_legend.png", sep = ""), width = 5, height = 8, p1)
  }, error=function(e){})
 }

motif.tf.exp("WRKY46", "WRKY46")

##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## genes linked with motifs...
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
motif.link.genes <- function(target.motif){
  DefaultAssay(data) <- "peaks"
  motifID <- dm.table$motif[dm.table$motif.name==target.motif]
  
  motif.peak.table <- Motifs(data)@data
  target.peaks <- motif.peak.table[, motifID]
  target.peaks <- target.peaks[target.peaks==1]
  
  target.links <- links.comb[links.comb$peak %in% names(target.peaks) ,]
  target.links.filter <- target.links[target.links$score > 0.1 & target.links$width < 5000 ,]
  write.table(target.links.filter, file=paste0("motif_linked_genes/", target.motif,".txt"), row.names=T, col.names=NA, sep="\t", quote=F)
  target.links.filter$gene %>% table() %>% sort(decreasing = T)
}