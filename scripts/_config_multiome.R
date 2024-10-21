
# Check if Seurat is installed
if (!requireNamespace("Seurat", quietly = TRUE)) {
  # If not installed, check for remotes package and install if missing
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  # Install Seurat version 5 from GitHub
  remotes::install_github("satijalab/seurat", ref = "seurat5")
} else {
  # If Seurat is installed, check if the installed version is below 5
  installed_version <- packageVersion("seurat")
  target_version <- as.numeric(sub("\\..*$", "", as.character(installed_version))) # Extract major version part
  if (target_version < 5) {
    # If the major version of the installed Seurat is less than 5, install version 5
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github("satijalab/seurat", ref = "seurat5")
  } else {
    message("Seurat version 5 or greater is already installed.")
  }
}


# Define a function to install CRAN packages only if they are not already installed
install_if_missing <- function(packages) {
  missing_packages <- packages[!packages %in% installed.packages()[,"Package"]]
  message(paste0("missing packages:  "))
  if (length(missing_packages)) install.packages(missing_packages)
}

# Install CRAN packages if missing
install_if_missing(c("dplyr", "ggplot2", "cowplot", "ggrepel", "reshape", 
                     "tictoc", "ggpubr", "tidyverse", "pheatmap", "RColorBrewer", "scales", 
                     "viridis", "magrittr",  "vegan", "labdsv", "splitstackshape", 
                     "rlist", "patchwork", "egg", "harmony", "hdf5r", "ggseqlogo"))

# Install Bioconductor packages if BiocManager is not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install() # Update BiocManager to ensure it installs the latest packages

# Define a function to install Bioconductor packages only if they are not already installed
install_bioc_if_missing <- function(packages) {
  missing_packages <- packages[!packages %in% rownames(installed.packages())]
  if (length(missing_packages)) BiocManager::install(missing_packages)
}

# Install Bioconductor packages if missing
install_bioc_if_missing(c("Signac", "BRGenomics", "fastmatch", "GenomicFeatures", "BSgenome.Athaliana.TAIR.TAIR9",
                          "org.Athaliana.eg.db", "clusterProfiler", "JASPAR2020", "TFBSTools", 
                          "SingleCellExperiment", "qvalue","chromVAR", "DESeq2", "motifmatchr", "edgeR", "glmGamPoi"))

# Install packages from GitHub if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

install_github_if_missing <- function(repos) {
  for (repo in repos) {
    package_name <- unlist(strsplit(repo, "/"))[2]
    if (!package_name %in% installed.packages()[, "Package"]) {
      remotes::install_github(repo)
    }
  }
}

# Install GitHub packages if missing
install_github_if_missing(c('chris-mcginnis-ucsf/DoubletFinder'))

# Custom updated GO database
if(!requireNamespace("org.Athaliana.eg.db", quietly = TRUE)){
  # URL of the package
  package_url <- "http://neomorph.salk.edu/download/Nobori_etal_merfish/multiome/go_database/org.Athaliana.eg.db.tar.gz"
  
  # Destination file path (change the path to where you want to save the file)
  dest_file <- tempfile(fileext = ".tar.gz")
  
  # Download the package
  download.file(package_url, destfile = dest_file, mode = "wb")
  
  # Install the package from the downloaded file
  remotes::install_local(dest_file)

} else{
  print("Custom GO database is installed")
}

# Load or Install CRAN packages
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}
# Install or load devtools
install_if_missing("devtools")

# Function to check and install GitHub packages
install_github_if_missing <- function(repo) {
  package <- unlist(strsplit(repo, "/"))[2]
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    devtools::install_github(repo)
    library(package, character.only = TRUE)
  }
}

# Install 'qlcMatrix' from GitHub if it's not installed
install_github_if_missing("cysouw/qlcMatrix")


# Check if Signac is installed and meets the minimum version requirement
desired_version <- "1.12.9007"
if (!requireNamespace("Signac", quietly = TRUE) ||
    packageVersion("Signac") < as.package_version(desired_version)) {
  # If Signac is not installed or the version is less than the desired version,
  # check for the 'remotes' package and install it if it's missing
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  # Install the development version of Signac from GitHub
  remotes::install_github("stuart-lab/signac", ref = "develop")
} else {
  message(sprintf("Signac version %s or greater is already installed.", desired_version))
}


# List of all package names to load
packages_to_load <- c("Seurat", "Signac", "BRGenomics", "fastmatch", "dplyr", "ggplot2", 
                      "cowplot", "ggrepel", "reshape", "tictoc", "ggpubr", "tidyverse", 
                      "pheatmap", "RColorBrewer", "scales", "viridis", "magrittr", "edgeR", 
                      "vegan", "labdsv", "splitstackshape", "rlist", "patchwork", 
                      "egg", "SingleCellExperiment", "chromVAR", "DoubletFinder", "GenomicFeatures", 
                      "BSgenome.Athaliana.TAIR.TAIR9", "org.Athaliana.eg.db", "clusterProfiler", 
                      "JASPAR2020", "TFBSTools", "harmony")

# Loop to load each package
for (package in packages_to_load) {
  try(library(package, character.only = TRUE))
}

# color schemes
cols <- c("0" = "#721c28",
          "1" = "#004242",
          "2" = "#B2DF8A",
          "3" = "#3c1655",
          "4" = "#418557",
          "5" = "#18392b",
          "6" = "#3815ef",
          "7" = "#E78AC3",
          "8" = "#33A02C",
          "9" = "#D95F02",
          "10" = "#A65628",
          "11" = "#984EA3",
          "12" = "#FB8072",
          "13" = "#009088",
          "14" = "#eacd01",
          "15" = "#3a451c",
          "16" = "#008fff",
          "17" = "#08ff00",
          "18" = "#F4CAE4",
          "19" = "#ff007f",
          "20" = "#1B9E77",
          "21" = "#E41A1C",
          "22" = "#6A3D9A",
          "23" = "#8DD3C7",
          "24" = "#bd3434",
          "25" = "#B3DE69",
          "26" = "#666666",
          "27" = "#A6761D",
          "28" = "#FB9A99",
          "29" = "#CCCCCC"
          
)

col.sample <- c("Mock_rep1" = "#005b96",
                "Mock_rep2" = "#005b96",
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


cols.celltype <- c("epidermis" = "#721c28",
                   "mesophyll" = "#94ac78",
                   "undifferentiated" = "#009088",
                   "vasculature" = "#ff9f22")



##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==
## Common functions
##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==##==

genes_to_id <- function(genes){
  if (length(grep(genes, annotation$ID))==0 ){
    tryCatch({id <- annotation$ID[annotation$gene_name==genes] %>% na.omit
    id <- gsub("gene:", "", id)
    id[1] }, error=function(e){})
  }else{
    genes
  }
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

kplot <- function(input){
  wss <- function(k) {
    kmeans(input, k, nstart = 10 )$tot.withinss
  }
  # Compute and plot wss for k = 1 to k = 15
  k.values <- 1:15
  # extract wss for 2-15 clusters
  wss_values <- map_dbl(k.values, wss)
  plot(k.values, wss_values,
       type="b", pch = 19, frame = FALSE,
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
}
