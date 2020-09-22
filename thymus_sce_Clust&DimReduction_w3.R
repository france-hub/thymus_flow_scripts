rm(list = ls())

###############################
## INSTALL REQUIRED PACKAGES ##
###############################

if(!require('rstudioapi')) {
  install.packages('rstudioapi')
}
if(!require('devtools')){
  install.packages("devtools")
}
if(!require('flowCore')){
  install.packages("flowCore")
}
if(!require('cytofCore')){
  devtools::install_github("nolanlab/cytofCore")
}
if(!require('cytofkit2')){
  #do not update hexbin
  devtools::install_github("JinmiaoChenLab/cytofkit2")
}
if(!require('FlowSOM')){
  install.packages("FlowSOM")
}
if(!require('cluster')){
  install.packages("cluster")
}
if(!require('Rtsne')){
  install.packages("Rtsne")
}
if(!require('ggplot2')){
  install.packages("ggplot2")
}
if(!require('dplyr')){
  install.packages("dplyr")
}
if(!require('ggthemes')){
  install.packages("ggthemes")
}
if(!require('RColorBrewer')){
  install.packages('RColorBrewer')
}
if(!require("uwot")){
  install.packages("uwot")
}
if(!require("CATALYST")){
  BiocManager::install("CATALYST")
}
if(!require("diffcyt")){
  BiocManager::install("diffcyt")
}
if(!require("stringr")){
  BiocManager::install("stringr")
}
if(!require("scran")){
  BiocManager::install("scran")
}
if(!require("scater")){
  BiocManager::install("scater")
}
if(!require("ggcyto")){
  BiocManager::install("ggcyto")
}
if(!require("SingleCellExperiment")){
  BiocManager::install("SingleCellExperiment")
}
if(!require("flowWorkspace")){
  BiocManager::install("flowWorkspace")
}
if(!require("reshape2")){
  BiocManager::install("reshape2")
}
if(!require("ggrepel")){
  BiocManager::install("ggrepel")
}
if(!require("knn.covertree")){
  devtools::install_github('flying-sheep/knn.covertree')
}
if(!require("theislab/destiny")){
  devtools::install_github('theislab/destiny')
}

###############################
### LOAD REQUIRED PACKAGES ####
###############################

suppressPackageStartupMessages({
  library(rstudioapi)
  library(flowCore)
  library(cytofCore)
  library(cytofkit2)
  library(FlowSOM)
  library(cluster)
  library(Rtsne)
  library(ggplot2)
  library(dplyr)
  library(ggthemes)
  library(RColorBrewer)
  library(uwot)
  library(CATALYST)
  library(diffcyt)
  library(stringr)
  library(scran)
  library(scater)
  library(ggcyto)
  library(SingleCellExperiment)
  library(flowWorkspace)
  library(reshape2)
  library(ggrepel)
  library(slingshot)
  library(knn.covertree)
  library(destiny)
  library(readxl)
})

############################
######## LOAD DATA #########
############################

### Load the FCS files 
# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

#load("thymus_step1.rds")

############################
######BUILD SCE OBJECT######
############################

# Keyword ($CYT = "FACS")
names(keyword(fs[[1]]))[15] <- "$CYT"
ds <- keyword(fs[[1]])
l <- list(cyt = "\\$CYT$")
keep <- lapply(l, grep, names(ds))
ds[[keep$cyt]] <- "FACS"
keyword(fs[[1]])[[keep$cyt]] <- "FACS"

##Building panel dataframe
# Define channels of interest and marker_class
fcs_colname <- colnames(fs)[7:15]

marker_class <- c("type", "state", rep("type", 7))
antigen <- fcs_colname
length(marker_class) == length(fcs_colname)

#Panel
panel <- cbind(fcs_colname, antigen, marker_class)
panel <- as.data.frame(panel)
all(panel$fcs_colname %in% colnames(fs))

##Building metadata dataframe
condition <- FCSfiles
condition <- word(condition, 1, sep = "_")
condition[grepl("C", condition)] <- "Ctrl"
condition[grepl("Z", condition)] <- "Zd"


patient_id <- FCSfiles
patient_id <- word(patient_id, 2, sep = "_")

sample_id <- paste(patient_id, condition, sep = "_")

file_name <- FCSfiles

md <- cbind(file_name, sample_id, condition, patient_id)
md <- data.frame(md)

#Check if ids and md$file_name are the same
ids <- c(fsApply(fs, identifier))
ids%in%md$file_name

##SCE object
sce <- CATALYST::prepData(fs, panel = panel, md = md, transform = FALSE, features = panel$fcs_colname)
sce@assays@data$exprs <- sce@assays@data$counts

## QC
# Density plots
p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 4
p

# Counts
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")

# Multidimensional scaling (MDS)
CATALYST::pbMDS(sce, color_by = "condition", label = "patient_id")

# Non Redundancy Score (NRS)
plotNRS(sce, features = type_markers(sce), color_by = "condition")

#Delete outlier 011
sce <- filterSCE(sce, patient_id != "011")

## Clustering
#FlowSOM
set.seed(1234)
sce <- cluster(sce,
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)
delta_area(sce)
plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "meta16", 
                bars = TRUE, perc = TRUE, col_anno = TRUE, scale = "last")

#Filter very small clusters (< 0.4%)
sce <- filterSCE(sce, cluster_id %in% c(1:15), k = "meta16")

#Run UMAP algorithm with a downsampling of 3000
set.seed(1234)
n_cells <- 3000
exaggeration_factor <- 12.0
eta <- n_cells/exaggeration_factor
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type")

#Plot UMAP
plotDR(sce, "UMAP", color_by = "condition")
plotDR(sce, "UMAP", color_by = type_markers(sce)) #markers distribution
plotDR(sce, "UMAP", color_by = "meta16") + facet_wrap("condition") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))) #different conditions
plotDR(sce, "UMAP", color_by = "meta16", facet_by = "sample_id") #by sample_id

#Plot UMAP + contour
plot <- plotDR(sce, "UMAP", color_by = "meta16") + facet_wrap("condition")
plot + geom_density2d(binwidth = 0.006, colour = "black")

#Abundances
#Differential abundance
plotAbundances(sce, k = "meta6", by = "sample_id", group_by = "condition")
plotAbundances(sce, k = "meta16", by = "cluster_id", group_by = "condition")

##Add annotations
#Read annotation excel file
annotation_table <- readxl::read_excel(file.choose())
annotation_table

#convert to factor with merged clusters in desired order
annotation_table$new_cluster <- factor(annotation_table$new_cluster, 
                                       levels = c("DN2a", "DN1", "ETP", "DN2b", "DN3",  "Non lymphoid cells", "DP CD90 low",
                                                  "CD4+", "CD8+", "CD4+ immature", "ISP", "DN4", "DP CD90 intermediate", 
                                                  "DP CD90 high"))
# apply manual annotation
sce <- mergeClusters(sce, k = "meta16", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)

#Heatmap with annotations
plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta16", m = "cluster_annotation", scale = "last")

plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "cluster_annotation", scale = "last")

#UMAP with annotations
dev.off()
levels(cluster_id) == levels("cluster")
plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("condition")
plot <- plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("condition")
plot + geom_density2d(binwidth = 0.006, colour = "black")

#Abundancies with the annotations
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", shape_by = NULL)

# Save workspace 
save(list = ls(), file = "CD45+_timo_2nd.rds")
