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
if(!require("Rphenograph")){
  BiocManager::install("Rphenograph")
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
if(!require("slingshot")){
  BiocManager::install("slingshot")
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
library(Rphenograph)
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

### Load the (transformed, normalized) FCS files 

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

#load("timo_CD45w8.rds")

fs
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
condition <- word(condition, 4, sep = "_")


patient_id <- FCSfiles
patient_id <- word(patient_id, 5, sep = "_")

sample_id <- paste(patient_id, condition, sep = "_")

file_name <- FCSfiles

md <- cbind(file_name, sample_id, condition, patient_id)
md <- data.frame(md)

ids <- c(fsApply(fs, identifier))
ids%in%md$file_name

##SCE object
sce <- CATALYST::prepData(fs, panel = panel, md = md, transform = FALSE, features = panel$fcs_colname)
sce@assays@data$exprs <- sce@assays@data$counts

## QC
# Density
p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 4
p

# Counts
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")


# Multidimensional scaling (MDS)
CATALYST::pbMDS(sce, color_by = "condition", label = "sample_id")

# Non Redundancy Score (NRS)
plotNRS(sce, features = type_markers(sce), color_by = "condition")

#delete outlier
sce <- filterSCE(sce, sample_id != "005_Zd")

## Clustering
#FlowSOM
set.seed(1234)
sce <- cluster(sce,
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)

delta_area(sce)
plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "meta16", 
                bars = TRUE, perc = TRUE, col_anno = TRUE, scale = "last", bin_anno = TRUE)
plotExprHeatmap(sce, features = "type",
                by = "sample_id", k = "meta20", 
                bars = TRUE, perc = TRUE, col_anno = TRUE, scale = "first", row_anno = "condition")

#delete clusters under 1%
sce <- filterSCE(sce, cluster_id %in% c(1:15), k = "meta16")

set.seed(1234)
n_cells <- 2200
exaggeration_factor <- 12.0
eta <- n_cells/exaggeration_factor
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type")

# UMAP
plotDR(sce, "UMAP", color_by = "condition")
plotDR(sce, "UMAP", color_by = type_markers(sce)) #markers distribution
plotDR(sce, "UMAP", color_by = "meta15") + facet_wrap("condition") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))) #different conditions
plotDR(sce, "UMAP", color_by = "meta15", facet_by = "sample_id") #by sample_id

save(list = ls(), file =("timo_w8_CD45_def.rds"))


#Abundances
# Differential abundance
plotAbundances(sce, k = "meta6", by = "sample_id", group_by = "condition")
plotAbundances(sce, k = "meta15", by = "cluster_id", group_by = "condition")


##Add annotations
#Read annotation file
annotation_table <- readxl::read_excel(file.choose())
annotation_table

levels <- factor(c(annotation_table$new_cluster))

# convert to factor with merged clusters in desired order
annotation_table$new_cluster <- factor(annotation_table$new_cluster, 
                                       levels = levels)

# apply manual annotation
sce <- mergeClusters(sce, k = "meta15", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)

#Heatmap with annotations
plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta15", m = "cluster_annotation", scale = "last")

plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "cluster_annotation", scale = "last", bin_anno = TRUE)

#UMAP with annotations
dev.off()
plot <- plotDR(sce, "UMAP", color_by = "meta15") + facet_wrap("condition")
plot <- plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("condition")
plot + geom_density2d(binwidth = 0.006, colour = "black")

#tSNE with annotations
plotDR(sce, "TSNE", color_by = "cluster_annotation") + facet_wrap("condition")

#Abundancies
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", shape_by = NULL)

save(list = ls(), file = "timo_CD45w8.rds")
load("timo_CD45w8.rds")

# Statistical Analysis
ei <- sce@metadata$experiment_info

design <- createDesignMatrix(
  ei, cols_design = c("condition", "patient_id")
)

design <- design[,-c(20)] # not a full rank matrix, drop column 19
#which(colnames(design) == "patient_id151")
#DG NR vs DG CR
contrast <- createContrast(c(0, 1, 0, 0, rep(0, 17)))
out_DA <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE)

topTable(out_DA, format_vals = TRUE)
da_1 <- rowData(out_DA$res)
sub_1 <- filterSCE(sce, condition == c("DG_CR", "DG_NR"))
CR_ind <- grep("CR", levels(sub_1$sample_id))
NR_ind <- grep("NR", levels(sub_1$sample_id))
levels(sub_1$sample_id) <- levels(sub_1$sample_id)[c(CR_ind, NR_ind)]
plotDiffHeatmap(sub_1, da_1, top_n = 12, all = TRUE, fdr = 0.05, sort_by = "padj", fun = "mean")


#post CR vs post NR
contrast <- createContrast(c(0, 0, -1, 1, rep(0, 17)))

out_DA <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "meta7", verbose = TRUE, subsampling = TRUE, transform = FALSE)
topTable(out_DA, format_vals = TRUE)
da <- rowData(out_DA$res)
sub_2 <- filterSCE(sce, condition == c("NR_BM_DG", "HC"))
NR_ind <- grep("NR", levels(sub_2$sample_id))
HC_ind <- grep("HC", levels(sub_2$sample_id))
levels(sub_2$sample_id) <- levels(sub_2$sample_id)[c(NR_ind, HC_ind)]
plotDiffHeatmap(sub_2, da, top_n = 12, all = TRUE, fdr = 0.05, sort_by = "padj", fun = "mean")

#DG CR post_CR
contrast <- createContrast(c(0, 0, 1, 0, rep(0, 17)))

out_DA <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE)
topTable(out_DA, format_vals = TRUE)
da_2 <- rowData(out_DA$res)
sub_3 <- filterSCE(sce, condition == c("DG_CR", "post_CR"))
DG_ind <- grep("DG", levels(sub_3$sample_id))
post_ind <- grep("post", levels(sub_3$sample_id))
levels(sub_3$sample_id) <- levels(sub_3$sample_id)[c(DG_ind, post_ind)]
plotDiffHeatmap(sub_3, da_2, top_n = 12, all = TRUE, fdr = 0.05, sort_by = "padj", fun = "mean")

#NR_DG vs NR_post
contrast <- createContrast(c(0, -1, 0, 1, rep(0, 17)))
out_DA <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE)
topTable(out_DA, format_vals = TRUE)
da_3 <- rowData(out_DA$res)
sub_4 <- filterSCE(sce, condition %in% c("DG_NR", "post_NR"))
DG_ind <- grep("DG", levels(sub_4$sample_id))
post_ind <- grep("post", levels(sub_4$sample_id))
levels(sub_4$sample_id) <- levels(sub_4$sample_id)[c(DG_ind, post_ind)]
plotDiffHeatmap(sub_4, da_3, top_n = 12, all = TRUE, fdr = 0.05, sort_by = "padj", fun = "mean")

library(ggpubr)

#DG CR vs DG NR
stat.test_1 <- as_tibble(da_1)
p.adj.signif <- c("ns", "1.9e-02", rep("ns",3),"1.9e-02","ns")
y.position <- c(58, 35, 60, 20, 22, 52, 10)
group1 <- (rep("DG_CR",7))
group2 <- (rep("DG_NR", 7))
stat.test <- cbind(stat.test, group1, group2, p.adj.signif, y.position)
bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", shape_by = "patient_id")
bxp <- bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, size = 2.5) 
bxp


#pre_CR vs post_CR
stat.test_2 <- as_tibble(da_2)
p.adj.signif <- c(rep("ns", 3), "6.8e-02", rep("ns", 3))
y.position <- c(50, 30,45,24,22,40,13)
group1 <- (rep("DG_CR",7))
group2 <- (rep("post_CR", 7))
stat.test <- cbind(stat.test_2, group1, group2, p.adj.signif, y.position)
bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", shape_by = "patient_id")
bxp <- bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, size = 2.5) 
bxp

#pre_NR vs post_NR
stat.test <- as_tibble(da_3)
p.adj.signif <- c(rep("ns",2), "8.2e-02", "5.53-03", '4.9e-02', "ns", "ns")
y.position <- c(40, 32, 38,35, 25,40,14)
group1 <- (rep("DG_NR",7))
group2 <- (rep("post_NR", 7))
stat.test <- cbind(stat.test, group1, group2, p.adj.signif, y.position)
bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", shape_by = "patient_id")
bxp <- bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, size = 2.5) 
bxp
  
# save workspace 
