rm(list = ls())

###############################
### LOAD REQUIRED PACKAGES ####
###############################

suppressPackageStartupMessages({
  library(rstudioapi)
  library(devtools)
  library(flowCore)
  library(cytofCore)
  library(FlowSOM)
  library(cluster)
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

load("FlowThymus_CD45+_w3_part1.rds")

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
antigen <- fcs_colname  
marker_class <- c("type", "state", rep("type", 7))
length(marker_class) == length(fcs_colname) #check 

# Create panel
panel <- cbind(fcs_colname, antigen, marker_class)
panel <- as.data.frame(panel)
all(panel$fcs_colname %in% colnames(fs)) #check

##Building metadata dataframe
#Define condition
condition <- FCSfiles
condition <- word(condition, 1, sep = "_")
condition[grepl("C", condition)] <- "Ctrl"
condition[grepl("Z", condition)] <- "Zd"

#Define patient_id
patient_id <- FCSfiles
patient_id <- word(patient_id, 2, sep = "_")

#Define sample_id
sample_id <- paste(patient_id, condition, sep = "_")

#Define file_name
file_name <- FCSfiles

#Create metadata dataframe
md <- as.data.frame(cbind(file_name, sample_id, condition, patient_id))


#Check if ids and md$file_name are the same
ids <- c(fsApply(fs, identifier))
ids%in%md$file_name

# create SingleCellExperiment object
sce <- CATALYST::prepData(fs, panel = panel, md = md, transform = FALSE, features = panel$fcs_colname)
sce@assays@data$exprs <- sce@assays@data$counts
assays(sce)

## QC
# Density plots
p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 4
p

# Check number of cells for each sample
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")

# Pseudo-bulk multidimensional scaling (pbMDS)
CATALYST::pbMDS(sce, color_by = "condition", label = "patient_id")

# Non Redundancy Score (NRS)
plotNRS(sce, features = type_markers(sce), color_by = "condition")

#Delete outlier 011
sce <- filterSCE(sce, patient_id != "011")

# Run FlowSOM and ConsensusClusterPlus
set.seed(1234)
sce <- cluster(sce,
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)
delta_area(sce)

#Check Heatmap
plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "meta16", 
                bars = TRUE, perc = TRUE, col_anno = TRUE, scale = "last", bin_anno = TRUE)

#Remove C16 - 0.17%
sce <- filterSCE(sce, cluster_id %in% c(1:15), k = "meta16")

# Run dimensionality reduction - UMAP
n_cells <- 3000
n_events <- min(n_cells(sce))
if(!(n_cells < n_events))
  n_cells <- n_events
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type")

# Plot UMAP
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
annotation_table <- readxl::read_excel("annotation.xlsx")
annotation_table

#Convert to factor with merged clusters in the desired order
annotation_table$new_cluster <- factor(annotation_table$new_cluster, 
                                       levels = c("DN2a", "DN1", "ETP", "DN2b", "DN3",  "Non lymphoid cells", "DP CD90 low",
                                                  "CD4+", "CD8+", "CD4+ immature", "ISP", "DN4", "DP CD90 intermediate", 
                                                  "DP CD90 high"))
#Apply manual annotation
sce <- mergeClusters(sce, k = "meta16", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)

#Heatmap with annotations
plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta16", m = "cluster_annotation", scale = "last", bin_anno = TRUE)

plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "cluster_annotation", scale = "last", bin_anno = TRUE)

#UMAP with annotations
dev.off()
plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("condition")
plot <- plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("condition")
plot + geom_density2d(binwidth = 0.006, colour = "black")

#Plot abundancies 
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", shape_by = NULL)

# Save workspace 
save(list = ls(), file = "FlowThymus_sce_CD45+_w3_part2.rds")
