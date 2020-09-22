rm(list = ls())

###############################
## INSTALL REQUIRED PACKAGES ##
###############################

if(!require('rstudioapi')) {
  install.packages('rstudioapi')
}
if(!require('flowCore')){
  install.packages("flowCore")
}
if(!require('ggplot2')){
  install.packages("ggplot2")
}
if(!require('dplyr')){
  install.packages('dplyr')
}
if(!require('cytofCore')){
  BiocManager::install('cytofCore')
}

###############################
### LOAD REQUIRED PACKAGES ####
###############################

suppressPackageStartupMessages({
  library(rstudioapi)
  library(ggplot2) 
  library(flowCore) 
  library(cytofCore)
  library(dplyr)
})

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# Define csv directory
csv <- "csv"
csvDirectory <- paste(PrimaryDirectory, csv, sep = "/")
dir.create(csvDirectory)

# Define workingDirectory
wdName <- "WorkingDirectory"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")
dir.create(workingDirectory)

#Define csv2fcsDirectory
csv2fcsDirectory <- "csv2fcsDirectory"
csv2fcsDirectory <- paste(PrimaryDirectory, csv2fcsDirectory, sep = "/")
dir.create(csv2fcsDirectory)

#List CSV files and add .fcs extension
CSVfiles <- list.files(csvDirectory, pattern = ".csv$", full = FALSE)
fileName <- gsub(".csv", ".fcs", CSVfiles)

#Obtain fcs files from csv files
for(i in c(1:length(CSVfiles))){
  data <- read.csv(paste(csvDirectory, CSVfiles[i], sep = "/"))
  print(CSVfiles[i])
  print(fileName[i])
  cytofCore.write.FCS(as.matrix(data), 
                      filename = paste(csv2fcsDirectory, fileName[i], sep = "/"),
                      what = "numeric")
}

# Create flowSet from FCSfiles
FCSfiles <- list.files(csv2fcsDirectory, pattern = ".fcs$", full = FALSE)
fs <- read.flowSet(files = FCSfiles, path = csv2fcsDirectory, truncate_max_range = FALSE)

#Look at the colnames(fs) and rename those of interest with the specific markers
colnames(fs)
colnames(fs)[7:15] <- c("c-kit", "ki-67", "CD44", "CD45", "CD4", "CD8", "Thy1", "CD3", 
                        "CD25")

save(list = ls(), file =("thymus_step1.rds"))

