rm(list = ls())

###############################
## INSTALL REQUIRED PACKAGES ##
###############################

if(!require('rstudioapi')) {
  install.packages('rstudioapi')
}
if(!require('Biobase')){
  install.packages('Biobase')
}
if(!require('flowCore')){
  install.packages("flowCore")
}
if(!require('FlowSOM')){
  install.packages("FlowSOM")
}
if(!require('ggplot2')){
  install.packages("ggplot2")
}
if(!require('dplyr')){
  install.packages('dplyr')
}
if(!require("flowVS")){
  install.packages(file.choose(), repos = NULL, type = "source")
}
if(!require('flowStats')){
  BiocManager::install('flowStats')
}
if(!require('flowSOMworkshop')){
  BiocManager::install('flowSOMworkshop')
}
if(!require('stringr')){
  BiocManager::install('stringr')
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
  library(readxl) 
  library(flowCore) 
  library(flowDensity)
  library(flowAI) 
  library(FlowSOM) 
  library(FlowSOMworkshop) 
  library(flowStats)
  library(CytoNorm)
  library(PeacoQC) 
  library(cytofCore)
  library(flowVS)
})

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# Define fcs_directory
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

#List CSV files
CSVfiles <- list.files(csvDirectory, pattern = ".csv$", full = FALSE)
fileName <- gsub(".csv", ".fcs", CSVfiles)

for(i in c(1:length(CSVfiles))){
  data <- read.csv(paste(csvDirectory, CSVfiles[i], sep = "/"))
  print(CSVfiles[i])
  print(fileName[i])
  cytofCore.write.FCS(as.matrix(data), 
                      filename = paste(csv2fcsDirectory, fileName[i], sep = "/"),
                      what = "numeric")
}


# read flowSet
# Create flowSet from FCSfiles
FCSfiles <- list.files(csv2fcsDirectory, pattern = ".fcs$", full = FALSE)
fs <- read.flowSet(files = FCSfiles, path = csv2fcsDirectory, truncate_max_range = FALSE)
colnames(fs)

colnames(fs)[7:15] <- c('c-kit', 'ki-67', 'CD44', 'CD45', 'CD4', 'CD8', 
                        'Thy1', 'CD3', 'CD25')

channels <- colnames(fs)
flowViz.par.set(theme =  trellis.par.get(), reset = TRUE)
as.vector(channels)
densityplot(~KLRG1+CD45RA+CD27+CD69+CD8+CD57+L_D+TIGIT+CD28+DNAM1+CD56+CCR7+PD1+CD25+Tim3+CD3, fs[c(21:24, 27:30)])

save(list = ls(), file =("timo_w8_step1.rds"))

