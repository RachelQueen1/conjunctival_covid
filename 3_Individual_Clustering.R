library(Seurat)
library(cowplot)
library(SCFunctionsV3)
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(DoubletFinder)
library(pheatmap)
library(ggplot2)

print("Individual clustering ...")
# Functions
howManyCells <- function(sObj){dim(sObj)[2]}

clusterCells <- function(sObj){
  sObj <- ScaleData(object = sObj, features = rownames(sObj))
  sObj <- RunPCA(object = sObj)
  sObj <- FindNeighbors(object = sObj)
  sObj <- FindClusters(object = sObj, resolution = 1.0)
  sObj <- RunUMAP(object = sObj, dims = 1:20)
  return(sObj)
}


# Cluster and save
clusterAndSave <- function(runNo){
  sObj_Filtered <- readRDS(paste0("rObjects/sObj_Filtered_", runNo, ".rds"))
  sampleInfo <- readRDS(paste0("rObjects/sampleInfo_", runNo, ".rds"))
  sObj_Filtered <- lapply(sObj_Filtered, clusterCells)
  saveRDS(sObj_Filtered, paste0("rObjects/sObj_Filtered_", runNo, ".rds"))
  return(sObj_Filtered)
}

# Find Markers
saveFindMarkers <- function(sObj_Filtered, runNo){
  for (i in 1:length(sObj_Filtered)){
    sampleInfo <- readRDS(paste0("rObjects/sampleInfo_", runNo, ".rds"))
    
    markers <- FindAllMarkers(sObj_Filtered[[i]], only.pos = TRUE, logfc.threshold = 0.5 )
    markers <- markers %>% arrange(cluster, desc(avg_log2FC))
    
    fName <- paste0("csvFiles/", sampleInfo$fName[i], runNo, ".csv")
    write.csv(markers, fName)
    fName <- gsub("csvFiles", "rObjects", fName) %>% gsub(".csv", ".rds", .)
    saveRDS(markers, fName)
    
  }
}

runNo <- "run1"
sObj_Filtered <- clusterAndSave(runNo)
saveFindMarkers(sObj_Filtered, runNo)

runNo <- "run2"
sObj_Filtered <- clusterAndSave(runNo)
saveFindMarkers(sObj_Filtered, runNo)
