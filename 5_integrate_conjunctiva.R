library(harmony)
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

## read data
sObj_Filtered_run1 <- readRDS("rObjects/sObj_Filtered_run1.rds")
sObj_Filtered_run2 <- readRDS("rObjects/sObj_Filtered_run2.rds")

sampleInfo_run1 <- readRDS("rObjects/sampleInfo_run1.rds")
sampleInfo_run2 <- readRDS("rObjects/sampleInfo_run2.rds")


## select samples
samples_1 <- sampleInfo_run1$fName[sampleInfo_run1$Tissue == "Conjunctiva"]
samples_2 <- sampleInfo_run2$fName[sampleInfo_run2$Cell.line == "C10"]

sObj_Use <- c(sObj_Filtered_run1[samples_1],
              sObj_Filtered_run2[samples_2])


## cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

## create batch corrected object
seuratObj <- merge(sObj_Use[[1]], sObj_Use[-1]) %>% 
  CellCycleScoring(s.features = s.genes, 
                   g2m.features = g2m.genes, 
                   set.ident = FALSE) %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE, 
            vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")) %>% 
  RunPCA(pc.genes = seuratObj_ds@var.genes, npcs = 20, verbose = FALSE) %>%
  RunHarmony("orig.ident", plot_convergence = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:10) %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = c(0.2,1,2.2))


#find_markers.R
markers <- FindAllMarkers(seuratObj, logfc.threshold = 0.7)## DE

markers <- markers %>% dplyr::arrange(cluster, desc(avg_log2FC))

fName <- paste0("csvFiles/", "C16_C10_C15_D33", ".csv")
write.csv(markers, fName)
fName <- gsub("csvFiles", "rObjects", fName) %>% gsub(".csv", ".rds", .)
saveRDS(markers, fName)

saveRDS(seuratObj, "rObjects/seuratObj_C16_C10_C15_D33.rds")





