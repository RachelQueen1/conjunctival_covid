library(Seurat)
library(harmony)

Combined <- readRDS("/data/rachel/Linda_Lako/Cornea/Adult_Cornea/rObjects/Combined.rds")
Combined$tissue <- "Adult Cornea"
Combined$cellType <- Combined$integrated_snn_res.0.6
Combined@active.assay <- "RNA"

seuratObj_C16_C10_C15_D33 <- readRDS("/data/rachel/Linda_Lako/Covid_SC/Repeat/rObjects/seuratObj_C16_C10_C15_D33.rds")
seuratObj_C16_C10_C15_D33$tissue <- "Organoid"
seuratObj_C16_C10_C15_D33$cellType <- seuratObj_C16_C10_C15_D33$annotation
seuratObj_C16_C10_C15_D33 <- SetIdent(seuratObj_C16_C10_C15_D33, value = "annotation")

# ## cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

compareObjects <- function(clustersUse, 
                           outName, 
                           clustersUse2 = levels(seuratObj_C16_C10_C15_D33@active.ident)){
  adult <- Combined[,Combined@active.ident %in% clustersUse]
  ALI <- seuratObj_C16_C10_C15_D33[, seuratObj_C16_C10_C15_D33@active.ident %in% clustersUse2]
  ## create batch corrected object
  seuratObj <- merge(adult, ALI ) %>%
    CellCycleScoring(s.features = s.genes,
                     g2m.features = g2m.genes,
                     set.ident = FALSE) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE,
              vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")) %>%
    RunPCA(pc.genes = seuratObj_ds@var.genes, npcs = 20, verbose = FALSE) %>%
    RunHarmony("orig.ident", plot_convergence = FALSE) 
  
  seuratObj <-seuratObj %>% RunUMAP(reduction = "harmony", dims = 1:10) 
  saveRDS(seuratObj, paste0("rObjects/cornea_comparison_", outName, ".rds"))
  
  return(seuratObj)
}
getwd()

seurat_all <- compareObjects(0:20, "all")
seurat_0_6 <- compareObjects(c(0,6), "0_6")
seurat_0_6_9 <- compareObjects(c(0,6,9), "0_6_9")

seurat_basal <- compareObjects(0, 
                               "basal_0", 
                               "basal conjunctival epithelium" )


seurat_basal_0_6 <- compareObjects(c(0,6), 
                               "basal_0_6", 
                               "basal conjunctival epithelium" )

seurat_superficial <- compareObjects(6, 
                               "superficial_6", 
                               "superficial conjunctival epithelium" )

seurat_superficial_0_6 <- compareObjects(c(0,6), 
                                     "superficial_0_6", 
                                     "superficial conjunctival epithelium" )

DimPlot(seurat_all,  group.by = "cellType", split.by = "tissue", label = TRUE)
DimPlot(seurat_0_6,  group.by = "cellType", split.by = "tissue", label = TRUE)
DimPlot(seurat_0_6_9,  group.by = "cellType", split.by = "tissue", label = TRUE)

DimPlot(seurat_basal,  group.by = "cellType", label = TRUE)

DimPlot(seurat_basal_0_6,  group.by = "cellType", label = TRUE)
DimPlot(seurat_basal_0_6,  group.by = "cellType", label = TRUE)
DimPlot(seurat_superficial_0_6,  group.by = "cellType", label = TRUE)
DimPlot(seurat_superficial,  group.by = "cellType", label = TRUE)
