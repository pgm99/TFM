library("beyondcell")
library("Seurat")
library("clustree")
library("patchwork")

set.seed(1)
# Read single-cell experiment.
A375S_data <- readRDS(file="/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/A375S.rds")
bc_recomputed <- readRDS(file="/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/bc_A375S.rds")

SSc <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"), )

bc <- bcScore(A375S_data, SSc, expr.thres = 0.1) 
bc_filtered <- bcSubset(bc, nan.cells = 0.95)
bc_filtered@normalized[is.na(bc_filtered@normalized)] <- 0                                                                                                                                                                       
bc_recomputed <- bcRecompute(bc_filtered, slot = "normalized") 

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
saveRDS(bc_recomputed, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/bc_A375S.rds")
saveRDS(A375S_data, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/A375S.rds")
bc_recomputed <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/A375S_bc.rds")

clustree(bc_recomputed@meta.data, prefix ="bc_clusters_res.")
clustree(x = bc@meta.data, prefix = "bc_clusters_res.")

bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)
bc_recomputed <- bcRegressOut(bc_recomputed, vars.to.regress = "nFeature_RNA", k.neighbors = 40)

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)


bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.6", label.size = 5, pt.size = 1.5)
A375S_data <- AddMetaData(object = A375S_data, metadata = bc_recomputed@meta.data)

# Expression UMAP.
bcClusters(bc_recomputed, UMAP = "Seurat", idents = "RNA_snn_res.0.6", pt.size = 1.5, label = TRUE)


DimPlot(
  A375S_data,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  group.by = "bc_clusters_res.0.6",
  #shape.by = "trametinib",
  raster = FALSE
)



DimPlot(
  A375S_data,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  group.by = "bc_clusters_res.0.4",
  shape.by = "subclone",
  raster = FALSE
)




# Obtain general statistics.
bc_recomputed <- bcRanks(bc_recomputed)
# Obtain condition-based statistics.
bc_recomputed <- bcRanks(bc_recomputed, idents = "bc_clusters_res.0.6")
# Obtain unextended therapeutic cluster-based statistics.
#bc <- bcRanks(bc, idents = "bc_clusters_res.0.4", extended = FALSE)

#Búsqueda de firmas de fármacos
vemurafenib_IDs <- FindDrugs(bc_recomputed, "vemurafenib")$IDs
vemurafenib_IDs <- c("sig-20547", "sig-3758", "sig-6291")

dabrafenib_IDs <- FindDrugs(bc_recomputed, "dabrafenib")$IDs
trametinib_IDs <- FindDrugs(bc_recomputed, "trametinib")$IDs
cobimetinib_IDs <- FindDrugs(bc_recomputed, "cobimetinib")$IDs
PLX <- FindDrugs(bc_recomputed, "PLX_4032")$IDs


# Patchwork object with all bortezomib plots.
dabrafenib <- bcSignatures(bc_recomputed, UMAP = "Seurat", 
                           signatures = list(values = "sig-21060"), pt.size = 1.5)


trametinib <- bcSignatures(bc_recomputed, UMAP = "Seurat", 
                           signatures = list(values = "sig-20908"), pt.size = 1.5)

cobimetinib <- bcSignatures(bc_recomputed, UMAP = "Seurat", 
                           signatures = list(values = cobimetinib_IDs), pt.size = 1.5)

prueba <- bcSignatures(bc_recomputed, UMAP = "Seurat", 
                            signatures = list(values = "PLX-4720"), pt.size = 1.5)


# Plot all three plots in one image.
wrap_plots(dabrafenib, ncol = 1)
wrap_plots(trametinib, ncol = 1)
wrap_plots(cobimetinib, ncol = 1)












bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 0 , top = 5)
ggsave(filename = "bc_recomputed_drugs0_A375s.jpg", path = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/plots_reg/bc_recomputed/", width = 16, height = 9)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 1 , top = 5)
ggsave(filename = "bc_recomputed_drugs1_A375s.jpg", path = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/plots_reg/bc_recomputed/", width = 16, height = 9)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 2 , top = 5)
ggsave(filename = "bc_recomputed_drugs2_A375s.jpg", path = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/plots_reg/bc_recomputed/", width = 16, height = 9)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 3 , top = 5)
ggsave(filename = "bc_recomputed_drugs3_A375s.jpg", path = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/plots_reg/bc_recomputed/", width = 16, height = 9)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 4 , top = 5)
ggsave(filename = "bc_recomputed_drugs4_A375s.jpg", path = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/plots_reg/bc_recomputed/", width = 16, height = 9)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 5 , top = 5)
ggsave(filename = "bc_recomputed_drugs5_A375s.jpg", path = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/plots_reg/bc_recomputed/", width = 16, height = 9)





#Sacar clusters terapeuticos con ssc, pensar si los clusters son adecuados en cuanto a número y composición, crear tabla con célula y a qué cluster pertenece




A375S_data <- AddMetaData(A375S_data, metadata =bc_recomputed@meta.data)


###########VEMURAFENIB#############
SSc <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"), )
vemurafenib_gs <- GenerateGenesets(x = "/home/lpgonzalezm/Downloads/vemurafenib_signatures.gmt", )

bc <- bcScore(A375S_data, SSc, expr.thres = 0.05) 
bc_filtered <- bcSubset(bc, nan.cells = 0.95)
bc_filtered@normalized[is.na(bc_filtered@normalized)] <- 0                                                                                                                                                                       
bc_recomputed <- bcRecompute(bc_filtered, slot = "normalized") 

bc_vemurafenib <- bcScore(A375S_data, vemurafenib_gs, expr.thres = 0.1) 
bc_vemurafenib_filtered <- bcSubset(bc_vemurafenib, nan.cells = 0.95)
bc_vemurafenib_filtered@normalized[is.na(bc_vemurafenib_filtered@normalized)] <- 0                                                                                                                                                                       
bc_vemurafenib_recomputed <- bcRecompute(bc_vemurafenib_filtered, slot = "normalized") 

bc_recomputed_merged <- bcMerge(bc1 = bc_vemurafenib_recomputed, bc2 = bc_recomputed)

bc_recomputed_merged <- bcRegressOut(bc_recomputed_merged, vars.to.regress = "nFeature_RNA")
bc_recomputed_merged <- bcUMAP(bc_recomputed_merged, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
clustree(bc_recomputed_merged@meta.data, prefix ="bc_clusters_res.")

bcClusters(bc_recomputed_merged, UMAP = "beyondcell", idents = "bc_clusters_res.0.8", label.size = 5, pt.size = 1.5)

vemurafenib_skin <- bcSignatures(bc_recomputed_merged, UMAP = "beyondcell", 
                                 signatures = list(values = "vemurafenib_skin"), pt.size = 1.5)
vemurafenib_pancancer <- bcSignatures(bc_recomputed_merged, UMAP = "beyondcell", 
                                      signatures = list(values = "vemurafenib_pancancer"), pt.size = 1.5)
wrap_plots(vemurafenib_skin, ncol = 1)
wrap_plots(vemurafenib_pancancer, ncol = 1)

saveRDS(bc_recomputed_merged, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/bc_merged_A375S.rds")



###############FUNCTIONAL##################

SSc_functional <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"), include.pathways = TRUE)
bc_functional <- bcScore(A375S_data, SSc_functional, expr.thres = 0.1) 
bc_filtered_functional <- bcSubset(bc_functional, nan.cells = 0.95)
bc_filtered_functional@normalized[is.na(bc_filtered_functional@normalized)] <- 0                                                                                                                                                                       
bc_recomputed_functional <- bcRecompute(bc_filtered_functional, slot = "normalized") 


bcClusters(bc_recomputed_functional, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bc_recomputed_functional <- bcRegressOut(bc_recomputed_functional, vars.to.regress = "nFeature_RNA", k.neighbors = 40)
bc_recomputed_functional <- bcUMAP(bc_recomputed_functional, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
clustree(bc_recomputed_functional@meta.data, prefix ="bc_clusters_res.")
bcClusters(bc_recomputed_functional, UMAP = "beyondcell", idents = "bc_clusters_res.0.8", label.size = 5, pt.size = 1.5)



