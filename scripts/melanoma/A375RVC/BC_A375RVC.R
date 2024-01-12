library("beyondcell")
library("Seurat")
library("clustree")
set.seed(1)
# Read single-cell experiment.
A375RVC_data <- readRDS(file="/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/A375RVC.rds")
bc_recomputed <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/bc_A375RVC.rds")

SSc <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"))
bc <- bcScore(A375RVC_data, SSc, expr.thres = 0.1) 
bc_filtered <- bcSubset(bc, nan.cells = 0.95)
bc_filtered@normalized[is.na(bc_filtered@normalized)] <- 0                                                                                                                                                                       
bc_recomputed <- bcRecompute(bc_filtered, slot = "normalized") 

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
saveRDS(bc_recomputed, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/bc_A375RVC.rds")
saveRDS(A375RVC_data, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/A375RVC.rds")

clustree(bc_recomputed@meta.data, prefix ="bc_clusters_res.")

bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)
bc_recomputed <- bcRegressOut(bc_recomputed, vars.to.regress = "nFeature_RNA", k.neighbors = 40)

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)
clustree(bc_recomputed@meta.data, prefix ="bc_clusters_res.")


bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.8", label.size = 5, pt.size = 1.5)
A375RVC_data <- AddMetaData(object = A375RVC_data, metadata = bc_recomputed@meta.data)

# Expression UMAP.
bcClusters(bc_recomputed, UMAP = "Seurat", idents = "RNA_snn_res.0.4", pt.size = 1.5, label = TRUE)


DimPlot(
  A375RVC_data,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  group.by = "bc_clusters_res.0.8",
  #shape.by = "RNA_snn_res.0.4",
  raster = FALSE
)



DimPlot(
  A375RVC_data,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  group.by = "bc_clusters_res.0.6",
  shape.by = "subclone",
  raster = FALSE
)




# Obtain general statistics.
bc_recomputed <- bcRanks(bc_recomputed)
# Obtain condition-based statistics.
bc_recomputed <- bcRanks(bc_recomputed, idents = "bc_clusters_res.0.4")
# Obtain unextended therapeutic cluster-based statistics.
#bc <- bcRanks(bc, idents = "bc_clusters_res.0.4", extended = FALSE)

#Búsqueda de firmas de fármacos
dabrafenib_IDs <- FindDrugs(bc_recomputed, "dabrafenib")$IDs
#vemurafenib_IDs <- c("sig-20547", "sig-3758", "sig-6291")
trametinib_IDs <- FindDrugs(bc_recomputed, "trametinib")$IDs
cobimetinib_IDs <- FindDrugs(bc_recomputed, "cobimetinib")$IDs

# Patchwork object with all bortezomib plots.
dabrafenib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = dabrafenib_IDs), pt.size = 1.5)


trametinib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                           signatures = list(values = trametinib_IDs), pt.size = 1.5)

cobimetinib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = cobimetinib_IDs), pt.size = 1.5)
# Plot all three plots in one image.
wrap_plots(dabrafenib, ncol = 2, widths = 7)
wrap_plots(trametinib, ncol = 3, widths = 7)
wrap_plots(cobimetinib, ncol = 1)







bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 0 , top = 5)
ggsave(filename = "bc_recomputed_drugs0_A375RVC.jpg", path = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/plots_reg/bc_recomputed/", width = 16, height = 9)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 1 , top = 5)
ggsave(filename = "bc_recomputed_drugs1_A375RVC.jpg", path = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/plots_reg/bc_recomputed/", width = 16, height = 9)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 2 , top = 5)
ggsave(filename = "bc_recomputed_drugs2_A375RVC.jpg", path = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/plots_reg/bc_recomputed/", width = 16, height = 9)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 3 , top = 5)
ggsave(filename = "bc_recomputed_drugs3_A375RVC.jpg", path = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/plots_reg/bc_recomputed/", width = 16, height = 9)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 4 , top = 5)
ggsave(filename = "bc_recomputed_drugs4_A375RVC.jpg", path = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/plots_reg/bc_recomputed/", width = 16, height = 9)





#Sacar clusters terapeuticos con ssc, pensar si los clusters son adecuados en cuanto a número y composición, crear tabla con célula y a qué cluster pertenece




A375RVC_data <- AddMetaData(A375RVC_data, metadata =bc_recomputed@meta.data)

###########VEMURAFENIB#############
SSc <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"), )
vemurafenib_gs <- GenerateGenesets(x = "/home/lpgonzalezm/Downloads/vemurafenib_signatures.gmt", )

bc <- bcScore(A375RVC_data, SSc, expr.thres = 0.05) 
bc_filtered <- bcSubset(bc, nan.cells = 0.95)
bc_filtered@normalized[is.na(bc_filtered@normalized)] <- 0                                                                                                                                                                       
bc_recomputed <- bcRecompute(bc_filtered, slot = "normalized") 

bc_vemurafenib <- bcScore(A375RVC_data, vemurafenib_gs, expr.thres = 0.05) 
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

saveRDS(bc_recomputed_merged, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/bc_merged_A375RVC.rds")











