library("beyondcell")
library("Seurat")
library("clustree")
set.seed(1)
# Read single-cell experiment.
A375RV_data <- readRDS(file="/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RV/A375RV.rds")
bc_recomputed <- readRDS(file="/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RV/bc_A375RV.rds")

SSc <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"), )
vemurafenib_gs <- GenerateGenesets(x = "/home/lpgonzalezm/Downloads/vemurafenib_signatures.gmt", )
# ListFilters(entry = "MoAs")

bc <- bcScore(A375RV_data, SSc, expr.thres = 0.05) 
bc_filtered <- bcSubset(bc, nan.cells = 0.95)
bc_filtered@normalized[is.na(bc_filtered@normalized)] <- 0                                                                                                                                                                       
bc_recomputed <- bcRecompute(bc_filtered, slot = "normalized") 

bc_vemurafenib <- bcScore(A375RV_data, vemurafenib_gs, expr.thres = 0.05) 
bc_vemurafenib_filtered <- bcSubset(bc_vemurafenib, nan.cells = 0.95)
bc_vemurafenib_filtered@normalized[is.na(bc_vemurafenib_filtered@normalized)] <- 0                                                                                                                                                                       
bc_vemurafenib_recomputed <- bcRecompute(bc_vemurafenib_filtered, slot = "normalized") 

bc_recomputed_merged <- bcMerge(bc1 = bc_vemurafenib_recomputed, bc2 = bc_recomputed)

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
saveRDS(bc_recomputed, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RV/bc_A375RV.rds")
saveRDS(A375RV_data, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RV/A375RV.rds")
bc_recomputed <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RV/A375RV_bc.rds")

clustree(bc_recomputed@meta.data, prefix ="bc_clusters_res.")
clustree(x = bc@meta.data, prefix = "bc_clusters_res.")

bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)
bc_recomputed <- bcRegressOut(bc_recomputed, vars.to.regress = "nFeature_RNA")


bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.8", label.size = 5, pt.size = 1.5)
A375RV_data <- AddMetaData(object = A375RV_data, metadata = bc_recomputed@meta.data)

# Expression UMAP.
bcClusters(bc_recomputed, UMAP = "Seurat", idents = "RNA_snn_res.0.4", pt.size = 1.5, label = TRUE)


DimPlot(
  A375RV_data,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  group.by = "bc_clusters_res.0.8",
  #shape.by = "RNA_snn_res.0.4",
  raster = FALSE
)



DimPlot(
  A375RV_data,
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
bc <- bcRanks(bc_recomputed, idents = "bc_clusters_res.0.6", extended = FALSE)

#Búsqueda de firmas de fármacos
vemurafenib_IDs <- FindDrugs(bc_recomputed, "vemurafenib")$IDs
vemurafenib_IDs <- c("sig-20547", "sig-3758", "sig-6291")
dabrafenib_IDs <- FindDrugs(bc_recomputed, "dabrafenib")$IDs

trametinib_IDs <- FindDrugs(bc_recomputed, "trametinib")$IDs
cobimetinib_IDs <- FindDrugs(bc_recomputed, "cobimetinib")$IDs

# Patchwork object with all bortezomib plots.
dabrafenib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                           signatures = list(values = dabrafenib_IDs), pt.size = 1.5, idents)

trametinib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                           signatures = list(values = trametinib_IDs), pt.size = 1.5)

cobimetinib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = cobimetinib_IDs), pt.size = 1.5)


functional <- bcSignatures(bc_recomputed, UMAP = "beyondcell", genes = list(values=A375RV_data.markers$gene), pt.size = 1.5)
# Plot all three plots in one image.
wrap_plots(dabrafenib, ncol = 3, widths = 7)
wrap_plots(trametinib, ncol = 3)
wrap_plots(cobimetinib, ncol = 1)





bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 0 , top = 5)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 1 , top = 5)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 2 , top = 5)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 3 , top = 5)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.6", lvl = 4 , top = 5)




#Sacar clusters terapeuticos con ssc, pensar si los clusters son adecuados en cuanto a número y composición, crear tabla con célula y a qué cluster pertenece




A375RV_data <- AddMetaData(A375RV_data, metadata =bc@meta.data)



#########################################
bc_recomputed_merged <- bcRegressOut(bc_recomputed_merged, vars.to.regress = "nFeature_RNA")
bc_recomputed_merged <- bcUMAP(bc_recomputed_merged, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
clustree(bc_recomputed_merged@meta.data, prefix ="bc_clusters_res.")

bcClusters(bc_recomputed_merged, UMAP = "beyondcell", idents = "bc_clusters_res.1", label.size = 5, pt.size = 1.5)

vemurafenib_skin <- bcSignatures(bc_recomputed_merged, UMAP = "beyondcell", 
                           signatures = list(values = "vemurafenib_skin"), pt.size = 1.5)
vemurafenib_pancancer <- bcSignatures(bc_recomputed_merged, UMAP = "beyondcell", 
                                 signatures = list(values = "vemurafenib_pancancer"), pt.size = 1.5)
wrap_plots(vemurafenib_skin, ncol = 1)
wrap_plots(vemurafenib_pancancer, ncol = 1)





