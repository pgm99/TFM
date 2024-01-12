###Script para integrar objetos Seurat y realizar análisis de expresión single cell, análisis funcional y análisis de susceptibilidad a fármacos###
library("beyondcell")
library("Seurat")
library("clustree")
library("patchwork")
library("pheatmap")
library("dplyr")
library("tidyr")
library("viridis")
library("RColorBrewer")
library("ggplot2")
library("ggalluvial")
library("forcats")
library("VennDetail")
library("UCell")
library("GSEABase")
library("ggpubr")
library("rstatix")
library("corrplot")
library("Hmisc")

##Cargar objetos Seurat a integrar si fuera necesario##
A375S_data <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/A375S.rds")
A375RV_data <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RV/A375RV.rds")
A375RVC_data <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/A375RVC.rds")
A375RVT_data <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVT/A375RVT.rds")


##En el caso de integración, añadir metadata para diferenciar muestras/tratamientos/etc
A375S_data$model = "A375S"
A375RV_data$model = "A375RV"
A375RVC_data$model = "A375RVC"
A375RVT_data$model = "A375RVT"


##Integración de los objetos Seurat, añadiendo un id por modelo
A375_merged <- merge(A375S_data, y = c(A375RV_data, A375RVC_data, A375RVT_data), 
                     add.cell.ids = c("A375S", "A375RV","A375RVC","A375RVT"), project = "A375")
##Guardar y cargar el objeto Seurat integrado

saveRDS(A375_merged_raw, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/seurat_A375_raw.rds")
A375_merged_raw <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/seurat_A375_raw.rds")

saveRDS(A375_merged, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/seurat_A375.rds")
A375_merged <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/seurat_A375.rds")

 A375_merged$model <- A375_merged$type

##Preprocesado para normalizar y escalar 
A375_merged <- NormalizeData(A375_merged, normalization.method = "LogNormalize", scale.factor = 10000)
A375_merged <- ScaleData(A375_merged)
A375_merged <- FindVariableFeatures(A375_merged, selection.method = "vst", nfeatures = 2000)

##Control de calidad de las muestras: Comprobar y filtrar en función de contenido mitocondrial y nFeature##
A375_merged[["percent.mt"]] <- PercentageFeatureSet(A375_merged, pattern = "^MT-")
A375_merged[["percent.rb"]] <- PercentageFeatureSet(A375_merged, pattern = "^RP[SL]")

Idents(object = A375_merged) <- "type"
VlnPlot(A375_merged, features="nCount_RNA")
VlnPlot(A375_merged, features="nFeature_RNA")

plot1 <- FeatureScatter(A375_merged, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(A375_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 4)
plot1 + plot2

A375_merged <- subset(A375_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


##Reducción PCA##
A375_merged <- RunPCA(A375_merged, features = VariableFeatures(object = A375_merged))

##Clustering. Elbowplot para decidir dimensiones, resolución en los intervalos óptimos para estudios single cell, elegir con clustree
ElbowPlot(A375_merged)
A375_merged <- FindNeighbors(A375_merged, dims = 1:15)
A375_merged <- FindClusters(A375_merged, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
clustree(A375_merged)

##Reducción UMAP por modelo y por clusters de expresión. Los parámetros del UMAP pueden ajustarse para mayor resolución o mayor contexto general
A375_merged <- RunUMAP(A375_merged, dims = 1:15, min.dist = 0, n.neighbors = 40)

Idents(object = A375_merged) <- "model"
plot_model <- DimPlot(A375_merged, reduction = "umap", pt.size = 1.2,  label = TRUE, label.size = 5)+
  labs(color = "Cell line model")+theme(legend.text = element_text(size=18), plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                                        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2))
                                       

Idents(object = A375_merged) <- "RNA_snn_res.0.4"
seurat_umap <- DimPlot(A375_merged, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)+
  labs(color = "Expression clusters")

##Comprobación y regresión del efecto en el análisis de la fase celular( )
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
A375_merged<- CellCycleScoring(A375_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
A375_merged <- RunPCA(A375_merged, features = c(s.genes, g2m.genes))
Idents(object = A375_merged) <- "Phase"
umap_seurat <- DimPlot(A375_merged, pt.size = 1.1)
A375_merged <- ScaleData(A375_merged, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(A375_merged))

Idents(object = A375_merged) <- "model"

A375_markers <- FindAllMarkers(A375_merged, only.pos = TRUE, logfc.threshold = 1)
A375_markers_T <- FindAllMarkers(A375_merged, only.pos = TRUE, logfc.threshold = 1, test.use = "t")

DoHeatmap(A375_merged, features = A375_markers$gene)+ guides(color="none")
A375_markerst <- t(A375_markers)
write.csv(A375_markerst, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/markers_clusters.csv", sep = "\t", quote = F)
write.csv(A375_markers, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/markers_clusters.csv", sep = "\t", quote = F)



A375S.markers <- FindMarkers(A375_merged, ident.1 = "A375S", min.pct = 0.25, logfc.threshold = 1)
A375RV.markers <- FindMarkers(A375_merged, ident.1 = "A375RV", min.pct = 0.25, logfc.threshold = 1)
A375RVC.markers <- FindMarkers(A375_merged, ident.1 = "A375RVC", min.pct = 0.25, logfc.threshold = 1)
A375RVT.markers <- FindMarkers(A375_merged, ident.1 = "A375RVT", min.pct = 0.25, logfc.threshold = 1)

DoHeatmap(A375_merged, features = c("MARCKS", "WNT5A", "UTS2", "COX2", "NTS", "MMP1", "EGR1", "MALAT1", "SUCNR1", "NEAT1", "OCT4",
                                    "UBE2S", "KRT81", "HMGN2", "FOSL1", "KRT18", "RANBP1", "STMN1", "CDK4"), 
          group.colors = c("#440154", "#21908c", "#fde725", "#e7298a")) + guides(color="none")
DoHeatmap(A375_merged, features = c("JUN", "WNT5A", "PDGFRB", "EGFR", "NRG1", "FGFR1", "AXL"), 
          group.colors = c("#440154", "#21908c", "#fde725", "#e7298a")) + guides(color="none")

DoHeatmap(A375_merged, features = genes_of_interest)


DoHeatmap(A375S.markers)

DimPlot(
  A375_merged,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  group.by = c("RNA_snn_res.0.6", "subclone")
)

A375RVC_data <- AddModuleScore(A375RVC_data, features = list(c("JUN", "WNT5A", "PDGFRB", "EGFR", "NRG1", "FGFR1", "AXL")), name = "resistance_markers", nbin = 1)

genes_of_interest <- c("JUN", "WNT5A", "PDGFRB", "EGFR", "NRG1", "FGFR1", "AXL")
BRAF_counts <- colSums(A375S_data@assays[["RNA"]]@data[genes_of_interest, ])
A375_merged@meta.data$BRAF_inhibitor_resistant_biomarkers <- BRAF_counts
FeaturePlot(A375_merged, features = "MAPK1", pt.size = 1.5, cols = c("grey", "yellow"))
FeaturePlot(A375_merged, features= "resistance_markers1", pt.size = 1.3, cols = c("grey", "yellow"),  label = TRUE, repel = TRUE)
bcClusters(A375_merged, UMAP = "beyondcell", idents = "resistance_markers1", pt.size = 1.5, label = TRUE)

FeaturePlot(A375_merged, features = "MAP2K2", pt.size = 2, cols = c("black", "yellow"))

# feature_values <- FetchData(A375_merged, vars = c("MARCKS", "WNT5A"))
# feature_values_long <- pivot_longer(data = feature_values, cols = c("MARCKS", "WNT5A"), names_to = "biomarker", values_to = "value")
# 
# ggplot(feature_values_long, aes(x = biomarker, y = value)) +
#   geom_boxplot(show.legend = FALSE, outlier.shape = NA, coef=0) +
#   labs(title = "Boxplot of gene_of_interest")


###CORRELATION###
# WNT5A <- matrix(A375_merged@assays$RNA@scale.data["WNT5A", ], )

gs_functional <- GenerateGenesets(x = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/data/genesets/gs_functional2.gmt")
bc_functional <- bcScore(A375RVC_data, gs_functional, expr.thres = 0.05)
bc_filtered_functional <- bcSubset(bc_functional, nan.cells = 0.95)
bc_filtered_functional@normalized[is.na(bc_filtered_functional@normalized)] <- 0                                                                                                                                                                       
bc_recomputed_functional <- bcRecompute(bc_filtered_functional, slot = "normalized") 


bc_recomputed_functional <- bcUMAP(bc_recomputed_functional, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
bc_recomputed_functional <- bcRegressOut(bc_recomputed_functional, vars.to.regress = "nFeature_RNA", k.neighbors = 40)
bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
             signatures = list(values = "vemurafenib_skin"), pt.size = 1.5)+theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                                                                                  panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
                                                                                  plot.subtitle = element_text(size=15))

resistance_markers <- as.data.frame(A375RVC_data@meta.data$resistance_markers1)
resistance_markers$BCscore <- (bc_functional@normalized["GOBP_MAPK_CASCADE", ])
colnames(resistance_markers) <- c("resistance_markers", "BCscore")

resistance_markers$resistance_markers[resistance_markers$resistance_markers < 0] <- NA
cor_resistance <- cor(x = resistance_markers)
matrix_resistance_markers <- as.matrix(resistance_markers)
cor_resistanceP <- rcorr(x = matrix_resistance_markers, type = "pearson")

ggplot(resistance_markers, aes(x=resistance_markers, y=BCscore)) + 
  geom_point()+ 
  stat_cor(method = "spearman", label.y = 30, size=7)+ geom_smooth(method = "glm", se = FALSE)+ 
  ggtitle("Correlation between GOBP MAPK cascade and BRAF resistance markers in RVC model")+theme(plot.title = element_text(size = 20), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"), 
                                                                          panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1.7), 
                                                                          plot.subtitle = element_text(size=15))
  


# normalized <- as.data.frame(A375S_data@assays$RNA@data[genes_of_interest, ])
# normalized <- t(normalized)
# normalized <- as.data.frame(normalized)
# normalized$BCscore <- (bc_recomputed_functional@normalized["vemurafenib_skin", ])
# res_all_normalized <- cor(x = normalized, method = "pearson")

# corrplot(res_all_normalized, type = "upper", order = "hclust", 
#          tl.col = "black", tl.srt = 45)
# col<- colorRampPalette(c("blue", "white", "red"))(20)
# heatmap(x = res_all_normalized, col = col, symm = TRUE)
# matrix_normalized <- as.matrix(normalized)
# res_all_normalizedP <- rcorr(x = matrix_normalized, type = "pearson")
# corrplot(res_all_normalizedP$r, p.mat = res_all_normalizedP$P, sig.level = 0.01, type="upper")
# 
# AXL <- ggplot(normalized, aes(x=AXL, y=BCscore)) + 
#   geom_point()
# JUN <- ggplot(normalized, aes(x=JUN, y=BCscore)) + 
#   geom_point()
# WNT5A <- ggplot(normalized, aes(x=WNT5A, y=BCscore)) + 
#   geom_point()
# PDGFRB <- ggplot(normalized, aes(x=PDGFRB, y=BCscore)) + 
#   geom_point()
# EGFR <- ggplot(normalized, aes(x=EGFR, y=BCscore)) + 
#   geom_point()
# NRG1 <- ggplot(normalized, aes(x=NRG1, y=BCscore)) + 
#   geom_point()
# FGFR1 <- ggplot(normalized, aes(x=FGFR1, y=BCscore)) + 
#   geom_point()
# 
# AXL+JUN+WNT5A+PDGFRB+EGFR+NRG1+FGFR1

################BEYONDCELL#################
SSc <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"), )

bc <- bcScore(A375_merged, SSc, expr.thres = 0.1)

saveRDS(bc, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_A375_0.1RAW.rds")
bc <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_A375_0.1RAW.rds")

bc_filtered <- bcSubset(bc, nan.cells = 0.95)
bc_filtered@normalized[is.na(bc_filtered@normalized)] <- 0                                                                                                                                                                       
bc_recomputed <- bcRecompute(bc_filtered, slot = "normalized") 


bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = NULL, k.neighbors = 40)
bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)

saveRDS(bc_recomputed, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_A375.rds")
bc_recomputed <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_A375.rds")
bc_recomputed <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_global_0.1.rds")

bc_recomputeda <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_A375a.rds") ##Este es sin filtrar con nan.cells=0.95

bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)
bc_recomputed <- bcRegressOut(bc_recomputed, vars.to.regress = "nFeature_RNA", k.neighbors = 40)

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)

clustree(bc_recomputed@meta.data, prefix ="bc_clusters_res.")
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.4", label.size = 5, pt.size = 1.5)
A375_merged <- AddMetaData(object = A375_merged, metadata = bc_recomputed@meta.data)

bcClusters(bc_recomputed, UMAP = "Seurat", idents = "RNA_snn_res.0.4", pt.size = 1.5, label = TRUE)
DimPlot(object = A375_merged)

DimPlot(
  A375_merged,
  reduction = "beyondcell",
  label = TRUE,
  label.size = 5,
  group.by = "bc_clusters_res.0.4",
  shape.by = "type",
  raster = FALSE,
  pt.size = 0.8
)

bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "type", pt.size = 1.5)






vemurafenib_IDs <- FindDrugs(bc_recomputed, "vemurafenib")$IDs
dabrafenib_IDs <- FindDrugs(bc_recomputed, "dabrafenib")$IDs
trametinib_IDs <- FindDrugs(bc_recomputed, "trametinib")$IDs
cobimetinib_IDs <- FindDrugs(bc_recomputed, "cobimetinib")$IDs


dabrafenib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                           signatures = list(values = "sig-21060"), pt.size = 1.5)
dabrafenib +

trametinib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                           signatures = list(values = "sig-21070"), pt.size = 1.5)

cobimetinib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = cobimetinib_IDs), pt.size = 1.5)
wrap_plots(dabrafenib, ncol = 1)
wrap_plots(cobitinib, ncol = 1)
wrap_plots(cobimetinib, ncol = 1)


bcSignatures(bc_recomputed, UMAP = "beyondcell", genes = list(values=c("JUN", "WNT5A", "PDGFRB", "EGFR", "NRG1", "FGFR1", "AXL")), 
             signatures = "all", pt.size = 1.5, blend = TRUE)




# Obtain condition-based statistics.
bc_recomputed <- bcRanks(bc_recomputed, idents = "type")
#Obtain unextended therapeutic cluster-based statistics.
bc_recomputed <- bcRanks(bc_recomputed, idents = "bc_clusters_res.0.2", extended = FALSE)


bc4Squares(bc_recomputed, idents = "type", lvl = NULL, top = 3)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.2", lvl = NULL, top = 3)



#####SANKEY#####
samples <- colnames(A375_merged)
clusters <- A375_merged@meta.data$RNA_snn_res.0.4
tc <- A375_merged@meta.data$bc_clusters_res.0.4
models <- A375_merged@meta.data$model
gg_data <- A375_merged@meta.data %>%
  group_by(colnames(A375_merged), clusters, tc, models) %>%
  tally() 


ggplot(data = gg_data,
       aes(axis1 = models, axis2 = clusters)) +
  geom_alluvium(aes(fill = models), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("models", "clusters"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Distribución de models en los clusters")

ggplot(data = gg_data,
       aes(axis1 = type, axis2 = tc)) +
  geom_alluvium(aes(fill = type), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("type", "TCs"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Distribución de type en los TCs")

ggplot(data = gg_data,
       aes(axis1 = models, axis2 = tc)) +
  geom_alluvium(aes(fill = models), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("models", "TCs"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Distribución de type en los TCs")





################FUNCIONAL###################
gs_functional <- GenerateGenesets(x = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/data/genesets/gs_functional2.gmt")
bc_functional <- bcScore(A375_merged, gs_functional, expr.thres = 0.05)
bc_filtered_functional <- bcSubset(bc_functional, nan.cells = 0.95)
bc_filtered_functional@normalized[is.na(bc_filtered_functional@normalized)] <- 0                                                                                                                                                                       
bc_recomputed_functional <- bcRecompute(bc_filtered_functional, slot = "normalized") 
saveRDS(bc_recomputed_functional, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_functional_A375.rds")

bc_recomputed_functional <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_functional_A375.rds")

bc_recomputed_functional <- bcUMAP(bc_recomputed_functional, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
bc_recomputed_functional <- bcRegressOut(bc_recomputed_functional, vars.to.regress = "nFeature_RNA", k.neighbors = 40)
bc_recomputed_functional <- bcUMAP(bc_recomputed_functional, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)


bcClusters(bc_recomputed_functional, UMAP = "beyondcell", idents = "type", label.size = 5, pt.size = 1.5)+ labs(title = "UMAP by sample")

data <- data.frame(bc_recomputed_functional@normalized)
clusters <- data.frame(A375_merged$bc_clusters_res.0.6, A375_merged$RNA_snn_res.0.6)
model$model <- data.frame(A375_merged$model)
clusters <- setNames(clusters, nm = c("BC_clusters", "clusters"))
model <- setNames(type, nm= "model")
data <- t(data)
data <- data[order(model$A375_merged.model), , drop= FALSE]
data <- t(data)
data <- as.matrix(data)


color = colorRampPalette(c("blue", "white", "red"))(200)

pheatmap(data, scale = "row", 
         show_colnames = F, cluster_rows = F, cluster_cols = F, annotation_col = model, treeheight_col = 0, 
         legend = T,  annotation_names_col = T, annotation_legend = T, color = color)


data_scaled <- data.frame(bc_recomputed_functional@scaled)
data_scaled <- t(data_scaled)
data_scaled <- data_scaled[order(model$A375_merged.model), , drop= FALSE]
data_scaled <- t(data_scaled)
data_scaled <- as.matrix(data_scaled)
data_scaled_breaks <- seq(min(data_scaled), max(data_scaled))

pheatmap(data_scaled, scale = "none",
         show_colnames = F, cluster_rows = F, cluster_cols = F, annotation_col = model, treeheight_col = 0, 
         legend = T , annotation_names_col = T, annotation_legend = T, 
         color = color)

plot1 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "GOBP_NEGATIVE_REGULATION_OF_MAPK_CASCADE"), pt.size = 0.8)

plot2 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "GOBP_POSITIVE_REGULATION_OF_MAPK_CASCADE"), pt.size = 0.8)

plot3 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "HALLMARK_ALLOGRAFT_REJECTION"), pt.size = 0.8)

plot4 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "HALLMARK_APOPTOSIS"), pt.size = 1)

plot5 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_DNA_REPAIR"), pt.size = 0.8)

plot6 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "HALLMARK_FATTY_ACID_METABOLISM"), pt.size = 0.8)

plot7 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "HALLMARK_MYC_TARGETS_V1"), pt.size = 0.8)

plot8 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"), pt.size = 0.8)

plot9 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_P53_PATHWAY"), pt.size = 0.8)

resumen_plots <- wrap_plots(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9)

plot11 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "HALLMARK_PEROXISOME"), pt.size = 0.8)

plot12 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"), pt.size = 0.8)

plot13 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "HALLMARK_TGF_BETA_SIGNALING"), pt.size = 0.8)

plot14 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "HALLMARK_UV_RESPONSE"), pt.size = 0.8)

plot15 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "HALLMARK_XENOBIOTIC_METABOLISM"), pt.size = 0.8)

plot16 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "HALLMARK_HYPOXIA"), pt.size = 0.8)

plot17 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "MY_GENESET_S"), pt.size = 0.8)

plot18 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "MY_GENESET_RV"), pt.size = 0.8)

plot19 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "MY_GENESET_RVC"), pt.size = 0.8)

plot20 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "MY_GENESET_RVT"), pt.size = 0.8)

plot21 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"), pt.size = 0.8)

plot30 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"), pt.size = 0.8)

plot31 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "HALLMARK_MTORC1_SIGNALING"), pt.size = 0.8)




resumen_plots2 <-wrap_plots(plot11, plot12, plot13, plot14, plot15, plot16 )
resumen_my_genesets <- wrap_plots(plot17, plot18, plot19, plot20)


plot21 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "GOBP_DEDIFFERENTIATION"), pt.size = 0.8)

plot22 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "MIKKELSEN_DEDIFFERENTIATED_STATE"), pt.size = 0.8)

plot23 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "MIKKELSEN_PLURIPOTENT_STATE"), pt.size = 0.8)

plot24 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "MIKKELSEN_PARTIALLY_REPROGRAMMED_TO_PLURIPOTENCY"), pt.size = 0.8)

plot25 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "PID_PI3KCI_AKT_PATHWAY"), pt.size = 0.8)

plot26 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "BIOCARTA_PTEN_PATHWAY"), pt.size = 0.8)

plot27 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                      signatures = list(values = "GOBP_MAPK_CASCADE"), pt.size = 0.8)

plot28 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "BIOCARTA_MAPK_PATHWAY"), pt.size = 0.8)

plot29 <- bcSignatures(bc_recomputed_functional, UMAP = "Seurat", 
                       signatures = list(values = "PID_PI3KCI_PATHWAY"), pt.size = 0.8)


resumen_plots3 <- wrap_plots(plot21, plot22, plot23, plot24, plot25, plot26)
resumen_plots_MAPK <- wrap_plots(plot25, plot26, plot27, plot28, plot29)

####Plots Vemurafenib####

plot32 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "vemurafenib_skin"), pt.size = 1.5)

plot33 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "vemurafenib_pancancer"), pt.size = 0.8)

plot34 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "vemurafenib_skin"), pt.size = 1.5)

plot35 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "vemurafenib_pancancer"), pt.size = 0.8)

resumen_plots_vemura <- wrap_plots(plot32, plot33)

plot32 +theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
              panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
              plot.subtitle = element_text(size=15))


#####BOXPLOTS#########
boxplot_bc <- CreatebcObject(bc = bc_recomputed_functional)
MAP2K1_bc <- bcSubset(bc = boxplot_bc, signatures = c("REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR"))
MAP2K1_frame <- data.frame(MAP2K1_bc@normalized)
rownames(MAP2K1_frame) <- c("REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR")
MAP2K1_frame <- data.frame(t(MAP2K1_frame))
model <- data.frame(A375_merged$model)
colnames(model) <- "model"
TCs <- data.frame(A375_merged$bc_clusters_res.0.4)
colnames(TCs) <- "TCs"
MAP2K1_frame_model <- merge(x = MAP2K1_frame, y = model, by= 0)
MAP2K1_frame_TCs <- merge(x = MAP2K1_frame_model, y = TCs$TCs, by= 0)
MAP2K1_frame_TCs <- MAP2K1_frame_TCs[,-1]
colnames(MAP2K1_frame_TCs) <- c("cells", "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR", "model", "TCs")
MAP2K1_frame_long <- pivot_longer(data = MAP2K1_frame_TCs, cols = c("REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR"), names_to = "sig", values_to = "enrichment_score")

MAP2K1_frame_long$model <- factor(MAP2K1_frame_long$model, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
TCs <- factor(MAP2K1_frame_long$TCs, levels=c("0", "1", "2", "3", "4", "5"))



MAP2K1_frame_long$model <- as.factor(MAP2K1_frame_long$model)
anova.test <- MAP2K1_frame_long %>% anova_test(formula = enrichment_score~model)
# group_by(MAP2K1_frame_long$sig)%>%


pwc <- MAP2K1_frame_long %>% tukey_hsd(enrichment_score ~ model)
pwc
# options("scipen" = 999, "digits" = 22)
# options("scipen" = 0, "digits" = 7)



boxplot_funcional_A375 <- ggplot(MAP2K1_frame_long) +
  geom_boxplot((aes(x=model, y=enrichment_score, fill=model)), show.legend = FALSE, outlier.shape = NA, coef=0) + 
  stat_pvalue_manual(pwc, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE, y.position = c(9, 9.5, 10, 10.5, 11, 11.5), size = 5) +
  ggtitle("Comparison of UPR signature between models")+
  theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
        panel.grid = element_line(colour = "grey", linewidth = 0.25), plot.margin = margin(t=10, r = 10, b = 5, l = 5)) + ylim(0, 14)

boxplot_funcional_A375 + scale_fill_manual(values = c("#440154", "#21908c", "#fde725", "#e7298a"))


###BOXPLOT TRAMETINIB###

trametinib <- bcSubset(bc = bc, signatures = c("sig-21070", "sig-20908", "sig-21071"))
trame_frame <- data.frame(trametinib@normalized)
rownames(trame_frame) <- c("trametinib_CTRP", "trametinib_GDSC", "trametinib_PRISM")
trame_frame <- data.frame(t(trame_frame))

type <- data.frame(A375_merged$type)
trame_frame_type <- merge(x = trame_frame, y = type, by= 0)
trame_frame_long <- pivot_longer(data = trame_frame_type, cols = c("trametinib_CTRP","trametinib_GDSC","trametinib_PRISM"), names_to = "sig", values_to = "enrichment_score")

sample_trame <- factor(trame_frame_long$A375_merged.type, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
ggplot(trame_frame_long) +
  geom_boxplot((aes(x=sig, y=enrichment_score, fill=sample_trame)))

colSums(is.na(trame_frame))

###BOXPLOT DABRAFENIB###

dabrafenib <- bcSubset(bc = bc, signatures = c("sig-20909", "sig-21059", "sig-21060"))
dabra_frame <- data.frame(dabrafenib@normalized)
rownames(dabra_frame) <- c("dabrafenib_GDSC", "dabrafenib_CTRP", "dabrafenib_PRISM")
dabra_frame <- data.frame(t(dabra_frame))

type <- data.frame(A375_merged$type)
dabra_frame_type <- merge(x = dabra_frame, y = type, by= 0)
dabra_frame_long <- pivot_longer(data = dabra_frame_type, cols = c("dabrafenib_GDSC", "dabrafenib_CTRP", "dabrafenib_PRISM"), names_to = "sig", values_to = "enrichment_score")

sample_dabra <- factor(dabra_frame_long$A375_merged.type, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
ggplot(dabra_frame_long) +
  geom_boxplot((aes(x=sig, y=enrichment_score, fill=sample_dabra)))


###BOXPLOT COBIMETINIB###

cobimetinib <- bcSubset(bc = bc, signatures = c("sig-21037"))
cobi_frame <- data.frame(cobimetinib@normalized)
rownames(cobi_frame) <- c("cobimetinib_PRISM")
cobi_frame <- data.frame(t(cobi_frame))

type <- data.frame(A375_merged$type)
cobi_frame_type <- merge(x = cobi_frame, y = type, by= 0)
cobi_frame_long <- pivot_longer(data = cobi_frame_type, cols = c("cobimetinib_PRISM"), names_to = "sig", values_to = "enrichment_score")

sample_cobi <- factor(cobi_frame_long$A375_merged.type, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
ggplot(cobi_frame_long) +
  geom_boxplot((aes(x=sig, y=enrichment_score, fill=sample_cobi)))

colSums(is.na(cobi_frame))





NEG_REG_MAPK <- read.delim(file = "/home/lpgonzalezm/Downloads/GOBP_NEGATIVE_REGULATION_OF_MAPK_CASCADE.v2023.1.Hs.tsv", header = F, sep = ",", quote = "")
NEG_REG_MAPK <- t(NEG_REG_MAPK)
DoHeatmap(A375_merged, features = NEG_REG_MAPK, slot = "scale.data")

MAPK_CASCADE <- read.delim(file = "/home/lpgonzalezm/Downloads/GOBP_MAPK_CASCADE.v2023.1.Hs.tsv", header = F, sep = ",", quote = "")
MAPK_CASCADE <- t(MAPK_CASCADE)
DoHeatmap(A375_merged, features = MAPK_CASCADE, slot = "scale.data")


MAPK_CASCADE_2 <- read.delim(file = "/home/lpgonzalezm/Downloads/BIOCARTA_MAPK_PATHWAY.v2023.1.Hs.tsv", header = F, sep = ",", quote = "")
MAPK_CASCADE_2 <- t(MAPK_CASCADE_2)
DoHeatmap(A375_merged, features = MAPK_CASCADE_2)




A375S_RVT.markers <- FindMarkers(A375_merged, ident.1 = "A375S", ident.2 = "A375RVT", min.pct = 0.25, logfc.threshold = 1)
DoHeatmap(A375_merged, features = rownames(A375S_RVT.markers))

A375RVC_RVT.markers <- FindMarkers(A375_merged, ident.1 = "A375RVC", ident.2 = "A375RVT", min.pct = 0.25, logfc.threshold = 1)
DoHeatmap(A375_merged, features = rownames(A375RVC_RVT.markers))


NEG_REG_MAPK <- as.data.frame(NEG_REG_MAPK)
MAPK_S_RVT <- intersect(rownames(A375S_RVT.markers), NEG_REG_MAPK)

######firmas dabrafenib skin#########

gs_dabrafenib <- GenerateGenesets(x = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/data/genesets/dabrafenib_signatures.gmt")
bc_dabra <- bcScore(A375_merged, gs_dabrafenib, expr.thres = 0.05) 
#bc_filtered_functional <- bcSubset(bc_functional, nan.cells = 0.95)
bc_dabra@normalized[is.na(bc_dabra@normalized)] <- 0                                                                                                                                                                       
bc_recomputed_dabra <- bcRecompute(bc_dabra, slot = "normalized") 

bc_recomputed_dabra <- bcUMAP(bc_recomputed_dabra, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40, add.DSS = TRUE)
bc_recomputed_dabra <- bcRegressOut(bc_recomputed_dabra, vars.to.regress = "nFeature_RNA", k.neighbors = 40)
bc_recomputed_dabra <- bcUMAP(bc_recomputed_dabra, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)


FeaturePlot(A375_merged, features = "TRIM51")
x <- drugInfo(bc = bc_recomputed, "trametinib")

####DIAGRAMAS VENN####
trametinib_GDSC <- SSc@genelist$`sig-20908`
trametinib_CTRP <- SSc@genelist$`sig-21070`
trametinib_PRISM <- SSc@genelist$`sig-21071`

ven_trameUP <- venndetail(list(GDSC=trametinib_GDSC$up, CTRP=sapply(trametinib_CTRP$up, head, 50), PRISM=sapply(trametinib_PRISM$up, head, 50)))
plot(ven_trameUP)
ven_trameDOWN <- venndetail(list(GDSC=trametinib_GDSC$down, CTRP=trametinib_CTRP$down, PRISM=trametinib_PRISM$down))
plot(ven_trameDOWN)


trametinib_GDSC_UP50 <- head(trametinib_GDSC$up, 50)
trametinib_GDSC_DOWN50 <- head(trametinib_GDSC$down, 50)
trametinib_CTRP_UP50 <- head(trametinib_CTRP$up, 50)
trametinib_CTRP_DOWN50 <- head(trametinib_CTRP$down, 50)
trametinib_PRISM_UP50 <- head(trametinib_PRISM$up, 50)
trametinib_PRISM_DOWN50 <- head(trametinib_PRISM$down, 50)

ven_trameUP50 <- venndetail(list(GDSC=trametinib_GDSC_UP50, CTRP=trametinib_CTRP_UP50, PRISM=trametinib_PRISM_UP50))
plot(ven_trameUP50)
ven_trameDOWN50 <- venndetail(list(GDSC=trametinib_GDSC_DOWN50, CTRP=trametinib_CTRP_DOWN50, PRISM=trametinib_PRISM_DOWN50))
plot(ven_trameDOWN50)


vemurafenib_skin <- gs_functional@genelist$vemurafenib_skin
vemurafenib_pancancer <- gs_functional@genelist$vemurafenib_pancancer

ven_vemuraUP <- venndetail(list(skin=vemurafenib_skin$up, pancancer=vemurafenib_pancancer$up))
plot(ven_vemuraUP)

ven_vemuraDOWN <- venndetail(list(skin=vemurafenib_skin$down, pancancer=vemurafenib_pancancer$down))
plot(ven_vemuraDOWN)

dabrafenib_GDSC <- SSc@genelist$`sig-20909`
dabrafenib_CTRP <- SSc@genelist$`sig-21059`
dabrafenib_PRISM <- SSc@genelist$`sig-21060`

ven_dabraUP <- venndetail(list(GDSC=dabrafenib_GDSC$up, CTRP=dabrafenib_CTRP$up, PRISM=dabrafenib_PRISM$up))
plot(ven_dabraUP)

ven_dabraDOWN <- venndetail(list(GDSC=dabrafenib_GDSC$down, CTRP=dabrafenib_CTRP$down, PRISM=dabrafenib_PRISM$down))
plot(ven_dabraDOWN)


#######UCELL#######

gs_UCELL <- list("MAP2K1" = gs_depend@genelist$MAP2K1_DepMap_22Q4$up,
                 "BRAF" = gs_depend@genelist$BRAF_DepMap_22Q4$up,
                 "MAP2K2" = gs_depend@genelist$MAP2K2_DepMap_22Q4$up,
                 "RAF1" = gs_depend@genelist$RAF1_DepMap_22Q4$up
                 )

BRAF_ucell_set_up <- paste0(gs_depend@genelist$BRAF_DepMap_22Q4$up, "+")
BRAF_ucell_set_dn <- paste0(gs_depend@genelist$BRAF_DepMap_22Q4$down, "-")

MAP2K1_ucell_set_up <- paste0(gs_depend@genelist$MAP2K1_DepMap_22Q4$up, "+")
MAP2K1_ucell_set_dn <- paste0(gs_depend@genelist$MAP2K1_DepMap_22Q4$down, "-")

MAP2K2_ucell_set_up <- paste0(gs_depend@genelist$MAP2K2_DepMap_22Q4$up, "+")
MAP2K2_ucell_set_dn <- paste0(gs_depend@genelist$MAP2K2_DepMap_22Q4$down, "-")

RAF1_ucell_set_up <- paste0(gs_depend@genelist$RAF1_DepMap_22Q4$up, "+")
RAF1_ucell_set_dn <- paste0(gs_depend@genelist$RAF1_DepMap_22Q4$down, "-")


gs_ucell_test <- list("BRAF_bidirect" = c(braf_ucell_set_up, braf_ucell_set_dn), "map2k1")

A375_merged <- AddModuleScore_UCell(A375_merged, features = gs_UCELL)
A375_merged <- AddModuleScore_UCell(A375_merged, features = gs_ucell_test, maxRank = 3000)

featnames <- c("MAP2K1_UCell", "MAP2K2_UCell", "BRAF_UCell", "RAF1_UCell")
FeaturePlot(A375_merged, features = featnames, pt.size = 0.1, ncol = 2, cols = c("blue", "red"))
FeaturePlot(A375_merged, features = "BRAF_bidirect_UCell", pt.size = 0.1)

VlnPlot(A375_merged, features = featnames, pt.size = 0, split.by = "type", ncol = 2)

gs_trametinib <- GenerateGenesets(x = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/data/genesets/trametinib.gmt")
trame_ucell_set_up <- paste0(gs_trametinib@genelist$TRAMETINIB_GS$up, "+")
trame_ucell_set_dn <- paste0(gs_trametinib@genelist$TRAMETINIB_GS$down, "-")

#UCELL TRAMETINIB

gs_trametinib_UCELL <- list("trametinib" = c(trame_ucell_set_up, trame_ucell_set_dn ))

A375_merged <- AddModuleScore_UCell(A375_merged, features = gs_trametinib_UCELL)
FeaturePlot(A375_merged, features = "trametinib_UCell", pt.size = 0.1, cols = c("blue", "red"))

trame_frame_UCELL <- data.frame(A375_merged$trametinib_UCell)
trame_frame_UCELL <- data.frame(t(trame_frame_UCELL))
rownames(trame_frame_UCELL) <- c("trametinib_ucell")

type <- data.frame(A375_merged$type)
trame_frame_UCELL <- merge(x = trame_frame_UCELL, y = type, by= 0)
trame_frame_UCELL_long <- pivot_longer(data = trame_frame_UCELL, cols = c("A375_merged.trametinib_UCell"), names_to = "sig", values_to = "enrichment_score")

sample_trame <- factor(trame_frame_UCELL$A375_merged.type, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
box_trame <-ggplot(trame_frame_UCELL_long) +
  geom_boxplot((aes(x=sig, y=enrichment_score, fill=sample_trame)))

box_trame + stat_compare_means(comparisons = comparisons)

ggboxplot(trame_frame_UCELL_long, x = trame_frame_UCELL_long$sig, y = trame_frame_UCELL_long$enrichment_score, fill = sample_trame)

#UCELL VEMURAFENIB
gs_vemurafenib <- GenerateGenesets(x = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/data/genesets/vemurafenib.gmt")
gs_vemurafenib_UCELL <- list("vemurafenib" = gs_vemurafenib@genelist$VEMURAFENIB_GS)

A375_merged <- AddModuleScore_UCell(A375_merged, features = gs_vemurafenib_UCELL)
FeaturePlot(A375_merged, features = "vemurafenib_UCell", pt.size = 0.1, cols = c("blue", "red"))

vemura_frame_UCELL <- data.frame(A375_merged$vemurafenib_UCell)
vemura_frame_UCELL <- data.frame(t(vemura_frame_UCELL))
rownames(vemura_frame_UCELL) <- c("vemurafenib_ucell")

type <- data.frame(A375_merged$type)
vemura_frame_UCELL <- merge(x = vemura_frame_UCELL, y = type, by= 0)
vemura_frame_UCELL_long <- pivot_longer(data = vemura_frame_UCELL, cols = c("vemurafenib_ucell"), names_to = "sig", values_to = "enrichment_score")

sample_vemura <- factor(vemura_frame_UCELL$A375_merged.type, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
ggplot(vemura_frame_UCELL_long) +
  geom_boxplot((aes(x=sig, y=enrichment_score, fill=sample_vemura)))


##########DEPENDENCIAS DRUGS###################
gs_depend <- GenerateGenesets(x = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/data/genesets/pablo_signatures_2.gmt")
bc_depend <- bcScore(A375_merged, gs_depend, expr.thres = 0.1) 
bc_filtered_depend <- bcSubset(bc_depend, nan.cells = 0.95)
bc_filtered_depend@normalized[is.na(bc_filtered_depend@normalized)] <- 0                                                                                                                                                                       
bc_recomputed_depend <- bcRecompute(bc_filtered_depend, slot = "normalized") 
saveRDS(bc_recomputed_depend, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_depend_A375.rds")
bc_recomputed_depend <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_depend_A375.rds")


bc_recomputed_depend <- bcUMAP(bc_recomputed_depend, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40, add.DSS = TRUE)
bc_recomputed_depend <- bcRegressOut(bc_recomputed_depend, vars.to.regress = "nFeature_RNA", k.neighbors = 40)





plot35 <- bcSignatures(bc_recomputed_depend, UMAP = "Seurat", 
                      signatures = list(values = "RAF1_DepMap_22Q4"), pt.size = 0.8)
plot36 <- bcSignatures(bc_recomputed_depend, UMAP = "Seurat", 
                       signatures = list(values = "MAP2K1_DepMap_22Q4"), pt.size = 0.8)
plot37 <- bcSignatures(bc_recomputed_depend, UMAP = "Seurat", 
                       signatures = list(values = "BRAF_DepMap_22Q4"), pt.size = 0.8)
plot38 <- bcSignatures(bc_recomputed_depend, UMAP = "Seurat", 
                       signatures = list(values = "MAP2K2_DepMap_22Q4"), pt.size = 0.8)

plot39 <- bcSignatures(bc_recomputed_depend, UMAP = "Seurat", 
                       signatures = list(values = "PI3K_depmap"), pt.size = 0.8)


plot_type <- plot_type + labs(title = "UMAP_sample")

resumen_plots_dep <- wrap_plots(plot36, plot37, plot38, plot39)

dep_frame <- data.frame(bc_depend@normalized)
rownames(dep_frame) <- c("RAF1","MAP2K1","BRAF","MAP2K2")
dep_frame <- data.frame(t(dep_frame))

type <- data.frame(A375_merged$type)
dep_frame <- merge(x = dep_frame, y = type, by= 0)

dep_frame_longer <- pivot_longer(data = dep_frame, cols = c("RAF1","MAP2K1","BRAF","MAP2K2"), names_to = "sig", values_to = "enrichment_score")

sample <- factor(dep_frame_longer$A375_merged.type, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
box_dep <- ggplot(dep_frame_longer) + 
  geom_boxplot(aes(x=sig, y=enrichment_score, fill=sample))

comparisons <- list( c("A375S", "A375RV"), c("A375S", "A375RVC"), c("A375S", "A375RVT"))
stat_compare_means(comparisons = comparisons)



##############BC GLOBAL##############

bc_functional_depend <- bcMerge(bc1 =bc_recomputed_depend , bc2 = bc_recomputed_functional)

bc_recomputed_global <- bcMerge(bc1 = bc_recomputed, bc2 = bc_recomputed_functional, keep.bc.clusters = TRUE)
bcClusters(bc_recomputed_global, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)

bc_recomputed_global <- bcRegressOut(bc_recomputed_global, vars.to.regress = "nFeature_RNA", k.neighbors = 40)
bc_recomputed_global <- bcUMAP(bc_recomputed_global, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
#umaps generales
bc_seurat <- bcClusters(bc_recomputed_global, UMAP = "Seurat", idents = "type", label.size = 5, pt.size = 1.5)
bc_type <- bcClusters(bc_recomputed_global, UMAP = "beyondcell", idents = "type", label.size = 5, pt.size = 1.5)
bc_umap <- bcClusters(bc_recomputed_global, UMAP = "beyondcell", idents = "bc_clusters_res.0.4", label.size = 5, pt.size = 1.5)

saveRDS(bc_recomputed_global, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_global_0.1.rds")
bc_recomputed_global <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_global_0.1.rds")

#DRUGS
dabrafenib_IDs <- FindDrugs(bc_recomputed_global, "dabrafenib")$IDs
trametinib_IDs <- FindDrugs(bc_recomputed_global, "trametinib")$IDs
cobimetinib_IDs <- FindDrugs(bc_recomputed_global, "cobimetinib")$IDs
FindDrugs(bc_recomputed_global, "olaparib")$IDs

dabrafenib <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                           signatures = list(values = "sig-21060"), pt.size = 2)

dabrafenib +theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                  panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
                  plot.subtitle = element_text(size=15))

trametinib <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                           signatures = list(values = trametinib_IDs), pt.size = 1.5)

cobimetinib <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                            signatures = list(values = cobimetinib_IDs), pt.size = 1.5)
wrap_plots(dabrafenib, ncol = 3)
wrap_plots(trametinib, ncol = 3)
wrap_plots(cobimetinib, ncol = 1)

dabrafenib+trametinib+cobimetinib+bc_umap
dabrafenib+trametinib+cobimetinib+seurat_umap

#Boxplot Trametinib
trametinib_bc <- CreatebcObject(bc = bc_recomputed_global)
trametinib_bc <- bcSubset(bc = trametinib_bc, signatures = c("sig-21060"))
trame_frame <- data.frame(trametinib_bc@normalized)
rownames(trame_frame) <- c("dabrafenib_PRISM")
trame_frame <- data.frame(t(trame_frame))

model <- data.frame(A375_merged$model)
trame_frame_type <- merge(x = trame_frame, y = model, by= 0)
trame_frame_long <- pivot_longer(data = trame_frame_type, cols = c("dabrafenib_PRISM"), names_to = "sig", values_to = "enrichment_score")
colnames(trame_frame_long) <- c("Row.names", "model", "sig", "enrichment_score")


sample <- factor(trame_frame_long$model, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
ggplot(trame_frame_long) +
  geom_boxplot((aes(x=sample, y=enrichment_score, fill=sample)), show.legend = FALSE, outlier.shape = NA)+ 
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, hide.ns = FALSE, y.position = c(8,9,10,11,12,13)) +
  ggtitle("Comparison of dabrafenib signature between cell lines models")

stat.test <- trame_frame_long %>% t_test((enrichment_score~model))
stat.test <- stat.test %>% add_xy_position(x = "model")

#Boxplot Dabrafenib
dabrafenib_bc <- CreatebcObject(bc = bc_recomputed_global)
dabrafenib_bc <- bcSubset(bc = dabrafenib_bc, signatures = c("sig-21060"))
dabra_frame <- data.frame(dabrafenib_bc@normalized)
rownames(dabra_frame) <- c("dabrafenib_PRISM")
dabra_frame <- data.frame(t(dabra_frame))

model <- data.frame(A375_merged$model)
dabra_frame_type <- merge(x = dabra_frame, y = model, by= 0)
dabra_frame_long <- pivot_longer(data = dabra_frame_type, cols = c("dabrafenib_PRISM"), names_to = "sig", values_to = "enrichment_score")
colnames(dabra_frame_long) <- c("Row.names", "model", "sig", "enrichment_score")


sample <- factor(dabra_frame_long$model, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
ggplot(dabra_frame_long) +
  geom_boxplot((aes(x=sample, y=enrichment_score, fill=sample)), show.legend = FALSE, outlier.shape = NA)+ 
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, hide.ns = FALSE, y.position = c(8,9,10,11,12,13)) +
  ggtitle("Comparison of dabrafenib signature between cell lines models")

stat.test <- dabra_frame_long %>% t_test((enrichment_score~model))
stat.test <- stat.test %>% add_xy_position(x = "model")



#Boxplot Cobimetinib
cobimetinib_bc <- CreatebcObject(bc = bc_recomputed_global)
cobimetinib_bc <- bcSubset(bc = cobimetinib_bc, signatures = "sig-21037")
cobi_frame <- data.frame(cobimetinib_bc@normalized)
rownames(cobi_frame) <- "cobimetinib_PRISM"
cobi_frame <- data.frame(t(cobi_frame))

model <- data.frame(A375_merged$model)
cobi_frame_type <- merge(x = cobi_frame, y = model, by= 0)
cobi_frame_long <- pivot_longer(data = cobi_frame_type, cols = c("cobimetinib_PRISM"), names_to = "sig", values_to = "enrichment_score")
colnames(cobi_frame_long) <- c("Row.names", "model", "sig", "enrichment_score")


sample <- factor(cobi_frame_long$model, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
ggplot(cobi_frame_long) +
  geom_boxplot((aes(x=sample, y=enrichment_score, fill=sample)), show.legend = FALSE, outlier.shape = NA)+ 
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, hide.ns = FALSE, ) +
  ggtitle("Comparison of cobimetinib signature between cell lines models")

stat.test <- cobi_frame_long %>% t_test((enrichment_score~model))
stat.test <- stat.test %>% add_xy_position(x = "model")

#Boxplot Global
drugs_bc <- CreatebcObject(bc = bc_recomputed_global)
drugs_bc <- bcSubset(bc = drugs_bc, signatures = c("sig-21060","sig-20908", "sig-21037"))
drugs_frame <- data.frame(drugs_bc@normalized)
rownames(drugs_frame) <- c("dabrafenib_PRISM","trametinib_GDSC", "cobimetinib_PRISM")
drugs_frame <- data.frame(t(drugs_frame))

type <- data.frame(A375_merged$type)
drugs_frame_type <- merge(x = drugs_frame, y = type, by= 0)
drugs_frame_long <- pivot_longer(data = drugs_frame_type, cols = c("dabrafenib_PRISM","trametinib_GDSC", "cobimetinib_PRISM"), names_to = "sig", values_to = "enrichment_score")

sample <- factor(drugs_frame_long$A375_merged.type, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
ggplot(drugs_frame_long) +
  geom_boxplot((aes(x=sig, y=enrichment_score, fill=sample)))

#FUNCIONAL


plot1 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                      signatures = list(values = "GOBP_NEGATIVE_REGULATION_OF_MAPK_CASCADE"), pt.size = 0.8)

plot2 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                      signatures = list(values = "GOBP_POSITIVE_REGULATION_OF_MAPK_CASCADE"), pt.size = 0.8)

plot3 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_ALLOGRAFT_REJECTION"), pt.size = 0.8)

plot4 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_APOPTOSIS"), pt.size = 0.8)

plot5 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_DNA_REPAIR"), pt.size = 0.8)

plot6 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_FATTY_ACID_METABOLISM"), pt.size = 0.8)

plot7 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_MYC_TARGETS_V1"), pt.size = 0.8)

plot8 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"), pt.size = 0.8)

plot9 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_P53_PATHWAY"), pt.size = 0.8)

resumen_plots <- wrap_plots(plot1, plot2, plot3, plot4, plot6, plot5, plot7, plot8, plot9)
`plots_diapos` <- wrap_plots(plot1, plot3, plot4, plot9, plot14, plot15)

plot11 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_PEROXISOME"), pt.size = 0.8)

plot12 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"), pt.size = 0.8)

plot13 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_TGF_BETA_SIGNALING"), pt.size = 0.8)

plot14 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_UV_RESPONSE"), pt.size = 0.8)

plot15 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_XENOBIOTIC_METABOLISM"), pt.size = 0.8)

plot16 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_HYPOXIA"), pt.size = 0.8)

plot17 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "MY_GENESET_S"), pt.size = 0.8)

plot18 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "MY_GENESET_RV"), pt.size = 0.8)

plot19 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "MY_GENESET_RVC"), pt.size = 0.8)

plot20 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "MY_GENESET_RVT"), pt.size = 0.8)

plot21 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"), pt.size = 0.8)

plot30 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"), pt.size = 0.8)

plot31 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_MTORC1_SIGNALING"), pt.size = 0.8)




resumen_plots2 <-wrap_plots(plot11, plot12, plot13, plot14, bc_umap, plot15, plot16, plot21, plot22)
resumen_my_genesets <- wrap_plots(plot17, plot18, plot19, plot20)


plot21 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "GOBP_DEDIFFERENTIATION"), pt.size = 0.8)

plot22 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "MIKKELSEN_DEDIFFERENTIATED_STATE"), pt.size = 0.8)

plot23 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "MIKKELSEN_PLURIPOTENT_STATE"), pt.size = 0.8)

plot24 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "MIKKELSEN_PARTIALLY_REPROGRAMMED_TO_PLURIPOTENCY"), pt.size = 0.8)

plot25 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "PID_PI3KCI_AKT_PATHWAY"), pt.size = 0.8)

plot26 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "BIOCARTA_PTEN_PATHWAY"), pt.size = 0.8)

plot27 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "GOBP_MAPK_CASCADE"), pt.size = 0.8)

plot28 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "BIOCARTA_MAPK_PATHWAY"), pt.size = 0.8)

plot29 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "PID_PI3KCI_PATHWAY"), pt.size = 0.8)


resumen_plots3 <- wrap_plots(plot23, plot24, bc_umap)
resumen_plots_MAPK <- wrap_plots(plot26, plot27, plot28, plot29)



#Dependencias

plot40 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "RAF1_DepMap_22Q4"), pt.size = 0.8)
plot41 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "MAP2K1_DepMap_22Q4"), pt.size = 1.8)
plot42 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "BRAF_DepMap_22Q4"), pt.size = 0.8)
plot43 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "MAP2K2_DepMap_22Q4"), pt.size = 1.8)

resumen_plots_dep_umap <- wrap_plots(plot40, plot41, plot42, plot43)
plot41 + bc_type
plot41 + bc_umap

plot44 <- bcSignatures(bc_recomputed_global, UMAP = "beyondcell", 
                       signatures = list(values = "vemurafenib_skin"), pt.size = 0.8)











boxplot_bc <- CreatebcObject(bc = bc_recomputed_global)
boxplot_bc <- bcSubset(bc = boxplot_bc, signatures = c("MAP2K1_DepMap_22Q4", "GOBP_NEGATIVE_REGULATION_OF_MAPK_CASCADE",
                                                       "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_APOPTOSIS", 
                                                       "HALLMARK_P53_PATHWAY", "HALLMARK_UV_RESPONSE", "HALLMARK_XENOBIOTIC_METABOLISM"))
boxplot_frame <- data.frame(boxplot_bc@normalized)
rownames(boxplot_frame) <- c("MAP2K1_DepMap_22Q4", "NEG_REGULATION_MAPK_CASCADE",
                             "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_APOPTOSIS", 
                             "HALLMARK_P53_PATHWAY", "HALLMARK_UV_RESPONSE", "HALLMARK_XENOBIOTIC_METABOLISM")
boxplot_frame <- data.frame(t(boxplot_frame))
model <- data.frame(A375_merged$model)
colnames(model) <- "model"
boxplot_frame_model <- merge(x = boxplot_frame, y = model, by= 0)
boxplot_frame_long <- pivot_longer(data = boxplot_frame_model, cols = c("MAP2K1_DepMap_22Q4", "NEG_REGULATION_MAPK_CASCADE",
                                                                        "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_APOPTOSIS", 
                                                                        "HALLMARK_P53_PATHWAY", "HALLMARK_UV_RESPONSE", "HALLMARK_XENOBIOTIC_METABOLISM"), names_to = "sig", values_to = "enrichment_score")

model <- factor(boxplot_frame_long$model, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
ggplot(boxplot_frame_long) +
  geom_boxplot((aes(x=sig, y=enrichment_score, fill=model)))

###Boxplots funcionales individuales###

boxplot_bc <- CreatebcObject(bc = bc_recomputed_global)
MAP2K1_bc <- bcSubset(bc = boxplot_bc, signatures = c("MAP2K1_DepMap_22Q4"))
MAP2K1_frame <- data.frame(MAP2K1_bc@normalized)
rownames(MAP2K1_frame) <- c("MAP2K1_DepMap_22Q4")
MAP2K1_frame <- data.frame(t(MAP2K1_frame))
model <- data.frame(A375_merged$model)
colnames(model) <- "model"
TCs <- data.frame(A375_merged$bc_clusters_res.0.4)
colnames(TCs) <- "TCs"
MAP2K1_frame_model <- merge(x = MAP2K1_frame, y = model, by= 0)
MAP2K1_frame_TCs <- merge(x = MAP2K1_frame_model, y = TCs$TCs, by= 0)
MAP2K1_frame_TCs <- MAP2K1_frame_TCs[,-1]
colnames(MAP2K1_frame_TCs) <- c("cells", "MAP2K1_DepMap_22Q4", "model", "TCs")
MAP2K1_frame_long <- pivot_longer(data = MAP2K1_frame_TCs, cols = c("MAP2K1_DepMap_22Q4"), names_to = "sig", values_to = "enrichment_score")

MAP2K1_frame_long$model <- factor(MAP2K1_frame_long$model, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
TCs <- factor(MAP2K1_frame_long$TCs, levels=c("0", "1", "2", "3", "4", "5"))



MAP2K1_frame_long$model <- as.factor(MAP2K1_frame_long$model)
anova.test <- MAP2K1_frame_long %>% anova_test(formula = enrichment_score~model)
# group_by(MAP2K1_frame_long$sig)%>%
   

pwc <- MAP2K1_frame_long %>% tukey_hsd(enrichment_score ~ model)
pwc
# options("scipen" = 999, "digits" = 22)
# options("scipen" = 0, "digits" = 7)



boxplot_funcional_A375 <- ggplot(MAP2K1_frame_long) +
  geom_boxplot((aes(x=model, y=enrichment_score, fill=model)), show.legend = FALSE, outlier.shape = NA, coef=0) + 
  stat_pvalue_manual(pwc, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE, y.position = c(12, 13.5, 15, 16.5, 18, 19.5), size = 5) +
  ggtitle("Comparison of MAP2K1 depmap signature between models")+
  theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
        panel.grid = element_line(colour = "grey", linewidth = 0.25), plot.margin = margin(t=10, r = 10, b = 5, l = 5)) + ylim(-10, 20)
  
boxplot_funcional_A375 + scale_fill_manual(values = c("#440154", "#21908c", "#fde725", "#e7298a"))


# # stat_pvalue_manual(pwc, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE, y.position = c(19, 22, 25)) +
# labs(subtitle = get_test_label(anova.test, detailed = TRUE))+



boxplot_bc <- CreatebcObject(bc = bc_prueba_global)
MAP2K1_bc <- bcSubset(bc = boxplot_bc, signatures = c("PID_PI3KCI_AKT_PATHWAY"))
MAP2K1_frame <- data.frame(MAP2K1_bc@normalized)
rownames(MAP2K1_frame) <- c("PID_PI3KCI_AKT_PATHWAY")
MAP2K1_frame <- data.frame(t(MAP2K1_frame))
model <- data.frame(A375_merged$model)
colnames(model) <- "model"
TCs <- data.frame(A375_merged$bc_clusters_res.0.4)
colnames(TCs) <- "TCs"
MAP2K1_frame_model <- merge(x = MAP2K1_frame, y = model, by= 0)
MAP2K1_frame_TCs <- merge(x = MAP2K1_frame_model, y = TCs$TCs, by= 0)
MAP2K1_frame_TCs <- MAP2K1_frame_TCs[,-1]
colnames(MAP2K1_frame_TCs) <- c("cells", "PID_PI3KCI_AKT_PATHWAY", "model", "TCs")
MAP2K1_frame_long <- pivot_longer(data = MAP2K1_frame_TCs, cols = c("PID_PI3KCI_AKT_PATHWAY"), names_to = "sig", values_to = "enrichment_score")

MAP2K1_frame_long$model <- factor(MAP2K1_frame_long$model, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
TCs <- factor(MAP2K1_frame_long$TCs, levels=c("0", "1", "2", "3", "4", "5"))



MAP2K1_frame_long$model <- as.factor(MAP2K1_frame_long$model)
anova.test <- MAP2K1_frame_long %>% anova_test(formula = enrichment_score~model)
# group_by(MAP2K1_frame_long$sig)%>%


pwc <- MAP2K1_frame_long %>% tukey_hsd(enrichment_score ~ model)
pwc
# options("scipen" = 999, "digits" = 22)
# options("scipen" = 0, "digits" = 7)



boxplot_funcional_A375 <- ggplot(MAP2K1_frame_long) +
  geom_boxplot((aes(x=model, y=enrichment_score, fill=model)), show.legend = FALSE, outlier.shape = NA) + 
  stat_pvalue_manual(pwc, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE, y.position = c(13, 16, 19, 22, 25, 28), size = 5) +
  ggtitle("Comparison of PI3KCI_AKT signature between models")+
  theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
        panel.grid = element_line(colour = "grey", linewidth = 0.25), plot.margin = margin(t=10, r = 10, b = 5, l = 5)) + ylim(-3, 5)

boxplot_funcional_A375 + scale_fill_manual(values = c("#440154", "#21908c", "#fde725", "#e7298a"))


# # stat_pvalue_manual(pwc, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE, y.position = c(19, 22, 25)) +
# labs(subtitle = get_test_label(anova.test, detailed = TRUE))+

boxplot_bc_filtered <- bcSubset(boxplot_bc, nan.cells = 0.95)
boxplot_bc_filtered@normalized[is.na(boxplot_bc_filtered@normalized)] <- 0                                                                                                                                                                       
boxplot_bc_recomputed <- bcRecompute(boxplot_bc_filtered, slot = "normalized") 

boxplot_bc_recomputed <- bcUMAP(boxplot_bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
boxplot_bc_recomputed <- bcRegressOut(boxplot_bc_recomputed, vars.to.regress = "nFeature_RNA", k.neighbors = 40)

bcSignatures(boxplot_bc_recomputed, UMAP = "beyondcell", 
             signatures = list(values = "BRAF_resistance"), pt.size = 1.5)
