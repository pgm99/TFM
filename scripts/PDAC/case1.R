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
library("SingleR")
library("rstatix")
set.seed(123)

case1YF_mtx <- ReadMtx(mtx = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case1/Case1_YF/Case1_YF_matrix.mtx", skip.cell = 1, 
                       cells = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case1/Case1_YF/Case1_YF_barcodes.tsv", 
                       features = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case1/Case1_YF/Case1_YF_features.tsv"
)
case1YF <- CreateSeuratObject(counts = case1YF_mtx)
case1YF <- SetIdent(case1YF, value = "case1FY")
####QC####
case1YF[["percent.mt"]] <- PercentageFeatureSet(case1YF, pattern = "^MT-")
case1YF[["percent.rb"]] <- PercentageFeatureSet(case1YF, pattern = "^RP[SL]")
head(case1YF@meta.data, 5)
VlnPlot(case1YF, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(case1YF, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(case1YF, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
case1YF <- subset(case1YF, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 10)


case1YF <- NormalizeData(case1YF, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
case1YF <- FindVariableFeatures(case1YF, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(case1YF), 10)
plot1 <- VariableFeaturePlot(case1YF)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(case1YF)
case1YF <- ScaleData(case1YF, features = all.genes)

#Regresar cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
case1YF<- CellCycleScoring(case1YF, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
case1YF <- RunPCA(case1YF, features = c(s.genes, g2m.genes))
Idents(object = case1YF) <- "Phase"
DimPlot(case1YF, pt.size = 1.1)
case1YF <- ScaleData(case1YF, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(case1YF))



case1YF <- RunPCA(case1YF, features = VariableFeatures(object = case1YF))
ElbowPlot(case1YF)

case1YF <- FindNeighbors(case1YF, dims = 1:15)
case1YF <- FindClusters(case1YF, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
clustree(case1YF)
Idents(object = case1YF) <- "RNA_snn_res.0.2"
case1YF <- RunUMAP(case1YF, dims = 1:15, min.dist = 0.5, n.neighbors = 30)
seurat_Case1_YF_umap <- DimPlot(case1YF, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)


saveRDS(case1YF, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/case1_YF.rds")
case1YF <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/case1_YF.rds")


library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
library(scRNAseq)
library(scuttle)
sceM <- MuraroPancreasData()
sceM <- sceM[,!is.na(sceM$label)]
sceM <- logNormCounts(sceM)
rownames(sceM) <- sub("[_].*", "", rownames(sceM))


Idents(object = case1YF) <- "RNA_snn_res.0.2"
case1YF_celltype <- SingleR(GetAssayData(case1YF, assay = "RNA", slot = "data"),
                                 clusters = Idents(case1YF), ref = hpca.se,  labels =hpca.se$label.main, de.method = "wilcox") 

case1YF[["cell_type"]] <- case1YF_celltype$labels[match(case1YF[[]][["RNA_snn_res.0.2"]], rownames(case1YF_celltype))]
Idents(object = case1YF) <- "cell_type"
case1YF <- RunUMAP(case1YF, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
cell_type_case1_YF_umap <- DimPlot(case1YF, reduction = "umap", label = TRUE, label.size = 4, pt.size = 1.2)
annotated_celltype <- case1YF$cell_type
write.table(x = annotated_celltype, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case1/Case1_YF/annotated_celltype.tsv", sep = "\t")















###############ZY################

case1ZY_mtx <- ReadMtx(mtx = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case1/Case1_ZY/Case1_ZY_matrix.mtx",
                   cells = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case1/Case1_ZY/Case1_ZY_barcodes.tsv", 
                   features = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case1/Case1_ZY/Case1_ZY_features.tsv"
                     )
case1ZY <- CreateSeuratObject(counts = case1ZY_mtx)
case1ZY <- SetIdent(case1ZY, value = "case1ZY")
####QC####
case1ZY[["percent.mt"]] <- PercentageFeatureSet(case1ZY, pattern = "^MT-")
case1ZY[["percent.rb"]] <- PercentageFeatureSet(case1ZY, pattern = "^RP[SL]")
head(case1ZY@meta.data, 5)
VlnPlot(case1ZY, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(case1ZY, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(case1ZY, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
case1ZY <- subset(case1ZY, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 10)

case1ZY <- NormalizeData(case1ZY, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
case1ZY <- FindVariableFeatures(case1ZY, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(case1ZY), 10)
plot1 <- VariableFeaturePlot(case1ZY)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(case1ZY)
case1ZY <- ScaleData(case1ZY, features = all.genes)

#Regresar cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
case1ZY<- CellCycleScoring(case1ZY, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
case1ZY <- RunPCA(case1ZY, features = c(s.genes, g2m.genes))
Idents(object = case1ZY) <- "Phase"
DimPlot(case1ZY, pt.size = 1.1)
case1ZY <- ScaleData(case1ZY, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(case1ZY))



case1ZY <- RunPCA(case1ZY, features = VariableFeatures(object = case1ZY))
ElbowPlot(case1ZY)

case1ZY <- FindNeighbors(case1ZY, dims = 1:15)
case1ZY <- FindClusters(case1ZY, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
clustree(case1ZY)
Idents(object = case1ZY) <- "RNA_snn_res.0.2"
case1ZY <- RunUMAP(case1ZY, dims = 1:15, min.dist = 0.5, n.neighbors = 30)
seurat_Case1_ZY_umap <- DimPlot(case1ZY, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)

saveRDS(case1ZY, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/case1_ZY.rds")
case1ZY <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/case1_ZY.rds")

library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()

Idents(object = case1ZY) <- "RNA_snn_res.0.2"
case1ZY_celltype <- SingleR(GetAssayData(case1ZY, assay = "RNA", slot = "data"),
                            clusters = Idents(case1ZY), ref = hpca.se,  labels =hpca.se$label.main, de.method = "wilcox") 

case1ZY[["cell_type"]] <- case1ZY_celltype$labels[match(case1ZY[[]][["RNA_snn_res.0.2"]], rownames(case1ZY_celltype))]
Idents(object = case1ZY) <- "cell_type"
case1ZY <- RunUMAP(case1ZY, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
cell_type_case1_ZY_umap <- DimPlot(case1ZY, reduction = "umap", label = TRUE, label.size = 4, pt.size = 1.2)
annotated_celltype <- case1ZY$cell_type
write.table(x = annotated_celltype, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case1/Case1_ZY/annotated_celltype_case1_ZY.tsv", sep = "\t")

###MERGED###
case1YF$sample = "YF"
case1ZY$sample = "ZY"
case1_merged <- merge(case1YF, y = case1ZY, project = "case1_PDAC")

saveRDS(case1_merged, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/case1_merged.rds")
case1_merged <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/case1_merged.rds")


Idents(object = case1_merged) <- "sample"

case1_merged[["percent.mt"]] <- PercentageFeatureSet(case1_merged, pattern = "^MT-")
case1_merged[["percent.rb"]] <- PercentageFeatureSet(case1_merged, pattern = "^RP[SL]")
head(case1_merged@meta.data, 5)
VlnPlot(case1_merged, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(case1_merged, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(case1_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
case1_merged <- subset(case1_merged, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 10)


case1_merged <- NormalizeData(case1_merged, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
case1_merged <- FindVariableFeatures(case1_merged, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(case1_merged), 10)
plot1 <- VariableFeaturePlot(case1_merged)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(case1_merged)
case1_merged <- ScaleData(case1_merged, features = all.genes)

#Regresar cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
case1_merged<- CellCycleScoring(case1_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
case1_merged <- RunPCA(case1_merged, features = c(s.genes, g2m.genes))
Idents(object = case1_merged) <- "Phase"
DimPlot(case1_merged, pt.size = 1.1)
case1_merged <- ScaleData(case1_merged, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(case1_merged))



case1_merged <- RunPCA(case1_merged, features = VariableFeatures(object = case1_merged))
ElbowPlot(case1_merged)

case1_merged <- FindNeighbors(case1_merged, dims = 1:15)
case1_merged <- FindClusters(case1_merged, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
clustree(case1_merged)
Idents(object = case1_merged) <- "RNA_snn_res.0.2"


case1_merged <- RunUMAP(case1_merged, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
seurat_case1_merged_umap30 <- DimPlot(case1_merged, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)

Idents(object = case1_merged) <- "sample"
sample_case1_merged_umap <- DimPlot(case1_merged, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)


case1_markers <- FindAllMarkers(case1_merged, only.pos = TRUE, logfc.threshold = 1, min.pct = 0.15)
case1_YF_markers <- FindMarkers(case1_merged, ident.1 = "YF", only.pos = TRUE, logfc.threshold = 1, min.pct = 0.15)
DoHeatmap(case1_merged, features = case1_markers$gene) + guides(color="none")
DoHeatmap(case1_merged, features = rownames(case1_YF_markers))



###CELL TYPE###
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
library(scRNAseq)
library(scuttle)
sceM <- MuraroPancreasData()
sceM <- sceM[,!is.na(sceM$label)]
sceM <- logNormCounts(sceM)
rownames(sceM) <- sub("[_].*", "", rownames(sceM))



Idents(object = case1_merged) <- "RNA_snn_res.0.2"
case1_merged_celltype <- SingleR(GetAssayData(case1_merged, assay = "RNA", slot = "data"),
                                 clusters = Idents(case1_merged), ref = hpca.se,  labels =hpca.se$label.main, de.method = "wilcox") 

case1_merged[["cell_type"]] <- case1_merged_celltype$labels[match(case1_merged[[]][["RNA_snn_res.0.2"]], rownames(case1_merged_celltype))]
Idents(object = case1_merged) <- "cell_type"
case1_merged <- RunUMAP(case1_merged, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
cell_type_case1_merged_umap <- DimPlot(case1_merged, reduction = "umap", label = TRUE, label.size = 4, pt.size = 1.2)

Idents(object = case1_merged) <- "RNA_snn_res.0.2"
case1_merged_celltype2 <- SingleR(GetAssayData(case1_merged, assay = "RNA", slot = "data"),
                                  clusters = Idents(case1_merged), ref = sceM,  labels =sceM$label) 
case1_merged[["cell_type2"]] <- case1_merged_celltype2$labels[match(case1_merged[[]][["RNA_snn_res.0.2"]], rownames(case1_merged_celltype2))]
Idents(object = case1_merged) <- "cell_type2"
case1_merged <- RunUMAP(case1_merged, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
cell_type2_case2_ZY_umap <- DimPlot(case1_merged, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2, )
 



################BEYONDCELL#################
SSc <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"), )

bc_case1 <- bcScore(case1_merged, SSc, expr.thres = 0.1)

saveRDS(bc_case1, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/bc_case1.rds")
bc_case1 <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/bc_case1.rds")

bc_filtered <- bcSubset(bc_case1, nan.cells = 0.95)
bc_filtered@normalized[is.na(bc_filtered@normalized)] <- 0                                                                                                                                                                       
bc_recomputed <- bcRecompute(bc_filtered, slot = "normalized") 


bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = NULL, k.neighbors = 40)
bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)

saveRDS(bc_recomputed, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/bc_case1_recomputed.rds")
bc_recomputed <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/bc_case1_recomputed.rds")

bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)
bc_recomputed <- bcRegressOut(bc_recomputed, vars.to.regress = "nFeature_RNA", k.neighbors = 40)

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)

clustree(bc_recomputed@meta.data, prefix ="bc_clusters_res.")
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.4", label.size = 5, pt.size = 1.5)
case1_merged <- AddMetaData(object = case1_merged, metadata = bc_recomputed@meta.data)

bcClusters(bc_recomputed, UMAP = "Seurat", idents = "RNA_snn_res.0.4", pt.size = 1.5, label = TRUE)
DimPlot(object = case1_merged)

DimPlot(
  case1_merged,
  reduction = "beyondcell",
  label = TRUE,
  label.size = 5,
  group.by = "bc_clusters_res.0.4",
  shape.by = "type",
  raster = FALSE,
  pt.size = 0.8
)

bc_samples <- bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "sample", pt.size = 1.5)






dabrafenib_info <- FindDrugs(bc_recomputed, "dabrafenib")

dabrafenib_IDs <- FindDrugs(bc_recomputed, "dabrafenib")$IDs
trametinib_IDs <- FindDrugs(bc_recomputed, "trametinib")$IDs
cobimetinib_IDs <- FindDrugs(bc_recomputed, "cobimetinib")$IDs
gemcitabine_IDs <- FindDrugs(bc_recomputed, "gemcitabine")$IDs


dabrafenib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                           signatures = list(values = "sig-20909"), pt.size = 2)


trametinib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                           signatures = list(values = "sig-21071"), pt.size = 2)

cobimetinib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = cobimetinib_IDs), pt.size = 2)

gemcitabine <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = "sig-21081"), pt.size = 2)




wrap_plots(dabrafenib, ncol = 1)
wrap_plots(cobitinib, ncol = 1)
wrap_plots(cobimetinib, ncol = 1)
wrap_plots(gemcitabine, ncol = 3)

# Obtain condition-based statistics.
bc_recomputed <- bcRanks(bc_recomputed, idents = "sample")
#Obtain unextended therapeutic cluster-based statistics.
bc_recomputed <- bcRanks(bc_recomputed, idents = "sample", extended = FALSE)


bc4Squares(bc_recomputed, idents = "sample", lvl = NULL, top = 5)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.2", lvl = NULL, top = 3)


YF1 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20948"), pt.size = 1.5)
YF2 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20908"), pt.size = 1.5)
YF3 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20925"), pt.size = 1.5)
YF4 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20931"), pt.size = 1.5)
YF5 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20917"), pt.size = 1.5)

top_diff_YF <- bc_samples+YF1+YF2+YF3+YF4+YF5


YF6 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21011"), pt.size = 1.5)
YF7 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21053"), pt.size = 1.5)
YF8 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21322"), pt.size = 1.5)
YF9 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21119"), pt.size = 1.5)
YF10 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                     signatures = list(values = "sig-21055"), pt.size = 1.5)

low_diff_YF <- bc_samples+YF6+YF7+YF8+YF9+YF10






ZY1 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
             signatures = list(values = "sig-20948"), pt.size = 1.5)
ZY2 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
             signatures = list(values = "sig-20908"), pt.size = 1.5)
ZY3 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
             signatures = list(values = "sig-20925"), pt.size = 1.5)
ZY4 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
             signatures = list(values = "sig-20931"), pt.size = 1.5)
ZY5 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
             signatures = list(values = "sig-20917"), pt.size = 1.5)

top_diff_ZY <- bc_samples+ZY1+ZY2+ZY3+ZY4+ZY5


ZY6 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21011"), pt.size = 1.5)
ZY7 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21053"), pt.size = 1.5)
ZY8 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21322"), pt.size = 1.5)
ZY9 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21119"), pt.size = 1.5)
ZY10 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21055"), pt.size = 1.5)

low_diff_ZY <- bc_samples+ZY6+ZY7+ZY8+ZY9+ZY10


bcHistogram(bc_recomputed, signatures = "sig-21081", idents = "bc_clusters_res.0.4")

#Boxplots
trametinib_bc <- CreatebcObject(bc = bc_recomputed)
trametinib_bc <- bcSubset(bc = trametinib_bc, signatures = c("sig-21037"))
trame_frame <- data.frame(trametinib_bc@normalized)
rownames(trame_frame) <- c("cobimetinib_PRISM")
trame_frame <- data.frame(t(trame_frame))
rownames(trame_frame) <- gsub("\\.", "-", rownames(trame_frame))


sample <- data.frame(case1_merged$sample)
sample$case1_merged.sample <- gsub("YF", "tumor", sample$case1_merged.sample)
sample$case1_merged.sample <- gsub("ZY", "meta", sample$case1_merged.sample)
trame_frame_type <- merge(x = trame_frame, y = sample, by= 0)
trame_frame_long <- pivot_longer(data = trame_frame_type, cols = c("cobimetinib_PRISM"), names_to = "sig", values_to = "enrichment_score")
colnames(trame_frame_long) <- c("Cells", "sample", "sig", "enrichment_score")


trame_frame_long$sample <- factor(trame_frame_long$sample, levels=c("tumor", "meta"), ordered = TRUE)
ggplot(trame_frame_long) +
  geom_boxplot((aes(x=sample, y=enrichment_score, fill=sample)), show.legend = FALSE, outlier.shape = NA)+ 
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, hide.ns = FALSE, size = 8) + scale_fill_brewer(palette="Dark2")+
  ggtitle("Comparison of cobimetinib signature between samples") + theme(plot.title = element_text(size = 30))


stat.test <- trame_frame_long %>% t_test((enrichment_score~sample))
stat.test <- stat.test %>% add_xy_position(x = "sample")

# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.6141




####BARPLOT CON %########
data2 <- data %>%
  group_by(TCs, sample) %>%
  summarise(count = n(), .groups = "keep") %>%
  group_by(TCs) %>%
  mutate(percentage = (count / sum(count)) * 100)

# Create the barplot
ggplot(data2, aes(x = as.factor(TCs), y = percentage, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Percentage (%)", x = "Values in TCs") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = c("YF" = "blue", "ZY" = "red")) +
  theme_minimal()
















####Funcional####

# hallmarks <- read.table(file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/data/genesets/h.all.v2023.1.Hs.symbols.gmt", header = F, sep = "\t", fill = T)
# hallmarks$V1=paste(hallmarks$V1,"_UP", sep = "")
# write.table(hallmarks, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/data/genesets/h.all.v2023.1.Hs.symbols.gmt", sep = "\t", row.names = F, col.names = F, quote = F)





gs_functional <- GenerateGenesets(x = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/data/genesets/h.all.v2023.1.Hs.symbols.gmt")
bc_functional <- bcScore(case1_merged, gs_functional, expr.thres = 0.1)
bc_filtered_functional <- bcSubset(bc_functional, nan.cells = 0.95)
bc_filtered_functional@normalized[is.na(bc_filtered_functional@normalized)] <- 0                                                                                                                                                                       
bc_recomputed_functional <- bcRecompute(bc_filtered_functional, slot = "normalized") 

bc_recomputed_functional <- bcUMAP(bc_recomputed_functional, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40, add.DSS = FALSE)
bc_recomputed_functional <- bcRegressOut(bc_recomputed_functional, vars.to.regress = "nFeature_RNA", k.neighbors = 40)
bc_recomputed_functional <- bcUMAP(bc_recomputed_functional, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)

clustree(bc_recomputed_functional@meta.data, prefix ="bc_clusters_res.")
bcClusters(bc_recomputed_functional, UMAP = "beyondcell", idents = "celltype", label.size = 5, pt.size = 1.3)
bc_sample_umap <- bcClusters(bc_recomputed_functional, UMAP = "beyondcell", idents = "sample", label.size = 5, pt.size = 1.2)

case1_merged <- AddMetaData(object = case1_merged, metadata = bc_recomputed@meta.data)



saveRDS(bc_recomputed_functional, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/case1_merged_funcional.rds")
bc_recomputed_functional <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/case1_merged_funcional.rds")



# data <- data.frame(bc_recomputed_functional@normalized)
# clusters <- data.frame(case1_merged$RNA_snn_res.0.2)
# sample <- data.frame(case1_merged$sample)
# clusters <- setNames(clusters, nm = "clusters")
# sample <- setNames(sample, nm= "sample")
# data <- t(data)
# data <- data[order(sample$sample), , drop= FALSE]
# data <- t(data)
# data <- as.matrix(data)
# data <- as.data.frame(data)
# 
# color = colorRampPalette(c("navy", "white", "firebrick3"))(200)
# pheatmap(data, scale = "row", 
#          show_colnames = F, cluster_rows = F, cluster_cols = F, annotation_col = type, treeheight_col = 0, 
#          legend = T,  annotation_names_col = T, annotation_legend = T, color = color)
# 
# DoHeatmap(data)


# plot1 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
#                       signatures = list(values = "GOBP_NEGATIVE_REGULATION_OF_MAPK_CASCADE"), pt.size = 0.8)
# 
# plot2 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
#                       signatures = list(values = "GOBP_POSITIVE_REGULATION_OF_MAPK_CASCADE"), pt.size = 0.8)

plot1 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_PEROXISOME"), pt.size = 2)

plot2 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_HYPOXIA"), pt.size = 2)

plot3 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_ALLOGRAFT_REJECTION"), pt.size = 2)

plot4 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_APOPTOSIS"), pt.size = 2)

plot5 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_DNA_REPAIR"), pt.size = 2)

plot6 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_FATTY_ACID_METABOLISM"), pt.size = 2)

plot7 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_MYC_TARGETS_V1"), pt.size = 2)

plot8 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"), pt.size = 2)

plot9 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_P53_PATHWAY"), pt.size = 2)

resumen_plots <- wrap_plots(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9)


plot12 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"), pt.size = 2)

plot13 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_TGF_BETA_SIGNALING"), pt.size = 2)

plot14 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_UV_RESPONSE"), pt.size = 2)

plot15 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_XENOBIOTIC_METABOLISM"), pt.size = 2)

plot11 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"), pt.size = 2)

plot16 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"), pt.size = 2)

plot17 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_MTORC1_SIGNALING"), pt.size = 2)

resumen_plots2 <-wrap_plots(plot11, plot12, plot13, plot15, plot16, plot17)



boxplot_bc <- CreatebcObject(bc = bc_recomputed_functional)
MAP2K1_bc <- bcSubset(bc = boxplot_bc, signatures = c("HALLMARK_HYPOXIA"))
MAP2K1_frame <- data.frame(MAP2K1_bc@normalized)
rownames(MAP2K1_frame) <- c("HALLMARK_HYPOXIA")
MAP2K1_frame <- data.frame(t(MAP2K1_frame))
sample <- data.frame(case1_merged$sample)
colnames(sample) <- "sample"
# TCs <- data.frame(case1_merged$bc_clusters_res.0.4)
# colnames(TCs) <- "TCs"
rownames(MAP2K1_frame) <- gsub("\\.", "-", rownames(MAP2K1_frame))
sample$sample <- gsub("YF", "tumor", sample$sample)
sample$sample <- gsub("ZY", "meta", sample$sample)
sample$sample <- gsub("ZC", "normal", sample$sample)

MAP2K1_frame_model <- merge(x = MAP2K1_frame, y = sample, by= 0)
# MAP2K1_frame_TCs <- merge(x = MAP2K1_frame_model, y = TCs$TCs, by= 0)
# MAP2K1_frame_TCs <- MAP2K1_frame_TCs[,-1]
colnames(MAP2K1_frame_model) <- c("Cells", "HALLMARK_HYPOXIA", "sample")
MAP2K1_frame_long <- pivot_longer(data = MAP2K1_frame_model, cols = c("HALLMARK_HYPOXIA"), names_to = "sig", values_to = "enrichment_score")

MAP2K1_frame_long$sample <- factor(MAP2K1_frame_long$sample, levels=c("tumor", "meta", "normal"), ordered = TRUE)
# TCs <- factor(MAP2K1_frame_long$TCs, levels=c("0", "1", "2", "3", "4", "5"))



stat.test <- MAP2K1_frame_long %>% t_test((enrichment_score~sample))
stat.test <- stat.test %>% add_xy_position(x = "sample")

ggplot(MAP2K1_frame_long) +
  geom_boxplot((aes(x=sample, y=enrichment_score, fill=sample)), show.legend = FALSE, outlier.shape = NA) + 
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, hide.ns = FALSE, y.position = c(20, 25, 30)) + scale_fill_brewer(palette="Dark2")+
  ggtitle("Comparison of hypoxia signature between samples") + theme(plot.title = element_text(size = 30))



ggplot(MAP2K1_frame_long) +
  geom_dotplot(aes(x = sample, y = enrichment_score, fill = sample), binaxis = "y", stackdir = "center", dotsize = 0.1) +
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, hide.ns = FALSE) +
  scale_fill_brewer(palette = "Dark2") +
  ggtitle("Comparison of hypoxia signature between samples") +
  theme(plot.title = element_text(size = 30))



?geom_dotplot







