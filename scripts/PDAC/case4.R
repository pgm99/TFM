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
library("GSEABase")
library("ggpubr")
library("SingleR")
set.seed(123)

case4ZY_mtx <- ReadMtx(mtx = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case4/Case4_ZY_matrix.mtx", 
                       cells = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case4/Case4_ZY_barcodes.tsv", 
                       features = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case4/Case4_ZY_features.tsv"
)
case4_meta <- CreateSeuratObject(counts = case4ZY_mtx)
case4_meta <- SetIdent(case4_meta, value = "case4_meta")

case4_meta[["percent.mt"]] <- PercentageFeatureSet(case4_meta, pattern = "^MT-")
case4_meta[["percent.rb"]] <- PercentageFeatureSet(case4_meta, pattern = "^RP[SL]")
head(case4_meta@meta.data, 5)
VlnPlot(case4_meta, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(case4_meta, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(case4_meta, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
case4_meta <- subset(case4_meta, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 10)


case4_meta <- NormalizeData(case4_meta, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
case4_meta <- FindVariableFeatures(case4_meta, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(case4_meta), 10)
plot1 <- VariableFeaturePlot(case4_meta)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(case4_meta)
case4_meta <- ScaleData(case4_meta, features = all.genes)



#Regresar cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
case4_meta <- CellCycleScoring(case4_meta, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
case4_meta <- RunPCA(case4_meta, features = c(s.genes, g2m.genes))
Idents(object = case4_meta) <- "Phase"
DimPlot(case4_meta, pt.size = 1.1)
case4_meta <- ScaleData(case4_meta, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(case4_meta))




case4_meta <- RunPCA(case4_meta, features = VariableFeatures(object = case4_meta))
ElbowPlot(case4_meta)

case4_meta <- FindNeighbors(case4_meta, dims = 1:15)
case4_meta <- FindClusters(case4_meta, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
clustree(case4_meta)
Idents(object = case4_meta) <- "RNA_snn_res.0.4"
case4_meta <- RunUMAP(case4_meta, dims = 1:15, min.dist = 0.5, n.neighbors = 30)
seurat_case4_meta_umap <- DimPlot(case4_meta, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)


saveRDS(case4_meta, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case4/case4_meta.rds")
case4_meta <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case1/case4_meta.rds")


###Cell type###
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
library(scRNAseq)
library(scuttle)
sceM <- MuraroPancreasData()
sceM <- sceM[,!is.na(sceM$label)]
sceM <- logNormCounts(sceM)
rownames(sceM) <- sub("[_].*", "", rownames(sceM))


Idents(object = case4_meta) <- "RNA_snn_res.0.4"
case4_meta_celltype <- SingleR(GetAssayData(case4_meta, assay = "RNA", slot = "data"),
                                clusters = Idents(case4_meta), ref = hpca.se,  labels =hpca.se$label.main, de.method = "wilcox") 

case4_meta[["cell_type"]] <- case4_meta_celltype$labels[match(case4_meta[[]][["RNA_snn_res.0.4"]], rownames(case4_meta_celltype))]
Idents(object = case4_meta) <- "cell_type"
case4_meta <- RunUMAP(case4_meta, dims = 1:15, min.dist = 0.5, n.neighbors = 30)
cell_type_case4_meta_umap <- DimPlot(case4_meta, reduction = "umap", label = TRUE, label.size = 4, pt.size = 1.2)
annotated_celltype <- case4_meta$cell_type
write.table(x = annotated_celltype, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case4/annotated_celltype.tsv", sep = "\t")






















