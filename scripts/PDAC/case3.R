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

###YF/TUMOR###


case3YF_mtx <- ReadMtx(mtx = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case3/Case3_YF/Case3_YF_matrix.mtx", 
                       cells = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case3/Case3_YF/Case3_YF_barcodes.tsv", 
                       features = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case3/Case3_YF/Case3_YF_features.tsv"
)
case3_tumor <- CreateSeuratObject(counts = case3YF_mtx)
case3_tumor <- SetIdent(case3_tumor, value = "case3_tumor")

case3_tumor[["percent.mt"]] <- PercentageFeatureSet(case3_tumor, pattern = "^MT-")
case3_tumor[["percent.rb"]] <- PercentageFeatureSet(case3_tumor, pattern = "^RP[SL]")
head(case3_tumor@meta.data, 5)
VlnPlot(case3_tumor, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(case3_tumor, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(case3_tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
case3_tumor <- subset(case3_tumor, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 10)


case3_tumor <- NormalizeData(case3_tumor, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
case3_tumor <- FindVariableFeatures(case3_tumor, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(case3_tumor), 10)
plot1 <- VariableFeaturePlot(case3_tumor)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(case3_tumor)
case3_tumor <- ScaleData(case3_tumor, features = all.genes)



#Regresar cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
case3_tumor <- CellCycleScoring(case3_tumor, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
case3_tumor <- RunPCA(case3_tumor, features = c(s.genes, g2m.genes))
Idents(object = case3_tumor) <- "Phase"
DimPlot(case3_tumor, pt.size = 1.1)
case3_tumor <- ScaleData(case3_tumor, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(case3_tumor))




case3_tumor <- RunPCA(case3_tumor, features = VariableFeatures(object = case3_tumor))
ElbowPlot(case3_tumor)

case3_tumor <- FindNeighbors(case3_tumor, dims = 1:15)
case3_tumor <- FindClusters(case3_tumor, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
clustree(case3_tumor)
Idents(object = case3_tumor) <- "RNA_snn_res.0.2"
case3_tumor <- RunUMAP(case3_tumor, dims = 1:15, min.dist = 0.5, n.neighbors = 30)
seurat_case3_tumor_umap <- DimPlot(case3_tumor, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)


saveRDS(case3_tumor, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case3/case3_tumor.rds")
case3_tumor <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case3/case3_tumor.rds")


###Cell type###
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
library(scRNAseq)
library(scuttle)
sceM <- MuraroPancreasData()
sceM <- sceM[,!is.na(sceM$label)]
sceM <- logNormCounts(sceM)
rownames(sceM) <- sub("[_].*", "", rownames(sceM))


Idents(object = case3_tumor) <- "RNA_snn_res.0.2"
case3_tumor_celltype <- SingleR(GetAssayData(case3_tumor, assay = "RNA", slot = "data"),
                            clusters = Idents(case3_tumor), ref = hpca.se,  labels =hpca.se$label.main, de.method = "wilcox") 

case3_tumor[["cell_type"]] <- case3_tumor_celltype$labels[match(case3_tumor[[]][["RNA_snn_res.0.2"]], rownames(case3_tumor_celltype))]
Idents(object = case3_tumor) <- "cell_type"
case3_tumor <- RunUMAP(case3_tumor, dims = 1:15, min.dist = 0.5, n.neighbors = 30)
cell_type_case3_tumor_umap <- DimPlot(case3_tumor, reduction = "umap", label = TRUE, label.size = 4, pt.size = 1.2)
annotated_celltype <- case3_tumor$cell_type
write.table(x = annotated_celltype, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case3/Case3_tumor/annotated_celltype_case3_tumor.tsv", sep = "\t")





#########ZY/META#########

case3ZY_mtx <- ReadMtx(mtx = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case3/Case3_ZY/Case3_ZY_matrix.mtx", 
                       cells = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case3/Case3_ZY/Case3_ZY_barcodes.tsv", 
                       features = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case3/Case3_ZY/Case3_ZY_features.tsv"
)
case3_meta <- CreateSeuratObject(counts = case3ZY_mtx)
case3_meta <- SetIdent(case3_meta, value = "case3_meta")

case3_meta[["percent.mt"]] <- PercentageFeatureSet(case3_meta, pattern = "^MT-")
case3_meta[["percent.rb"]] <- PercentageFeatureSet(case3_meta, pattern = "^RP[SL]")
head(case3_meta@meta.data, 5)
VlnPlot(case3_meta, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(case3_meta, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(case3_meta, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
case3_meta <- subset(case3_meta, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 10)


case3_meta <- NormalizeData(case3_meta, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
case3_meta <- FindVariableFeatures(case3_meta, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(case3_meta), 10)
plot1 <- VariableFeaturePlot(case3_meta)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(case3_meta)
case3_meta <- ScaleData(case3_meta, features = all.genes)



#Regresar cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
case3_meta <- CellCycleScoring(case3_meta, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
case3_meta <- RunPCA(case3_meta, features = c(s.genes, g2m.genes))
Idents(object = case3_meta) <- "Phase"
DimPlot(case3_meta, pt.size = 1.1)
case3_meta <- ScaleData(case3_meta, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(case3_meta))




case3_meta <- RunPCA(case3_meta, features = VariableFeatures(object = case3_meta))
ElbowPlot(case3_meta)

case3_meta <- FindNeighbors(case3_meta, dims = 1:15)
case3_meta <- FindClusters(case3_meta, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
clustree(case3_meta)
Idents(object = case3_meta) <- "RNA_snn_res.0.4"
case3_meta <- RunUMAP(case3_meta, dims = 1:15, min.dist = 0.5, n.neighbors = 30)
seurat_case3_meta_umap <- DimPlot(case3_meta, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)


saveRDS(case3_meta, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case3/case3_meta.rds")
case3_meta <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case3/case3_meta.rds")


###Cell type###
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
library(scRNAseq)
library(scuttle)
sceM <- MuraroPancreasData()
sceM <- sceM[,!is.na(sceM$label)]
sceM <- logNormCounts(sceM)
rownames(sceM) <- sub("[_].*", "", rownames(sceM))


Idents(object = case3_meta) <- "RNA_snn_res.0.4"
case3_meta_celltype <- SingleR(GetAssayData(case3_meta, assay = "RNA", slot = "data"),
                                clusters = Idents(case3_meta), ref = hpca.se,  labels =hpca.se$label.main, de.method = "wilcox") 

case3_meta[["cell_type"]] <- case3_meta_celltype$labels[match(case3_meta[[]][["RNA_snn_res.0.4"]], rownames(case3_meta_celltype))]
Idents(object = case3_meta) <- "cell_type"
case3_meta <- RunUMAP(case3_meta, dims = 1:15, min.dist = 0.5, n.neighbors = 30)
cell_type_case3_meta_umap <- DimPlot(case3_meta, reduction = "umap", label = TRUE, label.size = 4, pt.size = 1.2)
annotated_celltype <- case3_meta$cell_type
write.table(x = annotated_celltype, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case3/Case3_meta/annotated_celltype_case3_meta.tsv", sep = "\t")




########MERGED########
case3_tumor$sample = "tumor"
case3_meta$sample = "meta"
case3_merged <- merge(case3_tumor, y = case3_meta, project = "case3_PDAC")

saveRDS(case3_merged, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case3/case3_merged.rds")
case3_merged <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case3/case3_merged.rds")


Idents(object = case3_merged) <- "sample"

case3_merged[["percent.mt"]] <- PercentageFeatureSet(case3_merged, pattern = "^MT-")
case3_merged[["percent.rb"]] <- PercentageFeatureSet(case3_merged, pattern = "^RP[SL]")
head(case3_merged@meta.data, 5)
VlnPlot(case3_merged, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =2)

plot1 <- FeatureScatter(case3_merged, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 2)
plot2 <- FeatureScatter(case3_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 2)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
case3_merged <- subset(case3_merged, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 10)


case3_merged <- NormalizeData(case3_merged, normalization.method = "LogNormalize", scale.factor = 10000)

case3_merged <- FindVariableFeatures(case3_merged, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(case3_merged), 10)
plot1 <- VariableFeaturePlot(case3_merged)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(case3_merged)
case3_merged <- ScaleData(case3_merged, features = all.genes)

#Regresar cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
case3_merged<- CellCycleScoring(case3_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
case3_merged <- RunPCA(case3_merged, features = c(s.genes, g2m.genes))
Idents(object = case3_merged) <- "Phase"
DimPlot(case3_merged, pt.size = 1.1)
case3_merged <- ScaleData(case3_merged, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(case3_merged))



case3_merged <- RunPCA(case3_merged, features = VariableFeatures(object = case3_merged))
ElbowPlot(case3_merged)

case3_merged <- FindNeighbors(case3_merged, dims = 1:15)
case3_merged <- FindClusters(case3_merged, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
clustree(case3_merged)
Idents(object = case3_merged) <- "RNA_snn_res.0.2"


case3_merged <- RunUMAP(case3_merged, dims = 1:15, min.dist = 0.5, n.neighbors = 40)
seurat_case3_merged_umap30 <- DimPlot(case3_merged, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)

Idents(object = case3_merged) <- "sample"
sample_case3_merged_umap <- DimPlot(case3_merged, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)

ccase3_markers <- FindAllMarkers(case3_merged, only.pos = TRUE, logfc.threshold = 1, min.pct = 0.15)
DoHeatmap(case3_merged, features = ccase3_markers$gene) + guides(color="none")


################BEYONDCELL#################
SSc <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"), include.pathways = FALSE)

bc_case3 <- bcScore(case3_merged, SSc, expr.thres = 0.1)

saveRDS(bc_case3, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case3/bc_case3.rds")
bc_case3 <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case3/bc_case3.rds")

bc_filtered <- bcSubset(bc_case3, nan.cells = 0.95)
bc_filtered@normalized[is.na(bc_filtered@normalized)] <- 0                                                                                                                                                                       
bc_recomputed <- bcRecompute(bc_filtered, slot = "normalized") 


bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = NULL, k.neighbors = 40)
bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)

saveRDS(bc_recomputed, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case3/bc_case3_recomputed.rds")
bc_recomputed <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case3/bc_case3_recomputed.rds")

bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)
bc_recomputed <- bcRegressOut(bc_recomputed, vars.to.regress = "nFeature_RNA", k.neighbors = 40)

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)

clustree(bc_recomputed@meta.data, prefix ="bc_clusters_res.")
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.4", label.size = 5, pt.size = 1.5)
case3_merged <- AddMetaData(object = case3_merged, metadata = bc_recomputed@meta.data)

bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "sample", pt.size = 1.5, label = TRUE)
DimPlot(object = case3_merged)


#Drugs
dabrafenib_IDs <- FindDrugs(bc_recomputed, "dabrafenib")$IDs
trametinib_IDs <- FindDrugs(bc_recomputed, "trametinib")$IDs
cobimetinib_IDs <- FindDrugs(bc_recomputed, "cobimetinib")$IDs
gemcitabine_IDs <- FindDrugs(bc_recomputed, "gemcitabine")$IDs



dabrafenib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                           signatures = list(values = dabrafenib_IDs), pt.size = 1.5)


trametinib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                           signatures = list(values = trametinib_IDs), pt.size = 2)

cobimetinib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = cobimetinib_IDs), pt.size = 1.5)

gemcitabine <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = gemcitabine_IDs), pt.size = 1.8)




wrap_plots(dabrafenib, ncol = 3)
wrap_plots(trametinib, ncol = 3)
wrap_plots(cobimetinib, ncol = 1)
wrap_plots(gemcitabine, ncol = 3)



#Boxplot
gemcitabine_bc <- CreatebcObject(bc = bc_recomputed)
gemcitabine_bc <- bcSubset(bc = gemcitabine_bc, signatures = c("sig-20902"))
gemci_frame <- data.frame(gemcitabine_bc@normalized)
rownames(gemci_frame) <- c("gemcitabine_GDSC")
gemci_frame <- data.frame(t(gemci_frame))
rownames(gemci_frame) <- gsub("tumor_","", rownames(gemci_frame))
rownames(gemci_frame) <- gsub("meta_","", rownames(gemci_frame))

model <- data.frame(case3_merged$sample)
gemci_frame_type <- merge(x = gemci_frame, y = model, by= 0)
gemci_frame_long <- pivot_longer(data = gemci_frame_type, cols = c("gemcitabine_GDSC"), names_to = "sig", values_to = "enrichment_score")
colnames(gemci_frame_long) <- c("Row.names", "model", "sig", "enrichment_score")


sample <- factor(gemci_frame_long$model, levels=c("A375S", "A375RV", "A375RVC", "A375RVT"))
ggplot(gemci_frame_long) +
  geom_boxplot((aes(x=sample, y=enrichment_score, fill=sample)), show.legend = FALSE, outlier.shape = NA)+ 
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, hide.ns = FALSE, y.position = c(8,9,10,11,12,13)) +
  ggtitle("Comparison of gemcitabine signature between cell lines models")

stat.test <- gemci_frame_long %>% t_test((enrichment_score~model))
stat.test <- stat.test %>% add_xy_position(x = "model")








bc_recomputed <- bcRanks(bc_recomputed, idents = "sample")
#Obtain unextended therapeutic cluster-based statistics.
bc_recomputed <- bcRanks(bc_recomputed, idents = "sample", extended = FALSE)


bc4Squares(bc_recomputed, idents = "sample", lvl = NULL, top = 5)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.2", lvl = NULL, top = 3)









