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
library("GSEABase")
library("ggpubr")
library("SingleR")
library("biomaRt")
library("curl")
library("SingleR")
library("rstatix")
library("ComplexHeatmap")
set.seed(123)

case2YF_mtx <- ReadMtx(mtx = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case2/Case2_YF/Case2_YF_matrix.mtx", 
                       cells = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case2/Case2_YF/Case2_YF_barcodes.tsv", 
                       features = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case2/Case2_YF/Case2_YF_features.tsv"
)
case2YF <- CreateSeuratObject(counts = case2YF_mtx)
case2YF <- SetIdent(case2YF, value = "case2FY")
####QC####
case2YF[["percent.mt"]] <- PercentageFeatureSet(case2YF, pattern = "^MT-")
case2YF[["percent.rb"]] <- PercentageFeatureSet(case2YF, pattern = "^RP[SL]")
head(case2YF@meta.data, 5)
VlnPlot(case2YF, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(case2YF, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(case2YF, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
case2YF <- subset(case2YF, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 10)


case2YF <- NormalizeData(case2YF, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
case2YF <- FindVariableFeatures(case2YF, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(case2YF), 10)
plot1 <- VariableFeaturePlot(case2YF)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(case2YF)
case2YF <- ScaleData(case2YF, features = all.genes)

#Regresar cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
case2YF<- CellCycleScoring(case2YF, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
case2YF <- RunPCA(case2YF, features = c(s.genes, g2m.genes))
Idents(object = case2YF) <- "Phase"
DimPlot(case2YF, pt.size = 1.1)
case2YF <- ScaleData(case2YF, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(case2YF))



case2YF <- RunPCA(case2YF, features = VariableFeatures(object = case2YF))
ElbowPlot(case2YF)

case2YF <- FindNeighbors(case2YF, dims = 1:15)
case2YF <- FindClusters(case2YF, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
clustree(case2YF)
Idents(object = case2YF) <- "RNA_snn_res.0.2"
case2YF <- RunUMAP(case2YF, dims = 1:15, min.dist = 0.2, n.neighbors = 30)
seurat_case2_YF_umap <- DimPlot(case2YF, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)




saveRDS(case2YF, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_YF.rds")
case2YF <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_YF.rds")
case2YF <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_YF_bcmeta.rds")

Idents(object = case2YF) <- "subclone"

case2_YF_markers <- FindAllMarkers(case2YF, only.pos = TRUE, logfc.threshold = 1, min.pct = 0.15)
# case2_YF_markers$subclone <- case2YF@meta.data$subclone
# case2YF$subclone <- as.factor(case2YF$subclone)
# DoHeatmap(case2YF, features = case2_YF_markers$gene, label = FALSE, group.by = c("RNA_snn_res.0.2", "subclone"))

gene_expression_matrix <- case2YF@assays$RNA@scale.data
selected_genes <- case2_YF_markers$gene
gene_expression_matrix_selected <- gene_expression_matrix[selected_genes, ]
case2YF@meta.data$subclone <- as.factor(case2YF@meta.data$subclone)

df <- data.frame("subclone" = case2YF@meta.data$subclone)
df[is.na(df)] <- 0
df <- df[complete.cases(df),]


subclone_colors <- c("1" = "green", "2" = "purple", "3" = "orange", "4" = "blue", "5" = "yellow", "6" = "red")
clusters_colors <- c("0" = "brown","1" = "green", "2" = "purple", "3" = "orange", "4" = "blue", "5" = "yellow", "6" = "red", "7" = "pink", "8"="gold")
groups <- c("Cluster", "subclone")
heatmap_object <- Heatmap(
  gene_expression_matrix_selected,
  name = "Gene Expression",
  cluster_columns = FALSE,  # Puedes ajustar esto según tus preferencias
  cluster_rows = FALSE,     # Puedes ajustar esto según tus preferencias
  show_row_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_column_names = FALSE,
  top_annotation = HeatmapAnnotation(
    df = df,
    col = list("subclone" = subclone_colors)
  ),
  column_title = "Samples",
  column_title_side = "bottom",
  row_title = "Genes",
  row_title_side = "left"
)
ordered_columns <- order(case2YF@meta.data$subclone)
heatmap_object <- heatmap_object[, ordered_columns]
draw(heatmap_object)

#subclones###
column_to_extract <- "subclone"
matching_cells <- intersect(colnames(case2YF), colnames(pdac_shu_zhang[[4]]))
metadata_from_object2 <- pdac_shu_zhang[[4]]@meta.data[matching_cells, column_to_extract]
case2YF@meta.data[matching_cells, column_to_extract] <- metadata_from_object2

Idents(object = case2YF) <- "subclone"
DimPlot(case2YF, pt.size = 1.1)

DimPlot(
  case2YF,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  group.by = c("RNA_snn_res.0.2", "subclone")
  #  shape.by = "subclone",
)


###CELL TYPE###

library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
library(scRNAseq)
library(scuttle)
sceM <- MuraroPancreasData(location = TRUE)
sceM <- sceM[,!is.na(sceM$label)]
sceM <- logNormCounts(sceM)
rownames(sceM) <- sub("[_].*", "", rownames(sceM))


Idents(object = case2YF) <- "RNA_snn_res.0.2"
case2YF_celltype <- SingleR(GetAssayData(case2YF, assay = "RNA", slot = "data"),
                             clusters = Idents(case2YF), ref = hpca.se,  labels =hpca.se$label.main, de.method = "wilcox") 

case2YF[["cell_type"]] <- case2YF_celltype$labels[match(case2YF[[]][["RNA_snn_res.0.2"]], rownames(case2YF_celltype))]
Idents(object = case2YF) <- "cell_type"
case2YF <- RunUMAP(case2YF, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
cell_type_case2_YF_umap <- DimPlot(case2YF, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)

annotated_celltype <- case2YF$cell_type
write.table(x = annotated_celltype, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case2/Case2_YF/annotated_celltype.tsv", sep = "\t")




# Idents(object = case2YF) <- "RNA_snn_res.0.2"
# case2YF_celltype2 <- SingleR(GetAssayData(case2YF, assay = "RNA", slot = "data"),
#                              clusters = Idents(case2YF), ref = sceM,  labels =sceM$label, de.method = "wilcox") 
# case2YF[["cell_type2"]] <- case2YF_celltype2$labels[match(case2YF[[]][["RNA_snn_res.0.2"]], rownames(case2YF_celltype2))]
# Idents(object = case2YF) <- "cell_type2"
# case2YF <- RunUMAP(case2YF, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
# cell_type2_case2_YF_umap <- DimPlot(case2YF, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)



# new.cluster.ids <- c("T cell", "Fibroblasts", "", "Monocyte", "Mast cell", "", "", "", "B cell", "Ductal cell", "Plasma cell")


SSc <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"), )

bc_case2YF <- bcScore(case2YF, SSc, expr.thres = 0.1)
bc_filtered <- bcSubset(bc_case2YF, nan.cells = 0.95)
bc_filtered@normalized[is.na(bc_filtered@normalized)] <- 0                                                                                                                                                                       
bc_recomputed <- bcRecompute(bc_filtered, slot = "normalized") 

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)

bc_recomputed <- bcRegressOut(bc_recomputed, vars.to.regress = "nFeature_RNA", k.neighbors = 40)

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
gemcitabine <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = "sig-20902"), pt.size = 2.5)
gemcitabine+theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                  panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
                  plot.subtitle = element_text(size=15))

gs_functional <- GenerateGenesets(x = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/data/genesets/h.all.v2023.1.Hs.symbols.gmt")
bc_functional <- bcScore(case2YF, gs_functional, expr.thres = 0.1)
resistance_markers <- as.data.frame(case2_merged@meta.data$resistance_markers1)
resistance_markers <- as.data.frame(bc_functional@normalized["HALLMARK_XENOBIOTIC_METABOLISM", ])
resistance_markers$BCscore <- (bc_case2YF@normalized["sig-20902", ])
colnames(resistance_markers) <- c("xenobiotic_metabolism", "BCscore")

resistance_markers$xenobiotic_metabolism[resistance_markers$xenobiotic_metabolism < 0] <- NA

ggplot(resistance_markers, aes(x=xenobiotic_metabolism, y=BCscore)) + 
  geom_point()+ 
  stat_cor(method = "spearman", label.y = 60, size=7)+ geom_smooth(method = "glm", se = FALSE)+ 
  ggtitle("Correlation between xenobiotic metabolism and gemcitabine BCscore in tumor sample")+theme(plot.title = element_text(size = 20), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"), 
                                                                                               panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1.7), 
                                                                                               plot.subtitle = element_text(size=15))





###############ZY################

case2ZY_mtx <- ReadMtx(mtx = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case2/Case2_ZY/Case2_ZY_matrix.mtx",
                       cells = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case2/Case2_ZY/Case2_ZY_barcodes.tsv", 
                       features = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case2/Case2_ZY/Case2_ZY_features.tsv"
)
case2ZY <- CreateSeuratObject(counts = case2ZY_mtx)
case2ZY <- SetIdent(case2ZY, value = "case2ZY")
####QC####
case2ZY[["percent.mt"]] <- PercentageFeatureSet(case2ZY, pattern = "^MT-")
case2ZY[["percent.rb"]] <- PercentageFeatureSet(case2ZY, pattern = "^RP[SL]")
head(case2ZY@meta.data, 5)
VlnPlot(case2ZY, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(case2ZY, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(case2ZY, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
case2ZY <- subset(case2ZY, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 10)

case2ZY <- NormalizeData(case2ZY, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
case2ZY <- FindVariableFeatures(case2ZY, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(case2ZY), 10)
plot1 <- VariableFeaturePlot(case2ZY)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(case2ZY)
case2ZY <- ScaleData(case2ZY, features = all.genes)

#Regresar cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
case2ZY<- CellCycleScoring(case2ZY, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
case2ZY <- RunPCA(case2ZY, features = c(s.genes, g2m.genes))
Idents(object = case2ZY) <- "Phase"
DimPlot(case2ZY, pt.size = 1.1)
case2ZY <- ScaleData(case2ZY, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(case2ZY))


case2ZY <- RunPCA(case2ZY, features = VariableFeatures(object = case2ZY))
ElbowPlot(case2ZY)

case2ZY <- FindNeighbors(case2ZY, dims = 1:15)
case2ZY <- FindClusters(case2ZY, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
clustree(case2ZY)
Idents(object = case2ZY) <- "RNA_snn_res.0.2"
case2ZY <- RunUMAP(case2ZY, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
seurat_case2_ZY_umap <- DimPlot(case2ZY, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)

saveRDS(case2ZY, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_ZY.rds")
case2ZY <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_ZY.rds")

library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
library(scRNAseq)
library(scuttle)
sceM <- MuraroPancreasData()
sceM <- sceM[,!is.na(sceM$label)]
sceM <- logNormCounts(sceM)
rownames(sceM) <- sub("[_].*", "", rownames(sceM))



Idents(object = case2ZY) <- "RNA_snn_res.0.2"
case2ZY_celltype <- SingleR(GetAssayData(case2ZY, assay = "RNA", slot = "data"),
                            clusters = Idents(case2ZY), ref = hpca.se,  labels =hpca.se$label.main, de.method = "wilcox") 

case2ZY[["cell_type"]] <- case2ZY_celltype$labels[match(case2ZY[[]][["RNA_snn_res.0.2"]], rownames(case2ZY_celltype))]
Idents(object = case2ZY) <- "cell_type"
case2ZY <- RunUMAP(case2ZY, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
cell_type_case2_ZY_umap <- DimPlot(case2ZY, reduction = "umap", label = TRUE, label.size = 4, pt.size = 1.2)

annotated_celltype <- case2ZY$cell_type
write.table(x = annotated_celltype, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case2/Case2_ZY/annotated_celltype_case2ZY.tsv", sep = "\t")




# Idents(object = case2ZY) <- "RNA_snn_res.0.2"
# case2ZY_celltype2 <- SingleR(GetAssayData(case2ZY, assay = "RNA", slot = "data"),
#                              clusters = Idents(case2ZY), ref = sceM,  labels =sceM$label, de.method = "wilcox") 
# case2ZY[["cell_type2"]] <- case2ZY_celltype2$labels[match(case2ZY[[]][["RNA_snn_res.0.2"]], rownames(case2ZY_celltype2))]
# Idents(object = case2ZY) <- "cell_type2"
# case2ZY <- RunUMAP(case2ZY, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
# cell_type2_case2_ZY_umap <- DimPlot(case2ZY, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)

SSc <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"), )

bc_case2ZY <- bcScore(case2ZY, SSc, expr.thres = 0.1)
bc_filtered <- bcSubset(bc_case2ZY, nan.cells = 0.95)
bc_filtered@normalized[is.na(bc_filtered@normalized)] <- 0                                                                                                                                                                       
bc_recomputed <- bcRecompute(bc_filtered, slot = "normalized") 

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)

bc_recomputed <- bcRegressOut(bc_recomputed, vars.to.regress = "nFeature_RNA", k.neighbors = 40)

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
gemcitabine <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = "sig-20902"), pt.size = 2.5)
gemcitabine+theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                  panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
                  plot.subtitle = element_text(size=15))

gs_functional <- GenerateGenesets(x = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/data/genesets/h.all.v2023.1.Hs.symbols.gmt")
bc_functional_ZY <- bcScore(case2ZY, gs_functional, expr.thres = 0.1)
resistance_markers <- as.data.frame(bc_case2ZY@meta.data$resistance_markers1)
resistance_markers <- as.data.frame(bc_functional_ZY@normalized["HALLMARK_XENOBIOTIC_METABOLISM", ])
resistance_markers$BCscore <- (bc_case2ZY@normalized["sig-20902", ])
colnames(resistance_markers) <- c("xenobiotic_metabolism", "BCscore")

resistance_markers$xenobiotic_metabolism[resistance_markers$xenobiotic_metabolism < 0] <- NA

ggplot(resistance_markers, aes(x=xenobiotic_metabolism, y=BCscore)) + 
  geom_point()+ 
  stat_cor(method = "spearman", label.y = 60, size=7)+ geom_smooth(method = "glm", se = FALSE)+ 
  ggtitle("Correlation between xenobiotic metabolism and gemcitabine BCscore in metastasis sample")+theme(plot.title = element_text(size = 20), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"), 
                                                                                                     panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1.7), 
                                                                                                     plot.subtitle = element_text(size=15))











###############ZC################

case2ZC_mtx <- ReadMtx(mtx = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case2/Case2_ZC/Case2_ZC_matrix.mtx",
                       cells = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case2/Case2_ZC/Case2_ZC_barcodes.tsv", 
                       features = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case2/Case2_ZC/Case2_ZC_features.tsv"
)
case2ZC <- CreateSeuratObject(counts = case2ZC_mtx)
case2ZC <- SetIdent(case2ZC, value = "case2ZC")
####QC####
case2ZC[["percent.mt"]] <- PercentageFeatureSet(case2ZC, pattern = "^MT-")
case2ZC[["percent.rb"]] <- PercentageFeatureSet(case2ZC, pattern = "^RP[SL]")
head(case2ZC@meta.data, 5)
VlnPlot(case2ZC, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(case2ZC, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(case2ZC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
case2ZC <- subset(case2ZC, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 10)

case2ZC <- NormalizeData(case2ZC, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
case2ZC <- FindVariableFeatures(case2ZC, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(case2ZC), 10)
plot1 <- VariableFeaturePlot(case2ZC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(case2ZC)
case2ZC <- ScaleData(case2ZC, features = all.genes)

#Regresar cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
case2ZC<- CellCycleScoring(case2ZC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
case2ZC <- RunPCA(case2ZC, features = c(s.genes, g2m.genes))
Idents(object = case2ZC) <- "Phase"
DimPlot(case2ZC, pt.size = 1.1)
case2ZC <- ScaleData(case2ZC, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(case2ZC))



case2ZC <- RunPCA(case2ZC, features = VariableFeatures(object = case2ZC))
ElbowPlot(case2ZC)

case2ZC <- FindNeighbors(case2ZC, dims = 1:15)
case2ZC <- FindClusters(case2ZC, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
clustree(case2ZC)
Idents(object = case2ZC) <- "RNA_snn_res.0.2"
case2ZC <- RunUMAP(case2ZC, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
seurat_case2_ZC_umap <- DimPlot(case2ZC, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)

saveRDS(case2ZC, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_ZC.rds")
case2ZC <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_ZC.rds")




SSc <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"), )

bc_case2ZC <- bcScore(case2ZC, SSc, expr.thres = 0.1)
bc_filtered <- bcSubset(bc_case2ZC, nan.cells = 0.95)
bc_filtered@normalized[is.na(bc_filtered@normalized)] <- 0                                                                                                                                                                       
bc_recomputed <- bcRecompute(bc_filtered, slot = "normalized") 

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)

bc_recomputed <- bcRegressOut(bc_recomputed, vars.to.regress = "nFeature_RNA", k.neighbors = 40)

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
gemcitabine <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = "sig-20902"), pt.size = 2.5)
gemcitabine+theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                  panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
                  plot.subtitle = element_text(size=15))

resistance_markers <- as.data.frame(bc_functional_ZC@normalized["HALLMARK_XENOBIOTIC_METABOLISM", ])
resistance_markers$BCscore <- (bc_case2ZC@normalized["sig-20902", ])
colnames(resistance_markers) <- c("xenobiotic_metabolism", "BCscore")

resistance_markers$xenobiotic_metabolism[resistance_markers$xenobiotic_metabolism < 0] <- NA

ggplot(resistance_markers, aes(x=xenobiotic_metabolism, y=BCscore)) + 
  geom_point()+ 
  stat_cor(method = "spearman", label.y = 60, size=7)+ geom_smooth(method = "glm", se = FALSE)+ 
  ggtitle("Correlation between xenobiotic metabolism and gemcitabine BCscore in metastasis sample")+theme(plot.title = element_text(size = 20), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"), 
                                                                                                          panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1.7), 
                                                                                                          plot.subtitle = element_text(size=15))











###CELL TYPE###

library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
library(scRNAseq)
library(scuttle)
sceM <- MuraroPancreasData()
sceM <- sceM[,!is.na(sceM$label)]
sceM <- logNormCounts(sceM)

# ensemble_id <- row.names(sceM)
# ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# attributes = listAttributes(ensembl)
# listDatasets(ensembl)
# listFilters(ensembl)
# hgnc_symbol <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), values = ensemble_id, mart = ensembl, uniqueRows = TRUE, filters = "ensembl_gene_id")
# hgnc_symbol <- arrange(hgnc_symbol, hgnc_symbol)
# hgnc_symbol <- hgnc_symbol[!duplicated(hgnc_symbol$hgnc_symbol),]
# 
# nuevo_dataframe <- merge(hgnc_symbol, row.names(sceM), by.x = "ensembl_gene_id", by.y = "row.names", all.x = TRUE)
# # Eliminar la columna adicional creada por el merge
# nuevo_dataframe <- nuevo_dataframe[, -3]
# 
# rownames(nuevo_dataframe) <- c("hgnc_symbol", sceM$label)
# nuevo_dataframe[nuevo_dataframe == ""] <- NA
# sceM_final <- na.omit(nuevo_dataframe)
# row.names(sceM) <- nuevo_dataframe$hgnc_symbol


Idents(object = case2ZC) <- "RNA_snn_res.0.2"
case2ZC_celltype <- SingleR(GetAssayData(case2ZC, assay = "RNA", slot = "data"),
                            clusters = Idents(case2ZC), ref = hpca.se,  labels =hpca.se$label.fine, de.method = "wilcox") 

case2ZC[["cell_type"]] <- case2ZC_celltype$labels[match(case2ZC[[]][["RNA_snn_res.0.2"]], rownames(case2ZC_celltype))]
Idents(object = case2ZC) <- "cell_type"
case2ZC <- RunUMAP(case2ZC, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
cell_type_case2_ZC_umap <- DimPlot(case2ZC, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)

# Idents(object = case2ZC) <- "RNA_snn_res.0.2"
# case2ZC_celltype2 <- SingleR(GetAssayData(case2ZC, assay = "RNA", slot = "data"),
#                              clusters = Idents(case2ZC), ref = sceM,  labels =sceM$label, de.method = "wilcox") 
# case2ZC[["cell_type2"]] <- case2ZC_celltype2$labels[match(case2ZC[[]][["RNA_snn_res.0.2"]], rownames(case2ZC_celltype2))]
# Idents(object = case2ZC) <- "cell_type2"
# case2ZC <- RunUMAP(case2ZC, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
# cell_type2_case2_ZC_umap <- DimPlot(case2ZC, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)

annotated_celltype <- case2ZC$cell_type
write.table(x = annotated_celltype, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/data/Case2/Case2_ZC/annotated_celltype_case2ZC.tsv", sep = "\t")

gs_functional <- GenerateGenesets(x = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/data/genesets/h.all.v2023.1.Hs.symbols.gmt")
bc_functional_ZC <- bcScore(case2ZC, gs_functional, expr.thres = 0.1)
resistance_markers <- as.data.frame(bc_functional_ZC@normalized["HALLMARK_XENOBIOTIC_METABOLISM", ])
resistance_markers$BCscore <- (bc_case2ZC@normalized["sig-20902", ])
colnames(resistance_markers) <- c("xenobiotic_metabolism", "BCscore")

resistance_markers$xenobiotic_metabolism[resistance_markers$xenobiotic_metabolism < 0] <- NA

ggplot(resistance_markers, aes(x=xenobiotic_metabolism, y=BCscore)) + 
  geom_point()+ 
  stat_cor(method = "spearman", label.y = 60, size=7)+ geom_smooth(method = "glm", se = FALSE)+ 
  ggtitle("Correlation between xenobiotic metabolism and gemcitabine BCscore in normal sample")+theme(plot.title = element_text(size = 20), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"), 
                                                                                                          panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1.7), 
                                                                                                          plot.subtitle = element_text(size=15))











###MERGED###
case2YF$sample = "tumor"
case2ZY$sample = "meta"
case2ZC$sample = "normal"
case2_merged <- merge(case2YF, y = c(case2ZY, case2ZC), project = "case2_PDAC")
case2_merged <- RenameIdents(case2_merged, 'YF' = 'tumor', 'ZY' = 'meta', 'ZC' = 'normal')
saveRDS(case2_merged, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_merged_raw.rds")
case2_merged <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_merged_raw.rds")

Idents(object = case2_merged) <- "sample"

case2_merged[["percent.mt"]] <- PercentageFeatureSet(case2_merged, pattern = "^MT-")
case2_merged[["percent.rb"]] <- PercentageFeatureSet(case2_merged, pattern = "^RP[SL]")
head(case2_merged@meta.data, 5)
palette <- rep(c("white", "red", "blue"), length.out = 3)
VlnPlot(case2_merged, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(case2_merged, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(case2_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
case2_merged <- subset(case2_merged, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 10)


case2_merged <- NormalizeData(case2_merged, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
case2_merged <- FindVariableFeatures(case2_merged, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(case2_merged), 10)
plot1 <- VariableFeaturePlot(case2_merged)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(case2_merged)
case2_merged <- ScaleData(case2_merged, features = all.genes)

#Regresar cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
case2_merged<- CellCycleScoring(case2_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
case2_merged <- RunPCA(case2_merged, features = c(s.genes, g2m.genes))
Idents(object = case2_merged) <- "Phase"
DimPlot(case2_merged, pt.size = 1.1)
case2_merged <- ScaleData(case2_merged, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(case2_merged))



case2_merged <- RunPCA(case2_merged, features = VariableFeatures(object = case2_merged))
ElbowPlot(case2_merged)

case2_merged <- FindNeighbors(case2_merged, dims = 1:17)
case2_merged <- FindClusters(case2_merged, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
clustree(case2_merged)
Idents(object = case2_merged) <- "RNA_snn_res.0.4"


case2_merged <- RunUMAP(case2_merged, dims = 1:17, min.dist = 0.15, n.neighbors = 30)
seurat_case2_merged_umap30 <- DimPlot(case2_merged, reduction = "umap", label = FALSE, label.size = 5, pt.size = 1.1)+
  theme(legend.text = element_text(size=21), plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                                        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2))

Idents(object = case2_merged) <- "sample"
sample_case2_merged_umap <- DimPlot(case2_merged, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1)

saveRDS(case2_merged, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_merged.rds")
case2_merged <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_merged.rds")


case2_markers <- FindAllMarkers(case2_merged, only.pos = TRUE, logfc.threshold = 1, min.pct = 0.15)
case2_YF_markers <- FindMarkers(case2_merged, ident.1 = "YF", only.pos = TRUE, logfc.threshold = 1, min.pct = 0.15)
DoHeatmap(case2_merged, features = case2_markers$gene, draw.lines = TRUE, 
          group.colors = c("#440154", "#21908c", "#fde725")) + guides(color="none")
DoHeatmap(case2_merged, features = rownames(case2_YF_markers))


CD44_umap <- FeaturePlot(case2_merged, features = "MUC1", pt.size = 1.2)
HMGA1_umap <- FeaturePlot(case2_merged, features = "HMGA1", pt.size = 1.2)

HMGA1_umap+theme(legend.text = element_text(size=21), plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2))
CD44_umap+theme(legend.text = element_text(size=21), plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                 panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2))


DoHeatmap(case2_merged, features = c("MUC1", "DCK", "CD44","RRM1", "RRM2","HMGA1", "NOTCH1", "NOTCH2","NOTCH3", 
                                     "BRCA2", "ROCK2", "SLC29A1", "PTK2", "BNIP3"))

DoHeatmap(case2_merged, features = c("MAPK1", "MAP2K1", "MAP2K2","DCK", "CDA", "RRM1", "RRM2", "HMGA1"))

DoHeatmap(case2_merged, features = c("CCDC148", "SH3RF2", "CACNA1D", "POLD3", "PARP1", "AP1M2", "C4ORF19", "ANO1", "VGLL1", "SCEL", "INPP4B", "NET1", "INSIG2", "BVES"))

library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
library(scRNAseq)
library(scuttle)
sceM <- MuraroPancreasData()
sceM <- sceM[,!is.na(sceM$label)]
sceM <- logNormCounts(sceM)
rownames(sceM) <- sub("[_].*", "", rownames(sceM))



Idents(object = case2_merged) <- "RNA_snn_res.0.4"
case2_merged_celltype <- SingleR(GetAssayData(case2_merged, assay = "RNA", slot = "data"),
                            clusters = Idents(case2_merged), ref = hpca.se,  labels =hpca.se$label.main, de.method = "wilcox") 

case2_merged[["cell_type"]] <- case2_merged_celltype$labels[match(case2_merged[[]][["RNA_snn_res.0.2"]], rownames(case2_merged_celltype))]
Idents(object = case2_merged) <- "cell_type"
case2_merged <- RunUMAP(case2_merged, dims = 1:17, min.dist = 0.25, n.neighbors = 30)
cell_type_case2_merged_umap <- DimPlot(case2_merged, reduction = "umap", label = TRUE, label.size = 4, pt.size = 1.2)

Idents(object = case2_merged) <- "RNA_snn_res.0.2"
case2_merged_celltype2 <- SingleR(GetAssayData(case2_merged, assay = "RNA", slot = "data"),
                             clusters = Idents(case2_merged), ref = sceM,  labels =sceM$label, de.method = "wilcox") 
case2_merged[["cell_type2"]] <- case2_merged_celltype2$labels[match(case2_merged[[]][["RNA_snn_res.0.2"]], rownames(case2_merged_celltype2))]
Idents(object = case2_merged) <- "cell_type2"
case2_merged <- RunUMAP(case2_merged, dims = 1:17, min.dist = 0.25, n.neighbors = 30)
cell_type2_case2_ZY_umap <- DimPlot(case2_merged, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2, )

###CORRELATION###
case2_merged <- AddModuleScore(case2_merged, features = list(c("MUC1")), name = "resistance_markers", nbin = 1)

resistance_markers <- as.data.frame(bc_functional@normalized["HALLMARK_XENOBIOTIC_METABOLISM", ])
resistance_markers <- as.data.frame(case2_merged@meta.data$resistance_markers1)
resistance_markers$BCscore <- (bc_case2@normalized["sig-20902", ])
colnames(resistance_markers) <- c("MUC1", "BCscore")
resistance_markers$MUC1[resistance_markers$MUC1 < 0] <- NA

cor_resistance <- cor(x = resistance_markers)
matrix_resistance_markers <- as.matrix(resistance_markers)
cor_resistanceP <- rcorr(x = matrix_resistance_markers, type = "pearson")

ggplot(resistance_markers, aes(x=MUC1, y=BCscore)) + 
  geom_point()+ 
  stat_cor(method = "spearman", label.y = 60, size=7)+ geom_smooth(method = "glm", se = FALSE)+ 
  ggtitle("Correlation between MUC1 expression and gemcitabine BCscore in PDAC samples")+theme(plot.title = element_text(size = 24), axis.text=element_text(size=15), axis.title=element_text(size=18,face="bold"), 
                                                                                                  panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1.7), 
                                                                                                  plot.subtitle = element_text(size=15))











################BEYONDCELL#################
SSc <- GetCollection(SSc, n.genes = 250, mode = c("up", "down"), )

bc_case2 <- bcScore(case2_merged, SSc, expr.thres = 0.1)

saveRDS(bc_case2, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/bc_case2.rds")
bc_case2 <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/bc_case2.rds")

bc_filtered <- bcSubset(bc_case2, nan.cells = 0.95)
bc_filtered@normalized[is.na(bc_filtered@normalized)] <- 0                                                                                                                                                                       
bc_recomputed <- bcRecompute(bc_filtered, slot = "normalized") 


bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = NULL, k.neighbors = 40)
bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)


bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)
bc_recomputed <- bcRegressOut(bc_recomputed, vars.to.regress = "nFeature_RNA", k.neighbors = 40)

bc_recomputed <- bcUMAP(bc_recomputed, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "nCount_RNA", factor.col = FALSE)

saveRDS(bc_recomputed, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/bc_case2_recomputed.rds")
bc_recomputed <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/bc_case2_recomputed.rds")


clustree(bc_recomputed@meta.data, prefix ="bc_clusters_res.")
bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "sample", label.size = 5, pt.size = 1.5)
case2_merged <- AddMetaData(object = case2_merged, metadata = bc_recomputed@meta.data)

bcClusters(bc_recomputed, UMAP = "Seurat", idents = "RNA_snn_res.0.4", pt.size = 1.5, label = TRUE)
DimPlot(object = case2_merged)

DimPlot(
  case2_merged,
  reduction = "beyondcell",
  label = TRUE,
  label.size = 5,
  group.by = "bc_clusters_res.0.4",
  shape.by = "type",
  raster = FALSE,
  pt.size = 0.8
)

bc_samples <- bcClusters(bc_recomputed, UMAP = "beyondcell", idents = "sample", pt.size = 2)

bc_samples+theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                 panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
                 plot.subtitle = element_text(size=15))




dabrafenib_info <- FindDrugs(bc_recomputed, "dabrafenib")

dabrafenib_IDs <- FindDrugs(bc_recomputed, "dabrafenib")$IDs
trametinib_IDs <- FindDrugs(bc_recomputed, "trametinib")$IDs
cobimetinib_IDs <- FindDrugs(bc_recomputed, "cobimetinib")$IDs


dabrafenib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                           signatures = list(values = "sig-21060"), pt.size = 2)
dabrafenib+theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                 panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
                 plot.subtitle = element_text(size=15))

trametinib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                           signatures = list(values = "sig-20908"), pt.size = 2.5)
trametinib+theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                 panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
                 plot.subtitle = element_text(size=15))

cobimetinib <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = cobimetinib_IDs), pt.size = 2)
cobimetinib+theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                 panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
                 plot.subtitle = element_text(size=15))

gemcitabine <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                            signatures = list(values = "sig-20902"), pt.size = 2.5)
gemcitabine+theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
                 panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
                 plot.subtitle = element_text(size=15))

wrap_plots(dabrafenib, ncol = 1)
wrap_plots(cobitinib, ncol = 1)
wrap_plots(cobimetinib, ncol = 1)

##Boxplots
trametinib_bc <- CreatebcObject(bc = bc_recomputed)
trametinib_bc <- bcSubset(bc = trametinib_bc, signatures = c("sig-21037"))
trame_frame <- data.frame(trametinib_bc@normalized)
trame_frame <- as.data.frame(case2_merged$resistance_markers1)
colnames(trame_frame) <- c("MUC1")
trame_frame <- data.frame(t(trame_frame))
rownames(trame_frame) <- gsub("\\.", "-", rownames(trame_frame))


sample <- data.frame(case2_merged$sample)
sample$case2_merged.sample <- gsub("YF", "tumor", sample$case2_merged.sample)
sample$case2_merged.sample <- gsub("ZY", "meta", sample$case2_merged.sample)
sample$case2_merged.sample <- gsub("ZC", "normal", sample$case2_merged.sample)

trame_frame_type <- merge(x = trame_frame, y = sample, by= 0)
trame_frame_long <- pivot_longer(data = trame_frame_type, cols = c("MUC1"), names_to = "sig", values_to = "enrichment_score")
colnames(trame_frame_long) <- c("Cells", "sample", "sig", "enrichment_score")


trame_frame_long$sample <- factor(trame_frame_long$sample, levels=c("tumor", "meta", "normal"), ordered = TRUE)

anova.test <- trame_frame_long %>% anova_test(formula = enrichment_score~sample)

pwc <- trame_frame_long %>% tukey_hsd(enrichment_score ~ sample)
pwc

boxplot_case2 <- ggplot(trame_frame_long) +
  geom_boxplot((aes(x=sample, y=enrichment_score, fill=sample)), show.legend = FALSE, outlier.shape = NA) + 
  stat_pvalue_manual(pwc, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE, y.position = c(2, 3, 3.4), size = 7)+
  ggtitle("Comparison of MUC1 expression between samples") + 
  theme(plot.title = element_text(size = 30), axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), 
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 2), 
        panel.grid = element_line(colour = "grey", linewidth = 0.25), plot.margin = margin(t=10, r = 10, b = 5, l = 5)) + ylim(-1, 4)
boxplot_case2 + scale_fill_manual(values = c("#440154", "#21908c", "#fde725", "#e7298a"))
#1b9e77


ggplot(trame_frame_long) +
  geom_boxplot((aes(x=sample, y=enrichment_score, fill=sample)), show.legend = FALSE, outlier.shape = NA)+ 
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, hide.ns = FALSE, size = 8, y.position =c(5, 6.5, 8) ) + scale_fill_brewer(palette="Dark2")+
  ggtitle("Comparison of gemcitabine signature between samples") + theme(plot.title = element_text(size = 30))

















# Obtain condition-based statistics.
bc_recomputed <- bcRanks(bc_recomputed, idents = "sample")
#Obtain unextended therapeutic cluster-based statistics.
bc_recomputed <- bcRanks(bc_recomputed, idents = "sample", extended = FALSE)


bc4Squares(bc_recomputed, idents = "sample", lvl = NULL, top = 5)
bc4Squares(bc_recomputed, idents = "bc_clusters_res.0.2", lvl = NULL, top = 3)


YF1 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20923"), pt.size = 1.5)
YF2 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21232"), pt.size = 1.5)
YF3 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20912"), pt.size = 1.5)
YF4 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20937"), pt.size = 1.5)
YF5 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21414"), pt.size = 1.5)

top_diff_YF <- bc_samples+YF1+YF2+YF3+YF4+YF5


YF6 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21213"), pt.size = 1.5)
YF7 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21456"), pt.size = 1.5)
YF8 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21288"), pt.size = 1.5)
YF9 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21105"), pt.size = 1.5)
YF10 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                     signatures = list(values = "sig-21058"), pt.size = 1.5)

low_diff_YF <- bc_samples+YF6+YF7+YF8+YF9+YF10






ZC1 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21213"), pt.size = 1.5)
ZC2 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21105"), pt.size = 1.5)
ZC3 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21456"), pt.size = 1.5)
ZC4 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21288"), pt.size = 1.5)
ZC5 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21058"), pt.size = 1.5)

top_diff_ZC <- bc_samples+ZC1+ZC2+ZC3+ZC4+ZC5


ZC6 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20923"), pt.size = 1.5)
ZC7 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21232"), pt.size = 1.5)
ZC8 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20937"), pt.size = 1.5)
ZC9 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21414"), pt.size = 1.5)
ZC10 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                     signatures = list(values = "sig-20912"), pt.size = 1.5)

low_diff_ZC <- bc_samples+ZC6+ZC7+ZC8+ZC9+ZC10


ZY1 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20922"), pt.size = 1.5)
ZY2 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20923"), pt.size = 1.5)
ZY3 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20912"), pt.size = 1.5)
ZY4 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20937"), pt.size = 1.5)
ZY5 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-20905"), pt.size = 1.5)

top_diff_ZY <- bc_samples+ZY1+ZY2+ZY3+ZY4+ZY5


ZY6 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21456"), pt.size = 1.5)
ZY7 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21195"), pt.size = 1.5)
ZY8 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21058"), pt.size = 1.5)
ZY9 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                    signatures = list(values = "sig-21288"), pt.size = 1.5)
ZY10 <- bcSignatures(bc_recomputed, UMAP = "beyondcell", 
                     signatures = list(values = "sig-21322"), pt.size = 1.5)

low_diff_ZY <- bc_samples+ZY6+ZY7+ZY8+ZY9+ZY10












####Funcional####

# hallmarks <- read.table(file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/data/genesets/h.all.v2023.1.Hs.symbols.gmt", header = F, sep = "\t", fill = T)
# hallmarks$V1=paste(hallmarks$V1,"_UP", sep = "")
# write.table(hallmarks, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/data/genesets/h.all.v2023.1.Hs.symbols.gmt", sep = "\t", row.names = F, col.names = F, quote = F)


gs_functional <- GenerateGenesets(x = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/data/genesets/h.all.v2023.1.Hs.symbols.gmt")
bc_functional <- bcScore(case2_merged, gs_functional, expr.thres = 0.1)
bc_filtered_functional <- bcSubset(bc_functional, nan.cells = 0.95)
bc_filtered_functional@normalized[is.na(bc_filtered_functional@normalized)] <- 0                                                                                                                                                                       
bc_recomputed_functional <- bcRecompute(bc_filtered_functional, slot = "normalized") 

bc_recomputed_functional <- bcUMAP(bc_recomputed_functional, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40, add.DSS = FALSE)
bc_recomputed_functional <- bcRegressOut(bc_recomputed_functional, vars.to.regress = "nFeature_RNA", k.neighbors = 40)
bc_recomputed_functional <- bcUMAP(bc_recomputed_functional, res = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4), pc = 10, k.neighbors = 40)

clustree(bc_recomputed_functional@meta.data, prefix ="bc_clusters_res.")
bc_umap <- bcClusters(bc_recomputed_functional, UMAP = "beyondcell", idents = "bc_clusters_res.0.4", label.size = 5, pt.size = 1.3)
bc_sample_umap <- bcClusters(bc_recomputed_functional, UMAP = "beyondcell", idents = "sample", label.size = 5, pt.size = 1.2)

saveRDS(bc_functional, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_merged_funcionalRAW.rds")
bc_functional <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_merged_funcionalRAW.rds")

saveRDS(bc_recomputed_functional, file = "/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_merged_funcional.rds")
bc_recomputed_functional <- readRDS("/home/lpgonzalezm/TFM/beyondcell/PDAC/seurat/Case2/case2_merged_funcional.rds")



# data <- data.frame(bc_recomputed_functional@normalized)
# clusters <- data.frame(case2_merged$RNA_snn_res.0.2)
# sample <- data.frame(case2_merged$sample)
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
                      signatures = list(values = "HALLMARK_MITOTIC_SPINDLE"), pt.size = 0.8)

plot2 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_HYPOXIA"), pt.size = 0.8)

plot3 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_ALLOGRAFT_REJECTION"), pt.size = 0.8)

plot4 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_APOPTOSIS"), pt.size = 1)

plot5 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_DNA_REPAIR"), pt.size = 0.8)

plot6 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_FATTY_ACID_METABOLISM"), pt.size = 0.8)

plot7 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_MYC_TARGETS_V1"), pt.size = 0.8)

plot8 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"), pt.size = 0.8)

plot9 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                      signatures = list(values = "HALLMARK_P53_PATHWAY"), pt.size = 0.8)

resumen_plots <- wrap_plots(plot1, plot2, plot3, plot4)
resumen_plots2 <- wrap_plots(plot5, plot7, plot8, plot9)


plot12 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"), pt.size = 0.8)

plot13 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_TGF_BETA_SIGNALING"), pt.size = 0.8)

plot14 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_UV_RESPONSE"), pt.size = 0.8)

plot15 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_XENOBIOTIC_METABOLISM"), pt.size = 0.8)

plot11 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"), pt.size = 0.8)

plot16 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"), pt.size = 0.8)

plot17 <- bcSignatures(bc_recomputed_functional, UMAP = "beyondcell", 
                       signatures = list(values = "HALLMARK_MTORC1_SIGNALING"), pt.size = 0.8)

resumen_plots2 <-wrap_plots(plot11, plot12, plot13, plot15, plot16, plot17)

##Boxplots funcionales##

boxplot_bc <- CreatebcObject(bc = bc_recomputed_functional)
MAP2K1_bc <- bcSubset(bc = boxplot_bc, signatures = c("HALLMARK_XENOBIOTIC_METABOLISM"))
MAP2K1_frame <- data.frame(MAP2K1_bc@normalized)
rownames(MAP2K1_frame) <- c("HALLMARK_XENOBIOTIC_METABOLISM")
MAP2K1_frame <- data.frame(t(MAP2K1_frame))
sample <- data.frame(case2_merged$sample)
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
colnames(MAP2K1_frame_model) <- c("Cells", "HALLMARK_XENOBIOTIC_METABOLISM", "sample")
MAP2K1_frame_long <- pivot_longer(data = MAP2K1_frame_model, cols = c("HALLMARK_XENOBIOTIC_METABOLISM"), names_to = "sig", values_to = "enrichment_score")

MAP2K1_frame_long$sample <- factor(MAP2K1_frame_long$sample, levels=c("tumor", "meta", "normal"), ordered = TRUE)
# TCs <- factor(MAP2K1_frame_long$TCs, levels=c("0", "1", "2", "3", "4", "5"))



# stat.test <- MAP2K1_frame_long %>% t_test((enrichment_score~sample))
# stat.test <- stat.test %>% add_xy_position(x = "sample")

MAP2K1_frame_long$sample <- as.factor(MAP2K1_frame_long$sample)
anova.test <- MAP2K1_frame_long %>% anova_test(formula = enrichment_score~sample)

pwc <- MAP2K1_frame_long %>% tukey_hsd(enrichment_score ~ sample)
pwc
options("scipen" = 999, "digits" = 22)
options("scipen" = 0, "digits" = 7)


ggplot(MAP2K1_frame_long) +
  geom_boxplot((aes(x=sample, y=enrichment_score, fill=sample)), show.legend = FALSE, outlier.shape = NA) + 
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, hide.ns = FALSE, y.position = c(20, 25, 30)) + scale_fill_brewer(palette="Dark2") +
  ggtitle("Comparison of mitotic spindle signature between samples") + 
  theme(plot.title = element_text(size = 30), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))


boxplot_funcional_case2 <- ggplot(MAP2K1_frame_long) +
  geom_boxplot((aes(x=sample, y=enrichment_score, fill=sample)), show.legend = FALSE, outlier.shape = NA) + 
  stat_pvalue_manual(pwc, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE, y.position = c(11, 14, 17)) +
  labs(subtitle = get_test_label(anova.test, detailed = TRUE))+
  ggtitle("Comparison of xenobiotic metabolism signature between samples") + 
  theme(plot.title = element_text(size = 30), axis.text=element_text(size=16), axis.title=element_text(size=19,face="bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", linewidth = 2, fill = NA),
        panel.grid = element_line(colour = "grey", linewidth = 0.25)) + ylim(-10, 20)
boxplot_funcional_case2 + scale_fill_manual(values = c("#440154", "#21908c", "#fde725", "#e7298a"))


weight_aov = aov(enrichment_score ~ sample, data = MAP2K1_frame_long)
summary(weight_aov)
#Homogeneity of variances
plot(weight_aov, 1)
#Normality
plot(weight_aov, 2)


####OBJETO BCMETA####

case2YF <- pdac_shu_zhang[[4]]
case2ZY <- pdac_shu_zhang[[5]]
case2ZC <- pdac_shu_zhang[[3]]





