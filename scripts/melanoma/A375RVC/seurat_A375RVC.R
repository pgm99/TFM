library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(biomaRt)
library(curl)
library(tidyverse)
library(ggplot2)
library(ggalluvial)
A375RVC <- read.table("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/data/GSM5022597_raw_counts_probe5.txt.gz", header = TRUE)
saveRDS(A375RVC_data, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/A375RVC.rds")
A375RVC_data <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/A375RVC.rds")


# Anotar hgnc symbol en vez de ensembl
ensemble_id <- row.names(A375RVC)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl)
listDatasets(ensembl)
listFilters(ensembl)
hgnc_symbol <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), values = ensemble_id, mart = ensembl, uniqueRows = TRUE, filters = "ensembl_gene_id")
hgnc_symbol <- arrange(hgnc_symbol, hgnc_symbol)
hgnc_symbol <- hgnc_symbol[!duplicated(hgnc_symbol$hgnc_symbol),]

nuevo_dataframe <- merge(hgnc_symbol, A375RVC, by.x = "ensembl_gene_id", by.y = "row.names", all.x = TRUE)
# Eliminar la columna adicional creada por el merge
nuevo_dataframe <- nuevo_dataframe[, -1]

colnames(nuevo_dataframe) <- c("hgnc_symbol", colnames(A375RVC))
nuevo_dataframe[nuevo_dataframe == ""] <- NA
A375RVC_final <- na.omit(nuevo_dataframe)
row.names(A375RVC_final) <- A375RVC_final$hgnc_symbol
A375RVC_final <- A375RVC_final[,-1]

write.table(A375RVC_final, file= "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/data/GSM5022597_raw_counts_probe5_anotado.txt.gz", row.names = TRUE, col.names = TRUE)

#objeto Seurat
A375RVC_data <- CreateSeuratObject(counts= A375RVC_final)
A375RVC_data <- SetIdent(A375RVC_data, value = "A375RVC")
####QC####
A375RVC_data[["percent.mt"]] <- PercentageFeatureSet(A375RVC_data, pattern = "^MT-")
A375RVC_data[["percent.rb"]] <- PercentageFeatureSet(A375RVC_data, pattern = "^RP[SL]")
head(A375RVC_data@meta.data, 5)
VlnPlot(A375RVC_data, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(A375RVC_data, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(A375RVC_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
A375RVC_data <- subset(A375RVC_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)




#Normalización
A375RVC_data <- NormalizeData(A375RVC_data, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
A375RVC_data <- FindVariableFeatures(A375RVC_data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(A375RVC_data), 10)
plot1 <- VariableFeaturePlot(A375RVC_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(A375RVC_data)
A375RVC_data <- ScaleData(A375RVC_data, features = all.genes)

#Cell cycling
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
A375RVC_data<- CellCycleScoring(A375RVC_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
A375RVC_data <- RunPCA(A375RVC_data, features = c(s.genes, g2m.genes))
Idents(object = A375RVC_data) <- "Phase"
DimPlot(A375RVC_data)

A375RVC_data <- RunUMAP(A375RVC_data, dims = 1:15, n.neighbors = 40, min.dist = 0)
DimPlot(A375RVC_data, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.2)

A375RVC_data <- ScaleData(A375RVC_data, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(A375RVC_data))

A375RVC_data <- SetIdent(A375RVC_data, value = "A375RVC")




#PCA y formas de visualizarlo
A375RVC_data <- RunPCA(A375RVC_data, features = VariableFeatures(object = A375RVC_data))
#print(A375RVC_data[["pca"]], dims = 1:10, nfeatures = 5)
#VizDimLoadings(A375RVC_data, dims = 1:2, reduction = "pca")
#DimPlot(A375RVC_data, reduction = "pca")

#VizDimLoadings(A375RVC_data, dims = 1:5, reduction = "pca")
#DimPlot(A375RVC_data, reduction = "pca", dims = 1:2)

#Heatmaps
DimHeatmap(A375RVC_data, dims = 1:10, cells = 500, balanced = TRUE)

#Determinar dimensionalidad de PC con JackStraw y con ElbowPlot
A375RVC_data <- JackStraw(A375RVC_data, num.replicate = 1000)
A375RVC_data <- ScoreJackStraw(A375RVC_data, dims = 1:20)
JackStrawPlot(A375RVC_data, dims = 1:20)
ElbowPlot(A375RVC_data)

#Clustering
A375RVC_data <- FindNeighbors(A375RVC_data, dims = 1:15)
A375RVC_data <- FindClusters(A375RVC_data, resolution = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4))

clustree(A375RVC_data)
Idents(object = A375RVC_data) <- "RNA_snn_res.0.6"

#UMAP
A375RVC_data <- RunUMAP(A375RVC_data, dims = 1:15, n.neighbors = 40, min.dist = 0)
DimPlot(A375RVC_data, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.5)

#UMAP para los clusters de BC
Idents(object = A375RVC_data) <- "bc_clusters_res.0.4"
A375RVC_data <- RunUMAP(A375RVC_data, dims = 1:15)
DimPlot(A375RVC_data, reduction = "umap", label = TRUE, label.size = 5)

DimPlot(
  A375RVC_data,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  group.by = c("RNA_snn_res.0.4", "bc_clusters_res.0.6"),
)



n_cells <- table(A375RVC_data@meta.data$RNA_snn_res.0.4, A375RVC_data@meta.data$orig.ident)

#Biomarkers
A375RVC_data.markers <- FindAllMarkers(A375RVC_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, return.thresh = 0.0000025)
A375RVC_data.allmarkers <- FindAllMarkers(A375RVC_data, only.pos = TRUE)

A375RVC_top_markers <- A375RVC_data.markers %>% group_by(cluster) %>% slice_max(n = 12, order_by = avg_log2FC)

DoHeatmap(A375RVC_data, features = A375RVC_data.markers$gene)

A375RVC_annotated_markers <- getBM(attributes = c("hgnc_symbol","description", "start_position", "end_position", "chromosome_name", "gene_biotype"), 
                                  values = A375RVC_data.markers$gene, mart = ensembl, uniqueRows = TRUE, filters = "hgnc_symbol")
A375RVC_annotated_markers <- A375RVC_annotated_markers[-c(29,31,35), ]

markers_list <- data.frame(A375RVC_data.allmarkers$gene, A375RVC_data.allmarkers$avg_log2FC)
write.csv(x = A375RVC_annotated_markers, file = "A375RVC_annotated_markers.csv", row.names = TRUE)
write.table(markers_list, file = "A375RVC_markers_list.rnk", row.names = F, quote = F, col.names = F, sep = "\t")


#Subclones SCEVAN
scevan <- readRDS(file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RVC/scevan_A375RVC.rds")
scevan_subclone1 <- read.table(file = "/home/lpgonzalezm/output_old/A375RVC_subclone1_CN.seg", header = TRUE, sep = "\t")
scevan$subclone <- as.factor(scevan$subclone)
rownames(scevan) <- gsub(pattern = "-", replacement = ".", x = rownames(scevan))
scevan_filtered <- scevan[colnames(A375RVC_data), ]
A375RVC_data <- AddMetaData(object = A375RVC_data, metadata = scevan_filtered)

DimPlot(
  A375RVC_data,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  group.by = c("RNA_snn_res.0.6", "subclone")
  #  shape.by = "subclone",
)




#SANKEY
samples <- colnames(A375RVC_data)
clusters <- A375RVC_data@meta.data$RNA_snn_res.0.6
subclones <- A375RVC_data@meta.data$subclone
gg_data <- A375RVC_data@meta.data %>%
  group_by(colnames(A375RVC_data), clusters, subclones) %>%
  tally() 


ggplot(data = gg_data,
       aes(axis1 = subclones, axis2 = clusters)) +
  geom_alluvium(aes(fill = subclones), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("subclone", "RNA_snn_res.0.6"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Distribución de subclones en los clusters")
# ggplot(data = gg_data,
#        aes(y = "subclone", axis = "RNA_snn_res.0.6"))

gg_data_filtrada <- subset(gg_data, !(subclone == 3 & clusters == 1), )



ggplot(data = gg_data_filtrada,
       aes(axis1 = subclones, axis2 = clusters)) +
  geom_alluvium(aes(fill = subclones), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("subclone", "RNA_snn_res.0.6"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Distribución de subclones en los clusters")


