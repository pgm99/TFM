library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(biomaRt)
library(curl)
library(tidyverse)
library(ggplot2)
library(ggalluvial)
A375RV <- read.table("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/data/GSM5022595_raw_counts_probe3.txt.gz", header = TRUE)
saveRDS(A375RV_data, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RV/A375RV.rds")
A375RV_data <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375RV/A375RV.rds")


# Anotar hgnc symbol en vez de ensembl
ensemble_id <- row.names(A375RV)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl)
listDatasets(ensembl)
listFilters(ensembl)
hgnc_symbol <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), values = ensemble_id, mart = ensembl, uniqueRows = TRUE, filters = "ensembl_gene_id")
hgnc_symbol <- arrange(hgnc_symbol, hgnc_symbol)
hgnc_symbol <- hgnc_symbol[!duplicated(hgnc_symbol$hgnc_symbol),]

nuevo_dataframe <- merge(hgnc_symbol, A375RV, by.x = "ensembl_gene_id", by.y = "row.names", all.x = TRUE)
# Eliminar la columna adicional creada por el merge
nuevo_dataframe <- nuevo_dataframe[, -1]

colnames(nuevo_dataframe) <- c("hgnc_symbol", colnames(A375RV))
nuevo_dataframe[nuevo_dataframe == ""] <- NA
A375RV_final <- na.omit(nuevo_dataframe)
row.names(A375RV_final) <- A375RV_final$hgnc_symbol
A375RV_final <- A375RV_final[,-1]

write.table(A375RV_final, file= "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/data/GSM5022595_raw_counts_probe3_anotado.txt.gz", row.names = TRUE, col.names = TRUE)

#objeto Seurat
A375RV_data <- CreateSeuratObject(counts= A375RV_final)
A375RV_data <- SetIdent(A375RV_data, value = "A375RV")
####QC####
A375RV_data[["percent.mt"]] <- PercentageFeatureSet(A375RV_data, pattern = "^MT-")
A375RV_data[["percent.rb"]] <- PercentageFeatureSet(A375RV_data, pattern = "^RP[SL]")
head(A375RV_data@meta.data, 5)
VlnPlot(A375RV_data, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(A375RV_data, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(A375RV_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
A375RV_data <- subset(A375RV_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)




#Normalización
A375RV_data <- NormalizeData(A375RV_data, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
A375RV_data <- FindVariableFeatures(A375RV_data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(A375RV_data), 10)
plot1 <- VariableFeaturePlot(A375RV_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(A375RV_data)
A375RV_data <- ScaleData(A375RV_data, features = all.genes)

#Cell cycling
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
A375RV_data<- CellCycleScoring(A375RV_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
A375RV_data <- RunPCA(A375RV_data, features = c(s.genes, g2m.genes))
Idents(object = A375RV_data) <- "Phase"
DimPlot(A375RV_data)

A375RV_data <- ScaleData(A375RV_data, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(A375RV_data))

A375RV_data <- SetIdent(A375RV_data, value = "A375RV")




#PCA y formas de visualizarlo
A375RV_data <- RunPCA(A375RV_data, features = VariableFeatures(object = A375RV_data))
#print(A375RV_data[["pca"]], dims = 1:10, nfeatures = 5)
#VizDimLoadings(A375RV_data, dims = 1:2, reduction = "pca")
#DimPlot(A375RV_data, reduction = "pca")

#VizDimLoadings(A375RV_data, dims = 1:5, reduction = "pca")
#DimPlot(A375RV_data, reduction = "pca", dims = 1:2)

#Heatmaps
DimHeatmap(A375RV_data, dims = 1:10, cells = 500, balanced = TRUE)

#Determinar dimensionalidad de PC con JackStraw y con ElbowPlot
A375RV_data <- JackStraw(A375RV_data, num.replicate = 1000)
A375RV_data <- ScoreJackStraw(A375RV_data, dims = 1:20)
JackStrawPlot(A375RV_data, dims = 1:20)
ElbowPlot(A375RV_data)

#Clustering
A375RV_data <- FindNeighbors(A375RV_data, dims = 1:15)
A375RV_data <- FindClusters(A375RV_data, resolution = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4))

clustree(A375RV_data)
Idents(object = A375RV_data) <- "RNA_snn_res.0.6"

#UMAP
A375RV_data <- RunUMAP(A375RV_data, dims = 1:15, n.neighbors = 40, min.dist = 0)
DimPlot(A375RV_data, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.5)

#UMAP para los clusters de BC
Idents(object = A375RV_data) <- "bc_clusters_res.0.4"
A375RV_data <- RunUMAP(A375RV_data, dims = 1:15)
DimPlot(A375RV_data, reduction = "umap", label = TRUE, label.size = 5)

DimPlot(
  A375RV_data,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  group.by = c("RNA_snn_res.0.6", "bc_clusters_res.0.6"),
)



n_cells <- table(A375RV_data@meta.data$RNA_snn_res.0.6, A375RV_data@meta.data$orig.ident)

#Biomarkers
A375RV_data.markers <- FindAllMarkers(A375RV_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, return.thresh = 0.0000025 )
A375RV_data.allmarkers <- FindAllMarkers(A375RV_data, only.pos = TRUE)

A375RV_top_markers <- A375RV_data.markers %>% group_by(cluster) %>% slice_max(n = 12, order_by = avg_log2FC)

DoHeatmap(A375RV_data, features = A375RV_data.markers$gene)

annotated_markers <- getBM(attributes = c("hgnc_symbol","description", "start_position", "end_position", "chromosome_name", "gene_biotype", "entrezgene_description"), 
                           values = A375RV_data.markers$gene, mart = ensembl, uniqueRows = TRUE, filters = "hgnc_symbol")
A375RV_annotated_markers <- annotated_markers[-c(19,20,21,26,34), ]
markers_list <- data.frame(A375RV_data.allmarkers$gene, A375RV_data.allmarkers$avg_log2FC)
write.csv(x = A375RV_annotated_markers, file = "A375RV_annotated_markers.csv", row.names = TRUE)
write.table(markers_list, file = "A375RV_markers_list.rnk", row.names = F, quote = F, col.names = F, sep = "\t")


#Subclones SCEVAN
scevan <- readRDS(file = "/home/lpgonzalezm/TFM/A375RV_scevan.rds")
scevan_subclone1 <- read.table(file = "/home/lpgonzalezm/output_old/A375RV_subclone1_CN.seg", header = TRUE, sep = "\t")
scevan$subclone <- as.factor(scevan$subclone)
rownames(scevan) <- gsub(pattern = "-", replacement = ".", x = rownames(scevan))
scevan_filtered <- scevan[colnames(A375RV_data), ]
A375RV_data <- AddMetaData(object = A375RV_data, metadata = scevan)

DimPlot(
  A375RV_data,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  group.by = c("RNA_snn_res.0.6", "subclone")
  #  shape.by = "subclone",
)




#SANKEY
samples <- colnames(A375RV_data)
clusters <- A375RV_data@meta.data$RNA_snn_res.0.4
subclones <- A375RV_data@meta.data$subclone
gg_data <- A375RV_data@meta.data %>%
  group_by(colnames(A375RV_data), clusters, subclones) %>%
  tally() 


ggplot(data = gg_data,
       aes(axis1 = subclones, axis2 = clusters)) +
  geom_alluvium(aes(fill = subclones), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("subclone", "RNA_snn_res.0.4"), expand = c(.05, .05)) +
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


