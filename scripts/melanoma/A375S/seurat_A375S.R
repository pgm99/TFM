library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(biomaRt)
library(curl)
library(tidyverse)
library(ggplot2)
library(ggalluvial)
A375S <- read.table("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/data/GSM5022596_raw_counts_probe4.txt.gz", header = TRUE)
saveRDS(A375S_data, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/A375S.rds")
A375S_data <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375S/A375S.rds")


# Anotar hgnc symbol en vez de ensembl
ensemble_id <- row.names(A375S)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl)
listDatasets(ensembl)
listFilters(ensembl)
hgnc_symbol <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), values = ensemble_id, mart = ensembl, uniqueRows = TRUE, filters = "ensembl_gene_id")
hgnc_symbol <- arrange(hgnc_symbol, hgnc_symbol)
hgnc_symbol <- hgnc_symbol[!duplicated(hgnc_symbol$hgnc_symbol),]

nuevo_dataframe <- merge(hgnc_symbol, A375S, by.x = "ensembl_gene_id", by.y = "row.names", all.x = TRUE)
# Eliminar la columna adicional creada por el merge
nuevo_dataframe <- nuevo_dataframe[, -1]

colnames(nuevo_dataframe) <- c("hgnc_symbol", colnames(A375S))
nuevo_dataframe[nuevo_dataframe == ""] <- NA
A375S_final <- na.omit(nuevo_dataframe)
row.names(A375S_final) <- A375S_final$hgnc_symbol
A375S_final <- A375S_final[,-1]


#objeto Seurat
A375S_data <- CreateSeuratObject(counts= A375S_final)
A375S_data <- SetIdent(A375S_data, value = "A375S")
####QC####
A375S_data[["percent.mt"]] <- PercentageFeatureSet(A375S_data, pattern = "^MT-")
A375S_data[["percent.rb"]] <- PercentageFeatureSet(A375S_data, pattern = "^RP[SL]")
head(A375S_data@meta.data, 5)
VlnPlot(A375S_data, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), cols =3)

plot1 <- FeatureScatter(A375S_data, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 3)
plot2 <- FeatureScatter(A375S_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
#Eliminamos células con % mitocondrial elevado (células muertas) y con pocas o demasiadas cuentas(doublets)
A375S_data <- subset(A375S_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# assay <- GetAssay(A375S_data, assay = "RNA")
# assay <- as.data.frame(assay)
# write.table(assay, file= "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/data/GSM5022596_probe4_anotado.txt.gz", row.names = TRUE, col.names = TRUE)



#Normalización
A375S_data <- NormalizeData(A375S_data, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
A375S_data <- FindVariableFeatures(A375S_data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(A375S_data), 10)
plot1 <- VariableFeaturePlot(A375S_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Escalado linear prePCA
all.genes <- rownames(A375S_data)
A375S_data <- ScaleData(A375S_data, features = all.genes)

#PCA y formas de visualizarlo
A375S_data <- RunPCA(A375S_data, features = VariableFeatures(object = A375S_data))
#print(A375S_data[["pca"]], dims = 1:10, nfeatures = 5)
#VizDimLoadings(A375S_data, dims = 1:2, reduction = "pca")
#DimPlot(A375S_data, reduction = "pca")

#VizDimLoadings(A375S_data, dims = 1:5, reduction = "pca")
#DimPlot(A375S_data, reduction = "pca", dims = 1:2)

#Heatmaps
DimHeatmap(A375S_data, dims = 1:10, cells = 500, balanced = TRUE)

#Determinar dimensionalidad de PC con JackStraw y con ElbowPlot
A375S_data <- JackStraw(A375S_data, num.replicate = 1000)
A375S_data <- ScoreJackStraw(A375S_data, dims = 1:20)
JackStrawPlot(A375S_data, dims = 1:20)
ElbowPlot(A375S_data)

#Clustering
A375S_data <- FindNeighbors(A375S_data, dims = 1:15)
A375S_data <- FindClusters(A375S_data, resolution = c(0.4, 0.6, 0.8, 1.0,1.2, 1.4))

clustree(A375S_data)
Idents(object = A375S_data) <- "RNA_snn_res.0.6"

#UMAP
A375S_data <- RunUMAP(A375S_data, dims = 1:15, min.dist = 0, n.neighbors = 40)
DimPlot(A375S_data, reduction = "umap", label = TRUE, label.size = 5)

#UMAP para los clusters de BC
Idents(object = A375S_data) <- "bc_clusters_res.0.4"
A375S_data <- RunUMAP(A375S_data, dims = 1:15)
DimPlot(A375S_data, reduction = "umap", label = TRUE, label.size = 5)


#Cell cycling
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
A375S_data<- CellCycleScoring(A375S_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
A375S_data <- RunPCA(A375S_data, features = c(s.genes, g2m.genes))
Idents(object = A375S_data) <- "Phase"
DimPlot(A375S_data)

A375S_data <- ScaleData(A375S_data, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(A375S_data))



#Biomarkers
A375S_data.allmarkers <- FindAllMarkers(A375S_data, only.pos = TRUE)
# A375S_top_markers <- A375S_data.markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC)

DoHeatmap(A375S_data, features = A375S_data.markers$gene)

annotated_markers <- getBM(attributes = c("hgnc_symbol","description", "start_position", "end_position", "chromosome_name", "gene_biotype", "entrezgene_description"), 
                           values = A375S_data.markers$gene, mart = ensembl, uniqueRows = TRUE, filters = "hgnc_symbol")
annotated_markers <- annotated_markers[-c(44,45,46,47,48,49,50), ]
markers_list <- data.frame(A375S_data.allmarkers$gene, A375S_data.allmarkers$avg_log2FC)
write.csv(x = annotated_markers, file = "A375S_annotated_markers.csv", row.names = TRUE)
write.table(markers_list, file = "markers_list.rnk", row.names = F, quote = F, col.names = F, sep = "\t")


#Subclones SCEVAN
scevan <- readRDS(file = "/home/lpgonzalezm/TFM/A375S_scevan.rds")
scevan_subclone1 <- read.table(file = "/home/lpgonzalezm/output_old/A375S_subclone1_CN.seg", header = TRUE, sep = "\t")
scevan$subclone <- as.factor(scevan$subclone)
rownames(scevan) <- gsub(pattern = "-", replacement = ".", x = rownames(scevan))
scevan_filtered <- scevan[colnames(A375S_data), ]
A375S_data <- AddMetaData(object = A375S_data, metadata = scevan)

DimPlot(
  A375S_data,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  group.by = c("RNA_snn_res.0.6", "subclone")
#  shape.by = "subclone",
  )




#SANKEY
samples <- colnames(A375S_data)
clusters <- A375S_data@meta.data$RNA_snn_res.0.6
subclones <- A375S_data@meta.data$subclone
gg_data <- A375S_data@meta.data %>%
  group_by(colnames(A375S_data), clusters, subclones) %>%
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








library(SingleR)
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
library(scRNAseq)
library(scuttle)
sceM <- MuraroPancreasData()
sceM <- sceM[,!is.na(sceM$label)]
sceM <- logNormCounts(sceM)
rownames(sceM) <- sub("[_].*", "", rownames(sceM))


Idents(object = A375S_data) <- "RNA_snn_res.0.4"
A375S_data_celltype <- SingleR(GetAssayData(A375S_data, assay = "RNA", slot = "data"),
                            clusters = Idents(A375S_data), ref = hpca.se,  labels =hpca.se$label.main, de.method = "wilcox") 

A375S_data[["cell_type"]] <- A375S_data_celltype$labels[match(A375S_data[[]][["RNA_snn_res.0.2"]], rownames(A375S_data_celltype))]
Idents(object = A375S_data) <- "cell_type"
A375S_data <- RunUMAP(A375S_data, dims = 1:15, min.dist = 0.25, n.neighbors = 30)
cell_type_case1_YF_umap <- DimPlot(A375S_data, reduction = "umap", label = TRUE, label.size = 4, pt.size = 1.2)
annotated_celltype <- A375S_data$class
write.table(x = A375S_data$class, file = "/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/SComatic/A375S_annotated_celltype.tsv", sep = "\t")









