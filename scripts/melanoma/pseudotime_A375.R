library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(beyondcell)
A375_merged <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/seurat_A375.rds")
cds <- as.cell_data_set(A375_merged)
bc_recomputed_global <- readRDS("/home/lpgonzalezm/TFM/beyondcell/prueba_melanoma/A375_merged/bc_global_0.1.rds")
cds_bc <- as.cell_data_set(bc_recomputed_global)


recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
Idents(object = A375_data) <- "RNA_snn_res.0.4"
Idents(object = A375_merged) <- "bc_clusters_res.0.4"

list.cluster <- A375_merged@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- A375_merged@reductions$umap@cell.embeddings

cluster.before.traj <-plot_cells(cds, color_cells_by = "type", graph_label_size = 3, cell_size = 0.8, label_groups_by_cluster = F, 
                                 group_label_size = 6) + theme(legend.position = "right")
cluster.before.traj

cds <- learn_graph(cds, use_partition = F)

plot_cells(cds, color_cells_by = "type", cell_size = 0.8,
           graph_label_size = 3,
           label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = T, label_cell_groups = F,
           group_label_size = 7)

cds_pt <- order_cells(cds, reduction_method = "UMAP")

plot_cells(cds_pt, color_cells_by = "pseudotime", label_groups_by_cluster = F,
           cell_size = 0.8, group_label_size = 7,
           graph_label_size = 3, label_branch_points = T,
           label_roots = T, label_leaves = T)

cds_pt$monocle3_pseudotime <- pseudotime(cds_pt)
data.pseudo <- as.data.frame(colData(cds_pt))

ggplot(data.pseudo, aes(monocle3_pseudotime, RNA_snn_res.0.6, fill = RNA_snn_res.0.6)) + geom_boxplot()


deg <- graph_test(cds_pt, neighbor_graph = "principal_graph")
deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()

# FeaturePlot(A375_data, features = c("EPS8L3", "RASSF9", "LRRTM3", "PLXNB3"))


# a helper function to identify the root principal points:


get_earliest_principal_node <- function(cds, type="A375S"){
  cell_ids <- which(colData(cds)[, "type"] == type)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=T,
           label_branch_points=F,
           graph_label_size=4,
           cell_size = 0.8
           )











