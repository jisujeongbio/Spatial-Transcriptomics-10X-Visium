library(Seurat)
library(dplyr)

Idents(brain.harmony) <- brain.harmony$region_cluster

common_groups <- c("WT", "G2_3")
deg_by_cluster <- list()

for (cluster in unique(brain.harmony$region_cluster)) {
  message("Running DEG analysis for cluster: ", cluster)
  
  cluster_subset <- subset(brain.harmony, subset = region_cluster == !!cluster)
  groups_present <- unique(cluster_subset$orig.ident)
  
  if (all(common_groups %in% groups_present)) {
    deg_res <- FindMarkers(
      cluster_subset,
      ident.1 = "G2_3",
      ident.2 = "WT",
      group.by = "orig.ident",
      test.use = "wilcox",
      logfc.threshold = 1,
      min.pct = 0.1,
      min.cells.group = 3,
      only.pos = FALSE,
      recorrect_umi = FALSE
    )
    
    deg_res$gene <- rownames(deg_res)
    deg_res$cluster <- cluster
    deg_by_cluster[[cluster]] <- deg_res
  }
}

deg_by_cluster_filter <- lapply(deg_by_cluster, function(df) {
  df %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC))
})

deg_combined <- bind_rows(deg_by_cluster_filter)

write.csv(deg_combined, "DEG_all_clusters.csv", row.names = FALSE)

deg_count_df <- data.frame(
  Cluster = names(deg_by_cluster_filter),
  DEG_Count = sapply(deg_by_cluster_filter, nrow)
)

write.csv(deg_count_df, "DEG_count_by_cluster.csv", row.names = FALSE)