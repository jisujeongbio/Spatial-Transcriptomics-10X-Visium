# 02_enrichment_analysis.R
# Perform GO Biological Process and KEGG enrichment analysis

library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gprofiler2)
library(openxlsx)

# 1. GO Biological Process enrichment
all_go_results <- list()

for (cluster in names(deg_by_cluster_filter)) {
  message("Running GO analysis for cluster: ", cluster)
  
  gene_symbols <- deg_by_cluster_filter[[cluster]]$gene
  
  gene_entrez <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  entrez_ids <- unique(gene_entrez$ENTREZID)
  
  if (length(entrez_ids) < 10) {
    message("Skipping ", cluster, ": too few mapped genes")
    next
  }
  
  go_result <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  go_df <- as.data.frame(go_result)
  
  if (nrow(go_df) == 0) next
  
  go_df$cluster <- cluster
  all_go_results[[cluster]] <- go_df
}

all_go_df <- bind_rows(all_go_results)
write.xlsx(all_go_df, "all_GO_BP_results.xlsx", rowNames = FALSE)

# 2. KEGG enrichment using g:Profiler
all_kegg_results <- list()

for (cluster in names(deg_by_cluster_filter)) {
  message("Running KEGG analysis for cluster: ", cluster)
  
  gene_symbols <- deg_by_cluster_filter[[cluster]]$gene
  
  if (length(gene_symbols) < 10) {
    message("Skipping ", cluster, ": too few genes")
    next
  }
  
  gostres <- gost(
    query = gene_symbols,
    organism = "mmusculus",
    sources = "KEGG",
    correction_method = "g_SCS",
    evcodes = TRUE,
    significant = TRUE
  )
  
  if (is.null(gostres$result) || nrow(gostres$result) == 0) next
  
  kegg_df <- gostres$result
  kegg_df$cluster <- cluster
  all_kegg_results[[cluster]] <- kegg_df
}

all_kegg_df <- bind_rows(all_kegg_results)
write.xlsx(all_kegg_df, "all_KEGG_results.xlsx", rowNames = FALSE)