# 00_region_annotation_summary.R
# Purpose:
# Summarize region annotation workflow for integrated spatial transcriptomics data.
# This script assigns anatomical region labels to clusters using:
# 1) region-specific marker reference,
# 2) overlap between reference and cluster markers,
# 3) anatomical context.

library(Seurat)
library(dplyr)
library(readxl)

# 1. Load region reference markers
region_ref <- read_excel("data/Region_type_marker.xlsx") |> 
  as.data.frame()

region_ref_list <- split(region_ref$genes, region_ref$region)

region_ref_list <- lapply(region_ref_list, function(x) {
  genes <- unlist(strsplit(unlist(x), ",\\s*"))
  unique(genes)
})

genes_in_object <- rownames(brain.harmony)

region_ref_list_filtered <- lapply(region_ref_list, function(x) {
  x[x %in% genes_in_object]
})

# use top markers only
region_ref_list_top50 <- lapply(region_ref_list_filtered, function(x) {
  head(x, 50)
})

# 2. Identify cluster markers
Idents(brain.harmony) <- brain.harmony$seurat_clusters

cluster_markers <- FindAllMarkers(
  brain.harmony,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

cluster_marker_list <- split(cluster_markers$gene, cluster_markers$cluster)

cluster_marker_list_top50 <- lapply(cluster_marker_list, function(x) {
  unique(head(x, 50))
})

# 3. Define overlap score function
calc_overlap_jaccard <- function(ref_genes, cluster_genes) {
  ref_genes <- na.omit(unique(ref_genes))
  cluster_genes <- na.omit(unique(cluster_genes))
  
  inter <- intersect(ref_genes, cluster_genes)
  union_set <- union(ref_genes, cluster_genes)
  
  data.frame(
    overlap_n = length(inter),
    jaccard = ifelse(length(union_set) == 0, NA, length(inter) / length(union_set))
  )
}

# 4. Compare each cluster against each region reference
annotation_scores <- lapply(names(cluster_marker_list_top50), function(cl) {
  cluster_genes <- cluster_marker_list_top50[[cl]]
  
  res <- lapply(names(region_ref_list_top50), function(region) {
    score <- calc_overlap_jaccard(region_ref_list_top50[[region]], cluster_genes)
    score$cluster <- cl
    score$region <- region
    score
  })
  
  bind_rows(res)
}) |> bind_rows()

# 5. Assign best-matching region per cluster
cluster_region_assignment <- annotation_scores |>
  group_by(cluster) |>
  arrange(desc(overlap_n), desc(jaccard), .by_group = TRUE) |>
  slice(1) |>
  ungroup()

cluster_region_assignment

# 6. Add region labels to Seurat object
cluster_to_region <- setNames(
  cluster_region_assignment$region,
  cluster_region_assignment$cluster
)

brain.harmony$region_cluster <- cluster_to_region[
  as.character(brain.harmony$seurat_clusters)
]

write.csv(annotation_scores, "results/region_annotation_scores.csv", row.names = FALSE)
write.csv(cluster_region_assignment, "results/cluster_region_assignment.csv", row.names = FALSE)
