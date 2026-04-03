# 03_CARD_deconvolution.R
# Region-matched CARD deconvolution using multiple ABA/WMB-10X references

library(Seurat)
library(CARD)
library(zellkonverter)
library(dplyr)
library(SummarizedExperiment)

prepare_reference_expr <- function(se_obj) {
  expr <- assay(se_obj, "X")
  gene_symbol <- rowData(se_obj)$gene_symbol
  
  valid <- !is.na(gene_symbol) & gene_symbol != ""
  expr <- expr[valid, ]
  gene_symbol <- gene_symbol[valid]
  
  rownames(expr) <- gene_symbol
  expr <- expr[!duplicated(rownames(expr)), ]
  
  expr
}

run_card_for_reference <- function(ref_path, ref_meta, spatial_count, spatial_location) {
  ref_obj <- zellkonverter::readH5AD(ref_path, reader = "R", verbose = TRUE)
  ref_expr <- prepare_reference_expr(ref_obj)
  
  ref_meta$custom <- recode_celltype(ref_meta$subclass)
  
  card_obj <- createCARDObject(
    sc_count = ref_expr,
    sc_meta = ref_meta,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    ct.varname = "custom",
    ct.select = unique(ref_meta$custom),
    sample.varname = NULL,
    minCountGene = 100,
    minCountSpot = 5
  )
  
  CARD_deconvolution(CARD_object = card_obj)
}

extract_region_proportion <- function(card_obj, coord_df) {
  prop <- as.data.frame(card_obj@Proportion_CARD)
  prop[rownames(prop) %in% coord_df$Coordinate, , drop = FALSE]
}

WMB10Xv3_metadata <- read.csv("REFERENCE_METADATA.csv")
WMB10Xv3_metadata <- WMB10Xv3_metadata[WMB10Xv3_metadata$library_method == "10Xv3", ]

reference_config <- list(
  HPF = list(
    ref_path = "WMB-10Xv3-HPF-log2.h5ad",
    meta = WMB10Xv3_metadata[
      WMB10Xv3_metadata$anatomical_division_label == "HPF", ],
    region_labels = c("CA1_CA2", "DG-mo", "DG-sg", "CA3", "SUB"),
    swap_coord = TRUE
  ),
  OLF = list(
    ref_path = "WMB-10Xv3-OLF-log2.h5ad",
    meta = WMB10Xv3_metadata[
      WMB10Xv3_metadata$anatomical_division_label == "OLF", ],
    region_labels = c("OLF"),
    swap_coord = FALSE
  ),
  P = list(
    ref_path = "WMB-10Xv3-P-log2.h5ad",
    meta = WMB10Xv3_metadata[
      WMB10Xv3_metadata$anatomical_division_label == "P", ],
    region_labels = c("P"),
    swap_coord = FALSE
  ),
  TH = list(
    ref_path = "WMB-10Xv3-TH-log2.h5ad",
    meta = WMB10Xv3_metadata[
      WMB10Xv3_metadata$anatomical_division_label == "TH", ],
    region_labels = c("TH"),
    swap_coord = FALSE
  ),
  Isocortex1 = list(
    ref_path = "WMB-10Xv3-Isocortex-1-log2.h5ad",
    meta = WMB10Xv3_metadata[
      WMB10Xv3_metadata$matrix_label == "WMB-10Xv3-Isocortex-1", ],
    region_labels = c("Isocortex1"),
    swap_coord = FALSE
  ),
  Isocortex2 = list(
    ref_path = "WMB-10Xv3-Isocortex-2-log2.h5ad",
    meta = WMB10Xv3_metadata[
      WMB10Xv3_metadata$matrix_label == "WMB-10Xv3-Isocortex-2", ],
    region_labels = c("Isocortex2"),
    swap_coord = FALSE
  ),
  MB = list(
    ref_path = "WMB-10Xv3-MB-log2.h5ad",
    meta = WMB10Xv3_metadata[
      WMB10Xv3_metadata$anatomical_division_label == "MB", ],
    region_labels = c("MB", "SN"),
    swap_coord = FALSE
  )
)

card_results <- list()
region_proportion_list <- list()

for (ref_name in names(reference_config)) {
  message("Running CARD for: ", ref_name)
  
  cfg <- reference_config[[ref_name]]
  
  card_obj <- run_card_for_reference(
    ref_path = cfg$ref_path,
    ref_meta = cfg$meta,
    spatial_count = spatial_counts_data,
    spatial_location = spatial_coords
  )
  
  prop <- as.data.frame(card_obj@Proportion_CARD)
  
  if (isTRUE(cfg$swap_coord)) {
    rownames(prop) <- swap_xy(rownames(prop))
  }
  
  coord_sub <- coord_barcode[coord_barcode$brain_name %in% cfg$region_labels, , drop = FALSE]
  prop_real <- prop[rownames(prop) %in% coord_sub$Coordinate, , drop = FALSE]
  
  card_results[[ref_name]] <- card_obj
  region_proportion_list[[ref_name]] <- prop_real
}

merged_proportion <- bind_rows(region_proportion_list)
merged_proportion[is.na(merged_proportion)] <- 0

matched_indices <- match(rownames(merged_proportion), coord_barcode$Coordinate)
merged_proportion$brain_region <- coord_barcode$brain_name[matched_indices]
merged_proportion$ident <- coord_barcode$ident[matched_indices]

write.csv(merged_proportion, "merged_CARD_proportions.csv", row.names = TRUE)
saveRDS(card_results, file = "CARD_objects_by_reference.rds")