# 5XFAD;PS19 snRNA-seq analysis
# Step 5: Cell-type annotation and marker validation
#
# Run from the project root directory.
#
# Inputs:
#   ./publication_output/Step4_cluster_markers/
#       5XFAD_PS19_res0.6_cluster_markers.csv
#       5XFAD_PS19_cluster_markers_prepared.rds
#
# Outputs:
#   ./publication_output/Step5_annotation/
#
# Cell-type assignments were informed by SAHA (https://github.com/neurogenetics/SAHA) marker-based annotation and
# manual inspection of canonical lineage markers.
#
# Abbreviations:
#   ExN = excitatory neuron
#   InN = inhibitory neuron
#   O   = oligodendrocyte
#   M   = microglia
#   A   = astrocyte
#   OPC = oligodendrocyte precursor cell
#   EC  = endothelial cell
#   FB  = fibroblast
#
# R and package versions are documented at the project level (renv.lock).

library(Seurat)
library(SAHA)
library(SAHAdata)
library(ggplot2)

# -----------------------------------------------------------------------------
# Paths and annotation parameters
# -----------------------------------------------------------------------------

input_dir <- file.path(
  "publication_output",
  "Step4_cluster_markers"
)

marker_file <- file.path(
  input_dir,
  "5XFAD_PS19_res0.6_cluster_markers.csv"
)

input_object_file <- file.path(
  input_dir,
  "5XFAD_PS19_cluster_markers_prepared.rds"
)

output_dir <- file.path(
  "publication_output",
  "Step5_annotation"
)

plot_dir <- file.path(
  output_dir,
  "plots"
)

# Prevent accidental replacement of existing publication outputs.
# Change to TRUE only when intentional overwriting is desired.
overwrite <- FALSE

cluster_column <- "SCT_snn_res.0.6"
umap_reduction <- "umap"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

check_output <- function(file) {
  if (file.exists(file) && !overwrite) {
    stop(
      "Output already exists: ", file,
      "\nSet overwrite <- TRUE only if you intend to replace it."
    )
  }
}


# -----------------------------------------------------------------------------
# Define output files
# -----------------------------------------------------------------------------

saha_self_similarity_file <- file.path(
  plot_dir,
  "Step5_SAHA_self_similarity.pdf"
)

saha_marker_file <- file.path(
  plot_dir,
  "Step5_SAHA_marker_based.pdf"
)

saha_object_file <- file.path(
  output_dir,
  "Step5_SAHA_annotation.rds"
)

manual_qc_file <- file.path(
  plot_dir,
  "Step5_manual_marker_sanity_checks.pdf"
)

annotation_table_file <- file.path(
  output_dir,
  "Step5_cluster_annotation_table.csv"
)

cluster_umap_file <- file.path(
  plot_dir,
  "Step5_UMAP_annotated_clusters.pdf"
)

cell_type_umap_file <- file.path(
  plot_dir,
  "Step5_UMAP_cell_types.pdf"
)

dotplot_file <- file.path(
  plot_dir,
  "Step5_cell_type_marker_dotplot.pdf"
)

representative_vln_file <- file.path(
  plot_dir,
  "Step5_representative_cluster_markers.pdf"
)

cell_type_vln_file <- file.path(
  plot_dir,
  "Step5_cell_type_markers.pdf"
)

annotated_object_file <- file.path(
  output_dir,
  "5XFAD_PS19_Harmony_annotated.rds"
)

output_files <- c(
  saha_self_similarity_file,
  saha_marker_file,
  saha_object_file,
  manual_qc_file,
  annotation_table_file,
  cluster_umap_file,
  cell_type_umap_file,
  dotplot_file,
  representative_vln_file,
  cell_type_vln_file,
  annotated_object_file
)

invisible(lapply(output_files, check_output))

# -----------------------------------------------------------------------------
# Check inputs
# -----------------------------------------------------------------------------

if (!file.exists(marker_file)) {
  stop("Marker file not found: ", marker_file)
}

if (!file.exists(input_object_file)) {
  stop("Input Seurat object not found: ", input_object_file)
}

# =============================================================================
# SAHA marker-based annotation
# =============================================================================

data("ISOCTX_Markers", package = "SAHAdata")
data("ABC_meta", package = "SAHAdata")

cluster_markers <- read.csv(
  marker_file,
  stringsAsFactors = FALSE
)

saha_meta <- ABC_meta

ann <- Create_SAHA_object(
  query = cluster_markers,
  db = ISOCTX_Markers,
  data_type = "Markers"
)

ann <- Initialize_Markers(ann)

ann <- Initialize_Self_Similiarity(
  ann,
  slot = "Markers"
)

ann <- Create_SelfSimilarity_Viz(
  ann,
  slot = "Markers"
)

pdf(saha_self_similarity_file)

call_SAHA_plots(
  ann,
  plot_type = "self-similarity",
  data_type = "Markers"
)

dev.off()

ann <- Run_Marker_Based(ann)

ann <- Create_MarkerBased_Viz(
  ann,
  meta = saha_meta
)

pdf(saha_marker_file)

call_SAHA_plots(
  ann,
  plot_type = "Marker-based",
  data_type = "Markers"
)

dev.off()

saveRDS(
  ann,
  saha_object_file
)

# =============================================================================
# Manual annotation and marker validation
# =============================================================================

obj <- readRDS(input_object_file)

DefaultAssay(obj) <- "SCT"

if (!cluster_column %in% colnames(obj[[]])) {
  stop(
    "Cluster metadata column not found: ",
    cluster_column
  )
}

if (!umap_reduction %in% Reductions(obj)) {
  stop(
    "UMAP reduction '", umap_reduction, "' was not found. ",
    "Available reductions: ",
    paste(Reductions(obj), collapse = ", ")
  )
}

Idents(obj) <- cluster_column

# -----------------------------------------------------------------------------
# Manual marker sanity checks
# -----------------------------------------------------------------------------

# Neuronal markers
neuronal_markers <- c(
  "Slc17a7", "Slc17a6", "Gad1", "Gad2"
)

# Astrocyte markers
astrocyte_markers <- c(
  "Slc1a3", "Aqp4", "Gja1", "Gfap", "Aldoc", "F3"
)

# Oligodendrocyte markers
oligodendrocyte_markers <- c(
  "Mog", "Mbp"
)

# Microglial markers
microglia_markers <- c(
  "P2ry12", "Inpp5d", "Csf1r", "Hexb"
)

# Endothelial markers
endothelial_markers <- c(
  "Flt1", "Id1", "Foxf2", "Foxq1"
)

marker_panels <- list(
  Neuronal = neuronal_markers,
  Astrocyte = astrocyte_markers,
  Oligodendrocyte = oligodendrocyte_markers,
  Microglia = microglia_markers,
  Endothelial = endothelial_markers
)

markers_to_check <- unique(unlist(marker_panels))


# -----------------------------------------------------------------------------
# Plot canonical markers across clusters
# -----------------------------------------------------------------------------

cluster_ids <- levels(Idents(obj))

if (all(grepl("^[0-9]+$", cluster_ids))) {
  cluster_ids <- as.character(sort(as.integer(cluster_ids)))
}

cluster_groups <- split(
  cluster_ids,
  ceiling(seq_along(cluster_ids) / 10)
)

pdf(
  manual_qc_file,
  width = 12,
  height = 7
)

for (cluster_group in cluster_groups) {
  for (panel_name in names(marker_panels)) {
    print(
      VlnPlot(
        obj,
        features = marker_panels[[panel_name]],
        idents = cluster_group,
        raster = FALSE
      ) +
        NoLegend() +
        ggtitle(
          paste0(
            panel_name,
            " markers | clusters ",
            paste(cluster_group, collapse = ", ")
          )
        )
    )
  }
}

dev.off()

# =============================================================================
# Final annotation maps
# =============================================================================

cell_type_map <- c(
  "0"  = "InN",
  "1"  = "O",
  "2"  = "M",
  "3"  = "ExN",
  "4"  = "ExN",
  "5"  = "ExN",
  "6"  = "InN",
  "7"  = "A",
  "8"  = "ExN",
  "9"  = "InN",
  "10" = "O",
  "11" = "ExN",
  "12" = "InN",
  "13" = "InN",
  "14" = "OPC",
  "15" = "InN",
  "16" = "ExN",
  "17" = "ExN",
  "18" = "InN",
  "19" = "M",
  "20" = "InN",
  "21" = "EC",
  "22" = "ExN",
  "23" = "ExN",
  "24" = "InN",
  "25" = "InN",
  "26" = "ExN",
  "27" = "ExN",
  "28" = "InN",
  "29" = "ExN",
  "30" = "A",
  "31" = "A",
  "32" = "ExN",
  "33" = "InN",
  "34" = "O",
  "35" = "FB",
  "36" = "EC",
  "37" = "InN"
)

clust_cell_map <- c(
  "0"  = "InN_1",
  "1"  = "O_1",
  "2"  = "M_1",
  "3"  = "ExN_1",
  "4"  = "ExN_2",
  "5"  = "ExN_3",
  "6"  = "InN_2",
  "7"  = "A_1",
  "8"  = "ExN_4",
  "9"  = "InN_3",
  "10" = "O_2",
  "11" = "ExN_5",
  "12" = "InN_4",
  "13" = "InN_5",
  "14" = "OPC",
  "15" = "InN_6",
  "16" = "ExN_6",
  "17" = "ExN_7",
  "18" = "InN_7",
  "19" = "M_2",
  "20" = "InN_8",
  "21" = "EC_1",
  "22" = "ExN_8",
  "23" = "ExN_9",
  "24" = "InN_9",
  "25" = "InN_10",
  "26" = "ExN_10",
  "27" = "ExN_11",
  "28" = "InN_11",
  "29" = "ExN_12",
  "30" = "A_2",
  "31" = "A_3",
  "32" = "ExN_13",
  "33" = "InN_12",
  "34" = "O_3",
  "35" = "FB",
  "36" = "EC_2",
  "37" = "InN_13"
)

cluster_ids <- as.character(
  obj[[]][[cluster_column]]
)

observed_clusters <- sort(
  unique(cluster_ids)
)


obj$cell_type <- unname(
  cell_type_map[cluster_ids]
)

obj$clust_cell <- unname(
  clust_cell_map[cluster_ids]
)

annotation_table <- data.frame(
  cluster = names(cell_type_map),
  cell_type = unname(cell_type_map),
  clust_cell = unname(clust_cell_map[names(cell_type_map)]),
  stringsAsFactors = FALSE
)

annotation_table$cluster <- as.integer(
  annotation_table$cluster
)

annotation_table <- annotation_table[
  order(annotation_table$cluster),
]

write.csv(
  annotation_table,
  annotation_table_file,
  row.names = FALSE
)

# =============================================================================
# Annotation visualization and validation
# =============================================================================

# -----------------------------------------------------------------------------
# Cluster-level annotation UMAP
# -----------------------------------------------------------------------------

Idents(obj) <- "clust_cell"

pdf(
  cluster_umap_file,
  width = 10,
  height = 8
)

print(
  DimPlot(
    obj,
    reduction = umap_reduction,
    label = TRUE,
    pt.size = 0.5,
    repel = TRUE,
    raster = FALSE,
    label.size = 4
  ) +
    NoLegend()
)

dev.off()

# -----------------------------------------------------------------------------
# Broad cell-type UMAP
# -----------------------------------------------------------------------------

Idents(obj) <- "cell_type"

pdf(
  cell_type_umap_file,
  width = 10,
  height = 8
)

print(
  DimPlot(
    obj,
    reduction = umap_reduction,
    label = TRUE,
    pt.size = 0.5,
    repel = TRUE,
    raster = FALSE
  ) +
    NoLegend()
)

dev.off()

# -----------------------------------------------------------------------------
# Canonical marker dot plot by broad cell type
# -----------------------------------------------------------------------------

annotation_markers <- c(
  "Gad1", "Gad2",
  "Mog", "Mbp",
  "P2ry12", "Inpp5d",
  "Slc17a7", "Satb2",
  "Gja1", "Slc1a3",
  "Olig1", "Vcan",
  "Slc6a20a", "Flt1", "Id1"
)

pdf(
  dotplot_file,
  width = 11,
  height = 6
)

print(
  DotPlot(
    obj,
    features = annotation_markers,
    group.by = "cell_type"
  ) +
    ggtitle("5XFAD;PS19 cell-type markers") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 9
      ),
      axis.text.y = element_text(size = 11)
    )
)

dev.off()

# -----------------------------------------------------------------------------
# Stacked violin plot using one representative cluster per broad cell type
# -----------------------------------------------------------------------------

representative_features <- c(
  "Slc17a7", "Satb2",
  "Gad1", "Gad2",
  "Mag", "Mog",
  "P2ry12", "Csf1r",
  "Aqp4", "Aldoc",
  "Vcan", "Slc6a20a",
  "Id1"
)

representative_clusters <- c(
  "ExN_1",
  "InN_1",
  "O_1",
  "M_1",
  "A_1",
  "OPC",
  "EC_1",
  "FB"
)

Idents(obj) <- "clust_cell"

pdf(
  representative_vln_file,
  width = 8,
  height = 9
)

print(
  VlnPlot(
    obj,
    features = representative_features,
    stack = TRUE,
    sort = FALSE,
    flip = TRUE,
    idents = representative_clusters
  ) +
    NoLegend() +
    theme(
      text = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    )
)

dev.off()

# -----------------------------------------------------------------------------
# Stacked violin plot by broad cell type
# -----------------------------------------------------------------------------

Idents(obj) <- "cell_type"

pdf(
  cell_type_vln_file,
  width = 8,
  height = 9
)

print(
  VlnPlot(
    obj,
    features = representative_features,
    stack = TRUE,
    sort = FALSE,
    flip = TRUE,
    idents = cell_type_order
  ) +
    NoLegend() +
    theme(
      text = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 18)
    )
)

dev.off()

# -----------------------------------------------------------------------------
# Save final annotated object
# -----------------------------------------------------------------------------

Idents(obj) <- "clust_cell"

saveRDS(
  obj,
  annotated_object_file
)
