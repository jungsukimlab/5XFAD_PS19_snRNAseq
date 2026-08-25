# 5XFAD;PS19 snRNA-seq analysis
# Step 3: Clustering and resolution selection
#
# Run from the project root directory.
#
# Input:
#   ./publication_output/Step2_integration/SCTintegrated_Harmony.rds
#
# Output:
#   ./publication_output/Step3_clustering/
#
# The SNN graph used here was constructed from Harmony embeddings in Step 2.
# Multiple clustering resolutions are evaluated with clustree, and resolution
# 0.6 is retained for downstream annotation.

library(Seurat)
library(clustree)

set.seed(4444)

# -----------------------------------------------------------------------------
# Paths and clustering parameters
# -----------------------------------------------------------------------------

input_file <- file.path(
  "publication_output",
  "Step2_integration",
  "SCTintegrated_Harmony.rds"
)

output_dir <- file.path(
  "publication_output",
  "Step3_clustering"
)

plot_dir <- file.path(
  output_dir,
  "plots"
)

# Prevent accidental replacement of existing publication outputs.
# Change to TRUE only when intentional overwriting is desired.
overwrite <- FALSE

cluster_resolutions <- c(0.1, 0.4, 0.6, 0.8, 1.2)
selected_resolution <- 0.6

# Graph generated in Step 2 by FindNeighbors(..., reduction = "harmony").
snn_graph <- "SCT_snn"

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
# Load Harmony-integrated object
# -----------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

obj <- readRDS(input_file)

DefaultAssay(obj) <- "SCT"

if (!snn_graph %in% Graphs(obj)) {
  stop(
    "Expected SNN graph '", snn_graph, "' was not found. ",
    "Available graphs: ", paste(Graphs(obj), collapse = ", ")
  )
}

if (!"umap" %in% Reductions(obj)) {
  stop(
    "Expected UMAP reduction 'umap' was not found. ",
    "Available reductions: ", paste(Reductions(obj), collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Evaluate clustering across multiple resolutions
# -----------------------------------------------------------------------------

obj <- FindClusters(
  object = obj,
  graph.name = snn_graph,
  resolution = cluster_resolutions
)

resolution_columns <- paste0(
  snn_graph,
  "_res.",
  cluster_resolutions
)

missing_columns <- setdiff(
  resolution_columns,
  colnames(obj[[]])
)

if (length(missing_columns) > 0) {
  stop(
    "Expected clustering metadata columns were not created: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Summarize the number of clusters produced at each tested resolution.
resolution_summary <- data.frame(
  resolution = cluster_resolutions,
  n_clusters = vapply(
    resolution_columns,
    function(column) {
      length(unique(obj[[]][[column]]))
    },
    integer(1)
  )
)

resolution_summary_file <- file.path(
  output_dir,
  "Step3_resolution_summary.csv"
)

check_output(resolution_summary_file)

write.csv(
  resolution_summary,
  resolution_summary_file,
  row.names = FALSE
)

# Save a checkpoint containing all tested clustering resolutions.
before_clustree_file <- file.path(
  output_dir,
  "Step3_before_resolution_selection.rds"
)

check_output(before_clustree_file)

saveRDS(
  obj,
  before_clustree_file
)

# -----------------------------------------------------------------------------
# Clustering tree
# -----------------------------------------------------------------------------

clustree_file <- file.path(
  plot_dir,
  "Step3_clustree.pdf"
)

check_output(clustree_file)

pdf(clustree_file)

print(
  clustree::clustree(
    obj,
    prefix = paste0(snn_graph, "_res.")
  )
)

dev.off()

# -----------------------------------------------------------------------------
# Select final clustering resolution
# -----------------------------------------------------------------------------

selected_cluster_column <- paste0(
  snn_graph,
  "_res.",
  selected_resolution
)

if (!selected_cluster_column %in% colnames(obj[[]])) {
  stop(
    "Selected clustering column not found: ",
    selected_cluster_column
  )
}

# Set resolution 0.6 as the active cluster identity.
Idents(obj) <- selected_cluster_column

# Keep Seurat's standard cluster metadata synchronized with the selected
# resolution to avoid ambiguity in downstream analyses.
obj$seurat_clusters <- obj[[]][[selected_cluster_column]]

# -----------------------------------------------------------------------------
# UMAP of selected clusters
# -----------------------------------------------------------------------------

umap_cluster_file <- file.path(
  plot_dir,
  paste0(
    "Step3_UMAP_Harmony_res",
    selected_resolution,
    ".pdf"
  )
)

check_output(umap_cluster_file)

pdf(umap_cluster_file)

print(
  DimPlot(
    obj,
    reduction = "umap",
    label = TRUE,
    raster = FALSE
  )
)

dev.off()

# -----------------------------------------------------------------------------
# Check cluster representation across biological replicates
# -----------------------------------------------------------------------------

if (!"ID" %in% colnames(obj[[]])) {
  stop("Sample metadata column 'ID' was not found.")
}

cluster_sample_table <- table(
  Cluster = obj[[]][[selected_cluster_column]],
  Sample = obj$ID
)

cluster_sample_file <- file.path(
  output_dir,
  "Step3_cluster_by_sample_counts.csv"
)

check_output(cluster_sample_file)

write.csv(
  as.data.frame.matrix(cluster_sample_table),
  cluster_sample_file,
  row.names = TRUE
)

print(cluster_sample_table)

# -----------------------------------------------------------------------------
# Remove clustering resolutions not selected for downstream analysis
# -----------------------------------------------------------------------------

columns_to_remove <- setdiff(
  resolution_columns,
  selected_cluster_column
)

for (column in columns_to_remove) {
  obj[[column]] <- NULL
}

# Reset identities after metadata cleanup.
Idents(obj) <- selected_cluster_column

# -----------------------------------------------------------------------------
# Save final unannotated clustered object
# -----------------------------------------------------------------------------

clustered_file <- file.path(
  output_dir,
  "5XFAD_PS19_Harmony_unannotated.rds"
)

check_output(clustered_file)

saveRDS(
  obj,
  clustered_file
)
