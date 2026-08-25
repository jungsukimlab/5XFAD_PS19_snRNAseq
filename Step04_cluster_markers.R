# 5XFAD;PS19 snRNA-seq analysis
# Step 4: Cluster marker identification
#
# Run from the project root directory.
#
# Input:
#   ./publication_output/Step3_clustering/5XFAD_PS19_Harmony_unannotated.rds
#
# Output:
#   ./publication_output/Step4_cluster_markers/
#
# Cluster markers identified here are used for cell-type annotation.
# This is distinct from downstream genotype differential expression analysis.

library(Seurat)

# -----------------------------------------------------------------------------
# Paths and marker-analysis parameters
# -----------------------------------------------------------------------------

input_file <- file.path(
  "publication_output",
  "Step3_clustering",
  "5XFAD_PS19_Harmony_unannotated.rds"
)

output_dir <- file.path(
  "publication_output",
  "Step4_cluster_markers"
)

# Prevent accidental replacement of existing publication outputs.
# Change to TRUE only when intentional overwriting is desired.
overwrite <- FALSE

cluster_column <- "SCT_snn_res.0.6"

min_pct <- 0.25
logfc_threshold <- 0.25

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

check_output <- function(file) {
  if (file.exists(file) && !overwrite) {
    stop(
      "Output already exists: ", file,
      "\nSet overwrite <- TRUE only if you intend to replace it."
    )
  }
}

# -----------------------------------------------------------------------------
# Load clustered object
# -----------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

obj <- readRDS(input_file)

DefaultAssay(obj) <- "SCT"

if (!cluster_column %in% colnames(obj[[]])) {
  stop(
    "Cluster metadata column not found: ",
    cluster_column
  )
}

Idents(obj) <- cluster_column

message(
  "Identifying markers for ",
  length(levels(Idents(obj))),
  " clusters."
)

# -----------------------------------------------------------------------------
# Prepare SCT assay for marker testing
# -----------------------------------------------------------------------------

# The merged object contains multiple SCT models. PrepSCTFindMarkers()
# recorrects the SCT counts using a common sequencing-depth covariate so that
# marker testing can be performed across cells originating from different
# sample-specific SCT models.

obj <- PrepSCTFindMarkers(
  object = obj,
  assay = "SCT",
  verbose = TRUE
)

# -----------------------------------------------------------------------------
# Identify positive markers for each cluster
# -----------------------------------------------------------------------------

# Wilcoxon rank-sum testing is used here to identify cluster-enriched markers
# for cell-type annotation. This marker analysis is separate from downstream
# genotype differential expression analysis using MAST.

cluster_markers <- FindAllMarkers(
  object = obj,
  assay = "SCT",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold,
  verbose = TRUE
)

# -----------------------------------------------------------------------------
# Save marker table and prepared Seurat object
# -----------------------------------------------------------------------------

markers_file <- file.path(
  output_dir,
  "5XFAD_PS19_res0.6_cluster_markers.csv"
)

prepared_object_file <- file.path(
  output_dir,
  "5XFAD_PS19_cluster_markers_prepared.rds"
)

check_output(markers_file)
check_output(prepared_object_file)

write.csv(
  cluster_markers,
  markers_file,
  row.names = FALSE
)

saveRDS(
  obj,
  prepared_object_file
)
