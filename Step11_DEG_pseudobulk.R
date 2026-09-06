# =============================================================================
# 5XFAD;PS19 snRNA-seq analysis
# Step 11: Pseudobulk differential expression
# =============================================================================
#
# Run from the project root directory.
#
# Input:
#   publication_output/Step5_annotation/
#       5XFAD_PS19_Harmony_annotated.rds
#
# Output:
#   publication_output/Step11_DEG_pseudobulk/
#       pseudo_ID_geno_DEG.csv
#
# Analysis:
#   1. Aggregate SCT counts by biological sample
#   2. Differential expression using DESeq2
#
# This pseudobulk analysis provides the mouse-side input for Step 12.

library(Seurat)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

input_file <- file.path(
  "publication_output",
  "Step5_annotation",
  "5XFAD_PS19_Harmony_annotated.rds"
)

output_dir <- file.path(
  "publication_output",
  "Step11_DEG_pseudobulk"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_file <- file.path(
  output_dir,
  "pseudo_ID_geno_DEG.csv"
)

overwrite <- FALSE

check_output <- function(file) {
  if (file.exists(file) && !overwrite) {
    stop(
      "Output already exists: ", file,
      "\nSet overwrite <- TRUE only if you intend to replace it."
    )
  }
}

if (!file.exists(input_file)) {
  stop(
    "Input Seurat object not found: ",
    input_file
  )
}

check_output(
  output_file
)

# =============================================================================
# 1. Aggregate SCT counts by biological sample
# =============================================================================

obj <- readRDS(
  input_file
)

required_metadata <- c(
  "ID",
  "Genotype"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(obj[[]])
)

if (length(missing_metadata) > 0) {
  stop(
    "Required metadata columns are missing: ",
    paste(
      missing_metadata,
      collapse = ", "
    )
  )
}

if (!"SCT" %in% Assays(obj)) {
  stop(
    "SCT assay not found in the Seurat object."
  )
}


DefaultAssay(obj) <- "SCT"

pseudo_obj <- AggregateExpression(
  obj,
  group.by = c(
    "ID",
    "Genotype"
  ),
  assays = "SCT",
  slot = "counts",
  return.seurat = TRUE
)

Idents(pseudo_obj) <- "Genotype"

# =============================================================================
# 2. Differential expression using DESeq2
# =============================================================================

conditions <- c(
  "5XFAD;PS19",
  "5XFAD",
  "PS19"
)

combined_results <- list()

for (condition in conditions) {

  if (
    !condition %in% levels(Idents(pseudo_obj)) ||
      !"WT" %in% levels(Idents(pseudo_obj))
  ) {
    stop(
      "Required genotype comparison not found: ",
      condition,
      " vs WT"
    )
  }

  markers <- FindMarkers(
    pseudo_obj,
    ident.1 = condition,
    ident.2 = "WT",
    verbose = TRUE,
    min.cells.group = 1,
    min.cells.feature = 1,
    min.pct = 0,
    logfc.threshold = 0,
    only.pos = FALSE,
    test.use = "DESeq2"
  )

  markers$Gene <- rownames(
    markers
  )

  markers$Condition <- condition

  combined_results[[condition]] <- markers
}

combined_df <- do.call(
  rbind,
  combined_results
)

write.csv(
  combined_df,
  output_file,
  row.names = FALSE
)

message(
  "Step 11 complete: pseudobulk differential expression."
)
