# =============================================================================
# 5XFAD;PS19 snRNA-seq analysis
# Step 7: Differential expression analysis by annotated cluster
# =============================================================================
#
# Run from the project root directory.
#
# Input:
#   publication_output/Step5_annotation/
#       5XFAD_PS19_Harmony_annotated.rds
#
# Outputs:
#   publication_output/Step7_DEG/
#       5XFAD_PS19_RNA_normalized.rds
#       clust_DEGs/
#
# Analysis:
#   MAST differential expression within each annotated cluster (clust_cell):
#       5XFAD vs WT
#       PS19 vs WT
#       5XFAD;PS19 vs WT
#
# R and package versions will be documented at the project level (renv.lock).

library(Seurat)
library(MAST)

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
  "Step7_DEG"
)

deg_dir <- file.path(
  output_dir,
  "clust_DEGs"
)

normalized_rna_file <- file.path(
  output_dir,
  "5XFAD_PS19_RNA_normalized.rds"
)

overwrite <- FALSE

dir.create(
  deg_dir,
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

if (!file.exists(input_file)) {
  stop(
    "Input Seurat object not found: ",
    input_file
  )
}

check_output(
  normalized_rna_file
)

# -----------------------------------------------------------------------------
# Load annotated object
# -----------------------------------------------------------------------------

obj <- readRDS(
  input_file
)

required_metadata <- c(
  "Genotype",
  "clust_cell"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(obj[[]])
)

if (length(missing_metadata) > 0) {
  stop(
    "Required metadata columns are missing: ",
    paste(missing_metadata, collapse = ", ")
  )
}

if (!"RNA" %in% Assays(obj)) {
  stop(
    "RNA assay not found in the Seurat object."
  )
}

# =============================================================================
# 1. Prepare RNA assay
# =============================================================================

DefaultAssay(obj) <- "RNA"

# Join sample-specific RNA layers before differential expression.
obj[["RNA"]] <- JoinLayers(
  obj[["RNA"]]
)

# Log-normalize RNA counts.
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

# Save the joined, log-normalized RNA object for Step 8.
saveRDS(
  obj,
  normalized_rna_file
)


# =============================================================================
# 2. Set cluster-by-genotype identities
# =============================================================================

obj$DEG <- paste0(
  obj$clust_cell,
  "_",
  obj$Genotype
)

Idents(obj) <- "DEG"

conditions <- c(
  "5XFAD",
  "PS19",
  "5XFAD;PS19"
)

clusters <- unique(
  as.character(obj$clust_cell)
)

# =============================================================================
# 3. MAST differential expression
# =============================================================================

for (clust in clusters) {
  
  for (condition in conditions) {
    
    ident_1 <- paste0(
      clust,
      "_",
      condition
    )
    
    ident_2 <- paste0(
      clust,
      "_WT"
    )
    
    if (
      ident_1 %in% levels(Idents(obj)) &&
      ident_2 %in% levels(Idents(obj))
    ) {
      
      message(
        "Running MAST: ",
        clust,
        " | ",
        condition,
        " vs WT"
      )
      
      markers <- FindMarkers(
        object = obj,
        ident.1 = ident_1,
        ident.2 = ident_2,
        assay = "RNA",
        test.use = "MAST",
        min.cells.group = 1,
        min.cells.feature = 1,
        min.pct = 0,
        logfc.threshold = 0,
        only.pos = FALSE,
        verbose = TRUE
      )
      
      result_name <- paste0(
        "clust_DEG_",
        clust,
        "_",
        condition,
        "_vs_WT"
      )
      
      result_file <- file.path(
        deg_dir,
        paste0(
          result_name,
          ".csv"
        )
      )
      
      check_output(
        result_file
      )
      
      write.csv(
        markers,
        result_file,
        row.names = TRUE
      )
      
    } else {
      
      message(
        "Skipping: ",
        ident_1,
        " or ",
        ident_2,
        " not found."
      )
    }
  }
}