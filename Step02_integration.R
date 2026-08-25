# 5XFAD;PS19 snRNA-seq analysis
# Step 2: SCT-based sample integration
#
# Run from the project root directory.
#
# Input:
#   ./publication_output/RDS/postSXandDF_<sample_id>_20250111.rds
#
# Output:
#   ./publication_output/Step2_integration/

library(Seurat)
library(harmony)
set.seed(4444)

# -----------------------------------------------------------------------------
# Paths and analysis parameters
# -----------------------------------------------------------------------------

input_dir <- file.path("publication_output", "RDS")
output_dir <- file.path("publication_output", "Step2_integration")
plot_dir <- file.path(output_dir, "plots")

analysis_tag <- "20250111"

# Prevent accidental replacement of existing publication outputs.
# Change to TRUE only when intentional overwriting is desired.
overwrite <- FALSE

# Main integration parameters
n_integration_features <- 3000
n_pcs <- 50
integration_dims <- 1:12
cluster_resolution <- 0.1

# Hard-doublet thresholds retained from the original analysis.
hard_doublet_feature_threshold <- 8000
hard_doublet_count_threshold <- 40000

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

samples <- c(
  "1608", "1609", "1617", "1619", "1620", "1624", "1625", "1626",
  "1627", "1629", "1636", "1643", "1646", "1699", "1702", "1705"
)

# -----------------------------------------------------------------------------
# Load Step 1 objects and define combined doublet calls
# -----------------------------------------------------------------------------

sample_objects <- vector("list", length(samples))
names(sample_objects) <- samples

doublet_summary <- vector("list", length(samples))

for (i in seq_along(samples)) {
  sample_id <- samples[i]
  
  input_file <- file.path(
    input_dir,
    paste0("postSXandDF_", sample_id, "_", analysis_tag, ".rds")
  )
  
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }
  
  message("Loading sample ", sample_id)
  
  obj <- readRDS(input_file)
  
  Idents(obj) <- "orig.ident"
  DefaultAssay(obj) <- "SCT"
  
  # Hard-doublet call 
  # nCount_RNA >= 40,000 was based on the sequencing/loading design.
  # nFeature_RNA >= 8,000 was used as a high-feature hard threshold.
  obj$HardDoublet <- ifelse(
    obj$nFeature_RNA >= hard_doublet_feature_threshold |
      obj$nCount_RNA >= hard_doublet_count_threshold,
    "Doublet",
    "Singlet"
  )
  
  # Combined doublet call: DoubletFinder (DF) OR hard-threshold doublet (HardDoublet).
  obj$Doublet <- ifelse(
    obj$DF == "Doublet" | obj$HardDoublet == "Doublet",
    "Doublet",
    "Singlet"
  )
  
  # Diagnostic plot for the hard and combined doublet calls.
  doublet_plot_file <- file.path(
    plot_dir,
    paste0(analysis_tag, "_DoubletComparison_", sample_id, ".pdf")
  )
  check_output(doublet_plot_file)
  
  pdf(doublet_plot_file)
  print(
    DimPlot(
      obj,
      reduction = "umap",
      group.by = c("HardDoublet", "DF", "Doublet")
    )
  )
  dev.off()
  
  doublet_summary[[i]] <- data.frame(
    sample_id = sample_id,
    n_nuclei = ncol(obj),
    DF_doublets = sum(obj$DF == "Doublet", na.rm = TRUE),
    HardDoublet_doublets = sum(obj$HardDoublet == "Doublet", na.rm = TRUE),
    Combined_doublets = sum(obj$Doublet == "Doublet", na.rm = TRUE),
    Combined_doublet_percent = 100 *
      sum(obj$Doublet == "Doublet", na.rm = TRUE) / ncol(obj),
    stringsAsFactors = FALSE
  )
  
  # Doublets are flagged but retained during this integration step
  sample_objects[[sample_id]] <- obj
  
  rm(obj)
}

doublet_summary <- do.call(rbind, doublet_summary)

doublet_summary_file <- file.path(
  output_dir,
  paste0(analysis_tag, "_Step2_doublet_summary.csv")
)
check_output(doublet_summary_file)

write.csv(
  doublet_summary,
  file = doublet_summary_file,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Select integration features across samples
# -----------------------------------------------------------------------------

integration_features <- SelectIntegrationFeatures(
  object.list = sample_objects,
  nfeatures = n_integration_features
)

# -----------------------------------------------------------------------------
# Merge SCT-normalized samples
# -----------------------------------------------------------------------------

merged_seurat <- merge(
  x = sample_objects[[1]],
  y = sample_objects[2:length(sample_objects)],
  merge.data = TRUE
)

DefaultAssay(merged_seurat) <- "SCT"

# Manually set the consensus variable features selected across samples.
VariableFeatures(merged_seurat) <- integration_features

# PCA on the merged SCT assay.
merged_seurat <- RunPCA(
  merged_seurat,
  assay = "SCT",
  npcs = n_pcs
)

preintegration_file <- file.path(output_dir, "preIntegration.rds")
check_output(preintegration_file)
saveRDS(merged_seurat, preintegration_file)


# -----------------------------------------------------------------------------
# Harmony integration
# -----------------------------------------------------------------------------

# Harmony was used to obtain a shared low-dimensional representation in which
# corresponding cell populations align across experimental groups.

obj <- merged_seurat
rm(merged_seurat)
gc()

DefaultAssay(obj) <- "SCT"

obj <- harmony::RunHarmony(
  object = obj,
  group.by.vars = c("Genotype", "Sex"),
  reduction.use = "pca",
  dims.use = seq_len(n_pcs),
  reduction.save = "harmony",
  verbose = TRUE
)

harmony_checkpoint_file <- file.path(
  output_dir,
  "SCTintegrated_Harmony_before_clustering.rds"
)
check_output(harmony_checkpoint_file)

saveRDS(
  obj,
  harmony_checkpoint_file
)

# -----------------------------------------------------------------------------
# Dimensionality reduction and clustering using Harmony embeddings
# -----------------------------------------------------------------------------

obj <- FindNeighbors(
  obj,
  dims = integration_dims,
  reduction = "harmony"
)

obj <- FindClusters(
  obj,
  resolution = cluster_resolution
)

obj <- RunUMAP(
  obj,
  dims = integration_dims,
  reduction = "harmony"
)

# -----------------------------------------------------------------------------
# Harmony UMAP diagnostics
# -----------------------------------------------------------------------------

umap_genotype_file <- file.path(
  plot_dir,
  paste0(analysis_tag, "_UMAP_Harmony_Genotype.pdf")
)
check_output(umap_genotype_file)

pdf(umap_genotype_file)
print(
  DimPlot(
    obj,
    reduction = "umap",
    group.by = "Genotype",
    raster = FALSE,
    alpha = 0.1,
    label = FALSE
  )
)
dev.off()

umap_sex_file <- file.path(
  plot_dir,
  paste0(analysis_tag, "_UMAP_Harmony_Sex.pdf")
)
check_output(umap_sex_file)

pdf(umap_sex_file)
print(
  DimPlot(
    obj,
    reduction = "umap",
    group.by = "Sex",
    raster = FALSE,
    alpha = 0.1,
    label = FALSE
  )
)
dev.off()

umap_sample_file <- file.path(
  plot_dir,
  paste0(analysis_tag, "_UMAP_Harmony_Sample.pdf")
)
check_output(umap_sample_file)

pdf(umap_sample_file)
print(
  DimPlot(
    obj,
    reduction = "umap",
    group.by = "ID",
    raster = FALSE,
    alpha = 0.1,
    label = FALSE
  )
)
dev.off()

# -----------------------------------------------------------------------------
# Save final Harmony-integrated object
# -----------------------------------------------------------------------------

integrated_file <- file.path(
  output_dir,
  "SCTintegrated_Harmony.rds"
)
check_output(integrated_file)

saveRDS(
  obj,
  integrated_file
)
