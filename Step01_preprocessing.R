# 5XFAD;PS19 snRNA-seq analysis
# Step 1: Ambient RNA correction, quality control, and doublet detection
#
# Run from the project root directory.
# Input:  Cell Ranger output directories in ./data/<sample_id>/
# Output: Publication rerun outputs in ./publication_output/

library(Seurat)
library(SoupX)
library(DoubletFinder)

# -----------------------------------------------------------------------------
# Paths and analysis parameters
# -----------------------------------------------------------------------------

data_dir <- "data"
output_dir <- "publication_output"
qc_dir <- file.path(output_dir, "QC_plots")
rds_dir <- file.path(output_dir, "RDS")

# Identifier retained from the original preprocessing analysis.
analysis_tag <- "20250111"

min_features <- 200
max_percent_mt <- 1
n_pcs <- 10
cluster_resolution <- 0.2
pN <- 0.25

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)


samples <- c(
  "1608", "1609", "1617", "1619", "1620", "1624", "1625", "1626",
  "1627", "1629", "1636", "1643", "1646", "1699", "1702", "1705"
)

# Sample-level metadata.
sample_metadata <- data.frame(
  sample_id = c(
    "1608", "1609", "1617", "1702",
    "1636", "1699", "1624", "1626",
    "1620", "1643", "1629", "1705",
    "1619", "1646", "1625", "1627"
  ),
  ID = paste0("S", 1:16),
  Genotype = c(
    rep("WT", 4),
    rep("5XFAD", 4),
    rep("PS19", 4),
    rep("5XFAD;PS19", 4)
  ),
  Sex = c(
    "M", "M", "F", "F",
    "M", "M", "F", "F",
    "M", "M", "F", "F",
    "M", "M", "F", "F"
  ),
  Age = c(
    "9.4", "9.4", "9.4", "9.4",
    "9.5", "9.4", "9.4", "9.4",
    "9.4", "9.5", "9.5", "9.4",
    "9.5", "9.5", "8.8", "9.5"
  ),
  stringsAsFactors = FALSE
)

stopifnot(setequal(samples, sample_metadata$sample_id))

qc_summary <- vector("list", length(samples))

# -----------------------------------------------------------------------------
# Process each sample independently
# -----------------------------------------------------------------------------

for (i in seq_along(samples)) {
  sample_id <- samples[i]
  sample_info <- sample_metadata[sample_metadata$sample_id == sample_id, ]
  sample_dir <- file.path(data_dir, sample_id)
  
  if (!dir.exists(sample_dir)) {
    stop("Input directory not found: ", sample_dir)
  }
  
  message("Processing sample ", sample_id, " (", sample_info$ID, ")")
  
  # ---------------------------------------------------------------------------
  # Define output files 
  # ---------------------------------------------------------------------------
  
  contam_pdf <- file.path(
    qc_dir,
    paste0(analysis_tag, "_Contam", sample_id, ".pdf")
  )
  
  qc_pdf <- file.path(
    qc_dir,
    paste0(analysis_tag, "_QCplots", sample_id, ".pdf")
  )
  
  pk_csv <- file.path(
    qc_dir,
    paste0("SoupXpK", sample_id, ".csv")
  )
  
  doublet_qc_pdf <- file.path(
    qc_dir,
    paste0(analysis_tag, "_DoubletQCplot", sample_id, ".pdf")
  )
  
  doublet_umap_pdf <- file.path(
    qc_dir,
    paste0(analysis_tag, "_DoubletUMAP", sample_id, ".pdf")
  )
  
  rds_file <- file.path(
    rds_dir,
    paste0("postSXandDF_", sample_id, "_", analysis_tag, ".rds")
  )
  
  output_files <- c(
    contam_pdf,
    qc_pdf,
    pk_csv,
    doublet_qc_pdf,
    doublet_umap_pdf,
    rds_file
  )
  
  invisible(lapply(output_files, check_output))
  
  # ---------------------------------------------------------------------------
  # SoupX ambient RNA correction
  # ---------------------------------------------------------------------------
  
  sample_soup <- load10X(sample_dir)
  
  pdf(contam_pdf)
  sample_soup <- autoEstCont(sample_soup)
  print(sample_soup)
  dev.off()
  
  soup_rho <- sample_soup$fit$rhoEst
  corrected_counts <- adjustCounts(sample_soup)
  
  # ---------------------------------------------------------------------------
  # Create Seurat object and perform nucleus-level QC
  # ---------------------------------------------------------------------------
  
  sample_seurat <- CreateSeuratObject(
    counts = corrected_counts,
    project = "5XFAD_PS19_10X",
    min.cells = 3,
    min.features = min_features
  )
  
  sample_seurat[["percent.mt"]] <- PercentageFeatureSet(
    sample_seurat,
    pattern = "^mt-"
  )
  
  n_nuclei_pre_qc <- ncol(sample_seurat)
  
  pdf(qc_pdf)
  print(
    VlnPlot(
      sample_seurat,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      ncol = 3
    )
  )
  dev.off()
  
  # QC thresholds used in this study for snRNA-seq nuclei.
  sample_seurat <- subset(
    sample_seurat,
    subset = nFeature_RNA > min_features & percent.mt < max_percent_mt
  )
  
  n_nuclei_post_qc <- ncol(sample_seurat)
  
  # ---------------------------------------------------------------------------
  # Preprocessing for DoubletFinder
  # ---------------------------------------------------------------------------
  
  sample_seurat <- SCTransform(sample_seurat)
  sample_seurat <- RunPCA(sample_seurat)
  sample_seurat <- FindNeighbors(sample_seurat, dims = seq_len(n_pcs))
  
  # A low clustering resolution was used so that homotypic doublets were
  # estimated from broad cell populations rather than fine subclusters.
  sample_seurat <- FindClusters(
    sample_seurat,
    resolution = cluster_resolution
  )
  
  sample_seurat <- RunUMAP(
    sample_seurat,
    dims = seq_len(n_pcs)
  )
  
  # ---------------------------------------------------------------------------
  # DoubletFinder parameter sweep and doublet classification
  # ---------------------------------------------------------------------------
  
  sweep_res <- paramSweep(
    sample_seurat,
    PCs = seq_len(n_pcs),
    sct = TRUE
  )
  
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  
  write.csv(
    bcmvn,
    file = pk_csv,
    row.names = FALSE
  )
  
  best_pk_row <- which.max(bcmvn$BCmetric)
  pK <- as.numeric(as.character(bcmvn$pK[best_pk_row]))
  
  # Expected multiplet-rate calculation retained from the original analysis.
  multiplet_rate <- (ncol(sample_seurat) / 0.57) * 4.6e-6
  
  homotypic_prop <- modelHomotypic(sample_seurat$seurat_clusters)
  n_exp_poi <- round(multiplet_rate * ncol(sample_seurat))
  n_exp_poi_adj <- round(n_exp_poi * (1 - homotypic_prop))
  
  # reuse.pANN = FALSE is retained from the original analysis. The compatible
  # DoubletFinder version will be recorded in the project dependency lockfile.
  sample_seurat <- doubletFinder(
    sample_seurat,
    PCs = seq_len(n_pcs),
    pN = pN,
    pK = pK,
    nExp = n_exp_poi_adj,
    reuse.pANN = NULL,
    sct = TRUE
  )
  
  # Standardize DoubletFinder metadata column names for downstream analyses.
  metadata_names <- colnames(sample_seurat[[]])
  pann_col <- grep("^pANN_", metadata_names, value = TRUE)
  df_col <- grep("^DF\\.classifications_", metadata_names, value = TRUE)
  
  if (length(pann_col) != 1 || length(df_col) != 1) {
    stop(
      "Unexpected DoubletFinder metadata columns for sample ", sample_id,
      ". pANN columns: ", paste(pann_col, collapse = ", "),
      "; classification columns: ", paste(df_col, collapse = ", ")
    )
  }
  
  sample_seurat$pANN <- sample_seurat[[]][[pann_col]]
  sample_seurat$DF <- sample_seurat[[]][[df_col]]
  sample_seurat[[pann_col]] <- NULL
  sample_seurat[[df_col]] <- NULL
  
  n_doublets <- sum(sample_seurat$DF == "Doublet", na.rm = TRUE)
  doublet_percent <- 100 * n_doublets / ncol(sample_seurat)
  
  # ---------------------------------------------------------------------------
  # DoubletFinder QC plots
  # ---------------------------------------------------------------------------
  
  pdf(doublet_qc_pdf)
  print(
    VlnPlot(
      sample_seurat,
      features = c("nFeature_SCT", "nCount_SCT", "percent.mt"),
      group.by = "DF"
    )
  )
  dev.off()
  
  pdf(doublet_umap_pdf)
  print(
    DimPlot(
      sample_seurat,
      reduction = "umap",
      group.by = "DF"
    )
  )
  dev.off()
  
  # Doublets are retained at this stage and removed later in the pipeline,
  # matching the original analysis workflow.
  
  # ---------------------------------------------------------------------------
  # Add sample metadata and save processed object
  # ---------------------------------------------------------------------------
  
  sample_seurat$ID <- sample_info$ID
  sample_seurat$Genotype <- sample_info$Genotype
  sample_seurat$Sex <- sample_info$Sex
  sample_seurat$Age <- sample_info$Age
  
  saveRDS(sample_seurat, file = rds_file)
  
  # Save key QC metrics for reporting and reproducibility.
  qc_summary[[i]] <- data.frame(
    sample_id = sample_id,
    ID = sample_info$ID,
    Genotype = sample_info$Genotype,
    Sex = sample_info$Sex,
    Age = sample_info$Age,
    SoupX_rho = soup_rho,
    nuclei_pre_qc = n_nuclei_pre_qc,
    nuclei_post_qc = n_nuclei_post_qc,
    pK = pK,
    multiplet_rate = multiplet_rate,
    homotypic_proportion = homotypic_prop,
    expected_doublets_unadjusted = n_exp_poi,
    expected_doublets_adjusted = n_exp_poi_adj,
    predicted_doublets = n_doublets,
    predicted_doublet_percent = doublet_percent,
    stringsAsFactors = FALSE
  )
  
  rm(sample_soup, corrected_counts, sample_seurat)
  gc()
}

# -----------------------------------------------------------------------------
# Save Step 1.1 QC summary
# -----------------------------------------------------------------------------

qc_summary <- do.call(rbind, qc_summary)

qc_summary_file <- file.path(
  qc_dir,
  paste0(analysis_tag, "_Step1.1_QC_summary.csv")
)

check_output(qc_summary_file)

write.csv(
  qc_summary,
  file = qc_summary_file,
  row.names = FALSE
)
