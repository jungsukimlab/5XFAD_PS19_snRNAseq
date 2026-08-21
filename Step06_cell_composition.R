# =============================================================================
# 5XFAD;PS19 snRNA-seq analysis
# Step 6: Cell-type and cluster proportion analysis by genotype
# =============================================================================
#
# Run from the project root directory.
#
# Input:
#   publication_output/Step5_annotation/
#       5XFAD_PS19_Harmony_annotated.rds
#
# Outputs:
#   publication_output/Step6_cell_composition/
#
# Analyses:
#   1. Annotated cluster proportions by genotype
#      (clust_cell: ExN_1, InN_1, O_1, M_1, etc.)
#   2. Broad cell-type proportions by genotype
#      (cell_type: ExN, InN, O, M, A, OPC, EC, FB)
#
# R and package versions will be documented at the project level (renv.lock).

library(Seurat)
library(speckle)
library(limma)
library(openxlsx)
library(ggplot2)
library(dplyr)

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
  "Step6_cell_composition"
)

plot_dir <- file.path(
  output_dir,
  "plots"
)

results_file <- file.path(
  output_dir,
  "Step6_propeller_results.xlsx"
)

cluster_summary_file <- file.path(
  output_dir,
  "Step6_cluster_proportion_summary.csv"
)

cell_type_summary_file <- file.path(
  output_dir,
  "Step6_cell_type_proportion_summary.csv"
)

cluster_plot_file <- file.path(
  plot_dir,
  "Step6_cluster_proportions_by_genotype.pdf"
)

cell_type_plot_file <- file.path(
  plot_dir,
  "Step6_cell_type_proportions_by_genotype.pdf"
)

overwrite <- FALSE

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plot_dir,
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

output_files <- c(
  results_file,
  results_rds_file,
  cluster_summary_file,
  cell_type_summary_file,
  cluster_plot_file,
  cell_type_plot_file
)

invisible(
  lapply(
    output_files,
    check_output
  )
)

if (!file.exists(input_file)) {
  stop(
    "Input Seurat object not found: ",
    input_file
  )
}

# -----------------------------------------------------------------------------
# Load annotated object
# -----------------------------------------------------------------------------

obj <- readRDS(
  input_file
)

required_metadata <- c(
  "ID",
  "Genotype",
  "cell_type",
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

genotype_levels <- c(
  "WT",
  "5XFAD",
  "PS19",
  "5XFAD;PS19"
)

genotype_colors <- c(
  "WT" = "seagreen4",
  "5XFAD" = "darkgoldenrod1",
  "PS19" = "darkslategray3",
  "5XFAD;PS19" = "palevioletred"
)


# =============================================================================
# 1. Annotated cluster proportions by genotype
# =============================================================================

props_cluster <- getTransformedProps(
  clusters = obj$clust_cell,
  sample = obj$ID,
  transform = "logit"
)

# Match genotype labels to the exact sample order returned by Propeller.
cluster_sample_ids <- colnames(
  props_cluster$Proportions
)

group_cluster <- obj$Genotype[
  match(
    cluster_sample_ids,
    as.character(obj$ID)
  )
]

group_cluster <- factor(
  as.character(group_cluster),
  levels = genotype_levels
)

if (any(is.na(group_cluster))) {
  stop(
    "Could not assign genotype to one or more samples in the cluster ",
    "proportion analysis."
  )
}

# Design matrix.
design_cluster <- model.matrix(
  ~ 0 + group_cluster
)

colnames(design_cluster) <- c(
  "WT",
  "g5XFAD",
  "gPS19",
  "g5XFAD_PS19"
)

# -----------------------------------------------------------------------------
# Cluster proportion contrasts
# -----------------------------------------------------------------------------

contr_cluster_5XFAD <- makeContrasts(
  g5XFAD - WT,
  levels = design_cluster
)

contr_cluster_PS19 <- makeContrasts(
  gPS19 - WT,
  levels = design_cluster
)

contr_cluster_5XFAD_PS19 <- makeContrasts(
  g5XFAD_PS19 - WT,
  levels = design_cluster
)

# -----------------------------------------------------------------------------
# Cluster proportion tests
# -----------------------------------------------------------------------------

result_cluster_5XFAD <- as.data.frame(
  propeller.ttest(
    props_cluster,
    design = design_cluster,
    contrasts = contr_cluster_5XFAD,
    robust = TRUE,
    trend = FALSE,
    sort = TRUE
  )
)

result_cluster_PS19 <- as.data.frame(
  propeller.ttest(
    props_cluster,
    design = design_cluster,
    contrasts = contr_cluster_PS19,
    robust = TRUE,
    trend = FALSE,
    sort = TRUE
  )
)

result_cluster_5XFAD_PS19 <- as.data.frame(
  propeller.ttest(
    props_cluster,
    design = design_cluster,
    contrasts = contr_cluster_5XFAD_PS19,
    robust = TRUE,
    trend = FALSE,
    sort = TRUE
  )
)

# -----------------------------------------------------------------------------
# Cluster proportion summary for plotting
# -----------------------------------------------------------------------------

cluster_long <- as.data.frame(
  as.table(
    props_cluster$Proportions
  ),
  stringsAsFactors = FALSE
)

colnames(cluster_long) <- c(
  "Cluster",
  "ID",
  "Proportion"
)

cluster_long$Genotype <- group_cluster[
  match(
    cluster_long$ID,
    cluster_sample_ids
  )
]

# Preserve the cluster order established during Step 5 annotation.
if (is.factor(obj$clust_cell)) {
  cluster_order <- levels(
    obj$clust_cell
  )
} else {
  cluster_order <- rownames(
    props_cluster$Proportions
  )
}

cluster_long$Cluster <- factor(
  cluster_long$Cluster,
  levels = cluster_order
)

cluster_long$Genotype <- factor(
  cluster_long$Genotype,
  levels = genotype_levels
)

cluster_summary <- cluster_long %>%
  group_by(
    Cluster,
    Genotype
  ) %>%
  summarise(
    mean_proportion = mean(Proportion),
    sd_proportion = sd(Proportion),
    n = n(),
    se_proportion = sd_proportion / sqrt(n),
    .groups = "drop"
  )

write.csv(
  cluster_summary,
  cluster_summary_file,
  row.names = FALSE
)

pdf(
  cluster_plot_file,
  width = 14,
  height = 8
)

print(
  ggplot(
    cluster_summary,
    aes(
      x = Cluster,
      y = mean_proportion,
      fill = Genotype
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.75
    ) +
    geom_errorbar(
      aes(
        ymin = pmax(
          mean_proportion - se_proportion,
          0
        ),
        ymax = mean_proportion + se_proportion
      ),
      position = position_dodge(width = 0.8),
      width = 0.25
    ) +
    scale_fill_manual(
      values = genotype_colors
    ) +
    labs(
      x = "Cluster",
      y = "Proportion of nuclei",
      fill = "Genotype"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
)

dev.off()

# =============================================================================
# 2. Broad cell-type proportions by genotype
# =============================================================================

props_celltype <- getTransformedProps(
  clusters = obj$cell_type,
  sample = obj$ID,
  transform = "logit"
)

# Match genotype labels to the exact sample order returned by Propeller.
celltype_sample_ids <- colnames(
  props_celltype$Proportions
)

group_celltype <- obj$Genotype[
  match(
    celltype_sample_ids,
    as.character(obj$ID)
  )
]

group_celltype <- factor(
  as.character(group_celltype),
  levels = genotype_levels
)

if (any(is.na(group_celltype))) {
  stop(
    "Could not assign genotype to one or more samples in the cell-type ",
    "proportion analysis."
  )
}

# Design matrix.
design_celltype <- model.matrix(
  ~ 0 + group_celltype
)

colnames(design_celltype) <- c(
  "WT",
  "g5XFAD",
  "gPS19",
  "g5XFAD_PS19"
)

# -----------------------------------------------------------------------------
# Cell-type proportion contrasts
# -----------------------------------------------------------------------------

contr_celltype_5XFAD <- makeContrasts(
  g5XFAD - WT,
  levels = design_celltype
)

contr_celltype_PS19 <- makeContrasts(
  gPS19 - WT,
  levels = design_celltype
)

contr_celltype_5XFAD_PS19 <- makeContrasts(
  g5XFAD_PS19 - WT,
  levels = design_celltype
)

contr_celltype_5XFAD_PS19_vs_5XFAD <- makeContrasts(
  g5XFAD_PS19 - g5XFAD,
  levels = design_celltype
)

contr_celltype_5XFAD_PS19_vs_PS19 <- makeContrasts(
  g5XFAD_PS19 - gPS19,
  levels = design_celltype
)

# -----------------------------------------------------------------------------
# Cell-type proportion tests
# -----------------------------------------------------------------------------

result_celltype_5XFAD <- as.data.frame(
  propeller.ttest(
    props_celltype,
    design = design_celltype,
    contrasts = contr_celltype_5XFAD,
    robust = TRUE,
    trend = FALSE,
    sort = TRUE
  )
)

result_celltype_PS19 <- as.data.frame(
  propeller.ttest(
    props_celltype,
    design = design_celltype,
    contrasts = contr_celltype_PS19,
    robust = TRUE,
    trend = FALSE,
    sort = TRUE
  )
)

result_celltype_5XFAD_PS19 <- as.data.frame(
  propeller.ttest(
    props_celltype,
    design = design_celltype,
    contrasts = contr_celltype_5XFAD_PS19,
    robust = TRUE,
    trend = FALSE,
    sort = TRUE
  )
)

result_celltype_5XFAD_PS19_vs_5XFAD <- as.data.frame(
  propeller.ttest(
    props_celltype,
    design = design_celltype,
    contrasts = contr_celltype_5XFAD_PS19_vs_5XFAD,
    robust = TRUE,
    trend = FALSE,
    sort = TRUE
  )
)

result_celltype_5XFAD_PS19_vs_PS19 <- as.data.frame(
  propeller.ttest(
    props_celltype,
    design = design_celltype,
    contrasts = contr_celltype_5XFAD_PS19_vs_PS19,
    robust = TRUE,
    trend = FALSE,
    sort = TRUE
  )
)

# -----------------------------------------------------------------------------
# Cell-type proportion summary for plotting
# -----------------------------------------------------------------------------

celltype_long <- as.data.frame(
  as.table(
    props_celltype$Proportions
  ),
  stringsAsFactors = FALSE
)

colnames(celltype_long) <- c(
  "CellType",
  "ID",
  "Proportion"
)

celltype_long$Genotype <- group_celltype[
  match(
    celltype_long$ID,
    celltype_sample_ids
  )
]

cell_type_order <- c(
  "ExN",
  "InN",
  "O",
  "M",
  "A",
  "OPC",
  "EC",
  "FB"
)

cell_type_order <- intersect(
  cell_type_order,
  rownames(
    props_celltype$Proportions
  )
)

celltype_long$CellType <- factor(
  celltype_long$CellType,
  levels = cell_type_order
)

celltype_long$Genotype <- factor(
  celltype_long$Genotype,
  levels = genotype_levels
)

cell_type_summary <- celltype_long %>%
  group_by(
    CellType,
    Genotype
  ) %>%
  summarise(
    mean_proportion = mean(Proportion),
    sd_proportion = sd(Proportion),
    n = n(),
    se_proportion = sd_proportion / sqrt(n),
    .groups = "drop"
  )

write.csv(
  cell_type_summary,
  cell_type_summary_file,
  row.names = FALSE
)

pdf(
  cell_type_plot_file,
  width = 9,
  height = 7
)

print(
  ggplot(
    cell_type_summary,
    aes(
      x = CellType,
      y = mean_proportion,
      fill = Genotype
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.75
    ) +
    geom_errorbar(
      aes(
        ymin = pmax(
          mean_proportion - se_proportion,
          0
        ),
        ymax = mean_proportion + se_proportion
      ),
      position = position_dodge(width = 0.8),
      width = 0.25
    ) +
    scale_fill_manual(
      values = genotype_colors
    ) +
    labs(
      x = "Cell type",
      y = "Proportion of nuclei",
      fill = "Genotype"
    ) +
    theme_classic()
)

dev.off()

# =============================================================================
# Export statistical results
# =============================================================================

wb <- createWorkbook()

# Cluster comparisons
addWorksheet(
  wb,
  "Cluster_5XFAD_vs_WT"
)

writeData(
  wb,
  "Cluster_5XFAD_vs_WT",
  result_cluster_5XFAD,
  rowNames = TRUE
)

addWorksheet(
  wb,
  "Cluster_PS19_vs_WT"
)

writeData(
  wb,
  "Cluster_PS19_vs_WT",
  result_cluster_PS19,
  rowNames = TRUE
)

addWorksheet(
  wb,
  "Cluster_Double_vs_WT"
)

writeData(
  wb,
  "Cluster_Double_vs_WT",
  result_cluster_5XFAD_PS19,
  rowNames = TRUE
)

# Cell-type comparisons
addWorksheet(
  wb,
  "CellType_5XFAD_vs_WT"
)

writeData(
  wb,
  "CellType_5XFAD_vs_WT",
  result_celltype_5XFAD,
  rowNames = TRUE
)

addWorksheet(
  wb,
  "CellType_PS19_vs_WT"
)

writeData(
  wb,
  "CellType_PS19_vs_WT",
  result_celltype_PS19,
  rowNames = TRUE
)

addWorksheet(
  wb,
  "CellType_Double_vs_WT"
)

writeData(
  wb,
  "CellType_Double_vs_WT",
  result_celltype_5XFAD_PS19,
  rowNames = TRUE
)

addWorksheet(
  wb,
  "CellType_Double_vs_5XFAD"
)

writeData(
  wb,
  "CellType_Double_vs_5XFAD",
  result_celltype_5XFAD_PS19_vs_5XFAD,
  rowNames = TRUE
)

addWorksheet(
  wb,
  "CellType_Double_vs_PS19"
)

writeData(
  wb,
  "CellType_Double_vs_PS19",
  result_celltype_5XFAD_PS19_vs_PS19,
  rowNames = TRUE
)

saveWorkbook(
  wb,
  results_file,
  overwrite = overwrite
)
