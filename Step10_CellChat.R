# =============================================================================
# 5XFAD;PS19 snRNA-seq analysis
# Step 10: CellChat analysis
# =============================================================================
#
# Run from the project root directory.
#
# Analysis:
#   1. Create CellChat objects for each genotype and merge them
#   2. Interaction analysis
#   3. Communication pattern analysis
#   4. Pathway-specific ligand-receptor analysis
#
# R and package versions will be documented at the project level (renv.lock).

library(Seurat)
library(CellChat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(NMF)
library(ggalluvial)
library(future)

plan("sequential")

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

input_file <- file.path(
  "publication_output",
  "Step7_DEG",
  "5XFAD_PS19_RNA_normalized.rds"
)

output_dir <- file.path(
  "publication_output",
  "Step10_CellChat"
)

object_dir <- file.path(output_dir, "objects")
interaction_dir <- file.path(output_dir, "interactions")
pattern_dir <- file.path(output_dir, "communication_patterns")
signaling_dir <- file.path(output_dir, "PROS_VISTA")

dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(interaction_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pattern_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(signaling_dir, recursive = TRUE, showWarnings = FALSE)

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
  stop("Input Seurat object not found: ", input_file)
}

# -----------------------------------------------------------------------------
# Cluster and genotype order
# -----------------------------------------------------------------------------

cluster_order <- c(
  "ExN_1", "ExN_2", "ExN_3", "ExN_4", "ExN_5", "ExN_6", "ExN_7",
  "ExN_8", "ExN_9", "ExN_10", "ExN_11", "ExN_12", "ExN_13",
  "InN_1", "InN_2", "InN_3", "InN_4", "InN_5", "InN_6", "InN_7",
  "InN_8", "InN_9", "InN_10", "InN_11", "InN_12", "InN_13",
  "O_1", "O_2", "O_3",
  "M_1", "M_2",
  "A_1", "A_2", "A_3",
  "OPC",
  "EC_1", "EC_2",
  "FB"
)

genotypes <- c(
  "WT",
  "5XFAD",
  "PS19",
  "5XFAD;PS19"
)

genotype_file_names <- c(
  "WT" = "WT",
  "5XFAD" = "5XFAD",
  "PS19" = "PS19",
  "5XFAD;PS19" = "5XFAD_PS19"
)

genotype_colors <- c(
  "seagreen4",
  "darkgoldenrod1",
  "darkslategray3",
  "palevioletred"
)

# =============================================================================
# 1. CellChat analysis
# =============================================================================

obj <- readRDS(input_file)

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

DefaultAssay(obj) <- "RNA"

obj$Genotype <- factor(
  as.character(obj$Genotype),
  levels = genotypes
)

obj$clust_cell <- factor(
  as.character(obj$clust_cell),
  levels = cluster_order
)

options(
  future.globals.maxSize = 4 * 1024^3
)

object.list <- list()

for (geno in genotypes) {
  
  cellchat_file <- file.path(
    object_dir,
    paste0(
      "cellchat_",
      genotype_file_names[[geno]],
      ".rds"
    )
  )
  
  check_output(cellchat_file)
  
  sub_obj <- subset(
    obj,
    cells = colnames(obj)[
      as.character(obj$Genotype) == geno
    ]
  )
  
  sub_obj$clust_cell <- droplevels(
    sub_obj$clust_cell
  )
  
  Idents(sub_obj) <- "clust_cell"
  
  cellchat_obj <- createCellChat(
    object = sub_obj,
    group.by = "clust_cell",
    assay = "RNA"
  )
  
  cellchat_obj@DB <- CellChatDB.mouse
  
  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  
  cellchat_obj <- computeCommunProb(
    cellchat_obj,
    type = "triMean"
  )
  
  cellchat_obj <- filterCommunication(
    cellchat_obj,
    min.cells = 10
  )
  
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  
  cellchat_obj <- netAnalysis_computeCentrality(
    cellchat_obj,
    slot.name = "netP"
  )
  
  saveRDS(
    cellchat_obj,
    cellchat_file
  )
  
  object.list[[geno]] <- cellchat_obj
}

# -----------------------------------------------------------------------------
# Lift up CellChat objects and merge them together
# -----------------------------------------------------------------------------
#
# WT and 5XFAD;PS19 contained all 38 clusters,
# while 5XFAD and PS19 lacked ExN_13.

group.new <- levels(
  object.list[["WT"]]@idents
)

for (geno in genotypes[-1]) {
  
  current_groups <- levels(
    object.list[[geno]]@idents
  )
  
  if (!identical(current_groups, group.new)) {
    object.list[[geno]] <- liftCellChat(
      object.list[[geno]],
      group.new
    )
  }
}

for (geno in genotypes) {
  object.list[[geno]]@idents <- factor(
    as.character(object.list[[geno]]@idents),
    levels = cluster_order
  )
}

merged_file <- file.path(
  object_dir,
  "cellchat_merged.rds"
)

object_list_file <- file.path(
  object_dir,
  "cellchat_object_list.rds"
)

check_output(merged_file)
check_output(object_list_file)

cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)

saveRDS(
  object.list,
  object_list_file
)

saveRDS(
  cellchat,
  merged_file
)

# -----------------------------------------------------------------------------
# Total number of interactions and interaction strength
# -----------------------------------------------------------------------------

interaction_summary_file <- file.path(
  interaction_dir,
  "interaction_number_and_strength.pdf"
)

check_output(interaction_summary_file)

custom_theme <- theme(
  text = element_text(size = 14),
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 16),
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.text.y = element_text(size = 14)
)

p_number <- compareInteractions(
  cellchat,
  show.legend = FALSE,
  group = 1:4
) +
  scale_fill_manual(values = genotype_colors) +
  custom_theme

p_strength <- compareInteractions(
  cellchat,
  show.legend = FALSE,
  group = 1:4,
  measure = "weight"
) +
  scale_fill_manual(values = genotype_colors) +
  custom_theme

pdf(
  interaction_summary_file,
  width = 8,
  height = 4
)

print(
  p_number + p_strength
)

dev.off()

# -----------------------------------------------------------------------------
# Cell-type interaction network in 5XFAD;PS19
# -----------------------------------------------------------------------------

clusters_use <- cluster_order[
  grepl(
    "^(ExN_|InN_|O_|M_|A_)",
    cluster_order
  )
]

cellchat_sub <- subsetCellChat(
  object.list[["5XFAD;PS19"]],
  idents.use = clusters_use
)

group.cellType <- c(
  rep("ExN", 13),
  rep("InN", 13),
  rep("O", 3),
  rep("M", 2),
  rep("A", 3)
)

group.cellType <- factor(
  group.cellType,
  levels = c(
    "ExN",
    "InN",
    "O",
    "M",
    "A"
  )
)

cellchat_sub <- mergeInteractions(
  cellchat_sub,
  group.cellType
)

weight.max <- max(
  cellchat_sub@net$count.merged,
  na.rm = TRUE
)

coarse_network_file <- file.path(
  interaction_dir,
  "cell_type_interactions_5XFAD_PS19.pdf"
)

check_output(
  coarse_network_file
)

pdf(
  coarse_network_file,
  width = 7,
  height = 7
)

netVisual_circle(
  cellchat_sub@net$count.merged,
  weight.scale = TRUE,
  label.edge = TRUE,
  edge.weight.max = weight.max,
  edge.width.max = 12,
  edge.label.cex = 1.4,
  vertex.label.cex = 1.6
)

dev.off()

# =============================================================================
# 2. Communication pattern analysis
# =============================================================================

# -----------------------------------------------------------------------------
# Overall information flow: WT versus 5XFAD;PS19
# -----------------------------------------------------------------------------

information_flow_file <- file.path(
  pattern_dir,
  "information_flow_WT_vs_5XFAD_PS19.pdf"
)

check_output(information_flow_file)

p_flow <- rankNet(
  cellchat,
  mode = "comparison",
  comparison = c(1, 4),
  measure = "weight",
  cutoff.pvalue = 0.05,
  thresh = 0.05,
  sources.use = NULL,
  targets.use = NULL,
  stacked = FALSE,
  do.stat = TRUE
)

ggsave(
  information_flow_file,
  plot = p_flow,
  width = 7,
  height = 12
)

# -----------------------------------------------------------------------------
# Outgoing communication patterns in 5XFAD;PS19
# -----------------------------------------------------------------------------

object.list[["5XFAD;PS19"]] <- identifyCommunicationPatterns(
  object.list[["5XFAD;PS19"]],
  pattern = "outgoing",
  width = 6,
  height = 16,
  k = 9
)

outgoing_file <- file.path(
  pattern_dir,
  "outgoing_communication_patterns.pdf"
)

check_output(outgoing_file)

pdf(
  outgoing_file,
  width = 8,
  height = 12
)

netAnalysis_river(
  object.list[["5XFAD;PS19"]],
  pattern = "outgoing",
  cutoff = 0.5,
  font.size = 2.7,
  color.use.signaling = "white"
)

dev.off()

# oligodendrocyte/microglia outgoing patterns
patterns_to_keep <- c(
  "Pattern 2",
  "Pattern 8"
)

cellchat_filtered <- object.list[["5XFAD;PS19"]]

cellchat_filtered@netP$pattern$outgoing$pattern$signaling <-
  cellchat_filtered@netP$pattern$outgoing$pattern$signaling[
    cellchat_filtered@netP$pattern$outgoing$pattern$signaling$Pattern %in%
      patterns_to_keep,
    ,
    drop = FALSE
  ]

outgoing_om_file <- file.path(
  pattern_dir,
  "outgoing_patterns_oligodendrocyte_microglia.pdf"
)

check_output(outgoing_om_file)

pdf(
  outgoing_om_file,
  width = 8,
  height = 7
)

netAnalysis_river(
  cellchat_filtered,
  pattern = "outgoing",
  sources.use = c(
    "O_1",
    "M_1",
    "O_2",
    "M_2",
    "O_3"
  ),
  cutoff = 0.5,
  font.size = 2.7,
  color.use.signaling = "white"
)

dev.off()

# -----------------------------------------------------------------------------
# Incoming communication patterns in 5XFAD;PS19
# -----------------------------------------------------------------------------

object.list[["5XFAD;PS19"]] <- identifyCommunicationPatterns(
  object.list[["5XFAD;PS19"]],
  pattern = "incoming",
  width = 6,
  height = 16,
  k = 3
)

incoming_file <- file.path(
  pattern_dir,
  "incoming_communication_patterns.pdf"
)

check_output(incoming_file)

pdf(
  incoming_file,
  width = 8,
  height = 12
)

netAnalysis_river(
  object.list[["5XFAD;PS19"]],
  pattern = "incoming",
  cutoff = 0.5,
  font.size = 2.7,
  color.use.signaling = "white"
)

dev.off()

# =============================================================================
# 3. PROS and VISTA signaling
# =============================================================================

# -----------------------------------------------------------------------------
# Custom signaling heatmap with shared WT and 5XFAD;PS19 barplot limits
# -----------------------------------------------------------------------------

netVisual_heatmap_edit <- function(
    object,
    signaling,
    ordered_levels,
    row.limit,
    col.limit,
    color.heatmap = "Reds",
    title.name = NULL,
    font.size = 8,
    font.size.title = 10
) {
  
  mat <- object@netP$prob[, , signaling]
  
  mat <- mat[
    ordered_levels,
    ordered_levels,
    drop = FALSE
  ]
  
  mat[is.na(mat)] <- 0
  
  color.use <- scPalette(
    ncol(mat)
  )
  
  names(color.use) <- colnames(mat)
  
  heatmap_colors <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(
      n = 9,
      name = color.heatmap
    )
  )(100)
  
  df <- data.frame(
    group = colnames(mat)
  )
  
  rownames(df) <- colnames(mat)
  
  col_annotation <- HeatmapAnnotation(
    df = df,
    col = list(group = color.use),
    which = "column",
    show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(0.2, "cm")
  )
  
  row_annotation <- HeatmapAnnotation(
    df = df,
    col = list(group = color.use),
    which = "row",
    show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(0.2, "cm")
  )
  
  right_annotation <- rowAnnotation(
    Strength = anno_barplot(
      rowSums(abs(mat)),
      border = FALSE,
      gp = grid::gpar(
        fill = color.use,
        col = color.use
      ),
      ylim = c(0, row.limit)
    ),
    show_annotation_name = FALSE
  )
  
  top_annotation <- HeatmapAnnotation(
    Strength = anno_barplot(
      colSums(abs(mat)),
      border = FALSE,
      gp = grid::gpar(
        fill = color.use,
        col = color.use
      ),
      ylim = c(0, col.limit)
    ),
    show_annotation_name = FALSE
  )
  
  mat[mat == 0] <- NA
  
  Heatmap(
    mat,
    col = heatmap_colors,
    na_col = "white",
    name = "Communication Prob.",
    bottom_annotation = col_annotation,
    left_annotation = row_annotation,
    top_annotation = top_annotation,
    right_annotation = right_annotation,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_rot = 0,
    row_names_gp = grid::gpar(fontsize = font.size),
    column_names_gp = grid::gpar(fontsize = font.size),
    column_title = title.name,
    column_title_gp = grid::gpar(fontsize = font.size.title),
    column_names_rot = 90,
    row_title = "Sources (Sender)",
    row_title_gp = grid::gpar(fontsize = font.size.title),
    row_title_rot = 90,
    heatmap_legend_param = list(
      title_gp = grid::gpar(
        fontsize = 8,
        fontface = "plain"
      ),
      title_position = "leftcenter-rot",
      border = NA,
      legend_height = grid::unit(20, "mm"),
      labels_gp = grid::gpar(fontsize = 8),
      grid_width = grid::unit(2, "mm")
    )
  )
}

# =============================================================================
# PROS signaling
# =============================================================================

pathways.show <- "PROS"

# -----------------------------------------------------------------------------
# Ligand-receptor contribution
# -----------------------------------------------------------------------------

pros_contribution_file <- file.path(
  signaling_dir,
  "PROS_LR_contribution_5XFAD_PS19.pdf"
)

check_output(pros_contribution_file)

p <- netAnalysis_contribution(
  object.list[["5XFAD;PS19"]],
  signaling = pathways.show
)

ggsave(
  pros_contribution_file,
  plot = p,
  width = 6,
  height = 3
)

# -----------------------------------------------------------------------------
# Pros1-Mertk ligand-receptor interaction
# -----------------------------------------------------------------------------

weight.max <- getMaxWeight(
  object.list,
  slot.name = "netP",
  attribute = pathways.show
)

pairLR.PROS <- extractEnrichedLR(
  cellchat,
  signaling = pathways.show,
  geneLR.return = FALSE
)

# Row 2 corresponds to Pros1-Mertk.
LR.show <- pairLR.PROS[2, ]

print(LR.show)

pros_pair_file <- file.path(
  signaling_dir,
  "pair_PROS_WT_5XFAD_PS19.pdf"
)

check_output(pros_pair_file)

pdf(
  pros_pair_file,
  width = 8,
  height = 8
)

netVisual_individual(
  object.list[["WT"]],
  signaling = pathways.show,
  pairLR.use = LR.show,
  layout = "circle",
  arrow.size = 0.8,
  arrow.width = 2,
  remove.isolate = TRUE,
  edge.label.cex = 1.4,
  vertex.label.cex = 1.4,
  edge.weight.max = weight.max[1],
  edge.width.max = 10,
  signaling.name = "PROS WT"
)

netVisual_individual(
  object.list[["5XFAD;PS19"]],
  signaling = pathways.show,
  pairLR.use = LR.show,
  layout = "circle",
  arrow.size = 0.8,
  arrow.width = 2,
  remove.isolate = TRUE,
  edge.label.cex = 1.4,
  vertex.label.cex = 1.4,
  edge.weight.max = weight.max[1],
  edge.width.max = 10,
  signaling.name = "PROS 5XFAD;PS19"
)

dev.off()

# -----------------------------------------------------------------------------
# PROS signaling heatmap
# -----------------------------------------------------------------------------

mat_wt <- object.list[["WT"]]@netP$prob[
  ,
  ,
  pathways.show
]

mat_double <- object.list[["5XFAD;PS19"]]@netP$prob[
  ,
  ,
  pathways.show
]

mat_wt <- mat_wt[
  cluster_order,
  cluster_order
]

mat_double <- mat_double[
  cluster_order,
  cluster_order
]

# Use the same sender and receiver barplot scales for both genotypes.
shared_row_limit <- max(
  c(
    rowSums(abs(mat_wt)),
    rowSums(abs(mat_double))
  ),
  na.rm = TRUE
)

shared_col_limit <- max(
  c(
    colSums(abs(mat_wt)),
    colSums(abs(mat_double))
  ),
  na.rm = TRUE
)

ht_wt <- netVisual_heatmap_edit(
  object.list[["WT"]],
  signaling = pathways.show,
  ordered_levels = cluster_order,
  row.limit = shared_row_limit,
  col.limit = shared_col_limit,
  color.heatmap = "Reds",
  title.name = "PROS signaling WT"
)

ht_double <- netVisual_heatmap_edit(
  object.list[["5XFAD;PS19"]],
  signaling = pathways.show,
  ordered_levels = cluster_order,
  row.limit = shared_row_limit,
  col.limit = shared_col_limit,
  color.heatmap = "Reds",
  title.name = "PROS signaling 5XFAD;PS19"
)

pros_heatmap_file <- file.path(
  signaling_dir,
  "heat_PROS_WT_5XFAD_PS19.pdf"
)

check_output(pros_heatmap_file)

pdf(
  pros_heatmap_file,
  width = 12,
  height = 10
)

ComplexHeatmap::draw(
  ht_wt + ht_double,
  ht_gap = grid::unit(
    0.5,
    "cm"
  )
)

dev.off()


# =============================================================================
# VISTA signaling
# =============================================================================

pathways.show <- "VISTA"

# -----------------------------------------------------------------------------
# Ligand-receptor contribution
# -----------------------------------------------------------------------------

vista_contribution_file <- file.path(
  signaling_dir,
  "VISTA_LR_contribution_5XFAD_PS19.pdf"
)

check_output(vista_contribution_file)

p <- netAnalysis_contribution(
  object.list[["5XFAD;PS19"]],
  signaling = pathways.show
)

ggsave(
  vista_contribution_file,
  plot = p,
  width = 6,
  height = 3
)

# -----------------------------------------------------------------------------
# Vsir-Igsf11 ligand-receptor interaction
# -----------------------------------------------------------------------------

weight.max <- getMaxWeight(
  object.list,
  slot.name = "netP",
  attribute = pathways.show
)

pairLR.VISTA <- extractEnrichedLR(
  cellchat,
  signaling = pathways.show,
  geneLR.return = FALSE
)

# Row 1 corresponds to Vsir-Igsf11.
LR.show <- pairLR.VISTA[1, ]

print(LR.show)

vista_pair_file <- file.path(
  signaling_dir,
  "pair_VISTA_WT_5XFAD_PS19.pdf"
)

check_output(vista_pair_file)

pdf(
  vista_pair_file,
  width = 8,
  height = 8
)

netVisual_individual(
  object.list[["WT"]],
  signaling = pathways.show,
  pairLR.use = LR.show,
  layout = "circle",
  arrow.size = 0.8,
  arrow.width = 2,
  remove.isolate = TRUE,
  edge.label.cex = 1.4,
  vertex.label.cex = 1.4,
  edge.weight.max = weight.max[1],
  edge.width.max = 10,
  signaling.name = "VISTA WT"
)

netVisual_individual(
  object.list[["5XFAD;PS19"]],
  signaling = pathways.show,
  pairLR.use = LR.show,
  layout = "circle",
  arrow.size = 0.8,
  arrow.width = 2,
  remove.isolate = TRUE,
  edge.label.cex = 1.4,
  vertex.label.cex = 1.4,
  edge.weight.max = weight.max[1],
  edge.width.max = 10,
  signaling.name = "VISTA 5XFAD;PS19"
)

dev.off()

# -----------------------------------------------------------------------------
# VISTA signaling heatmap
# -----------------------------------------------------------------------------

mat_wt <- object.list[["WT"]]@netP$prob[
  ,
  ,
  pathways.show
]

mat_double <- object.list[["5XFAD;PS19"]]@netP$prob[
  ,
  ,
  pathways.show
]

mat_wt <- mat_wt[
  cluster_order,
  cluster_order
]

mat_double <- mat_double[
  cluster_order,
  cluster_order
]

# Use the same sender and receiver barplot scales for both genotypes.
shared_row_limit <- max(
  c(
    rowSums(abs(mat_wt)),
    rowSums(abs(mat_double))
  ),
  na.rm = TRUE
)

shared_col_limit <- max(
  c(
    colSums(abs(mat_wt)),
    colSums(abs(mat_double))
  ),
  na.rm = TRUE
)

ht_wt <- netVisual_heatmap_edit(
  object.list[["WT"]],
  signaling = pathways.show,
  ordered_levels = cluster_order,
  row.limit = shared_row_limit,
  col.limit = shared_col_limit,
  color.heatmap = "Reds",
  title.name = "VISTA signaling WT"
)

ht_double <- netVisual_heatmap_edit(
  object.list[["5XFAD;PS19"]],
  signaling = pathways.show,
  ordered_levels = cluster_order,
  row.limit = shared_row_limit,
  col.limit = shared_col_limit,
  color.heatmap = "Reds",
  title.name = "VISTA signaling 5XFAD;PS19"
)

vista_heatmap_file <- file.path(
  signaling_dir,
  "heat_VISTA_WT_5XFAD_PS19.pdf"
)

check_output(vista_heatmap_file)

pdf(
  vista_heatmap_file,
  width = 12,
  height = 10
)

ComplexHeatmap::draw(
  ht_wt + ht_double,
  ht_gap = grid::unit(
    0.5,
    "cm"
  )
)

dev.off()

# -----------------------------------------------------------------------------
# Save CellChat objects with communication-pattern results
# -----------------------------------------------------------------------------

saveRDS(
  object.list,
  object_list_file
)

saveRDS(
  cellchat,
  merged_file
)