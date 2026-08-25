# =============================================================================
# 5XFAD;PS19 snRNA-seq analysis
# Step 9: hdWGCNA co-expression network analysis
# =============================================================================
#
# Run from the project root directory.
#
# Input:
#   publication_output/Step7_DEG/
#       5XFAD_PS19_RNA_normalized.rds
#
# Outputs:
#   publication_output/Step9_hdWGCNA/
#
# Analysis:
#   1. Set up hdWGCNA
#   2. Module eigengenes and module connectivity
#   3. Enrichr analysis
#   4. Differential module eigengene (DME) analysis

library(Seurat)
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)

theme_set(
  theme_cowplot()
)

set.seed(
  12345
)

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
  "Step9_hdWGCNA"
)

network_dir <- file.path(
  output_dir,
  "network"
)

enrichment_dir <- file.path(
  output_dir,
  "enrichment"
)

enrichment_plot_dir <- file.path(
  enrichment_dir,
  "plots"
)

dme_dir <- file.path(
  output_dir,
  "DME"
)

tom_dir <- file.path(
  output_dir,
  "TOM"
)

output_rds <- file.path(
  output_dir,
  "5XFAD_PS19_hdWGCNA.rds"
)

soft_power_plot_file <- file.path(
  network_dir,
  "soft_power.pdf"
)

soft_power_table_file <- file.path(
  network_dir,
  "soft_power_table.csv"
)

dendrogram_file <- file.path(
  network_dir,
  "dendrogram.pdf"
)

hmes_file <- file.path(
  network_dir,
  "harmonized_module_eigengenes.csv"
)

modules_file <- file.path(
  network_dir,
  "module_assignments.csv"
)

hub_genes_file <- file.path(
  network_dir,
  "hub_genes.csv"
)

hme_featureplot_file <- file.path(
  network_dir,
  "module_eigengene_featureplots.pdf"
)

cell_type_umap_file <- file.path(
  network_dir,
  "cell_type_umap.pdf"
)

module_dotplot_file <- file.path(
  network_dir,
  "module_eigengene_dotplot_cell_type.pdf"
)

enrichment_table_file <- file.path(
  enrichment_dir,
  "module_Enrichr_results.csv"
)

dme_file <- file.path(
  dme_dir,
  "DME_results.csv"
)

dme_filtered_file <- file.path(
  dme_dir,
  "DME_results_filtered.csv"
)

dme_plot_file <- file.path(
  dme_dir,
  "DME_heatmap.pdf"
)

overwrite <- FALSE

dir.create(
  network_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  enrichment_plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dme_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tom_dir,
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
  output_rds,
  soft_power_plot_file,
  soft_power_table_file,
  dendrogram_file,
  hmes_file,
  modules_file,
  hub_genes_file,
  hme_featureplot_file,
  cell_type_umap_file,
  module_dotplot_file,
  enrichment_table_file,
  dme_file,
  dme_filtered_file,
  dme_plot_file
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

# =============================================================================
# 1. Load object and set up hdWGCNA
# =============================================================================

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
    paste(
      missing_metadata,
      collapse = ", "
    )
  )
}

if (!"RNA" %in% Assays(obj)) {
  stop(
    "RNA assay not found in the Seurat object."
  )
}

if (!"harmony" %in% Reductions(obj)) {
  stop(
    "Harmony reduction not found in the Seurat object."
  )
}

if (!"umap" %in% Reductions(obj)) {
  stop(
    "UMAP reduction not found in the Seurat object."
  )
}

DefaultAssay(obj) <- "RNA"

obj$Genotype <- factor(
  as.character(obj$Genotype),
  levels = c(
    "WT",
    "5XFAD",
    "PS19",
    "5XFAD;PS19"
  )
)

Idents(obj) <- "clust_cell"

wgcna_name <- "5XFAD_PS19_clust"

# Use genes expressed in at least the specified fraction of nuclei.
obj <- SetupForWGCNA(
  obj,
  gene_select = "fraction",
  fraction = 0.0000001,
  wgcna_name = wgcna_name
)

# -----------------------------------------------------------------------------
# Construct metacells within annotated cluster x biological sample
# -----------------------------------------------------------------------------

obj <- MetacellsByGroups(
  obj,
  group.by = c(
    "clust_cell",
    "ID"
  ),
  reduction = "harmony",
  k = 35, # The number of cells to be aggregated (20~75). 45 & 55 got errors due to low cell number in EC.
  min_cells = 100,
  max_shared = 10,
  ident.group = "clust_cell",
  assay = "RNA",
  layer = "counts",
  wgcna_name = wgcna_name
)

obj <- NormalizeMetacells(
  obj,
  wgcna_name = wgcna_name
)

# Use all annotated-cluster metacells to construct one global network.
obj <- SetDatExpr(
  obj,
  group_name = NULL,
  group.by = NULL, # If NULL (default), hdWGCNA uses the Seurat Idents as the group.
  use_metacells = TRUE,
  assay = "RNA",
  layer = "data",
  wgcna_name = wgcna_name
)

# -----------------------------------------------------------------------------
# Soft-power selection
# -----------------------------------------------------------------------------

obj <- TestSoftPowers(
  obj,
  networkType = "signed",
  wgcna_name = wgcna_name
)

soft_power_plots <- PlotSoftPowers(
  obj,
  wgcna_name = wgcna_name
)

pdf(
  soft_power_plot_file,
  width = 10,
  height = 8
)

print(
  wrap_plots(
    soft_power_plots,
    ncol = 2
  )
)

dev.off()

power_table <- GetPowerTable(
  obj,
  wgcna_name = wgcna_name
)

write.csv(
  power_table,
  soft_power_table_file,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Construct co-expression network
# -----------------------------------------------------------------------------

obj <- ConstructNetwork(
  obj,
  tom_outdir = tom_dir,
  tom_name = wgcna_name,
  overwrite_tom = TRUE
)

pdf(
  dendrogram_file,
  width = 10,
  height = 6
)

PlotDendrogram(
  obj,
  main = "5XFAD;PS19 hdWGCNA dendrogram",
  wgcna_name = wgcna_name
)

dev.off()

# =============================================================================
# 2. Module eigengenes and connectivity
# =============================================================================
#
# Module eigengenes are harmonized by biological sample ID.

obj <- ModuleEigengenes(
  obj,
  group.by.vars = "ID",
  assay = "RNA",
  wgcna_name = wgcna_name
)

# Compute eigengene-based connectivity across the same annotated populations
# used to build the network.
obj <- ModuleConnectivity(
  obj,
  group.by = "clust_cell",
  harmonized = TRUE,
  assay = "RNA",
  layer = "data",
  wgcna_name = wgcna_name
)

# Rename modules Mod1, Mod2, ...
obj <- ResetModuleNames(
  obj,
  new_name = "Mod",
  wgcna_name = wgcna_name
)

hMEs <- GetMEs(
  obj,
  harmonized = TRUE,
  wgcna_name = wgcna_name
)

write.csv(
  hMEs,
  hmes_file,
  row.names = TRUE
)

modules <- GetModules(
  obj,
  wgcna_name = wgcna_name
)

write.csv(
  modules,
  modules_file,
  row.names = FALSE
)

mods <- levels(
  modules$module
)

mods <- mods[
  mods != "grey"
]

mods <- gtools::mixedsort(
  mods
)

hub_df <- GetHubGenes(
  obj,
  n_hubs = 50,
  wgcna_name = wgcna_name
)

write.csv(
  hub_df,
  hub_genes_file,
  row.names = FALSE
)

# Add harmonized module eigengenes to Seurat metadata for plotting.
obj <- AddMetaData(
  obj,
  metadata = hMEs
)

# -----------------------------------------------------------------------------
# Module eigengenes on the Harmony-derived UMAP
# -----------------------------------------------------------------------------

hme_plots <- ModuleFeaturePlot(
  obj,
  features = "hMEs",
  module_names = mods,
  order_points = TRUE,
  reduction = "umap",
  raster = FALSE,
  wgcna_name = wgcna_name
)

hme_plots <- lapply(
  hme_plots,
  function(p) {
    p +
      theme(
        plot.title = element_text(
          size = 9
        )
      )
  }
)

pdf(
  hme_featureplot_file,
  width = 12,
  height = 10
)

print(
  wrap_plots(
    hme_plots,
    ncol = 4
  )
)

dev.off()

# -----------------------------------------------------------------------------
# Broad cell-type UMAP
# -----------------------------------------------------------------------------

pdf(
  cell_type_umap_file,
  width = 7,
  height = 6
)

print(
  DimPlot(
    obj,
    reduction = "umap",
    group.by = "cell_type",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  )
)

dev.off()

# -----------------------------------------------------------------------------
# Module eigengene dot plot by cluster and cell type
# -----------------------------------------------------------------------------

# Module eigengene dot plot by annotated cluster
obj$clust_cell <- factor(
  obj$clust_cell,
  levels = c(
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
)

pdf(
  file.path(
    network_dir,
    "module_eigengene_dotplot_cluster.pdf"
  ),
  width = 11,
  height = 8
)

p <- DotPlot(
  obj,
  features = mods,
  group.by = "clust_cell",
  dot.scale = 5
) +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(
    high = "red",
    mid = "grey95",
    low = "blue"
  ) +
  theme(
    axis.text.x = element_text(size = 8),
    panel.grid.major = element_line(
      color = "grey80",
      linewidth = 0.5
    ),
    panel.grid.minor = element_line(
      color = "grey90",
      linewidth = 0.2
    )
  )

print(p)

dev.off()


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

obj$cell_type <- factor(
  as.character(obj$cell_type),
  levels = cell_type_order
)

pdf(
  module_dotplot_file,
  width = 9,
  height = 7
)

p <- DotPlot(
  obj,
  features = mods,
  group.by = "cell_type",
  dot.scale = 6
) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(
    x = "Cell type",
    y = "Module",
    color = "Average expression",
    size = "Percent expressed"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.4
    )
  )

print(
  p
)

dev.off()

# =============================================================================
# 3. hdWGCNA Enrichr analysis
# =============================================================================
#
# Enrichment is performed with the hdWGCNA RunEnrichr pipeline using the top
# 100 genes in each module ranked by module connectivity.

setEnrichrSite(
  "Enrichr"
)

dbs <- c(
  "GO_Biological_Process_2023",
  "KEGG_2019_Mouse",
  "WikiPathways_2024_Mouse"
)

obj <- RunEnrichr(
  obj,
  dbs = dbs,
  max_genes = 100,
  wgcna_name = wgcna_name
)

enrichr_df <- GetEnrichrTable(
  obj,
  wgcna_name = wgcna_name
)

write.csv(
  enrichr_df,
  enrichment_table_file,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Enrichment plots for each module
# -----------------------------------------------------------------------------
#
# Top terms are ranked by adjusted P value. Bar length shows the Enrichr
# Combined Score and fill shows -log10(adjusted P value).

EnrichrBarPlot_edit<- function (obj, outdir = "enrichr_plots", n_terms = 10,
                                plot_size = c(10, 10), logscale = FALSE, plot_bar_color = NULL,
                                plot_text_color = NULL, wgcna_name = NULL, ...)
{
  if (is.null(wgcna_name)) {
    wgcna_name <- obj@misc$active_wgcna
  }
  modules <- GetModules(obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != "grey"]
  enrichr_df <- GetEnrichrTable(obj, wgcna_name)
  wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"),
           USE.NAMES = FALSE)
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  for (i in 1:length(mods)) {
    cur_mod <- mods[i]
    cur_terms <- subset(enrichr_df, module == cur_mod)
    print(cur_mod)
    cur_color <- modules %>% subset(module == cur_mod) %>%
      .$color %>% unique %>% as.character
    if (!is.null(plot_bar_color)) {
      cur_color <- plot_bar_color
    }
    if (nrow(cur_terms) == 0) {
      next
    }
    cur_terms$wrap <- wrapText(cur_terms$Term, 45)
    plot_list <- list()
    for (cur_db in dbs) {
      plot_df <- subset(cur_terms, db == cur_db) %>% 
        arrange(Adjusted.P.value) %>%  # Sort by Adjusted.P.value in ascending order
        slice_head(n = n_terms)         # Get the top n_terms with the lowest Adjusted.P.value
      if (is.null(plot_text_color)) {
        if (cur_color %in% c("black","brown","blue")) {
          text_color = "grey"
        }
        else {
          text_color = "black"
        }
      }
      else {
        text_color <- plot_text_color
      }
      if (logscale) {
        plot_df$Adjusted.P.value <- -log10(plot_df$Adjusted.P.value)
        lab <- "-log10(p.adjust)"
        
      }
      else {
        lab <- "p.adjust"
        
      }
      plot_list[[cur_db]] <- ggplot(plot_df, aes(x = Adjusted.P.value,
                                                 y = reorder(wrap, Adjusted.P.value))) + geom_bar(stat = "identity",
                                                                                                  position = "identity", color = "white", fill = cur_color, width = 0.8) +
        geom_text(aes(label = wrap, x = 0), color = text_color,
                  size = 5, hjust = 0) + ylab(cur_mod) +
        xlab(lab) + ggtitle(cur_db) + theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), legend.title = element_blank(),
                                            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                                            plot.title = element_text(hjust = 0.5))
    }
    pdf(paste0(outdir, "/", cur_mod, ".pdf"), width = plot_size[1],
        height = plot_size[2])
    for (plot in plot_list) {
      print(plot)
    }
    dev.off()
  }
}

EnrichrBarPlot_edit(
  obj,
  outdir = enrichment_plot_dir,
  n_terms = 10,
  plot_size = c(7, 10),
  logscale = TRUE,
  wgcna_name = wgcna_name
)


# =============================================================================
# 4. Differential module eigengene analysis
# =============================================================================
#
# Compare each disease genotype with WT within every annotated cluster.

if (is.factor(obj$clust_cell)) {
  cluster_order <- levels(
    obj$clust_cell
  )
} else {
  cluster_order <- unique(
    as.character(obj$clust_cell)
  )
}

dme_genotypes <- c(
  "5XFAD",
  "PS19",
  "5XFAD;PS19"
)

DMEs <- data.frame()

for (cur_genotype in dme_genotypes) {
  
  for (cur_cluster in cluster_order) {
    
    group1 <- rownames(
      obj[[]]
    )[
      as.character(obj$clust_cell) == cur_cluster &
        as.character(obj$Genotype) == cur_genotype
    ]
    
    group2 <- rownames(
      obj[[]]
    )[
      as.character(obj$clust_cell) == cur_cluster &
        as.character(obj$Genotype) == "WT"
    ]
    
    if (
      length(group1) == 0 ||
      length(group2) == 0
    ) {
      message(
        "Skipping DME: ",
        cur_cluster,
        " | ",
        cur_genotype,
        " vs WT because one group has no nuclei."
      )
      next
    }
    
    message(
      "Running DME: ",
      cur_cluster,
      " | ",
      cur_genotype,
      " vs WT"
    )
    
    cur_DMEs <- FindDMEs(
      obj,
      barcodes1 = group1,
      barcodes2 = group2,
      test.use = "wilcox",
      wgcna_name = wgcna_name
    )
    
    cur_DMEs$Cluster <- cur_cluster
    cur_DMEs$Genotype <- cur_genotype
    
    rownames(
      cur_DMEs
    ) <- NULL
    
    DMEs <- rbind(
      DMEs,
      cur_DMEs
    )
  }
}

write.csv(
  DMEs,
  dme_file,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Filter modules with positive average expression in each annotated cluster
# -----------------------------------------------------------------------------

dot_data <- DotPlot(
  obj,
  features = mods,
  group.by = "clust_cell",
  dot.scale = 5
)

avg_expr_per_module <- dot_data$data %>%
  select(
    id,
    features.plot,
    avg.exp.scaled
  )

DMEs_filtered <- DMEs %>%
  left_join(
    avg_expr_per_module,
    by = c(
      "Cluster" = "id",
      "module" = "features.plot"
    )
  ) %>%
  filter(
    avg.exp.scaled > 0,
    pct.1 > 0
  )

write.csv(
  DMEs_filtered,
  dme_filtered_file,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# DME heatmaps
# -----------------------------------------------------------------------------

plot_df <- DMEs_filtered

plot_df$module <- factor(
  as.character(plot_df$module),
  levels = mods
)

plot_df$Cluster <- factor(
  as.character(plot_df$Cluster),
  levels = cluster_order
)

plot_df$Genotype <- factor(
  as.character(plot_df$Genotype),
  levels = dme_genotypes
)

maxval <- 10
minval <- -10

plot_df$avg_log2FC <- pmax(
  pmin(
    plot_df$avg_log2FC,
    maxval
  ),
  minval
)

plot_df$Significance <- gtools::stars.pval(
  plot_df$p_val_adj
)

dme_plots <- list()

for (cur_genotype in dme_genotypes) {
  
  genotype_plot_df <- plot_df %>%
    filter(
      Genotype == cur_genotype
    ) %>%
    mutate(
      text_color = ifelse(
        abs(avg_log2FC) >= 1.5,
        "white",
        "black"
      ),
      avg_log2FC_plot = ifelse(
        p_val_adj <= 0.05,
        avg_log2FC,
        NA
      )
    )
  
  p <- ggplot(
    genotype_plot_df,
    aes(
      y = Cluster,
      x = module,
      fill = avg_log2FC_plot
    )
  ) +
    geom_tile() +
    geom_text(
      aes(
        label = ifelse(
          !is.na(avg_log2FC_plot),
          Significance,
          ""
        ),
        color = text_color
      ),
      size = 3
    ) +
    scale_fill_gradientn(
      colors = c(
        "blue",
        "white",
        "red"
      ),
      values = scales::rescale(
        c(
          minval,
          0,
          maxval
        )
      ),
      limits = c(
        minval,
        maxval
      ),
      na.value = "white",
      breaks = c(
        minval,
        0,
        maxval
      ),
      labels = c(
        paste0(
          "<=",
          minval
        ),
        "0",
        paste0(
          ">=",
          maxval
        )
      )
    ) +
    scale_color_identity() +
    labs(
      x = NULL,
      y = NULL,
      title = paste0(
        cur_genotype,
        " vs WT"
      ),
      fill = "Average log2FC\n(capped)"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      panel.border = element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.5
      )
    ) +
    coord_equal()
  
  dme_plots[[cur_genotype]] <- p
}

combined_dme_plot <- (
  dme_plots[["5XFAD"]] +
    theme(
      legend.position = "none"
    )
) +
  (
    dme_plots[["PS19"]] +
      theme(
        legend.position = "none"
      )
  ) +
  dme_plots[["5XFAD;PS19"]]

pdf(
  dme_plot_file,
  width = 20,
  height = 12
)

print(
  combined_dme_plot
)

dev.off()

# =============================================================================
# 5. Save final hdWGCNA object
# =============================================================================

saveRDS(
  obj,
  output_rds
)
