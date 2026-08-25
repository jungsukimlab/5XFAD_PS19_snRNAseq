# =============================================================================
# 5XFAD;PS19 snRNA-seq analysis
# Step 8: DEG downstream analysis 
# =============================================================================
#
# Run from the project root directory.
#
# Inputs:
#   1. Step 7 MAST DEG tables:
#      publication_output/Step7_DEG/clust_DEGs/
#
#   2. Step 7 joined/log-normalized RNA object:
#      publication_output/Step7_DEG/
#          5XFAD_PS19_RNA_normalized.rds
#
# Outputs:
#   publication_output/Step8_DEG_downstream/
#       significant_DEGs/
#       overlap/
#       enrichment/
#
# Analyses:
#   1. Filter significant DEGs:
#          |avg_log2FC| >= 1.5
#          adjusted P value <= 0.05
#
#   2. Compare significant DEG sets among:
#          5XFAD vs WT
#          PS19 vs WT
#          5XFAD;PS19 vs WT
#      and generate UpSet plots.
#
#   3. Run Enrichr 

library(Seurat)
library(ComplexUpset)
library(enrichR)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

deg_input_dir <- file.path(
  "publication_output",
  "Step7_DEG",
  "clust_DEGs"
)

seurat_input_file <- file.path(
  "publication_output",
  "Step7_DEG",
  "5XFAD_PS19_RNA_normalized.rds"
)

output_dir <- file.path(
  "publication_output",
  "Step8_DEG_downstream"
)

sig_deg_dir <- file.path(
  output_dir,
  "significant_DEGs"
)

overlap_dir <- file.path(
  output_dir,
  "overlap"
)

enrichment_dir <- file.path(
  output_dir,
  "enrichment"
)

overlap_summary_file <- file.path(
  overlap_dir,
  "Step8_DEG_overlap_summary.csv"
)

gene_sets_file <- file.path(
  overlap_dir,
  "Step8_DEG_gene_sets.csv"
)

overwrite <- FALSE

dir.create(sig_deg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(overlap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)

check_output <- function(file) {
  if (file.exists(file) && !overwrite) {
    stop(
      "Output already exists: ", file,
      "\nSet overwrite <- TRUE only if you intend to replace it."
    )
  }
}


check_output(overlap_summary_file)
check_output(gene_sets_file)

# =============================================================================
# 1. Significant DEG filtering and genotype overlap
# =============================================================================

log2fc_cutoff <- 1.5
padj_cutoff <- 0.05

deg_files <- list.files(
  deg_input_dir,
  pattern = "^clust_DEG_.*_(5XFAD;PS19|5XFAD|PS19)_vs_WT\\.csv$",
  full.names = FALSE
)

if (length(deg_files) == 0) {
  stop(
    "No Step 7 DEG CSV files were found in: ",
    deg_input_dir
  )
}

# Extract annotated cluster names from the Step 7 filenames.
clusters <- sub(
  "^clust_DEG_",
  "",
  deg_files
)

clusters <- sub(
  "_(5XFAD;PS19|5XFAD|PS19)_vs_WT\\.csv$",
  "",
  clusters
)

clusters <- unique(
  clusters
)

overlap_summary <- data.frame()
deg_gene_sets <- data.frame(
  Cluster = character(),
  Set = character(),
  Gene = character(),
  stringsAsFactors = FALSE
)

# Store DEG gene lists for enrichment later in this script.
deg_sets <- list()

for (clust in clusters) {
  
  message(
    "Processing DEG overlap: ",
    clust
  )
  
  file_5XFAD <- file.path(
    deg_input_dir,
    paste0(
      "clust_DEG_",
      clust,
      "_5XFAD_vs_WT.csv"
    )
  )
  
  file_PS19 <- file.path(
    deg_input_dir,
    paste0(
      "clust_DEG_",
      clust,
      "_PS19_vs_WT.csv"
    )
  )
  
  file_double <- file.path(
    deg_input_dir,
    paste0(
      "clust_DEG_",
      clust,
      "_5XFAD;PS19_vs_WT.csv"
    )
  )
  
  if (
    !file.exists(file_5XFAD) ||
    !file.exists(file_PS19) ||
    !file.exists(file_double)
  ) {
    message(
      "Skipping ",
      clust,
      ": one or more Step 7 DEG files are missing."
    )
    next
  }
  
  data_5XFAD <- read.csv(
    file_5XFAD,
    row.names = 1,
    check.names = FALSE
  )
  
  data_PS19 <- read.csv(
    file_PS19,
    row.names = 1,
    check.names = FALSE
  )
  
  data_double <- read.csv(
    file_double,
    row.names = 1,
    check.names = FALSE
  )
  
  required_columns <- c(
    "avg_log2FC",
    "p_val_adj"
  )
  
  if (
    !all(required_columns %in% colnames(data_5XFAD)) ||
    !all(required_columns %in% colnames(data_PS19)) ||
    !all(required_columns %in% colnames(data_double))
  ) {
    stop(
      "Required DEG columns are missing for cluster ",
      clust,
      ". Expected avg_log2FC and p_val_adj."
    )
  }
  
  # ---------------------------------------------------------------------------
  # Significant DEGs
  # ---------------------------------------------------------------------------
  
  sig_5XFAD <- data_5XFAD[
    abs(data_5XFAD$avg_log2FC) >= log2fc_cutoff &
      data_5XFAD$p_val_adj <= padj_cutoff,
    ,
    drop = FALSE
  ]
  
  sig_PS19 <- data_PS19[
    abs(data_PS19$avg_log2FC) >= log2fc_cutoff &
      data_PS19$p_val_adj <= padj_cutoff,
    ,
    drop = FALSE
  ]
  
  sig_double <- data_double[
    abs(data_double$avg_log2FC) >= log2fc_cutoff &
      data_double$p_val_adj <= padj_cutoff,
    ,
    drop = FALSE
  ]
  
  sig_file_5XFAD <- file.path(
    sig_deg_dir,
    paste0(
      clust,
      "_5XFAD_vs_WT_significant_DEGs.csv"
    )
  )
  
  sig_file_PS19 <- file.path(
    sig_deg_dir,
    paste0(
      clust,
      "_PS19_vs_WT_significant_DEGs.csv"
    )
  )
  
  sig_file_double <- file.path(
    sig_deg_dir,
    paste0(
      clust,
      "_5XFAD_PS19_vs_WT_significant_DEGs.csv"
    )
  )
  
  check_output(sig_file_5XFAD)
  check_output(sig_file_PS19)
  check_output(sig_file_double)
  
  write.csv(
    sig_5XFAD,
    sig_file_5XFAD,
    row.names = TRUE
  )
  
  write.csv(
    sig_PS19,
    sig_file_PS19,
    row.names = TRUE
  )
  
  write.csv(
    sig_double,
    sig_file_double,
    row.names = TRUE
  )
  
  genes_5XFAD <- rownames(sig_5XFAD)
  genes_PS19 <- rownames(sig_PS19)
  genes_double <- rownames(sig_double)
  
  all_genes <- unique(
    c(
      genes_5XFAD,
      genes_PS19,
      genes_double
    )
  )
  
  upset_df <- data.frame(
    Gene = all_genes,
    `5XFAD` = all_genes %in% genes_5XFAD,
    `PS19` = all_genes %in% genes_PS19,
    `5XFAD;PS19` = all_genes %in% genes_double,
    check.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # Exclusive/shared DEG sets
  # ---------------------------------------------------------------------------
  
  genes_5XFAD_unique <- upset_df$Gene[
    upset_df$`5XFAD` &
      !upset_df$PS19 &
      !upset_df$`5XFAD;PS19`
  ]
  
  genes_PS19_unique <- upset_df$Gene[
    !upset_df$`5XFAD` &
      upset_df$PS19 &
      !upset_df$`5XFAD;PS19`
  ]
  
  genes_double_unique <- upset_df$Gene[
    !upset_df$`5XFAD` &
      !upset_df$PS19 &
      upset_df$`5XFAD;PS19`
  ]
  
  genes_5XFAD_double <- upset_df$Gene[
    upset_df$`5XFAD` &
      !upset_df$PS19 &
      upset_df$`5XFAD;PS19`
  ]
  
  genes_PS19_double <- upset_df$Gene[
    !upset_df$`5XFAD` &
      upset_df$PS19 &
      upset_df$`5XFAD;PS19`
  ]
  
  genes_all_three <- upset_df$Gene[
    upset_df$`5XFAD` &
      upset_df$PS19 &
      upset_df$`5XFAD;PS19`
  ]
  
  genes_5XFAD_PS19_only <- upset_df$Gene[
    upset_df$`5XFAD` &
      upset_df$PS19 &
      !upset_df$`5XFAD;PS19`
  ]
  
  # Store the DEG sets.
  deg_sets[[paste0(clust, "_5XFAD")]] <- genes_5XFAD
  deg_sets[[paste0(clust, "_PS19")]] <- genes_PS19
  deg_sets[[paste0(clust, "_5XFAD;PS19")]] <- genes_double
  deg_sets[[paste0(clust, "_5XFAD_unique")]] <- genes_5XFAD_unique
  deg_sets[[paste0(clust, "_PS19_unique")]] <- genes_PS19_unique
  deg_sets[[paste0(clust, "_5XFAD;PS19_unique")]] <- genes_double_unique
  deg_sets[[paste0(clust, "_5XFAD_and_5XFAD;PS19")]] <- genes_5XFAD_double
  deg_sets[[paste0(clust, "_PS19_and_5XFAD;PS19")]] <- genes_PS19_double
  
  # Export UpSet membership table.
  membership_file <- file.path(
    overlap_dir,
    paste0(
      clust,
      "_DEG_overlap_membership.csv"
    )
  )
  
  check_output(membership_file)
  
  write.csv(
    upset_df,
    membership_file,
    row.names = FALSE
  )
  
  # Long-format table of the DEG sets used for enrichment.
  current_sets <- list(
    `5XFAD;PS19` = genes_double,
    `5XFAD;PS19_unique` = genes_double_unique
  )
  
  for (set_name in names(current_sets)) {
    
    current_genes <- current_sets[[set_name]]
    
    if (length(current_genes) > 0) {
      deg_gene_sets <- rbind(
        deg_gene_sets,
        data.frame(
          Cluster = clust,
          Set = set_name,
          Gene = current_genes,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  overlap_summary <- rbind(
    overlap_summary,
    data.frame(
      Cluster = clust,
      n_5XFAD = length(genes_5XFAD),
      n_PS19 = length(genes_PS19),
      n_5XFAD_PS19 = length(genes_double),
      n_5XFAD_unique = length(genes_5XFAD_unique),
      n_PS19_unique = length(genes_PS19_unique),
      n_5XFAD_PS19_unique = length(genes_double_unique),
      n_5XFAD_and_double = length(genes_5XFAD_double),
      n_PS19_and_double = length(genes_PS19_double),
      n_all_three = length(genes_all_three),
      n_5XFAD_and_PS19_only = length(genes_5XFAD_PS19_only)
    )
  )
  
  # ---------------------------------------------------------------------------
  # UpSet plot
  # ---------------------------------------------------------------------------
  
  if (nrow(upset_df) > 0) {
    
    upset_file <- file.path(
      overlap_dir,
      paste0(
        clust,
        "_DEG_UpSet.pdf"
      )
    )
    
    check_output(upset_file)
    
    pdf(
      upset_file,
      width = 5,
      height = 6
    )
    
    print(
      upset(
        upset_df,
        c(
          "5XFAD;PS19",
          "5XFAD",
          "PS19"
        ),
        sort_intersections = FALSE,
        intersections = list(
          "5XFAD;PS19",
          "5XFAD",
          "PS19",
          c(
            "5XFAD",
            "5XFAD;PS19"
          ),
          c(
            "PS19",
            "5XFAD;PS19"
          ),
          c(
            "5XFAD",
            "PS19",
            "5XFAD;PS19"
          )
        )
      ) +
        theme(
          text = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 14)
        )
    )
    
    dev.off()
  }
}

write.csv(
  overlap_summary,
  overlap_summary_file,
  row.names = FALSE
)

write.csv(
  deg_gene_sets,
  gene_sets_file,
  row.names = FALSE
)

# =============================================================================
# 2. Enrichr analysis
# =============================================================================

setEnrichrSite(
  "Enrichr"
)

dbs <- c(
  "GO_Biological_Process_2023",
  "KEGG_2019_Mouse",
  "WikiPathways_2024_Mouse"
)

enrichment_sets <- c(
  "5XFAD;PS19",
  "5XFAD;PS19_unique"
)

# -----------------------------------------------------------------------------
# Edited plotEnrich function 
# -----------------------------------------------------------------------------

.enrichment_prep_df <- function(df, showTerms, orderBy) {
  
  if(is.null(showTerms)) {
    showTerms = nrow(df)
  } else if(!is.numeric(showTerms)) {
    stop(paste0("showTerms '", showTerms, "' is invalid."))
  }
  
  Annotated <- as.numeric(sub("^\\d+/", "", as.character(df$Overlap)))
  Significant <- as.numeric(sub("/\\d+$", "", as.character(df$Overlap)))
  
  # Build data frame
  df <- cbind(df, data.frame(Annotated = Annotated, Significant = Significant,
                             stringsAsFactors = FALSE))
  
  # Order data frame (P.value or Combined.Score)
  if(orderBy == "Combined.Score") {
    idx <- order(df$Combined.Score, decreasing = TRUE)
  } else {
    idx <- order(df$Adjusted.P.value, decreasing = FALSE)
  }
  df <- df[idx,]
  
  # Subset to selected number of terms
  if(showTerms <= nrow(df)) {
    df <- df[1:showTerms,]
  }
  
  return(df)
}

# Define the wrapText function to wrap long text
wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

plotEnrich_edit <- function(df, showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",
                            xlab = NULL, ylab = NULL, title = NULL) {
  if(!is.data.frame(df)) {
    stop("Input df is malformed - must be a data.frame object.")
  }
  if(nrow(df) == 0 | ncol(df) == 0) {
    stop("Input df is empty.")
  }
  if(!is.numeric(numChar)) {
    stop(paste0("numChar '", numChar, "' is invalid."))
  }
  
  df <- .enrichment_prep_df(df, showTerms, orderBy)
  # Transform Adjusted.P.value for y-axis
  df$Adjusted.P.value <- -log10(df$Adjusted.P.value)
  # Clean the Term names to remove unwanted parts like "(GO: xxxx)" or "WP xxx"
  df$Term <- gsub("\\s*\\(GO:\\s*\\d+\\)|\\s*WP\\d+", "", df$Term)
  # Wrap the Term column for better display
  df$wrap <- wrapText(df$Term, numChar)
  # Combine Term and Genes for a single label
  df$Label <- paste0(df$wrap, " (", df$Genes, ")")
  df$wrap2 <- wrapText(df$Label, numChar)
  
  df$Ratio <- df$Significant/df$Annotated
  
  df$wrap <- factor(df$wrap, levels = rev(df$wrap))
  
  # Define y variable (Count or Ratio)
  if(y != "Adjusted.P.value") {
    y <- "Significant"
  }
  
  # Define variable mapping
  map <- aes_string(x = "wrap", y = "Adjusted.P.value")
  
  # Define labels
  if(is.null(xlab)) {
    xlab <- "Enriched terms"
  }
  
  if(is.null(ylab)) {
    if(y == "Adjusted.P.value") {
      ylab <- "-log10(p.adjust)"
    } else {
      ylab <- "Gene count"
    }
  }
  
  if(is.null(title)) {
    title <- "Enrichment analysis by Enrichr"
  }
  
  # Make the ggplot
  p <- ggplot(df, map) + geom_bar(stat = "identity", fill = "lightgrey", width = .9) + 
    coord_flip() + theme_bw() 
  
  # Adjust theme components
  p <- p + theme(axis.text.x = element_text(size=16,colour = "black", vjust = 1),
                 axis.text.y = element_blank(),
                 axis.title = element_text(size=16,color = "black", margin = margin(10, 5, 0, 0)),
                 axis.title.y = element_text(size=16,angle = 90),
                 plot.title = element_text(size = 20)
  ) + 
    xlab(xlab) + ylab(ylab) + ggtitle(title) + 
    geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
    annotate("text", x = Inf, y = 1.3, label = "p.adjust = 0.05",
             hjust = -0.1, vjust = -0.5, color = "red") +
    geom_text(aes(label = wrap, y = 0),
              hjust = 0, 
              size = 8, 
              color = "black")  
  
  return(p)
}

# -----------------------------------------------------------------------------
# Run enrichment for every cluster x DEG-set 
# -----------------------------------------------------------------------------

for (clust in clusters) {
  
  for (set_name in enrichment_sets) {
    
    condition <- paste0(
      clust,
      "_",
      set_name
    )

    gene_list <- deg_sets[[condition]]
    
    if (length(gene_list) == 0) {
      message(
        "Skipping enrichment: ",
        condition,
        " has no genes."
      )
      next
    }
    
    message(
      "Running Enrichr: ",
      condition,
      " (",
      length(gene_list),
      " genes)"
    )
    
    enriched <- enrichr(
      gene_list,
      dbs
    )
    
    # One PDF per cluster/set combination, with one page per database.
    plot_file <- file.path(
      enrichment_dir,
      paste0(
        "EnrichR_",
        condition,
        ".pdf"
      )
    )
    
    check_output(
      plot_file
    )
    
    pdf(
      plot_file,
      width = 9,
      height = 11
    )
    
    for (db in names(enriched)) {
      
      enrichment_table <- as.data.frame(
        enriched[[db]]
      )
      
      result_file <- file.path(
        enrichment_dir,
        paste0(
          "EnrichR_",
          condition,
          "_",
          db,
          ".csv"
        )
      )
      
      check_output(
        result_file
      )
      
      write.csv(
        enrichment_table,
        result_file,
        row.names = FALSE
      )
      
      if (nrow(enrichment_table) > 0) {
        
        plot_obj <- plotEnrich_edit(
          enrichment_table,
          showTerms = 10,
          numChar = 50,
          y = "Adjusted.P.value",
          orderBy = "Adjusted.P.value",
          xlab = db,
          ylab = "-log10(p.adjust)",
          title = condition
        )
        
        print(
          plot_obj
        )
        
      } else {
        
        message(
          "No enrichment terms for ",
          condition,
          ": ",
          db
        )
      }
    }
    
    dev.off()
  }
}
