# =============================================================================
# 5XFAD;PS19 snRNA-seq analysis
# Step 12: Cross-species transcriptomic alignment
# =============================================================================
#
# Run from the project root directory.
#
# Inputs:
#   1. Step 11 pseudobulk DEG table:
#      publication_output/Step11_DEG_pseudobulk/
#          pseudo_ID_geno_DEG.csv
#
#   2. AMP-AD reference data (downloaded from Synapse on first run, cached
#      thereafter under publication_output/Step12_cross_species/source_data/):
#          syn11932957  consensus co-expression modules (Wan et al. 2020)
#          syn17010253  human-mouse ortholog table
#          syn14237651  AMP-AD differential expression (all 7 brain regions)
#          syn25428992  AD biological domain/subdomain GO term definitions
#          syn26856828  biological domain display labels/colors
#      Requires a Synapse personal access token with "view" + "download"
#      scope, set as the SYNAPSE_AUTH_TOKEN environment variable -- never
#      commit this token to git.
#
# Outputs:
#   publication_output/Step12_cross_species/
#       source_data/               cached Synapse downloads (see Inputs #2)
#       AMPAD_module_correlation/  mouse vs. AMP-AD module log2FC correlation
#       GSEA/                      mouse and human GO GSEA results
#       biodomain/                 NES-by-biodomain summaries and concordance
#
# Analysis:
#   1. AMP-AD module correlation
#   2. Biological domain enrichment
#   3. Cross-species functional alignment
#
# The cross-species workflow was adapted from teaching materials by
# Gregory Cary at the Jackson Workshop on Computational Techniques and Resources
# for Effective Translational Research in Alzheimer's Disease.

library(reticulate)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(forcats)
library(readr)
library(ggplot2)
library(ggrepel)
library(fgsea)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(GO.db)

theme_set(
  theme_bw()
)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

mouse_deg_file <- file.path(
  "publication_output",
  "Step11_DEG_pseudobulk",
  "pseudo_ID_geno_DEG.csv"
)

output_dir <- file.path(
  "publication_output",
  "Step12_cross_species"
)

source_data_dir <- file.path(
  output_dir,
  "source_data"
)

correlation_dir <- file.path(
  output_dir,
  "AMPAD_module_correlation"
)

gsea_dir <- file.path(
  output_dir,
  "GSEA"
)

biodomain_dir <- file.path(
  output_dir,
  "biodomain"
)

dir.create(
  source_data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  correlation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  gsea_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  biodomain_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(mouse_deg_file)) {
  stop(
    "Step 11 DEG file not found: ",
    mouse_deg_file
  )
}

overwrite <- FALSE

check_output <- function(file) {
  if (file.exists(file) && !overwrite) {
    stop(
      "Output already exists: ", file,
      "\nSet overwrite <- TRUE only if you intend to replace it."
    )
  }
}

# -----------------------------------------------------------------------------
# Output files
# -----------------------------------------------------------------------------

module_source_file <- file.path(
  source_data_dir,
  "AMPAD_consensus_modules_syn11932957.csv"
)

ortholog_source_file <- file.path(
  source_data_dir,
  "human_mouse_orthologs_syn17010253.tsv"
)

ampad_deg_source_file <- file.path(
  source_data_dir,
  "AMPAD_modules_raw_syn14237651.tsv"
)

biodomain_source_file <- file.path(
  source_data_dir,
  "AD_biodomain_definitions_syn25428992.rds"
)

biodomain_label_source_file <- file.path(
  source_data_dir,
  "AD_biodomain_labels_syn26856828.csv"
)

correlation_table_file <- file.path(
  correlation_dir,
  "AMPAD_module_correlations.csv"
)

correlation_plot_file <- file.path(
  correlation_dir,
  "AMPAD_module_correlations.pdf"
)

ifgturquoise_file <- file.path(
  correlation_dir,
  "IFGturquoise_5XFAD_PS19_correlation.pdf"
)

ifgbrown_file <- file.path(
  correlation_dir,
  "IFGbrown_5XFAD_PS19_correlation.pdf"
)

phggreen_file <- file.path(
  correlation_dir,
  "PHGgreen_5XFAD_PS19_correlation.pdf"
)

mouse_gsea_file <- file.path(
  gsea_dir,
  "mouse_GO_GSEA.rds"
)

human_gsea_file <- file.path(
  gsea_dir,
  "human_AMPAD_GO_GSEA.rds"
)

combined_biodomain_file <- file.path(
  biodomain_dir,
  "mouse_human_biodomain_NES.pdf"
)

human_region_biodomain_file <- file.path(
  biodomain_dir,
  "human_region_biodomain_NES.pdf"
)

concordant_terms_file <- file.path(
  biodomain_dir,
  "concordant_terms_by_biodomain.pdf"
)

concordant_terms_table_file <- file.path(
  biodomain_dir,
  "concordant_terms_by_biodomain.csv"
)

analysis_outputs <- c(
  module_source_file,
  ortholog_source_file,
  ampad_deg_source_file,
  biodomain_source_file,
  biodomain_label_source_file,
  correlation_table_file,
  correlation_plot_file,
  ifgturquoise_file,
  ifgbrown_file,
  phggreen_file,
  mouse_gsea_file,
  human_gsea_file,
  combined_biodomain_file,
  human_region_biodomain_file,
  concordant_terms_file,
  concordant_terms_table_file
)

invisible(
  lapply(
    analysis_outputs,
    check_output
  )
)

# -----------------------------------------------------------------------------
# Synapse login
# -----------------------------------------------------------------------------

synapse_token <- Sys.getenv(
  "SYNAPSE_AUTH_TOKEN"
)

if (identical(synapse_token, "")) {
  stop(
    "SYNAPSE_AUTH_TOKEN is not set. ",
    "Set a Synapse personal access token as an environment variable ",
    "before running Step 12."
  )
}

syn.client <- reticulate::import(
  "synapseclient"
)

syn <- syn.client$Synapse()

syn$login(
  authToken = synapse_token
)

# =============================================================================
# 1. AMP-AD module correlation
# =============================================================================

# -----------------------------------------------------------------------------
# AMP-AD consensus modules and mouse orthologs
# -----------------------------------------------------------------------------

query <- syn$tableQuery(
  "SELECT * FROM syn11932957"
)

module_table <- read_csv(
  query$filepath,
  show_col_types = FALSE
)

write_csv(
  module_table,
  module_source_file
)

mouse_human_ortho <- read_tsv(
  syn$get("syn17010253")$path,
  show_col_types = FALSE
)

write_tsv(
  mouse_human_ortho,
  ortholog_source_file
)

module_table$Mouse_gene_symbol <-
  mouse_human_ortho$mouse_symbol[
    match(
      module_table$GeneID,
      mouse_human_ortho$human_ensembl_gene
    )
  ]

ampad_modules <- module_table %>%
  distinct(
    tissue = brainRegion,
    module = Module,
    gene = GeneID,
    Mouse_gene_symbol
  ) %>%
  filter(
    !is.na(Mouse_gene_symbol),
    Mouse_gene_symbol != ""
  )

# -----------------------------------------------------------------------------
# AMP-AD AD-versus-control differential expression
# -----------------------------------------------------------------------------

ampad_modules_raw <- read_tsv(
  syn$get("syn14237651")$path,
  show_col_types = FALSE
)

write_tsv(
  ampad_modules_raw,
  ampad_deg_source_file
)

ampad_fc <- ampad_modules_raw %>%
  filter(
    Model == "Diagnosis",
    Comparison == "AD-CONTROL"
  ) %>%
  transmute(
    tissue = Tissue,
    gene = ensembl_gene_id,
    ampad_fc = logFC,
    ampad_adj.P = adj.P.Val
  ) %>%
  distinct()

ampad_modules_fc <- inner_join(
  ampad_modules,
  ampad_fc,
  by = c(
    "gene",
    "tissue"
  )
) %>%
  select(
    symbol = Mouse_gene_symbol,
    module,
    ampad_fc,
    ampad_adj.P
  )

# -----------------------------------------------------------------------------
# AMP-AD consensus-cluster assignments
# -----------------------------------------------------------------------------

cluster_a <- tibble(
  module = c(
    "TCXblue",
    "PHGyellow",
    "IFGyellow"
  ),
  cluster = "Consensus Cluster A (ECM organization)",
  cluster_label = "Consensus Cluster A\n(ECM organization)"
)

cluster_b <- tibble(
  module = c(
    "DLPFCblue",
    "CBEturquoise",
    "STGblue",
    "PHGturquoise",
    "IFGturquoise",
    "TCXturquoise",
    "FPturquoise"
  ),
  cluster = "Consensus Cluster B (Immune system)",
  cluster_label = "Consensus Cluster B\n(Immune system)"
)

cluster_c <- tibble(
  module = c(
    "IFGbrown",
    "STGbrown",
    "DLPFCyellow",
    "TCXgreen",
    "FPyellow",
    "CBEyellow",
    "PHGbrown"
  ),
  cluster = "Consensus Cluster C (Neuronal system)",
  cluster_label = "Consensus Cluster C\n(Neuronal system)"
)

cluster_d <- tibble(
  module = c(
    "DLPFCbrown",
    "STGyellow",
    "PHGgreen",
    "CBEbrown",
    "TCXyellow",
    "IFGblue",
    "FPblue"
  ),
  cluster = "Consensus Cluster D (Cell Cycle, NMD)",
  cluster_label = "Consensus Cluster D\n(Cell Cycle, NMD)"
)

cluster_e <- tibble(
  module = c(
    "FPbrown",
    "CBEblue",
    "DLPFCturquoise",
    "TCXbrown",
    "STGturquoise",
    "PHGblue"
  ),
  cluster = "Consensus Cluster E (Organelle Biogenesis, Cellular stress response)",
  cluster_label = "Consensus Cluster E\n(Organelle Biogenesis,\nCellular stress response)"
)

module_clusters <- bind_rows(
  cluster_a,
  cluster_b,
  cluster_c,
  cluster_d,
  cluster_e
) %>%
  mutate(
    cluster_label = fct_inorder(
      cluster_label
    )
  )

module_order <- module_clusters$module

# -----------------------------------------------------------------------------
# Pearson correlation of mouse and human log2 fold changes
# -----------------------------------------------------------------------------

fad.deg <- read.csv(
  mouse_deg_file,
  check.names = FALSE
) %>%
  rename(
    symbol = Gene
  )

required_mouse_columns <- c(
  "symbol",
  "Condition",
  "avg_log2FC"
)

missing_mouse_columns <- setdiff(
  required_mouse_columns,
  colnames(fad.deg)
)

if (length(missing_mouse_columns) > 0) {
  stop(
    "Step 11 DEG table is missing required columns: ",
    paste(
      missing_mouse_columns,
      collapse = ", "
    )
  )
}

model_vs_ampad <- inner_join(
  fad.deg,
  ampad_modules_fc,
  by = "symbol",
  multiple = "all"
)

cor_df <- model_vs_ampad %>%
  select(
    module,
    Condition,
    symbol,
    avg_log2FC,
    ampad_fc
  ) %>%
  group_by(
    module,
    Condition
  ) %>%
  nest(
    data = c(
      symbol,
      avg_log2FC,
      ampad_fc
    )
  ) %>%
  mutate(
    cor_test = map(
      data,
      ~ cor.test(
        .x$avg_log2FC,
        .x$ampad_fc,
        method = "pearson"
      )
    ),
    estimate = map_dbl(
      cor_test,
      "estimate"
    ),
    p_value = map_dbl(
      cor_test,
      "p.value"
    )
  ) %>%
  ungroup() %>%
  select(
    -cor_test
  )

correlation_for_plot <- cor_df %>%
  mutate(
    significant = p_value < 0.05
  ) %>%
  left_join(
    module_clusters,
    by = "module"
  ) %>%
  select(
    cluster,
    cluster_label,
    module,
    Condition,
    correlation = estimate,
    p_value,
    significant
  ) %>%
  arrange(
    cluster
  ) %>%
  mutate(
    module = factor(
      module,
      levels = module_order
    ),
    Condition = factor(
      Condition,
      levels = c(
        "PS19",
        "5XFAD",
        "5XFAD;PS19"
      )
    )
  )

write.csv(
  correlation_for_plot,
  correlation_table_file,
  row.names = FALSE
)

p_cor <- ggplot() +
  geom_tile(
    data = correlation_for_plot,
    aes(
      x = module,
      y = Condition
    ),
    colour = "black",
    fill = "white"
  ) +
  geom_point(
    data = correlation_for_plot,
    aes(
      x = module,
      y = Condition,
      colour = correlation,
      size = abs(correlation)
    )
  ) +
  geom_point(
    data = filter(
      correlation_for_plot,
      significant
    ),
    aes(
      x = module,
      y = Condition
    ),
    color = "black",
    shape = 0,
    size = 9
  ) +
  scale_x_discrete(
    position = "top"
  ) +
  scale_size_continuous(
    name = "Correlation strength",
    limits = c(
      0,
      max(
        abs(
          correlation_for_plot$correlation
        ),
        na.rm = TRUE
      )
    )
  ) +
  scale_color_gradient2(
    low = "#164B6E",
    mid = "white",
    high = "#85070C",
    midpoint = 0,
    name = "Correlation",
    guide = guide_colorbar(
      ticks = FALSE
    )
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  facet_grid(
    cols = vars(
      cluster_label
    ),
    scales = "free",
    space = "free"
  ) +
  theme(
    strip.text.x = element_text(
      size = 14,
      colour = "black"
    ),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 0,
      size = 14
    ),
    axis.text.y = element_text(
      size = 16
    ),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(
      size = 14
    ),
    legend.title = element_text(
      size = 15
    )
  )

ggsave(
  correlation_plot_file,
  p_cor,
  width = 18,
  height = 5
)

# -----------------------------------------------------------------------------
# Selected AMP-AD module correlations
# -----------------------------------------------------------------------------

# IFGturquoise
indiv_corr <- cor_df %>%
  filter(
    module == "IFGturquoise",
    Condition == "5XFAD;PS19"
  ) %>%
  unnest(
    data
  ) %>%
  mutate(
    facet = str_c(
      "r = ",
      signif(
        estimate,
        3
      ),
      " ; p = ",
      signif(
        p_value,
        3
      )
    )
  )

top_pos <- indiv_corr %>%
  filter(
    avg_log2FC > 0,
    ampad_fc > 0
  ) %>%
  arrange(
    desc(
      abs(
        avg_log2FC
      )
    )
  ) %>%
  slice_head(
    n = 5
  )

top_neg <- indiv_corr %>%
  filter(
    avg_log2FC < 0,
    ampad_fc < 0
  ) %>%
  arrange(
    desc(
      abs(
        avg_log2FC
      )
    )
  ) %>%
  slice_head(
    n = 1
  )

goi <- c(
  "Ly86",
  "Apbb1ip",
  "Inpp5d",
  "Cd84",
  "Ptprc",
  "Selplg",
  "Klhl6",
  "Ctsd",
  "Trem2",
  "Tgfbr2",
  "Chrna7"
)

label_df <- bind_rows(
  top_pos,
  top_neg,
  filter(
    indiv_corr,
    symbol %in% goi
  )
) %>%
  distinct(
    symbol,
    .keep_all = TRUE
  )

p <- ggplot(
  indiv_corr,
  aes(
    avg_log2FC,
    ampad_fc
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.1
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.1
  ) +
  geom_point(
    size = 0.4,
    color = "darkred"
  ) +
  geom_smooth(
    method = "lm",
    linewidth = 0.5
  ) +
  geom_text_repel(
    data = label_df,
    aes(
      label = symbol
    ),
    size = 4.5,
    min.segment.length = 0
  ) +
  labs(
    x = "5XFAD;PS19 Log2FC",
    y = "AMP-AD Log2FC",
    title = "IFGturquoise"
  ) +
  facet_wrap(
    ~ facet
  ) +
  theme_bw(
    base_size = 12
  )

ggsave(
  ifgturquoise_file,
  p,
  width = 5,
  height = 5
)

# IFGbrown
indiv_corr <- cor_df %>%
  filter(
    module == "IFGbrown",
    Condition == "5XFAD;PS19"
  ) %>%
  unnest(
    data
  ) %>%
  mutate(
    facet = str_c(
      "r = ",
      signif(
        estimate,
        3
      ),
      " ; p = ",
      signif(
        p_value,
        3
      )
    )
  )

top_pos <- indiv_corr %>%
  filter(
    avg_log2FC > 0,
    ampad_fc > 0
  ) %>%
  arrange(
    desc(
      abs(
        avg_log2FC
      )
    )
  ) %>%
  slice_head(
    n = 2
  )

top_neg <- indiv_corr %>%
  filter(
    avg_log2FC < 0,
    ampad_fc < 0
  ) %>%
  arrange(
    desc(
      abs(
        avg_log2FC
      )
    )
  ) %>%
  slice_head(
    n = 3
  )

goi <- c(
  "Arhgap24",
  "St3gal6",
  "Elk3",
  "Vgf",
  "Bdnf",
  "Neurod6",
  "Pcsk1",
  "Rph3a"
)

label_df <- bind_rows(
  top_pos,
  top_neg,
  filter(
    indiv_corr,
    symbol %in% goi
  )
) %>%
  distinct(
    symbol,
    .keep_all = TRUE
  )

p <- ggplot(
  indiv_corr,
  aes(
    avg_log2FC,
    ampad_fc
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.1
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.1
  ) +
  geom_point(
    size = 0.4,
    color = "darkred"
  ) +
  geom_smooth(
    method = "lm",
    linewidth = 0.5
  ) +
  geom_text_repel(
    data = label_df,
    aes(
      label = symbol
    ),
    size = 4.5,
    min.segment.length = 0
  ) +
  labs(
    x = "5XFAD;PS19 Log2FC",
    y = "AMP-AD Log2FC",
    title = "IFGbrown"
  ) +
  facet_wrap(
    ~ facet
  ) +
  theme_bw(
    base_size = 12
  )

ggsave(
  ifgbrown_file,
  p,
  width = 5,
  height = 5
)

# PHGgreen
indiv_corr <- cor_df %>%
  filter(
    module == "PHGgreen",
    Condition == "5XFAD;PS19"
  ) %>%
  unnest(
    data
  ) %>%
  mutate(
    facet = str_c(
      "r = ",
      signif(
        estimate,
        3
      ),
      " ; p = ",
      signif(
        p_value,
        3
      )
    )
  )

top_pos <- indiv_corr %>%
  filter(
    avg_log2FC > 0,
    ampad_fc > 0
  ) %>%
  arrange(
    desc(
      abs(
        avg_log2FC
      )
    )
  ) %>%
  slice_head(
    n = 5
  )

top_neg <- indiv_corr %>%
  filter(
    avg_log2FC < 0,
    ampad_fc < 0
  ) %>%
  arrange(
    desc(
      abs(
        avg_log2FC
      )
    )
  ) %>%
  slice_head(
    n = 3
  )

goi <- c(
  "Rasgrp3",
  "Cd9",
  "Myo1e",
  "Mog",
  "Mag",
  "Cldn11",
  "Trf",
  "Ugt8a",
  "Plp1",
  "Cdk19",
  "Nek7"
)

label_df <- bind_rows(
  top_pos,
  top_neg,
  filter(
    indiv_corr,
    symbol %in% goi
  )
) %>%
  distinct(
    symbol,
    .keep_all = TRUE
  )

p <- ggplot(
  indiv_corr,
  aes(
    avg_log2FC,
    ampad_fc
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.1
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.1
  ) +
  geom_point(
    size = 0.4,
    color = "darkred"
  ) +
  geom_smooth(
    method = "lm",
    linewidth = 0.5
  ) +
  geom_text_repel(
    data = label_df,
    aes(
      label = symbol
    ),
    size = 4.5,
    min.segment.length = 0
  ) +
  labs(
    x = "5XFAD;PS19 Log2FC",
    y = "AMP-AD Log2FC",
    title = "PHGgreen"
  ) +
  facet_wrap(
    ~ facet
  ) +
  theme_bw(
    base_size = 12
  )

ggsave(
  phggreen_file,
  p,
  width = 5,
  height = 5
)

# =============================================================================
# 2. Biological domain enrichment
# =============================================================================

# -----------------------------------------------------------------------------
# GO annotations
# -----------------------------------------------------------------------------

min_size <- 15
max_size <- 500

mouse_go_terms <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = keys(
    org.Mm.eg.db,
    keytype = "SYMBOL"
  ),
  columns = c(
    "GO",
    "SYMBOL"
  ),
  keytype = "SYMBOL"
)

mouse_go_sets <- split(
  mouse_go_terms$SYMBOL,
  mouse_go_terms$GO
)

mouse_go_sets <- mouse_go_sets[
  sapply(
    mouse_go_sets,
    function(x) {
      length(x) >= min_size &&
        length(x) <= max_size
    }
  )
]

human_go_terms <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(
    org.Hs.eg.db,
    keytype = "SYMBOL"
  ),
  columns = c(
    "GO",
    "SYMBOL"
  ),
  keytype = "SYMBOL"
)

human_go_sets <- split(
  human_go_terms$SYMBOL,
  human_go_terms$GO
)

human_go_sets <- human_go_sets[
  sapply(
    human_go_sets,
    function(x) {
      length(x) >= min_size &&
        length(x) <= max_size
    }
  )
]

# -----------------------------------------------------------------------------
# Mouse GSEA
# -----------------------------------------------------------------------------

fad.enr <- fad.deg %>%
  group_by(
    Condition
  ) %>%
  summarise(
    gl = list(
      sort(
        setNames(
          avg_log2FC,
          symbol
        ),
        decreasing = TRUE
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(
    gse = map(
      gl,
      ~ fgseaMultilevel(
        pathways = mouse_go_sets,
        stats = .x,
        minSize = min_size,
        maxSize = max_size,
        nproc = 1,
        nPermSimple = 100000
      )
    ),
    res = map(
      gse,
      ~ {
        go_term_map <- AnnotationDbi::select(
          GO.db,
          keys = .x$pathway,
          columns = c(
            "GOID",
            "TERM"
          ),
          keytype = "GOID"
        )

        inner_join(
          .x,
          go_term_map %>%
            select(
              pathway = GOID,
              TERM
            ),
          by = "pathway"
        ) %>%
          relocate(
            TERM,
            .after = pathway
          )
      }
    )
  )

saveRDS(
  fad.enr,
  mouse_gsea_file
)

# -----------------------------------------------------------------------------
# Human AMP-AD GSEA
# -----------------------------------------------------------------------------

hs.gsea <- ampad_modules_raw %>%
  filter(
    Model == "Diagnosis",
    Comparison == "AD-CONTROL",
    !is.na(hgnc_symbol)
  ) %>%
  group_by(
    Study,
    Tissue
  ) %>%
  summarise(
    gl = list(
      sort(
        setNames(
          logFC,
          hgnc_symbol
        ),
        decreasing = TRUE
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(
    gse = map(
      gl,
      ~ fgseaMultilevel(
        pathways = human_go_sets,
        stats = .x,
        minSize = min_size,
        maxSize = max_size,
        nproc = 1,
        nPermSimple = 100000
      )
    ),
    res = map(
      gse,
      ~ {
        go_term_map <- AnnotationDbi::select(
          GO.db,
          keys = .x$pathway,
          columns = c(
            "GOID",
            "TERM"
          ),
          keytype = "GOID"
        )

        inner_join(
          .x,
          go_term_map %>%
            select(
              pathway = GOID,
              TERM
            ),
          by = "pathway"
        ) %>%
          relocate(
            TERM,
            .after = pathway
          )
      }
    )
  )

saveRDS(
  hs.gsea,
  human_gsea_file
)

# -----------------------------------------------------------------------------
# AD biological-domain annotations
# -----------------------------------------------------------------------------

biodom <- readRDS(
  syn$get("syn25428992")$path
)

saveRDS(
  biodom,
  biodomain_source_file
)

dom.lab <- read_csv(
  syn$get("syn26856828")$path,
  show_col_types = FALSE
)

write_csv(
  dom.lab,
  biodomain_label_source_file
)

# =============================================================================
# 3. Cross-species functional alignment
# =============================================================================

# -----------------------------------------------------------------------------
# Mouse models and human AMP-AD studies
# -----------------------------------------------------------------------------

combined_nes <- bind_rows(
  fad.enr %>%
    mutate(
      model = as.character(
        Condition
      )
    ) %>%
    select(
      model,
      res
    ) %>%
    unnest(
      res
    ),
  hs.gsea %>%
    mutate(
      model = Study
    ) %>%
    select(
      model,
      res
    ) %>%
    unnest(
      res
    )
) %>%
  left_join(
    biodom %>%
      select(
        Biodomain,
        Subdomain,
        pathway = GO_ID
      ),
    by = "pathway"
  ) %>%
  full_join(
    dom.lab,
    by = c(
      "Biodomain" = "domain"
    )
  ) %>%
  filter(
    !is.na(model)
  ) %>%
  mutate(
    Biodomain = if_else(
      is.na(Biodomain),
      "none",
      Biodomain
    ),
    Subdomain = if_else(
      is.na(Subdomain),
      "none",
      Subdomain
    )
  ) %>%
  mutate(
    n_sig = length(
      unique(
        pathway
      )
    ),
    .by = Biodomain
  ) %>%
  mutate(
    Biodomain = fct_reorder(
      Biodomain,
      n_sig
    )
  ) %>%
  arrange(
    Biodomain,
    padj
  ) %>%
  mutate(
    model = factor(
      model,
      levels = c(
        "5XFAD",
        "PS19",
        "5XFAD;PS19",
        "MAYO",
        "MSSM",
        "ROSMAP"
      )
    )
  )

p <- ggplot(
  combined_nes,
  aes(
    NES,
    Biodomain
  )
) +
  facet_wrap(
    ~ model,
    nrow = 1
  ) +
  geom_violin(
    data = subset(
      combined_nes,
      NES > 0
    ),
    aes(
      color = color
    ),
    scale = "width"
  ) +
  geom_violin(
    data = subset(
      combined_nes,
      NES < 0
    ),
    aes(
      color = color
    ),
    scale = "width"
  ) +
  geom_jitter(
    aes(
      size = -log10(
        padj
      ),
      fill = color
    ),
    color = "grey20",
    shape = 21,
    alpha = 0.3
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.1
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(
    limits = c(
      -2.5,
      2.5
    ),
    oob = scales::oob_squish
  ) +
  scale_size_continuous(
    name = expression(
      -log[10](
        padj
      )
    ),
    limits = c(
      0,
      20
    ),
    range = c(
      1,
      6
    )
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  combined_biodomain_file,
  p,
  width = 11.5,
  height = 4.5
)

# -----------------------------------------------------------------------------
# Human AMP-AD brain regions
# -----------------------------------------------------------------------------

human_region_levels <- c(
  "ROSMAP, DLPFC",
  "MSSM, FP",
  "MSSM, STG",
  "MSSM, PHG",
  "MSSM, IFG",
  "MAYO, CBE",
  "MAYO, TCX"
)

human_region_nes <- hs.gsea %>%
  mutate(
    model = str_c(
      Study,
      ", ",
      Tissue
    )
  ) %>%
  select(
    model,
    res
  ) %>%
  unnest(
    res
  ) %>%
  left_join(
    biodom %>%
      select(
        Biodomain,
        Subdomain,
        pathway = GO_ID
      ),
    by = "pathway"
  ) %>%
  full_join(
    dom.lab,
    by = c(
      "Biodomain" = "domain"
    )
  ) %>%
  filter(
    !is.na(model)
  ) %>%
  mutate(
    Biodomain = if_else(
      is.na(Biodomain),
      "none",
      Biodomain
    ),
    Subdomain = if_else(
      is.na(Subdomain),
      "none",
      Subdomain
    )
  ) %>%
  mutate(
    n_sig = length(
      unique(
        pathway
      )
    ),
    .by = Biodomain
  ) %>%
  mutate(
    Biodomain = fct_reorder(
      Biodomain,
      n_sig
    )
  ) %>%
  arrange(
    Biodomain,
    padj
  ) %>%
  mutate(
    model = factor(
      model,
      levels = human_region_levels
    )
  )

p <- ggplot(
  human_region_nes,
  aes(
    NES,
    Biodomain
  )
) +
  facet_wrap(
    ~ model,
    nrow = 1
  ) +
  geom_violin(
    data = subset(
      human_region_nes,
      NES > 0
    ),
    aes(
      color = color
    ),
    scale = "width"
  ) +
  geom_violin(
    data = subset(
      human_region_nes,
      NES < 0
    ),
    aes(
      color = color
    ),
    scale = "width"
  ) +
  geom_jitter(
    aes(
      size = -log10(
        padj
      ),
      fill = color
    ),
    color = "grey20",
    shape = 21,
    alpha = 0.3
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.1
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(
    limits = c(
      -2.5,
      2.5
    ),
    oob = scales::oob_squish
  ) +
  scale_size_continuous(
    name = expression(
      -log[10](
        padj
      )
    ),
    limits = c(
      0,
      20
    ),
    range = c(
      1,
      6
    )
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  human_region_biodomain_file,
  p,
  width = 11.5,
  height = 4.5
)

# -----------------------------------------------------------------------------
# Concordant GO terms between mouse models and human AMP-AD
# -----------------------------------------------------------------------------

nes_tol <- 0.05

terms_tbl <- combined_nes %>%
  select(
    Biodomain,
    model,
    TERM,
    NES
  ) %>%
  mutate(
    dir_NES = case_when(
      is.na(NES) ~ NA_integer_,
      NES > nes_tol ~ 1L,
      NES < -nes_tol ~ -1L,
      TRUE ~ 0L
    )
  )

human_summary <- terms_tbl %>%
  filter(
    model %in% c(
      "ROSMAP",
      "MSSM",
      "MAYO"
    )
  ) %>%
  group_by(
    Biodomain,
    TERM
  ) %>%
  summarise(
    dir_mean = mean(
      dir_NES,
      na.rm = TRUE
    ),
    dir_human = case_when(
      dir_mean > 0.05 ~ 1L,
      dir_mean < -0.05 ~ -1L,
      TRUE ~ 0L
    ),
    .groups = "drop"
  )

mouse_summary <- terms_tbl %>%
  filter(
    model %in% c(
      "5XFAD",
      "PS19",
      "5XFAD;PS19"
    )
  ) %>%
  left_join(
    human_summary,
    by = c(
      "Biodomain",
      "TERM"
    )
  ) %>%
  mutate(
    same_dir = case_when(
      dir_human == 1L &
        dir_NES == 1L ~ TRUE,
      dir_human == -1L &
        dir_NES == -1L ~ TRUE,
      TRUE ~ FALSE
    )
  )

same_dir_true_bd <- mouse_summary %>%
  filter(
    same_dir,
    Biodomain != "none"
  ) %>%
  group_by(
    Biodomain,
    model
  ) %>%
  summarise(
    n_true = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Biodomain = factor(
      Biodomain,
      levels = levels(
        combined_nes$Biodomain
      )
    ),
    model = factor(
      model,
      levels = c(
        "5XFAD",
        "PS19",
        "5XFAD;PS19"
      )
    )
  )

write.csv(
  same_dir_true_bd,
  concordant_terms_table_file,
  row.names = FALSE
)

p <- ggplot(
  same_dir_true_bd,
  aes(
    x = Biodomain,
    y = n_true,
    fill = model
  )
) +
  geom_col(
    position = "dodge"
  ) +
  scale_fill_manual(
    values = c(
      "5XFAD" = "darkgoldenrod1",
      "PS19" = "darkslategray3",
      "5XFAD;PS19" = "palevioletred"
    )
  ) +
  labs(
    x = "Biodomain",
    y = "Number of concordant terms",
    fill = "Genotype"
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  concordant_terms_file,
  p,
  width = 8,
  height = 5
)
