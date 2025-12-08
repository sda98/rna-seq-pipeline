#!/usr/bin/env Rscript
# 05A_pca.R
# PCA and sample QC plots from a DESeq2 vst-transformed object (vsd)
#
# Usage (default paths):
#   Rscript R/05A_pca.R
#
# Or with custom paths:
#   Rscript R/05A_pca.R results/deseq2/vsd.rds results/figures
#
# Expects:
#   - vsd.rds: DESeqTransform object with assay(vsd) = transformed counts
#   - colData(vsd): sample metadata (e.g., condition, batch, etc.)
#
# Outputs:
#   - results/figures/pca_PC1_PC2.png
#   - results/figures/pca_PC3_PC4.png
#   - results/figures/pca_PC1_PC2_PC3_3D.html
#   - results/figures/sample_correlation_heatmap.png
#   - results/deseq2/pca_top_genes_PC1_PC2.csv

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(matrixStats)
  library(plotly)
  library(htmlwidgets)
})

# ----------------------------
# 1. Config & input arguments
# ----------------------------

args <- commandArgs(trailingOnly = TRUE)

vsd_file <- ifelse(length(args) >= 1, args[1], "results/deseq2/vsd.rds")
fig_dir  <- ifelse(length(args) >= 2, args[2], "results/figures")
out_dir  <- "results/deseq2"

if (!file.exists(vsd_file)) {
  stop("vsd file not found: ", vsd_file,
       "\nMake sure 01_deseq2_analysis.R saves vsd to this path.")
}

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Loading vsd from: ", vsd_file)
vsd <- readRDS(vsd_file)

# ----------------------------
# 2. Prepare data for PCA
# ----------------------------

mat <- assay(vsd)
if (!is.matrix(mat)) {
  stop("assay(vsd) is not a matrix.")
}

# Use top variable genes (up to 500)
rv <- matrixStats::rowVars(mat)
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
mat_top <- mat[select, , drop = FALSE]

# PCA
pca <- prcomp(t(mat_top), scale. = FALSE)
percent_var <- (pca$sdev^2) / sum(pca$sdev^2) * 100

# Sample metadata
md <- as.data.frame(colData(vsd))
md$sample_id <- rownames(md)

# PCA scores for all PCs
pca_scores <- as.data.frame(pca$x)
pca_scores$sample_id <- rownames(pca_scores)

# Merge scores + metadata
pca_df <- pca_scores |>
  left_join(md, by = "sample_id")

# Decide aesthetics (best-effort)
color_col <- intersect(c("condition", "Condition", "group", "Group"), colnames(pca_df))
shape_col <- intersect(c("batch", "Batch"), colnames(pca_df))

color_col <- if (length(color_col) > 0) color_col[1] else NA
shape_col <- if (length(shape_col) > 0) shape_col[1] else NA

# ----------------------------
# 3. PCA plot (PC1 vs PC2)
# ----------------------------

message("Creating PCA plot (PC1 vs PC2)...")

p_pc12 <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.9)

if (!is.na(color_col)) {
  p_pc12 <- p_pc12 + aes_string(color = color_col)
}
if (!is.na(shape_col)) {
  p_pc12 <- p_pc12 + aes_string(shape = shape_col)
}

p_pc12 <- p_pc12 +
  xlab(sprintf("PC1 (%.1f%% variance)", percent_var[1])) +
  ylab(sprintf("PC2 (%.1f%% variance)", percent_var[2])) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text = element_text(colour = "black")
  )

pca_pc12_file <- file.path(fig_dir, "pca_PC1_PC2.png")
ggsave(pca_pc12_file, plot = p_pc12, width = 6, height = 5, dpi = 300)
message("Saved PCA PC1–PC2 plot to: ", pca_pc12_file)

# ----------------------------
# 4. PCA plot (PC3 vs PC4)
# ----------------------------

if (ncol(pca$x) >= 4) {
  message("Creating PCA plot (PC3 vs PC4)...")

  p_pc34 <- ggplot(pca_df, aes(x = PC3, y = PC4)) +
    geom_point(size = 3, alpha = 0.9)

  if (!is.na(color_col)) {
    p_pc34 <- p_pc34 + aes_string(color = color_col)
  }
  if (!is.na(shape_col)) {
    p_pc34 <- p_pc34 + aes_string(shape = shape_col)
  }

  p_pc34 <- p_pc34 +
    xlab(sprintf("PC3 (%.1f%% variance)", percent_var[3])) +
    ylab(sprintf("PC4 (%.1f%% variance)", percent_var[4])) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text = element_text(colour = "black")
    )

  pca_pc34_file <- file.path(fig_dir, "pca_PC3_PC4.png")
  ggsave(pca_pc34_file, plot = p_pc34, width = 6, height = 5, dpi = 300)
  message("Saved PCA PC3–PC4 plot to: ", pca_pc34_file)
} else {
  message("Not enough PCs for PC3/PC4 plot (need at least 4 PCs).")
}

# ----------------------------
# 5. 3D PCA plot (PC1, PC2, PC3)
# ----------------------------

if (ncol(pca$x) >= 3) {
  message("Creating 3D PCA plot (PC1, PC2, PC3)...")

  # Choose color for 3D plot (if available)
  color_vec <- if (!is.na(color_col)) pca_df[[color_col]] else NULL

  plt3d <- plot_ly(
    data = pca_df,
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 6,
      line = list(color = "black", width = 1)
    )
  )

  if (!is.null(color_vec)) {
    plt3d <- plt3d %>%
      add_markers(color = color_vec)
  }

  plt3d <- plt3d %>%
    layout(
      scene = list(
        xaxis = list(title = sprintf("PC1 (%.1f%%)", percent_var[1])),
        yaxis = list(title = sprintf("PC2 (%.1f%%)", percent_var[2])),
        zaxis = list(title = sprintf("PC3 (%.1f%%)", percent_var[3]))
      ),
      title = "3D PCA (PC1, PC2, PC3)"
    )

  pca_3d_file <- file.path(fig_dir, "pca_PC1_PC2_PC3_3D.html")
  saveWidget(plt3d, file = pca_3d_file, selfcontained = TRUE)
  message("Saved 3D PCA HTML plot to: ", pca_3d_file)
} else {
  message("Not enough PCs for 3D plot (need at least 3 PCs).")
}

# ----------------------------
# 6. Sample correlation heatmap
# ----------------------------

message("Creating sample correlation heatmap...")

sample_cor <- cor(mat_top)

# Use at most one annotation column for clarity
ann_cols <- intersect(c("condition", "Condition", "group", "Group"), colnames(md))
ann_df <- if (length(ann_cols) > 0) {
  md[, ann_cols[1, drop = FALSE]]
} else {
  NULL
}

heatmap_file <- file.path(fig_dir, "sample_correlation_heatmap.png")

pheatmap(
  sample_cor,
  annotation_col = ann_df,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 10,
  filename = heatmap_file
)

message("Saved correlation heatmap to: ", heatmap_file)

# ----------------------------
# 7. Top genes contributing to PC1 & PC2
# ----------------------------

message("Extracting top genes contributing to PC1 and PC2...")

loadings <- pca$rotation

get_top_genes <- function(pc_name, top_n = 50) {
  pc_load <- loadings[, pc_name]
  idx <- order(abs(pc_load), decreasing = TRUE)[seq_len(min(top_n, length(pc_load)))]
  data.frame(
    gene = rownames(loadings)[idx],
    loading = pc_load[idx],
    PC = pc_name,
    stringsAsFactors = FALSE
  )
}

top_pc1 <- get_top_genes("PC1", top_n = 50)
top_pc2 <- get_top_genes("PC2", top_n = 50)

top_genes_df <- bind_rows(top_pc1, top_pc2)

top_genes_file <- file.path(out_dir, "pca_top_genes_PC1_PC2.csv")
write.csv(top_genes_df, top_genes_file, row.names = FALSE, quote = FALSE)

message("Saved top PC genes to: ", top_genes_file)
message("PCA QC finished.")
