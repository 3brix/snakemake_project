# pca_visualization.R
library(Seurat)
library(ggplot2)

# Parse input/output arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]

# Read input
pbmc <- readRDS(input_file)

# Plot PCA loadings
pca_loadings <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave(filename = file.path(output_dir, "pca_loading_plot.png"), plot = pca_loadings, height = 6, width = 6, dpi = 300)

# PCA DimPlot
pca_dimplot <- DimPlot(pbmc, reduction = "pca") + NoLegend()
ggsave(filename = file.path(output_dir, "pca_dimplot.png"), plot = pca_dimplot, height = 6, width = 6, dpi = 300)

# Save single-dimension heatmap
png(filename = file.path(output_dir, "pca_heatmap1.png"), height = 1200, width = 1200, res = 300)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
dev.off()

# Save multi-dimension heatmap
png(filename = file.path(output_dir, "pca_heatmap15.png"), height = 2000, width = 1200, res = 300)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# Elbow plot
pca_elbowplot <- ElbowPlot(pbmc)
ggsave(filename = file.path(output_dir, "elbow_plot.png"), plot = pca_elbowplot, height = 6, width = 6, dpi = 300)
