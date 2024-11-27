# pca_visualization.R
library(Seurat)
library(ggplot2)

# Parse input/output arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_pca_loadings <- args[2]
output_pca_dimplot <- args[3]
output_pca_heatmap1 <- args[4]
output_pca_heatmap15 <- args[5]
output_elbow_plot <- args[6]

# Read input Seurat object
pbmc <- readRDS(input_file)

# PCA Loadings Plot
pca_loadings <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave(filename = output_pca_loadings, plot = pca_loadings, height = 6, width = 6, dpi = 300)

# PCA DimPlot
pca_dimplot <- DimPlot(pbmc, reduction = "pca") + NoLegend()
ggsave(filename = output_pca_dimplot, plot = pca_dimplot, height = 6, width = 6, dpi = 300)

# PCA Heatmap (1st Dimension)
png(filename = output_pca_heatmap1, height = 1200, width = 1200, res = 300)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
dev.off()

# PCA Heatmap (1-15 Dimensions)
png(filename = output_pca_heatmap15, height = 2000, width = 1200, res = 300)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# PCA Elbow Plot
pca_elbowplot <- ElbowPlot(pbmc)
ggsave(filename = output_elbow_plot, plot = pca_elbowplot, height = 6, width = 6, dpi = 300)
