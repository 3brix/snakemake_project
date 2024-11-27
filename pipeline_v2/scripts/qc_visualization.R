# qc_visualization.R
library(Seurat)
library(ggplot2)

# Parse input/output arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_violin <- args[2]
output_scatter <- args[3]

# Read input
pbmc <- readRDS(input_file)

# Create QC violin plot
QC_violin <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = output_violin, plot = QC_violin, height = 6, width = 12, dpi = 300)

# Create feature scatter plots
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
feature_scatter <- plot1 + plot2
ggsave(filename = output_scatter, plot = feature_scatter, height = 6, width = 12, dpi = 300)
