# variable_features_visualization.R
library(Seurat)
library(ggplot2)

# Parse input/output arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_plot <- args[2]

# Read input
pbmc <- readRDS(input_file)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
variable_features <- plot1 + plot2

# Save the plot
ggsave(filename = output_plot, plot = variable_features, height = 6, width = 12, dpi = 300)
