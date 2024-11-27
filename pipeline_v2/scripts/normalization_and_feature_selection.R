# normalization_and_feature_selection.R
library(Seurat)

# Parse input/output arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read input
pbmc <- readRDS(input_file)

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Save normalized data with variable features
saveRDS(pbmc, file = output_file)


