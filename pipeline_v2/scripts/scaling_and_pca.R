# scaling_and_pca.R
library(Seurat)

# Parse input/output arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read input
pbmc <- readRDS(input_file)

# Scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Save processed data
saveRDS(pbmc, file = output_file)

