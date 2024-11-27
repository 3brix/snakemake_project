# qc_preprocessing.R
library(Seurat)

# Parse input/output arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read input
pbmc <- readRDS(input_file)

# Calculate mitochondrial gene percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Apply filtering based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Save filtered data
saveRDS(pbmc, file = output_file)



