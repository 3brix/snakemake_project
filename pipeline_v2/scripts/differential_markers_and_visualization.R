# differential_markers_and_visualization.R
library(Seurat)
library(ggplot2)
library(dplyr)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
markers_output <- args[2]
plot_output <- args[3]

pbmc <- readRDS(input_file)

# Find markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

# Violin plot
marker_violin <- VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
ggsave("results/plots/marker_violin.png", marker_violin, height = 6, width = 6, dpi = 300)

# Violin plot (log scale)
marker_violin_raw <- VlnPlot(pbmc, features = c("NKG7", "PF4"), log = TRUE) 
ggsave("results/plots/marker_violin_raw.png", marker_violin_raw, height = 6, width = 6, dpi = 300)

# Feature plot
marker_feature <- FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
ggsave("results/plots/marker_feature.png", marker_feature, height = 8, width = 10, dpi = 300)

# Heatmap of top 10 markers for each cluster
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
marker_heatmap <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
ggsave("results/plots/marker_heatmap.png", marker_heatmap, height = 8, width = 10, dpi = 300)


# UMAP plot
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = plot_output, height = 7, width = 12, plot = plot, dpi = 300)


# Save marker data
saveRDS(pbmc.markers, file = markers_output)
