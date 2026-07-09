###############################################################
# Script 05
# MSC/stromal subclustering and transfer labels to main object
###############################################################

library(Seurat)
library(ggplot2)

dir.create("figures/healthy_msc", recursive = TRUE, showWarnings = FALSE)
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Parameters
msc_labels <- c(
  "Tissue_stem_cells",
  "Endothelial_cells",
  "Osteoblasts",
  "Chondrocytes",
  "BM & Prog.",
  "Fibroblasts"
)

npcs <- 30
dims_use <- 1:15
resolution <- 0.3

marker_genes <- c(
  "PLIN1", "CEBPA",
  "WNT4", "FRZB",
  "RUNX2", "BGLAP",
  "PECAM1", "VWF"
)

# Load main healthy object
healthy <- readRDS("results/objects/04_Healthy_integrated_singler.rds")

# Subset MSC/stromal candidate cells based on SingleR labels
msc <- subset(
  healthy,
  subset = SingleR_label %in% msc_labels
)

# Recluster MSC/stromal cells using SCT
DefaultAssay(msc) <- "RNA"

msc <- SCTransform(
  msc,
  verbose = FALSE
)

DefaultAssay(msc) <- "SCT"

msc <- RunPCA(
  msc,
  npcs = npcs,
  verbose = FALSE
)

msc <- FindNeighbors(
  msc,
  dims = dims_use
)

msc <- FindClusters(
  msc,
  resolution = resolution
)

msc <- RunUMAP(
  msc,
  dims = dims_use
)

# Marker plots for deciding MSC identities
marker_genes <- marker_genes[marker_genes %in% rownames(msc)]

pdf("figures/healthy_msc/05_msc_umap_clusters.pdf", width = 7, height = 6)
print(
  DimPlot(
    msc,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  )
)
dev.off()

pdf("figures/healthy_msc/05_msc_marker_featureplots.pdf", width = 10, height = 8)
for(gene in marker_genes){
  print(
    FeaturePlot(
      msc,
      features = gene,
      reduction = "umap"
    ) +
      ggtitle(gene)
  )
}
dev.off()

pdf("figures/healthy_msc/05_msc_marker_dotplot.pdf", width = 10, height = 5)
print(
  DotPlot(
    msc,
    features = marker_genes,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis()
)
dev.off()

# Export cluster table for manual annotation
msc_cluster_table <- as.data.frame(table(msc$seurat_clusters))
colnames(msc_cluster_table) <- c("msc_cluster", "n_cells")

write.csv(
  msc_cluster_table,
  "results/tables/05_msc_cluster_table_for_manual_annotation.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# MANUAL ANNOTATION SECTION
# Edit this after inspecting the marker plots.
# Replace the example cluster numbers with your real clusters.
# ------------------------------------------------------------

msc$msc_annotation <- NA


msc$msc_annotation[msc$seurat_clusters %in% c("1")] <- "Preadipocyte"
msc$msc_annotation[msc$seurat_clusters %in% c("0")] <- "Adipocyte"
msc$msc_annotation[msc$seurat_clusters %in% c("3")] <- "Osteoblast"
msc$msc_annotation[msc$seurat_clusters %in% c("2")] <- "Endothelial"

# Transfer MSC labels back to main healthy object
healthy$final_annotation <- healthy$SingleR_label

healthy$final_annotation[colnames(msc)] <- msc$msc_annotation

# Save annotated MSC object and final healthy object
saveRDS(
  msc,
  "results/objects/05_Healthy_MSC_reclustered_annotated.rds"
)

saveRDS(
  healthy,
  "results/objects/05_Healthy_final_annotated.rds"
)

# Save final annotation table
final_annotation_table <- data.frame(
  barcode = colnames(healthy),
  sample = healthy$sample,
  seurat_cluster = healthy$seurat_clusters,
  SingleR_label = healthy$SingleR_label,
  final_annotation = healthy$final_annotation
)

write.csv(
  final_annotation_table,
  "results/tables/05_healthy_final_annotations.csv",
  row.names = FALSE
)

cat("Finished Script 05\n")