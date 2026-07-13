###############################################################
# Script 12
# ALL sample MSC/stromal subclustering and manual annotation
###############################################################

library(Seurat)
library(ggplot2)

dir.create("figures/ALL", recursive = TRUE, showWarnings = FALSE)
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ==========================
# Parameters
# ==========================

input_object <- "results/objects/11_ALL_preprocessed_SingleR.rds"

npcs <- 30
dims_use <- 1:15
resolution <- 0.2

# Same stromal-like SingleR labels used for healthy analysis
msc_labels <- c(
  "Tissue_stem_cells",
  "Endothelial_cells",
  "Osteoblasts",
  "Chondrocytes",
  "Fibroblasts",
  "Smooth_muscle_cells"
)

marker_genes <- c(
  "PLIN1", "CEBPA",
  "WNT4", "FRZB",
  "RUNX2", "BGLAP",
  "PECAM1", "VWF",
  "ACTA2", "TAGLN"
)

# ==========================
# Load ALL object
# ==========================

all_obj <- readRDS(input_object)

DefaultAssay(all_obj) <- "RNA"

if(inherits(all_obj[["RNA"]], "Assay5")){
  all_obj[["RNA"]] <- JoinLayers(all_obj[["RNA"]])
}

# ==========================
# Subset stromal-like cells
# ==========================

all_msc <- subset(
  all_obj,
  subset = SingleR_label %in% msc_labels
)

write.csv(
  as.data.frame(table(all_msc$SingleR_label)),
  "results/tables/12_ALL_MSC_subset_SingleR_labels.csv",
  row.names = FALSE
)

cat("ALL MSC/stromal subset nuclei:", ncol(all_msc), "\n")

# ==========================
# SCT reclustering
# ==========================

all_msc <- SCTransform(
  all_msc,
  assay = "RNA",
  verbose = FALSE
)

DefaultAssay(all_msc) <- "SCT"

all_msc <- RunPCA(
  all_msc,
  npcs = npcs,
  verbose = FALSE
)

all_msc <- FindNeighbors(
  all_msc,
  dims = dims_use,
  verbose = FALSE
)

all_msc <- FindClusters(
  all_msc,
  resolution = resolution,
  verbose = FALSE
)

all_msc <- RunUMAP(
  all_msc,
  dims = dims_use,
  reduction = "pca",
  verbose = FALSE
)

# ==========================
# Marker inspection plots
# ==========================

pdf(
  "figures/ALL/12_ALL_MSC_UMAP_clusters.pdf",
  width = 7,
  height = 6
)

print(
  DimPlot(
    all_msc,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE
  ) +
    theme_classic() +
    ggtitle("ALL MSC/stromal subset clusters")
)

dev.off()

available_markers <- marker_genes[
  marker_genes %in% rownames(all_msc)
]

pdf(
  "figures/ALL/12_ALL_MSC_marker_featureplots.pdf",
  width = 10,
  height = 8
)

print(
  FeaturePlot(
    all_msc,
    features = available_markers,
    reduction = "umap",
    ncol = 4
  )
)

dev.off()

pdf(
  "figures/ALL/12_ALL_MSC_marker_dotplot.pdf",
  width = 8,
  height = 5
)

print(
  DotPlot(
    all_msc,
    features = available_markers,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    theme_classic() +
    ggtitle("ALL MSC/stromal marker expression")
)

dev.off()

# ==========================
# Export tables for manual annotation
# ==========================

cluster_table <- as.data.frame(
  table(
    seurat_cluster = all_msc$seurat_clusters,
    SingleR_label = all_msc$SingleR_label
  )
)

write.csv(
  cluster_table,
  "results/tables/12_ALL_MSC_cluster_SingleR_table_for_manual_annotation.csv",
  row.names = FALSE
)

# Average marker expression by cluster
DefaultAssay(all_msc) <- "SCT"

avg_markers <- AverageExpression(
  all_msc,
  assays = "SCT",
  features = available_markers,
  group.by = "seurat_clusters"
)$SCT

write.csv(
  as.data.frame(avg_markers),
  "results/tables/12_ALL_MSC_average_marker_expression_by_cluster.csv"
)

# ==========================
# Manual annotation section
# ==========================

all_msc$msc_annotation <- NA

all_msc$msc_annotation[all_msc$seurat_clusters %in% c("2")] <- "Preadipocyte"
all_msc$msc_annotation[all_msc$seurat_clusters %in% c("4")] <- "Adipocyte"
all_msc$msc_annotation[all_msc$seurat_clusters %in% c("0")] <- "Endothelial"
all_msc$msc_annotation[all_msc$seurat_clusters %in% c("3", "1")] <- "Fibroblast"

# Stop if any MSC clusters are still unannotated
if(any(is.na(all_msc$msc_annotation))){
  
  unannotated_clusters <- unique(
    all_msc$seurat_clusters[is.na(all_msc$msc_annotation)]
  )
  
  saveRDS(
    all_msc,
    "results/objects/12_ALL_MSC_reclustered_partly_annotated.rds"
  )
  
  stop(
    paste(
      "Some MSC clusters are still unannotated:",
      paste(unannotated_clusters, collapse = ", ")
    )
  )
}

# ==========================
# Transfer manual MSC labels back to ALL object
# ==========================

all_obj$final_annotation <- all_obj$SingleR_label

all_obj$final_annotation[colnames(all_msc)] <- all_msc$msc_annotation

# Add compartment labels
msc_populations <- c(
  "Adipocyte",
  "Preadipocyte",
  "Osteoblast",
  "Endothelial",
  "Fibroblast"
)

all_obj$compartment <- ifelse(
  all_obj$final_annotation %in% msc_populations,
  "MSC_stromal",
  "Hematopoietic"
)

# ==========================
# Save annotated objects and tables
# ==========================

write.csv(
  as.data.frame(table(all_msc$seurat_clusters, all_msc$msc_annotation)),
  "results/tables/12_ALL_MSC_cluster_manual_annotation_table.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(table(all_obj$final_annotation)),
  "results/tables/12_ALL_final_annotation_summary.csv",
  row.names = FALSE
)

saveRDS(
  all_msc,
  "results/objects/12_ALL_MSC_reclustered_annotated.rds"
)

saveRDS(
  all_obj,
  "results/objects/12_ALL_final_annotated.rds"
)

cat("Finished Script 12\n")
cat("Output object: results/objects/12_ALL_final_annotated.rds\n")