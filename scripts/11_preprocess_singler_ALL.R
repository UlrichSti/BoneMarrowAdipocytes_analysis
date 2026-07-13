###############################################################
# Script 11
# ALL sample preprocessing, clustering, UMAP, and SingleR
###############################################################

library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(ggplot2)

dir.create("figures/ALL", recursive = TRUE, showWarnings = FALSE)
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ==========================
# Parameters
# ==========================

input_object <- "results/objects/ALL_Patient_SN3_D33_clean.rds"

nfeatures <- 3000
npcs <- 30
dims_use <- 1:20
resolution <- 0.5

sample_name <- "ALL_Patient_SN3_D33"

# ==========================
# Load ALL clean object
# ==========================

all_obj <- readRDS(input_object)

all_obj$sample <- sample_name
all_obj$condition <- "ALL"

DefaultAssay(all_obj) <- "RNA"

if(inherits(all_obj[["RNA"]], "Assay5")){
  all_obj[["RNA"]] <- JoinLayers(all_obj[["RNA"]])
}

# ==========================
# Standard Seurat preprocessing
# ==========================

all_obj <- NormalizeData(
  all_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

all_obj <- FindVariableFeatures(
  all_obj,
  selection.method = "vst",
  nfeatures = nfeatures,
  verbose = FALSE
)

all_obj <- ScaleData(
  all_obj,
  verbose = FALSE
)

all_obj <- RunPCA(
  all_obj,
  npcs = npcs,
  verbose = FALSE
)

all_obj <- FindNeighbors(
  all_obj,
  dims = dims_use,
  verbose = FALSE
)

all_obj <- FindClusters(
  all_obj,
  resolution = resolution,
  verbose = FALSE
)

all_obj <- RunUMAP(
  all_obj,
  dims = dims_use,
  reduction = "pca",
  verbose = FALSE
)

# ==========================
# SingleR annotation
# ==========================

ref <- HumanPrimaryCellAtlasData()

sce <- as.SingleCellExperiment(
  all_obj,
  assay = "RNA"
)

singler <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

all_obj$SingleR_label <- singler$labels
all_obj$SingleR_pruned_label <- singler$pruned.labels

# ==========================
# Export annotation tables
# ==========================

singler_summary <- as.data.frame(
  table(all_obj$SingleR_label)
)

colnames(singler_summary) <- c("SingleR_label", "n_cells")

singler_summary <- singler_summary[
  order(singler_summary$n_cells, decreasing = TRUE),
]

write.csv(
  singler_summary,
  "results/tables/11_ALL_SingleR_label_summary.csv",
  row.names = FALSE
)

cluster_singler_table <- as.data.frame(
  table(
    seurat_cluster = all_obj$seurat_clusters,
    SingleR_label = all_obj$SingleR_label
  )
)

write.csv(
  cluster_singler_table,
  "results/tables/11_ALL_cluster_SingleR_composition.csv",
  row.names = FALSE
)

# ==========================
# UMAPs for inspection
# ==========================

pdf(
  "figures/ALL/11_ALL_UMAP_clusters.pdf",
  width = 7,
  height = 6
)

print(
  DimPlot(
    all_obj,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE
  ) +
    theme_classic() +
    ggtitle("ALL sample clusters")
)

dev.off()

pdf(
  "figures/ALL/11_ALL_UMAP_SingleR.pdf",
  width = 9,
  height = 6
)

print(
  DimPlot(
    all_obj,
    reduction = "umap",
    group.by = "SingleR_label",
    label = FALSE
  ) +
    theme_classic() +
    ggtitle("ALL sample SingleR annotation")
)

dev.off()

# ==========================
# Save object
# ==========================

saveRDS(
  all_obj,
  "results/objects/11_ALL_preprocessed_SingleR.rds"
)

cat("Finished Script 11\n")
cat("ALL nuclei:", ncol(all_obj), "\n")
cat("Output object: results/objects/11_ALL_preprocessed_SingleR.rds\n")