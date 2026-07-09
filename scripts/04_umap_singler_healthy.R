###############################################################
# Script 04
# SingleR annotation for healthy integrated data
###############################################################

library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)

healthy <- readRDS("results/objects/03_Healthy_integrated.rds")

DefaultAssay(healthy) <- "RNA"

healthy[["RNA"]] <- JoinLayers(healthy[["RNA"]])

healthy <- NormalizeData(
  healthy,
  assay = "RNA",
  verbose = FALSE
)

ref <- HumanPrimaryCellAtlasData()

sce <- as.SingleCellExperiment(
  healthy,
  assay = "RNA"
)

singler_results <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

healthy$SingleR_label <- singler_results$labels
healthy$SingleR_pruned_label <- singler_results$pruned.labels

singler_table <- data.frame(
  barcode = colnames(healthy),
  sample = healthy$sample,
  seurat_cluster = healthy$seurat_clusters,
  SingleR_label = healthy$SingleR_label,
  SingleR_pruned_label = healthy$SingleR_pruned_label
)

write.csv(
  singler_table,
  "results/tables/04_healthy_singler_annotations.csv",
  row.names = FALSE
)

singler_summary <- as.data.frame(table(healthy$SingleR_label))
colnames(singler_summary) <- c("SingleR_label", "n_cells")
singler_summary <- singler_summary[order(-singler_summary$n_cells), ]

write.csv(
  singler_summary,
  "results/tables/04_healthy_singler_label_summary.csv",
  row.names = FALSE
)

saveRDS(
  healthy,
  "results/objects/04_Healthy_integrated_singler.rds"
)

cat("Finished Script 04\n")