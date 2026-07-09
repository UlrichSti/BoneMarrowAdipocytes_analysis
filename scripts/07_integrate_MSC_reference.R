###############################################################
# Script 07
# Integrate healthy BM MSC/stromal cells with public MSC atlas
# GSE253355
###############################################################

library(Seurat)
library(ggplot2)

dir.create("figures/MSC_reference", recursive = TRUE, showWarnings = FALSE)
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ==========================
# Parameters
# ==========================

reference_path <- "data/public_reference/GSE253355/GSE253355_MSC_reference.rds"
reference_label_column <- "cluster_anno_l2"

nfeatures <- 3000
npcs <- 30
dims_use <- 1:20

# ==========================
# Colors
# ==========================

bm_msc_colors <- c(
  Adipocyte = "#2D6A4F",
  Preadipocyte = "#52B788",
  Osteoblast = "#74C69D"
)

dataset_colors <- c(
  Bone_marrow = "#8DBBE8",
  GSE253355 = "#D0D0D0"
)

plin1_colors <- c(
  `BM Adipocyte` = "#2D6A4F",
  `BM Preadipocyte` = "#52B788",
  `Public Adipo-MSC` = "grey70"
)

# ==========================
# Load objects
# ==========================

healthy <- readRDS("results/objects/06_Healthy_final_annotated_with_compartment.rds")
reference <- readRDS(reference_path)

# ==========================
# Prepare bone marrow MSC/stromal cells
# Endothelial cells are excluded before integration
# ==========================

bm_msc_populations <- c(
  "Adipocyte",
  "Preadipocyte",
  "Osteoblast"
)

bm_msc <- subset(
  healthy,
  subset = final_annotation %in% bm_msc_populations
)

bm_msc$dataset <- "Bone_marrow"
bm_msc$reference_annotation <- bm_msc$final_annotation

bm_msc <- RenameCells(
  bm_msc,
  add.cell.id = "BM"
)

DefaultAssay(bm_msc) <- "RNA"

# ==========================
# Prepare public reference
# ==========================

if(!(reference_label_column %in% colnames(reference@meta.data))){
  stop(
    paste0(
      "The metadata column '", reference_label_column,
      "' was not found in the reference object. Available columns are:\n",
      paste(colnames(reference@meta.data), collapse = ", ")
    )
  )
}

reference$dataset <- "GSE253355"
reference$reference_annotation <- reference@meta.data[[reference_label_column]]

# Remove old reductions because they can cause errors during RenameCells
reference@reductions <- list()

reference <- RenameCells(
  reference,
  add.cell.id = "GSE253355"
)

DefaultAssay(reference) <- "RNA"

# ==========================
# Public reference UMAP alone
# No integration
# ==========================

reference_umap <- reference
reference_umap@reductions <- list()

reference_umap <- NormalizeData(
  reference_umap,
  verbose = FALSE
)

reference_umap <- FindVariableFeatures(
  reference_umap,
  selection.method = "vst",
  nfeatures = nfeatures,
  verbose = FALSE
)

reference_umap <- ScaleData(
  reference_umap,
  verbose = FALSE
)

reference_umap <- RunPCA(
  reference_umap,
  npcs = npcs,
  verbose = FALSE
)

reference_umap <- RunUMAP(
  reference_umap,
  dims = dims_use,
  reduction = "pca"
)

ref_labels <- sort(unique(reference_umap$reference_annotation))

ref_colors <- setNames(
  grDevices::hcl.colors(length(ref_labels), palette = "Set 3"),
  ref_labels
)

pdf("figures/MSC_reference/07_public_reference_umap_only.pdf", width = 8, height = 6)

print(
  DimPlot(
    reference_umap,
    reduction = "umap",
    group.by = "reference_annotation",
    label = TRUE,
    repel = TRUE,
    cols = ref_colors
  ) +
    ggtitle("Public MSC reference atlas GSE253355") +
    theme_classic()
)

dev.off()

# ==========================
# Normalize both datasets before integration
# ==========================

bm_msc <- NormalizeData(
  bm_msc,
  verbose = FALSE
)

bm_msc <- FindVariableFeatures(
  bm_msc,
  selection.method = "vst",
  nfeatures = nfeatures,
  verbose = FALSE
)

reference <- NormalizeData(
  reference,
  verbose = FALSE
)

reference <- FindVariableFeatures(
  reference,
  selection.method = "vst",
  nfeatures = nfeatures,
  verbose = FALSE
)

object_list <- list(
  Bone_marrow = bm_msc,
  GSE253355 = reference
)

# ==========================
# RPCA-based integration
# ==========================

features <- SelectIntegrationFeatures(
  object.list = object_list,
  nfeatures = nfeatures
)

object_list <- lapply(object_list, function(obj) {
  
  obj <- ScaleData(
    obj,
    features = features,
    verbose = FALSE
  )
  
  obj <- RunPCA(
    obj,
    features = features,
    npcs = npcs,
    verbose = FALSE
  )
  
  return(obj)
})

anchors <- FindIntegrationAnchors(
  object.list = object_list,
  anchor.features = features,
  reduction = "rpca",
  dims = dims_use
)

msc_integrated <- IntegrateData(
  anchorset = anchors,
  dims = dims_use
)

# ==========================
# PCA / UMAP after integration
# ==========================

DefaultAssay(msc_integrated) <- "integrated"

msc_integrated <- ScaleData(
  msc_integrated,
  verbose = FALSE
)

msc_integrated <- RunPCA(
  msc_integrated,
  npcs = npcs,
  verbose = FALSE
)

msc_integrated <- RunUMAP(
  msc_integrated,
  dims = dims_use,
  reduction = "pca"
)

# ==========================
# Integrated UMAP by dataset
# ==========================

pdf("figures/MSC_reference/07_integrated_umap_by_dataset.pdf", width = 7, height = 6)

print(
  DimPlot(
    msc_integrated,
    reduction = "umap",
    group.by = "dataset",
    cols = dataset_colors
  ) +
    ggtitle("Bone marrow MSCs integrated with GSE253355") +
    theme_classic()
)

dev.off()

# ==========================
# Integrated UMAP with BM cells in front
# ==========================

umap_df <- as.data.frame(Embeddings(msc_integrated, reduction = "umap"))

umap_df$cell <- rownames(umap_df)
umap_df$dataset <- msc_integrated$dataset
umap_df$reference_annotation <- msc_integrated$reference_annotation

colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

public_df <- subset(umap_df, dataset == "GSE253355")
bm_df <- subset(umap_df, dataset == "Bone_marrow")

p_integrated_front <- ggplot() +
  geom_point(
    data = public_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey80",
    size = 0.35,
    alpha = 0.6
  ) +
  geom_point(
    data = bm_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = reference_annotation
    ),
    size = 1.2,
    alpha = 0.95
  ) +
  scale_color_manual(values = bm_msc_colors) +
  theme_classic() +
  ggtitle("Bone marrow MSC populations integrated with public MSC reference") +
  labs(color = "")

pdf(
  "figures/MSC_reference/07_integrated_umap_bm_cells_in_front.pdf",
  width = 8,
  height = 6
)

print(p_integrated_front)

dev.off()

# ==========================
# Integration summary tables
# ==========================

integration_summary <- as.data.frame(
  table(
    dataset = msc_integrated$dataset,
    annotation = msc_integrated$reference_annotation
  )
)

write.csv(
  integration_summary,
  "results/tables/07_MSC_reference_integration_summary.csv",
  row.names = FALSE
)

adipo_msc_summary <- subset(
  integration_summary,
  annotation == "Adipo-MSC"
)

write.csv(
  adipo_msc_summary,
  "results/tables/07_public_Adipo_MSC_summary.csv",
  row.names = FALSE
)

# ==========================
# PLIN1 expression in PLIN1-positive adipogenic populations
# ==========================

DefaultAssay(msc_integrated) <- "RNA"

msc_integrated[["RNA"]] <- JoinLayers(msc_integrated[["RNA"]])

msc_integrated <- NormalizeData(
  msc_integrated,
  assay = "RNA",
  verbose = FALSE
)

plin1_df <- FetchData(
  msc_integrated,
  vars = c("PLIN1", "dataset", "reference_annotation")
)

plin1_df$population <- NA

plin1_df$population[
  plin1_df$dataset == "Bone_marrow" &
    plin1_df$reference_annotation == "Adipocyte"
] <- "BM Adipocyte"

plin1_df$population[
  plin1_df$dataset == "Bone_marrow" &
    plin1_df$reference_annotation == "Preadipocyte"
] <- "BM Preadipocyte"

plin1_df$population[
  plin1_df$dataset == "GSE253355" &
    plin1_df$reference_annotation == "Adipo-MSC"
] <- "Public Adipo-MSC"

plin1_df <- subset(
  plin1_df,
  !is.na(population) & PLIN1 > 0
)

plin1_df$population <- factor(
  plin1_df$population,
  levels = c(
    "BM Adipocyte",
    "BM Preadipocyte",
    "Public Adipo-MSC"
  )
)

write.csv(
  plin1_df,
  "results/tables/07_PLIN1_positive_adipogenic_populations.csv",
  row.names = FALSE
)

# Violin plot
p_plin1_violin <- ggplot(
  plin1_df,
  aes(
    x = population,
    y = PLIN1,
    fill = population
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    color = "black",
    linewidth = 0.3
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    color = "black",
    fill = "white",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = plin1_colors) +
  theme_classic() +
  xlab("") +
  ylab("PLIN1 expression") +
  ggtitle("PLIN1 expression in PLIN1-positive adipogenic cells") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

pdf(
  "figures/MSC_reference/07_PLIN1_positive_violin.pdf",
  width = 6,
  height = 5
)

print(p_plin1_violin)

dev.off()

# Density plot
p_plin1_density <- ggplot(
  plin1_df,
  aes(
    x = PLIN1,
    color = population,
    fill = population
  )
) +
  geom_density(
    alpha = 0.25,
    linewidth = 1
  ) +
  scale_color_manual(values = plin1_colors) +
  scale_fill_manual(values = plin1_colors) +
  theme_classic() +
  xlab("PLIN1 expression") +
  ylab("Density") +
  ggtitle("PLIN1-positive cells across adipogenic populations") +
  theme(
    legend.title = element_blank()
  )

pdf(
  "figures/MSC_reference/07_PLIN1_positive_density.pdf",
  width = 6,
  height = 5
)

print(p_plin1_density)

dev.off()

# ==========================
# Save object
# ==========================

saveRDS(
  msc_integrated,
  "results/objects/07_MSC_reference_integrated.rds"
)

cat("Finished Script 07\n")