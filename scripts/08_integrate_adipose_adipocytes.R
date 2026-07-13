###############################################################
# Script 08
# Integrate bone marrow adipocytes with public adipose tissue
# adipocyte reference dataset
###############################################################

library(Seurat)
library(ggplot2)

dir.create("figures/adipose_reference", recursive = TRUE, showWarnings = FALSE)
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ==========================
# Parameters
# ==========================

adipose_reference_path <- "data/public_reference/GSE225700/GSE225700_adipose_reference.rds"

nfeatures <- 3000
npcs <- 30
dims_use <- 1:20

set.seed(1234)

# Downsampling: keep all BM adipocytes, downsample public adipocytes
max_public_cells_per_source <- 150

# ==========================
# Colors
# ==========================

source_colors <- c(
  BM_Adipocyte = "#2D6A4F",
  Subcutaneous = "#BDBDBD"
)

# ==========================
# Load objects
# ==========================

healthy <- readRDS("results/objects/06_Healthy_final_annotated_with_compartment.rds")
adipose_ref <- readRDS(adipose_reference_path)

# ==========================
# Prepare BM adipocytes
# ==========================

bm_adipocytes <- subset(
  healthy,
  subset = final_annotation == "Adipocyte"
)

bm_adipocytes$adipocyte_source <- "BM_Adipocyte"
bm_adipocytes$dataset <- "Bone_marrow"

bm_adipocytes <- RenameCells(
  bm_adipocytes,
  add.cell.id = "BM"
)

DefaultAssay(bm_adipocytes) <- "RNA"

# ==========================
# Prepare public adipose adipocytes
# ==========================

adipose_ref$dataset <- "GSE225700"
DefaultAssay(adipose_ref) <- "RNA"

# Remove old reductions if present
adipose_ref@reductions <- list()

# Keep only subcutaneous adipocytes
adipose_ref_sc <- subset(
  adipose_ref,
  subset = adipocyte_source == "Subcutaneous"
)

# Downsample subcutaneous adipocytes
sc_cells <- colnames(adipose_ref_sc)

if(length(sc_cells) > max_public_cells_per_source){
  sc_cells_keep <- sample(sc_cells, max_public_cells_per_source)
} else {
  sc_cells_keep <- sc_cells
}

adipose_ref_ds <- subset(
  adipose_ref_sc,
  cells = sc_cells_keep
)

adipose_ref_ds <- RenameCells(
  adipose_ref_ds,
  add.cell.id = "GSE225700_SC"
)

# ==========================
# RPCA integration
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
  
  obj
})

anchors <- FindIntegrationAnchors(
  object.list = object_list,
  anchor.features = features,
  reduction = "rpca",
  dims = dims_use
)

adipocyte_integrated <- IntegrateData(
  anchorset = anchors,
  dims = dims_use,
  k.weight = 30
)

# ==========================
# PCA / UMAP after integration
# ==========================

DefaultAssay(adipocyte_integrated) <- "integrated"

adipocyte_integrated <- ScaleData(
  adipocyte_integrated,
  verbose = FALSE
)

adipocyte_integrated <- RunPCA(
  adipocyte_integrated,
  npcs = npcs,
  verbose = FALSE
)

adipocyte_integrated <- RunUMAP(
  adipocyte_integrated,
  dims = dims_use,
  reduction = "pca"
)

# ==========================
# Figures
# ==========================

pdf("figures/adipose_reference/08_umap_by_adipocyte_source.pdf", width = 8, height = 6)

print(
  DimPlot(
    adipocyte_integrated,
    reduction = "umap",
    group.by = "adipocyte_source",
    cols = source_colors
  ) +
    ggtitle("BM adipocytes integrated with adipose tissue adipocytes") +
    theme_classic()
)

dev.off()

# BM cells in front
umap_df <- as.data.frame(Embeddings(adipocyte_integrated, reduction = "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

umap_df$cell <- rownames(umap_df)
umap_df$dataset <- adipocyte_integrated$dataset
umap_df$adipocyte_source <- adipocyte_integrated$adipocyte_source

public_df <- subset(umap_df, dataset == "GSE225700")
bm_df <- subset(umap_df, dataset == "Bone_marrow")

p_front <- ggplot() +
  geom_point(
    data = public_df,
    aes(x = UMAP_1, y = UMAP_2, color = adipocyte_source),
    size = 0.4,
    alpha = 0.5
  ) +
  geom_point(
    data = bm_df,
    aes(x = UMAP_1, y = UMAP_2, color = adipocyte_source),
    size = 1.4,
    alpha = 0.95
  ) +
  scale_color_manual(values = source_colors) +
  theme_classic() +
  ggtitle("BM adipocytes overlaid on adipose tissue reference") +
  labs(color = "")

pdf("figures/adipose_reference/08_umap_bm_adipocytes_in_front.pdf", width = 8, height = 6)

print(p_front)

dev.off()

# ==========================
# Save object and tables
# ==========================

downsampling_summary <- data.frame(
  adipocyte_source = "Subcutaneous",
  original_n_cells = ncol(adipose_ref_sc),
  downsampled_n_cells = ncol(adipose_ref_ds)
)

write.csv(
  downsampling_summary,
  "results/tables/08_adipose_downsampling_summary.csv",
  row.names = FALSE
)


saveRDS(
  adipocyte_integrated,
  "results/objects/08_BM_adipose_adipocytes_integrated.rds"
)

cat("Finished Script 08\n")