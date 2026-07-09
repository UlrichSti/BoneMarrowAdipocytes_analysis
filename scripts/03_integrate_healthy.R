###############################################################
# Script 03
# Integrate healthy donor samples
###############################################################

library(Seurat)

# Parameters
nfeatures <- 3000
npcs <- 30
dims_use <- 1:20
resolution <- 0.5

# Load clean healthy samples
Healthy_1 <- readRDS("results/objects/Healthy_1_clean.rds")
Healthy_2 <- readRDS("results/objects/Healthy_2_clean.rds")
Healthy_3 <- readRDS("results/objects/Healthy_3_clean.rds")

healthy_list <- list(
  Healthy_1 = Healthy_1,
  Healthy_2 = Healthy_2,
  Healthy_3 = Healthy_3
)

# Normalize and find variable genes
healthy_list <- lapply(healthy_list, function(obj) {
  
  obj <- NormalizeData(obj)
  
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = nfeatures
  )
  
  return(obj)
})

# Select integration features
features <- SelectIntegrationFeatures(
  object.list = healthy_list,
  nfeatures = nfeatures
)

# Scale and run PCA on each sample
healthy_list <- lapply(healthy_list, function(obj) {
  
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

# Find integration anchors
anchors <- FindIntegrationAnchors(
  object.list = healthy_list,
  anchor.features = features,
  reduction = "cca",
  dims = dims_use
)

# Integrate data
healthy <- IntegrateData(
  anchorset = anchors,
  dims = dims_use
)

# Run dimensional reduction and clustering
DefaultAssay(healthy) <- "integrated"

healthy <- ScaleData(healthy, verbose = FALSE)

healthy <- RunPCA(
  healthy,
  npcs = npcs,
  verbose = FALSE
)

healthy <- RunUMAP(
  healthy,
  dims = dims_use
)

healthy <- FindNeighbors(
  healthy,
  dims = dims_use
)

healthy <- FindClusters(
  healthy,
  resolution = resolution
)

# Save integrated object
saveRDS(
  healthy,
  "results/objects/03_Healthy_integrated.rds"
)

cat("Finished Script 03\n")