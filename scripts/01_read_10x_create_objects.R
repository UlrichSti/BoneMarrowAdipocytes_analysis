###############################################################
# Script 01
# Read 10x matrices and create Seurat objects
###############################################################

library(Seurat)
library(Matrix)

# Create output folders if they do not exist
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

# Folder containing the extracted 10x matrices
base_dir <- "data/filtered_matrices"

samples <- c(
  "Healthy_1",
  "Healthy_2",
  "Healthy_3",
  "ALL_Patient_SN3_D33"
)

# Empty table for QC summary
qc_summary <- data.frame()

# Read one sample after another
for (sample in samples) {
  
  cat("Reading:", sample, "\n")
  
  counts <- Read10X(file.path(base_dir, sample))
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample,
    min.cells = 3,
    min.features = 200
  )
  
  obj$sample <- sample
  
  if(grepl("Healthy", sample)){
    obj$condition <- "Healthy"
  } else {
    obj$condition <- "ALL"
  }
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # Save Seurat object
  saveRDS(
    obj,
    file.path("results", "objects", paste0(sample, "_raw.rds"))
  )
  
  # QC summary
  qc_summary <- rbind(
    qc_summary,
    data.frame(
      sample = sample,
      nuclei = ncol(obj),
      median_nFeature_RNA = median(obj$nFeature_RNA),
      median_nCount_RNA = median(obj$nCount_RNA),
      median_percent_mt = median(obj$percent.mt)
    )
  )
  
}

write.csv(
  qc_summary,
  "results/qc/01_initial_qc_summary.csv",
  row.names = FALSE
)

cat("Finished Script 01\n")