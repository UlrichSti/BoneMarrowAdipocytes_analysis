###############################################################
# Script 02
# QC filtering and doublet removal
###############################################################

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)

dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

samples <- c(
  "Healthy_1",
  "Healthy_2",
  "Healthy_3",
  "ALL_Patient_SN3_D33"
)

qc_summary <- data.frame()

for(sample in samples){
  
  cat("Processing:", sample, "\n")
  
  obj <- readRDS(
    file.path("results", "objects", paste0(sample, "_raw.rds"))
  )
  
  initial_cells <- ncol(obj)
  
  # Light filtering
  obj <- subset(
    obj,
    subset =
      nFeature_RNA >= 300 &
      nCount_RNA >= 500
  )
  
  after_light_filter <- ncol(obj)
  
  # Doublet detection
  sce <- as.SingleCellExperiment(obj)
  
  sce <- scDblFinder(sce)
  
  obj$scDblFinder.class <- sce$scDblFinder.class
  obj$scDblFinder.score <- sce$scDblFinder.score
  
  predicted_doublets <- sum(obj$scDblFinder.class == "doublet")
  
  # Remove doublets + mitochondrial filter
  obj <- subset(
    obj,
    subset =
      scDblFinder.class == "singlet" &
      percent.mt < 10
  )
  
  final_cells <- ncol(obj)
  
  saveRDS(
    obj,
    file.path("results", "objects", paste0(sample, "_clean.rds"))
  )
  
  qc_summary <- rbind(
    qc_summary,
    data.frame(
      sample = sample,
      initial_nuclei = initial_cells,
      after_light_filter = after_light_filter,
      predicted_doublets = predicted_doublets,
      final_nuclei = final_cells,
      removed_total = initial_cells - final_cells,
      percent_removed = round((initial_cells-final_cells)/initial_cells*100,2)
    )
  )
  
}

write.csv(
  qc_summary,
  "results/qc/02_qc_cleaning_summary.csv",
  row.names = FALSE
)

cat("Finished Script 02\n")