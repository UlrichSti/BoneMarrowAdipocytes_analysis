###############################################################
# Script 10
# CellChat analysis of healthy BM snRNA-seq
# Focus: adipocyte-derived secreted signaling
###############################################################

library(Seurat)
library(CellChat)
library(Matrix)
library(ggplot2)

dir.create("figures/CellChat", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)

# ==========================
# Parameters
# ==========================

input_object <- "results/objects/06_Healthy_final_annotated_with_compartment.rds"

group_column <- "final_annotation"

source_celltype <- "Adipocyte"

min_cells_population <- 30
min_cells_interaction <- 10

cellchat_category <- "Secreted Signaling"

# ==========================
# Load healthy object
# ==========================

healthy <- readRDS(input_object)

DefaultAssay(healthy) <- "RNA"

if(inherits(healthy[["RNA"]], "Assay5")){
  healthy[["RNA"]] <- JoinLayers(healthy[["RNA"]])
}

healthy <- NormalizeData(
  healthy,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ==========================
# Filter small populations
# ==========================

population_counts_before <- as.data.frame(
  table(healthy@meta.data[[group_column]])
)

colnames(population_counts_before) <- c("celltype", "n_cells_before")

keep_celltypes <- population_counts_before$celltype[
  population_counts_before$n_cells_before >= min_cells_population
]

removed_celltypes <- population_counts_before$celltype[
  population_counts_before$n_cells_before < min_cells_population
]

healthy_cc <- subset(
  healthy,
  subset = final_annotation %in% keep_celltypes
)

healthy_cc@meta.data[[group_column]] <- droplevels(
  factor(healthy_cc@meta.data[[group_column]])
)

population_counts_after <- as.data.frame(
  table(healthy_cc@meta.data[[group_column]])
)

colnames(population_counts_after) <- c("celltype", "n_cells_after")

population_filter_summary <- merge(
  population_counts_before,
  population_counts_after,
  by = "celltype",
  all.x = TRUE
)

population_filter_summary$n_cells_after[
  is.na(population_filter_summary$n_cells_after)
] <- 0

population_filter_summary$kept_for_CellChat <- 
  population_filter_summary$n_cells_after >= min_cells_population

write.csv(
  population_filter_summary,
  "results/tables/10_CellChat_population_filter_summary.csv",
  row.names = FALSE
)

cat("Initial nuclei:", ncol(healthy), "\n")
cat("Nuclei after filtering:", ncol(healthy_cc), "\n")
cat("Removed populations:", paste(removed_celltypes, collapse = ", "), "\n")

# ==========================
# Prepare CellChat input
# ==========================

data_input <- GetAssayData(
  healthy_cc,
  assay = "RNA",
  layer = "data"
)

meta_input <- healthy_cc@meta.data

meta_input$cell <- rownames(meta_input)

meta_input <- meta_input[
  colnames(data_input),
]

stopifnot(
  all(rownames(meta_input) == colnames(data_input))
)

# ==========================
# Helper function to run CellChat
# ==========================

run_cellchat_analysis <- function(
    data_input,
    meta_input,
    group_column,
    population_size,
    suffix
){
  
  cellchat <- createCellChat(
    object = data_input,
    meta = meta_input,
    group.by = group_column
  )
  
  CellChatDB <- CellChatDB.human
  
  CellChatDB.use <- subsetDB(
    CellChatDB,
    search = cellchat_category
  )
  
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean",
    population.size = population_size
  )
  
  cellchat <- filterCommunication(
    cellchat,
    min.cells = min_cells_interaction
  )
  
  cellchat <- computeCommunProbPathway(cellchat)
  
  cellchat <- aggregateNet(cellchat)
  
  cellchat <- netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )
  
  saveRDS(
    cellchat,
    paste0("results/objects/10_CellChat_", suffix, ".rds")
  )
  
  return(cellchat)
}

# ==========================
# Run two CellChat analyses
# ==========================

# Primary analysis:
# abundance-independent signaling capacity
cellchat_intrinsic <- run_cellchat_analysis(
  data_input = data_input,
  meta_input = meta_input,
  group_column = group_column,
  population_size = FALSE,
  suffix = "intrinsic_population_size_FALSE"
)

# Secondary analysis:
# tissue-level abundance-weighted communication
cellchat_abundance <- run_cellchat_analysis(
  data_input = data_input,
  meta_input = meta_input,
  group_column = group_column,
  population_size = TRUE,
  suffix = "abundance_weighted_population_size_TRUE"
)

# ==========================
# Incoming and outgoing interaction strength for all cell types
# Intrinsic analysis
# ==========================

# cellchat_intrinsic@net$weight:
# rows = source / outgoing
# columns = target / incoming

weight_intrinsic <- cellchat_intrinsic@net$weight

interaction_strength_intrinsic <- data.frame(
  celltype = rownames(weight_intrinsic),
  outgoing_strength = rowSums(weight_intrinsic),
  incoming_strength = colSums(weight_intrinsic)
)

interaction_strength_intrinsic$total_strength <- 
  interaction_strength_intrinsic$outgoing_strength +
  interaction_strength_intrinsic$incoming_strength

interaction_strength_intrinsic <- interaction_strength_intrinsic[
  order(interaction_strength_intrinsic$total_strength, decreasing = TRUE),
]

write.csv(
  interaction_strength_intrinsic,
  "results/tables/10_CellChat_interaction_strength_all_celltypes_intrinsic.csv",
  row.names = FALSE
)

# Long format for plotting
interaction_strength_long_intrinsic <- rbind(
  data.frame(
    celltype = interaction_strength_intrinsic$celltype,
    direction = "Outgoing",
    strength = interaction_strength_intrinsic$outgoing_strength
  ),
  data.frame(
    celltype = interaction_strength_intrinsic$celltype,
    direction = "Incoming",
    strength = interaction_strength_intrinsic$incoming_strength
  )
)

interaction_strength_long_intrinsic$celltype <- factor(
  interaction_strength_long_intrinsic$celltype,
  levels = rev(interaction_strength_intrinsic$celltype)
)

interaction_strength_long_intrinsic$direction <- factor(
  interaction_strength_long_intrinsic$direction,
  levels = c("Outgoing", "Incoming")
)

write.csv(
  interaction_strength_long_intrinsic,
  "results/tables/10_CellChat_interaction_strength_all_celltypes_intrinsic_long.csv",
  row.names = FALSE
)

strength_colors <- c(
  Outgoing = "#D73027",
  Incoming = "#4575B4"
)

p_strength_intrinsic <- ggplot(
  interaction_strength_long_intrinsic,
  aes(
    x = strength,
    y = celltype,
    fill = direction
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  scale_fill_manual(values = strength_colors) +
  theme_classic() +
  xlab("Interaction strength") +
  ylab("") +
  ggtitle("CellChat interaction strength by cell type") +
  labs(fill = "")

pdf(
  "figures/CellChat/10_CellChat_interaction_strength_all_celltypes_intrinsic.pdf",
  width = 8,
  height = 7
)

print(p_strength_intrinsic)

dev.off()

# ==========================
# Incoming and outgoing interaction strength for all cell types
# Abundance-weighted analysis
# ==========================

weight_abundance <- cellchat_abundance@net$weight

interaction_strength_abundance <- data.frame(
  celltype = rownames(weight_abundance),
  outgoing_strength = rowSums(weight_abundance),
  incoming_strength = colSums(weight_abundance)
)

interaction_strength_abundance$total_strength <- 
  interaction_strength_abundance$outgoing_strength +
  interaction_strength_abundance$incoming_strength

interaction_strength_abundance <- interaction_strength_abundance[
  order(interaction_strength_abundance$total_strength, decreasing = TRUE),
]

write.csv(
  interaction_strength_abundance,
  "results/tables/10_CellChat_interaction_strength_all_celltypes_abundance_weighted.csv",
  row.names = FALSE
)

interaction_strength_long_abundance <- rbind(
  data.frame(
    celltype = interaction_strength_abundance$celltype,
    direction = "Outgoing",
    strength = interaction_strength_abundance$outgoing_strength
  ),
  data.frame(
    celltype = interaction_strength_abundance$celltype,
    direction = "Incoming",
    strength = interaction_strength_abundance$incoming_strength
  )
)

interaction_strength_long_abundance$celltype <- factor(
  interaction_strength_long_abundance$celltype,
  levels = rev(interaction_strength_abundance$celltype)
)

interaction_strength_long_abundance$direction <- factor(
  interaction_strength_long_abundance$direction,
  levels = c("Outgoing", "Incoming")
)

write.csv(
  interaction_strength_long_abundance,
  "results/tables/10_CellChat_interaction_strength_all_celltypes_abundance_weighted_long.csv",
  row.names = FALSE
)

p_strength_abundance <- ggplot(
  interaction_strength_long_abundance,
  aes(
    x = strength,
    y = celltype,
    fill = direction
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  scale_fill_manual(values = strength_colors) +
  theme_classic() +
  xlab("Interaction strength") +
  ylab("") +
  ggtitle("CellChat interaction strength by cell type, abundance-weighted") +
  labs(fill = "")

pdf(
  "figures/CellChat/10_CellChat_interaction_strength_all_celltypes_abundance_weighted.pdf",
  width = 8,
  height = 7
)

print(p_strength_abundance)

dev.off()

# ==========================
# Export global communication tables
# ==========================

comm_intrinsic <- subsetCommunication(
  cellchat_intrinsic
)

write.csv(
  comm_intrinsic,
  "results/tables/10_CellChat_all_interactions_intrinsic.csv",
  row.names = FALSE
)

comm_abundance <- subsetCommunication(
  cellchat_abundance
)

write.csv(
  comm_abundance,
  "results/tables/10_CellChat_all_interactions_abundance_weighted.csv",
  row.names = FALSE
)

# ==========================
# Extract adipocyte-derived signaling
# Primary: intrinsic analysis
# ==========================

source_celltype <- "Adipocyte"

adipo_comm <- subsetCommunication(
  cellchat_intrinsic,
  sources.use = source_celltype
)

adipo_comm <- adipo_comm[
  order(adipo_comm$prob, decreasing = TRUE),
]

write.csv(
  adipo_comm,
  "results/tables/10_CellChat_adipocyte_outgoing_interactions_intrinsic.csv",
  row.names = FALSE
)

cat("Adipocyte outgoing intrinsic interactions:", nrow(adipo_comm), "\n")

# ==========================
# Recipient summary: intrinsic analysis
# ==========================

if(nrow(adipo_comm) > 0){
  
  recipient_n <- aggregate(
    prob ~ target,
    adipo_comm,
    length
  )
  
  recipient_sum <- aggregate(
    prob ~ target,
    adipo_comm,
    sum
  )
  
  recipient_mean <- aggregate(
    prob ~ target,
    adipo_comm,
    mean
  )
  
  recipient_summary <- data.frame(
    target = recipient_n$target,
    n_interactions = recipient_n$prob,
    cumulative_probability = recipient_sum$prob,
    mean_probability = recipient_mean$prob
  )
  
  recipient_summary <- recipient_summary[
    order(recipient_summary$cumulative_probability, decreasing = TRUE),
  ]
  
} else {
  
  recipient_summary <- data.frame(
    target = character(),
    n_interactions = numeric(),
    cumulative_probability = numeric(),
    mean_probability = numeric()
  )
}

write.csv(
  recipient_summary,
  "results/tables/10_CellChat_adipocyte_recipient_summary_intrinsic.csv",
  row.names = FALSE
)

# ==========================
# Pathway summary: intrinsic analysis
# ==========================

if(nrow(adipo_comm) > 0){
  
  pathway_n <- aggregate(
    prob ~ pathway_name,
    adipo_comm,
    length
  )
  
  pathway_sum <- aggregate(
    prob ~ pathway_name,
    adipo_comm,
    sum
  )
  
  pathway_mean <- aggregate(
    prob ~ pathway_name,
    adipo_comm,
    mean
  )
  
  pathway_summary <- data.frame(
    pathway_name = pathway_n$pathway_name,
    n_interactions = pathway_n$prob,
    cumulative_probability = pathway_sum$prob,
    mean_probability = pathway_mean$prob
  )
  
  pathway_summary <- pathway_summary[
    order(pathway_summary$cumulative_probability, decreasing = TRUE),
  ]
  
} else {
  
  pathway_summary <- data.frame(
    pathway_name = character(),
    n_interactions = numeric(),
    cumulative_probability = numeric(),
    mean_probability = numeric()
  )
}

write.csv(
  pathway_summary,
  "results/tables/10_CellChat_adipocyte_pathway_summary_intrinsic.csv",
  row.names = FALSE
)

# ==========================
# Same adipocyte exports for abundance-weighted analysis
# ==========================

adipo_comm_abundance <- subsetCommunication(
  cellchat_abundance,
  sources.use = source_celltype
)

adipo_comm_abundance <- adipo_comm_abundance[
  order(adipo_comm_abundance$prob, decreasing = TRUE),
]

write.csv(
  adipo_comm_abundance,
  "results/tables/10_CellChat_adipocyte_outgoing_interactions_abundance_weighted.csv",
  row.names = FALSE
)

cat("Adipocyte outgoing abundance-weighted interactions:", nrow(adipo_comm_abundance), "\n")

# ==========================
# Circle plots: adipocyte outgoing signaling
# Intrinsic analysis
# ==========================

group_size <- as.numeric(table(cellchat_intrinsic@idents))
names(group_size) <- names(table(cellchat_intrinsic@idents))

# Number of interactions from adipocytes
if(source_celltype %in% rownames(cellchat_intrinsic@net$count)){
  
  mat_count <- matrix(
    0,
    nrow = nrow(cellchat_intrinsic@net$count),
    ncol = ncol(cellchat_intrinsic@net$count),
    dimnames = dimnames(cellchat_intrinsic@net$count)
  )
  
  mat_count[source_celltype, ] <- cellchat_intrinsic@net$count[source_celltype, ]
  
  pdf(
    "figures/CellChat/10_CellChat_adipocyte_outgoing_circle_number_intrinsic.pdf",
    width = 8,
    height = 8
  )
  
  netVisual_circle(
    mat_count,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Adipocyte outgoing interactions"
  )
  
  dev.off()
}

# Communication strength from adipocytes
if(source_celltype %in% rownames(cellchat_intrinsic@net$weight)){
  
  mat_weight <- matrix(
    0,
    nrow = nrow(cellchat_intrinsic@net$weight),
    ncol = ncol(cellchat_intrinsic@net$weight),
    dimnames = dimnames(cellchat_intrinsic@net$weight)
  )
  
  mat_weight[source_celltype, ] <- cellchat_intrinsic@net$weight[source_celltype, ]
  
  pdf(
    "figures/CellChat/10_CellChat_adipocyte_outgoing_circle_strength_intrinsic.pdf",
    width = 8,
    height = 8
  )
  
  netVisual_circle(
    mat_weight,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Adipocyte outgoing signaling strength"
  )
  
  dev.off()
}

# ==========================
# Bubble plot: adipocyte ligand-receptor interactions
# Intrinsic analysis
# ==========================

target_celltypes <- levels(cellchat_intrinsic@idents)

target_celltypes <- target_celltypes[
  target_celltypes != source_celltype
]

if(nrow(adipo_comm) > 0){
  
  pdf(
    "figures/CellChat/10_CellChat_adipocyte_outgoing_bubble_intrinsic.pdf",
    width = 12,
    height = 8
  )
  
  print(
    netVisual_bubble(
      cellchat_intrinsic,
      sources.use = source_celltype,
      targets.use = target_celltypes,
      remove.isolate = FALSE
    )
  )
  
  dev.off()
}

# ==========================
# Top adipocyte-derived interactions plot
# Intrinsic analysis
# ==========================

if(nrow(adipo_comm) > 0){
  
  top_n <- 25
  
  top_adipo_comm <- head(adipo_comm, top_n)
  
  top_adipo_comm$interaction_label <- paste0(
    top_adipo_comm$ligand,
    " - ",
    top_adipo_comm$receptor,
    " → ",
    top_adipo_comm$target
  )
  
  top_adipo_comm$interaction_label <- factor(
    top_adipo_comm$interaction_label,
    levels = rev(top_adipo_comm$interaction_label)
  )
  
  p_top <- ggplot(
    top_adipo_comm,
    aes(
      x = prob,
      y = interaction_label
    )
  ) +
    geom_col(fill = "#2D6A4F") +
    theme_classic() +
    xlab("Communication probability") +
    ylab("") +
    ggtitle("Top adipocyte-derived signaling interactions")
  
  pdf(
    "figures/CellChat/10_CellChat_top_adipocyte_outgoing_interactions_intrinsic.pdf",
    width = 9,
    height = 7
  )
  
  print(p_top)
  
  dev.off()
}

# ==========================
# Save filtered Seurat input
# ==========================

saveRDS(
  healthy_cc,
  "results/objects/10_Healthy_CellChat_input_filtered.rds"
)

cat("Finished Script 10 export and adipocyte visualization section\n")