###############################################################
# Script 13
# Final ALL figures/tables:
# UMAP, composition tables, PLIN1-positive violin and density
###############################################################

library(Seurat)
library(ggplot2)

dir.create("figures/ALL", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)

# ==========================
# Input objects
# ==========================

all_object_path <- "results/objects/12_ALL_final_annotated.rds"
healthy_object_path <- "results/objects/06_Healthy_final_annotated_with_compartment.rds"

all_obj <- readRDS(all_object_path)
healthy <- readRDS(healthy_object_path)

all_obj$condition <- "ALL"
healthy$condition <- "Healthy"

# ==========================
# Color schemes
# ==========================

celltype_colors <- c(
  Adipocyte = "#2D6A4F",
  Preadipocyte = "#52B788",
  Osteoblast = "#74C69D",
  Endothelial = "#95D5B2",
  Fibroblast = "#40916C",
  
  T_cells = "#A9D6E5",
  B_cell = "#89C2D9",
  NK_cell = "#61A5C2",
  `Pre-B_cell_CD34-` = "#BDE0FE",
  `Pro-B_cell_CD34+` = "#A2D2FF",
  
  Monocyte = "#F4A6A6",
  Macrophage = "#E5989B",
  Neutrophils = "#F28482",
  Myelocyte = "#D65A5A",
  `Pro-Myelocyte` = "#E07A5F",
  DC = "#FFB4A2",
  
  `HSC_-G-CSF` = "#CDB4DB",
  `HSC_CD34+` = "#B392AC",
  CMP = "#D0A2F7",
  GMP = "#C77DFF",
  MEP = "#B8A1D9",
  Erythroblast = "#9D4EDD",
  
  Platelets = "#CABBE9",
  BM = "#D8F3DC"
)

condition_colors <- c(
  Healthy = "#4575B4",
  ALL = "#D73027"
)

compartment_colors <- c(
  Hematopoietic = "#BDBDBD",
  MSC_stromal = "#2D6A4F"
)

# Add fallback colors for labels not listed above
all_labels <- sort(unique(as.character(all_obj$final_annotation)))
missing_labels <- setdiff(all_labels, names(celltype_colors))

if(length(missing_labels) > 0){
  extra_colors <- setNames(
    grDevices::hcl.colors(length(missing_labels), palette = "Dark 3"),
    missing_labels
  )
  celltype_colors <- c(celltype_colors, extra_colors)
}

all_obj$final_annotation <- factor(
  as.character(all_obj$final_annotation),
  levels = all_labels
)

used_celltype_colors <- celltype_colors[levels(all_obj$final_annotation)]

# ==========================
# Final ALL UMAP with labels/colors
# ==========================

pdf(
  "figures/ALL/13_ALL_UMAP_final_annotation.pdf",
  width = 9,
  height = 6
)

print(
  DimPlot(
    all_obj,
    reduction = "umap",
    group.by = "final_annotation",
    cols = used_celltype_colors,
    label = FALSE
  ) +
    theme_classic() +
    ggtitle("ALL sample final annotation")
)

dev.off()

pdf(
  "figures/ALL/13_ALL_UMAP_compartment.pdf",
  width = 7,
  height = 6
)

print(
  DimPlot(
    all_obj,
    reduction = "umap",
    group.by = "compartment",
    cols = compartment_colors,
    label = FALSE
  ) +
    theme_classic() +
    ggtitle("ALL sample compartments")
)

dev.off()

# ==========================
# Composition tables
# ==========================

celltype_counts <- as.data.frame(
  table(final_annotation = all_obj$final_annotation)
)

colnames(celltype_counts) <- c("final_annotation", "n_cells")

celltype_counts$percent_total <- 
  celltype_counts$n_cells / sum(celltype_counts$n_cells) * 100

celltype_counts <- celltype_counts[
  order(celltype_counts$n_cells, decreasing = TRUE),
]

write.csv(
  celltype_counts,
  "results/tables/13_ALL_celltype_counts_percent.csv",
  row.names = FALSE
)

compartment_counts <- as.data.frame(
  table(compartment = all_obj$compartment)
)

colnames(compartment_counts) <- c("compartment", "n_cells")

compartment_counts$percent_total <- 
  compartment_counts$n_cells / sum(compartment_counts$n_cells) * 100

write.csv(
  compartment_counts,
  "results/tables/13_ALL_compartment_counts_percent.csv",
  row.names = FALSE
)

msc_populations <- c(
  "Adipocyte",
  "Preadipocyte",
  "Osteoblast",
  "Endothelial",
  "Fibroblast"
)

msc_counts <- subset(
  celltype_counts,
  final_annotation %in% msc_populations
)

write.csv(
  msc_counts,
  "results/tables/13_ALL_MSC_stromal_population_counts_percent.csv",
  row.names = FALSE
)

# ==========================
# Composition barplots
# ==========================

p_celltype_bar <- ggplot(
  celltype_counts,
  aes(
    x = reorder(final_annotation, n_cells),
    y = percent_total,
    fill = final_annotation
  )
) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = celltype_colors) +
  theme_classic() +
  xlab("") +
  ylab("Percent of nuclei") +
  ggtitle("ALL sample cell-type composition") +
  labs(fill = "")

pdf(
  "figures/ALL/13_ALL_celltype_composition_percent_barplot.pdf",
  width = 8,
  height = 7
)

print(p_celltype_bar)

dev.off()

p_msc_bar <- ggplot(
  msc_counts,
  aes(
    x = final_annotation,
    y = percent_total,
    fill = final_annotation
  )
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = celltype_colors) +
  theme_classic() +
  xlab("") +
  ylab("Percent of total nuclei") +
  ggtitle("ALL stromal/adipocyte populations") +
  labs(fill = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

pdf(
  "figures/ALL/13_ALL_MSC_stromal_population_percent_barplot.pdf",
  width = 6,
  height = 5
)

print(p_msc_bar)

dev.off()

# ==========================
# Prepare RNA assay helper
# ==========================

prepare_rna <- function(obj){
  
  DefaultAssay(obj) <- "RNA"
  
  if(inherits(obj[["RNA"]], "Assay5")){
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  }
  
  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  return(obj)
}

all_obj <- prepare_rna(all_obj)
healthy <- prepare_rna(healthy)

# ==========================
# PLIN1-positive adipocytes and preadipocytes
# Healthy vs ALL
# ==========================

plin_celltypes <- c(
  "Adipocyte",
  "Preadipocyte"
)

extract_plin1 <- function(obj, dataset_name){
  
  df <- FetchData(
    obj,
    vars = c("PLIN1", "final_annotation")
  )
  
  df$condition <- dataset_name
  df$cell <- rownames(df)
  
  df <- df[
    df$final_annotation %in% plin_celltypes,
  ]
  
  df$PLIN1_positive <- df$PLIN1 > 0
  
  return(df)
}

plin_healthy <- extract_plin1(
  healthy,
  "Healthy"
)

plin_all <- extract_plin1(
  all_obj,
  "ALL"
)

plin_df <- rbind(
  plin_healthy,
  plin_all
)

plin_df$condition <- factor(
  plin_df$condition,
  levels = c("Healthy", "ALL")
)

plin_df$final_annotation <- factor(
  plin_df$final_annotation,
  levels = c("Adipocyte", "Preadipocyte")
)

write.csv(
  plin_df,
  "results/tables/13_Healthy_ALL_PLIN1_expression_adipocyte_preadipocyte_all_cells.csv",
  row.names = FALSE
)

# PLIN1-positive cells only
plin_pos_df <- plin_df[
  plin_df$PLIN1_positive == TRUE,
]

write.csv(
  plin_pos_df,
  "results/tables/13_Healthy_ALL_PLIN1_positive_adipocyte_preadipocyte_cells.csv",
  row.names = FALSE
)

# Summary table
plin_summary <- do.call(
  rbind,
  lapply(split(plin_df, list(plin_df$condition, plin_df$final_annotation), drop = TRUE), function(x) {
    
    data.frame(
      condition = unique(x$condition),
      final_annotation = unique(x$final_annotation),
      n_cells_total = nrow(x),
      n_PLIN1_positive = sum(x$PLIN1_positive),
      percent_PLIN1_positive = mean(x$PLIN1_positive) * 100,
      mean_PLIN1_all_cells = mean(x$PLIN1),
      median_PLIN1_all_cells = median(x$PLIN1),
      mean_PLIN1_positive_cells = ifelse(
        sum(x$PLIN1_positive) > 0,
        mean(x$PLIN1[x$PLIN1_positive]),
        NA
      ),
      median_PLIN1_positive_cells = ifelse(
        sum(x$PLIN1_positive) > 0,
        median(x$PLIN1[x$PLIN1_positive]),
        NA
      )
    )
  })
)

write.csv(
  plin_summary,
  "results/tables/13_Healthy_ALL_PLIN1_summary_adipocyte_preadipocyte.csv",
  row.names = FALSE
)

# ==========================
# PLIN1-positive violin plot
# Manual ggplot, NOT Seurat VlnPlot
# Same y-axis for adipocytes and preadipocytes
# ==========================

plin_pos_df$condition <- factor(
  plin_pos_df$condition,
  levels = c("Healthy", "ALL")
)

plin_pos_df$final_annotation <- factor(
  plin_pos_df$final_annotation,
  levels = c("Adipocyte", "Preadipocyte")
)

plin_pos_df$plot_group <- paste(
  plin_pos_df$condition,
  plin_pos_df$final_annotation,
  sep = " "
)

plin_pos_df$plot_group <- factor(
  plin_pos_df$plot_group,
  levels = c(
    "Healthy Adipocyte",
    "ALL Adipocyte",
    "Healthy Preadipocyte",
    "ALL Preadipocyte"
  )
)

violin_colors <- c(
  `Healthy Adipocyte` = "#4575B4",
  `ALL Adipocyte` = "#D73027",
  `Healthy Preadipocyte` = "#4575B4",
  `ALL Preadipocyte` = "#D73027"
)

p_plin_violin <- ggplot(
  plin_pos_df,
  aes(
    x = plot_group,
    y = PLIN1,
    fill = plot_group
  )
) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    alpha = 0.85
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    aes(color = plot_group),
    width = 0.12,
    size = 0.5,
    alpha = 0.35
  ) +
  scale_fill_manual(values = violin_colors) +
  scale_color_manual(values = violin_colors) +
  theme_classic() +
  xlab("") +
  ylab("Log-normalized PLIN1 expression") +
  ggtitle("PLIN1-positive adipocytes and preadipocytes") +
  labs(fill = "", color = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

pdf(
  "figures/ALL/13_Healthy_ALL_PLIN1_positive_violin_same_axis.pdf",
  width = 7,
  height = 5
)

print(p_plin_violin)

dev.off()

# ==========================
# PLIN1-positive violin group counts
# ==========================

violin_group_counts <- as.data.frame(
  table(
    plot_group = plin_pos_df$plot_group
  )
)

colnames(violin_group_counts) <- c(
  "plot_group",
  "n_PLIN1_positive"
)

write.csv(
  violin_group_counts,
  "results/tables/13_Healthy_ALL_PLIN1_positive_violin_group_counts.csv",
  row.names = FALSE
)

# ==========================
# Save updated object
# ==========================

saveRDS(
  all_obj,
  "results/objects/13_ALL_final_annotated_with_figures_tables.rds"
)

cat("Finished Script 13\n")
cat("Output object: results/objects/13_ALL_final_annotated_with_figures_tables.rds\n")