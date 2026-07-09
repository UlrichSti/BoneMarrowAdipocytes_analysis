###############################################################
# Script 06
# Final healthy figures and composition tables
###############################################################

library(Seurat)
library(ggplot2)

dir.create("figures/healthy_final", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ==========================
# Load final annotated object
# ==========================

healthy <- readRDS("results/objects/05_Healthy_final_annotated.rds")

DefaultAssay(healthy) <- "RNA"

# ==========================
# Colors
# ==========================

sample_colors <- c(
  Healthy_1 = "#F6D77A",
  Healthy_2 = "#F4A6C1",
  Healthy_3 = "#8DBBE8"
)

celltype_colors <- c(
  # MSC / stromal - greens
  Adipocyte = "#2D6A4F",
  Preadipocyte = "#52B788",
  Osteoblast = "#74C69D",
  Endothelial = "#95D5B2",
  
  # Lymphoid - blues
  T_cells = "#A9D6E5",
  B_cell = "#89C2D9",
  NK_cell = "#61A5C2",
  `Pre-B_cell_CD34-` = "#BDE0FE",
  `Pro-B_cell_CD34+` = "#A2D2FF",
  
  # Myeloid - reds
  Monocyte = "#F4A6A6",
  Macrophage = "#E5989B",
  Neutrophils = "#F28482",
  Myelocyte = "#D65A5A",
  `Pro-Myelocyte` = "#E07A5F",
  DC = "#FFB4A2",
  
  # HSC / progenitors - purples
  `HSC_-G-CSF` = "#CDB4DB",
  `HSC_CD34+` = "#B392AC",
  CMP = "#D0A2F7",
  GMP = "#C77DFF",
  MEP = "#B8A1D9",
  Erythroblast = "#9D4EDD",
  
  # Other
  Platelets = "#CABBE9",
  BM = "#D8F3DC"
)

celltype_colors <- celltype_colors[
  names(celltype_colors) %in% unique(healthy$final_annotation)
]

# ==========================
# Basic cell type counts
# ==========================

celltype_counts <- as.data.frame(table(healthy$final_annotation))
colnames(celltype_counts) <- c("cell_type", "n_cells")

write.csv(
  celltype_counts,
  "results/tables/06_final_celltype_counts.csv",
  row.names = FALSE
)

sample_counts <- as.data.frame(table(healthy$sample))
colnames(sample_counts) <- c("sample", "n_cells")

write.csv(
  sample_counts,
  "results/tables/06_sample_counts.csv",
  row.names = FALSE
)

# ==========================
# Final UMAPs
# ==========================

pdf("figures/healthy_final/06_umap_final_annotation.pdf", width = 9, height = 7)

print(
  DimPlot(
    healthy,
    reduction = "umap",
    group.by = "final_annotation",
    label = FALSE,
    cols = celltype_colors
  ) +
    ggtitle("Healthy bone marrow - final annotation") +
    theme_classic()
)

dev.off()

pdf("figures/healthy_final/06_umap_by_sample.pdf", width = 8, height = 6)

print(
  DimPlot(
    healthy,
    reduction = "umap",
    group.by = "sample",
    label = FALSE,
    cols = sample_colors
  ) +
    ggtitle("Healthy bone marrow - donors") +
    theme_classic()
)

dev.off()

# ==========================
# Composition tables
# ==========================

msc_populations <- c(
  "Adipocyte",
  "Preadipocyte",
  "Osteoblast",
  "Endothelial"
)

healthy$compartment <- ifelse(
  healthy$final_annotation %in% msc_populations,
  "MSC_stromal",
  "Hematopoietic"
)

# Total hematopoietic vs MSC/stromal
compartment_total <- as.data.frame(table(healthy$compartment))
colnames(compartment_total) <- c("compartment", "n_cells")

compartment_total$percent_total <- round(
  compartment_total$n_cells / sum(compartment_total$n_cells) * 100,
  2
)

write.csv(
  compartment_total,
  "results/tables/06_compartment_total.csv",
  row.names = FALSE
)

# Hematopoietic vs MSC/stromal per sample
compartment_by_sample <- as.data.frame(
  table(
    sample = healthy$sample,
    compartment = healthy$compartment
  )
)

colnames(compartment_by_sample) <- c(
  "sample",
  "compartment",
  "n_cells"
)

sample_totals <- as.data.frame(table(healthy$sample))
colnames(sample_totals) <- c("sample", "total_cells")

compartment_by_sample <- merge(
  compartment_by_sample,
  sample_totals,
  by = "sample"
)

compartment_by_sample$percent_of_sample <- round(
  compartment_by_sample$n_cells / compartment_by_sample$total_cells * 100,
  2
)

write.csv(
  compartment_by_sample,
  "results/tables/06_compartment_by_sample.csv",
  row.names = FALSE
)

# MSC population counts total
msc_meta <- healthy@meta.data[
  healthy$final_annotation %in% msc_populations,
]

msc_total <- as.data.frame(table(msc_meta$final_annotation))
colnames(msc_total) <- c("msc_population", "n_cells")

msc_total$percent_of_all_cells <- round(
  msc_total$n_cells / ncol(healthy) * 100,
  2
)

msc_total$percent_within_msc <- round(
  msc_total$n_cells / sum(msc_total$n_cells) * 100,
  2
)

write.csv(
  msc_total,
  "results/tables/06_msc_total.csv",
  row.names = FALSE
)

# MSC population counts per sample
msc_by_sample <- as.data.frame(
  table(
    sample = msc_meta$sample,
    msc_population = msc_meta$final_annotation
  )
)

colnames(msc_by_sample) <- c(
  "sample",
  "msc_population",
  "n_cells"
)

msc_by_sample <- merge(
  msc_by_sample,
  sample_totals,
  by = "sample"
)

msc_totals_by_sample <- aggregate(
  n_cells ~ sample,
  data = msc_by_sample,
  sum
)

colnames(msc_totals_by_sample) <- c(
  "sample",
  "total_msc_cells"
)

msc_by_sample <- merge(
  msc_by_sample,
  msc_totals_by_sample,
  by = "sample"
)

msc_by_sample$percent_of_sample <- round(
  msc_by_sample$n_cells / msc_by_sample$total_cells * 100,
  2
)

msc_by_sample$percent_within_sample_msc <- round(
  msc_by_sample$n_cells / msc_by_sample$total_msc_cells * 100,
  2
)

write.csv(
  msc_by_sample,
  "results/tables/06_msc_by_sample.csv",
  row.names = FALSE
)

# ==========================
# Marker DotPlots
# ==========================

healthy$compartment <- factor(
  healthy$compartment,
  levels = c("Hematopoietic", "MSC_stromal")
)

compartment_markers <- c(
  "PTPRC", "LAPTM5", "CORO1A",
  "COL1A1", "COL1A2", "CXCL12"
)

compartment_markers <- compartment_markers[
  compartment_markers %in% rownames(healthy)
]

p_compartment <- DotPlot(
  healthy,
  features = compartment_markers,
  group.by = "compartment",
  cols = c("grey90", "black"),
  dot.scale = 7
) +
  RotatedAxis() +
  ggtitle("Hematopoietic vs MSC/stromal marker expression") +
  theme_classic()

pdf(
  "figures/healthy_final/06_dotplot_hematopoietic_vs_msc.pdf",
  width = 7,
  height = 4
)

print(p_compartment)

dev.off()

msc <- subset(
  healthy,
  subset = final_annotation %in% msc_populations
)

msc$final_annotation <- factor(
  msc$final_annotation,
  levels = c(
    "Adipocyte",
    "Preadipocyte",
    "Osteoblast",
    "Endothelial"
  )
)

msc_markers <- c(
  "PLIN1", "CEBPA",
  "WNT4", "FRZB",
  "RUNX2", "BGLAP",
  "PECAM1", "VWF"
)

msc_markers <- msc_markers[
  msc_markers %in% rownames(msc)
]

p_msc_dotplot <- DotPlot(
  msc,
  features = msc_markers,
  group.by = "final_annotation",
  cols = c("grey90", "black"),
  dot.scale = 7
) +
  RotatedAxis() +
  ggtitle("MSC/stromal population marker expression") +
  theme_classic()

pdf(
  "figures/healthy_final/06_dotplot_msc_populations.pdf",
  width = 8,
  height = 4
)

print(p_msc_dotplot)

dev.off()

# ==========================
# Bar plots
# ==========================

# All cell types per donor
all_celltypes_by_sample <- as.data.frame(
  table(
    sample = healthy$sample,
    cell_type = healthy$final_annotation
  )
)

colnames(all_celltypes_by_sample) <- c(
  "sample",
  "cell_type",
  "n_cells"
)

all_celltypes_by_sample <- merge(
  all_celltypes_by_sample,
  sample_totals,
  by = "sample"
)

all_celltypes_by_sample$percent_of_sample <- round(
  all_celltypes_by_sample$n_cells / all_celltypes_by_sample$total_cells * 100,
  2
)

write.csv(
  all_celltypes_by_sample,
  "results/tables/06_all_celltypes_by_sample.csv",
  row.names = FALSE
)

p_all <- ggplot(
  all_celltypes_by_sample,
  aes(
    x = sample,
    y = percent_of_sample,
    fill = cell_type
  )
) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = celltype_colors) +
  ylab("Percent of donor nuclei") +
  xlab("") +
  ggtitle("Cell type composition per healthy donor") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

pdf(
  "figures/healthy_final/06_barplot_all_celltypes_by_donor.pdf",
  width = 10,
  height = 6
)

print(p_all)

dev.off()

# MSC populations per donor
p_msc_bar <- ggplot(
  msc_by_sample,
  aes(
    x = sample,
    y = percent_within_sample_msc,
    fill = msc_population
  )
) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = celltype_colors[msc_populations]) +
  ylab("Percent of MSC/stromal compartment") +
  xlab("") +
  ggtitle("MSC/stromal composition per healthy donor") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

pdf(
  "figures/healthy_final/06_barplot_msc_populations_by_donor.pdf",
  width = 7,
  height = 5
)

print(p_msc_bar)

dev.off()

# Total adipocyte number per donor
adipocyte_counts <- as.data.frame(
  table(
    sample = healthy$sample[
      healthy$final_annotation == "Adipocyte"
    ]
  )
)

colnames(adipocyte_counts) <- c(
  "sample",
  "n_adipocytes"
)

write.csv(
  adipocyte_counts,
  "results/tables/06_adipocyte_counts_by_sample.csv",
  row.names = FALSE
)

p_adipo <- ggplot(
  adipocyte_counts,
  aes(
    x = sample,
    y = n_adipocytes,
    fill = sample
  )
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = sample_colors) +
  ylab("Number of adipocytes") +
  xlab("") +
  ggtitle("Adipocyte number per healthy donor") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

pdf(
  "figures/healthy_final/06_barplot_adipocyte_counts_by_donor.pdf",
  width = 5,
  height = 5
)

print(p_adipo)

dev.off()

# ==========================
# Save updated object
# ==========================

saveRDS(
  healthy,
  "results/objects/06_Healthy_final_annotated_with_compartment.rds"
)

cat("Finished Script 06\n")