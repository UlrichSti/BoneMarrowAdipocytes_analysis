###############################################################
# Script 09
# RNA-based DEG, custom GSEA, and selected gene boxplots:
# BM adipocytes vs subcutaneous adipocytes
###############################################################

library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(Matrix)

dir.create("figures/adipose_reference", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)

# ==========================
# Parameters
# ==========================

input_object <- "results/objects/08_BM_adipose_adipocytes_integrated.rds"

group_column <- "adipocyte_source"
group_bm <- "BM_Adipocyte"
group_sc <- "Subcutaneous"

min_pct_both_groups <- 0.05

sig_padj_cutoff <- 0.001
sig_log2fc_cutoff <- 0.5

volcano_x_limit <- 8

condition_colors <- c(
  BM_Adipocyte = "#D73027",
  Subcutaneous = "#4575B4",
  not_significant = "grey80"
)

selected_genes <- c(
  "MGLL",
  "LIPE",
  "INSR",
  "IRS1",
  "CXCL12",
  "ACADVL"
)

# ==========================
# Custom gene sets
# ==========================

custom_sets <- list(
  
  LIPOLYSIS = c(
    "LIPE", "PNPLA2", "MGLL", "ABHD5", "PLIN1", "PLIN2",
    "PLIN3", "CIDEC", "CIDEA", "G0S2"
  ),
  
  FATTY_ACID_UPTAKE = c(
    "CD36", "LPL", "FABP4", "FABP5", "SLC27A1", "SLC27A4",
    "ACSL1", "ACSL3", "ACSL4", "ACSL5"
  ),
  
  FATTY_ACID_OXIDATION = c(
    "CPT1A", "CPT1B", "CPT2", "ACADM", "ACADL", "ACADVL",
    "HADHA", "HADHB", "ECHS1", "ETFA", "ETFB"
  ),
  
  INSULIN_RESPONSE = c(
    "INSR", "IRS1", "IRS2", "PIK3R1", "AKT2", "SLC2A4",
    "FOXO1", "RPTOR", "TSC2"
  ),
  
  BM_ADIPOCYTE = c(
    "CXCL12", "LEPR", "KITLG", "VCAM1", "ANGPT1", "SPP1",
    "COL1A1", "COL1A2", "DCN", "LUM"
  )
)

custom_term2gene <- do.call(
  rbind,
  lapply(names(custom_sets), function(term) {
    data.frame(
      term = term,
      gene = custom_sets[[term]]
    )
  })
)

write.csv(
  custom_term2gene,
  "results/tables/09_custom_gene_sets.csv",
  row.names = FALSE
)

# ==========================
# Load object and prepare RNA assay
# ==========================

adipo <- readRDS(input_object)

adipo <- subset(
  adipo,
  subset = adipocyte_source %in% c(group_bm, group_sc)
)

DefaultAssay(adipo) <- "RNA"

if(inherits(adipo[["RNA"]], "Assay5")){
  adipo[["RNA"]] <- JoinLayers(adipo[["RNA"]])
}

adipo <- NormalizeData(
  adipo,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

Idents(adipo) <- group_column

# ==========================
# Gene filter
# ==========================

counts <- GetAssayData(
  adipo,
  assay = "RNA",
  layer = "counts"
)

group_values <- adipo@meta.data[[group_column]]

cells_bm <- colnames(adipo)[group_values == group_bm]
cells_sc <- colnames(adipo)[group_values == group_sc]

pct_bm <- Matrix::rowMeans(counts[, cells_bm] > 0)
pct_sc <- Matrix::rowMeans(counts[, cells_sc] > 0)

genes_keep <- rownames(counts)[
  pct_bm >= min_pct_both_groups &
    pct_sc >= min_pct_both_groups
]

gene_filter_table <- data.frame(
  gene = rownames(counts),
  pct_BM = as.numeric(pct_bm),
  pct_SC = as.numeric(pct_sc),
  keep_for_DEG_GSEA = rownames(counts) %in% genes_keep
)

write.csv(
  gene_filter_table,
  "results/tables/09_gene_filter_BM_vs_SC_adipocytes.csv",
  row.names = FALSE
)

# ==========================
# DEG using RNA assay
# ==========================

deg <- FindMarkers(
  adipo,
  ident.1 = group_bm,
  ident.2 = group_sc,
  features = genes_keep,
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox",
  slot = "data"
)

deg$gene <- rownames(deg)

deg$significant <- ifelse(
  deg$p_val_adj < sig_padj_cutoff &
    abs(deg$avg_log2FC) > sig_log2fc_cutoff,
  "significant",
  "not_significant"
)

deg$enriched_in <- ifelse(
  deg$significant == "significant" &
    deg$avg_log2FC > sig_log2fc_cutoff,
  "BM_Adipocyte",
  ifelse(
    deg$significant == "significant" &
      deg$avg_log2FC < -sig_log2fc_cutoff,
    "Subcutaneous",
    "not_significant"
  )
)

write.csv(
  deg,
  "results/tables/09_DEG_RNA_BM_vs_SC_adipocytes.csv",
  row.names = FALSE
)

deg_sig <- deg[deg$significant == "significant", ]

write.csv(
  deg_sig,
  "results/tables/09_DEG_RNA_BM_vs_SC_adipocytes_significant_only.csv",
  row.names = FALSE
)

# ==========================
# Volcano plot with x-axis limits
# ==========================

deg$neg_log10_padj <- -log10(deg$p_val_adj)

max_finite <- max(
  deg$neg_log10_padj[is.finite(deg$neg_log10_padj)],
  na.rm = TRUE
)

deg$neg_log10_padj[is.infinite(deg$neg_log10_padj)] <- max_finite + 1

deg$enriched_in <- factor(
  deg$enriched_in,
  levels = c("not_significant", "Subcutaneous", "BM_Adipocyte")
)

deg_plot <- deg[order(deg$enriched_in), ]

volcano_colors <- c(
  not_significant = "grey80",
  Subcutaneous = "#4575B4",
  BM_Adipocyte = "#D73027"
)

p_volcano <- ggplot(
  deg_plot,
  aes(
    x = avg_log2FC,
    y = neg_log10_padj,
    color = enriched_in
  )
) +
  geom_point(size = 1.2, alpha = 0.8) +
  scale_color_manual(values = volcano_colors) +
  geom_vline(
    xintercept = c(-sig_log2fc_cutoff, sig_log2fc_cutoff),
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(sig_padj_cutoff),
    linetype = "dashed",
    color = "grey40"
  ) +
  coord_cartesian(
    xlim = c(-volcano_x_limit, volcano_x_limit)
  ) +
  theme_classic() +
  xlab("Average log2 fold-change") +
  ylab("-log10 adjusted P-value") +
  ggtitle("BM adipocytes vs subcutaneous adipocytes") +
  labs(color = "")

pdf(
  "figures/adipose_reference/09_volcano_RNA_BM_vs_SC_adipocytes_capped.pdf",
  width = 7,
  height = 6
)

print(p_volcano)

dev.off()

# ==========================
# GSEA using full ranked DEG table
# ==========================

deg_ranked <- deg
deg_ranked <- deg_ranked[!is.na(deg_ranked$avg_log2FC), ]
deg_ranked <- deg_ranked[!duplicated(deg_ranked$gene), ]

gene_list <- deg_ranked$avg_log2FC
names(gene_list) <- deg_ranked$gene
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_custom <- GSEA(
  geneList = gene_list,
  TERM2GENE = custom_term2gene,
  pAdjustMethod = "BH",
  minGSSize = 3,
  maxGSSize = 100,
  pvalueCutoff = 1,
  verbose = FALSE
)

gsea_table <- as.data.frame(gsea_custom)

gsea_table$direction <- ifelse(
  gsea_table$NES > 0,
  "BM adipocyte enriched",
  "Subcutaneous enriched"
)

gsea_table <- gsea_table[order(gsea_table$NES), ]

gsea_table$Description <- factor(
  gsea_table$Description,
  levels = gsea_table$Description
)

write.csv(
  gsea_table,
  "results/tables/09_GSEA_custom_gene_sets_BM_vs_SC_adipocytes.csv",
  row.names = FALSE
)

gsea_colors <- c(
  `BM adipocyte enriched` = "#D73027",
  `Subcutaneous enriched` = "#4575B4"
)

p_gsea <- ggplot(
  gsea_table,
  aes(
    x = NES,
    y = Description,
    fill = direction
  )
) +
  geom_col(width = 0.7) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  scale_fill_manual(values = gsea_colors) +
  theme_classic() +
  xlab("Normalized enrichment score") +
  ylab("") +
  ggtitle("Custom adipocyte gene set GSEA") +
  labs(fill = "")

pdf(
  "figures/adipose_reference/09_GSEA_custom_gene_sets_NES_BM_vs_SC_adipocytes.pdf",
  width = 7,
  height = 4.5
)

print(p_gsea)

dev.off()

# ==========================
# Selected gene expression summary
# ==========================

expr_selected <- FetchData(
  adipo,
  vars = c(selected_genes, group_column)
)

expr_summary <- do.call(
  rbind,
  lapply(selected_genes, function(gene) {
    
    data.frame(
      gene = gene,
      group = c(group_bm, group_sc),
      pct_expressing = c(
        mean(counts[gene, cells_bm] > 0),
        mean(counts[gene, cells_sc] > 0)
      ),
      mean_expression = c(
        mean(expr_selected[expr_selected[[group_column]] == group_bm, gene]),
        mean(expr_selected[expr_selected[[group_column]] == group_sc, gene])
      )
    )
  })
)

expr_summary <- merge(
  expr_summary,
  selected_gene_groups,
  by = "gene",
  all.x = TRUE
)

expr_summary$gene <- factor(
  expr_summary$gene,
  levels = selected_genes
)

expr_summary$group <- factor(
  expr_summary$group,
  levels = c(group_bm, group_sc)
)

write.csv(
  expr_summary,
  "results/tables/09_selected_genes_expression_summary.csv",
  row.names = FALSE
)

# ==========================
# Selected gene boxplots
# ==========================

expr_box <- FetchData(
  adipo,
  vars = c(selected_genes, group_column)
)

expr_box[[group_column]] <- factor(
  expr_box[[group_column]],
  levels = c(group_bm, group_sc)
)

expr_box_long <- do.call(
  rbind,
  lapply(selected_genes, function(gene) {
    data.frame(
      cell = rownames(expr_box),
      gene = gene,
      expression = expr_box[[gene]],
      adipocyte_source = expr_box[[group_column]]
    )
  })
)

expr_box_long$gene <- factor(
  expr_box_long$gene,
  levels = selected_genes
)

expr_box_long$adipocyte_source <- factor(
  expr_box_long$adipocyte_source,
  levels = c(group_bm, group_sc)
)

write.csv(
  expr_box_long,
  "results/tables/09_selected_genes_boxplot_expression_values.csv",
  row.names = FALSE
)

boxplot_colors <- c(
  BM_Adipocyte = "#D73027",
  Subcutaneous = "#4575B4"
)

p_box <- ggplot(
  expr_box_long,
  aes(
    x = adipocyte_source,
    y = expression,
    fill = adipocyte_source
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65,
    alpha = 0.85
  ) +
  geom_jitter(
    aes(color = adipocyte_source),
    width = 0.15,
    size = 0.5,
    alpha = 0.35
  ) +
  facet_wrap(
    ~ gene,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(values = boxplot_colors) +
  scale_color_manual(values = boxplot_colors) +
  theme_classic() +
  xlab("") +
  ylab("Log-normalized RNA expression") +
  ggtitle("Selected adipocyte genes") +
  labs(fill = "", color = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

pdf(
  "figures/adipose_reference/09_selected_genes_RNA_boxplots.pdf",
  width = 9,
  height = 6
)

print(p_box)

dev.off()

# ==========================
# PLIN1-positive cells: density plot
# BM adipocytes vs subcutaneous adipocytes
# ==========================

DefaultAssay(adipo) <- "RNA"

if(inherits(adipo[["RNA"]], "Assay5")){
  adipo[["RNA"]] <- JoinLayers(adipo[["RNA"]])
}

adipo <- NormalizeData(
  adipo,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

plin1_df <- FetchData(
  adipo,
  vars = c("PLIN1", group_column)
)

plin1_df <- plin1_df[
  plin1_df[[group_column]] %in% c(group_bm, group_sc) &
    plin1_df$PLIN1 > 0,
]

plin1_df[[group_column]] <- factor(
  plin1_df[[group_column]],
  levels = c(group_bm, group_sc)
)

write.csv(
  plin1_df,
  "results/tables/09_PLIN1_positive_expression_values_BM_vs_SC.csv",
  row.names = FALSE
)

plin1_summary <- aggregate(
  PLIN1 ~ adipocyte_source,
  plin1_df,
  function(x) {
    c(
      n_PLIN1_positive = length(x),
      mean_PLIN1 = mean(x),
      median_PLIN1 = median(x)
    )
  }
)

plin1_summary <- do.call(
  data.frame,
  plin1_summary
)

write.csv(
  plin1_summary,
  "results/tables/09_PLIN1_positive_expression_summary_BM_vs_SC.csv",
  row.names = FALSE
)

plin1_density_colors <- c(
  BM_Adipocyte = "#D73027",
  Subcutaneous = "#4575B4"
)

p_plin1_density <- ggplot(
  plin1_df,
  aes(
    x = PLIN1,
    color = adipocyte_source,
    fill = adipocyte_source
  )
) +
  geom_density(
    alpha = 0.25,
    linewidth = 1
  ) +
  scale_color_manual(values = plin1_density_colors) +
  scale_fill_manual(values = plin1_density_colors) +
  theme_classic() +
  xlab("Log-normalized PLIN1 expression") +
  ylab("Density") +
  ggtitle("PLIN1-positive adipocytes") +
  labs(color = "", fill = "")

pdf(
  "figures/adipose_reference/09_PLIN1_positive_density_BM_vs_SC_adipocytes.pdf",
  width = 7,
  height = 5
)

print(p_plin1_density)

dev.off()

# ==========================
# Save object
# ==========================

saveRDS(
  adipo,
  "results/objects/09_BM_SC_adipocytes_RNA_DEG_GSEA_input.rds"
)

cat("Finished Script 09\n")
cat("Positive log2FC / NES = BM adipocyte enriched\n")
cat("Negative log2FC / NES = subcutaneous adipocyte enriched\n")