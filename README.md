# BoneMarrowAdipocytes_analysis

Reproducible R analysis pipeline for single-nucleus RNA-seq analysis of human bone marrow adipocyte and stromal populations in healthy donor bone marrow and acute lymphoblastic leukemia (ALL) bone marrow.

## Project overview

This repository contains the analysis scripts used to process, annotate, and compare human bone marrow single-nucleus RNA-seq datasets generated using the 10x Genomics Chromium Next GEM Flex Gene Expression Singleplex workflow.

The analysis includes:

- preprocessing of Cell Ranger filtered 10x feature-barcode matrices
- doublet detection using `scDblFinder`
- quality control and filtering
- Seurat v5 normalization, clustering, and UMAP visualization
- SingleR-assisted cell type annotation
- manual annotation of mesenchymal/stromal populations
- focused reclustering of adipocyte-associated stromal populations
- integration with public MSC and adipose tissue reference datasets
- differential expression analysis between bone marrow and subcutaneous adipocytes
- custom adipocyte gene-set enrichment analysis
- CellChat-based analysis of adipocyte-derived secreted signaling
- comparison of healthy and ALL bone marrow adipocyte/stromal populations

## Data

Raw sequencing data, filtered 10x matrices, large Seurat objects, and generated output files are not included in this GitHub repository.

The expected local data structure is:

```text
BoneMarrowAdipocytes_analysis/
├── data/
│   ├── filtered_matrices/
│   │   ├── Healthy_1/
│   │   ├── Healthy_2/
│   │   ├── Healthy_3/
│   │   └── ALL_Patient_SN3_D33/
│   └── public_reference/
│       ├── GSE253355/
│       └── GSE225700/
├── scripts/
├── results/
└── figures/
cd "D:\2026\3. Projekte\Bone Marrow Adipocytes\BoneMarrowAdipocytes_analysis"
notepad README.md
notepad README.md
cd "D:\2026\3. Projekte\Bone Marrow Adipocytes\BoneMarrowAdipocytes_analysis"

@'
# BoneMarrowAdipocytes_analysis

Reproducible R analysis pipeline for single-nucleus RNA-seq analysis of human bone marrow adipocyte and stromal populations in healthy donor bone marrow and acute lymphoblastic leukemia bone marrow.

## Project overview

This repository contains the analysis scripts used to process, annotate, and compare human bone marrow single-nucleus RNA-seq datasets generated using the 10x Genomics Chromium Next GEM Flex Gene Expression Singleplex workflow.

The analysis includes:

- preprocessing of Cell Ranger filtered 10x feature-barcode matrices
- doublet detection using scDblFinder
- quality control and filtering
- Seurat v5 normalization, clustering, and UMAP visualization
- SingleR-assisted cell type annotation
- manual annotation of mesenchymal/stromal populations
- focused reclustering of adipocyte-associated stromal populations
- integration with public MSC and adipose tissue reference datasets
- differential expression analysis between bone marrow and subcutaneous adipocytes
- custom adipocyte gene-set enrichment analysis
- CellChat-based analysis of adipocyte-derived secreted signaling
- comparison of healthy and ALL bone marrow adipocyte/stromal populations

## Data

Raw sequencing data, filtered 10x matrices, large Seurat objects, and generated output files are not included in this GitHub repository.

The expected local data structure is:

    BoneMarrowAdipocytes_analysis/
    ├── data/
    │   ├── filtered_matrices/
    │   │   ├── Healthy_1/
    │   │   ├── Healthy_2/
    │   │   ├── Healthy_3/
    │   │   └── ALL_Patient_SN3_D33/
    │   └── public_reference/
    │       ├── GSE253355/
    │       └── GSE225700/
    ├── scripts/
    ├── results/
    └── figures/

Each sample folder in data/filtered_matrices/ should contain:

    barcodes.tsv.gz
    features.tsv.gz
    matrix.mtx.gz

## Analysis workflow

The scripts are numbered in the order in which they were run.

| Script | Purpose |
|---|---|
| 01_read_10x_create_objects.R | Import Cell Ranger filtered matrices and create initial Seurat objects |
| 02_qc_scDblFinder_filtering.R | Perform QC, doublet detection, and final sample-level filtering |
| 03_integrate_healthy.R | Integrate healthy donor bone marrow samples |
| 04_umap_singler_healthy.R | Annotate healthy integrated object using SingleR |
| 05_msc_subclustering_healthy.R | Recluster and manually annotate healthy mesenchymal/stromal populations |
| 06_final_healthy_figures_tables.R | Generate final healthy donor figures and composition tables |
| 07_integrate_MSC_reference.R | Integrate healthy bone marrow stromal cells with public MSC reference dataset GSE253355 |
| 08_integrate_adipose_reference.R | Integrate bone marrow adipocytes with public subcutaneous adipocytes from GSE225700 |
| 09_DEG_GSEA_BM_vs_SC_adipocytes.R | Perform RNA-based differential expression and custom adipocyte gene-set analysis |
| 10_CellChat_healthy_intrinsic.R | Analyze intrinsic adipocyte-derived secreted signaling using CellChat |
| 11_preprocess_singler_ALL.R | Process and annotate the ALL bone marrow sample |
| 12_msc_subclustering_ALL.R | Recluster and manually annotate ALL mesenchymal/stromal populations |
| 13_final_ALL_figures_tables_PLIN1.R | Generate final ALL figures, composition tables, and PLIN1 expression comparison |

## Software

The analysis was performed in R using the following main packages:

- Seurat v5
- scDblFinder
- SingleR
- celldex
- SingleCellExperiment
- clusterProfiler
- ggplot2
- CellChat

Cell Ranger Multi v9.0.1 was used for initial processing of raw sequencing data against the human GRCh38 2024-A reference genome.

## Differential expression

Differential expression between bone marrow adipocytes and subcutaneous adipocytes was performed on the log-normalized RNA assay rather than the integrated assay. The integrated assay was used for visualization and integration, whereas the RNA assay preserves the full gene set required for differential expression analysis.

## CellChat analysis

CellChat analysis was performed on the healthy bone marrow dataset using the manually curated final_annotation metadata column. The primary CellChat analysis used population.size = FALSE to estimate intrinsic signaling capacity independent of cell population abundance.

## Data availability

Raw FASTQ files, filtered Cell Ranger feature-barcode matrices, and processed metadata tables will be deposited in the Gene Expression Omnibus under accession number GSEXXXXX.

## Code availability

This repository contains the reproducible R analysis scripts for the manuscript-associated snRNA-seq analysis.
