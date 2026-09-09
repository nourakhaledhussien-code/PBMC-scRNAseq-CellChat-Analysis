# PBMC-scRNAseq-CellChat-Analysis

## Project Overview

This project presents a complete single-cell RNA sequencing (scRNA-seq) analysis of a human Peripheral Blood Mononuclear Cell (PBMC) dataset, including quality control, dimensionality reduction, clustering, cell-type annotation, annotation validation, and cell–cell communication analysis using CellChat.

The analysis was performed in **R** using the **Seurat** framework, with multiple approaches used for cell-type annotation and **CellChat** used to investigate intercellular communication networks.

---

## Biological Objective

The main objectives of this project are to:

* Perform quality control and preprocessing of PBMC scRNA-seq data.
* Identify distinct cellular populations through dimensionality reduction and clustering.
* Annotate cell populations using multiple complementary approaches.
* Compare and validate cell-type annotations.
* Characterize cell–cell communication networks using CellChat.
* Identify major signaling pathways and ligand–receptor interactions between immune cell populations.

---

## Analysis Workflow

The analysis follows this workflow:

1. Data loading and preprocessing
2. Quality control (QC)
3. Doublet/singlet assessment
4. Normalization
5. Identification of highly variable genes
6. Scaling and PCA
7. UMAP dimensionality reduction
8. Clustering
9. Cell-type annotation
10. Annotation validation
11. CellChat analysis
12. Identification of major signaling pathways
13. Ligand–receptor interaction analysis
14. Visualization and result organization

---

## Quality Control

Quality control was performed to remove low-quality cells and reduce technical noise.

The QC process included evaluation of:

* Number of detected genes
* Number of RNA features
* Mitochondrial gene percentage
* Droplet/cell quality
* Singlet versus doublet populations

Cells that did not meet the defined quality criteria were excluded before downstream analysis.

---

## Dimensionality Reduction and Clustering

After preprocessing, dimensionality reduction and clustering were performed using Seurat.

The workflow included:

* Highly variable gene identification
* Data scaling
* Principal Component Analysis (PCA)
* Elbow plot evaluation
* UMAP visualization
* Graph-based clustering

These steps were used to identify transcriptionally distinct cellular populations within the PBMC dataset.

---

## Cell-Type Annotation

Cell identities were assigned using multiple complementary strategies.

### 1. Canonical Marker-Based Annotation

Known canonical marker genes were used to identify major immune cell populations based on their characteristic expression profiles.

### 2. DEG-Based Annotation

Differentially expressed genes were identified for each cluster and used to support biological interpretation and cell-type assignment.

### 3. SingleR Annotation

Automated cell-type annotation was performed using **SingleR** and reference transcriptomic datasets.

### 4. Annotation Validation

The different annotation approaches were compared to evaluate their consistency and identify potential differences between:

* Canonical marker annotation
* DEG-based annotation
* SingleR annotation

This multi-method strategy provides a more robust interpretation of cell identities than relying on a single annotation method.

---

## CellChat Analysis

Cell–cell communication was investigated using **CellChat**.

The analysis was used to identify:

* Communication networks between cell populations
* Communication probability and interaction strength
* Major signaling pathways
* Ligand–receptor interactions
* Important sender and receiver cell populations

Several signaling pathways were investigated, including pathways related to:

* MHC-I
* MHC-II
* MIF
* TGF-β
* CLEC
* Galectin

Both overall communication networks and pathway-specific interactions were examined.

---

## Key Results

The analysis generated:

* QC and filtering summaries
* PCA and UMAP visualizations
* Cluster-level marker genes
* Cell-type annotation maps
* Annotation comparison results
* CellChat communication networks
* Signaling pathway summaries
* Ligand–receptor interaction tables
* Resolved cell–cell communication results

The generated figures and tables are organized under the `results/` directory.

---

## Repository Structure

```text
PBMC-scRNAseq-CellChat-Analysis/
│
├── README.md
│
├── .gitignore
│
├── scripts/
│   └── PBMC scRNA-seq Analysis + CellChat.R
│
└── results/
    │
    ├── README.md
    │
    ├── figures/
    │   │
    │   ├── CellChat/
    │   │
    │   ├── PCA_plot.png
    │   ├── QC_after_filtering.png
    │   ├── UMAP_clusters.png
    │   ├── UMAP_canonical_annotation.png
    │   ├── UMAP_DEG_annotation.png
    │   ├── UMAP_SingleR_annotation.png
    │   ├── annotation_comparison.png
    │   ├── canonical_markers.png
    │   ├── highly_variable_genes.png
    │   └── top_cluster_markers_heatmap.png
    │
    └── tables/
        │
        ├── CellChat/
        │
        ├── analysis_summary.csv
        ├── QC_cell_count.csv
        ├── cluster_cell_count.csv
        ├── cluster_markers.csv
        ├── highly_variable_genes.csv
        ├── canonical_vs_DEG_annotation.csv
        └── DEG_vs_SingleR_annotation.csv
```

---

## Tools and Packages

The analysis was performed using R and the following major packages:

* **Seurat**
* **SeuratObject**
* **SingleR**
* **celldex**
* **CellChat**
* **dplyr**
* **ggplot2**
* **patchwork**

---

## Reproducibility

The complete analysis workflow is provided in:

```text
scripts/PBMC scRNA-seq Analysis + CellChat.R
```

The script contains the main preprocessing, clustering, annotation, validation, and CellChat analysis steps used to generate the reported results.

Raw sequencing data and large intermediate objects are not included in this repository.

---

## Project Type

**Single-cell RNA-seq | Immunology | Cell-Type Annotation | Cell–Cell Communication | CellChat | R**

This project was developed as a portfolio/research re-analysis project to demonstrate practical skills in single-cell transcriptomics and downstream biological interpretation.
