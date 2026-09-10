# PBMC scRNA-seq & CellChat Analysis

## Project Overview

This project analyzes single-cell RNA-seq (scRNA-seq) data from a human peripheral blood mononuclear cell (PBMC) sample. It applies a standard Seurat-based preprocessing and clustering workflow, annotates the resulting clusters using three complementary approaches, and uses CellChat to infer potential cell-cell communication patterns between the identified populations.

## Biological Objective

What immune cell populations can be identified in a healthy human PBMC sample from single-cell transcriptomic data, and what ligand-receptor communication patterns can be inferred between them using CellChat?

## Dataset

- A single donor's PBMC sample, provided as a 10x Genomics filtered feature-barcode matrix (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`), loaded locally with `Seurat::Read10X()`.
- The data are treated as a locally supplied sample; no public GEO accession or DOI is associated with the dataset used here. Raw data are not included in this repository.
- 5,550 cells passed QC filtering; 5,152 cells remained after doublet removal and were carried through the rest of the analysis.

## Analysis Workflow

1. Data loading (`CreateSeuratObject` from a 10x filtered matrix)
2. Quality control and filtering
3. Doublet detection and removal (`scDblFinder`)
4. Normalization (log-normalization)
5. Highly variable feature selection (2,000 features)
6. Scaling and PCA
7. Neighbor graph construction and Louvain clustering
8. UMAP embedding
9. Cell-type annotation via canonical markers, cluster DEGs, and SingleR
10. Annotation comparison and refinement across the three approaches
11. CellChat cell-cell communication analysis (probability, pathway-level, and ligand-receptor analysis)
12. Export of figures and tables to `results/`

## Quality Control

Cells were filtered using standard PBMC QC thresholds:

- `nFeature_RNA > 200` (removes empty droplets / low-complexity cells)
- `nFeature_RNA < 6000` (removes likely doublets/multiplets)
- `percent.mt < 10` (removes low-quality/dying cells)

These thresholds are consistent with common Seurat/10x PBMC workflows. QC distributions were visualized both before and after filtering so the effect of these cutoffs on the cell population is directly inspectable.

- Cells after QC filtering: 5,550
- Cells after doublet removal (final dataset used downstream): 5,152

## Dimensionality Reduction and Clustering

- PCA was run on the 2,000 highly variable features, with an elbow plot used to guide dimensionality choice.
- Clustering used the first 10 principal components (`FindNeighbors`, `FindClusters`, resolution = 0.5), yielding 16 clusters (0-15).
- UMAP was computed on the same 10 principal components for 2D visualization of clusters and annotations.

## Cell-Type Annotation

Three complementary annotation strategies were applied and compared:

- **Canonical marker-based annotation**: broad lineage assignment (T cells, B cells, NK cells, Monocytes, Dendritic cells, Platelets, or Unresolved) based on well-established PBMC lineage markers (e.g., `CD3D`/`CD3E` for T cells, `MS4A1`/`CD79A` for B cells, `NKG7`/`GNLY` for NK cells, `LYZ`/`S100A8` for monocytes, `CD1C`/`FCER1A` for dendritic cells, `PPBP`/`PF4` for platelets).
- **Cluster DEG-based annotation**: a more granular label per cluster, derived from each cluster's top differentially expressed genes.
- **SingleR**: automated, reference-based annotation against the Human Primary Cell Atlas.

The three approaches were cross-tabulated to assess agreement and refine the final labels. Final annotation, cluster size, and representative supporting evidence:

| Cluster | Cells | Final label | Key supporting markers | SingleR concordance |
|---------|-------|--------------|--------------------------|----------------------|
| 0 | 965 | Naive CD4 T cells | `LEF1`, `NRCAM` | 965/965 T_cells |
| 1 | 683 | Classical Monocytes | `S100A8`, `S100A12`, `RBP7` | 683/683 Monocyte |
| 2 | 619 | Monocytes | `LGALS2`, `TMEM176A/B` | 619/619 Monocyte |
| 3 | 538 | Myeloid cells (generic) | non-canonical top DEGs; grouped with cluster 13 | mostly Monocyte |
| 4 | 468 | Gamma-delta Cytotoxic T cells | `KLRC2`, `GZMH`, `TRDC`, `TRGC2` | split T_cells/NK_cell (240/227) |
| 5 | 459 | Naive B cells | `TCL1A`, `IGHD` | 459/459 B_cell |
| 6 | 397 | Regulatory / Activated CD4 T cells | `CCR4`, `IL2RA`, `CTLA4` | 395/397 T_cells |
| 7 | 268 | NK cells | `KIR3DL1`, `XCL2`, `CCL3`, `CD160` | 268/268 NK_cell |
| 8 | 161 | NK-like cells | `NCAM1`, `BNC2` | 156/161 NK_cell |
| 9 | 141 | Plasma cells | `IGHA1`, `IGHG1` | 141/141 B_cell (no plasma-cell category in reference) |
| 10 | 138 | C1Q+ Monocytes | `C1QA` | 138/138 Monocyte |
| 11 | 130 | MAIT / CD161-like T cells | `SLC4A10`, `RORC` | 124/130 T_cells |
| 12 | 87 | Conventional Dendritic Cells | `FCER1A`, `CD1C`, `CD1E`, `CLEC10A` | 86/87 Monocyte (no DC-specific match) |
| 13 | 45 | Myeloid cells (generic) | non-canonical top DEGs; grouped with cluster 3 | mostly Monocyte |
| 14 | 38 | Unresolved | non-canonical/low-specificity DEGs | mixed, no dominant reference label |
| 15 | 15 | Platelets | `ITGB3`, `TUBB1`, `CAVIN2`, `CMTM5` | 14/15 Platelets |

Clusters 3 and 13 retain the generic "Myeloid cells" label: their top differentially expressed genes were not canonical monocyte/DC markers, so a more specific identity was not assigned, even though SingleR maps the majority of this combined group to `Monocyte`. Cluster 14 remains labeled "Unresolved" for the same reason, no consistent marker or reference signal supported a confident call.

## CellChat Analysis

CellChat was built from the final resolved cluster identities. The workflow used the human CellChat ligand-receptor database and covered:

- Identification of overexpressed genes and ligand-receptor interactions per group
- Communication probability inference (`computeCommunProb`)
- Pathway-level communication probability (`computeCommunProbPathway`)
- Aggregate communication network (interaction count and strength)
- Pathway centrality (sender/receiver roles)
- Source and target rankings on the resolved network (Unresolved cluster removed)

68 signaling pathways were detected in total. Signaling pathways were examined in greater detail based on their relevance to the observed immune populations and inferred communication patterns: **MHC-I, MHC-II, MIF, CLEC, GALECTIN, and TGFb**. This does not imply that pathways outside this set are biologically unimportant, only that they were not examined in further depth here.

### Platelet Population

The platelet cluster contains only 15 cells, well below the other resolved groups. A uniform CellChat minimum-cell filter (`filterCommunication(min.cells = 20)`) was applied across all groups, which excludes any population below this threshold from communication probability estimation. In this dataset, only the platelet population falls below the threshold; it remains a labeled population in the metadata but contributes no communication estimates in the final network.

## Key Results

- 14 resolved cell populations (16 clusters, 15 annotated types, 14 after removing "Unresolved") were used in the final CellChat network.
- By outgoing signal (source ranking), Conventional Dendritic Cells, Monocytes, C1Q+ Monocytes, and Classical Monocytes rank highest. Myeloid populations are the strongest overall senders in this network.
- By incoming signal (target ranking), NK cells, MAIT/CD161-like T cells, Monocytes, and Gamma-delta Cytotoxic T cells rank highest.
- Among the top resolved ligand-receptor interactions, the CLEC pathway (`CD69`-`KLRB1`, `CLEC2B`-`KLRB1`, `CLEC2D`-`KLRB1`) dominates, largely directed toward MAIT/CD161-like T cells and NK cells from multiple lymphoid sources.
- MHC-I signaling (`HLA-E`-`KLRK1`) is prominent between NK cells, Gamma-delta Cytotoxic T cells, and several monocyte populations, consistent with NK/cytotoxic-T inhibitory receptor engagement.
- Myeloid populations (Classical Monocytes, Myeloid cells) show notable APP-`CD74` and VISFATIN (`NAMPT`-`ITGA5_ITGB1`) signaling among themselves and toward Conventional Dendritic Cells and Naive B cells.
- The platelet population contributes no interactions to the final resolved network under the applied filtering threshold.

## Limitations

- Single donor/sample; findings are not generalizable across individuals without additional biological replicates.
- CellChat output represents predicted, ligand-receptor-expression-based communication probability, not experimentally validated signaling.
- The Human Primary Cell Atlas reference used by SingleR does not include fine-grained categories for some populations (e.g., plasma cells, gamma-delta/MAIT T cells, dendritic-cell subsets), so SingleR concordance for these groups is approximate.
- Clusters 3 and 13 remain under the generic "Myeloid cells" label, and cluster 14 remains "Unresolved," reflecting conservative calls where marker evidence was insufficient for a more specific identity.
- The 15-cell platelet population is small enough that communication probability estimates would be unreliable; it was excluded from CellChat's communication estimation via a uniform minimum cell-count filter.

## Reproducibility

- Analysis performed in R using Seurat, dplyr, ggplot2, patchwork, SingleCellExperiment, scDblFinder, SingleR, celldex, and CellChat.
- Random seed fixed with `set.seed(100)`.
- Full R session and package version information is written to `results/sessionInfo.txt`.
- The 10x filtered feature-barcode matrix is expected under `data/` and is not included in this repository; users must supply their own matrix at the path referenced in the script.
- The full analysis script is provided as a single, runnable `.R` file.

## Repository Structure

```
PBMC-scRNAseq-CellChat-Analysis/
├── README.md
├── .gitignore
├── scripts/
│   └── PBMC scRNA-seq Analysis + CellChat.R
└── results/
    ├── sessionInfo.txt
    ├── figures/
    │   ├── CellChat/
    │   │   └── ...
    │   └── ...
    └── tables/
        ├── CellChat/
        │   └── ...
        └── ...
```

## .gitignore

```
data/
*.rds
.Rhistory
.RData
```

Raw/local input data and large R data objects are excluded from version control.
This project was developed as a portfolio/research re-analysis project to demonstrate practical skills in single-cell transcriptomics and downstream biological interpretation.
