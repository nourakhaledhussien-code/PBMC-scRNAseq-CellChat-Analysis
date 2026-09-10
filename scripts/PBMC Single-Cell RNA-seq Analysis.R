# ============================================================
# PBMC Single-Cell RNA-seq Analysis
# Donor 1
# ============================================================
#
# Workflow:
# QC Filtering
# → Doublet Detection
# → Normalization
# → HVG Selection
# → Scaling
# → PCA
# → Clustering
# → UMAP
# → Canonical Annotation
# → DEG/Marker Analysis
# → Annotation Validation
# → SingleR
# → CellChat
# → Pathway Analysis
# → Sender/Receiver Analysis
# → Ligand-Receptor Analysis
#
# Dataset:
# Healthy PBMC — Donor 1
#
# Dataset source (GEO accession / 10x Genomics dataset URL / DOI):
# [TODO: INSERT THE OFFICIAL DATASET SOURCE IDENTIFIER HERE]
# (Not fabricated — this must be filled in with the actual
#  accession/link/citation for the data used, before the README
#  or repository is finalized.)
#
# Biological Question:
# What immune cell populations can be identified in a healthy
# PBMC sample, and what potential cell-cell communication
# patterns can be detected between these populations?
#
# ============================================================


# ============================================================
# 1. Load Required Packages
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

library(SingleCellExperiment)
library(scDblFinder)

library(SingleR)
library(celldex)

library(CellChat)


# ============================================================
# 2. Reproducibility
# ============================================================

set.seed(100)


# ============================================================
# 3. Create Output Directories
# ============================================================
#
# Directory layout matches the GitHub repository structure:
#
# results/
# ├── figures/
# │   ├── CellChat/      <- CellChat-specific figures
# │   └── (general figures)
# └── tables/
#     ├── CellChat/      <- CellChat-specific tables
#     └── (general tables)
# ============================================================

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "results/figures",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "results/figures/CellChat",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "results/tables",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "results/tables/CellChat",
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 4. Load PBMC Donor 1
# ============================================================

# Place the 10X filtered matrix inside the data/ directory.
# See dataset source note at the top of this script — the
# official accession/link must be inserted there.

data_dir <- "/Volumes/NGSـprojects/Abd_alfataah_academy/Course_project/sample_filtered_feature_bc_matrix"
pbmc_d1 <- CreateSeuratObject(
  counts = Read10X(data_dir),
  project = "PBMC_Donor1"
)


cat("\n============================================\n")
cat("Initial number of cells:", ncol(pbmc_d1), "\n")
cat("============================================\n")


# ============================================================
# 5. Quality Control
# ============================================================

pbmc_d1[["percent.mt"]] <- PercentageFeatureSet(
  pbmc_d1,
  pattern = "^MT-"
)


# ------------------------------------------------------------
# 5.0 QC Visualization — Before Filtering
# ------------------------------------------------------------
#
# Shown for transparency/documentation purposes only. This
# reflects the FULL pre-filtering distribution and is not a
# filtering step itself. Compare against QC_after_filtering.png
# to see the effect of the thresholds below.
# ------------------------------------------------------------

p_qc_before <- VlnPlot(
  pbmc_d1,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
) +
  plot_annotation(
    title = "PBMC Donor 1 — QC Before Filtering"
  )

ggsave(
  "results/figures/QC_before_filtering.png",
  p_qc_before,
  width = 10,
  height = 4
)


# ------------------------------------------------------------
# 5.1 QC Filtering
# ------------------------------------------------------------
#
# Thresholds used:
#   nFeature_RNA > 200   — removes empty droplets / low-complexity cells
#   nFeature_RNA < 6000  — removes likely doublets / multiplets
#   percent.mt   < 10    — removes dying / low-quality cells
#
# These are standard, widely used PBMC QC thresholds (consistent
# with common Seurat/10x Genomics PBMC workflows). They are a
# reasonable, data-dependent default rather than an arbitrary or
# broken cap: the upper nFeature_RNA bound is not a visualization
# limit, it is an actual filtering threshold intended to remove
# likely multiplets, and it is documented here explicitly rather
# than left unexplained. It was not changed as part of this
# correction pass because no evidence was provided that it is
# inappropriate for this dataset; QC_before_filtering.png and
# QC_after_filtering.png together let this be checked visually
# after rerunning.
# ------------------------------------------------------------

pbmc_d1 <- subset(
  pbmc_d1,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)


# ------------------------------------------------------------
# 5.2 QC Summary
# ------------------------------------------------------------

qc_summary <- data.frame(
  stage = "After QC filtering",
  cells = ncol(pbmc_d1)
)

print(qc_summary)

write.csv(
  qc_summary,
  "results/tables/QC_cell_count.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 5.3 QC Visualization — After Filtering
# ------------------------------------------------------------

p_qc_after <- VlnPlot(
  pbmc_d1,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
) +
  plot_annotation(
    title = "PBMC Donor 1 — QC After Filtering"
  )

ggsave(
  "results/figures/QC_after_filtering.png",
  p_qc_after,
  width = 10,
  height = 4
)


# ============================================================
# 6. Doublet Detection
# ============================================================

sce_d1 <- as.SingleCellExperiment(pbmc_d1)

sce_d1 <- scDblFinder(
  sce_d1
)


# ------------------------------------------------------------
# 6.1 Transfer Doublet Information
# ------------------------------------------------------------

pbmc_d1$doublet_score <-
  colData(sce_d1)$scDblFinder.score

pbmc_d1$doublet_class <-
  colData(sce_d1)$scDblFinder.class


# ------------------------------------------------------------
# 6.2 Doublet Summary
# ------------------------------------------------------------

doublet_summary <- as.data.frame(
  table(pbmc_d1$doublet_class)
)

colnames(doublet_summary) <- c(
  "classification",
  "cells"
)

print(doublet_summary)

write.csv(
  doublet_summary,
  "results/tables/doublet_summary.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 6.3 Doublet Visualization
# ------------------------------------------------------------

p_doublet <- ggplot(
  doublet_summary,
  aes(
    x = classification,
    y = cells
  )
) +
  geom_col() +
  labs(
    title = "Doublet Detection Summary",
    x = "Classification",
    y = "Number of Cells"
  ) +
  theme_classic()

ggsave(
  "results/figures/doublet_summary.png",
  p_doublet,
  width = 7,
  height = 5
)


# ------------------------------------------------------------
# 6.4 Remove Doublets
# ------------------------------------------------------------

pbmc_d1 <- subset(
  pbmc_d1,
  subset = doublet_class == "singlet"
)


# ------------------------------------------------------------
# 6.5 Singlet Summary
# ------------------------------------------------------------

singlet_summary <- data.frame(
  stage = "After doublet removal",
  cells = ncol(pbmc_d1)
)

print(singlet_summary)

write.csv(
  singlet_summary,
  "results/tables/singlet_cell_count.csv",
  row.names = FALSE
)


# ============================================================
# 7. Normalization
# ============================================================

pbmc_d1 <- NormalizeData(
  pbmc_d1,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# ============================================================
# 8. Highly Variable Genes
# ============================================================

pbmc_d1 <- FindVariableFeatures(
  pbmc_d1,
  selection.method = "vst",
  nfeatures = 2000
)


# ------------------------------------------------------------
# 8.1 Save HVGs
# ------------------------------------------------------------

hvg_table <- data.frame(
  gene = VariableFeatures(pbmc_d1)
)

write.csv(
  hvg_table,
  "results/tables/highly_variable_genes.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 8.2 HVG Visualization
# ------------------------------------------------------------

p_hvg <- VariableFeaturePlot(
  pbmc_d1
) +
  ggtitle(
    "Highly Variable Genes"
  )

ggsave(
  "results/figures/highly_variable_genes.png",
  p_hvg,
  width = 8,
  height = 6
)


# ============================================================
# 9. Scaling
# ============================================================

pbmc_d1 <- ScaleData(
  pbmc_d1,
  features = VariableFeatures(pbmc_d1)
)


# ============================================================
# 10. Principal Component Analysis
# ============================================================

pbmc_d1 <- RunPCA(
  pbmc_d1,
  features = VariableFeatures(pbmc_d1)
)


# ------------------------------------------------------------
# 10.1 PCA Elbow Plot
# ------------------------------------------------------------

p_elbow <- ElbowPlot(
  pbmc_d1,
  ndims = 30
) +
  ggtitle(
    "PCA Elbow Plot"
  )

ggsave(
  "results/figures/PCA_elbow_plot.png",
  p_elbow,
  width = 7,
  height = 5
)


# ------------------------------------------------------------
# 10.2 PCA Visualization
# ------------------------------------------------------------

p_pca <- DimPlot(
  pbmc_d1,
  reduction = "pca"
) +
  ggtitle(
    "PCA — PBMC Donor 1"
  )

ggsave(
  "results/figures/PCA_plot.png",
  p_pca,
  width = 8,
  height = 6
)


# ============================================================
# 11. Neighbor Graph and Clustering
# ============================================================

pbmc_d1 <- FindNeighbors(
  pbmc_d1,
  dims = 1:10
)

pbmc_d1 <- FindClusters(
  pbmc_d1,
  resolution = 0.5
)


# ------------------------------------------------------------
# 11.1 Cluster Cell Counts
# ------------------------------------------------------------

cluster_counts <- as.data.frame(
  table(pbmc_d1$seurat_clusters)
)

colnames(cluster_counts) <- c(
  "cluster",
  "cells"
)

print(cluster_counts)

write.csv(
  cluster_counts,
  "results/tables/cluster_cell_counts.csv",
  row.names = FALSE
)


# ============================================================
# 12. UMAP
# ============================================================

pbmc_d1 <- RunUMAP(
  pbmc_d1,
  dims = 1:10
)


# ------------------------------------------------------------
# 12.1 UMAP — Unsupervised Clustering
# ------------------------------------------------------------

p_umap_clusters <- DimPlot(
  pbmc_d1,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "PBMC Donor 1 — Unsupervised Clustering"
  )

ggsave(
  "results/figures/UMAP_clusters.png",
  p_umap_clusters,
  width = 8,
  height = 6
)


# ============================================================
# 13. Canonical Marker-Based Annotation
# ============================================================

canonical_markers <- c(
  
  # T cells
  "CD3D",
  "CD3E",
  "TRBC1",
  "TRBC2",
  
  # NK cells
  "NKG7",
  "GNLY",
  "KLRD1",
  
  # B cells
  "MS4A1",
  "CD79A",
  "CD74",
  "HLA-DRA",
  
  # Monocytes
  "LYZ",
  "S100A8",
  "S100A9",
  "LST1",
  
  # Dendritic cells
  "CD1C",
  "FCER1A",
  "CLEC10A",
  
  # Platelets
  "PPBP",
  "PF4",
  "GNG11"
)


# ------------------------------------------------------------
# 13.1 Canonical Marker Visualization
# ------------------------------------------------------------

canonical_features <- intersect(
  c(
    "CD3D",
    "MS4A1",
    "NKG7",
    "LYZ",
    "CD14",
    "FCGR3A",
    "CD1C",
    "PPBP"
  ),
  rownames(pbmc_d1)
)

p_canonical_markers <- FeaturePlot(
  pbmc_d1,
  features = canonical_features,
  ncol = 4
)

ggsave(
  "results/figures/canonical_markers.png",
  p_canonical_markers,
  width = 12,
  height = 8
)


# ------------------------------------------------------------
# 13.2 Assign Major Cell Types
# ------------------------------------------------------------
#
# CORRECTED: clusters 4 and 5 were previously swapped.
#   Cluster 4 marker evidence (KLRC2, GZMH, TRDC) is
#   cytotoxic/T-cell-associated -> "T cells".
#   Cluster 5 marker evidence (MS4A1, IGHD, TCL1A), together
#   with SingleR classifying 459/459 of this cluster's cells
#   as "B_cell" -> "B cells".
# This is fixed here at the source (the recode mapping), so
# it propagates through the annotation comparison, SingleR
# comparison, and the CellChat object built from these labels.
# ------------------------------------------------------------

pbmc_d1$canonical_annotation <- dplyr::recode(
  as.character(pbmc_d1$seurat_clusters),
  
  `0`  = "T cells",
  `1`  = "Monocytes",
  `2`  = "Monocytes",
  `3`  = "Myeloid cells",
  `4`  = "T cells",
  `5`  = "B cells",
  `6`  = "T cells",
  `7`  = "NK cells",
  `8`  = "NK cells",
  `9`  = "B cells",
  `10` = "Myeloid cells",
  `11` = "T cells",
  `12` = "Dendritic cells",
  `13` = "Myeloid cells",
  `14` = "Unresolved",
  `15` = "Platelets"
)


# ------------------------------------------------------------
# 13.3 Canonical Annotation UMAP
# ------------------------------------------------------------

p_canonical <- DimPlot(
  pbmc_d1,
  reduction = "umap",
  group.by = "canonical_annotation",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "Canonical Marker-Based Annotation"
  )

ggsave(
  "results/figures/UMAP_canonical_annotation.png",
  p_canonical,
  width = 9,
  height = 6
)


# ============================================================
# 14. Differentially Expressed Marker Analysis
# ============================================================

markers_d1 <- FindAllMarkers(
  pbmc_d1,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


# ------------------------------------------------------------
# 14.1 Save All Markers
# ------------------------------------------------------------

write.csv(
  markers_d1,
  "results/tables/cluster_markers.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 14.2 Top 5 Markers per Cluster
# ------------------------------------------------------------

top5_markers_d1 <- markers_d1 %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = 5
  ) %>%
  dplyr::arrange(
    cluster,
    dplyr::desc(avg_log2FC)
  )

print(top5_markers_d1)

write.csv(
  top5_markers_d1,
  "results/tables/top5_markers_per_cluster.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 14.3 Marker Heatmap
# ------------------------------------------------------------

heatmap_features <- unique(
  top5_markers_d1$gene
)

heatmap_features <- intersect(
  heatmap_features,
  rownames(pbmc_d1)
)

if (length(heatmap_features) > 0) {
  
  p_marker_heatmap <- DoHeatmap(
    pbmc_d1,
    features = heatmap_features,
    size = 3
  ) +
    ggtitle(
      "Top Differentially Expressed Marker Genes"
    )
  
  ggsave(
    "results/figures/top_cluster_markers_heatmap.png",
    p_marker_heatmap,
    width = 12,
    height = 10
  )
}


# ============================================================
# 15. DEG-Based Cell-Type Annotation
# ============================================================
#
# CORRECTED: clusters 4 and 5 were previously swapped (same
# issue as section 13.2). Cluster 4 -> "Gamma-delta Cytotoxic
# T cells" (KLRC2/GZMH/TRDC), cluster 5 -> "Naive B cells"
# (MS4A1/IGHD/TCL1A; SingleR: 459/459 "B_cell"). This is the
# identity table used downstream by CellChat (group.by =
# "DEG_annotation"), so correcting it here propagates the fix
# into CellChat automatically.
#
# NOTE ON OTHER CLUSTERS: clusters 3, 10, and 13 are currently
# labeled generically ("Myeloid cells" / "C1Q+ Monocytes").
# These were left as-is because no additional marker evidence
# for a more specific identity was available in this pass.
# After rerunning, check cluster_markers.csv / 
# top5_markers_per_cluster.csv and SingleR_results.csv for
# clusters 3 and 13 specifically — if they show clear,
# consistent monocyte/DC subtype markers, they can be made more
# specific at that point; otherwise they should remain generic
# rather than be assigned an unsupported identity.
# ------------------------------------------------------------

pbmc_d1$DEG_annotation <- dplyr::recode(
  as.character(pbmc_d1$seurat_clusters),
  
  `0`  = "Naive CD4 T cells",
  `1`  = "Classical Monocytes",
  `2`  = "Monocytes",
  `3`  = "Myeloid cells",
  `4`  = "Gamma-delta Cytotoxic T cells",
  `5`  = "Naive B cells",
  `6`  = "Regulatory / Activated CD4 T cells",
  `7`  = "NK cells",
  `8`  = "NK-like cells",
  `9`  = "Plasma cells",
  `10` = "C1Q+ Monocytes",
  `11` = "MAIT / CD161-like T cells",
  `12` = "Conventional Dendritic Cells",
  `13` = "Myeloid cells",
  `14` = "Unresolved",
  `15` = "Platelets"
)


# ------------------------------------------------------------
# 15.1 DEG Annotation UMAP
# ------------------------------------------------------------

p_deg <- DimPlot(
  pbmc_d1,
  reduction = "umap",
  group.by = "DEG_annotation",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "Differentially Expressed Marker-Based Annotation"
  )

ggsave(
  "results/figures/UMAP_DEG_annotation.png",
  p_deg,
  width = 10,
  height = 7
)


# ============================================================
# 16. Compare Annotation Approaches
# ============================================================

canonical_vs_deg <- table(
  Canonical = pbmc_d1$canonical_annotation,
  DEG = pbmc_d1$DEG_annotation
)

print(canonical_vs_deg)

write.csv(
  as.data.frame(canonical_vs_deg),
  "results/tables/canonical_vs_DEG_annotation.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 16.1 Visual Comparison
# ------------------------------------------------------------

p_annotation_comparison <-
  p_canonical + p_deg

ggsave(
  "results/figures/annotation_comparison.png",
  p_annotation_comparison,
  width = 14,
  height = 6
)


# ============================================================
# 17. Automated Annotation Using SingleR
# ============================================================

hpca_ref <-
  celldex::HumanPrimaryCellAtlasData()

SCE_d1 <-
  as.SingleCellExperiment(pbmc_d1)

singler_d1 <- SingleR(
  test = SCE_d1,
  ref = hpca_ref,
  labels = hpca_ref$label.main
)

pbmc_d1$SingleR_annotation <-
  singler_d1$labels


# ------------------------------------------------------------
# 17.1 SingleR UMAP
# ------------------------------------------------------------

p_singler <- DimPlot(
  pbmc_d1,
  reduction = "umap",
  group.by = "SingleR_annotation",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "SingleR Automated Annotation"
  )

ggsave(
  "results/figures/UMAP_SingleR_annotation.png",
  p_singler,
  width = 10,
  height = 7
)


# ------------------------------------------------------------
# 17.2 DEG vs SingleR Comparison
# ------------------------------------------------------------

deg_vs_singler <- table(
  DEG = pbmc_d1$DEG_annotation,
  SingleR = pbmc_d1$SingleR_annotation
)

print(deg_vs_singler)

write.csv(
  as.data.frame(deg_vs_singler),
  "results/tables/DEG_vs_SingleR_annotation.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 17.3 SingleR Scores
# ------------------------------------------------------------

singler_scores <- as.data.frame(
  singler_d1
)

write.csv(
  singler_scores,
  "results/tables/SingleR_results.csv",
  row.names = FALSE
)


# ============================================================
# 18. CellChat Analysis
# ============================================================
#
# Built directly from the CORRECTED "DEG_annotation" metadata,
# so the corrected cluster 4/5 identities (and all other
# labels) flow into every downstream CellChat step.
# ============================================================

cellchat <- createCellChat(
  object = pbmc_d1,
  group.by = "DEG_annotation"
)


# ------------------------------------------------------------
# 18.1 Load Human CellChat Database
# ------------------------------------------------------------

cellchat@DB <- CellChatDB.human


# ------------------------------------------------------------
# 18.2 Prepare Signaling Data
# ------------------------------------------------------------

cellchat <- subsetData(
  cellchat
)


# ------------------------------------------------------------
# 18.3 Identify Overexpressed Signaling Genes
# ------------------------------------------------------------

cellchat <- identifyOverExpressedGenes(
  cellchat
)


# ------------------------------------------------------------
# 18.4 Identify Ligand-Receptor Interactions
# ------------------------------------------------------------

cellchat <- identifyOverExpressedInteractions(
  cellchat
)


# ------------------------------------------------------------
# 18.5 Infer Communication Probability
# ------------------------------------------------------------

cellchat <- computeCommunProb(
  cellchat
)


# ------------------------------------------------------------
# 18.6 Remove Very Small Cell Groups
# ------------------------------------------------------------
#
# The platelet population (~15 cells) is small enough that its
# communication probability estimates are statistically
# unstable and can disproportionately dominate top
# ligand-receptor results.
#
# DECISION: the platelet population is EXCLUDED from CellChat's
# statistical inference (its interactions are zeroed out) by
# raising min.cells from CellChat's default of 10 to 20. This
# threshold is applied uniformly to ALL cell groups (it is not
# a rule written specifically to target platelets), so any
# group below 20 cells would be excluded the same way. The
# platelet cluster itself is NOT deleted from the metadata or
# from cellchat@idents — only its contribution to the
# communication network is removed due to insufficient cell
# number for reliable estimation. This limitation should be
# stated explicitly in the final README.
# ------------------------------------------------------------

cellchat <- filterCommunication(
  cellchat,
  min.cells = 20
)


# ------------------------------------------------------------
# 18.7 Compute Pathway-Level Communication
# ------------------------------------------------------------

cellchat <- computeCommunProbPathway(
  cellchat
)


# ------------------------------------------------------------
# 18.8 Aggregate Communication Network
# ------------------------------------------------------------

cellchat <- aggregateNet(
  cellchat
)


# ============================================================
# 19. CellChat Overall Communication Network
# ============================================================

# ------------------------------------------------------------
# 19.1 Number of Inferred Interactions
# ------------------------------------------------------------

png(
  "results/figures/CellChat/CellChat_communication_count_circle.png",
  width = 8,
  height = 8,
  units = "in",
  res = 300
)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = as.numeric(
    table(cellchat@idents)
  ),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of inferred interactions"
)

dev.off()


# ------------------------------------------------------------
# 19.2 Strength of Inferred Interactions
# ------------------------------------------------------------

png(
  "results/figures/CellChat/CellChat_communication_strength_circle.png",
  width = 8,
  height = 8,
  units = "in",
  res = 300
)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = as.numeric(
    table(cellchat@idents)
  ),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Strength of inferred interactions"
)

dev.off()


# ------------------------------------------------------------
# 19.3 Save Communication Matrices
# ------------------------------------------------------------

communication_count <-
  cellchat@net$count

write.csv(
  communication_count,
  "results/tables/CellChat/CellChat_communication_count.csv"
)


communication_weight <-
  cellchat@net$weight

write.csv(
  communication_weight,
  "results/tables/CellChat/CellChat_communication_weight.csv"
)


# ============================================================
# 20. Ligand-Receptor Results
# ============================================================

lr_results <- subsetCommunication(
  cellchat
)


# ------------------------------------------------------------
# 20.1 Save All Ligand-Receptor Results
# ------------------------------------------------------------

write.csv(
  lr_results,
  "results/tables/CellChat/CellChat_ligand_receptor_results.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 20.2 Top Ligand-Receptor Interactions
# ------------------------------------------------------------

if (nrow(lr_results) > 0) {
  
  top_lr <- lr_results %>%
    dplyr::arrange(
      dplyr::desc(prob)
    ) %>%
    dplyr::select(
      source,
      target,
      ligand,
      receptor,
      prob,
      pval,
      pathway_name
    ) %>%
    dplyr::slice_head(
      n = 20
    )
  
  print(top_lr)
  
  write.csv(
    top_lr,
    "results/tables/CellChat/CellChat_top20_ligand_receptor_interactions.csv",
    row.names = FALSE
  )
}


# ============================================================
# 21. Available Signaling Pathways
# ============================================================

available_pathways <-
  cellchat@netP$pathways

cat(
  "\nNumber of CellChat signaling pathways:",
  length(available_pathways),
  "\n"
)

print(available_pathways)

write.csv(
  data.frame(
    pathway = available_pathways
  ),
  "results/tables/CellChat/CellChat_signaling_pathways.csv",
  row.names = FALSE
)


# ============================================================
# 22. Pathway-Level Probability Analysis
# ============================================================

if (!is.null(cellchat@netP$prob)) {
  
  pathway_strength <- apply(
    cellchat@netP$prob,
    3,
    sum,
    na.rm = TRUE
  )
  
  pathway_strength <- sort(
    pathway_strength,
    decreasing = TRUE
  )
  
  pathway_summary <- data.frame(
    pathway = names(pathway_strength),
    total_probability = as.numeric(
      pathway_strength
    )
  )
  
  pathway_summary$rank <-
    seq_len(nrow(pathway_summary))
  
  
  write.csv(
    pathway_summary,
    "results/tables/CellChat/CellChat_pathway_strength_summary.csv",
    row.names = FALSE
  )
  
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "TOP CELLCHAT SIGNALING PATHWAYS\n"
  )
  
  cat(
    "============================================\n"
  )
  
  print(
    head(
      pathway_summary,
      20
    )
  )
  
} else {
  
  warning(
    "cellchat@netP$prob is NULL. Pathway analysis skipped."
  )
}


# ============================================================
# 23. Compute Pathway Centrality
# ============================================================

cellchat <- netAnalysis_computeCentrality(
  cellchat,
  slot.name = "netP"
)


saveRDS(
  cellchat@netP$centr,
  "results/cellchat_pathway_centrality.rds"
)


# ============================================================
# 24. Define Main Biological Pathways
# ============================================================

main_pathways <- c(
  "MHC-I",
  "MHC-II",
  "MIF",
  "CLEC"
)

secondary_pathways <- c(
  "GALECTIN",
  "TGFb"
)


main_pathways_present <- intersect(
  main_pathways,
  cellchat@netP$pathways
)

secondary_pathways_present <- intersect(
  secondary_pathways,
  cellchat@netP$pathways
)


cat(
  "\nMain pathways present:\n"
)

print(
  main_pathways_present
)


cat(
  "\nSecondary pathways present:\n"
)

print(
  secondary_pathways_present
)


write.csv(
  data.frame(
    pathway = main_pathways_present
  ),
  "results/tables/CellChat/main_CellChat_pathways.csv",
  row.names = FALSE
)


write.csv(
  data.frame(
    pathway = secondary_pathways_present
  ),
  "results/tables/CellChat/secondary_CellChat_pathways.csv",
  row.names = FALSE
)


# ============================================================
# 25. Pathway Sender / Receiver Analysis
# ============================================================

sender_receiver_results <- list()


for (pathway in main_pathways_present) {
  
  cat(
    "\n\n====================================\n"
  )
  
  cat(
    "PATHWAY:",
    pathway,
    "\n"
  )
  
  cat(
    "====================================\n"
  )
  
  
  centrality <- cellchat@netP$centr[[pathway]]
  
  
  # Sender strength
  sender_values <- sort(
    centrality$outdeg,
    decreasing = TRUE
  )
  
  
  # Receiver strength
  receiver_values <- sort(
    centrality$indeg,
    decreasing = TRUE
  )
  
  
  cat(
    "\nTop Senders:\n"
  )
  
  print(
    head(
      sender_values,
      15
    )
  )
  
  
  cat(
    "\nTop Receivers:\n"
  )
  
  print(
    head(
      receiver_values,
      15
    )
  )
  
  
  sender_df <- data.frame(
    cell_type = names(sender_values),
    sender_strength = as.numeric(
      sender_values
    ),
    pathway = pathway,
    row.names = NULL
  )
  
  
  receiver_df <- data.frame(
    cell_type = names(receiver_values),
    receiver_strength = as.numeric(
      receiver_values
    ),
    pathway = pathway,
    row.names = NULL
  )
  
  
  sender_receiver_results[[pathway]] <- list(
    senders = sender_df,
    receivers = receiver_df
  )
  
  
  write.csv(
    sender_df,
    paste0(
      "results/tables/CellChat/",
      "CellChat_",
      pathway,
      "_senders.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    receiver_df,
    paste0(
      "results/tables/CellChat/",
      "CellChat_",
      pathway,
      "_receivers.csv"
    ),
    row.names = FALSE
  )
}


saveRDS(
  sender_receiver_results,
  "results/cellchat_sender_receiver_results.rds"
)


# ============================================================
# 26. Resolved Ligand-Receptor Analysis
# ============================================================

lr_resolved <- lr_results %>%
  dplyr::filter(
    source != "Unresolved",
    target != "Unresolved"
  ) %>%
  dplyr::arrange(
    dplyr::desc(prob)
  )


write.csv(
  lr_resolved,
  "results/tables/CellChat/CellChat_resolved_ligand_receptor_results.csv",
  row.names = FALSE
)


# ============================================================
# 27. Top Resolved Ligand-Receptor Interactions
# ============================================================

if (nrow(lr_resolved) > 0) {
  
  top_lr_resolved <- lr_resolved %>%
    dplyr::select(
      source,
      target,
      ligand,
      receptor,
      prob,
      pval,
      pathway_name
    ) %>%
    dplyr::slice_head(
      n = 40
    )
  
  
  print(
    top_lr_resolved
  )
  
  
  write.csv(
    top_lr_resolved,
    "results/tables/CellChat/CellChat_top40_resolved_ligand_receptor_interactions.csv",
    row.names = FALSE
  )
}


# ============================================================
# 28. Top Ligand-Receptor Interactions
#     for Main Pathways
# ============================================================

for (pathway in main_pathways_present) {
  
  pathway_lr <- lr_resolved %>%
    dplyr::filter(
      pathway_name == pathway
    ) %>%
    dplyr::arrange(
      dplyr::desc(prob)
    )
  
  
  if (nrow(pathway_lr) > 0) {
    
    pathway_lr_top <- pathway_lr %>%
      dplyr::select(
        source,
        target,
        ligand,
        receptor,
        prob,
        pval,
        pathway_name
      ) %>%
      dplyr::slice_head(
        n = 20
      )
    
    
    write.csv(
      pathway_lr_top,
      paste0(
        "results/tables/CellChat/",
        "CellChat_",
        pathway,
        "_top20_LR.csv"
      ),
      row.names = FALSE
    )
    
    
    cat(
      "\n\n====================================\n"
    )
    
    cat(
      "TOP LR INTERACTIONS:",
      pathway,
      "\n"
    )
    
    cat(
      "====================================\n"
    )
    
    print(
      pathway_lr_top
    )
  }
}


# ============================================================
# 29. Secondary Pathway Ligand-Receptor Results
# ============================================================

for (pathway in secondary_pathways_present) {
  
  pathway_lr <- lr_resolved %>%
    dplyr::filter(
      pathway_name == pathway
    ) %>%
    dplyr::arrange(
      dplyr::desc(prob)
    )
  
  
  if (nrow(pathway_lr) > 0) {
    
    pathway_lr_top <- pathway_lr %>%
      dplyr::select(
        source,
        target,
        ligand,
        receptor,
        prob,
        pval,
        pathway_name
      ) %>%
      dplyr::slice_head(
        n = 20
      )
    
    
    write.csv(
      pathway_lr_top,
      paste0(
        "results/tables/CellChat/",
        "CellChat_",
        pathway,
        "_top20_LR.csv"
      ),
      row.names = FALSE
    )
  }
}


# ============================================================
# 30. Selected Pathway Visualizations
# ============================================================

all_selected_pathways <- unique(
  c(
    main_pathways_present,
    secondary_pathways_present
  )
)


for (pathway in all_selected_pathways) {
  
  output_file <- paste0(
    "results/figures/CellChat/CellChat_",
    pathway,
    "_circle.png"
  )
  
  
  png(
    output_file,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  
  
  try(
    netVisual_aggregate(
      cellchat,
      signaling = pathway,
      layout = "circle"
    ),
    silent = TRUE
  )
  
  
  dev.off()
}


# ============================================================
# 31. RESOLVED CELLCHAT NETWORK
# ============================================================
#
# Remove "Unresolved" while preserving the original
# CellChat object and its communication results.
#
# IMPORTANT:
# This is a group-level reconstruction of the communication
# network. The original CellChat object remains unchanged.
# ============================================================


# ------------------------------------------------------------
# 31.1 Define Resolved Cell Types
# ------------------------------------------------------------

old_levels <- levels(
  cellchat@idents
)

remove_type <- "Unresolved"

keep_levels <- setdiff(
  old_levels,
  remove_type
)

keep_pos <- match(
  keep_levels,
  old_levels
)

n_groups_old <- length(
  old_levels
)


stopifnot(
  length(keep_levels) == 14,
  length(keep_pos) == 14
)


cat(
  "Original CellChat cell types:",
  n_groups_old,
  "\n"
)

cat(
  "Resolved CellChat cell types:",
  length(keep_levels),
  "\n"
)

cat(
  "Removed cell type:",
  remove_type,
  "\n"
)

cat(
  "Resolved cell types:\n"
)

print(
  keep_levels
)


# ------------------------------------------------------------
# 31.2 Create Resolved CellChat Object
# ------------------------------------------------------------

cellchat_resolved <- cellchat


# ------------------------------------------------------------
# 31.3 Subset Ligand-Receptor Communication Results
# ------------------------------------------------------------

cellchat_resolved@net$prob <-
  cellchat@net$prob[
    keep_pos,
    keep_pos,
    ,
    drop = FALSE
  ]

cellchat_resolved@net$pval <-
  cellchat@net$pval[
    keep_pos,
    keep_pos,
    ,
    drop = FALSE
  ]

cellchat_resolved@net$count <-
  cellchat@net$count[
    keep_pos,
    keep_pos,
    drop = FALSE
  ]

cellchat_resolved@net$weight <-
  cellchat@net$weight[
    keep_pos,
    keep_pos,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 31.4 Subset Pathway-Level Communication Results
# ------------------------------------------------------------

cellchat_resolved@netP$prob <-
  cellchat@netP$prob[
    keep_pos,
    keep_pos,
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 31.5 Subset Pathway Centrality Information
# ------------------------------------------------------------

cellchat_resolved@netP$centr <-
  lapply(
    cellchat@netP$centr,
    function(path_centr) {
      
      lapply(
        path_centr,
        function(vec) {
          
          if (
            !is.null(names(vec)) &&
            all(keep_levels %in% names(vec))
          ) {
            
            vec[keep_levels]
            
          } else if (
            length(vec) == n_groups_old
          ) {
            
            vec[keep_pos]
            
          } else {
            
            vec
          }
        }
      )
    }
  )


# ------------------------------------------------------------
# 31.6 Update CellChat Identities
# ------------------------------------------------------------

cellchat_resolved@idents <-
  droplevels(
    cellchat@idents[
      cellchat@idents %in% keep_levels
    ]
  )


# ------------------------------------------------------------
# 31.7 Validate Resolved CellChat Object
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "RESOLVED CELLCHAT VALIDATION\n"
)

cat(
  "============================================================\n"
)

cat(
  "Class:",
  class(cellchat_resolved),
  "\n"
)

cat(
  "Number of resolved cell types:",
  length(
    levels(
      cellchat_resolved@idents
    )
  ),
  "\n"
)

cat(
  "Number of pathways:",
  length(
    cellchat_resolved@netP$pathways
  ),
  "\n"
)

cat(
  "LR probability dimensions:",
  paste(
    dim(
      cellchat_resolved@net$prob
    ),
    collapse = " × "
  ),
  "\n"
)

cat(
  "LR p-value dimensions:",
  paste(
    dim(
      cellchat_resolved@net$pval
    ),
    collapse = " × "
  ),
  "\n"
)

cat(
  "Interaction count dimensions:",
  paste(
    dim(
      cellchat_resolved@net$count
    ),
    collapse = " × "
  ),
  "\n"
)

cat(
  "Interaction weight dimensions:",
  paste(
    dim(
      cellchat_resolved@net$weight
    ),
    collapse = " × "
  ),
  "\n"
)

cat(
  "Pathway probability dimensions:",
  paste(
    dim(
      cellchat_resolved@netP$prob
    ),
    collapse = " × "
  ),
  "\n"
)

cat(
  "============================================================\n\n"
)


# ------------------------------------------------------------
# 31.8 Resolved Network — Interaction Count
# ------------------------------------------------------------

png(
  "results/figures/CellChat/CellChat_resolved_count_circle.png",
  width = 8,
  height = 8,
  units = "in",
  res = 300
)

netVisual_circle(
  cellchat_resolved@net$count,
  vertex.weight = as.numeric(
    table(
      cellchat_resolved@idents
    )
  ),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Resolved CellChat interactions — count"
)

dev.off()


# ------------------------------------------------------------
# 31.9 Resolved Network — Interaction Strength
# ------------------------------------------------------------

png(
  "results/figures/CellChat/CellChat_resolved_strength_circle.png",
  width = 8,
  height = 8,
  units = "in",
  res = 300
)

netVisual_circle(
  cellchat_resolved@net$weight,
  vertex.weight = as.numeric(
    table(
      cellchat_resolved@idents
    )
  ),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Resolved CellChat interactions — strength"
)

dev.off()


# ============================================================
# 32. Overall Source Ranking — Resolved Network
# ============================================================

source_ranking <- data.frame(
  Cell_Type = rownames(
    cellchat_resolved@net$weight
  ),
  Outgoing_Communication =
    rowSums(
      cellchat_resolved@net$weight,
      na.rm = TRUE
    ),
  row.names = NULL
) %>%
  dplyr::arrange(
    dplyr::desc(
      Outgoing_Communication
    )
  )

print(
  source_ranking
)

write.csv(
  source_ranking,
  "results/tables/CellChat/CellChat_resolved_source_ranking.csv",
  row.names = FALSE
)


# ============================================================
# 33. Overall Target Ranking — Resolved Network
# ============================================================

target_ranking <- data.frame(
  Cell_Type = colnames(
    cellchat_resolved@net$weight
  ),
  Incoming_Communication =
    colSums(
      cellchat_resolved@net$weight,
      na.rm = TRUE
    ),
  row.names = NULL
) %>%
  dplyr::arrange(
    dplyr::desc(
      Incoming_Communication
    )
  )

print(
  target_ranking
)

write.csv(
  target_ranking,
  "results/tables/CellChat/CellChat_resolved_target_ranking.csv",
  row.names = FALSE
)


# ============================================================
# 34. Resolved Cell-Type Composition
# ============================================================

celltype_counts <- as.data.frame(
  table(
    cellchat_resolved@idents
  )
)

colnames(
  celltype_counts
) <- c(
  "Cell_Type",
  "Cells"
)

celltype_counts <- celltype_counts %>%
  dplyr::arrange(
    dplyr::desc(
      Cells
    )
  )

print(
  celltype_counts
)

write.csv(
  celltype_counts,
  "results/tables/CellChat/CellChat_resolved_cell_type_composition.csv",
  row.names = FALSE
)


# ============================================================
# 35. Final Analysis Summary
# ============================================================

analysis_summary <- data.frame(
  
  metric = c(
    "Initial cells",
    "Cells after QC",
    "Cells after doublet removal",
    "Number of clusters",
    "Number of HVGs",
    "Number of annotated cell types",
    "Number of resolved CellChat cell types",
    "Number of CellChat pathways"
  ),
  
  value = c(
    NA,
    qc_summary$cells,
    singlet_summary$cells,
    length(
      levels(
        pbmc_d1$seurat_clusters
      )
    ),
    length(
      VariableFeatures(
        pbmc_d1
      )
    ),
    length(
      unique(
        pbmc_d1$DEG_annotation
      )
    ),
    length(
      levels(
        cellchat_resolved@idents
      )
    ),
    length(
      cellchat_resolved@netP$pathways
    )
  )
)

print(
  analysis_summary
)

write.csv(
  analysis_summary,
  "results/tables/analysis_summary.csv",
  row.names = FALSE
)


# ============================================================
# 36. Save Main Objects
# ============================================================

saveRDS(
  cellchat,
  "results/cellchat_donor1.rds"
)

saveRDS(
  cellchat_resolved,
  "results/cellchat_donor1_resolved.rds"
)

saveRDS(
  pbmc_d1,
  "results/pbmc_donor1_annotated.rds"
)


# ============================================================
# 37. Save Session Information
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  "results/sessionInfo.txt"
)


# ============================================================
# 38. Final Analysis Message
# ============================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "PBMC DONOR 1 ANALYSIS COMPLETED\n"
)

cat(
  "============================================================\n"
)

cat(
  "Cells after QC:",
  qc_summary$cells,
  "\n"
)

cat(
  "Cells after doublet removal:",
  singlet_summary$cells,
  "\n"
)

cat(
  "Clusters:",
  length(
    levels(
      pbmc_d1$seurat_clusters
    )
  ),
  "\n"
)

cat(
  "Annotated cell types:",
  length(
    unique(
      pbmc_d1$DEG_annotation
    )
  ),
  "\n"
)

cat(
  "Resolved CellChat cell types:",
  length(
    levels(
      cellchat_resolved@idents
    )
  ),
  "\n"
)

cat(
  "CellChat pathways:",
  length(
    cellchat_resolved@netP$pathways
  ),
  "\n"
)

cat(
  "Main pathways analyzed:",
  paste(
    main_pathways_present,
    collapse = ", "
  ),
  "\n"
)

cat(
  "\nNote: platelet population (~15 cells) was excluded from",
  "CellChat statistical inference (min.cells = 20) due to",
  "insufficient cell number for reliable estimation. See",
  "section 18.6 for details.\n"
)

cat(
  "\nAll results saved under:\n"
)

cat(
  "results/\n"
)

cat(
  "============================================================\n\n"
)