#!/usr/bin/env Rscript

# =========================================================
# File: 01_scRNA_CAF_core.R
# Purpose:
#   Curated core script for single-cell CAF analysis in prostate cancer.
#   This script performs single-cell preprocessing, Harmony integration,
#   major-cell-type annotation (SingleR-assisted with marker-based review),
#   fibroblast extraction, CAF reclustering, marker identification,
#   CopyKAT inference, and CAF-related pathway activity analysis based on
#   GSE185344.
#
# Main input files:
#   1. GSE185344_clean_seurat.rds
#   2. pmid_29625050_pathway.txt
#
# Main output files:
#   1. results/01_scRNA_CAF_core/CAF_marker_gene.tsv
#   2. results/01_scRNA_CAF_core/CAF_top10_markers.tsv
#   3. results/01_scRNA_CAF_core/CAF_copykat_prediction.tsv
#   4. results/01_scRNA_CAF_core/CAF_pathway_scores.tsv
#   5. results/01_scRNA_CAF_core/CAF_core_object.rds
#
# Notes:
#   - This is a curated submission-ready script.
#   - Historical intermediate objects, duplicated plotting code, and
#     exploratory auxiliary analyses have been removed.
#   - Harmony / annotation steps are explicitly included to better match the
#     manuscript. If a SingleR reference cannot be loaded, the script will
#     fall back to marker-guided broad annotation and write a warning file.
#   - The script assumes the input RDS already contains a Seurat object built
#     from GSE185344.
# =========================================================

options(stringsAsFactors = FALSE, check.names = FALSE)
set.seed(123)

## -----------------------------
## User-configurable parameters
## -----------------------------
INPUT_RDS <- "GSE185344_clean_seurat.rds"
PATHWAY_FILE <- "pmid_29625050_pathway.txt"
OUTDIR <- file.path("results", "01_scRNA_CAF_core")
N_CORES <- 6

# QC thresholds (matching Supplementary Table S4)
MIN_FEATURES <- 300
MAX_FEATURES <- 8000
MIN_COUNTS <- 1000
MAX_COUNTS <- 8000
MAX_PERCENT_MT <- 10
MAX_PERCENT_RIBO <- 40

# Clustering / reduction
N_HVG <- 2000
N_PCS <- 30
MAIN_RESOLUTION <- 0.2
CAF_RESOLUTION <- 0.3

# Marker calling
MARKER_LOGFC <- 1
MARKER_MIN_PCT <- 0.35
MARKER_PADJ <- 0.05

# CopyKAT settings
RUN_COPYKAT <- TRUE
COPYKAT_ID_TYPE <- "S"
COPYKAT_CELL_LINE <- "no"
COPYKAT_NGENE_CHR <- 5
COPYKAT_WIN_SIZE <- 25
COPYKAT_KS_CUT <- 0.05
COPYKAT_DISTANCE <- "euclidean"
COPYKAT_SAMPLE_NAME <- "PRAD"

# Major-cell-type marker sets used for broad manual review
BROAD_MARKERS <- list(
  Fibroblasts = c("ACTA2", "FAP", "POSTN", "ASPN", "MFAP5", "PDPN", "ITGA11", "PDGFRA", "PDGFRB", "COL11A1"),
  Immune = c("CD68", "CD163", "CD14", "CSF1R", "CD3E", "CD3D", "FGFBP2", "CCR7", "CX3CR1", "CXCR1", "CD19", "MS4A1", "IGHG1", "MZB1", "CD79A", "KIT"),
  Endothelial = c("PECAM1", "CDH5", "VWF", "SELE"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT5", "SOX9")
)

## -----------------------------
## Helpers
## -----------------------------
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")

ensure_dir <- function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE)

check_and_load <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }
  suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
}

write_tsv <- function(x, file) {
  data.table::fwrite(x, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
}

assay_counts <- function(obj, assay = "RNA") {
  if ("layer" %in% names(formals(SeuratObject::GetAssayData))) {
    SeuratObject::GetAssayData(obj, assay = assay, layer = "counts")
  } else {
    SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
  }
}

assay_data <- function(obj, assay = "RNA") {
  if ("layer" %in% names(formals(SeuratObject::GetAssayData))) {
    SeuratObject::GetAssayData(obj, assay = assay, layer = "data")
  } else {
    SeuratObject::GetAssayData(obj, assay = assay, slot = "data")
  }
}

save_pdf <- function(plot_obj, file, width = 8, height = 6) {
  ggplot2::ggsave(filename = file, plot = plot_obj, width = width, height = height)
}

safe_gene_vector <- function(x, universe) unique(x[x %in% universe])

ssgsea_matrix <- function(expr_mat, pathway_df) {
  if (!is.matrix(expr_mat)) expr_mat <- as.matrix(expr_mat)
  pathway_df <- as.data.frame(pathway_df)
  colnames(pathway_df)[1:2] <- c("Gene", "Pathway")
  pathway_df <- pathway_df[!is.na(pathway_df$Gene) & !is.na(pathway_df$Pathway), , drop = FALSE]

  pathway_score <- data.frame()
  for (pw in unique(pathway_df$Pathway)) {
    genes <- safe_gene_vector(pathway_df$Gene[pathway_df$Pathway == pw], rownames(expr_mat))
    if (length(genes) < 2) next
    gs <- GSEABase::GeneSet(
      setName = pw,
      geneIds = genes,
      geneIdType = GSEABase::SymbolIdentifier()
    )
    gsc <- GSEABase::GeneSetCollection(list(gs))
    param <- GSVA::ssgseaParam(
      exprData = expr_mat,
      geneSets = gsc,
      minSize = 2,
      maxSize = 5000,
      normalize = TRUE
    )
    score <- GSVA::gsva(param, verbose = FALSE)
    rownames(score) <- pw
    pathway_score <- rbind(pathway_score, score)
  }
  t(pathway_score)
}

infer_broad_label_from_singleR <- function(labels) {
  x <- tolower(labels)
  out <- rep("Unassigned", length(x))
  out[grepl("fibro|mesench|stromal|smooth", x)] <- "Fibroblasts"
  out[grepl("endo", x)] <- "Endothelial"
  out[grepl("epithelial|keratin|prostate|basal|luminal", x)] <- "Epithelial"
  out[grepl("t cell|b cell|nk|macroph|monocyte|dendritic|mast|immune|plasma|myeloid", x)] <- "Immune"
  out
}

cluster_majority <- function(cluster_vec, label_vec) {
  stopifnot(length(cluster_vec) == length(label_vec))
  tab <- data.frame(cluster = cluster_vec, label = label_vec, stringsAsFactors = FALSE)
  res <- tab |>
    dplyr::count(cluster, label, name = "n") |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  stats::setNames(res$label, res$cluster)
}

compute_cluster_marker_scores <- function(obj, markers) {
  expr <- as.matrix(assay_data(obj, assay = DefaultAssay(obj)))
  meta <- obj@meta.data
  clus <- as.character(Seurat::Idents(obj))
  out <- list()
  for (nm in names(markers)) {
    genes <- safe_gene_vector(markers[[nm]], rownames(expr))
    if (length(genes) == 0) next
    score <- colMeans(t(expr[genes, , drop = FALSE]))
    tmp <- data.frame(cluster = clus, score = score)
    out[[nm]] <- tmp |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(score = mean(score), .groups = "drop") |>
      dplyr::mutate(label = nm)
  }
  dplyr::bind_rows(out)
}

broad_annotation_from_markers <- function(obj, markers) {
  marker_scores <- compute_cluster_marker_scores(obj, markers)
  best <- marker_scores |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(order_by = score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  stats::setNames(best$label, best$cluster)
}

rename_caf_clusters <- function(obj) {
  cur <- levels(obj)
  new <- paste0("CAF_", seq_along(cur) - 1)
  names(new) <- cur
  obj <- Seurat::RenameIdents(obj, new)
  obj$CAF_subcluster <- as.character(Seurat::Idents(obj))
  obj
}

wilcox_by_subcluster <- function(df, pathways) {
  out <- list()
  for (sb in sort(unique(df$CAF_subcluster))) {
    sub_df <- df[df$CAF_subcluster == sb, , drop = FALSE]
    if (length(unique(sub_df$copykat_pred_binary)) < 2) next
    for (pw in pathways) {
      if (!pw %in% colnames(sub_df)) next
      wt <- tryCatch(
        wilcox.test(sub_df[[pw]] ~ sub_df$copykat_pred_binary),
        error = function(e) NULL
      )
      if (is.null(wt)) next
      out[[length(out) + 1]] <- data.frame(
        CAF_subcluster = sb,
        Pathway = pw,
        p_value = wt$p.value,
        median_malignant = median(sub_df[sub_df$copykat_pred_binary == "malignant", pw], na.rm = TRUE),
        median_non_malignant = median(sub_df[sub_df$copykat_pred_binary == "non_malignant", pw], na.rm = TRUE)
      )
    }
  }
  if (length(out) == 0) return(data.frame())
  dplyr::bind_rows(out) |>
    dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"))
}

## -----------------------------
## Main
## -----------------------------
ensure_dir(OUTDIR)
ensure_dir(file.path(OUTDIR, "plots"))

check_and_load(c(
  "Seurat", "SeuratObject", "dplyr", "ggplot2", "patchwork", "data.table",
  "Matrix", "harmony", "SingleR", "SummarizedExperiment",
  "copykat", "GSVA", "GSEABase", "clustree"
))

if (N_CORES > 1 && requireNamespace("doParallel", quietly = TRUE)) {
  doParallel::registerDoParallel(cores = N_CORES)
}

if (!file.exists(INPUT_RDS)) stop("Input RDS not found: ", INPUT_RDS)
if (!file.exists(PATHWAY_FILE)) stop("Pathway file not found: ", PATHWAY_FILE)

msg("Loading Seurat object: ", INPUT_RDS)
sce <- readRDS(INPUT_RDS)
if (!inherits(sce, "Seurat")) stop("Input RDS must contain a Seurat object.")

# Basic metadata / QC ---------------------------------------------------------
msg("Running QC and preprocessing")
sce[["percent.Ribo"]] <- Seurat::PercentageFeatureSet(sce, pattern = "^RP[SL]")

qc_before <- Seurat::VlnPlot(
  sce,
  group.by = "orig.ident",
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.Ribo"),
  pt.size = 0,
  ncol = 4
)
save_pdf(qc_before, file.path(OUTDIR, "plots", "qc_violin_before.pdf"), width = 15, height = 5)

raw_counts <- data.frame(sample = names(table(sce$orig.ident)), n_cells = as.integer(table(sce$orig.ident)))
write_tsv(raw_counts, file.path(OUTDIR, "raw_cell_counts.tsv"))

sce <- subset(
  sce,
  subset = nFeature_RNA > MIN_FEATURES & nFeature_RNA < MAX_FEATURES &
    nCount_RNA > MIN_COUNTS & nCount_RNA < MAX_COUNTS &
    percent.mt < MAX_PERCENT_MT & percent.Ribo < MAX_PERCENT_RIBO
)

qc_after <- Seurat::VlnPlot(
  sce,
  group.by = "orig.ident",
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.Ribo"),
  pt.size = 0,
  ncol = 4
)
save_pdf(qc_after, file.path(OUTDIR, "plots", "qc_violin_after.pdf"), width = 15, height = 5)

clean_counts <- data.frame(sample = names(table(sce$orig.ident)), n_cells = as.integer(table(sce$orig.ident)))
write_tsv(clean_counts, file.path(OUTDIR, "clean_cell_counts.tsv"))

# Normalize / HVG / Harmony / clustering -------------------------------------
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = N_HVG)
sce <- ScaleData(sce, vars.to.regress = c("percent.mt"))
sce <- RunPCA(sce, features = VariableFeatures(sce), npcs = max(50, N_PCS))

pca_plot <- Seurat::DimPlot(sce, reduction = "pca", group.by = "orig.ident")
save_pdf(pca_plot, file.path(OUTDIR, "plots", "pca_by_sample.pdf"), width = 8, height = 6)

elbow_plot <- Seurat::ElbowPlot(sce, ndims = 50, reduction = "pca")
save_pdf(elbow_plot, file.path(OUTDIR, "plots", "elbowplot.pdf"), width = 5, height = 5)

harmony_group <- if ("orig.ident" %in% colnames(sce@meta.data)) "orig.ident" else NULL
if (is.null(harmony_group)) stop("orig.ident metadata is required for Harmony integration.")

sce <- harmony::RunHarmony(sce, group.by.vars = harmony_group, reduction = "pca", assay.use = DefaultAssay(sce), plot_convergence = FALSE)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:N_PCS)
sce <- FindClusters(sce, resolution = seq(0.1, 1.0, by = 0.1))
sce <- RunUMAP(sce, dims = 1:N_PCS, reduction = "harmony")
sce <- RunTSNE(sce, dims = 1:N_PCS, reduction = "harmony")

clustree_plot <- clustree::clustree(sce)
save_pdf(clustree_plot, file.path(OUTDIR, "plots", "clustree_main.pdf"), width = 8, height = 15)

sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:N_PCS)
sce <- FindClusters(sce, resolution = MAIN_RESOLUTION)

main_umap <- (
  Seurat::DimPlot(sce, reduction = "umap", group.by = "orig.ident") + ggplot2::ggtitle("UMAP by sample")
) | (
  Seurat::DimPlot(sce, reduction = "umap", group.by = "type") + ggplot2::ggtitle("UMAP by type")
) | (
  Seurat::DimPlot(sce, reduction = "tsne", label = TRUE) + ggplot2::ggtitle("t-SNE by cluster")
) | (
  Seurat::DimPlot(sce, reduction = "umap", label = TRUE) + ggplot2::ggtitle("UMAP by cluster")
)
save_pdf(main_umap, file.path(OUTDIR, "plots", "main_umap_overview.pdf"), width = 20, height = 12)

# Annotation -----------------------------------------------------------------
msg("Running major-cell-type annotation")
DefaultAssay(sce) <- "RNA"
Idents(sce) <- "seurat_clusters"

marker_dot <- Seurat::DotPlot(sce, features = unlist(BROAD_MARKERS), cols = c("yellow", "red")) +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(marker_dot, file.path(OUTDIR, "plots", "major_celltype_marker_dotplot.pdf"), width = 16, height = 8)

singleR_warning <- NULL
singleR_labels <- rep(NA_character_, ncol(sce))
names(singleR_labels) <- colnames(sce)

if (requireNamespace("celldex", quietly = TRUE)) {
  msg("Trying SingleR annotation using celldex::HumanPrimaryCellAtlasData()")
  ref <- tryCatch(celldex::HumanPrimaryCellAtlasData(), error = function(e) NULL)
  if (!is.null(ref)) {
    singler_res <- SingleR::SingleR(
      test = SummarizedExperiment::SummarizedExperiment(list(logcounts = as.matrix(assay_data(sce, assay = DefaultAssay(sce))))),
      ref = ref,
      labels = ref$label.main
    )
    singleR_labels <- as.character(singler_res$labels)
    names(singleR_labels) <- rownames(singler_res)
  } else {
    singleR_warning <- "SingleR reference could not be loaded; falling back to marker-guided broad annotation only."
  }
} else {
  singleR_warning <- "Package 'celldex' not available; falling back to marker-guided broad annotation only."
}

sce$SingleR_label <- singleR_labels[colnames(sce)]
cluster_singler <- cluster_majority(as.character(Idents(sce)), infer_broad_label_from_singleR(sce$SingleR_label))
cluster_marker_based <- broad_annotation_from_markers(sce, BROAD_MARKERS)

all_clusters <- sort(unique(as.character(Idents(sce))))
cluster_annotation_tbl <- data.frame(
  cluster = all_clusters,
  broad_label_singleR = unname(cluster_singler[all_clusters]),
  broad_label_marker = unname(cluster_marker_based[all_clusters]),
  stringsAsFactors = FALSE
)
cluster_annotation_tbl$final_broad_label <- ifelse(
  !is.na(cluster_annotation_tbl$broad_label_singleR) & cluster_annotation_tbl$broad_label_singleR != "Unassigned",
  cluster_annotation_tbl$broad_label_singleR,
  cluster_annotation_tbl$broad_label_marker
)
cluster_annotation_tbl$final_broad_label[is.na(cluster_annotation_tbl$final_broad_label)] <- "Unassigned"
write_tsv(cluster_annotation_tbl, file.path(OUTDIR, "cluster_annotation_review.tsv"))

cluster_to_label <- stats::setNames(cluster_annotation_tbl$final_broad_label, cluster_annotation_tbl$cluster)
sce$celltype <- unname(cluster_to_label[as.character(Idents(sce))])
sce$celltype <- factor(sce$celltype, levels = c("Fibroblasts", "Endothelial", "Epithelial", "Immune", "Unassigned"))

annot_colors <- c(Fibroblasts = "#1F77B4", Endothelial = "#FF7F0E", Epithelial = "lightgreen", Immune = "#D62728", Unassigned = "grey70")
annot_umap <- Seurat::DimPlot(
  sce,
  group.by = "celltype",
  cols = annot_colors,
  shuffle = TRUE
) + ggplot2::ggtitle("Cell type annotation")
save_pdf(annot_umap, file.path(OUTDIR, "plots", "celltype_annotation_umap.pdf"), width = 8, height = 6)

celltype_table <- data.frame(
  cell = colnames(sce),
  sample = sce$orig.ident,
  seurat_cluster = as.character(Idents(sce)),
  celltype = as.character(sce$celltype),
  singleR_label = sce$SingleR_label,
  stringsAsFactors = FALSE
)
write_tsv(celltype_table, file.path(OUTDIR, "celltype_annotation_table.tsv"))

if (!is.null(singleR_warning)) {
  writeLines(singleR_warning, con = file.path(OUTDIR, "annotation_warning.txt"))
}

# Extract fibroblasts ---------------------------------------------------------
fibro_clusters <- cluster_annotation_tbl$cluster[cluster_annotation_tbl$final_broad_label == "Fibroblasts"]
if (length(fibro_clusters) == 0) {
  stop("No fibroblast-related clusters were identified. Review cluster_annotation_review.tsv.")
}

msg("Extracting fibroblast-related clusters: ", paste(fibro_clusters, collapse = ", "))
CAF <- subset(sce, idents = fibro_clusters)
saveRDS(CAF, file = file.path(OUTDIR, "fibroblast_subset_raw.rds"))

# CAF reclustering ------------------------------------------------------------
msg("Running CAF reclustering")
DefaultAssay(CAF) <- "RNA"
CAF <- NormalizeData(CAF, normalization.method = "LogNormalize", scale.factor = 10000)
CAF <- FindVariableFeatures(CAF, selection.method = "vst", nfeatures = N_HVG)
CAF <- ScaleData(CAF, vars.to.regress = c("percent.mt"))
CAF <- RunPCA(CAF, features = VariableFeatures(CAF), npcs = max(50, N_PCS))
CAF <- FindNeighbors(CAF, reduction = "pca", dims = 1:N_PCS)
CAF <- FindClusters(CAF, resolution = seq(0.1, 1.0, by = 0.1))
CAF <- RunUMAP(CAF, dims = 1:N_PCS, reduction = "pca")
CAF <- RunTSNE(CAF, dims = 1:N_PCS, reduction = "pca")

clustree_caf <- clustree::clustree(CAF)
save_pdf(clustree_caf, file.path(OUTDIR, "plots", "clustree_caf.pdf"), width = 8, height = 15)

CAF <- FindNeighbors(CAF, reduction = "pca", dims = 1:N_PCS)
CAF <- FindClusters(CAF, resolution = CAF_RESOLUTION)
CAF <- rename_caf_clusters(CAF)
Idents(CAF) <- "CAF_subcluster"

caf_umap <- Seurat::DimPlot(CAF, reduction = "umap", label = TRUE) + ggplot2::ggtitle("Reclustered CAF subclusters")
save_pdf(caf_umap, file.path(OUTDIR, "plots", "caf_umap.pdf"), width = 8, height = 6)

# Marker calling --------------------------------------------------------------
msg("Calling CAF markers")
CAF.markers <- FindAllMarkers(
  object = CAF,
  logfc.threshold = MARKER_LOGFC,
  min.pct = MARKER_MIN_PCT,
  only.pos = TRUE
)
CAF.markers$pct.diff <- CAF.markers$pct.1 - CAF.markers$pct.2
CAF.markers <- CAF.markers[CAF.markers$p_val_adj < MARKER_PADJ, , drop = FALSE]
CAF.markers <- CAF.markers |>
  dplyr::mutate(cluster = as.character(cluster)) |>
  dplyr::arrange(cluster, dplyr::desc(avg_log2FC), p_val_adj)
write_tsv(CAF.markers, file.path(OUTDIR, "CAF_marker_gene.tsv"))

CAF_top10 <- CAF.markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) |>
  dplyr::ungroup()
write_tsv(CAF_top10, file.path(OUTDIR, "CAF_top10_markers.tsv"))

Top5 <- CAF.markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) |>
  dplyr::ungroup()
Top5_genes <- unique(Top5$gene[Top5$gene %in% rownames(CAF)])
marker_dot_caf <- Seurat::DotPlot(CAF, features = Top5_genes, cols = c("blue", "red"), scale = TRUE) +
  Seurat::RotatedAxis() +
  ggplot2::ggtitle("Top 5 marker genes for each CAF subcluster") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
save_pdf(marker_dot_caf, file.path(OUTDIR, "plots", "caf_marker_dotplot.pdf"), width = 10, height = 7)

sample_clust <- as.matrix(table(CAF$orig.ident, CAF$CAF_subcluster))
sample_clust <- apply(sample_clust, 1, function(x) x / sum(x))
sample_clust <- reshape2::melt(sample_clust)
colnames(sample_clust) <- c("CAF_subcluster", "sample", "proportion")
write_tsv(sample_clust, file.path(OUTDIR, "CAF_cluster_sample_fraction.tsv"))

# CopyKAT ---------------------------------------------------------------------
msg("Running CopyKAT")
copykat_table <- NULL
if (RUN_COPYKAT) {
  matrix_counts <- as.matrix(assay_counts(CAF, assay = "RNA"))
  copykat_res <- copykat::copykat(
    rawmat = matrix_counts,
    id.type = COPYKAT_ID_TYPE,
    cell.line = COPYKAT_CELL_LINE,
    ngene.chr = COPYKAT_NGENE_CHR,
    win.size = COPYKAT_WIN_SIZE,
    KS.cut = COPYKAT_KS_CUT,
    sam.name = COPYKAT_SAMPLE_NAME,
    distance = COPYKAT_DISTANCE
  )

  pred <- copykat_res$prediction
  if (is.null(pred)) stop("CopyKAT did not return a prediction table.")
  pred <- as.data.frame(pred)
  if (!"copykat.pred" %in% colnames(pred)) {
    pred_col <- grep("pred", colnames(pred), value = TRUE)[1]
    if (!is.na(pred_col)) colnames(pred)[colnames(pred) == pred_col] <- "copykat.pred"
  }
  if (!"cell.names" %in% colnames(pred)) {
    cell_col <- grep("cell", colnames(pred), value = TRUE)[1]
    if (!is.na(cell_col)) colnames(pred)[colnames(pred) == cell_col] <- "cell.names"
  }
  if (!all(c("cell.names", "copykat.pred") %in% colnames(pred))) {
    stop("CopyKAT prediction table does not contain identifiable cell-name and prediction columns.")
  }
  pred$cell.names <- as.character(pred$cell.names)
  pred <- pred[match(colnames(CAF), pred$cell.names), , drop = FALSE]
  rownames(pred) <- pred$cell.names
  copykat_table <- data.frame(
    cell = pred$cell.names,
    copykat.pred = pred$copykat.pred,
    stringsAsFactors = FALSE
  )
  copykat_table$copykat.pred[is.na(copykat_table$copykat.pred)] <- "Unknown"
  CAF <- AddMetaData(CAF, metadata = setNames(copykat_table$copykat.pred, copykat_table$cell), col.name = "copykat.pred")
} else {
  CAF$copykat.pred <- "Unknown"
  copykat_table <- data.frame(cell = colnames(CAF), copykat.pred = CAF$copykat.pred, stringsAsFactors = FALSE)
}

CAF$copykat_pred_binary <- ifelse(CAF$copykat.pred %in% c("aneuploid", "malignant"), "malignant", ifelse(CAF$copykat.pred %in% c("diploid", "no_malignant", "non-malignant", "non_malignant"), "non_malignant", "Unknown"))
copykat_table$copykat_pred_binary <- CAF$copykat_pred_binary[match(copykat_table$cell, colnames(CAF))]
copykat_table$CAF_subcluster <- CAF$CAF_subcluster[match(copykat_table$cell, colnames(CAF))]
write_tsv(copykat_table, file.path(OUTDIR, "CAF_copykat_prediction.tsv"))

copykat_umap <- Seurat::DimPlot(CAF, reduction = "umap", group.by = "copykat_pred_binary", cols = c(malignant = "red", non_malignant = "blue", Unknown = "grey70")) +
  ggplot2::ggtitle("CopyKAT-derived malignant-like inference")
save_pdf(copykat_umap, file.path(OUTDIR, "plots", "caf_copykat_umap.pdf"), width = 8, height = 6)

clust_malig <- as.data.frame.matrix(table(CAF$copykat_pred_binary, CAF$CAF_subcluster))
clust_malig$copykat_pred_binary <- rownames(clust_malig)
write_tsv(clust_malig, file.path(OUTDIR, "CAF_cluster_malignant_fraction.tsv"))

# Pathway activity ------------------------------------------------------------
msg("Scoring tumor-related pathways")
pathway_df <- data.table::fread(PATHWAY_FILE, data.table = FALSE)
if (ncol(pathway_df) < 2) stop("Pathway file must contain at least two columns.")
pathway_df <- pathway_df[, 1:2, drop = FALSE]
colnames(pathway_df) <- c("Gene", "Pathway")
pathway_df <- pathway_df[pathway_df$Pathway %in% c("CellCyle", "HIPPO", "MYC", "NOTCH", "NRF1", "PI3K", "TGF.Beta", "RAS", "TP53", "WNT") | TRUE, , drop = FALSE]

matrix_counts <- as.matrix(assay_counts(CAF, assay = "RNA"))
pathway_scores <- ssgsea_matrix(matrix_counts, pathway_df)
pathway_scores <- as.data.frame(pathway_scores)
pathway_scores$cell <- rownames(pathway_scores)

pathway_scores_group <- dplyr::left_join(
  data.frame(
    cell = colnames(CAF),
    sample = CAF$orig.ident,
    CAF_subcluster = CAF$CAF_subcluster,
    copykat.pred = CAF$copykat.pred,
    copykat_pred_binary = CAF$copykat_pred_binary,
    stringsAsFactors = FALSE
  ),
  pathway_scores,
  by = "cell"
)
write_tsv(pathway_scores_group, file.path(OUTDIR, "CAF_pathway_scores.tsv"))

pw_cols <- intersect(unique(pathway_df$Pathway), colnames(pathway_scores_group))
if ("TGF.Beta" %in% pw_cols) {
  colnames(pathway_scores_group)[match("TGF.Beta", colnames(pathway_scores_group))] <- "TGF-Beta"
  pw_cols[pw_cols == "TGF.Beta"] <- "TGF-Beta"
}

anno_col <- pathway_scores_group[, c("CAF_subcluster", "copykat_pred_binary"), drop = FALSE]
rownames(anno_col) <- pathway_scores_group$cell
mat <- as.matrix(pathway_scores_group[, pw_cols, drop = FALSE])
rownames(mat) <- pathway_scores_group$cell
ord <- order(anno_col$copykat_pred_binary, anno_col$CAF_subcluster)

pdf(file.path(OUTDIR, "plots", "caf_pathway_heatmap.pdf"), width = 12, height = 9)
pheatmap::pheatmap(
  t(mat[ord, , drop = FALSE]),
  scale = "row",
  show_colnames = FALSE,
  annotation_col = anno_col[ord, , drop = FALSE],
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  annotation_colors = list(copykat_pred_binary = c(malignant = "red", non_malignant = "blue", Unknown = "grey70"))
)
dev.off()

pathway_stats <- wilcox_by_subcluster(pathway_scores_group, pw_cols)
write_tsv(pathway_stats, file.path(OUTDIR, "CAF_pathway_by_copykat_within_subcluster.tsv"))

# Final save ------------------------------------------------------------------
msg("Saving final object")
saveRDS(CAF, file = file.path(OUTDIR, "CAF_core_object.rds"))

session_info <- capture.output(sessionInfo())
writeLines(session_info, con = file.path(OUTDIR, "sessionInfo.txt"))

msg("Done. Key outputs written to: ", OUTDIR)
