# =========================================================
# File: 02_CAF_score_candidate_genes.R
# Purpose:
#   Curated, submission-ready script for TCGA CAF subgroup-score analysis
#   and CAF-related candidate-gene derivation in prostate cancer.
#
# What this script does:
#   1. Reads CAF marker genes derived from the single-cell workflow.
#   2. Calculates TCGA CAF subgroup scores by ssGSEA.
#   3. Compares CAF subgroup scores between normal and tumor tissues.
#   4. Performs optimal-cutoff Kaplan-Meier analysis for each CAF subgroup
#      in the TCGA tumor cohort.
#   5. Identifies prognosis-related CAF states (default: log-rank P < 0.05).
#   6. Integrates marker genes from prognosis-related CAF states and exports
#      CAF_candidate_gene_list.txt for downstream modeling.
#
# Expected manuscript-consistent outputs:
#   - tcga_caf_score_matrix.tsv
#   - tcga_caf_score_with_group.tsv
#   - normal_tumor_comparison.tsv
#   - tcga_caf_score_violin.pdf / .png
#   - tcga_caf_cli.tsv
#   - caf_km_statistics.tsv
#   - tcga_caf_km.pdf / .png
#   - candidate_gene_membership.tsv
#   - CAF_candidate_gene_list.txt
#
# Notes:
#   - This script is intended for the curated GitHub submission package.
#   - Historical local paths, intermediate .RData files, and unrelated
#     exploratory clinical extension plots have been removed.
#   - The primary endpoint is handled as bRFS/BCR; legacy OS/OS.time column
#     names are supported only as input aliases.
# =========================================================

options(stringsAsFactors = FALSE, check.names = FALSE, warn = 1)

CFG <- list(
  outdir = file.path("results", "02_CAF_score_candidate_genes"),
  p_cutoff = 0.05,
  gsva_min_size = 2,
  gsva_max_size = 5000,
  km_palette = c("Low" = "#2F7FC1FF", "High" = "#D94738FF"),
  violin_palette = c("Normal" = "#3C5488FF", "Tumor" = "#E64B35FF"),
  output_candidate_file = "CAF_candidate_gene_list.txt",
  force_selected_states = NULL,   # e.g. c("CAF_0","CAF_1","CAF_3"); keep NULL to use P < 0.05
  save_png_dpi = 300
)

dir.create(CFG$outdir, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

check_and_load_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing_pkgs, collapse = ", "),
      "\nPlease install them before running this script."
    )
  }
  invisible(lapply(pkgs, function(x) suppressPackageStartupMessages(library(x, character.only = TRUE))))
}

check_and_load_packages(c(
  "data.table", "dplyr", "ggplot2", "ggpubr", "survival", "survminer",
  "GSVA", "GSEABase", "patchwork"
))

find_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_delim_auto <- function(file) {
  ext <- tolower(tools::file_ext(file))
  if (ext == "csv") {
    df <- data.table::fread(file, data.table = FALSE, check.names = FALSE)
  } else {
    df <- data.table::fread(file, data.table = FALSE, check.names = FALSE, sep = "\t")
  }
  df
}

average_duplicate_rows <- function(mat) {
  mat <- as.matrix(mat)
  if (is.null(rownames(mat))) stop("Expression matrix must have rownames.")
  rn <- rownames(mat)
  if (!anyDuplicated(rn)) return(mat)
  uniq <- unique(rn)
  out <- sapply(uniq, function(g) {
    idx <- which(rn == g)
    if (length(idx) == 1) {
      mat[idx, ]
    } else {
      colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
    }
  })
  out <- t(out)
  rownames(out) <- uniq
  out
}

to_numeric_matrix <- function(df_or_mat) {
  x <- as.matrix(df_or_mat)
  storage.mode(x) <- "numeric"
  x
}

normalize_cluster_label <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- ifelse(grepl("^CAF_", x), x, paste0("CAF_", x))
  x
}

safe_sample_col <- function(df) {
  cand <- c("Samples", "Sample", "sample", "ID", "Id", "id")
  hit <- cand[cand %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  colnames(df)[1]
}

safe_group_col <- function(df) {
  cand <- c("Type", "type", "Group", "group", "SampleType", "sample_type")
  hit <- cand[cand %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  stop("Could not identify the sample-type/group column in the TCGA group file.")
}

classify_sample_type <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[y %in% c("t", "tumor", "primary tumor", "tumour", "cancer")] <- "Tumor"
  out[y %in% c("n", "normal", "benign", "solid tissue normal", "adjacent normal")] <- "Normal"
  out[grepl("tumor|tumour|cancer", y)] <- "Tumor"
  out[grepl("normal|benign", y)] <- "Normal"
  out
}

standardize_tcga_prefix12 <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.", "-", x)
  ifelse(grepl("^TCGA-", x), substr(x, 1, 12), x)
}

merge_by_samples <- function(left, right, left_sample_col = "Samples", right_sample_col = "Samples", mode = c("direct", "auto")) {
  mode <- match.arg(mode)
  out <- merge(left, right, by.x = left_sample_col, by.y = right_sample_col)
  if (mode == "direct" || nrow(out) > 0) return(out)

  left2 <- left
  right2 <- right
  left2$.__merge_id__ <- standardize_tcga_prefix12(left2[[left_sample_col]])
  right2$.__merge_id__ <- standardize_tcga_prefix12(right2[[right_sample_col]])
  out2 <- merge(left2, right2, by = ".__merge_id__")
  out2
}

guess_time_status_cols <- function(df) {
  cn <- colnames(df)

  time_cands <- c("BCR.time", "bRFS.time", "bRFS_time", "BCR_time", "futime", "OS.time", "time", "Time")
  status_cands <- c("BCR", "bRFS", "BCR.status", "BCR_status", "status", "Status", "fustat", "OS")

  time_hit <- time_cands[time_cands %in% cn]
  status_hit <- status_cands[status_cands %in% cn]

  if (length(time_hit) == 0 || length(status_hit) == 0) {
    stop(
      "Could not identify time/status columns in the TCGA clinical file.\n",
      "Supported time aliases: ", paste(time_cands, collapse = ", "), "\n",
      "Supported status aliases: ", paste(status_cands, collapse = ", ")
    )
  }

  list(time = time_hit[1], status = status_hit[1])
}

build_gene_set_list <- function(marker_df) {
  stopifnot(all(c("gene", "cluster") %in% colnames(marker_df)))
  marker_df <- marker_df[, c("gene", "cluster"), drop = FALSE]
  marker_df$gene <- trimws(as.character(marker_df$gene))
  marker_df$cluster <- normalize_cluster_label(marker_df$cluster)
  marker_df <- marker_df[!is.na(marker_df$gene) & nzchar(marker_df$gene), , drop = FALSE]
  split(unique(marker_df$gene), marker_df$cluster)
}

run_ssgsea <- function(expr_mat, gene_sets, min_size = 2, max_size = 5000) {
  gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= min_size]
  if (length(gene_sets) == 0) stop("No gene set remains after size filtering.")

  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    gs_list <- lapply(names(gene_sets), function(nm) {
      GSEABase::GeneSet(
        setName = nm,
        geneIds = unique(gene_sets[[nm]][gene_sets[[nm]] %in% rownames(expr_mat)]),
        geneIdType = GSEABase::SymbolIdentifier()
      )
    })
    gs_coll <- GSEABase::GeneSetCollection(gs_list)
    param <- GSVA::ssgseaParam(
      exprData = expr_mat,
      geneSets = gs_coll,
      minSize = min_size,
      maxSize = max_size,
      normalize = TRUE
    )
    score <- GSVA::gsva(param, verbose = FALSE)
  } else {
    score <- GSVA::gsva(
      expr = expr_mat,
      gset.idx.list = gene_sets,
      method = "ssgsea",
      kcdf = "Gaussian",
      abs.ranking = TRUE,
      min.sz = min_size,
      max.sz = max_size,
      verbose = FALSE
    )
  }
  score
}

sig_violin <- function(df, group_col, value_col, ylab_text, palette) {
  dat <- df[, c(group_col, value_col), drop = FALSE]
  colnames(dat) <- c("Group", "Score")
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  dat$Group <- factor(dat$Group, levels = names(palette))
  dat <- dat[!is.na(dat$Group), , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)

  comps <- if (length(unique(dat$Group)) >= 2) {
    list(c(levels(dat$Group)[1], levels(dat$Group)[2]))
  } else {
    NULL
  }

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = Group, y = Score, color = Group, fill = Group)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.15) +
    ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA, color = "black") +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = NULL, y = ylab_text) +
    ggplot2::theme_classic(base_size = 12)

  if (!is.null(comps)) {
    p <- p + ggpubr::stat_compare_means(
      comparisons = comps,
      method = "wilcox.test",
      label = "p.signif"
    )
  }
  p
}

fit_km_one <- function(df, score_col, time_col, status_col, state_name, palette) {
  dat <- df[, c(score_col, time_col, status_col), drop = FALSE]
  colnames(dat) <- c("Score", "time", "status")
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  dat$status <- as.numeric(as.character(dat$status))
  dat$time <- as.numeric(as.character(dat$time))

  if (nrow(dat) < 10 || length(unique(dat$status)) < 2 || sd(dat$Score, na.rm = TRUE) == 0) {
    return(list(
      plot = NULL,
      stat = data.frame(
        CAF_state = state_name,
        cutoff = NA_real_,
        pvalue = NA_real_,
        n_low = NA_integer_,
        n_high = NA_integer_,
        selected = FALSE,
        stringsAsFactors = FALSE
      )
    ))
  }

  cut_res <- survminer::surv_cutpoint(
    dat,
    time = "time",
    event = "status",
    variables = "Score"
  )
  cutoff <- as.numeric(cut_res$cutpoint[1, 1])
  dat$RiskGroup <- ifelse(dat$Score >= cutoff, "High", "Low")
  dat$RiskGroup <- factor(dat$RiskGroup, levels = c("Low", "High"))

  surv_obj <- survival::survfit(survival::Surv(time, status) ~ RiskGroup, data = dat)
  pval <- tryCatch({
    sdiff <- survival::survdiff(survival::Surv(time, status) ~ RiskGroup, data = dat)
    1 - stats::pchisq(sdiff$chisq, df = length(sdiff$n) - 1)
  }, error = function(e) NA_real_)

  g <- survminer::ggsurvplot(
    surv_obj,
    data = dat,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = FALSE,
    palette = unname(palette[c("Low", "High")]),
    legend.title = state_name,
    legend.labs = c("Low", "High"),
    xlab = "Time",
    ylab = "bRFS probability"
  )

  stat <- data.frame(
    CAF_state = state_name,
    cutoff = cutoff,
    pvalue = pval,
    n_low = sum(dat$RiskGroup == "Low"),
    n_high = sum(dat$RiskGroup == "High"),
    selected = isTRUE(!is.na(pval) && pval < CFG$p_cutoff),
    stringsAsFactors = FALSE
  )

  list(plot = g$plot + ggplot2::theme_classic(base_size = 12), stat = stat)
}

save_plot_dual <- function(plot_obj, file_stub, width = 7, height = 5) {
  ggplot2::ggsave(
    filename = paste0(file_stub, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    filename = paste0(file_stub, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = CFG$save_png_dpi
  )
}

# -----------------------------
# 1. Locate input files
# -----------------------------
msg("Locating input files ...")

candidate_marker_files <- c(
  "CAF_marker_gene.tsv",
  "scRNA_marker_gene.tsv",
  file.path("results", "scRNA_marker_gene.txt"),
  file.path("results", "03_caf", "scRNA_marker_gene.tsv"),
  file.path("results", "03_caf", "CAF_marker_gene.tsv"),
  file.path("01.scRNA", "results", "scRNA_marker_gene.txt"),
  file.path("..", "01.scRNA", "results", "scRNA_marker_gene.txt")
)

candidate_tcga_expr_files <- c(
  "TCGA_expression_tumor_normal.txt",
  "tcga_dat.txt",
  "TCGA_expression.txt",
  file.path("results", "tcga_dat.txt"),
  file.path("03.data.pre", "TCGA", "results", "tcga_dat.txt"),
  file.path("..", "03.data.pre", "TCGA", "results", "tcga_dat.txt")
)

candidate_tcga_group_files <- c(
  "TCGA_group_info.txt",
  "tcga.group.txt",
  file.path("results", "tcga.group.txt"),
  file.path("03.data.pre", "TCGA", "results", "tcga.group.txt"),
  file.path("..", "03.data.pre", "TCGA", "results", "tcga.group.txt")
)

candidate_tcga_cli_files <- c(
  "TCGA_bcr_clinical.txt",
  "tcga_cli.txt",
  file.path("results", "tcga_cli.txt"),
  file.path("03.data.pre", "TCGA", "results", "tcga_cli.txt"),
  file.path("..", "03.data.pre", "TCGA", "results", "tcga_cli.txt")
)

marker_file <- find_first_existing(candidate_marker_files)
expr_file <- find_first_existing(candidate_tcga_expr_files)
group_file <- find_first_existing(candidate_tcga_group_files)
cli_file <- find_first_existing(candidate_tcga_cli_files)

if (is.na(marker_file)) stop("Could not find CAF marker file.")
if (is.na(expr_file)) stop("Could not find TCGA expression file.")
if (is.na(group_file)) stop("Could not find TCGA group file.")
if (is.na(cli_file)) stop("Could not find TCGA clinical file.")

msg("Using marker file: ", marker_file)
msg("Using TCGA expression file: ", expr_file)
msg("Using TCGA group file: ", group_file)
msg("Using TCGA clinical file: ", cli_file)

# -----------------------------
# 2. Read and standardize inputs
# -----------------------------
msg("Reading input data ...")

marker_df <- read_delim_auto(marker_file)
expr_df <- read_delim_auto(expr_file)
group_df <- read_delim_auto(group_file)
cli_df <- read_delim_auto(cli_file)

# Marker file
if (!all(c("gene", "cluster") %in% colnames(marker_df))) {
  stop("CAF marker file must contain at least two columns: gene and cluster.")
}
marker_df <- marker_df[, colnames(marker_df) %in% c("gene", "cluster", "avg_log2FC", "p_val_adj", "pct.1", "pct.2"), drop = FALSE]
marker_df$cluster <- normalize_cluster_label(marker_df$cluster)

# Expression file
expr_gene_col <- colnames(expr_df)[1]
rownames(expr_df) <- expr_df[[expr_gene_col]]
expr_df[[expr_gene_col]] <- NULL
expr_mat <- to_numeric_matrix(expr_df)
expr_mat <- average_duplicate_rows(expr_mat)
expr_mat <- expr_mat[rowMeans(expr_mat, na.rm = TRUE) > 0, , drop = FALSE]

# Group file
group_sample_col <- safe_sample_col(group_df)
group_type_col <- safe_group_col(group_df)
group_df$Samples <- as.character(group_df[[group_sample_col]])
group_df$SampleType <- classify_sample_type(group_df[[group_type_col]])
if (all(is.na(group_df$SampleType))) {
  stop("Could not classify tumor/normal labels from the TCGA group file.")
}
group_df <- group_df[, c("Samples", setdiff(colnames(group_df), c("Samples"))), drop = FALSE]

# Clinical file
cli_sample_col <- safe_sample_col(cli_df)
cli_df$Samples <- as.character(cli_df[[cli_sample_col]])
ts_cols <- guess_time_status_cols(cli_df)
cli_df$time_bRFS <- as.numeric(as.character(cli_df[[ts_cols$time]]))
cli_df$status_bRFS <- as.numeric(as.character(cli_df[[ts_cols$status]]))
cli_df <- cli_df[, c("Samples", setdiff(colnames(cli_df), "Samples")), drop = FALSE]

# -----------------------------
# 3. Calculate CAF subgroup scores by ssGSEA
# -----------------------------
msg("Calculating TCGA CAF subgroup scores by ssGSEA ...")

gene_sets <- build_gene_set_list(marker_df)
msg("Detected CAF states: ", paste(names(gene_sets), collapse = ", "))

score_mat <- run_ssgsea(
  expr_mat = expr_mat,
  gene_sets = gene_sets,
  min_size = CFG$gsva_min_size,
  max_size = CFG$gsva_max_size
)
score_df <- data.frame(Samples = colnames(score_mat), t(score_mat), check.names = FALSE)

write.table(
  score_df,
  file = file.path(CFG$outdir, "tcga_caf_score_matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------
# 4. Normal vs tumor comparison
# -----------------------------
msg("Running normal-vs-tumor comparison ...")

score_group <- merge_by_samples(group_df, score_df, left_sample_col = "Samples", right_sample_col = "Samples", mode = "auto")
if (!"SampleType" %in% colnames(score_group)) {
  stop("SampleType column missing after merging group and score tables.")
}

write.table(
  score_group,
  file = file.path(CFG$outdir, "tcga_caf_score_with_group.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

caf_states <- names(gene_sets)
nt_stats <- lapply(caf_states, function(state) {
  dat <- score_group[, c("SampleType", state), drop = FALSE]
  colnames(dat) <- c("SampleType", "Score")
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  dat <- dat[dat$SampleType %in% c("Normal", "Tumor"), , drop = FALSE]
  pval <- if (length(unique(dat$SampleType)) == 2) {
    tryCatch(stats::wilcox.test(Score ~ SampleType, data = dat)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }
  data.frame(
    CAF_state = state,
    n_normal = sum(dat$SampleType == "Normal"),
    n_tumor = sum(dat$SampleType == "Tumor"),
    median_normal = stats::median(dat$Score[dat$SampleType == "Normal"], na.rm = TRUE),
    median_tumor = stats::median(dat$Score[dat$SampleType == "Tumor"], na.rm = TRUE),
    pvalue = pval,
    stringsAsFactors = FALSE
  )
})
nt_stats <- do.call(rbind, nt_stats)
write.table(
  nt_stats,
  file = file.path(CFG$outdir, "normal_tumor_comparison.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

violin_list <- lapply(caf_states, function(state) {
  sig_violin(
    df = score_group[score_group$SampleType %in% c("Normal", "Tumor"), , drop = FALSE],
    group_col = "SampleType",
    value_col = state,
    ylab_text = paste0(state, " score"),
    palette = CFG$violin_palette
  )
})
violin_list <- violin_list[!vapply(violin_list, is.null, logical(1))]
if (length(violin_list) > 0) {
  violin_plot <- patchwork::wrap_plots(violin_list, nrow = 1) + patchwork::plot_annotation(tag_levels = "A")
  save_plot_dual(violin_plot, file.path(CFG$outdir, "tcga_caf_score_violin"), width = 15, height = 4.8)
}

# -----------------------------
# 5. Kaplan-Meier analysis in TCGA tumor cohort
# -----------------------------
msg("Running Kaplan-Meier analysis for CAF subgroup scores ...")

# Keep tumor samples only for survival analysis
score_cli <- merge_by_samples(cli_df, score_df, left_sample_col = "Samples", right_sample_col = "Samples", mode = "auto")
if (!all(c("time_bRFS", "status_bRFS") %in% colnames(score_cli))) {
  stop("Clinical merge failed: missing time_bRFS or status_bRFS.")
}

# If group labels are available, keep tumor samples only
score_cli2 <- merge_by_samples(score_cli, group_df[, c("Samples", "SampleType"), drop = FALSE],
                               left_sample_col = "Samples", right_sample_col = "Samples", mode = "auto")
if ("SampleType" %in% colnames(score_cli2) && any(!is.na(score_cli2$SampleType))) {
  tumor_cli <- score_cli2[score_cli2$SampleType == "Tumor", , drop = FALSE]
} else {
  tumor_cli <- score_cli2
  msg("Warning: no usable SampleType after clinical merge; proceeding without explicit tumor-only filter.")
}

tumor_cli <- tumor_cli[tumor_cli$time_bRFS > 30, , drop = FALSE]
write.table(
  tumor_cli,
  file = file.path(CFG$outdir, "tcga_caf_cli.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

km_results <- lapply(caf_states, function(state) {
  fit_km_one(
    df = tumor_cli,
    score_col = state,
    time_col = "time_bRFS",
    status_col = "status_bRFS",
    state_name = state,
    palette = CFG$km_palette
  )
})

km_stats <- do.call(rbind, lapply(km_results, `[[`, "stat"))
write.table(
  km_stats,
  file = file.path(CFG$outdir, "caf_km_statistics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

km_plots <- lapply(km_results, `[[`, "plot")
names(km_plots) <- caf_states
km_plots_nonnull <- km_plots[!vapply(km_plots, is.null, logical(1))]
if (length(km_plots_nonnull) > 0) {
  km_combined <- patchwork::wrap_plots(km_plots_nonnull, nrow = 1) + patchwork::plot_annotation(tag_levels = "A")
  save_plot_dual(km_combined, file.path(CFG$outdir, "tcga_caf_km"), width = 16, height = 4.8)
}

# -----------------------------
# 6. Derive CAF candidate gene list
# -----------------------------
msg("Deriving CAF candidate gene list ...")

if (!is.null(CFG$force_selected_states)) {
  selected_states <- intersect(normalize_cluster_label(CFG$force_selected_states), caf_states)
  if (length(selected_states) == 0) {
    stop("CFG$force_selected_states was set, but none matched the detected CAF states.")
  }
  km_stats$selected <- km_stats$CAF_state %in% selected_states
} else {
  selected_states <- km_stats$CAF_state[isTRUE(km_stats$selected) | km_stats$selected %in% TRUE]
  selected_states <- as.character(selected_states)
}

if (length(selected_states) == 0) {
  stop(
    "No prognosis-related CAF state was selected.\n",
    "You can either inspect caf_km_statistics.tsv or set CFG$force_selected_states manually."
  )
}

candidate_membership <- marker_df[marker_df$cluster %in% selected_states, c("gene", "cluster"), drop = FALSE]
candidate_membership <- unique(candidate_membership)
candidate_membership <- candidate_membership[order(candidate_membership$cluster, candidate_membership$gene), , drop = FALSE]

candidate_genes <- unique(candidate_membership$gene)

write.table(
  candidate_membership,
  file = file.path(CFG$outdir, "candidate_gene_membership.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  data.frame(gene = candidate_genes, stringsAsFactors = FALSE),
  file = file.path(CFG$outdir, "CAF_candidate_gene_list.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

# Also write a root-level copy for direct downstream compatibility
write.table(
  data.frame(gene = candidate_genes, stringsAsFactors = FALSE),
  file = CFG$output_candidate_file,
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

summary_tab <- data.frame(
  item = c(
    "n_total_marker_rows",
    "n_unique_marker_genes",
    "n_detected_caf_states",
    "selected_states",
    "n_selected_states",
    "n_candidate_genes"
  ),
  value = c(
    nrow(marker_df),
    length(unique(marker_df$gene)),
    length(caf_states),
    paste(selected_states, collapse = "; "),
    length(selected_states),
    length(candidate_genes)
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_tab,
  file = file.path(CFG$outdir, "candidate_gene_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

msg("Done.")
msg("Selected prognosis-related CAF states: ", paste(selected_states, collapse = ", "))
msg("Candidate gene count: ", length(candidate_genes))
msg("Key output: ", normalizePath(file.path(CFG$outdir, "CAF_candidate_gene_list.txt"), winslash = "/"))
msg("Root-level compatibility copy: ", normalizePath(CFG$output_candidate_file, winslash = "/"))
