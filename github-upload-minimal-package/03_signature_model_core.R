# =========================================================
# File: 03_signature_model_core.R
# Purpose:
#   Curated, submission-ready script for multicohort prognostic-model
#   development, validation, risk-score generation, clinicopathological
#   integration, and nomogram construction in prostate cancer.
#
# What this script does:
#   1. Reads the TCGA-PRAD training cohort and five external validation cohorts.
#   2. Uses CAF_candidate_gene_list.txt as the biologically constrained
#      candidate-gene input for model development.
#   3. Runs the 101-combination machine-learning survival framework
#      implemented in Mime1.
#   4. Fixes the final selected model to RSF+stepCOX[both] by default,
#      consistent with the current manuscript and supplementary materials.
#   5. Exports all-cohort riskscore tables, TCGA-specific risk files, and
#      TCGA expression matrices for downstream analyses.
#   6. Generates Kaplan-Meier and time-dependent ROC plots for the training
#      and validation cohorts.
#   7. Performs univariable/multivariable Cox analyses and builds a simplified
#      4-variable nomogram in the TCGA cohort.
#
# Main required input files:
#   1. TCGA_expression.txt
#   2. time.TCGA.txt
#   3. GSE70768_expression.txt
#   4. time.GSE70768.txt
#   5. GSE70769_expression.txt
#   6. time.GSE70769.txt
#   7. GSE46602_expression.txt
#   8. time.GSE46602.txt
#   9. MSKCC_expression.txt
#   10. time.MSKCC.txt
#   11. DKFZ_expression.txt
#   12. time.DKFZ.txt
#   13. CAF_candidate_gene_list.txt
#   14. clinical.txt
#
# Main output files:
#   1. results/00_tables/Cindex.res.txt
#   2. results/00_tables/unicox.txt
#   3. results/00_tables/best_model_genes.txt
#   4. results/00_tables/riskscore.txt
#   5. results/00_tables/riskscore.selected_model.allcohort.txt
#   6. results/00_tables/risk.TCGA.txt
#   7. results/00_tables/risk.TCGA.extended.txt
#   8. results/00_tables/tcga_dat_T.txt
#   9. results/03_clinical_nomogram/*
#   10. results/99_parameters/important_parameters.tsv
#
# Notes:
#   - This script intentionally excludes GO/KEGG enrichment, GSEA, TIDE,
#     immune subtype, immune checkpoint, and drug-response analyses.
#     Those downstream modules belong in 04_riskgroup_key_analyses.R.
#   - This script depends on the Mime1 package for the 101-combination
#     survival machine-learning workflow.
#   - The primary endpoint is handled as bRFS/BCR. Legacy OS/OS.time or
#     futime/fustat column names are supported as input aliases.
# =========================================================

options(stringsAsFactors = FALSE, check.names = FALSE)
set.seed(1234)

CFG <- list(
  ROOT_DIR = getwd(),
  RESULTS_DIR = file.path(getwd(), "results"),
  FINAL_MODEL = "RSF+stepCOX[both]",
  AUTO_DETECT_FINAL_MODEL = FALSE,
  ML_SEED = 1234,
  UNICOX_FILTER_FOR_CANDI = TRUE,
  UNICOX_P_CUTOFF = 0.05,
  DATASET_ORDER = c("TCGA", "GSE70768", "GSE70769", "GSE46602", "MSKCC", "DKFZ"),
  ROC_YEARS = c(1, 3, 5),
  DAYS_PER_YEAR = 365,
  KM_WIDTH = 7.5,
  KM_HEIGHT = 7,
  ROC_WIDTH = 7,
  ROC_HEIGHT = 7,
  PDF_WIDTH = 7,
  PDF_HEIGHT = 6,
  TIFF_DPI = 600,
  RISK_COL = c("Low" = "#2F7FC1FF", "High" = "#D94738FF"),
  ROC_COL = c("#2F7FC1FF", "#3BA272FF", "#D94738FF", "#845EC2", "#F28E2B", "#76B7B2", "#59A14F"),
  PUB_BASE_SIZE = 12,
  ENCODINGS_TO_TRY = c("UTF-8", "GB18030", "GBK", "CP936", "latin1")
)

DIRS <- list(
  tables = file.path(CFG$RESULTS_DIR, "00_tables"),
  ml = file.path(CFG$RESULTS_DIR, "01_machine_learning"),
  perf = file.path(CFG$RESULTS_DIR, "02_performance"),
  clinical = file.path(CFG$RESULTS_DIR, "03_clinical_nomogram"),
  params = file.path(CFG$RESULTS_DIR, "99_parameters")
)

invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

`%||%` <- function(a, b) if (!is.null(a)) a else b

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")

check_and_load_packages <- function(pkgs){
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing_pkgs) > 0){
    stop(
      "Missing required packages: ",
      paste(missing_pkgs, collapse = ", "),
      "\nPlease install them before running this script.\n",
      "CRAN packages: install.packages(...); Bioconductor packages: BiocManager::install(...)"
    )
  }
  suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
}

read_table_auto <- function(file, sep = "\t", header = TRUE, row.names = NULL, check.names = FALSE){
  last_err <- NULL
  for(enc in CFG$ENCODINGS_TO_TRY){
    x <- tryCatch(
      utils::read.table(
        file = file, header = header, sep = sep, check.names = check.names,
        quote = "\"", comment.char = "", fill = TRUE, fileEncoding = enc,
        row.names = row.names, stringsAsFactors = FALSE
      ),
      error = function(e) { last_err <<- e; NULL }
    )
    if(!is.null(x)) return(x)
  }
  stop("Failed to read file: ", file, "\nLast error: ", conditionMessage(last_err))
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

clean_relop <- function(x){
  x <- as.character(x)
  x <- gsub("＞", ">", x, fixed = TRUE)
  x <- gsub("＜", "<", x, fixed = TRUE)
  x <- gsub("≤", "<=", x, fixed = TRUE)
  x <- gsub("≥", ">=", x, fixed = TRUE)
  trimws(x)
}

clean_na_string <- function(x){
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "na", "Na", "nan", "unknown", "Unknown", "unknow", ".")] <- NA
  x
}

standardize_tcga_id <- function(x){
  x <- as.character(x)
  x <- gsub("\\.", "-", x)
  ifelse(grepl("^TCGA-", x), substr(x, 1, 12), x)
}

save_gg_dual <- function(plot_obj, filename_stub, width = CFG$PDF_WIDTH, height = CFG$PDF_HEIGHT, dpi = CFG$TIFF_DPI){
  pdf_file <- paste0(filename_stub, ".pdf")
  tiff_file <- paste0(filename_stub, ".tiff")
  ggplot2::ggsave(pdf_file, plot = plot_obj, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(tiff_file, plot = plot_obj, width = width, height = height, units = "in", dpi = dpi, compression = "lzw")
}

save_base_dual <- function(plot_fun, filename_stub, width = CFG$PDF_WIDTH, height = CFG$PDF_HEIGHT, dpi = CFG$TIFF_DPI){
  pdf_file <- paste0(filename_stub, ".pdf")
  tiff_file <- paste0(filename_stub, ".tiff")
  grDevices::pdf(pdf_file, width = width, height = height, family = "Helvetica")
  plot_fun()
  grDevices::dev.off()
  grDevices::tiff(tiff_file, width = width, height = height, units = "in", res = dpi, compression = "lzw", pointsize = 12)
  plot_fun()
  grDevices::dev.off()
}

get_surv_alias_df <- function(df){
  x <- df
  cn <- colnames(x)
  if("BCR" %in% cn && "BCR.time" %in% cn){
    x$fustat <- safe_numeric(x$BCR)
    x$futime <- safe_numeric(x$BCR.time)
    x$OS <- x$fustat
    x$OS.time <- x$futime
  } else if("bRFS" %in% cn && "bRFS.time" %in% cn){
    x$fustat <- safe_numeric(x$bRFS)
    x$futime <- safe_numeric(x$bRFS.time)
    x$OS <- x$fustat
    x$OS.time <- x$futime
  } else if("OS" %in% cn && "OS.time" %in% cn){
    x$fustat <- safe_numeric(x$OS)
    x$futime <- safe_numeric(x$OS.time)
    x$OS <- x$fustat
    x$OS.time <- x$futime
  } else if("fustat" %in% cn && "futime" %in% cn){
    x$OS <- safe_numeric(x$fustat)
    x$OS.time <- safe_numeric(x$futime)
    x$fustat <- x$OS
    x$futime <- x$OS.time
  } else if("status" %in% cn && "time" %in% cn){
    x$OS <- safe_numeric(x$status)
    x$OS.time <- safe_numeric(x$time)
    x$fustat <- x$OS
    x$futime <- x$OS.time
  } else {
    stop("No supported survival columns found. Expected one of: BCR/BCR.time, bRFS/bRFS.time, OS/OS.time, futime/fustat, time/status.")
  }
  x$OS <- safe_numeric(x$OS)
  x$OS.time <- safe_numeric(x$OS.time)
  x$fustat <- safe_numeric(x$fustat)
  x$futime <- safe_numeric(x$futime)
  x
}

select_clinical_features <- function(rt, preferred = c("Age", "PSA", "gleasonscore", "Gleason", "Tstage", "pathologic_N", "clinical_M", "residual_tumor")){
  vars <- intersect(preferred, colnames(rt))
  vars <- vars[vapply(vars, function(v){
    x <- rt[[v]]
    sum(!is.na(x)) >= 10 && length(unique(stats::na.omit(x))) >= 2
  }, logical(1))]
  vars
}

publish_theme <- function(){
  ggplot2::theme_bw(base_size = CFG$PUB_BASE_SIZE) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )
}

simple_forest_plot <- function(df, title = NULL){
  df$Variable <- factor(df$Variable, levels = rev(df$Variable))
  ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = HR, ymin = Lower95, ymax = Upper95)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    ggplot2::geom_errorbar(width = 0.2, colour = "#2F5597") +
    ggplot2::geom_point(size = 2.8, colour = "#2F5597") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = CFG$PUB_BASE_SIZE) +
    ggplot2::labs(x = NULL, y = "Hazard Ratio (95% CI)", title = title)
}

make_risk_group <- function(score, method = "median", surv_df = NULL){
  score <- safe_numeric(score)
  cutoff <- stats::median(score, na.rm = TRUE)
  if(identical(method, "surv_cutpoint") && !is.null(surv_df) && requireNamespace("survminer", quietly = TRUE)){
    tmp <- data.frame(OS = surv_df$OS, OS.time = surv_df$OS.time, riskScore = score)
    tmp <- tmp[stats::complete.cases(tmp), , drop = FALSE]
    if(nrow(tmp) >= 20 && length(unique(tmp$OS)) > 1){
      sc <- tryCatch(
        survminer::surv_cutpoint(tmp, time = "OS.time", event = "OS", variables = "riskScore"),
        error = function(e) NULL
      )
      if(!is.null(sc)) cutoff <- as.numeric(sc$cutpoint[1, 1])
    }
  }
  list(cutoff = cutoff, group = ifelse(score >= cutoff, "High", "Low"))
}

coerce_age_group <- function(x){
  if(is.factor(x)) x <- as.character(x)
  x0 <- clean_na_string(clean_relop(x))
  num <- safe_numeric(x0)
  out <- ifelse(!is.na(num), ifelse(num <= 65, "<=65", ">65"), x0)
  out <- ifelse(out %in% c("<=65", ">65"), out, NA)
  factor(out, levels = c("<=65", ">65"))
}

coerce_psa_group <- function(x){
  if(is.factor(x)) x <- as.character(x)
  x0 <- clean_na_string(clean_relop(x))
  num <- safe_numeric(x0)
  out <- ifelse(!is.na(num), ifelse(num <= 10, "<=10", ">10"), x0)
  out <- ifelse(out %in% c("<=10", ">10"), out, NA)
  factor(out, levels = c("<=10", ">10"))
}

coerce_gleason_group <- function(x){
  if(is.factor(x)) x <- as.character(x)
  x0 <- clean_na_string(clean_relop(x))
  x0 <- gsub("\\s+", "", x0)
  num <- safe_numeric(x0)
  out <- ifelse(!is.na(num), ifelse(num <= 7, "<=7", ">7"), x0)
  out <- ifelse(out %in% c("<=7", ">7"), out, NA)
  factor(out, levels = c("<=7", ">7"))
}

extract_model_coefficients <- function(res, final_model){
  obj <- tryCatch(res[["ml.res"]][[final_model]], error = function(e) NULL)
  if(is.null(obj)) return(NULL)
  candidates <- list(
    obj[["fit"]] %||% NULL,
    obj[["model"]] %||% NULL,
    obj[["final.model"]] %||% NULL,
    obj
  )
  for(one in candidates){
    if(is.null(one)) next
    cf <- tryCatch(stats::coef(one), error = function(e) NULL)
    if(!is.null(cf) && length(cf) > 0){
      tab <- data.frame(
        Gene = names(cf),
        Coefficient = as.numeric(cf),
        Direction = ifelse(as.numeric(cf) > 0, "Positive", "Negative"),
        stringsAsFactors = FALSE
      )
      return(tab)
    }
  }
  NULL
}

base_pkgs <- c(
  "Mime1", "limma", "survival", "survminer", "timeROC", "ggplot2", "ggpubr",
  "rms", "pec", "regplot", "survcomp", "data.table", "tidyverse"
)
check_and_load_packages(base_pkgs)

msg("Step 1/7: reading expression matrices, survival files, and candidate genes")

TCGA <- read.table("TCGA_expression.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
TCGA.time <- read.table("time.TCGA.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
GSE70768 <- read.table("GSE70768_expression.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
GSE70768.time <- read.table("time.GSE70768.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
GSE70769 <- read.table("GSE70769_expression.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
GSE70769.time <- read.table("time.GSE70769.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
GSE46602 <- read.table("GSE46602_expression.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
GSE46602.time <- read.table("time.GSE46602.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
MSKCC <- read.table("MSKCC_expression.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
MSKCC.time <- read.table("time.MSKCC.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
DKFZ <- read.table("DKFZ_expression.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
DKFZ.time <- read.table("time.DKFZ.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
genelist <- read.table("CAF_candidate_gene_list.txt", header = FALSE, sep = "\t", check.names = FALSE)

required_manifest <- data.frame(
  type = "base_required",
  file = c(
    "TCGA_expression.txt", "time.TCGA.txt",
    "GSE70768_expression.txt", "time.GSE70768.txt",
    "GSE70769_expression.txt", "time.GSE70769.txt",
    "GSE46602_expression.txt", "time.GSE46602.txt",
    "MSKCC_expression.txt", "time.MSKCC.txt",
    "DKFZ_expression.txt", "time.DKFZ.txt",
    "CAF_candidate_gene_list.txt", "clinical.txt"
  ),
  stringsAsFactors = FALSE
)
write.table(required_manifest, file.path(DIRS$params, "required_files_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

TCGA.gene <- rownames(TCGA)[rowMeans(TCGA) > 1]
GSE70768.gene <- rownames(GSE70768)[rowMeans(GSE70768) > 1]
GSE70769.gene <- rownames(GSE70769)[rowMeans(GSE70769) > 1]
GSE46602.gene <- rownames(GSE46602)[rowMeans(GSE46602) > 1]
MSKCC.gene <- rownames(MSKCC)[rowMeans(MSKCC) > 1]
DKFZ.gene <- rownames(DKFZ)[rowMeans(DKFZ) > 1]

samegene <- Reduce(intersect, list(
  TCGA.gene, GSE70768.gene, GSE70769.gene, GSE46602.gene, MSKCC.gene, DKFZ.gene, genelist[, 1]
))
msg("samegene count: ", length(samegene))
if(length(samegene) == 0) stop("No overlapping candidate genes were found across all cohorts.")

prepare_ml_dataset <- function(expr_mat, time_df, standardize_tcga = FALSE){
  expr_t <- t(expr_mat)
  if(standardize_tcga) rownames(expr_t) <- standardize_tcga_id(rownames(expr_t))
  if(standardize_tcga) rownames(time_df) <- standardize_tcga_id(rownames(time_df))
  sameSample <- intersect(rownames(expr_t), rownames(time_df))
  expr_t <- expr_t[sameSample, , drop = FALSE]
  time_df <- time_df[sameSample, , drop = FALSE]
  out <- cbind(rownames(expr_t), time_df, expr_t)
  colnames(out)[1] <- "ID"
  out
}

list_train_vali_Data <- list(
  TCGA = prepare_ml_dataset(TCGA, TCGA.time, standardize_tcga = TRUE),
  GSE70768 = prepare_ml_dataset(GSE70768, GSE70768.time),
  GSE70769 = prepare_ml_dataset(GSE70769, GSE70769.time),
  GSE46602 = prepare_ml_dataset(GSE46602, GSE46602.time),
  MSKCC = prepare_ml_dataset(MSKCC, MSKCC.time),
  DKFZ = prepare_ml_dataset(DKFZ, DKFZ.time)
)

msg("Step 2/7: running 101-combination machine-learning survival modeling")
res <- ML.Dev.Prog.Sig(
  train_data = list_train_vali_Data[[1]],
  list_train_vali_Data = list_train_vali_Data,
  unicox.filter.for.candi = CFG$UNICOX_FILTER_FOR_CANDI,
  unicox_p_cutoff = CFG$UNICOX_P_CUTOFF,
  candidate_genes = samegene,
  mode = "all",
  nodesize = 7,
  seed = CFG$ML_SEED
)
saveRDS(res, file.path(DIRS$tables, "res.rds"))

final_model <- CFG$FINAL_MODEL
if(CFG$AUTO_DETECT_FINAL_MODEL){
  msg("AUTO_DETECT_FINAL_MODEL = TRUE is not used in this curated version. Final model remains fixed as: ", final_model)
}
msg("Final model: ", final_model)

best_model_genes <- tryCatch(colnames(res[["ml.res"]][[final_model]][["dataX"]]), error = function(e) character(0))
if(length(best_model_genes) == 0) msg("Warning: best_model_genes could not be extracted from res[[\"ml.res\"]][[final_model]][[\"dataX\"]].")

write.table(res[["Cindex.res"]], file = file.path(DIRS$tables, "Cindex.res.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(res[["Sig.genes"]], file = file.path(DIRS$tables, "unicox.txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(data.frame(gene = best_model_genes), file = file.path(DIRS$tables, "best_model_genes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

coef_tab <- extract_model_coefficients(res, final_model)
if(!is.null(coef_tab)){
  write.table(coef_tab, file = file.path(DIRS$tables, "final_model_coefficients.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

msg("Step 3/7: exporting riskscore files and TCGA-specific downstream inputs")
RS_full <- data.frame()
for(i in seq_along(res[["riskscore"]])){
  model_name <- names(res[["riskscore"]])[i]
  for(l in seq_along(list_train_vali_Data)){
    one <- as.data.frame(res[["riskscore"]][[i]][[l]][, c(1, 4), drop = FALSE])
    colnames(one) <- c("ID", "RS")
    one$ID <- as.character(one$ID)
    one$ID <- ifelse(grepl("^TCGA-", one$ID), standardize_tcga_id(one$ID), one$ID)
    one$model <- model_name
    one$dataset <- names(list_train_vali_Data)[l]
    RS_full <- rbind(RS_full, one)
  }
}
RS_full$RS <- safe_numeric(RS_full$RS)
RS_full <- RS_full[order(RS_full$model, RS_full$dataset, RS_full$ID), ]
RS_legacy <- RS_full[, c("ID", "RS", "model")]
write.table(RS_legacy, file = file.path(DIRS$tables, "riskscore.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(RS_full, file = file.path(DIRS$tables, "riskscore.with_dataset.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

RS_selected <- subset(RS_full, model == final_model)
write.table(RS_selected, file = file.path(DIRS$tables, "riskscore.selected_model.allcohort.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

TCGA_time_std <- get_surv_alias_df(data.frame(ID = rownames(TCGA.time), TCGA.time, check.names = FALSE))
rownames(TCGA_time_std) <- standardize_tcga_id(TCGA_time_std$ID)
TCGA_time_std <- TCGA_time_std[!duplicated(rownames(TCGA_time_std)), , drop = FALSE]

RS_TCGA <- subset(RS_selected, dataset == "TCGA")
if(nrow(RS_TCGA) == 0) RS_TCGA <- subset(RS_selected, grepl("^TCGA-", ID))
RS_TCGA$ID <- standardize_tcga_id(RS_TCGA$ID)
RS_TCGA <- RS_TCGA[!duplicated(RS_TCGA$ID), , drop = FALSE]
rownames(RS_TCGA) <- RS_TCGA$ID

common_tcga <- intersect(rownames(TCGA_time_std), rownames(RS_TCGA))
if(length(common_tcga) == 0) stop("No overlap between TCGA riskscore and time.TCGA.txt sample IDs.")

risk_tcga_ext <- cbind(
  data.frame(ID = common_tcga, check.names = FALSE),
  RS_TCGA[common_tcga, c("RS", "model", "dataset"), drop = FALSE],
  TCGA_time_std[common_tcga, c("OS", "OS.time", "fustat", "futime"), drop = FALSE]
)
risk_tcga_ext$riskScore <- risk_tcga_ext$RS
risk_group_res <- make_risk_group(risk_tcga_ext$riskScore, method = "median", surv_df = risk_tcga_ext[, c("OS", "OS.time"), drop = FALSE])
risk_tcga_ext$Risk <- factor(risk_group_res$group, levels = c("Low", "High"))
risk_tcga_ext <- risk_tcga_ext[, c("ID", "OS", "OS.time", "fustat", "futime", "RS", "riskScore", "Risk", "model", "dataset")]

risk_tcga_legacy <- risk_tcga_ext[, c("ID", "futime", "fustat", "riskScore", "Risk")]
write.table(risk_tcga_legacy, file = file.path(DIRS$tables, "risk.TCGA.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(risk_tcga_ext, file = file.path(DIRS$tables, "risk.TCGA.extended.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(risk_tcga_ext[, c("ID", "RS", "model")], file = file.path(DIRS$tables, "riskscore.TCGA.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

cutoff_record <- data.frame(
  metric = c("final_model", "tcga_risk_cutoff_method", "tcga_risk_cutoff_value", "tcga_n_samples"),
  value = c(final_model, "median", as.character(risk_group_res$cutoff), as.character(nrow(risk_tcga_ext))),
  stringsAsFactors = FALSE
)
write.table(cutoff_record, file = file.path(DIRS$params, "tcga_risk_cutoff_record.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

tcga_expr_full <- data.frame(Gene = rownames(TCGA), TCGA, check.names = FALSE)
write.table(tcga_expr_full, file = file.path(DIRS$tables, "tcga_dat_T.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
if(length(best_model_genes) > 0){
  keep_genes <- intersect(rownames(TCGA), best_model_genes)
  tcga_expr_model <- data.frame(Gene = keep_genes, TCGA[keep_genes, , drop = FALSE], check.names = FALSE)
  write.table(tcga_expr_model, file = file.path(DIRS$tables, "tcga_dat_T.modelGenes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

important_parameters <- data.frame(
  Parameter = c(
    "Primary endpoint", "Training cohort", "Validation cohorts", "Candidate genes",
    "Univariable Cox filter", "Machine-learning algorithms", "Model combinations",
    "Final model", "Risk grouping strategy", "Time-dependent ROC years"
  ),
  Value = c(
    "bRFS", "TCGA-PRAD", "MSKCC; GSE46602; GSE70768; GSE70769; DKFZ",
    as.character(length(samegene)), "P < 0.05",
    "CoxBoost; Elastic Net; GBM; Lasso; plsRcox; Ridge; RSF; StepCox; SuperPC; survival-SVM",
    "101", final_model, "Cohort-specific median cutoff", paste(CFG$ROC_YEARS, collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write.table(important_parameters, file = file.path(DIRS$params, "important_parameters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

msg("Step 4/7: generating model-development and validation plots")

save_base_dual(
  function(){ cindex_dis_all(res, validate_set = names(list_train_vali_Data)[-1], order = names(list_train_vali_Data), width = 0.35) },
  file.path(DIRS$ml, "CAF_BCR_01_AllModels_Cindex_Overview"), width = 8, height = 15
)

save_base_dual(
  function(){ cindex_dis_select(res, model = final_model, order = names(list_train_vali_Data)) },
  file.path(DIRS$ml, "CAF_BCR_02_SelectedModel_Cindex_ByCohort"), width = 8, height = 5
)

for(i in seq_along(list_train_vali_Data)){
  dataset_i <- names(list_train_vali_Data)[i]
  save_base_dual(
    function(){
      rs_sur(
        res, model_name = final_model, dataset = dataset_i,
        median.line = "hv", cutoff = 0.5, conf.int = TRUE,
        xlab = "Day", pval.coord = c(1000, 0.9)
      )
    },
    file.path(DIRS$ml, paste0("CAF_BCR_03_KM_", dataset_i)),
    width = CFG$KM_WIDTH, height = CFG$KM_HEIGHT
  )
}

all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res, train_data = list_train_vali_Data[[1]], inputmatrix.list = list_train_vali_Data, mode = "all", AUC_time = 1, auc_cal_method = "KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res, train_data = list_train_vali_Data[[1]], inputmatrix.list = list_train_vali_Data, mode = "all", AUC_time = 3, auc_cal_method = "KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res, train_data = list_train_vali_Data[[1]], inputmatrix.list = list_train_vali_Data, mode = "all", AUC_time = 5, auc_cal_method = "KM")

for(yr in c(1, 3, 5)){
  auc_obj <- get(paste0("all.auc.", yr, "y"))
  save_base_dual(
    function(){
      auc_dis_all(auc_obj, dataset = names(list_train_vali_Data), validate_set = names(list_train_vali_Data)[-1], order = names(list_train_vali_Data), width = 0.35, year = yr)
    },
    file.path(DIRS$ml, paste0("CAF_BCR_04_AllModels_AUC_", yr, "year")), width = 8, height = 15
  )
  save_base_dual(
    function(){
      roc_vis(auc_obj, model_name = final_model, dataset = names(list_train_vali_Data), order = names(list_train_vali_Data), anno_position = c(0.65, 0.55), year = yr)
    },
    file.path(DIRS$ml, paste0("CAF_BCR_05_SelectedModel_ROC_", yr, "year")), width = 8, height = 8
  )
}

save_base_dual(
  function(){
    auc_dis_select(list(all.auc.1y, all.auc.3y, all.auc.5y), model_name = final_model, dataset = names(list_train_vali_Data), order = names(list_train_vali_Data), year = c(1, 3, 5))
  },
  file.path(DIRS$ml, "CAF_BCR_06_SelectedModel_AUC_1_3_5year"), width = 8, height = 4
)

msg("Step 5/7: TCGA risk-group performance and clinical integration")
risk <- read.table(file.path(DIRS$tables, "risk.TCGA.txt"), header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
risk <- get_surv_alias_df(risk)
rownames(risk) <- standardize_tcga_id(rownames(risk))
risk <- risk[!duplicated(rownames(risk)), , drop = FALSE]
risk$Risk <- factor(as.character(risk$Risk), levels = c("Low", "High"))

cli <- read_table_auto("clinical.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
rownames(cli) <- standardize_tcga_id(rownames(cli))
cli <- cli[!duplicated(rownames(cli)), , drop = FALSE]
cli[] <- lapply(cli, function(x) if(is.character(x)) clean_na_string(clean_relop(x)) else x)

if("Age" %in% colnames(cli)) cli$Age <- coerce_age_group(cli$Age)
if("PSA" %in% colnames(cli)) cli$PSA <- coerce_psa_group(cli$PSA)
if("gleasonscore" %in% colnames(cli)) cli$gleasonscore <- coerce_gleason_group(cli$gleasonscore)
if("Gleason" %in% colnames(cli)) cli$Gleason <- coerce_gleason_group(cli$Gleason)

fac_cols <- intersect(c("clinical_M", "pathologic_N", "residual_tumor", "Age", "gleasonscore", "Gleason", "PSA", "Tstage"), colnames(cli))
cli[fac_cols] <- lapply(cli[fac_cols], factor)

common_cli <- intersect(rownames(risk), rownames(cli))
if(length(common_cli) == 0) stop("No overlap between risk.TCGA.txt and clinical.txt sample IDs.")
rt_all <- cbind(risk[common_cli, , drop = FALSE], cli[common_cli, , drop = FALSE])

# 5.1 KM curve in TCGA
rt_km <- rt_all[, , drop = FALSE]
rt_km <- rt_km[complete.cases(rt_km[, c("futime", "fustat", "riskScore", "Risk"), drop = FALSE]), , drop = FALSE]
rt_km$Risk <- droplevels(factor(rt_km$Risk, levels = c("Low", "High")))

if(nrow(rt_km) >= 10 && length(unique(rt_km$Risk)) == 2){
  fit_km <- survival::survfit(survival::Surv(OS.time, OS) ~ Risk, data = rt_km)
  p_km <- survminer::ggsurvplot(
    fit_km,
    data = rt_km,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    palette = unname(CFG$RISK_COL[c("Low", "High")]),
    legend.title = "Risk group",
    legend.labs = c("Low", "High"),
    xlab = "Time (days)",
    ylab = "bRFS probability"
  )
  save_gg_dual(p_km$plot + publish_theme(), file.path(DIRS$perf, "CAF_BCR_TCGA_01_KM_RiskGroup"), width = 7, height = 6)
}

# 5.2 time-dependent ROC in TCGA
roc_times_full <- CFG$ROC_YEARS * CFG$DAYS_PER_YEAR
roc_keep <- roc_times_full < max(rt_km$futime, na.rm = TRUE)
roc_times <- roc_times_full[roc_keep]
roc_years <- CFG$ROC_YEARS[roc_keep]
if(nrow(rt_km) >= 20 && length(unique(rt_km$fustat)) >= 2 && length(roc_times) >= 1){
  ROC_rt <- timeROC::timeROC(
    T = rt_km$futime,
    delta = rt_km$fustat,
    marker = rt_km$riskScore,
    cause = 1,
    weighting = "aalen",
    times = roc_times,
    ROC = TRUE
  )
  roc_plot_fun <- function(){
    for(i in seq_along(roc_times)){
      plot(ROC_rt, time = roc_times[i], col = CFG$ROC_COL[i], add = i > 1, title = FALSE, lwd = 3)
    }
    legend("bottomright", legend = paste0("AUC at ", roc_years, ifelse(roc_years > 1, " years: ", " year: "), sprintf("%.3f", ROC_rt$AUC)), col = CFG$ROC_COL[seq_along(roc_times)], lwd = 3, bty = "n")
  }
  save_base_dual(roc_plot_fun, file.path(DIRS$perf, "CAF_BCR_TCGA_02_TimeROC_1_3_5year"), width = 6, height = 6)
}

# 5.3 Clinical integration / Cox
clinical_vars <- select_clinical_features(rt_all)
rt_cli <- rt_all[, unique(c("futime", "fustat", "riskScore", clinical_vars)), drop = FALSE]
rt_cli <- rt_cli[complete.cases(rt_cli), , drop = FALSE]
rt_cli <- droplevels(rt_cli)

uni_res <- data.frame()
all_uni_vars <- unique(c("riskScore", clinical_vars))
for(v in all_uni_vars){
  dat_i <- rt_all[, intersect(c("futime", "fustat", v), colnames(rt_all)), drop = FALSE]
  dat_i <- dat_i[stats::complete.cases(dat_i), , drop = FALSE]
  if(nrow(dat_i) < 10) next
  fit_i <- tryCatch(survival::coxph(stats::as.formula(paste0("Surv(futime, fustat) ~ `", v, "`")), data = dat_i), error = function(e) NULL)
  if(is.null(fit_i)) next
  s <- summary(fit_i)
  if(nrow(s$coefficients) == 0) next
  uni_res <- rbind(uni_res, data.frame(
    Variable = v,
    HR = s$coefficients[1, "exp(coef)"],
    Lower95 = s$conf.int[1, "lower .95"],
    Upper95 = s$conf.int[1, "upper .95"],
    pvalue = s$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  ))
}
write.table(uni_res, file = file.path(DIRS$clinical, "univariate_cox_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if(nrow(uni_res) > 0){
  save_gg_dual(simple_forest_plot(uni_res, title = "Univariable Cox regression") + publish_theme(), file.path(DIRS$clinical, "CAF_BCR_TCGA_03_UnivariateCox_Forest"), width = 7, height = 4.8)
}

multi_res <- NULL
res.cox <- NULL
if(length(clinical_vars) > 0 && nrow(rt_cli) >= 30){
  multi_formula <- stats::as.formula(paste0("Surv(futime, fustat) ~ riskScore + ", paste(paste0("`", clinical_vars, "`"), collapse = " + ")))
  res.cox <- survival::coxph(multi_formula, data = rt_cli, x = TRUE, y = TRUE, model = TRUE)
  s <- summary(res.cox)
  multi_res <- data.frame(
    Variable = rownames(s$coefficients),
    HR = s$coefficients[, "exp(coef)"],
    Lower95 = s$conf.int[, "lower .95"],
    Upper95 = s$conf.int[, "upper .95"],
    pvalue = s$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  write.table(multi_res, file = file.path(DIRS$clinical, "multivariable_cox_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  save_gg_dual(simple_forest_plot(multi_res, title = "Multivariable Cox regression") + publish_theme(), file.path(DIRS$clinical, "CAF_BCR_TCGA_04_MultivariableCox_Forest"), width = 8, height = 5.5)
}

msg("Step 6/7: building the simplified 4-variable nomogram")
if(!is.null(res.cox)){
  rt_nomo <- rt_cli[, unique(c("futime", "fustat", "riskScore", clinical_vars)), drop = FALSE]
  rt_nomo <- rt_nomo[complete.cases(rt_nomo), , drop = FALSE]
  if(nrow(rt_nomo) >= 30){
    failtime <- c(365, 1095, 1825)
    prefix_nomo <- file.path(DIRS$clinical, "CAF_BCR_TCGA_07_Nomogram")
    prefix_cal  <- file.path(DIRS$clinical, "CAF_BCR_TCGA_08_Calibration_1_3_5year")
    prefix_risk <- file.path(DIRS$clinical, "CAF_BCR_TCGA_07_Nomogram_RiskTable")
    prefix_cidx <- file.path(DIRS$clinical, "CAF_BCR_TCGA_08_Calibration_Cindex")

    save_base_dual(
      plot_fun = function() {
        regplot::regplot(
          res.cox,
          plots = c("density", "boxes"),
          dencol = "#9467BDFF",
          boxcol = "#98DF8AFF",
          clickable = FALSE,
          title = "",
          points = FALSE,
          droplines = FALSE,
          observation = FALSE,
          rank = "sd",
          failtime = failtime,
          prfail = FALSE
        )
      },
      filename_stub = prefix_nomo,
      width = 8,
      height = 6
    )

    nomoRisk <- stats::predict(res.cox, newdata = rt_nomo, type = "risk")
    rt_out <- cbind(ID = rownames(rt_nomo), rt_nomo, Nomogram = nomoRisk)
    write.table(rt_out, file = paste0(prefix_risk, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

    dd <- rms::datadist(rt_out)
    options(datadist = "dd")
    m_size <- max(20, floor(nrow(rt_out) / 3))

    save_base_dual(
      plot_fun = function() {
        f1 <- rms::cph(survival::Surv(futime, fustat) ~ Nomogram, x = TRUE, y = TRUE, surv = TRUE, data = rt_out, time.inc = failtime[1])
        cal1 <- rms::calibrate(f1, cmethod = "KM", method = "boot", u = failtime[1], m = m_size, B = 1000)
        plot(cal1, xlim = c(0, 1), ylim = c(0, 1), xlab = "Nomogram-predicted bRFS (%)", ylab = "Observed bRFS (%)", lwd = 1.5, col = "green", sub = FALSE)

        f3 <- rms::cph(survival::Surv(futime, fustat) ~ Nomogram, x = TRUE, y = TRUE, surv = TRUE, data = rt_out, time.inc = failtime[2])
        cal3 <- rms::calibrate(f3, cmethod = "KM", method = "boot", u = failtime[2], m = m_size, B = 1000)
        plot(cal3, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", lwd = 1.5, col = "blue", sub = FALSE, add = TRUE)

        f5 <- rms::cph(survival::Surv(futime, fustat) ~ Nomogram, x = TRUE, y = TRUE, surv = TRUE, data = rt_out, time.inc = failtime[3])
        cal5 <- rms::calibrate(f5, cmethod = "KM", method = "boot", u = failtime[3], m = m_size, B = 1000)
        plot(cal5, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", lwd = 1.5, col = "red", sub = FALSE, add = TRUE)

        legend("topleft", legend = c("1-year", "3-year", "5-year"), col = c("green", "blue", "red"), lwd = 1.5, bty = "n")

        cindex <- survcomp::concordance.index(x = nomoRisk, surv.time = rt_out$futime, surv.event = rt_out$fustat, method = "noether")
        cindexLabel <- paste0(sprintf("%.03f", cindex$c.index), " (95% CI: ", sprintf("%.03f", cindex$lower), "-", sprintf("%.03f", cindex$upper), ")")
        text(0.35, 0.12, "C-index:", cex = 0.95)
        text(0.62, 0.06, cindexLabel, cex = 0.9)
      },
      filename_stub = prefix_cal,
      width = 5,
      height = 5
    )

    cindex <- survcomp::concordance.index(x = nomoRisk, surv.time = rt_out$futime, surv.event = rt_out$fustat, method = "noether")
    cindex_tab <- data.frame(Model = "Nomogram", C_index = sprintf("%.03f", cindex$c.index), CI_lower = sprintf("%.03f", cindex$lower), CI_upper = sprintf("%.03f", cindex$upper), stringsAsFactors = FALSE)
    write.table(cindex_tab, file = paste0(prefix_cidx, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

msg("Step 7/7: finished. Core outputs are ready for downstream risk-group analyses.")
msg("Please check:")
msg(" - ", file.path(DIRS$tables, "risk.TCGA.txt"))
msg(" - ", file.path(DIRS$tables, "tcga_dat_T.txt"))
msg(" - ", file.path(DIRS$tables, "best_model_genes.txt"))
msg(" - ", DIRS$ml)
msg(" - ", DIRS$clinical)
