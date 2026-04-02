# Curated Core Code for Single-Cell-Informed CAF Analysis in Prostate Cancer

## Overview

This repository contains a simplified, submission-ready code package for the principal analytical workflow of our prostate cancer study. Only the core scripts directly supporting the main conclusions are included.

The workflow covers four major parts:

1. Single-cell CAF-state identification, marker discovery, CopyKAT inference, and pathway activity analysis.
2. TCGA CAF subgroup-score analysis and CAF-related candidate-gene derivation.
3. Multicohort prognostic model development, validation, risk-score generation, and clinicopathological integration.
4. Manuscript-relevant downstream analyses between the high-risk and low-risk groups.

This repository is a curated submission version and does not include all historical intermediate files, exploratory side analyses, or duplicated plotting scripts from the original project.

## Files Included

- `01_scRNA_CAF_core.R`  
  Core single-cell analysis script for GSE185344, including quality control, Harmony integration, annotation, fibroblast extraction, CAF reclustering, marker identification, CopyKAT inference, and CAF-related pathway activity analysis.

- `02_CAF_score_candidate_genes.R`  
  TCGA CAF subgroup-score analysis script. It calculates CAF subgroup scores, compares normal and tumor tissues, performs Kaplan-Meier screening of prognosis-related CAF states, and exports `CAF_candidate_gene_list.txt`.

- `03_signature_model_core.R`  
  Core prognostic modeling script. It performs multicohort model development and validation, exports risk-score files, and runs clinicopathological integration and simplified nomogram analysis.

- `04_riskgroup_key_analyses.R`  
  Downstream analysis script for the TCGA training cohort, including GO/KEGG enrichment, GSEA, immune-related analyses, immune checkpoint comparison, and optional transcriptome-based predicted drug-response analysis.

- `README.md`  
  Repository overview and usage instructions.

## Public Datasets Used

The code is based on publicly available datasets from the original repositories, including:

- GSE185344
- TCGA-PRAD
- MSKCC
- GSE46602
- GSE70768
- GSE70769
- DKFZ

Large public raw or processed expression matrices are not bundled in this repository. Users should obtain the original datasets from the corresponding public sources and prepare the required input files before running the scripts.

## Running Order

Run the scripts in the following order:

```r
source("01_scRNA_CAF_core.R")
source("02_CAF_score_candidate_genes.R")
source("03_signature_model_core.R")
source("04_riskgroup_key_analyses.R")
```

## Key Workflow Logic

- `01_scRNA_CAF_core.R` generates CAF marker results from single-cell data.
- `02_CAF_score_candidate_genes.R` uses CAF marker genes to derive `CAF_candidate_gene_list.txt`.
- `03_signature_model_core.R` uses `CAF_candidate_gene_list.txt` for multicohort model development and exports `risk.TCGA.txt` and `tcga_dat_T.txt`.
- `04_riskgroup_key_analyses.R` uses `risk.TCGA.txt` and `tcga_dat_T.txt` for downstream analyses.

## Software Environment

This workflow was developed in R 4.3.2.

Major packages used across the scripts include Seurat, harmony, SingleR, copykat, GSVA, clusterProfiler, survival, survminer, limma, timeROC, rms, Mime1, and oncoPredict.

It is recommended to also provide a `sessionInfo.txt` or `package_versions.txt` file when sharing the code publicly.

## Notes

- This repository is intended to improve transparency and reproducibility of the principal analytical framework of the study.
- The code package is intentionally simplified for manuscript submission.
- Optional analyses in `04_riskgroup_key_analyses.R` require additional auxiliary files if users want to run those modules.

## Citation

If you use this repository, please cite the associated manuscript and acknowledge the original public datasets used in the study.
