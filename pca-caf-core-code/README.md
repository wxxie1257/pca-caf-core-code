# pca-caf-core-code

Core R scripts for single-cell-informed CAF signature development, multicohort validation, clinicopathological integration, and downstream risk-group analyses in prostate cancer.

## Contents

This repository contains four core R scripts:

- `01_scRNA_CAF_core.R`  
  Single-cell RNA-seq processing, fibroblast extraction, and CAF subcluster identification.

- `02_CAF_score_candidate_genes.R`  
  Identification of prognosis-related CAF states and derivation of the CAF-related candidate-gene set.

- `03_signature_model_core.R`  
  Multicohort prognostic-model development and validation, including risk-score generation, TCGA clinicopathological integration, Cox regression, nomogram construction, calibration, and decision curve analysis (DCA).

- `04_riskgroup_key_analyses.R`  
  Downstream analyses of high-risk and low-risk groups, including functional, immune-related, and exploratory therapeutic-response analyses.

## Recommended order

Run the scripts in the following order:

1. `01_scRNA_CAF_core.R`
2. `02_CAF_score_candidate_genes.R`
3. `03_signature_model_core.R`
4. `04_riskgroup_key_analyses.R`

## Notes

- This repository is a compact, submission-oriented core-code package.
- The updated clinicopathological workflow, including nomogram, calibration, and DCA, has been integrated into `03_signature_model_core.R`.
- No separate `cli.R` file is included in this repository.
- Public datasets should be downloaded and prepared by users before running the scripts.
- File paths and input file names may need local adjustment before execution.

## Requirements

Developed in R (version 4.5.1).

Main packages used include:

Major packages used across the scripts include Seurat, harmony, SingleR, copykat, GSVA, clusterProfiler, survival, survminer, limma, timeROC, rms, Mime1, and oncoPredict.

## Contact

For questions regarding the workflow, please contact the corresponding author of the associated manuscript.
