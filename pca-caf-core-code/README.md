# pca-caf-core-code

Core R scripts for single-cell-informed CAF signature development, multicohort validation, clinicopathological integration, and downstream risk-group analyses in prostate cancer.

## Overview

This repository contains a compact, submission-oriented code package supporting the main analytical workflow of our prostate cancer study. The project integrates single-cell-informed cancer-associated fibroblast (CAF) biology with bulk transcriptomic model development, multicohort validation, clinicopathological integration, and downstream characterization of risk groups.

This repository is intentionally streamlined and includes **four core R scripts only**. The updated clinicopathological workflow, including **univariable Cox analysis, multivariable Cox analysis, nomogram construction, calibration analysis, and decision curve analysis (DCA)**, has been incorporated into `03_signature_model_core.R`.

## Repository contents

The repository contains the following four scripts:

- `01_scRNA_CAF_core.R`  
  Core single-cell RNA-seq processing, fibroblast extraction, CAF subcluster identification, and related upstream analyses.

- `02_CAF_score_candidate_genes.R`  
  Derivation of prognosis-related CAF states and generation of the CAF-related candidate-gene set.

- `03_signature_model_core.R`  
  Core multicohort prognostic-model workflow, including:
  - candidate-gene restriction
  - machine-learning survival model development
  - multicohort validation
  - risk-score generation
  - TCGA clinicopathological integration
  - Cox regression analyses
  - clinicomolecular nomogram construction
  - calibration analysis
  - decision curve analysis (DCA)

- `04_riskgroup_key_analyses.R`  
  Downstream analyses of the high-risk and low-risk groups, including selected biological and immune-related comparisons and exploratory therapeutic-response analyses.

## Analytical workflow

The recommended running order is:

1. `01_scRNA_CAF_core.R`
2. `02_CAF_score_candidate_genes.R`
3. `03_signature_model_core.R`
4. `04_riskgroup_key_analyses.R`

## Main study scope

The code package supports a study framework that:

1. resolves CAF heterogeneity from single-cell RNA-seq data,
2. identifies prognosis-related CAF states,
3. derives a CAF-related candidate-gene set,
4. develops and validates a multicohort prognostic signature for biochemical recurrence-free survival (bRFS),
5. integrates the molecular risk score with clinicopathological variables,
6. evaluates the resulting model using Cox regression, nomogram analysis, calibration, and DCA,
7. performs downstream functional, immune-related, and exploratory response analyses.

## Expected input files

The scripts assume that users provide appropriately formatted input files in the working directory. Depending on the script, expected inputs may include:

### For single-cell and candidate-gene steps
- single-cell expression matrix
- cell annotation / metadata files
- marker-gene or fibroblast-related intermediate files

### For multicohort prognostic modeling
- `TCGA_expression.txt`
- `time.TCGA.txt`
- `GSE70768_expression.txt`
- `time.GSE70768.txt`
- `GSE70769_expression.txt`
- `time.GSE70769.txt`
- `GSE46602_expression.txt`
- `time.GSE46602.txt`
- `MSKCC_expression.txt`
- `time.MSKCC.txt`
- `DKFZ_expression.txt`
- `time.DKFZ.txt`
- `CAF_candidate_gene_list.txt`
- `clinical.txt`

### For the updated clinicopathological module in `03_signature_model_core.R`
The integrated clinical workflow may additionally require a TCGA clinical integration file containing variables such as:

- `ID`
- `riskScore`
- `age_group`
- `gleason`
- `psa`
- `pathologic_t_stage`
- `pathologic_n_stage`
- `futime`
- `fustat`

Column names and formatting should be checked carefully before running the script.

## Main outputs

Depending on the input data and script settings, the repository can generate outputs such as:

- CAF subcluster results
- candidate-gene lists
- machine-learning model comparison tables
- selected model coefficients
- multicohort risk-score files
- Kaplan–Meier survival plots
- time-dependent ROC plots
- Cox regression result tables
- clinicomolecular nomogram plots
- calibration plots
- DCA plots
- downstream risk-group analysis outputs

The exact output structure depends on the current script settings and file organization in the local working directory.

## Software requirements

This code was developed in **R (version 4.3.x)**.

Required R packages include, but are not limited to:

- `survival`
- `survminer`
- `timeROC`
- `ggplot2`
- `ggpubr`
- `rms`
- `pec`
- `regplot`
- `survcomp`
- `data.table`
- `tidyverse`
- `Hmisc`
- `dplyr`
- `patchwork`
- `dcurves`
- `limma`
- `Mime1`

Some packages may need to be installed separately before running the scripts.

## Notes on reproducibility

This repository is intended as a **core-code package** for the main analytical framework. It does not include raw source datasets or all intermediate files. Users should obtain the original public datasets independently and prepare them in the expected format before running the scripts.

The scripts were organized for transparency and reproducibility of the main workflow, but users may still need to adjust:

- working directories
- file paths
- input file names
- plotting parameters
- local package installation settings

before successful execution in a new environment.

## Important repository note

This repository contains **four core scripts only**.

The updated clinicopathological integration workflow, including the revised nomogram, calibration, and DCA procedures, is already integrated into:

- `03_signature_model_core.R`

Therefore, no separate `cli.R` file is included in this repository.

## Study context

This code package corresponds to a retrospective, multicohort transcriptomic workflow in prostate cancer and should be interpreted as a research-oriented analytical framework rather than a locked clinical decision tool.

## Citation

If you use this repository or adapt parts of the code, please cite the associated study once published.

## Contact

For questions regarding the analytical workflow, please contact the corresponding author listed in the associated manuscript.
