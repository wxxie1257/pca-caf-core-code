# GitHub update package

This package is only for the GitHub repository code update.

## Files
- `03_signature_model_core.R`
  - Updated version integrating the revised `cli.R`-style clinical workflow.
  - New clinical module:
    - univariable Cox
    - multivariable Cox
    - 6-variable nomogram
    - calibration
    - decision curve analysis (DCA)
  - Recommended additional input:
    - `TCGA_clinic_for_dca.txt`

- `cli.R`
  - Standalone revised clinical workflow script retained for reference.

## How to use in your repo
1. Replace the old `03_signature_model_core.R` with this updated file.
2. Add `cli.R` to the repository root or to the same analysis folder.
3. Update the repository README to mention that the clinical integration workflow was revised and now includes DCA.

## Suggested README sentence
The updated `03_signature_model_core.R` integrates the revised clinicopathological workflow (`cli.R` logic), including univariable and multivariable Cox analyses, 6-variable nomogram construction, calibration, and decision curve analysis in the TCGA cohort.
