# Data directory

This folder stores the raw and intermediate datasets used throughout the DRD2 computational drug discovery workflow.

## Contents

- `DopamineD2_compounds.csv` — primary IC50 dataset with SMILES and activity values
- `DopamineD2_compounds_Kd.csv` — Kd dataset used to compare potency measures
- `DopamineD2_combined.csv` — merged IC50/Kd dataset aligned by SMILES
- `DopamineD2_maccs.csv`, `DopamineD2_morgan2.csv`, `DopamineD2_morgan3.csv` — fingerprint-encoded matrices for modeling
- `ml_results_single_split.csv` and `ml_results_crossvalidation.csv` — classification model results
- `reg_results_single_split.csv` and `reg_results_crossvalidation.csv` — regression model results
- `*_roc_*.png` and `pred_vs_measured_*.png` — performance visualizations

## Notes

- This folder is the main data source for the project pipeline.
- Keep raw datasets unchanged when possible; generate derived files in a reproducible way via the scripts in `src/`.
