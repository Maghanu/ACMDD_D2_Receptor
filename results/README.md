# Machine-learning results

This directory contains outputs from the current models trained on
experimentally measured DRD2 IC50-derived pIC50 data and molecular
fingerprints.

## Contents

- `ml_results_single_split.csv`, `ml_results_crossvalidation.csv` — active /
  inactive classification results
- `reg_results_single_split.csv`, `reg_results_crossvalidation.csv` —
  regression results for pIC50
- `roc_*.png` — classifier ROC plots
- `pred_vs_measured_*.png` — regression prediction plots

The classifier's active label is defined by the current code as
`pIC50 >= 6.3`. The regression target is the measured pIC50 value. The Kd table
is not used as a label by these models.

These results characterize the current dataset and evaluation splits; they do
not establish performance on novel chemical scaffolds, external candidate
libraries, or in vitro experiments. Before using model output to select
candidate SMILES, evaluate generalization and preserve uncertainty and
provenance. Docking results, stored separately under `docking/runs/`, are
computational evidence and must not be presented as measured activity.
