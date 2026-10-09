# Exploratory notebooks

These notebooks support interactive exploration and documentation of the
experimental DRD2 activity workflow.

## Notebooks

- `step1_ic50.ipynb` — IC50/pIC50 data retrieval and exploration
- `step1_kd.ipynb` — Kd data retrieval and exploration
- `step3_similar_scaffold.ipynb` — scaffold similarity and maximum common
  substructure analysis

The project uses experimentally measured compounds and assay endpoints as the
basis for its current ML work. Preserve endpoint and provenance information:
IC50-derived pIC50 and Kd are not interchangeable labels. For reproducible
preprocessing and model evaluation, use the scripts in `../src/`; future
screening of model-prioritized SMILES and docking are separate follow-up steps,
not experimental measurements.
