# Experimental data and derived features

This directory contains the experimental DRD2 activity data used to build and
evaluate the current ligand-based ML models, together with curated intermediate
tables and fingerprints.

## Experimental activity tables

- `DopamineD2_compounds.csv` — compound SMILES and experimentally measured
  IC50-derived `pIC50` values used by the current classifier and regressor
- `DopamineD2_compounds_Kd.csv` — experimental Kd measurements and units
- `DopamineD2_combined.csv` — outer merge of the IC50 and Kd tables by exact
  SMILES; Kd is normalized to nM, and an empty cell means that endpoint is
  unavailable for that compound

IC50 and Kd are distinct assay endpoints and should not be treated as
interchangeable measurements. The current ML scripts use the `pIC50` values in
the IC50 dataset; they do not currently train on the Kd values.

## Fingerprints and results

- `DopamineD2_maccs.csv`, `DopamineD2_morgan2.csv`, and
  `DopamineD2_morgan3.csv` — molecular fingerprints derived from the IC50
  dataset, retaining SMILES and activity columns alongside fingerprint bits
- `P14416_updated/` and `P14416_updated.tar.gz` — additional source/reference
  materials retained with the project

Regenerate the fingerprint files with
`python src/step2_convert_to_fingerprint.py`. The active ML scripts write
model metrics and plots to `results/`; older result files kept in `data/` are
not the output location used by those scripts.

## Data handling

Preserve the original experimental values and source information. Document
filtering, unit conversions, duplicate handling, and assay selection when
creating a new analysis table. The current merge matches exact SMILES and
rejects duplicates within each source table; it is not a chemical-identity
standardization or repeated-measurement aggregation pipeline.

Future candidate SMILES scored by an ML model are predictions, not experimental
labels. Keep them in a separate, clearly identified screening table and do not
append predictions to the experimental training data as if they were measured
activities.
