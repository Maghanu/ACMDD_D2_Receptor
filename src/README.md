# Source code: DRD2 data and ML workflow

The scripts here prepare experimentally measured DRD2 activity data, encode
experimental molecules as SMILES-derived fingerprints, and evaluate
ligand-based models. These models are intended to support future screening of
candidate SMILES; the current scripts do not generate molecules or connect
predictions automatically to docking.

## Scripts and sequence

1. `step1_ic50_kd.py` — combines the experimental IC50/pIC50 and Kd tables by
   exact SMILES; Kd values are converted to nM for the combined table.
2. `step2_convert_to_fingerprint.py` — builds MACCS, Morgan radius 2, and
   Morgan radius 3 fingerprints from the experimental IC50 dataset. It also
   assigns the active label using `pIC50 >= 6.3`.
3. `step2_machine_learning.py` — evaluates active/inactive classifiers.
4. `step2_machine_learning_pic50.py` — evaluates regression models for
   experimentally measured pIC50.
5. `step3_similar_scaffold.py` — explores scaffold similarity within the
   experimental DRD2 compound set.

Run from the repository root:

```bash
python src/step1_ic50_kd.py
python src/step2_convert_to_fingerprint.py
python src/step2_machine_learning.py
python src/step2_machine_learning_pic50.py
python src/step3_similar_scaffold.py
```

The fingerprint generator reads `data/DopamineD2_compounds.csv` and writes
derived fingerprint tables to `data/`. The classifier and regressor read those
tables and write metrics and plots to `results/`.

## Interpretation and future screening

Training targets are experimental measurements, not docking scores. IC50 and
Kd are different endpoints; the current models use pIC50 for classification
and regression, while Kd is retained separately for comparison.

The eventual screening workflow is to apply a validated model to external
candidate SMILES, rank candidates for follow-up, and use docking and
experimental assays as complementary evidence. Predictions are hypotheses for
triage, not confirmed DRD2 activity. Evaluation here uses random splits and
cross-validation; it is not by itself evidence of performance on novel
chemical scaffolds.

Install project dependencies from the repository root with:

```bash
pip install -r requirements.txt
```
