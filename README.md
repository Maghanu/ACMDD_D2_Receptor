# DRD2 ligand discovery: experimental data, machine learning, and docking

## Project purpose

This project uses experimentally measured dopamine D2 receptor (DRD2) activity
to train and evaluate ligand-based machine-learning models. The longer-term
goal is to use those models to prioritize candidate molecules represented as
SMILES, then investigate promising candidates further with structure-based
docking and ultimately in vitro experiments.

The current ML workflow is based on measured DRD2 `pIC50` data. It trains
classification models for an activity threshold and regression models for
`pIC50`. The Kd measurements are curated and combined separately for
comparison; they are not currently used as labels by the ML scripts. Docking
is a complementary prioritization tool, not an experimental measurement and
not proof that a compound binds DRD2.

## Repository layout

- `data/` — experimental activity tables and derived molecular fingerprints
- `src/` — data preparation, fingerprinting, ML evaluation, and scaffold
  analysis scripts
- `notebooks/` — exploratory and tutorial-style analyses
- `results/` — outputs from the current ML scripts
- `docking/` — shared-receptor Vina workflow, example SMILES, and per-ligand
  docking outputs
- `Workflow_steps/` — earlier workflow materials retained for reference

## ML workflow

Run these commands from the repository root in an environment with the
dependencies in `requirements.txt` installed:

```bash
python src/step1_ic50_kd.py
python src/step2_convert_to_fingerprint.py
python src/step2_machine_learning.py
python src/step2_machine_learning_pic50.py
python src/step3_similar_scaffold.py
```

The classification workflow labels the experimentally measured IC50 compounds
using the `pIC50 >= 6.3` cutoff. The regression workflow predicts measured
`pIC50`. Both use MACCS, Morgan radius 2, and Morgan radius 3 fingerprints.
Model outputs are written to `results/`.

The current scripts train and evaluate models on known compounds; they do not
generate new molecules or automatically score an external SMILES library.
Applying a trained model to candidate SMILES and selecting compounds for
follow-up is a future screening step. Keep candidate structures and their
predictions traceable to their source and do not treat predictions as
experimental activity.

## Docking workflow

Docking provides a structural follow-up for selected molecules. The Vina
workflow accepts `.smi`, `.sdf`, or `.mol` input, prepares each ligand against
a shared prepared DRD2 receptor, and writes poses and visualizations under
`docking/runs/<ligand-name>/`. On Windows, run Vina in Ubuntu under WSL; setup
and usage are documented in [docking/README.md](./docking/README.md).

Docking scores are approximate scoring-function outputs. They are useful, at
most, as one piece of evidence when prioritizing model-selected candidates.
Compare candidates using the same receptor, box, and settings; inspect the
poses and ligand preparation; then use experimental testing to establish
activity. Existing example SMILES are for workflow testing unless their
experimental DRD2 activity is explicitly documented.

## Biological context

DRD2 is a G-protein-coupled receptor involved in dopaminergic signaling and an
important target in several neurological and psychiatric disorders. Its
pharmacology is complex, and measured activity depends on assay conditions,
endpoint, and ligand context. The project therefore retains the distinction
between experimental measurements, ML predictions, and docking scores.
