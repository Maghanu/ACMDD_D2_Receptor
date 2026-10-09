# Source code directory

This folder contains the reproducible Python workflow for the D2 dopamine receptor study.

## Main scripts

- `step1_ic50_kd.py` — merges IC50 and Kd datasets by SMILES
- `step2_convert_to_fingerprint.py` — generates MACCS and Morgan fingerprints
- `step2_machine_learning.py` — binary classification workflow for active vs inactive compounds
- `step2_machine_learning_pic50.py` — regression workflow for predicting pIC50 values
- `step3_similar_scaffold.py` — scaffold similarity and maximum common substructure analysis

## Execution order

1. Run `step1_ic50_kd.py`
2. Run `step2_convert_to_fingerprint.py`
3. Run the machine learning scripts in `src/`
4. Use `step3_similar_scaffold.py` for scaffold-based analysis

## Requirements

Install project dependencies with:

```bash
pip install -r ../requirements.txt
```
