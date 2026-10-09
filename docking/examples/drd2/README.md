# DRD2 docking example

Use this directory for the receptor and ligand files for a DRD2 docking run.

## Expected files

- `drd2_clean.pdb` — cleaned DRD2 receptor structure before PDBQT preparation
- `drd2.pdbqt` — receptor prepared for Vina with Meeko
- `ligand.pdbqt` — prepared ligand to dock
- `results/` — Vina poses and logs

RCSB PDB 6CM4 is one possible starting structure: it contains DRD2 bound to risperidone. Inspect the structure before preparation, remove the co-crystallized ligand from the receptor, and retain its coordinates separately for setting the docking box and validating the workflow by redocking.

PDB 1IEP is an Abl kinase structure, not DRD2. Do not use its coordinates or docking box for this example.

From the repository root, prepare the cleaned receptor using the center and size measured for the DRD2 binding site:

```bash
mk_prepare_receptor.py -i docking/examples/drd2/drd2_clean.pdb \
  -o docking/examples/drd2/drd2 -p -v \
  --box_size SIZE_X SIZE_Y SIZE_Z \
  --box_center CENTER_X CENTER_Y CENTER_Z
```

Then run the docking wrapper:

```bash
python vina_autodock.py \
  --receptor docking/examples/drd2/drd2.pdbqt \
  --ligand docking/examples/drd2/ligand.pdbqt \
  --center_x CENTER_X --center_y CENTER_Y --center_z CENTER_Z \
  --size_x SIZE_X --size_y SIZE_Y --size_z SIZE_Z \
  --out docking/examples/drd2/results/docked_pose.pdbqt \
  --log docking/examples/drd2/results/docking.log
```

Replace the uppercase placeholders with the chosen numeric box values. The receptor and ligand files are not included; add prepared inputs locally before running.
