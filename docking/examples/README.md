# Docking example molecules

This directory contains small example inputs for testing the DRD2 docking
pipeline. They are not the experimentally measured DRD2 training dataset and
should not be interpreted as confirmed DRD2 ligands unless their experimental
activity and provenance are documented separately.

The project's current and future screening objective is to learn from
experimentally measured DRD2 activity, then apply validated models to
candidate SMILES. Docking examples are useful for checking file parsing,
ligand preparation, receptor preparation, and visualization—not for replacing
experimental labels or establishing target selectivity.

The workflow examples are in `drd2/`. To dock a molecule, follow
[`../README.md`](../README.md); outputs are written to
`docking/runs/<input-filename>/`.
