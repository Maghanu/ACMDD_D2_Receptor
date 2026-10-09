# DRD2 docking example inputs

These files are examples for exercising the docking workflow. They are
separate from `data/DopamineD2_compounds.csv`, the experimentally measured
DRD2 IC50/pIC50 dataset used by the current ML models. A sample SMILES file
must not be assumed to have measured DRD2 activity from its filename or from
its docking score.

| File | Purpose |
| --- | --- |
| `ligand.smi` | Risperidone-like positive/control example used during workflow development |
| `testligand.smi` | Ligand preparation and docking test input |
| `EGFRligandCHEMBL63786.smi`, `tailEGFRCHEMBL45068.smi`, `tail-5EGFRCHEMBL120564.smi` | External example structures used to test candidate inputs; not established DRD2 actives by virtue of being included here |
| `poor_candidate.smi` | Ethane (`CC`), a deliberately tiny, valid input for testing the workflow rather than a realistic hit |
| `6CM4.pdb` | Local structure-related file retained for reference; the main docking script downloads and processes PDB 6CM4 as mmCIF |

## Run an example

From the repository root in the configured WSL/Conda environment:

```bash
python vina_autodock.py
```

At the prompt, enter a path such as
`docking/examples/drd2/ligand.smi`. The run is saved under
`docking/runs/<input-filename-without-extension>/`.

The long-term use case is to pass experimentally grounded, ML-prioritized
candidate SMILES to docking for structural triage. The current example inputs
do not constitute a candidate library, and the docking script does not yet
automatically consume ML predictions. For comparisons, keep the shared
prepared receptor, box, and settings fixed, inspect the poses, and confirm
activity experimentally.
