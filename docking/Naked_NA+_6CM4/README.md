# Ligand 1 docking against the supplied 6CM4 receptor

This note records the docking workflow and parameters for preparing a
reproducible methods section. Vina scores are approximate scoring-function
estimates, not experimental binding affinities or free energies.

## Input and output files

- Receptor input: [`6cm4(1).m.mmcif`](./6cm4(1).m.mmcif)
- Final Moorhen-cleaned structure:
  [`Naked_6CM4_NA_lig1.mmcif`](./Naked_6CM4_NA_lig1.mmcif)
- Pre-cleanup docked structure:
  [`ligand1_docked_6CM4.mmcif`](./ligand1_docked_6CM4.mmcif)
- Search poses and run details: [`docking/runs/ligand_1_supplied_6CM4/`](../runs/ligand_1_supplied_6CM4/)
- Repeated-search score table: [`docking_summary.csv`](../runs/ligand_1_supplied_6CM4/docking_summary.csv)
- Preparation metadata: [`preparation.json`](../runs/ligand_1_supplied_6CM4/preparation.json)

## Reproducible method

1. **Ligand identity and protonation.** The submitted ligand was specified as
   the neutral SMILES
   `CC1=C(C(=O)N2CCCCC2=N1)CCN3CCN(CC3)C4=C(C(=CC=C4)F)F`. For docking, one
   hydrogen was added to the aliphatic piperazine nitrogen, giving the
   singly charged state
   `CC1=C(C(=O)N2CCCCC2=N1)CC[NH+]3CCN(CC3)C4=C(C(=CC=C4)F)F`. The input and
   docking-state SMILES are both retained in the run directory. This was a
   chemically plausible preparation choice; no pKa calculation or
   experimental protonation measurement was performed.
2. **Ligand coordinates and atom typing.** RDKit generated a 3D conformer
   (ETKDG random seed 22), then MMFF-optimized it. Meeko generated the
   AutoDock ligand input, including explicit hydrogens and atom types.
3. **Receptor preparation.** Protein heavy atoms were read from the supplied
   mmCIF `ATOM` records (2,131 atoms); non-protein `HETATM` records were not
   used to prepare the Vina receptor. PDBFixer added missing atoms and
   OpenMM added hydrogens at pH 7.0 and energy-minimized the prepared protein
   with positional restraints. Meeko prepared the receptor PDBQT. The
   sodium ion and waters in the supplied structure were retained in the
   output mmCIF but were **not** included in Vina scoring. Meeko also reported
   unmatched terminal-residue templates outside the docking-box neighborhood;
   those distant residues were omitted from its receptor PDBQT preparation.
4. **Binding-site box.** The supplied “naked” mmCIF has no co-crystallized
   ligand to define a search box. The box was therefore calculated from the
   coordinates of native risperidone (component `8NU`) in the cached 6CM4
   reference structure, with 8 Å added to each ligand bounding-box dimension:

   | Parameter | Value (Å) |
   | --- | ---: |
   | Center x | 9.1765 |
   | Center y | 5.4705 |
   | Center z | -8.6645 |
   | Size x | 14.757 |
   | Size y | 17.715 |
   | Size z | 22.815 |

5. **Docking and pose selection.** AutoDock Vina 1.2.7 (standard `vina`
   scoring function) was run five times with independent random seeds,
   exhaustiveness 32, and up to 20 output poses per search. The pose with the
   lowest Vina affinity score across the five searches was selected. The
   search summary records each seed and its best score:

   | Seed | Best Vina score (kcal/mol) |
   | ---: | ---: |
   | 42 | -10.418 |
   | 2026 | -10.443 |
   | 271828 | -10.455 |
   | 314159 | **-10.465** |
   | 8675309 | -10.451 |

   Thus, the selected pose is from seed 314159. The five best scores span
   0.047 kcal/mol. This close agreement is a search-repeatability observation
   for the score, not independent experimental validation of the pose.
6. **Pose assembly and post-docking cleanup.** The selected pose's 28 heavy
   atoms were added to a copy of the supplied mmCIF as component `LIG`, chain
   `L`, with ligand atom and bond definitions. This untrimmed pose is retained
   as `ligand1_docked_6CM4.mmcif`; at this stage, the original 2,148 atom
   records (protein, waters, and sodium) were preserved. Hydrogens used in
   docking preparation are not included in the appended ligand heavy-atom
   coordinates.

   The final file, `Naked_6CM4_NA_lig1.mmcif`, was subsequently edited in
   Moorhen. The user removed loose residues at the original T4-lysozyme
   connection: chain A author-numbered residues 1015 (LEU), 1016 (LYS),
   1017 (ILE), 1159 (ASP), 1160 (ALA), and 1161 (TYR), 46 protein atoms in
   total. This is a post-docking structural cleanup, not part of receptor
   preparation or the Vina scoring calculation; the ligand was not redocked
   or rescored after the cleanup. Comparison of the cleaned and untrimmed
   files confirmed that all 28 ligand atom coordinates are unchanged. The
   cleaned mmCIF parses with Gemmi and contains one `LIG` residue with 28
   atoms.

## Software and caveats

The run used Python 3.11.17, AutoDock Vina 1.2.7, RDKit 2026.03.1, Meeko
0.8.0, Biopython 1.88, PDBFixer 1.12.0, OpenMM 8.6.1, and Gemmi 0.7.5.

The receptor was treated as rigid during Vina search. The preparation does not
explicitly score the retained sodium ion or crystallographic waters, and the
selected protonation state, receptor preparation, and docking score introduce
uncertainty. The score is an approximate Vina ranking metric for this setup,
not a measured binding free energy, affinity, or pIC50.
