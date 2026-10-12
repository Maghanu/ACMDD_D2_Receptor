# DRD2 docking for ML-prioritized candidate investigation

## Scientific purpose

This docking workflow is intended as a **structure-based follow-up** to the
project's ligand-based machine-learning (ML) stage. The ML models learn from
known compounds with experimentally measured DRD2 activity, currently using
IC50-derived pIC50 values. In the eventual screening workflow, those models
will be applied to molecules without known DRD2 measurements, represented as
SMILES, to identify candidates worth further investigation.

Docking provides a complementary question: **can a candidate adopt a
chemically plausible pose in the DRD2 binding pocket, and are its predicted
contacts consistent with the known binding site?** It can help inspect and
triage ML-prioritized structures before committing resources to experimental
testing. It does not confirm binding, predict selectivity, or replace an
in-vitro assay.

The two computational stages estimate different things. The ML stage predicts
an activity endpoint from patterns in known active compounds; Vina searches
for poses in a prepared receptor and assigns an approximate score using its
scoring function. This workflow uses the standard Vina scoring function, which
approximates favorable steric/shape complementarity, hydrophobic contacts, and
hydrogen bonding, while penalizing steric clashes and ligand flexibility.
These terms are simplified and do not represent every physical contribution to
binding. In particular, this workflow docks against a rigid receptor; it does
not sample receptor conformational changes or explicitly model binding-site
water molecules.

A Vina score reported in kcal/mol is a **scoring-function estimate**, not an
experimentally measured binding free energy, affinity, or pIC50. It should not
be interpreted as a calibrated thermodynamic value or combined arithmetically
with predicted pIC50. Docking ligands separately is the usual screening
workflow; scores are cautiously comparable for relative triage when the same
receptor, box, ligand-preparation policy, and scoring/search settings are used.
They are not guaranteed to provide a reliable rank, and small differences may
not be meaningful. A favorable docking score alone is not evidence that an
otherwise unprioritized compound is a DRD2 hit.

## Intended ML-to-docking-to-experiment workflow

1. **Train and evaluate the ML model** using curated, experimentally measured
   DRD2 compounds. Keep assay endpoint, source, units, and structure
   standardization traceable.
2. **Score candidate SMILES** with a validated model. The candidate set should
   be clearly separated from the experimental training data; predicted values
   must not be written back as if they were measured activity.
3. **Select candidates for structural follow-up.** Consider ML score,
   prediction confidence/applicability domain, chemical diversity, and
   practical feasibility rather than selecting on a single score alone.
4. **Prepare and dock candidate structures** using the same receptor, binding
   box, protonation/tautomer policy, and Vina settings. Review poses and
   receptor interactions for chemical plausibility; do not rank compounds by
   small score differences without checking pose quality and search stability.
5. **Prioritize a manageable set for in-vitro testing.** Experimental DRD2
   measurements determine whether a candidate is active and provide new
   evidence for later model improvement.

The current scripts implement docking of **one supplied ligand per run**.
They do not yet read ML predictions, screen a SMILES library, or automatically
combine ML and docking rankings. The user must select a candidate and provide
its structure file manually. Example SMILES under `docking/examples/drd2/` are
for workflow testing; their presence does not establish experimental DRD2
activity.

## Scope and limitations

- The Vina score approximates shape/steric fit, hydrophobic contacts, hydrogen
  bonding, and a ligand-flexibility penalty. It omits or simplifies important
  contributions to binding, and the receptor is rigid in this workflow. Use
  scores for pose generation and cautious relative triage under a consistent
  protocol, not as calibrated experimental affinities.
- The workflow uses the DRD2 conformation in PDB 6CM4, an X-ray structure with
  co-crystallized risperidone. A single receptor conformation cannot represent
  all receptor states or conformational flexibility.
- Ligand protonation, tautomerism, stereochemistry, salt removal, and
  coordinate generation affect preparation and can affect the result. Check
  these explicitly for each candidate.
- Docking outcomes depend on receptor preparation, box placement, search
  settings, and stochastic search. Use matched settings for comparisons,
  consider repeated runs for close candidates, and inspect the poses rather
  than interpreting scores alone.
- Correct ligand protonation and atom typing affect the interactions evaluated
  by the scoring function. Check the prepared ligand and receptor before
  interpreting scores.
- A plausible pose does not establish cellular activity, functional
  pharmacology, selectivity, exposure, safety, or experimental binding.

Treat ML predictions and docking results as **complementary hypotheses** for
choosing experiments. Preserve the model version, input SMILES, prepared
structure/protonation state, receptor cache, docking settings, pose, and score
for every candidate so decisions can be reproduced and compared.

## Environment

AutoDock Vina's Python bindings are not supported on native Windows. On
Windows, run the workflow in Ubuntu under WSL. Install Miniforge inside
Ubuntu, not into a Windows Python environment. From the Ubuntu terminal:

```bash
cd ~
curl -L -o Miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh
source ~/miniforge3/etc/profile.d/conda.sh
conda create -n drd2-vina -c conda-forge python=3.11 vina meeko rdkit biopython pdbfixer openmm matplotlib -y
conda activate drd2-vina
cd "/mnt/c/Users/manup/Documents/Leiden Universiteit/LUMC/Advanced Computational Methods in Drug Discovery/Mini-Research/ACMDD_D2_Receptor"
```

For later terminal sessions, activate the environment and return to the
repository root:

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate drd2-vina
cd "/mnt/c/Users/manup/Documents/Leiden Universiteit/LUMC/Advanced Computational Methods in Drug Discovery/Mini-Research/ACMDD_D2_Receptor"
```

Install Vina from conda-forge rather than pip; pip may try to build it from
source. The environment belongs in Ubuntu's home directory. The repository can
remain on the Windows-mounted drive.

## Dock one ML-prioritized or control molecule

Provide the selected candidate's SMILES in a `.smi` file (one molecule per file
for this workflow), or provide an `.sdf`/`.mol` file. Prefer a standardized
structure and choose a protonation/tautomer state appropriate to the assay
conditions. Retain the exact input structure alongside its ML prediction so
that the docked molecule can be traced back to the prioritized candidate. For
example:

```bash
python Workflow_steps/Workflow_step7_vina_autodock.py
```

When prompted, enter a path such as
`docking/examples/drd2/ligand.smi`. The script downloads RCSB PDB 6CM4,
derives the docking box from co-crystallized risperidone, repairs and
minimizes the receptor, prepares the ligand, and runs Vina.

The prepared DRD2 receptor and docking box are cached in
`docking/prepared/6CM4/` and reused for every ligand. Per-ligand run folders
receive copies of the shared receptor files, allowing candidates to be
compared against identical receptor coordinates and box settings. Remove that
cache only when deliberately rebuilding the shared receptor; after doing so,
redock all candidates being compared. This controls receptor preparation
between candidates, but does not remove uncertainty from ligand preparation
or Vina's search. Consider repeated runs or a more systematic protocol when
ranking close candidates.

Salt inputs with exactly one carbon-containing component have disconnected
non-carbon ions removed with a warning. Multiple carbon-containing components
are rejected rather than guessed. Inspect protonation, stereochemistry, and
the resulting pose.

Outputs are stored under `docking/runs/<ligand-file-name>/`, including
`docked_poses.pdbqt`, `vina.log`, and prepared receptor/ligand files. Record
which ML-prioritized candidate and input structure each run represents. Treat
the scores as a computational prioritization signal only. Check pose
plausibility, compare against appropriate known ligands and controls, and
experimentally test promising candidates in vitro.

After every successful docking run, the workflow creates or updates
`docking/runs/docking_summary.csv`. It records the run name and directory,
input file, docked SMILES when available, optional prediction metadata, the
best Vina score in the explicitly named `vina_score_kcal_mol` column, and a
UTC timestamp. The score description in each row clarifies that this is an
approximate Vina score, not an experimental affinity. Re-running a ligand
updates that run's row rather than adding a duplicate.

## Dock the top predicted-pIC50 candidates

The prediction CSV can be ranked by predicted pIC50 and its ten highest-scoring
SMILES docked against the same cached receptor:

```bash
python Workflow_steps/Workflow_step7_vina_autodock.py --top-predictions
```

By default, this reads
`Workflow_steps/data/predicted_pIC50_svr.csv` and selects ten candidates in
descending predicted-pIC50 order. The CSV must contain `smiles` and either a
`predicted_pIC50` or `pred_pIC50` column; `compound_id` is used to identify
each run when present. To select another CSV or change the count, pass
`--predictions-csv path/to/predictions.csv --top 10`.

Each candidate gets an individual directory under `docking/runs/` by default.
Pass `--output-dir path/to/results-folder` to keep a batch and its
`top_N_docking_summary.csv` together in a separate results folder. Input CSVs
may use `smiles` or `SMILES` and may identify candidates with `compound_id` or
`zincid`.

The batch summary records the prediction rank, compound ID, input SMILES,
predicted pIC50, Vina score, and run directory.
Each successful candidate run is also recorded in the shared
`docking/runs/docking_summary.csv` file alongside results from single-ligand
runs.
Predicted pIC50 and Vina scores are separate computational estimates; neither
is an experimental activity measurement.

## Visualize docking results

Run the plotter and enter the run name (the ligand filename without its
extension):

```bash
python docking/plot_docking.py
```

It generates both `visualization.png` and `visualization_3d.html` in that
ligand's run directory. For example, entering `ligand` creates them under
`docking/runs/ligand/`. The HTML uses 3Dmol.js from its public CDN and requires
an internet connection. Open it in a browser, or from WSL run:

```bash
explorer.exe "$(wslpath -w docking/runs/ligand/visualization_3d.html)"
```

The 3D view shows crystallographic DRD2 helices, the docked pose, nearby
pocket residues, and element-colored pocket and ligand atoms. Its legend
identifies ligand elements. Gold connectors and listed distances mark close
nitrogen/oxygen pairs that could be polar contacts; green connectors mark
close carbon/carbon proximities that could reflect hydrophobic packing. These
are simple geometric cues, not confirmed hydrogen bonds, measured binding
forces, or an energy decomposition. Use the pose as a visual aid and assess
contacts in their chemical and protonation context.

## Scientific validation

PDB 6CM4 is DRD2 with co-crystallized risperidone. A useful workflow control is
to redock the native ligand and compare the predicted pose with its
crystallographic pose. This checks aspects of the docking setup but does not
validate predicted potency for new molecules. Include suitable known active
and inactive reference compounds when evaluating a prioritization protocol,
and assess whether docking adds value beyond the ML model alone. Ultimately,
experimental DRD2 measurements determine whether a prioritized candidate is
active and whether a prediction should inform future model training.
