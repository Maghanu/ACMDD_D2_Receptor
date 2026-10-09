# DRD2 docking for candidate prioritization

## How docking fits the project

The primary activity labels in this project are experimental DRD2 measurements
such as IC50-derived pIC50. Ligand-based ML models are trained on those known
experimental molecules. The longer-term goal is to apply a validated model to
candidate molecules represented as SMILES, prioritize candidates for further
investigation, and use docking as a structural follow-up before in vitro
testing.

Docking does not establish experimental binding or target selectivity. Vina
scores are approximate scoring-function outputs, not measured affinities, and
must not be used as substitutes for experimental training labels. The example
SMILES in `docking/examples/drd2/` are workflow inputs; their presence here
does not assert that they are experimentally confirmed DRD2 ligands.

The current docking script accepts one `.smi`, `.sdf`, or `.mol` ligand per
run. Model-to-SMILES-library screening and automated candidate ranking are
future integration steps; they are not currently performed by the docking
script.

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

## Dock one experimental or candidate molecule

Provide a SMILES string in a `.smi` file (one molecule per file for this
workflow), or an `.sdf`/`.mol` file. Prefer correctly standardized structures
and appropriate protonation/tautomer states for the assay context. For example:

```bash
python vina_autodock.py
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
redock all candidates being compared. Vina's search can still vary, so consider
repeated runs or a more systematic protocol when ranking close candidates.

Salt inputs with exactly one carbon-containing component have disconnected
non-carbon ions removed with a warning. Multiple carbon-containing components
are rejected rather than guessed. Inspect protonation, stereochemistry, and
the resulting pose.

Outputs are stored under `docking/runs/<ligand-file-name>/`, including
`docked_poses.pdbqt`, `vina.log`, and prepared receptor/ligand files. Treat the
scores as a computational prioritization signal only. Check pose plausibility,
compare against appropriate known ligands and controls, and experimentally
test promising candidates in vitro.

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
pocket residues, and close N/O contacts. These distance markers are not a
complete hydrogen-bond analysis.

## Scientific validation

PDB 6CM4 is DRD2 with co-crystallized risperidone. A useful workflow control is
to redock the native ligand and compare the predicted pose with its
crystallographic pose. This checks aspects of the docking setup but does not
validate predicted potency for new molecules. Ultimately, experimental
measurements determine whether a prioritized candidate is active.
