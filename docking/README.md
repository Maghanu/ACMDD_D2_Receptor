# DRD2 docking

AutoDock Vina's Python bindings are not supported on native Windows. On
Windows, run this workflow in Ubuntu under WSL. Install Miniforge in Ubuntu;
do not try to activate a Python virtual environment stored on `/mnt/c`.

From the Ubuntu terminal:

```bash
cd ~
curl -L -o Miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh
source ~/miniforge3/etc/profile.d/conda.sh
conda create -n drd2-vina -c conda-forge python=3.11 vina meeko rdkit biopython pdbfixer openmm matplotlib -y
conda activate drd2-vina
cd "/mnt/c/Users/manup/Documents/Leiden Universiteit/LUMC/Advanced Computational Methods in Drug Discovery/Mini-Research/ACMDD_D2_Receptor"
python vina_autodock.py
```

Install Vina using conda-forge, not pip: pip may try to compile it from source
and fail if Boost development libraries are unavailable. The conda environment
is located in Ubuntu's home directory; the project files can remain on the
Windows drive.

When prompted, enter the ligand file path. Supported inputs are `.sdf`, `.mol`,
and `.smi` (a SMILES string on the first line). The script downloads PDB 6CM4,
repairs missing receptor atoms with PDBFixer, performs a restrained OpenMM
minimization, prepares receptor and ligand PDBQT files, and centers the search
box on co-crystallized risperidone. It also corrects anomalous terminal-oxygen
coordinates before Meeko validates the receptor.

For salt inputs with one carbon-containing component and separate non-carbon
ions, the docking workflow ignores those ions and reports which ligand
component it uses. Inputs with multiple carbon-containing components are
rejected rather than choosing a ligand automatically.

Outputs are saved under `docking/runs/<ligand-file-name>/`. Docking scores are
computational estimates, not experimental binding affinities. Inspect poses
before interpreting them. PDB 1IEP is Abl kinase, not DRD2.

The minimization holds crystallographic heavy atoms close to their input
positions while allowing rebuilt atoms to relax and resolve severe clashes.

To visualize the best-ranked pose, nearby receptor atoms, polar contacts, and
the scores for all poses, run:

```bash
python docking/plot_docking.py
```

The plot is saved as `docking/runs/ligand/visualization.png`. Use
`--run-dir docking/runs/<ligand-file-name>` to plot another docking run.

For an interactive 3D view with crystallographic DRD2 helices, the docked
ligand, pocket residues, and close polar contacts, run:

```bash
python docking/visualize_docking_3d.py
```

Open `docking/runs/ligand/visualization_3d.html` in a browser with internet
access (the viewer loads 3Dmol.js from its public CDN). Drag to rotate, scroll
to zoom, and use **Focus pocket** to recenter the view.
