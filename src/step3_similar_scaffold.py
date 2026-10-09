from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, rdFMCS

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"

molecules = pd.read_csv(DATA_DIR / "DopamineD2_compounds.csv")[
    ["CHEMBLid", "smiles", "pIC50"]
].rename(columns={"CHEMBLid": "CHEMBL id"})
print(molecules.head(5).to_string(index=False))

mols = [Chem.MolFromSmiles(smiles) for smiles in molecules["smiles"]]
invalid_molecules = sum(mol is None for mol in mols)
if invalid_molecules:
    raise ValueError(f"Could not parse SMILES for {invalid_molecules} molecules")
print(f"Set with {len(mols)} molecules loaded.")

molecules["molecule_weight"] = molecules["smiles"].apply(
    lambda smiles: Descriptors.MolWt(Chem.MolFromSmiles(smiles))
)
molecules.sort_values(["molecule_weight"], ascending=False, inplace=True)

mcs1 = rdFMCS.FindMCS(mols)
print(f"MCS1 contains {mcs1.numAtoms} atoms and {mcs1.numBonds} bonds.")
print("MCS SMARTS string:", mcs1.smartsString)

mcs2 = rdFMCS.FindMCS(mols, threshold=0.8)
print(f"MCS2 contains {mcs2.numAtoms} atoms and {mcs2.numBonds} bonds.")
print("SMARTS string:", mcs2.smartsString)

m2 = Chem.MolFromSmarts(mcs2.smartsString)
if m2 is not None:
    Draw.MolsToGridImage([m2], legends=["MCS2: +threshold=0.8"])

mcs3 = rdFMCS.FindMCS(mols, threshold=0.8, ringMatchesRingOnly=True)
print(f"MCS3 contains {mcs3.numAtoms} atoms and {mcs3.numBonds} bonds.")
print("SMARTS string:", mcs3.smartsString)

m3 = Chem.MolFromSmarts(mcs3.smartsString)
if m3 is not None:
    Draw.MolsToGridImage([m3], legends=["MCS3: +ringmatch"])

reference_path = ROOT / "data" / "EGFR_compounds.csv"
if reference_path.exists():
    mol_df = pd.read_csv(reference_path, index_col=0)
    print("Total number of compounds:", mol_df.shape[0])
    mol_df = mol_df[mol_df.pIC50 > 9]
    print("Number of compounds with pIC50 > 9:", mol_df.shape[0])
else:
    print("Optional external reference data is not present; skipping comparison step.")
