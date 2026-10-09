from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    Draw,
    Descriptors,
    MACCSkeys,
    rdFingerprintGenerator,
)
from pathlib import Path
from copy import deepcopy
import random

from ipywidgets import interact, fixed, widgets
import pandas as pd
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS


HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

molecules = pd.read_csv(DATA / "DopamineD2_compounds.csv")[
    ["CHEMBLid", "smiles", "pIC50"]
].rename(columns={"CHEMBLid": "CHEMBL id"})
# Show the first 5 molecules
print(molecules.head(5).to_string(index=False))

mols = [Chem.MolFromSmiles(smiles) for smiles in molecules["smiles"]]

invalid_molecules = sum(mol is None for mol in mols)
if invalid_molecules:
    raise ValueError(f"Could not parse SMILES for {invalid_molecules} molecules")
print(f"Set with {len(mols)} molecules loaded.")

# Note -- we use pandas apply function to apply the MolWt function
molecules["molecule_weight"] = molecules["smiles"].apply(
    lambda smiles: Descriptors.MolWt(Chem.MolFromSmiles(smiles))
)
# Sort molecules by molecular weight
molecules.sort_values(["molecule_weight"], ascending=False, inplace=True)


# Find the maximum common substructure (MCS) of the molecules
mcs1 = rdFMCS.FindMCS(mols)
print(f"MCS1 contains {mcs1.numAtoms} atoms and {mcs1.numBonds} bonds.")
print("MCS SMARTS string:", mcs1.smartsString)
# NBVAL_CHECK_OUTPUT

mcs2 = rdFMCS.FindMCS(mols, threshold=0.8)
print(f"MCS2 contains {mcs2.numAtoms} atoms and {mcs2.numBonds} bonds.")
print("SMARTS string:", mcs2.smartsString)
# NBVAL_CHECK_OUTPUT

# Draw substructure
m2 = Chem.MolFromSmarts(mcs2.smartsString)
Draw.MolsToGridImage([m1, m2], legends=["MCS1", "MCS2: +threshold=0.8"])

mcs3 = rdFMCS.FindMCS(mols, threshold=0.8, ringMatchesRingOnly=True)
print(f"MCS3 contains {mcs3.numAtoms} atoms and {mcs3.numBonds} bonds.")
print("SMARTS string:", mcs3.smartsString)
# NBVAL_CHECK_OUTPUT

# Draw substructure
m3 = Chem.MolFromSmarts(mcs3.smartsString)
Draw.MolsToGridImage([m1, m2, m3], legends=["MCS1", "MCS2: +treshold=0.8", "mcs3: +ringmatch"])



# Find large database with potebtial new D2 receptor ligands
mol_df = pd.read_csv(HERE / "../T001_query_chembl/data/EGFR_compounds.csv", index_col=0)
print("Total number of compounds:", mol_df.shape[0])

# Only keep molecules with pIC50 > 9 (IC50 > 1nM)
mol_df = mol_df[mol_df.pIC50 > 9]
print("Number of compounds with pIC50 > 9:", mol_df.shape[0])
# NBVAL_CHECK_OUTPUT