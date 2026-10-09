"""Encode SMILES from DopamineD2_compounds.csv as fingerprints
(MACCS, Morgan radius 2, Morgan radius 3), following TeachOpenCADD T007."""

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdFingerprintGenerator
from tqdm.auto import tqdm

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
CSV_PATH = DATA_DIR / "DopamineD2_compounds.csv"
ACTIVITY_CUTOFF = 6.3  # pIC50 cutoff used in T007 (optional labeling)
N_BITS = 2048

# Generators are created once and reused
_MORGAN_GENERATORS = {
    "morgan2": rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=N_BITS),
    "morgan3": rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=N_BITS),
}


def smiles_to_fp(smiles, method="maccs"):
    """Encode a SMILES string as a numpy fingerprint array.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    method : str
        "maccs" (166 bits), "morgan2" or "morgan3" (2048 bits).

    Returns
    -------
    np.ndarray or None
        Fingerprint, or None if the SMILES cannot be parsed.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if method == "maccs":
        return np.array(MACCSkeys.GenMACCSKeys(mol))
    if method in _MORGAN_GENERATORS:
        return _MORGAN_GENERATORS[method].GetFingerprintAsNumPy(mol)
    raise ValueError(f"Unknown method: {method}")


def drop_invalid_smiles(df):
    """Remove rows whose SMILES RDKit cannot parse."""
    valid = df["smiles"].apply(lambda smi: Chem.MolFromSmiles(smi) is not None)
    if not valid.all():
        print(f"Dropped {(~valid).sum()} molecules with invalid SMILES")
    return df[valid].reset_index(drop=True)


def fingerprint_matrix(smiles_series, method):
    """Return an (n_molecules, n_bits) array for the given fingerprint method."""
    fps = [smiles_to_fp(smi, method) for smi in tqdm(smiles_series, desc=method)]
    return np.stack(fps)


if __name__ == "__main__":
    df = pd.read_csv(CSV_PATH)[["smiles", "pIC50"]]
    print("Shape:", df.shape)

    # Activity label (as in T007)
    df["active"] = (df["pIC50"] >= ACTIVITY_CUTOFF).astype(float)
    df = drop_invalid_smiles(df)

    # One CSV per fingerprint: bit columns, then smiles, pIC50, active
    for method in ("maccs", "morgan2", "morgan3"):
        X = fingerprint_matrix(df["smiles"], method)
        first_bit = 0
        if method == "maccs":
            X = X[:, 1:]  # RDKit MACCS bit 0 is unused (always 0)
            first_bit = 1  # keep original key numbering: bit_1 ... bit_166
        out = pd.DataFrame(
            X, columns=[f"bit_{i}" for i in range(first_bit, first_bit + X.shape[1])]
        )
        out["smiles"] = df["smiles"]
        out["pIC50"] = df["pIC50"]
        out["active"] = df["active"]
        out_path = DATA_DIR / f"DopamineD2_{method}.csv"
        out.to_csv(out_path, index=False)
        print(f"{method}: {X.shape} -> {out_path}")
