"""Fingerprint external SMILES and predict pIC50 with the trained Morgan2 SVR."""

from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from tqdm.auto import tqdm


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "random_smiles.csv"
DEFAULT_MODEL = ROOT / "data" / "final_svr_morgan2.joblib"
DEFAULT_FEATURES = ROOT / "data" / "random_smiles_morgan2.csv"
DEFAULT_PREDICTIONS = ROOT / "results" / "random_smiles_predicted_pic50.csv"
N_BITS = 2048
MORGAN_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=N_BITS,
)


def fingerprint_smiles(
    smiles: str,
) -> tuple[str | None, float | None, np.ndarray | None]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None, None, None
    canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
    molecular_weight = Descriptors.MolWt(molecule)
    fingerprint = MORGAN_GENERATOR.GetFingerprintAsNumPy(molecule)
    return canonical_smiles, molecular_weight, fingerprint


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Create a Morgan radius-2 feature map and predict pIC50 for "
            "SMILES using the saved experimental SVR model."
        )
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--features-output", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--predictions-output", type=Path, default=DEFAULT_PREDICTIONS)
    args = parser.parse_args()

    input_path = args.input.expanduser().resolve()
    model_path = args.model.expanduser().resolve()
    if not input_path.is_file():
        parser.error(f"Input CSV not found: {input_path}")
    if not model_path.is_file():
        parser.error(f"Trained model not found: {model_path}")

    source = pd.read_csv(input_path)
    if "SMILES" not in source.columns:
        parser.error(f"{input_path} must contain a 'SMILES' column.")

    feature_rows: list[dict] = []
    predictions = np.full(len(source), np.nan, dtype=float)
    valid_rows: list[int] = []
    fingerprints: list[np.ndarray] = []
    canonical_smiles: list[str | None] = []
    molecular_weights: list[float | None] = []

    for row_index, smiles_value in enumerate(
        tqdm(source["SMILES"], desc="Morgan2 fingerprints")
    ):
        smiles = "" if pd.isna(smiles_value) else str(smiles_value).strip()
        canonical, molecular_weight, fingerprint = fingerprint_smiles(smiles)
        canonical_smiles.append(canonical)
        molecular_weights.append(molecular_weight)
        if fingerprint is None:
            continue

        valid_rows.append(row_index)
        fingerprints.append(fingerprint)
        feature_row = source.iloc[row_index].to_dict()
        feature_row["canonical_smiles"] = canonical
        feature_row["molecular_weight"] = molecular_weight
        feature_row.update(
            {
                f"bit_{bit_index}": int(bit_value)
                for bit_index, bit_value in enumerate(fingerprint)
            }
        )
        feature_rows.append(feature_row)

    if not fingerprints:
        parser.error("No valid SMILES were found; no feature map or predictions written.")

    X = np.stack(fingerprints).astype(np.float32, copy=False)
    model = joblib.load(model_path)
    expected_bits = getattr(model, "n_features_in_", None)
    if expected_bits != N_BITS:
        parser.error(
            f"Model expects {expected_bits} features; Morgan2 produces {N_BITS}."
        )
    predictions[valid_rows] = model.predict(X)

    features_output = args.features_output.expanduser().resolve()
    predictions_output = args.predictions_output.expanduser().resolve()
    features_output.parent.mkdir(parents=True, exist_ok=True)
    predictions_output.parent.mkdir(parents=True, exist_ok=True)

    feature_columns = [
        *source.columns,
        "canonical_smiles",
        "molecular_weight",
        *(f"bit_{bit_index}" for bit_index in range(N_BITS)),
    ]
    pd.DataFrame(feature_rows, columns=feature_columns).to_csv(
        features_output,
        index=False,
    )

    prediction_table = source.copy()
    prediction_table["canonical_smiles"] = canonical_smiles
    prediction_table["molecular_weight"] = molecular_weights
    prediction_table["smiles_valid"] = [smiles is not None for smiles in canonical_smiles]
    prediction_table["predicted_pIC50"] = predictions
    prediction_table["predicted_IC50_nM"] = np.where(
        np.isfinite(predictions),
        10 ** (9 - predictions),
        np.nan,
    )
    prediction_table.sort_values(
        "predicted_pIC50",
        ascending=False,
        na_position="last",
    ).to_csv(predictions_output, index=False)

    print(f"Input rows: {len(source)}")
    print(f"Valid SMILES scored: {len(valid_rows)}")
    print(f"Invalid SMILES retained without predictions: {len(source) - len(valid_rows)}")
    print(f"Morgan2 feature map: {features_output}")
    print(f"SVR pIC50 predictions: {predictions_output}")
    print(
        "Predictions use the existing Morgan2 SVR trained on experimental DRD2 "
        "pIC50 data; random_smiles.csv was not used for model training."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
