"""Combine the DRD2 IC50 and Kd datasets using SMILES as the compound key.

The output is an outer merge: compounds present in only one source are kept,
and Kd values are reported in nM.
"""

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
IC50_PATH = DATA_DIR / "DopamineD2_compounds.csv"
KD_PATH = DATA_DIR / "DopamineD2_compounds_Kd.csv"
OUTPUT_PATH = DATA_DIR / "DopamineD2_combined.csv"

KD_TO_NM = {
    "M": 1_000_000_000,
    "mM": 1_000_000,
    "uM": 1_000,
    "µM": 1_000,
    "nM": 1,
    "pM": 0.001,
}


def read_unique_rows(path: Path, required_columns: set[str]) -> dict[str, dict[str, str]]:
    """Read CSV records indexed by SMILES and reject duplicate compounds."""
    with path.open("r", newline="", encoding="utf-8-sig") as csv_file:
        reader = csv.DictReader(csv_file)
        available_columns = set(reader.fieldnames or ())
        missing_columns = required_columns - available_columns
        if missing_columns:
            raise ValueError(
                f"{path} is missing required columns: {', '.join(sorted(missing_columns))}"
            )

        rows: dict[str, dict[str, str]] = {}
        for line_number, row in enumerate(reader, start=2):
            smiles = (row.get("smiles") or "").strip()
            if not smiles:
                raise ValueError(f"{path}, line {line_number}: SMILES is empty")
            if smiles in rows:
                raise ValueError(
                    f"{path}, line {line_number}: duplicate SMILES; "
                    "aggregate repeated measurements before combining"
                )
            rows[smiles] = row
    return rows


def main() -> None:
    ic50_rows = read_unique_rows(IC50_PATH, {"smiles", "pIC50"})
    kd_rows = read_unique_rows(KD_PATH, {"smiles", "Kd", "units"})

    combined_smiles = sorted(ic50_rows.keys() | kd_rows.keys())
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(
            csv_file,
            fieldnames=["smiles", "pIC50", "Kd_nM"],
        )
        writer.writeheader()

        for smiles in combined_smiles:
            ic50 = ic50_rows.get(smiles)
            kd = kd_rows.get(smiles)
            kd_nM = ""
            if kd is not None:
                unit = (kd.get("units") or "").strip()
                if unit not in KD_TO_NM:
                    raise ValueError(
                        f"Unsupported Kd unit {unit!r} for SMILES {smiles}"
                    )
                kd_nM = float(kd["Kd"]) * KD_TO_NM[unit]

            writer.writerow(
                {
                    "smiles": smiles,
                    "pIC50": ic50["pIC50"] if ic50 is not None else "",
                    "Kd_nM": kd_nM,
                }
            )

    print(
        f"Wrote {len(combined_smiles)} compounds to {OUTPUT_PATH} "
        f"({len(set(ic50_rows) & set(kd_rows))} present in both datasets)."
    )


if __name__ == "__main__":
    main()
