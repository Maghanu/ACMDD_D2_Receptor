"""Plot a Vina docking pose in its receptor pocket and compare pose scores."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNS_DIR = ROOT / "docking" / "runs"
DEFAULT_RUN_DIR = RUNS_DIR / "ligand"
ELEMENT_COLORS = {
    "C": "#59636e",
    "N": "#3977b8",
    "O": "#d34b4b",
    "S": "#d49a24",
    "P": "#e17a28",
    "CL": "#30965c",
    "F": "#55a868",
}
ATOM_TYPE_ELEMENTS = {
    "A": "C",
    "C": "C",
    "CL": "Cl",
    "F": "F",
    "N": "N",
    "NA": "N",
    "OA": "O",
    "O": "O",
    "P": "P",
    "S": "S",
    "SA": "S",
}


def parse_pose_models(path: Path) -> list[dict]:
    models = []
    current = None
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith("MODEL"):
            current = {"score": None, "atoms": [], "smiles_indices": {}}
        elif current is None:
            continue
        elif line.startswith("REMARK VINA RESULT:"):
            current["score"] = float(line.split()[3])
        elif line.startswith("REMARK SMILES IDX"):
            values = [int(value) for value in line.split()[3:]]
            current["smiles_indices"].update(
                {
                    smiles_index - 1: atom_serial - 1
                    for smiles_index, atom_serial in zip(values[::2], values[1::2])
                }
            )
        elif line.startswith(("ATOM  ", "HETATM")):
            atom_type = line.split()[-1].upper()
            current["atoms"].append(
                {
                    "serial": int(line[6:11]) - 1,
                    "element": ATOM_TYPE_ELEMENTS.get(atom_type, atom_type),
                    "position": tuple(
                        float(line[start:start + 8]) for start in (30, 38, 46)
                    ),
                }
            )
        elif line.startswith("ENDMDL"):
            if current["score"] is None:
                raise ValueError(f"Pose in {path} is missing its Vina score.")
            models.append(current)
            current = None

    if not models:
        raise ValueError(f"No scored poses were found in {path}.")
    return models


def parse_receptor(path: Path) -> list[dict]:
    atoms = []
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_type = line.split()[-1].upper()
        element = ATOM_TYPE_ELEMENTS.get(atom_type, atom_type)
        if element in {"H", "D"}:
            continue
        atoms.append(
            {
                "element": element,
                "residue": f"{line[17:20].strip()} {line[21].strip()}{line[22:26].strip()}",
                "position": tuple(
                    float(line[start:start + 8]) for start in (30, 38, 46)
                ),
            }
        )
    if not atoms:
        raise ValueError(f"No receptor atoms were found in {path}.")
    return atoms


def atom_positions_by_molecule_index(model: dict, molecule) -> dict[int, tuple[float, ...]]:
    pose_by_serial = {
        atom["serial"]: atom
        for atom in model["atoms"]
    }
    molecule_positions = {}
    for atom_index, pose_serial in model["smiles_indices"].items():
        if atom_index >= molecule.GetNumAtoms() or pose_serial not in pose_by_serial:
            continue
        molecule_atom = molecule.GetAtomWithIdx(atom_index)
        pose_atom = pose_by_serial[pose_serial]
        if molecule_atom.GetSymbol().upper() != pose_atom["element"].upper():
            raise ValueError(
                "Ligand atom mapping does not match the docked PDBQT pose. "
                "Use ligand_3d.sdf generated for this docking run."
            )
        molecule_positions[atom_index] = pose_atom["position"]

    heavy_atom_indices = [
        atom.GetIdx() for atom in molecule.GetAtoms()
        if atom.GetSymbol() not in {"H", "D"}
    ]
    missing = set(heavy_atom_indices) - molecule_positions.keys()
    if missing:
        raise ValueError(
            f"The docked pose does not map {len(missing)} ligand heavy atom(s) "
            "to ligand_3d.sdf."
        )
    return molecule_positions


def select_projection(positions: list[tuple[float, ...]]) -> tuple[int, int]:
    return max(
        itertools.combinations(range(3), 2),
        key=lambda axes: (
            max(point[axes[0]] for point in positions)
            - min(point[axes[0]] for point in positions)
        )
        * (
            max(point[axes[1]] for point in positions)
            - min(point[axes[1]] for point in positions)
        ),
    )


def make_plot(run_dir: Path, output_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from rdkit import Chem

    poses = parse_pose_models(run_dir / "docked_poses.pdbqt")
    receptor = parse_receptor(run_dir / "drd2.pdbqt")
    supplier = Chem.SDMolSupplier(str(run_dir / "ligand_3d.sdf"), removeHs=False)
    molecule = next((item for item in supplier if item is not None), None)
    if molecule is None:
        raise ValueError(f"Could not read a valid ligand from {run_dir / 'ligand_3d.sdf'}.")

    best_pose = poses[0]
    ligand_positions = atom_positions_by_molecule_index(best_pose, molecule)
    heavy_atom_indices = [
        atom.GetIdx() for atom in molecule.GetAtoms()
        if atom.GetSymbol() not in {"H", "D"}
    ]
    coordinate_by_index = {
        atom_index: np.asarray(ligand_positions[atom_index])
        for atom_index in heavy_atom_indices
    }
    axis_x, axis_y = select_projection(list(coordinate_by_index.values()))
    axes_labels = ("x", "y", "z")

    proximity_cutoff = 4.2
    closest_atom_by_residue = {}
    contact_by_residue = {}
    for receptor_atom in receptor:
        receptor_position = np.asarray(receptor_atom["position"])
        distances = [
            (float(np.linalg.norm(receptor_position - coordinate_by_index[index])), index)
            for index in heavy_atom_indices
        ]
        distance, ligand_index = min(distances)
        residue_key = receptor_atom["residue"]
        previous = contact_by_residue.get(residue_key)
        if previous is None or distance < previous[0]:
            contact_by_residue[residue_key] = (
                distance,
                receptor_position,
                coordinate_by_index[ligand_index],
                receptor_atom["element"],
                molecule.GetAtomWithIdx(ligand_index).GetSymbol(),
            )
        previous_atom = closest_atom_by_residue.get(residue_key)
        if distance <= proximity_cutoff and (
            previous_atom is None or distance < previous_atom[0]
        ):
            closest_atom_by_residue[residue_key] = (distance, receptor_atom)

    figure, (pose_axis, score_axis) = plt.subplots(
        1,
        2,
        figsize=(15, 8),
        gridspec_kw={"width_ratios": (1.6, 1)},
        layout="constrained",
    )
    figure.suptitle(
        f"DRD2 docking pose | Vina score {best_pose['score']:.3f} kcal/mol",
        fontsize=17,
        fontweight="bold",
    )

    for _, receptor_atom in closest_atom_by_residue.values():
        point = receptor_atom["position"]
        element = receptor_atom["element"].upper()
        pose_axis.scatter(
            point[axis_x],
            point[axis_y],
            s=54,
            color=ELEMENT_COLORS.get(element, "#9a9a9a"),
            alpha=0.86,
            edgecolors="white",
            linewidths=0.5,
            zorder=1,
        )

    for bond in molecule.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin not in coordinate_by_index or end not in coordinate_by_index:
            continue
        start_point = coordinate_by_index[begin]
        end_point = coordinate_by_index[end]
        pose_axis.plot(
            [start_point[axis_x], end_point[axis_x]],
            [start_point[axis_y], end_point[axis_y]],
            color="#263238",
            linewidth=2.6,
            zorder=3,
        )

    for element in sorted({molecule.GetAtomWithIdx(index).GetSymbol() for index in heavy_atom_indices}):
        indices = [
            index for index in heavy_atom_indices
            if molecule.GetAtomWithIdx(index).GetSymbol() == element
        ]
        pose_axis.scatter(
            [coordinate_by_index[index][axis_x] for index in indices],
            [coordinate_by_index[index][axis_y] for index in indices],
            s=115 if element != "C" else 64,
            color=ELEMENT_COLORS.get(element.upper(), "#9a9a9a"),
            edgecolors="white",
            linewidths=0.8,
            label=f"Ligand {element}",
            zorder=4,
        )

    polar_contacts = []
    for residue, (distance, receptor_point, ligand_point, receptor_element, ligand_element) in contact_by_residue.items():
        if (
            distance <= 3.5
            and receptor_element in {"N", "O"}
            and ligand_element in {"N", "O"}
        ):
            polar_contacts.append((distance, residue, receptor_point, ligand_point))
    for distance, residue, receptor_point, ligand_point in sorted(polar_contacts)[:6]:
        pose_axis.plot(
            [receptor_point[axis_x], ligand_point[axis_x]],
            [receptor_point[axis_y], ligand_point[axis_y]],
            color="#7b3294",
            linewidth=1.3,
            linestyle="--",
            alpha=0.8,
            zorder=2,
        )
        midpoint_x = (receptor_point[axis_x] + ligand_point[axis_x]) / 2
        midpoint_y = (receptor_point[axis_y] + ligand_point[axis_y]) / 2
        pose_axis.annotate(
            f"{distance:.1f} Å",
            (midpoint_x, midpoint_y),
            fontsize=8,
            color="#67217f",
            xytext=(3, 3),
            textcoords="offset points",
        )

    residue_labels = [
        (residue, values)
        for residue, values in contact_by_residue.items()
        if values[0] <= 4.0
    ]
    for residue, (distance, receptor_point, _, _, _) in sorted(
        residue_labels, key=lambda item: item[1][0]
    )[:10]:
        pose_axis.annotate(
            residue,
            (receptor_point[axis_x], receptor_point[axis_y]),
            fontsize=8,
            color="#39434b",
            xytext=(4, 5),
            textcoords="offset points",
        )

    pose_axis.set_title(
        f"Best pose and receptor atoms within {proximity_cutoff:.0f} Å\n"
        f"Projection: {axes_labels[axis_x]}–{axes_labels[axis_y]} plane",
        fontsize=12,
    )
    pose_axis.set_xlabel(f"{axes_labels[axis_x]} coordinate (Å)")
    pose_axis.set_ylabel(f"{axes_labels[axis_y]} coordinate (Å)")
    pose_axis.set_aspect("equal", adjustable="datalim")
    pose_axis.grid(alpha=0.16)
    pose_axis.legend(loc="upper left", fontsize=8, ncol=2)

    scores = [pose["score"] for pose in poses]
    ranks = list(range(1, len(scores) + 1))
    bars = score_axis.barh(
        ranks,
        scores,
        color=["#d1495b"] + ["#6c91bf"] * (len(scores) - 1),
        edgecolor="white",
    )
    score_axis.invert_yaxis()
    score_axis.set_yticks(ranks, [f"Pose {rank}" for rank in ranks])
    score_axis.set_xlabel("Vina affinity (kcal/mol; more negative is better)")
    score_axis.set_title("Docking scores", fontsize=12)
    score_axis.grid(axis="x", alpha=0.2)
    for bar, score in zip(bars, scores):
        score_axis.text(
            score + 0.08,
            bar.get_y() + bar.get_height() / 2,
            f"{score:.3f}",
            va="center",
            ha="left",
            color="white",
            fontsize=9,
            fontweight="bold",
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved docking visualization to: {output_path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot the top Vina pose and its receptor-pocket contacts."
    )
    parser.add_argument(
        "ligand_name",
        nargs="?",
        help="name of the ligand docking run under docking/runs/",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="output image path (default: docking/runs/<ligand-name>/visualization.png)",
    )
    args = parser.parse_args()
    required_files = (
        "docked_poses.pdbqt",
        "drd2.pdbqt",
        "ligand_3d.sdf",
        "6CM4.cif",
        "drd2_repaired.pdb",
    )
    prompted = args.ligand_name is None
    while True:
        ligand_name = args.ligand_name
        if prompted:
            ligand_name = input(
                "What ligand do you want to visualise? "
            ).strip()
        if not ligand_name:
            reason = "a ligand name is required"
        elif Path(ligand_name).name != ligand_name or ligand_name in {".", ".."}:
            reason = "enter a ligand run name, not a path"
        else:
            run_dir = RUNS_DIR / ligand_name
            missing_files = [
                filename for filename in required_files
                if not (run_dir / filename).is_file()
            ]
            if not missing_files:
                break
            reason = (
                f"no complete docking run found for {ligand_name!r}; "
                f"missing: {', '.join(missing_files)}"
            )

        if not prompted:
            parser.error(reason)
        print(f"Invalid ligand: {reason}. Please try again.")

    output_path = args.output or run_dir / "visualization.png"
    run_dir = run_dir.resolve()
    make_plot(run_dir, output_path.expanduser().resolve())

    from visualize_docking_3d import build_html

    html_path = run_dir / "visualization_3d.html"
    html_path.write_text(build_html(run_dir), encoding="utf-8")
    print(f"Saved interactive 3D visualization to: {html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
