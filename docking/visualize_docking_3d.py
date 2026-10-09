"""Create an interactive 3D view of the best Vina pose in its DRD2 pocket."""

from __future__ import annotations

import argparse
import json
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for import_path in (ROOT, Path(__file__).resolve().parent):
    if str(import_path) not in sys.path:
        sys.path.insert(0, str(import_path))

from plot_docking import (
    DEFAULT_RUN_DIR,
    atom_positions_by_molecule_index,
    parse_pose_models,
)


def parse_mmcif_loop(path: Path, category: str) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    prefix = f"_{category}."
    for line_index, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        headers = []
        cursor = line_index + 1
        while cursor < len(lines) and lines[cursor].strip().startswith(prefix):
            headers.append(lines[cursor].strip())
            cursor += 1
        if not headers:
            continue
        rows = []
        while cursor < len(lines):
            row_text = lines[cursor].strip()
            if not row_text:
                cursor += 1
                continue
            if row_text.startswith("#"):
                break
            if row_text.startswith(("loop_", "data_", "save_", "_")):
                break
            values = shlex.split(row_text)
            if len(values) != len(headers):
                raise ValueError(f"Could not parse {category} records in {path}.")
            rows.append(values)
            cursor += 1
        return headers, rows
    raise ValueError(f"The {category} records were not found in {path}.")


def pdb_helix_record(
    serial: int,
    helix_id: str,
    begin: tuple[str, str, str, str],
    end: tuple[str, str, str, str],
    helix_class: int,
    length: int,
) -> str:
    record = [" "] * 80
    record[0:6] = "HELIX "
    record[7:10] = f"{serial:3d}"
    record[11:14] = f"{helix_id:>3.3s}"
    record[15:18] = f"{begin[0]:>3.3s}"
    record[19] = begin[1]
    record[21:25] = f"{int(begin[2]):4d}"
    record[25] = " " if begin[3] in {".", "?"} else begin[3]
    record[27:30] = f"{end[0]:>3.3s}"
    record[31] = end[1]
    record[33:37] = f"{int(end[2]):4d}"
    record[37] = " " if end[3] in {".", "?"} else end[3]
    record[38:40] = f"{helix_class:2d}"
    record[71:76] = f"{length:5d}"
    return "".join(record) + "\n"


def receptor_pdb_with_helices(run_dir: Path) -> str:
    from vina_autodock import DRD2_SEQUENCE_RANGES

    cif_path = run_dir / "6CM4.cif"
    receptor_path = run_dir / "drd2_repaired.pdb"
    atom_headers, atom_rows = parse_mmcif_loop(cif_path, "atom_site")
    atom_index = {header: index for index, header in enumerate(atom_headers)}
    seq_to_auth = {}
    for row in atom_rows:
        if row[atom_index["_atom_site.group_PDB"]] != "ATOM":
            continue
        if row[atom_index["_atom_site.label_entity_id"]] != "1":
            continue
        if row[atom_index["_atom_site.label_asym_id"]] != "A":
            continue
        sequence = row[atom_index["_atom_site.label_seq_id"]]
        if sequence in {".", "?"}:
            continue
        seq_to_auth.setdefault(
            int(sequence),
            (
                row[atom_index["_atom_site.label_comp_id"]],
                row[atom_index["_atom_site.auth_seq_id"]],
                row[atom_index["_atom_site.pdbx_PDB_ins_code"]],
            ),
        )

    helix_headers, helix_rows = parse_mmcif_loop(cif_path, "struct_conf")
    helix_index = {header: index for index, header in enumerate(helix_headers)}
    helix_records = []
    serial = 1
    for row in helix_rows:
        if row[helix_index["_struct_conf.conf_type_id"]] not in {"HELX_P", "HELX_R"}:
            continue
        if row[helix_index["_struct_conf.beg_label_asym_id"]] != "A":
            continue
        first = int(row[helix_index["_struct_conf.beg_label_seq_id"]])
        last = int(row[helix_index["_struct_conf.end_label_seq_id"]])
        for segment_index, (segment_start, segment_end) in enumerate(DRD2_SEQUENCE_RANGES):
            clipped_start = max(first, segment_start)
            clipped_end = min(last, segment_end)
            if clipped_start > clipped_end:
                continue
            begin_info = seq_to_auth.get(clipped_start)
            end_info = seq_to_auth.get(clipped_end)
            if begin_info is None or end_info is None:
                continue
            output_chain = "A" if segment_index == 0 else "B"
            helix_id = row[helix_index["_struct_conf.pdbx_PDB_helix_id"]]
            helix_class = row[helix_index["_struct_conf.pdbx_PDB_helix_class"]]
            helix_records.append(
                pdb_helix_record(
                    serial,
                    helix_id,
                    (begin_info[0], output_chain, begin_info[1], begin_info[2]),
                    (end_info[0], output_chain, end_info[1], end_info[2]),
                    int(helix_class) if helix_class not in {".", "?"} else 1,
                    clipped_end - clipped_start + 1,
                )
            )
            serial += 1

    if not helix_records:
        raise ValueError(f"No DRD2 helix annotations were found in {cif_path}.")

    atom_lines = []
    for line in receptor_path.read_text(encoding="ascii").splitlines():
        if not line.startswith("ATOM  "):
            continue
        if line[76:78].strip().upper() in {"H", "D"}:
            continue
        atom_lines.append(line[:80].ljust(80) + "\n")
    if not atom_lines:
        raise ValueError(f"No receptor heavy atoms were found in {receptor_path}.")

    return "".join(helix_records + atom_lines + ["END\n"])


def receptor_atoms(path: Path) -> list[dict]:
    atoms = []
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.startswith("ATOM  ") or line[76:78].strip().upper() in {"H", "D"}:
            continue
        atoms.append(
            {
                "chain": line[21].strip(),
                "resi": int(line[22:26]),
                "residue": line[17:20].strip(),
                "atom": line[12:16].strip(),
                "element": line[76:78].strip().upper(),
                "position": tuple(
                    float(line[start:start + 8]) for start in (30, 38, 46)
                ),
            }
        )
    return atoms


def make_ligand_molblock(run_dir: Path, best_pose: dict) -> tuple[str, list[dict]]:
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    supplier = Chem.SDMolSupplier(
        str(run_dir / "ligand_3d.sdf"),
        removeHs=False,
    )
    source = next((molecule for molecule in supplier if molecule is not None), None)
    if source is None:
        raise ValueError(f"Could not read a ligand from {run_dir / 'ligand_3d.sdf'}.")
    positions = atom_positions_by_molecule_index(best_pose, source)

    ligand_atoms = []
    for atom_index, position in positions.items():
        atom = source.GetAtomWithIdx(atom_index)
        if atom.GetSymbol() in {"H", "D"}:
            continue
        ligand_atoms.append(
            {
                "element": atom.GetSymbol().upper(),
                "index": atom_index,
                "position": position,
            }
        )
    conformer = source.GetConformer()
    for atom_index, position in positions.items():
        conformer.SetAtomPosition(atom_index, Point3D(*position))
    molecule = Chem.RemoveHs(source)
    return Chem.MolToMolBlock(molecule), ligand_atoms


def contact_data(
    protein_atoms: list[dict],
    ligand_atoms: list[dict],
) -> tuple[list[dict], list[dict]]:
    contacts_by_residue = {}
    polar_pairs = []
    for protein_atom in protein_atoms:
        closest = min(
            (
                (
                    sum(
                        (
                            protein_atom["position"][axis]
                            - ligand_atom["position"][axis]
                        ) ** 2
                        for axis in range(3)
                    ) ** 0.5,
                    ligand_atom,
                )
                for ligand_atom in ligand_atoms
            ),
            key=lambda item: item[0],
        )
        distance, ligand_atom = closest
        residue_key = (protein_atom["chain"], protein_atom["resi"])
        current = contacts_by_residue.get(residue_key)
        if distance <= 4.5 and (current is None or distance < current["distance"]):
            contacts_by_residue[residue_key] = {
                "chain": protein_atom["chain"],
                "resi": protein_atom["resi"],
                "residue": protein_atom["residue"],
                "distance": distance,
                "position": protein_atom["position"],
            }
        if (
            distance <= 3.5
            and protein_atom["element"] in {"N", "O"}
            and ligand_atom["element"] in {"N", "O"}
        ):
            polar_pairs.append(
                {
                    "distance": distance,
                    "start": protein_atom["position"],
                    "end": ligand_atom["position"],
                    "residue": protein_atom["residue"],
                    "chain": protein_atom["chain"],
                    "resi": protein_atom["resi"],
                    "ligand_element": ligand_atom["element"],
                }
            )
    contacts = sorted(contacts_by_residue.values(), key=lambda item: item["distance"])
    polar_pairs.sort(key=lambda item: item["distance"])
    return contacts, polar_pairs[:6]


def build_html(run_dir: Path) -> str:
    poses = parse_pose_models(run_dir / "docked_poses.pdbqt")
    best_pose = poses[0]
    protein_pdb = receptor_pdb_with_helices(run_dir)
    atoms = receptor_atoms(run_dir / "drd2_repaired.pdb")
    ligand_molblock, ligand_atoms = make_ligand_molblock(run_dir, best_pose)
    contacts, polar_pairs = contact_data(atoms, ligand_atoms)
    pocket = [contact for contact in contacts if contact["distance"] <= 4.5]
    pocket_groups = {}
    for contact in pocket:
        pocket_groups.setdefault(contact["chain"], []).append(contact["resi"])
    labels = [
        {
            "text": f"{contact['residue']}{contact['resi']} ({contact['distance']:.1f} Å)",
            "chain": contact["chain"],
            "resi": contact["resi"],
        }
        for contact in contacts[:10]
    ]
    payload = json.dumps(
        {
            "protein": protein_pdb,
            "ligand": ligand_molblock,
            "score": best_pose["score"],
            "pocket": pocket_groups,
            "labels": labels,
            "polarPairs": polar_pairs,
            "poses": len(poses),
            "residueCount": len(contacts),
        }
    )
    viewer_script = r"""
const data = __DATA__;
const element = document.getElementById('viewer');
const viewer = $3Dmol.createViewer(element, {
  backgroundColor: '#101820',
  antialias: true
});
viewer.setBackgroundColor('#101820');
const protein = viewer.addModel(data.protein, 'pdb');
protein.setStyle({}, {});
viewer.addStyle(
  {model: protein, ss: 'h'},
  {cartoon: {color: '#5793b3', opacity: 0.48}}
);
for (const [chain, residues] of Object.entries(data.pocket)) {
  viewer.addStyle(
    {model: protein, chain: chain, resi: residues},
    {stick: {radius: 0.19, color: '#8fa3af'}}
  );
}
for (const label of data.labels) {
  viewer.addResLabels(
    {model: protein, chain: label.chain, resi: label.resi},
    {
      font: 'Arial',
      fontSize: 9,
      fontColor: '#f4f7fa',
      backgroundColor: 'rgba(16,24,32,0.82)',
      showBackground: true
    }
  );
}
const ligand = viewer.addModel(data.ligand, 'mol');
ligand.setStyle({}, {
  stick: {radius: 0.27, color: '#f0529c'},
  sphere: {scale: 0.27, color: '#f0529c'}
});
for (const [element, color] of Object.entries({
  N: '#2789e8',
  O: '#ed5149',
  S: '#d99b19',
  CL: '#31ad68'
})) {
  viewer.setStyle(
    {model: ligand, elem: element},
    {
      stick: {radius: 0.27, color: color},
      sphere: {scale: 0.27, color: color}
    },
    true
  );
}
for (const pair of data.polarPairs) {
  viewer.addCylinder({
    start: {x: pair.start[0], y: pair.start[1], z: pair.start[2]},
    end: {x: pair.end[0], y: pair.end[1], z: pair.end[2]},
    radius: 0.045,
    color: '#f0c75e',
    fromCap: 1,
    toCap: 1
  });
}
const focusSelections = [{model: ligand}];
for (const [chain, residues] of Object.entries(data.pocket)) {
  focusSelections.push({model: protein, chain: chain, resi: residues});
}
viewer.zoomTo({or: focusSelections});
viewer.render();
document.getElementById('score').textContent =
  `Best Vina score ${data.score.toFixed(3)} kcal/mol · ${data.poses} poses`;
document.getElementById('contacts').textContent =
  `${data.residueCount} residues within 4.5 Å of the ligand`;
document.getElementById('reset').addEventListener('click', () => {
  viewer.zoomTo({or: focusSelections});
  viewer.render();
});
"""
    viewer_script = viewer_script.replace("__DATA__", payload)
    legend_items = "".join(
        f"<li><span class='swatch {css_class}'></span>{text}</li>"
        for css_class, text in (
            ("helix", "DRD2 helix / backbone cartoon"),
            ("pocket", "Pocket residues within 4.5 Å"),
            ("ligand", "Docked ligand"),
            ("polar", "Close N/O contacts"),
        )
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>DRD2 docking pose · interactive 3D</title>
  <script src="https://3dmol.org/build/3Dmol-min.js"></script>
  <style>
    :root {{ color-scheme: dark; font-family: Inter, Segoe UI, sans-serif; }}
    body {{ margin: 0; background: #101820; color: #edf2f5; }}
    header {{ padding: 18px 24px 12px; display: flex; justify-content: space-between; gap: 16px; align-items: center; }}
    h1 {{ margin: 0 0 5px; font-size: 20px; }}
    .sub {{ color: #b7c4ce; font-size: 13px; }}
    #score {{ color: #f3d27a; font-weight: 700; }}
    main {{ display: grid; grid-template-columns: minmax(0, 1fr) 245px; gap: 12px; padding: 0 16px 16px; }}
    #viewer {{ height: min(78vh, 820px); min-height: 560px; border: 1px solid #40515d; border-radius: 10px; overflow: hidden; }}
    aside {{ background: #17232d; border: 1px solid #34444f; border-radius: 10px; padding: 16px; }}
    aside h2 {{ font-size: 14px; margin: 0 0 12px; }}
    ul {{ list-style: none; padding: 0; margin: 0 0 20px; }}
    li {{ display: flex; gap: 9px; align-items: center; margin: 12px 0; font-size: 13px; color: #d4dde3; }}
    .swatch {{ width: 13px; height: 13px; border-radius: 50%; flex: 0 0 auto; }}
    .helix {{ background: #6fa8dc; }} .pocket {{ background: #aebbc3; }}
    .ligand {{ background: #f0529c; }} .polar {{ background: #f0c75e; }}
    button {{ background: #263b49; color: white; border: 1px solid #537082; border-radius: 6px; padding: 8px 11px; cursor: pointer; }}
    .hint {{ color: #9daeb9; font-size: 12px; line-height: 1.5; margin-top: 14px; }}
    @media (max-width: 820px) {{ main {{ grid-template-columns: 1fr; }} #viewer {{ min-height: 430px; }} }}
  </style>
</head>
<body>
  <header>
    <div><h1>DRD2 · docked ligand in the binding pocket</h1>
    <div class="sub" id="score"></div><div class="sub" id="contacts"></div></div>
    <button id="reset" type="button">Focus pocket</button>
  </header>
  <main>
    <div id="viewer"></div>
    <aside>
      <h2>What you are seeing</h2>
      <ul>{legend_items}</ul>
      <div class="hint">The crystallographic helix annotations are used for the protein cartoon. Pocket-side chains are shown as sticks; the ligand is shown in ball-and-stick style.</div>
      <div class="hint">Drag to rotate · scroll to zoom · right-drag to pan. Gold connectors mark close ligand/receptor nitrogen–oxygen contacts; they are distance-based, not a complete hydrogen-bond analysis.</div>
      <div class="hint">Interactive viewer loads 3Dmol.js from 3dmol.org; an internet connection is required.</div>
    </aside>
  </main>
  <script>{viewer_script}</script>
</body>
</html>
"""


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create an interactive 3D view of the best-ranked Vina pose."
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=DEFAULT_RUN_DIR,
        help=f"folder containing docking outputs (default: {DEFAULT_RUN_DIR})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="output HTML path (default: <run-dir>/visualization_3d.html)",
    )
    args = parser.parse_args()
    run_dir = args.run_dir.expanduser().resolve()
    output = (args.output or run_dir / "visualization_3d.html").expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(build_html(run_dir), encoding="utf-8")
    print(f"Saved interactive 3D visualization to: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
