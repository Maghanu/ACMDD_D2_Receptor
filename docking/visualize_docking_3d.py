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
from Workflow_steps.Workflow_step7_vina_autodock import (
    RECEPTOR_STATES,
    get_receptor_state,
)

ELEMENT_COLORS = {
    "C": "#b8c2cc",
    "N": "#4da3ff",
    "O": "#ff5c5c",
    "S": "#f2c14e",
    "P": "#ff9f43",
    "F": "#62c370",
    "CL": "#62c370",
    "BR": "#a879d8",
    "I": "#a879d8",
}


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
    structure_ids = [
        structure_id
        for structure_id in RECEPTOR_STATES
        if (run_dir / f"{structure_id}.cif").is_file()
    ]
    if len(structure_ids) != 1:
        raise ValueError(
            f"Expected exactly one supported receptor CIF in {run_dir}; "
            f"found {structure_ids}."
        )
    receptor_state = get_receptor_state(structure_ids[0])
    cif_path = run_dir / f"{receptor_state.structure_id}.cif"
    receptor_path = run_dir / "drd2_repaired.pdb"
    atom_headers, atom_rows = parse_mmcif_loop(cif_path, "atom_site")
    atom_index = {header: index for index, header in enumerate(atom_headers)}
    seq_to_auth = {}
    entity_chains = set()
    for row in atom_rows:
        if row[atom_index["_atom_site.group_PDB"]] != "ATOM":
            continue
        if row[atom_index["_atom_site.label_entity_id"]] != "1":
            continue
        entity_chains.add(row[atom_index["_atom_site.label_asym_id"]])
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
    if len(entity_chains) != 1:
        raise ValueError(
            f"Expected one DRD2-containing polymer chain in {cif_path}; "
            f"found {sorted(entity_chains)}."
        )
    label_chain = next(iter(entity_chains))

    pdb_chain_by_residue = {}
    for line in receptor_path.read_text(encoding="ascii").splitlines():
        if line.startswith("ATOM  "):
            key = (
                line[17:20].strip(),
                line[22:26].strip(),
                line[26],
            )
            pdb_chain_by_residue[key] = line[21]
    seq_to_output_chain = {}
    for sequence, residue in seq_to_auth.items():
        insertion = " " if residue[2] in {".", "?"} else residue[2]
        key = (residue[0], residue[1].strip(), insertion)
        output_chain = pdb_chain_by_residue.get(key)
        if output_chain is not None:
            seq_to_output_chain[sequence] = output_chain

    helix_headers, helix_rows = parse_mmcif_loop(cif_path, "struct_conf")
    helix_index = {header: index for index, header in enumerate(helix_headers)}
    helix_records = []
    serial = 1
    for row in helix_rows:
        if row[helix_index["_struct_conf.conf_type_id"]] not in {"HELX_P", "HELX_R"}:
            continue
        if row[helix_index["_struct_conf.beg_label_asym_id"]] != label_chain:
            continue
        first = int(row[helix_index["_struct_conf.beg_label_seq_id"]])
        last = int(row[helix_index["_struct_conf.end_label_seq_id"]])
        for segment_start, segment_end in receptor_state.drd2_sequence_ranges:
            clipped_start = max(first, segment_start)
            clipped_end = min(last, segment_end)
            if clipped_start > clipped_end:
                continue
            helix_id = row[helix_index["_struct_conf.pdbx_PDB_helix_id"]]
            helix_class = row[helix_index["_struct_conf.pdbx_PDB_helix_class"]]
            blocks = []
            current_block = []
            for sequence in range(clipped_start, clipped_end + 1):
                residue = seq_to_auth.get(sequence)
                output_chain = seq_to_output_chain.get(sequence)
                if residue is None or output_chain is None:
                    if current_block:
                        blocks.append(current_block)
                        current_block = []
                    continue
                if current_block and (
                    sequence != current_block[-1][0] + 1
                    or output_chain != current_block[-1][2]
                ):
                    blocks.append(current_block)
                    current_block = []
                current_block.append((sequence, residue, output_chain))
            if current_block:
                blocks.append(current_block)

            for block in blocks:
                begin_info = block[0][1]
                end_info = block[-1][1]
                output_chain = block[0][2]
                helix_records.append(
                    pdb_helix_record(
                        serial,
                        helix_id,
                        (begin_info[0], output_chain, begin_info[1], begin_info[2]),
                        (end_info[0], output_chain, end_info[1], end_info[2]),
                        int(helix_class) if helix_class not in {".", "?"} else 1,
                        len(block),
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

    pose_smiles = best_pose.get("smiles")
    if pose_smiles:
        source = Chem.MolFromSmiles(pose_smiles)
        if source is None:
            raise ValueError(
                f"Could not parse the ligand SMILES stored in the docked pose "
                f"for {run_dir}."
            )
    else:
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
                "hydrophobic": atom.GetSymbol() == "C"
                and all(
                    neighbor.GetSymbol() == "C"
                    for neighbor in atom.GetNeighbors()
                    if neighbor.GetSymbol() not in {"H", "D"}
                ),
            }
        )
    conformer = Chem.Conformer(source.GetNumAtoms())
    for atom_index, position in positions.items():
        conformer.SetAtomPosition(atom_index, Point3D(*position))
    source.RemoveAllConformers()
    source.AddConformer(conformer, assignId=True)
    return Chem.MolToMolBlock(source), ligand_atoms


def contact_data(
    protein_atoms: list[dict],
    ligand_atoms: list[dict],
) -> tuple[list[dict], list[dict]]:
    hydrophobic_residues = {
        "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR", "CYS"
    }
    contacts_by_residue = {}
    interactions_by_residue: dict[tuple[str, int, str], dict] = {}
    for protein_atom in protein_atoms:
        residue_key = (protein_atom["chain"], protein_atom["resi"])
        for ligand_atom in ligand_atoms:
            distance = sum(
                (
                    protein_atom["position"][axis]
                    - ligand_atom["position"][axis]
                ) ** 2
                for axis in range(3)
            ) ** 0.5
            residue_contact = contacts_by_residue.get(residue_key)
            if distance <= 4.5 and (
                residue_contact is None or distance < residue_contact["distance"]
            ):
                contacts_by_residue[residue_key] = {
                    "chain": protein_atom["chain"],
                    "resi": protein_atom["resi"],
                    "residue": protein_atom["residue"],
                    "distance": distance,
                    "position": protein_atom["position"],
                }

            interaction_type = None
            cutoff = 0.0
            if (
                (protein_atom["element"] == "N" and ligand_atom["element"] == "O")
                or (protein_atom["element"] == "O" and ligand_atom["element"] == "N")
            ):
                interaction_type = "polar"
                cutoff = 3.5
            elif (
                protein_atom["residue"] in hydrophobic_residues
                and protein_atom["element"] == "C"
                and protein_atom["atom"] not in {"C", "CA"}
                and ligand_atom["hydrophobic"]
            ):
                interaction_type = "hydrophobic"
                cutoff = 4.0
            if interaction_type is None or distance > cutoff:
                continue

            interaction_key = (*residue_key, interaction_type)
            interaction = interactions_by_residue.get(interaction_key)
            if interaction is None or distance < interaction["distance"]:
                interactions_by_residue[interaction_key] = {
                    "type": interaction_type,
                    "distance": distance,
                    "start": protein_atom["position"],
                    "end": ligand_atom["position"],
                    "residue": protein_atom["residue"],
                    "chain": protein_atom["chain"],
                    "resi": protein_atom["resi"],
                    "protein_atom": protein_atom["atom"],
                    "ligand_element": ligand_atom["element"],
                }

    contacts = sorted(contacts_by_residue.values(), key=lambda item: item["distance"])
    interactions = []
    for interaction_type in ("polar", "hydrophobic"):
        interactions.extend(
            sorted(
                (
                    interaction
                    for interaction in interactions_by_residue.values()
                    if interaction["type"] == interaction_type
                ),
                key=lambda item: item["distance"],
            )[:6]
        )
    return contacts, interactions


def build_html(run_dir: Path) -> str:
    poses = parse_pose_models(run_dir / "docked_poses.pdbqt")
    best_pose = poses[0]
    protein_pdb = receptor_pdb_with_helices(run_dir)
    atoms = receptor_atoms(run_dir / "drd2_repaired.pdb")
    ligand_molblock, ligand_atoms = make_ligand_molblock(run_dir, best_pose)
    contacts, interactions = contact_data(atoms, ligand_atoms)
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
    ligand_elements = sorted({atom["element"] for atom in ligand_atoms})
    payload = json.dumps(
        {
            "protein": protein_pdb,
            "ligand": ligand_molblock,
            "score": best_pose["score"],
            "pocket": pocket_groups,
            "labels": labels,
            "interactions": interactions,
            "elementColors": ELEMENT_COLORS,
            "ligandElements": ligand_elements,
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
    {stick: {radius: 0.19, colorscheme: 'Jmol'}}
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
  stick: {radius: 0.27, color: '#b8c2cc'},
  sphere: {scale: 0.27, color: '#b8c2cc'}
});
for (const [element, color] of Object.entries(data.elementColors)) {
  viewer.setStyle(
    {model: ligand, elem: element},
    {
      stick: {radius: 0.27, color: color},
      sphere: {scale: 0.27, color: color}
    },
    true
  );
}
for (const interaction of data.interactions) {
  const isPolar = interaction.type === 'polar';
  viewer.addCylinder({
    start: {x: interaction.start[0], y: interaction.start[1], z: interaction.start[2]},
    end: {x: interaction.end[0], y: interaction.end[1], z: interaction.end[2]},
    radius: isPolar ? 0.045 : 0.03,
    color: isPolar ? '#f0c75e' : '#55c98a',
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
const interactionList = document.getElementById('interactions');
for (const interaction of data.interactions) {
  const item = document.createElement('li');
  const kind = interaction.type === 'polar'
    ? 'Candidate polar N/O contact'
    : 'Hydrophobic-like C/C proximity';
  item.textContent =
    `${kind}: ${interaction.residue}${interaction.resi} ` +
    `${interaction.protein_atom}···${interaction.ligand_element} ` +
    `${interaction.distance.toFixed(1)} Å`;
  interactionList.appendChild(item);
}
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
            ("pocket", "Pocket atoms (element colors)"),
            ("polar", "Candidate polar N/O contact"),
            ("hydrophobic", "Hydrophobic-like C/C proximity"),
        )
    )
    atom_legend_items = "".join(
        f"<li><span class='swatch atom' style='background:{ELEMENT_COLORS.get(element, '#d5dbe0')}'></span>"
        f"Ligand {element}</li>"
        for element in ligand_elements
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
    .polar {{ background: #f0c75e; }} .hydrophobic {{ background: #55c98a; }}
    .atom {{ border: 1px solid #45545e; }}
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
      <h2>Protein and interaction legend</h2>
      <ul>{legend_items}</ul>
      <h2>Ligand atom elements</h2>
      <ul>{atom_legend_items}</ul>
      <h2>Closest candidate contacts</h2>
      <ul id="interactions"></ul>
      <div class="hint">Pocket side chains use element colors. Candidate polar N/O contacts (gold) and hydrophobic-like C/C proximities (green) are geometric cues only; they do not confirm hydrogen bonds, hydrophobic free energy, or binding forces.</div>
      <div class="hint">Drag to rotate · scroll to zoom · right-drag to pan. Labels identify nearby DRD2 residues and their closest distance to the ligand.</div>
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
