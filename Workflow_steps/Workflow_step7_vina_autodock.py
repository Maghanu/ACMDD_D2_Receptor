"""Prepare and dock a ligand against the human D2 dopamine receptor (DRD2).

Activate the project's Vina environment, then run:
    conda install -c conda-forge vina meeko rdkit biopython pdbfixer openmm
    python Workflow_steps/Workflow_step7_vina_autodock.py

When prompted, enter the path to a ligand .sdf, .mol, or .smi file.
Use --top-predictions to dock the highest predicted-pIC50 candidates instead.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
from datetime import datetime, timezone
import math
import re
import shutil
import subprocess
import sys
import sysconfig
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / "docking" / "runs"
DOCKING_SUMMARY_PATH = RUN_DIR / "docking_summary.csv"
DEFAULT_PREDICTIONS_CSV = ROOT / "Workflow_steps" / "data" / "predicted_pIC50_svr.csv"
PDB_ID = "6CM4"
SHARED_RECEPTOR_DIR = ROOT / "docking" / "prepared" / PDB_ID
PDB_URL = f"https://files.rcsb.org/download/{PDB_ID}.cif"
RISPERIDONE_COMPONENT_ID = "8NU"

# 6CM4 contains DRD2 sequence positions 1-188 and 349-430 separated by a
# T4 lysozyme fusion. Keep only the DRD2 portions of the construct.
DRD2_SEQUENCE_RANGES = ((1, 188), (349, 430))
BOX_PADDING_ANGSTROM = 8.0
EXHAUSTIVENESS = 8
NUM_MODES = 9


def require_dependencies() -> dict[str, str]:
    if sys.platform == "win32":
        raise RuntimeError(
            "This docking script uses AutoDock Vina's Python bindings, which are "
            "supported on Linux and macOS, not native Windows. Run it in Ubuntu "
            "under WSL instead. In your WSL terminal, activate the drd2-vina "
            "environment and run:\n"
            "  conda install -c conda-forge vina meeko rdkit biopython "
            "pdbfixer openmm\n"
            "Then run this script from the WSL terminal."
        )

    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict  # noqa: F401
        from openmm.app import PDBFile  # noqa: F401
        from pdbfixer import PDBFixer  # noqa: F401
        from rdkit import Chem  # noqa: F401
        from rdkit.Chem import AllChem  # noqa: F401
        from vina import Vina  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "A required package is missing. Activate the Vina environment and run:\n"
            "  conda install -c conda-forge vina meeko rdkit biopython pdbfixer openmm\n"
            "Then run this script again with `python vina_autodock.py`."
        ) from exc

    scripts_dir = Path(sysconfig.get_path("scripts"))
    tools = {}
    for tool in ("mk_prepare_receptor.py", "mk_prepare_ligand.py"):
        executable = shutil.which(tool)
        if executable is None:
            candidate = scripts_dir / (
                f"{tool}.exe" if sys.platform == "win32" else tool
            )
            if candidate.is_file():
                executable = str(candidate)
        if executable is not None:
            tools[tool] = executable
    missing_tools = set(("mk_prepare_receptor.py", "mk_prepare_ligand.py")) - tools.keys()
    if missing_tools:
        raise RuntimeError(
            "Meeko preparation tools were not found. Activate the Vina environment and run:\n"
            "  conda install -c conda-forge meeko\n"
            "Then run this script again from that environment."
        )
    return tools


def cif_values(cif: dict, key: str) -> list[str]:
    value = cif.get(key, [])
    if isinstance(value, str):
        return [value]
    return list(value)


def download_structure(path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        PDB_URL,
        headers={"User-Agent": "DRD2-Vina-docking/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            path.write_bytes(response.read())
    except (urllib.error.URLError, TimeoutError) as exc:
        raise RuntimeError(
            f"Could not download DRD2 structure {PDB_ID} from RCSB PDB: {exc}"
        ) from exc


def make_receptor_and_box(cif_path: Path, receptor_pdb: Path) -> tuple[list[float], list[float]]:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    cif = MMCIF2Dict(str(cif_path))
    component_ids = set(cif_values(cif, "_chem_comp.id"))
    if RISPERIDONE_COMPONENT_ID not in component_ids:
        raise ValueError(
            f"Expected risperidone component {RISPERIDONE_COMPONENT_ID} was not "
            f"found in PDB {PDB_ID}; stopping rather than guessing a docking box."
        )

    fields = {
        "record": "_atom_site.group_PDB",
        "entity": "_atom_site.label_entity_id",
        "label_chain": "_atom_site.label_asym_id",
        "sequence": "_atom_site.label_seq_id",
        "component": "_atom_site.label_comp_id",
        "atom": "_atom_site.label_atom_id",
        "altloc": "_atom_site.label_alt_id",
        "residue": "_atom_site.auth_seq_id",
        "insertion": "_atom_site.pdbx_PDB_ins_code",
        "x": "_atom_site.Cartn_x",
        "y": "_atom_site.Cartn_y",
        "z": "_atom_site.Cartn_z",
        "occupancy": "_atom_site.occupancy",
        "b_factor": "_atom_site.B_iso_or_equiv",
        "element": "_atom_site.type_symbol",
    }
    atom_columns = {name: cif_values(cif, key) for name, key in fields.items()}
    atom_count = len(atom_columns["record"])
    if not atom_count or any(len(values) != atom_count for values in atom_columns.values()):
        raise ValueError(f"PDB {PDB_ID} has incomplete atom records.")

    alternate_occupancy: dict[tuple[str, str, str], dict[str, float]] = {}
    for index in range(atom_count):
        if atom_columns["record"][index] != "ATOM":
            continue
        altloc = atom_columns["altloc"][index]
        if altloc in (".", "?", " "):
            continue
        residue_key = (
            atom_columns["entity"][index],
            atom_columns["label_chain"][index],
            atom_columns["sequence"][index],
        )
        occupancies = alternate_occupancy.setdefault(residue_key, {})
        occupancies[altloc] = occupancies.get(altloc, 0.0) + float(
            atom_columns["occupancy"][index]
        )
    preferred_altloc = {
        residue_key: max(occupancies, key=occupancies.get)
        for residue_key, occupancies in alternate_occupancy.items()
    }

    native_ligand_points: list[tuple[float, float, float]] = []
    receptor_lines = ["REMARK DRD2 receptor prepared from RCSB PDB 6CM4\n"]
    serial = 1
    for index in range(atom_count):
        record = atom_columns["record"][index]
        component = atom_columns["component"][index]
        element = atom_columns["element"][index].strip().upper()
        point = tuple(
            float(atom_columns[axis][index]) for axis in ("x", "y", "z")
        )

        if (
            record == "HETATM"
            and component == RISPERIDONE_COMPONENT_ID
            and element != "H"
        ):
            native_ligand_points.append(point)

        if record != "ATOM" or atom_columns["entity"][index] != "1":
            continue
        try:
            sequence_number = int(atom_columns["sequence"][index])
        except ValueError:
            continue
        segment = next(
            (
                segment_index
                for segment_index, (start, end) in enumerate(DRD2_SEQUENCE_RANGES)
                if start <= sequence_number <= end
            ),
            None,
        )
        if segment is None:
            continue

        altloc = atom_columns["altloc"][index]
        if altloc not in (".", "?", " "):
            residue_key = (
                atom_columns["entity"][index],
                atom_columns["label_chain"][index],
                atom_columns["sequence"][index],
            )
            if altloc != preferred_altloc.get(residue_key):
                continue
        residue_id = int(atom_columns["residue"][index])
        insertion = atom_columns["insertion"][index]
        if insertion in (".", "?", " "):
            insertion = " "
        chain = "A" if segment == 0 else "B"
        atom_name = atom_columns["atom"][index]
        atom_field = (
            f"{atom_name:<4s}"
            if len(atom_name) == 4 or atom_name[:1].isdigit()
            else f" {atom_name:<3s}"
        )
        residue_name = component
        x, y, z = point
        occupancy = float(atom_columns["occupancy"][index])
        b_factor = float(atom_columns["b_factor"][index])
        receptor_lines.append(
            f"ATOM  {serial:5d} {atom_field} {residue_name:>3s} {chain}"
            f"{residue_id:4d}{insertion:1s}   {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{occupancy:6.2f}{b_factor:6.2f}          {element:>2s}\n"
        )
        serial += 1

    if not native_ligand_points:
        raise ValueError(f"No co-crystallized risperidone coordinates found in {PDB_ID}.")
    if serial <= 2:
        raise ValueError(f"No DRD2 protein atoms were extracted from PDB {PDB_ID}.")

    receptor_lines.extend(("TER\n", "END\n"))
    receptor_pdb.parent.mkdir(parents=True, exist_ok=True)
    receptor_pdb.write_text("".join(receptor_lines), encoding="ascii")

    axes = list(zip(*native_ligand_points))
    minimum = [min(axis) for axis in axes]
    maximum = [max(axis) for axis in axes]
    center = [(low + high) / 2 for low, high in zip(minimum, maximum)]
    size = [
        high - low + BOX_PADDING_ANGSTROM
        for low, high in zip(minimum, maximum)
    ]
    if not all(math.isfinite(value) for value in center) or not all(
        math.isfinite(value) and value > 0 for value in size
    ):
        raise ValueError("Could not calculate a valid docking box from native risperidone.")
    return center, size


def repair_missing_atoms(receptor_pdb: Path, repaired_pdb: Path) -> None:
    """Rebuild missing atoms and relax the structure with positional restraints."""
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(receptor_pdb))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    repaired_pdb.parent.mkdir(parents=True, exist_ok=True)
    with repaired_pdb.open("w", encoding="ascii") as output_file:
        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            output_file,
            keepIds=True,
        )
    minimize_repaired_receptor(receptor_pdb, repaired_pdb)
    corrected = correct_terminal_oxygen_geometry(repaired_pdb)
    if corrected:
        print(f"Corrected anomalous terminal oxygen geometry in {corrected} residue(s).")


def prepare_shared_receptor(tools: dict[str, str]) -> tuple[list[float], list[float]]:
    structure_path = SHARED_RECEPTOR_DIR / f"{PDB_ID}.cif"
    receptor_pdb = SHARED_RECEPTOR_DIR / "drd2_clean.pdb"
    repaired_receptor_pdb = SHARED_RECEPTOR_DIR / "drd2_repaired.pdb"
    receptor_pdbqt = SHARED_RECEPTOR_DIR / "drd2.pdbqt"

    download_structure(structure_path)
    center, size = make_receptor_and_box(structure_path, receptor_pdb)
    if all(path.is_file() for path in (repaired_receptor_pdb, receptor_pdbqt)):
        print(f"Reusing shared prepared receptor: {receptor_pdbqt}")
        return center, size

    print("Preparing shared DRD2 receptor for ligand screening...")
    repair_missing_atoms(receptor_pdb, repaired_receptor_pdb)
    run_preparation(
        [
            tools["mk_prepare_receptor.py"],
            "--read_pdb", str(repaired_receptor_pdb),
            "--write_pdbqt", str(receptor_pdbqt),
            "--box_center", *(f"{value:.3f}" for value in center),
            "--box_size", *(f"{value:.3f}" for value in size),
            "--delete_bad_res_from_box_radius", "8",
        ],
        "shared DRD2 receptor",
    )
    return center, size


def minimize_repaired_receptor(
    receptor_pdb: Path,
    repaired_pdb: Path,
) -> None:
    """Relax rebuilt atoms while strongly restraining experimental heavy atoms."""
    import openmm
    from openmm import unit
    from openmm.app import ForceField, Modeller, NoCutoff, PDBFile, Simulation

    original_atoms = set()
    for line in receptor_pdb.read_text(encoding="ascii").splitlines():
        if line.startswith("ATOM  "):
            original_atoms.add(
                (line[21], line[22:26].strip(), line[26], line[12:16].strip())
            )

    pdb = PDBFile(str(repaired_pdb))
    modeller = Modeller(pdb.topology, pdb.positions)
    forcefield = ForceField("amber14-all.xml")
    modeller.addHydrogens(forcefield, pH=7.0)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
    )
    restraints = openmm.CustomExternalForce(
        "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
    )
    for parameter in ("k", "x0", "y0", "z0"):
        restraints.addPerParticleParameter(parameter)

    positions = modeller.positions
    for atom in modeller.topology.atoms():
        if atom.element is None or atom.element.symbol in {"H", "D"}:
            continue
        residue = atom.residue
        atom_key = (
            residue.chain.id,
            residue.id.strip(),
            residue.insertionCode,
            atom.name,
        )
        position = positions[atom.index].value_in_unit(unit.nanometer)
        restraint_strength = 10000.0 if atom_key in original_atoms else 100.0
        restraints.addParticle(
            atom.index,
            [restraint_strength, position.x, position.y, position.z],
        )
    system.addForce(restraints)

    simulation = Simulation(
        modeller.topology,
        system,
        openmm.VerletIntegrator(0.001 * unit.picoseconds),
    )
    simulation.context.setPositions(positions)
    openmm.LocalEnergyMinimizer.minimize(
        simulation.context,
        tolerance=10.0 * unit.kilojoules_per_mole / unit.nanometer,
        maxIterations=1000,
    )
    minimized_positions = simulation.context.getState(getPositions=True).getPositions()
    with repaired_pdb.open("w", encoding="ascii") as output_file:
        PDBFile.writeFile(
            modeller.topology,
            minimized_positions,
            output_file,
            keepIds=True,
        )

    lines = repaired_pdb.read_text(encoding="ascii").splitlines(keepends=True)
    repaired_pdb.write_text(
        "".join(
            line for line in lines
            if not line.startswith(("ATOM  ", "HETATM"))
            or line[76:78].strip().upper() not in {"H", "D"}
        ),
        encoding="ascii",
    )


def correct_terminal_oxygen_geometry(repaired_pdb: Path) -> int:
    """Move mispositioned terminal oxygens away from false RDKit bond distances."""
    lines = repaired_pdb.read_text(encoding="ascii").splitlines(keepends=True)
    residues: dict[tuple[str, str], dict[str, int]] = {}
    for index, line in enumerate(lines):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue_key = (line[21], line[22:27])
        residues.setdefault(residue_key, {})[line[12:16].strip()] = index

    corrected = 0
    for atoms in residues.values():
        if not all(name in atoms for name in ("CA", "C", "O", "OXT")):
            continue

        def coordinates(atom_name: str) -> tuple[float, float, float]:
            line = lines[atoms[atom_name]]
            return tuple(float(line[start:start + 8]) for start in (30, 38, 46))

        ca, carbonyl_c, oxygen, terminal_oxygen = (
            coordinates(name) for name in ("CA", "C", "O", "OXT")
        )
        if math.dist(ca, terminal_oxygen) >= 1.8:
            continue

        ca_vector = tuple(ca[i] - carbonyl_c[i] for i in range(3))
        oxygen_vector = tuple(oxygen[i] - carbonyl_c[i] for i in range(3))

        def unit(vector: tuple[float, float, float]) -> tuple[float, float, float]:
            magnitude = math.sqrt(sum(value * value for value in vector))
            if magnitude < 1e-8:
                raise ValueError("Could not repair degenerate terminal oxygen geometry.")
            return tuple(value / magnitude for value in vector)

        ca_direction = unit(ca_vector)
        oxygen_direction = unit(oxygen_vector)
        opposite_bisector = tuple(
            -(ca_direction[i] + oxygen_direction[i]) for i in range(3)
        )
        terminal_direction = unit(opposite_bisector)
        repaired_position = tuple(
            carbonyl_c[i] + 1.25 * terminal_direction[i] for i in range(3)
        )

        index = atoms["OXT"]
        lines[index] = (
            lines[index][:30]
            + "".join(f"{value:8.3f}" for value in repaired_position)
            + lines[index][54:]
        )
        corrected += 1

    if corrected:
        repaired_pdb.write_text("".join(lines), encoding="ascii")
    return corrected


def create_ligand_sdf(source_path: Path, sdf_path: Path) -> None:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    suffix = source_path.suffix.lower()
    molecule = None
    if suffix == ".sdf":
        supplier = Chem.SDMolSupplier(str(source_path), removeHs=False)
        molecule = next((candidate for candidate in supplier if candidate is not None), None)
    elif suffix == ".mol":
        molecule = Chem.MolFromMolFile(str(source_path), removeHs=False)
    elif suffix in (".smi", ".smiles", ".txt"):
        smiles_tokens = source_path.read_text(encoding="utf-8-sig").strip().split()
        if not smiles_tokens:
            raise ValueError(f"No SMILES string found in {source_path}.")
        molecule = Chem.MolFromSmiles(smiles_tokens[0])
        if molecule is None:
            raise ValueError(
                f"Could not parse SMILES token {smiles_tokens[0]!r} from "
                f"{source_path}. Check the file contents and remove any "
                "unsupported characters."
            )
    else:
        raise ValueError("Use a ligand file ending in .sdf, .mol, .smi, .smiles, or .txt.")

    if molecule is None:
        raise ValueError(f"Could not read a valid molecule from {source_path}.")

    fragments = Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True)
    if len(fragments) > 1:
        organic_fragments = [
            fragment
            for fragment in fragments
            if any(atom.GetAtomicNum() == 6 for atom in fragment.GetAtoms())
        ]
        if len(organic_fragments) != 1:
            raise ValueError(
                f"Ligand input contains {len(fragments)} disconnected fragments "
                f"and {len(organic_fragments)} carbon-containing components. "
                "Provide a single ligand or a salt with exactly one "
                "carbon-containing component; refusing to guess which component "
                "to dock."
            )

        ignored_fragments = [
            Chem.MolToSmiles(fragment)
            for fragment in fragments
            if fragment is not organic_fragments[0]
        ]
        molecule = organic_fragments[0]
        print(
            "Ignoring disconnected non-carbon fragment(s) for docking: "
            f"{', '.join(ignored_fragments)}. "
            f"Using ligand component {Chem.MolToSmiles(molecule)}."
        )

    molecule = Chem.AddHs(molecule, addCoords=molecule.GetNumConformers() > 0)
    has_3d = (
        molecule.GetNumConformers() > 0
        and molecule.GetConformer().Is3D()
    )
    if not has_3d:
        molecule.RemoveAllConformers()
        if AllChem.EmbedMolecule(molecule, randomSeed=22) != 0:
            raise ValueError("Could not generate a 3D conformer for the ligand.")
        if AllChem.MMFFHasAllMoleculeParams(molecule):
            AllChem.MMFFOptimizeMolecule(molecule)

    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(molecule)
    writer.close()


def run_preparation(command: list[str], label: str) -> None:
    print(f"Preparing {label}...")
    subprocess.run(command, check=True)


def dock(receptor_pdbqt: Path, ligand_pdbqt: Path, output_path: Path,
         log_path: Path, center: list[float], size: list[float]) -> float:
    from vina import Vina

    vina = Vina(sf_name="vina", cpu=0)
    vina.set_receptor(str(receptor_pdbqt))
    vina.set_ligand_from_file(str(ligand_pdbqt))
    vina.compute_vina_maps(center=center, box_size=size)
    with log_path.open("w", encoding="utf-8") as log_file:
        with contextlib.redirect_stdout(log_file):
            vina.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=NUM_MODES)
            vina.write_poses(str(output_path), n_poses=NUM_MODES, overwrite=True)
            scores = vina.score()
    print(f"Best Vina score: {scores[0]:.2f} kcal/mol")
    return float(scores[0])


def load_top_predictions(csv_path: Path, count: int) -> list[dict[str, str]]:
    if count < 1:
        raise ValueError("The number of candidates to dock must be at least 1.")
    if not csv_path.is_file():
        raise FileNotFoundError(f"Prediction CSV not found: {csv_path}")

    with csv_path.open(encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file)
        smiles_column = next(
            (name for name in ("smiles", "SMILES") if name in (reader.fieldnames or [])),
            None,
        )
        if reader.fieldnames is None or smiles_column is None:
            raise ValueError(f"{csv_path} must contain a 'smiles' or 'SMILES' column.")
        score_column = next(
            (
                name for name in ("predicted_pIC50", "pred_pIC50")
                if name in reader.fieldnames
            ),
            None,
        )
        if score_column is None:
            raise ValueError(
                f"{csv_path} must contain 'predicted_pIC50' or 'pred_pIC50'."
            )

        candidates = []
        for row_number, row in enumerate(reader, start=2):
            smiles = (row.get(smiles_column) or "").strip()
            if not smiles:
                raise ValueError(f"{csv_path}, row {row_number}: SMILES is empty.")
            try:
                score = float(row[score_column])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"{csv_path}, row {row_number}: invalid predicted pIC50 "
                    f"{row.get(score_column)!r}."
                ) from exc
            if not math.isfinite(score):
                raise ValueError(
                    f"{csv_path}, row {row_number}: predicted pIC50 must be finite."
                )
            compound_id = (
                row.get("compound_id")
                or row.get("zincid")
                or row.get("zinc_id")
                or ""
            ).strip()
            candidates.append(
                {
                    "compound_id": compound_id or f"row_{row_number}",
                    "smiles": smiles,
                    "predicted_pIC50": str(score),
                }
            )

    if len(candidates) < count:
        raise ValueError(
            f"{csv_path} contains {len(candidates)} candidate(s), fewer than "
            f"the requested {count}."
        )
    candidates.sort(
        key=lambda candidate: float(candidate["predicted_pIC50"]),
        reverse=True,
    )
    return candidates[:count]


def record_docking_result(
    run_dir: Path,
    ligand_source: Path,
    vina_score: float,
    metadata: dict[str, str] | None = None,
) -> None:
    summary_fields = (
        "run_name",
        "run_directory",
        "input_file",
        "compound_id",
        "prediction_rank",
        "predicted_pIC50",
        "vina_score_kcal_mol",
        "score_description",
        "timestamp_utc",
        "smiles",
    )
    summary_path = DOCKING_SUMMARY_PATH
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    run_directory = str(run_dir.relative_to(ROOT))
    existing_rows = []
    if summary_path.exists():
        with summary_path.open(encoding="utf-8-sig", newline="") as summary_file:
            reader = csv.DictReader(summary_file)
            if reader.fieldnames != list(summary_fields):
                raise ValueError(
                    f"Cannot update docking summary with unexpected columns: "
                    f"{summary_path}"
                )
            existing_rows = list(reader)

    pose_path = run_dir / "docked_poses.pdbqt"
    pose_smiles = ""
    for line in pose_path.read_text(encoding="ascii").splitlines():
        if line.startswith("REMARK SMILES "):
            pose_smiles = line.partition("REMARK SMILES ")[2]
            break

    row = {
        "run_name": run_dir.name,
        "run_directory": run_directory,
        "input_file": str(ligand_source),
        "compound_id": "",
        "prediction_rank": "",
        "predicted_pIC50": "",
        "vina_score_kcal_mol": f"{vina_score:.3f}",
        "score_description": "AutoDock Vina score (approximate; not experimental affinity)",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "smiles": pose_smiles,
    }
    row.update(metadata or {})
    existing_rows = [
        existing
        for existing in existing_rows
        if existing["run_directory"] != run_directory
    ]
    existing_rows.append(row)

    temporary_path = summary_path.with_suffix(summary_path.suffix + ".tmp")
    with temporary_path.open("w", encoding="utf-8", newline="") as summary_file:
        writer = csv.DictWriter(summary_file, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(existing_rows)
    temporary_path.replace(summary_path)
    print(f"Docking summary updated: {summary_path}")


def dock_ligand(
    ligand_source: Path,
    run_name: str,
    tools: dict[str, str],
    center: list[float],
    size: list[float],
    summary_metadata: dict[str, str] | None = None,
    run_dir: Path | None = None,
) -> tuple[Path, float]:
    run_dir = run_dir or RUN_DIR / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    for filename in (
        f"{PDB_ID}.cif",
        "drd2_clean.pdb",
        "drd2_repaired.pdb",
        "drd2.pdbqt",
    ):
        shutil.copyfile(SHARED_RECEPTOR_DIR / filename, run_dir / filename)

    receptor_pdbqt = run_dir / "drd2.pdbqt"
    ligand_sdf = run_dir / "ligand_3d.sdf"
    ligand_pdbqt = run_dir / "ligand.pdbqt"
    output_path = run_dir / "docked_poses.pdbqt"
    log_path = run_dir / "vina.log"

    create_ligand_sdf(ligand_source, ligand_sdf)
    run_preparation(
        [
            tools["mk_prepare_ligand.py"],
            "-i", str(ligand_sdf),
            "-o", str(ligand_pdbqt),
        ],
        f"ligand {run_name}",
    )
    score = dock(receptor_pdbqt, ligand_pdbqt, output_path, log_path, center, size)
    record_docking_result(run_dir, ligand_source, score, summary_metadata)
    print(f"Docking poses: {output_path}")
    print(f"Vina log: {log_path}")
    return run_dir, score


def dock_top_predictions(
    csv_path: Path,
    count: int,
    tools: dict[str, str],
    output_dir: Path = RUN_DIR,
) -> Path:
    candidates = load_top_predictions(csv_path, count)
    center, size = prepare_shared_receptor(tools)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / f"top_{count}_docking_summary.csv"
    results = []

    for rank, candidate in enumerate(candidates, start=1):
        safe_id = re.sub(r"[^A-Za-z0-9._-]+", "_", candidate["compound_id"]).strip("._")
        if not safe_id:
            safe_id = f"candidate_{rank:02d}"
        run_name = f"top_{rank:02d}_{safe_id}"
        run_dir = output_dir / run_name
        ligand_source = run_dir / "candidate.smi"
        run_dir.mkdir(parents=True, exist_ok=True)
        ligand_source.write_text(candidate["smiles"] + "\n", encoding="utf-8")

        _, vina_score = dock_ligand(
            ligand_source,
            run_name,
            tools,
            center,
            size,
            {
                "compound_id": candidate["compound_id"],
                "prediction_rank": str(rank),
                "predicted_pIC50": candidate["predicted_pIC50"],
            },
            run_dir=run_dir,
        )
        results.append(
            {
                "rank": rank,
                "compound_id": candidate["compound_id"],
                "smiles": candidate["smiles"],
                "predicted_pIC50": candidate["predicted_pIC50"],
                "vina_score_kcal_mol": vina_score,
                "run_directory": str(run_dir.relative_to(ROOT)),
            }
        )
        with summary_path.open("w", encoding="utf-8", newline="") as output_file:
            writer = csv.DictWriter(output_file, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)

    print(f"Docking summary: {summary_path}")
    print(
        "Note: predicted pIC50 and Vina scores are distinct computational "
        "estimates, not experimental binding measurements."
    )
    return summary_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare and dock DRD2 ligands with AutoDock Vina."
    )
    parser.add_argument(
        "--top-predictions",
        action="store_true",
        help="dock the highest-pIC50 candidates from the prediction CSV",
    )
    parser.add_argument(
        "--predictions-csv",
        type=Path,
        default=DEFAULT_PREDICTIONS_CSV,
        help=f"prediction CSV (default: {DEFAULT_PREDICTIONS_CSV})",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="number of top predictions to dock (default: 10)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=RUN_DIR,
        help=f"directory for top-prediction runs and summary (default: {RUN_DIR})",
    )
    args = parser.parse_args()

    try:
        tools = require_dependencies()
        if args.top_predictions:
            dock_top_predictions(
                args.predictions_csv.expanduser().resolve(),
                args.top,
                tools,
                args.output_dir.expanduser().resolve(),
            )
        else:
            raw_path = input(
                "Enter the path to your ligand .sdf, .mol, or .smi file: "
            ).strip().strip('"')
            if not raw_path:
                raise ValueError("No ligand file was provided.")
            ligand_source = Path(raw_path).expanduser().resolve()
            if not ligand_source.is_file():
                raise FileNotFoundError(f"Ligand file not found: {ligand_source}")

            center, size = prepare_shared_receptor(tools)
            run_dir, _ = dock_ligand(
                ligand_source, ligand_source.stem, tools, center, size
            )
            print(
                "Note: this is a docking score, not an experimental binding affinity. "
                f"Review the receptor preparation and pose in {run_dir}."
            )
        return 0
    except (
        RuntimeError,
        ValueError,
        FileNotFoundError,
        OSError,
        csv.Error,
        subprocess.CalledProcessError,
    ) as exc:
        print(f"Docking could not be completed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
