#!/usr/bin/env python3
"""
Dopamine D2 receptor (CHEMBL217 / UniProt P14416): ligand and interaction analysis.

Part 1  ChEMBL: binding activities, curation, descriptors, scaffolds, action type.
Part 2  PDB:    residue-level interactions (PLIP), BW numbering (GPCRdb),
                PDB ligand -> ChEMBL link, interaction frequency by ligand class.

Install (in a virtual environment):
    python3 -m venv .venv
    .venv/bin/pip install -r requirements.txt
    .venv/bin/pip install --no-deps plip
    (PLIP needs Open Babel; requirements.txt uses the prebuilt openbabel-wheel)

Usage:
    .venv/bin/python CHEMBL_DRUG_D2.py
Outputs are written to ./d2_analysis/
"""

import os
import requests
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

TARGET = "CHEMBL217"
UNIPROT = "P14416"
OUT = "d2_analysis"
PDB_DIR = os.path.join(OUT, "pdb")
os.makedirs(PDB_DIR, exist_ok=True)

# Optional manual class override for PDB ligands, e.g. {"RIS": "ANTAGONIST"}
MANUAL_CLASS = {}

# Non-drug HETATM groups to ignore (lipids, detergents, solvents, ions)
EXCLUDE_HET = {
    "HOH", "OLC", "OLA", "OLB", "CLR", "PEG", "PGE", "PLM", "STE", "MYR",
    "LMT", "DMS", "GOL", "EDO", "SO4", "PO4", "NA", "CL", "ZN", "MG", "CA",
    "GDP", "GTP", "ACT", "NAG", "BMA", "PC1", "POV", "D10", "HEX",
}

# Aliphatic amine (approximate basic centre for the Asp3.32 salt bridge)
BASIC_N = Chem.MolFromSmarts(
    "[NX3;!$(N-[C,S,P]=[O,S,N]);!$(N-a);!$(N-[#7,#8]);!$(N-C=[C,N]);!$(N#*)]"
)


# ----------------------------------------------------------------------------
# Part 1: ChEMBL
# ----------------------------------------------------------------------------
def fetch_activities():
    """Binding-assay activities (Ki, Kd, IC50) with a pChEMBL value."""
    qs = new_client.activity.filter(
        target_chembl_id=TARGET,
        assay_type="B",
        standard_type__in=["Ki", "Kd", "IC50"],
        standard_relation="=",
        pchembl_value__isnull=False,
    ).only([
        "molecule_chembl_id", "canonical_smiles", "standard_type",
        "pchembl_value", "assay_chembl_id", "document_chembl_id",
    ])
    df = pd.DataFrame(list(qs))
    df["pchembl_value"] = df["pchembl_value"].astype(float)
    return df.dropna(subset=["canonical_smiles", "pchembl_value"])


def fetch_mechanisms():
    """Action type (agonist/antagonist/...) for annotated drugs."""
    m = pd.DataFrame(list(new_client.mechanism.filter(target_chembl_id=TARGET)))
    if m.empty:
        return pd.DataFrame(columns=["molecule_chembl_id", "action_type"])
    return (m[["molecule_chembl_id", "action_type", "mechanism_of_action"]]
            .drop_duplicates("molecule_chembl_id"))


def describe(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return pd.Series(dtype=float)
    try:
        scaf = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    except Exception:
        scaf = ""
    return pd.Series({
        "MW": Descriptors.MolWt(mol),
        "logP": Crippen.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "RotB": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "AromRings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "basic_amine": int(mol.HasSubstructMatch(BASIC_N)),
        "scaffold": scaf,
    })


def part1_chembl():
    print("[ChEMBL] fetching activities (this can take several minutes)...")
    acts = fetch_activities()
    acts.to_csv(os.path.join(OUT, "chembl_activities_raw.csv"), index=False)

    cpd = (acts.groupby("molecule_chembl_id")
           .agg(smiles=("canonical_smiles", "first"),
                pchembl_median=("pchembl_value", "median"),
                pchembl_sd=("pchembl_value", "std"),
                n_meas=("pchembl_value", "size"),
                n_assays=("assay_chembl_id", "nunique"))
           .reset_index())
    cpd["discordant"] = cpd["pchembl_sd"] > 1.0  # >10-fold spread between measurements

    cpd = pd.concat([cpd, cpd["smiles"].apply(describe)], axis=1)
    mech = fetch_mechanisms()
    cpd = cpd.merge(mech, on="molecule_chembl_id", how="left")
    cpd["action_type"] = cpd["action_type"].fillna("UNANNOTATED")
    cpd.to_csv(os.path.join(OUT, "chembl_compounds.csv"), index=False)

    summary = (cpd.groupby("action_type")
               .agg(n=("molecule_chembl_id", "size"),
                    pchembl_median=("pchembl_median", "median"),
                    MW=("MW", "median"), logP=("logP", "median"),
                    TPSA=("TPSA", "median"),
                    frac_basic_amine=("basic_amine", "mean"))
               .round(2))
    summary.to_csv(os.path.join(OUT, "chembl_summary_by_action.csv"))
    print(summary)

    scaf = (cpd[cpd.scaffold != ""].groupby("scaffold")
            .agg(n=("molecule_chembl_id", "size"),
                 pchembl_median=("pchembl_median", "median"))
            .sort_values("n", ascending=False).head(25).round(2))
    scaf.to_csv(os.path.join(OUT, "chembl_top_scaffolds.csv"))

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    ax[0].hist(cpd["pchembl_median"], bins=40, color="0.4")
    ax[0].set_xlabel("median pChEMBL")
    ax[0].set_ylabel("compounds")
    annotated = cpd[cpd.action_type != "UNANNOTATED"]
    groups = sorted(annotated.action_type.unique())
    for a, col in zip(ax[1:], ["MW", "logP"]):
        a.boxplot([annotated.loc[annotated.action_type == g, col] for g in groups])
        a.set_xticks(range(1, len(groups) + 1))
        a.set_xticklabels(groups)
        a.set_ylabel(col)
        a.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "chembl_overview.png"), dpi=200)
    plt.close(fig)
    return cpd


# ----------------------------------------------------------------------------
# Part 2: PDB interactions
# ----------------------------------------------------------------------------
def d2_pdb_entries():
    q = {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": UNIPROT}},
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    r = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=q, timeout=60)
    r.raise_for_status()
    return [h["identifier"] for h in r.json()["result_set"]]


def download_pdb(pdb_id):
    path = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
    if not os.path.exists(path):
        r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=120)
        if r.status_code != 200:
            return None
        with open(path, "w") as fh:
            fh.write(r.text)
    return path


def gpcrdb_map():
    """UniProt sequence number -> (amino acid, Ballesteros-Weinstein number)."""
    try:
        r = requests.get("https://gpcrdb.org/services/residues/extended/drd2_human/",
                         timeout=60)
        r.raise_for_status()
        out = {}
        for x in r.json():
            gn = x.get("display_generic_number") or ""
            bw = gn.split("x")[0] if gn else ""
            out[int(x["sequence_number"])] = (x["amino_acid"], bw)
        return out
    except Exception as e:
        print(f"[GPCRdb] mapping unavailable ({e}); BW numbers omitted")
        return {}


AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}

ITYPE_NAMES = {
    "hbond": "hbond", "hydroph_interaction": "hydrophobic",
    "saltbridge": "salt_bridge", "pistack": "pi_stack", "pication": "pi_cation",
    "halogenbond": "halogen", "waterbridge": "water_bridge",
    "metal_complex": "metal",
}


def itype_distance(i):
    for attr in ("distance", "distance_ad", "distance_aw"):
        if hasattr(i, attr):
            return float(getattr(i, attr))
    return np.nan


def het_to_chembl(het):
    """PDB ligand code -> ChEMBL ID via InChIKey."""
    try:
        j = requests.get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{het}",
                         timeout=30).json()
        d = j["rcsb_chem_comp_descriptor"]
        key = d.get("in_ch_i_key") or d.get("in_ch_ikey")
        hit = list(new_client.molecule.filter(
            molecule_structures__standard_inchi_key=key).only(["molecule_chembl_id"]))
        return hit[0]["molecule_chembl_id"] if hit else None
    except Exception:
        return None


def analyze_structure(pdb_id, path, bwmap):
    from plip.structure.preparation import PDBComplex
    mol = PDBComplex()
    mol.load_pdb(path)
    mol.analyze()
    rows = []
    for key, inter in mol.interaction_sets.items():
        het, lig_chain, lig_pos = key.split(":")
        if het in EXCLUDE_HET:
            continue
        hits = []
        for i in inter.all_itypes:
            hits.append({
                "interaction": ITYPE_NAMES.get(type(i).__name__, type(i).__name__),
                "restype": i.restype, "resnr": int(i.resnr),
                "reschain": i.reschain, "distance": itype_distance(i),
            })
        if not hits:
            continue
        h = pd.DataFrame(hits)
        # Keep the receptor chain: the chain with most contacts (excludes G protein etc.)
        h = h[h.reschain == h.reschain.value_counts().idxmax()]
        h["pdb"], h["het"], h["lig_site"] = pdb_id, het, key
        rows.append(h)
    if not rows:
        return pd.DataFrame()
    df = pd.concat(rows, ignore_index=True)
    # BW numbering only where the amino acid matches the UniProt sequence
    def bw(r):
        aa, b = bwmap.get(r.resnr, (None, ""))
        return b if aa == AA3.get(r.restype) else ""
    df["bw"] = df.apply(bw, axis=1) if bwmap else ""
    return df


def part2_pdb(cpd):
    ids = d2_pdb_entries()
    print(f"[PDB] {len(ids)} D2 entries: {ids}")
    bwmap = gpcrdb_map()
    frames = []
    for pid in ids:
        path = download_pdb(pid)
        if path is None:
            print(f"[PDB] {pid}: no .pdb file, skipped")
            continue
        try:
            df = analyze_structure(pid, path, bwmap)
        except Exception as e:
            print(f"[PDB] {pid}: PLIP failed ({e})")
            continue
        if not df.empty:
            frames.append(df)
    if not frames:
        print("[PDB] no interactions found")
        return
    inter = pd.concat(frames, ignore_index=True)

    # Link PDB ligands to ChEMBL potency and action type
    links = {h: het_to_chembl(h) for h in inter.het.unique()}
    inter["chembl_id"] = inter.het.map(links)
    mech = cpd.set_index("molecule_chembl_id")
    inter["action_type"] = inter.chembl_id.map(mech["action_type"])
    inter["pchembl_median"] = inter.chembl_id.map(mech["pchembl_median"])
    inter["action_type"] = (inter.het.map(MANUAL_CLASS)
                            .fillna(inter.action_type).fillna("UNASSIGNED"))
    inter["residue"] = (inter.restype + inter.resnr.astype(str)
                        + inter.bw.map(lambda s: f" ({s})" if s else ""))
    inter.to_csv(os.path.join(OUT, "pdb_interactions.csv"), index=False)

    # Interaction frequency per residue, by ligand class
    inter["ligand"] = inter.pdb + ":" + inter.het
    n_lig = inter.groupby("action_type")["ligand"].nunique()
    freq = (inter.drop_duplicates(["ligand", "residue", "interaction"])
            .groupby(["action_type", "residue", "interaction"])["ligand"]
            .nunique().rename("n_ligands").reset_index())
    freq["fraction"] = freq.apply(lambda r: r.n_ligands / n_lig[r.action_type], axis=1)
    freq.sort_values(["action_type", "fraction"], ascending=[True, False]).to_csv(
        os.path.join(OUT, "pdb_interaction_frequency.csv"), index=False)

    # Residue x ligand contact map
    m = (inter.drop_duplicates(["ligand", "residue"])
         .assign(v=1).pivot(index="residue", columns="ligand", values="v").fillna(0))
    m = m.loc[m.sum(axis=1).sort_values(ascending=False).head(40).index]
    fig, ax = plt.subplots(figsize=(0.45 * m.shape[1] + 3, 0.28 * m.shape[0] + 2))
    ax.imshow(m.values, aspect="auto", cmap="Greys")
    ax.set_xticks(range(m.shape[1]))
    ax.set_xticklabels(m.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(m.shape[0]))
    ax.set_yticklabels(m.index, fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "pdb_contact_map.png"), dpi=200)
    plt.close(fig)
    print(f"[PDB] {inter.ligand.nunique()} ligands, {len(inter)} interactions written")


if __name__ == "__main__":
    compounds = part1_chembl()
    part2_pdb(compounds)
    print(f"Done. Results in ./{OUT}/")
