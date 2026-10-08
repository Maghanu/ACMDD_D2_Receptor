import math
from pathlib import Path

import requests

import pandas as pd
from rdkit.Chem import PandasTools

HERE = Path(_dh[-1])
DATA = HERE / "data"

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

# Function to get JSON data from a the ChEMBL REST endpoint, handling pagination 
def get_json(endpoint, params=None):
    params = params.copy() if params else {}
    params.setdefault("limit", 1000)

    records = []
    offset = 0

    while True:
        params["offset"] = offset

        response = requests.get(
            f"{BASE_URL}/{endpoint}.json",
            params=params,
            headers={"Accept": "application/json"},
            timeout=30,
        )

        response.raise_for_status()
        data = response.json()
        collection_name = next(
            key
            for key in data.keys()
            if isinstance(data[key], list)
        )

        batch = data[collection_name]

        if not batch:
            break

        records.extend(batch)

        if len(batch) < params["limit"]:
            break

        offset += params["limit"]

    return records

    # Search for targets in chembl from uniprot_id and return a pandas DataFrame with the results
def get_targets(uniprot_id):
    targets = get_json(
        "target",
        {
            "target_components__accession": uniprot_id
        }
    )

    targets = pd.DataFrame(targets)

    # Add columns to the DataFrame if it is not empty
    if not targets.empty:
        columns = [
            "target_chembl_id",
            "organism",
            "pref_name",
            "target_type",
        ]

        targets = targets[columns]

    return targets

# Search for targets in chembl from uniprot_id and return a pandas DataFrame with the results
def get_targets(uniprot_id):
    targets = get_json(
        "target",
        {
            "target_components__accession": uniprot_id
        }
    )

    targets = pd.DataFrame(targets)

    # Add columns to the DataFrame if it is not empty
    if not targets.empty:
        columns = [
            "target_chembl_id",
            "organism",
            "pref_name",
            "target_type",
        ]

        targets = targets[columns]

    return targets

# Search for activities in chembl using a target_chembl_id and return a pandas DataFrame with the results
def get_activities(
    target_chembl_id,
    activity_type="IC50",
    relation="=",
    assay_type="B",
):

    activities = get_json(
        "activity",
        {
            "target_chembl_id": target_chembl_id,
            "type": activity_type,
            "relation": relation,
            "assay_type": assay_type,
        },
    )

    activities = pd.DataFrame(activities)

    if not activities.empty:
        columns = [
            "activity_id",
            "assay_chembl_id",
            "assay_description",
            "assay_type",
            "molecule_chembl_id",
            "type",
            "standard_units",
            "relation",
            "standard_value",
            "target_chembl_id",
            "target_organism",
        ]

        available = [
            c for c in columns
            if c in activities.columns
        ]

        activities = activities[available]

    return activities


# Search for compounds in chembl using a list of molecule_chembl_ids and return a pandas DataFrame with the results
def get_compounds(molecule_chembl_ids, batch_size=100):

    molecule_chembl_ids = list(set(molecule_chembl_ids))
    compounds = []

    for start in range(0, len(molecule_chembl_ids), batch_size):
        batch = molecule_chembl_ids[start:start + batch_size]
        ids = ";".join(batch)

        data = get_json(
            f"molecule/set/{ids}"
        )
        for molecule in data:
            compounds.append(
                {
                    "molecule_chembl_id": molecule["molecule_chembl_id"],
                    "molecule_structures": molecule["molecule_structures"],
                }
            )

        print(
            f"Retrieved {len(compounds)} compounds...",
            end="\r"
        )

    return pd.DataFrame(compounds)

# Input of protein of interst, in our case the uniprot_id of D2 receptor, which is P14416. This will be used to query the ChEMBL database for targets and activities related to this protein.
uniprot_id = "P14416"

# Get target information from ChEMBL but restrict it to specified values only
targets = get_targets(uniprot_id)

print(f"The type of targets is {type(targets)}")

# Let's have a look at the target information we retrieved from ChEMBL
targets.head()

target = targets.iloc[0]
target

chembl_id = target.target_chembl_id
print(f"The target ChEMBL ID is {chembl_id}")

bioactivities = get_activities(chembl_id)

print(
    f"Length and type of bioactivities object: "
    f"{len(bioactivities)}, {type(bioactivities)}"
)
