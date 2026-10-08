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

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

