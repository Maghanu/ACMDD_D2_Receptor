"""Ligand-based active vs inactive classification of dopamine D2 compounds (TeachOpenCADD T007).

Reads the fingerprint CSVs written by step2_convert_to_fingerprint.py,
trains RF, SVM and ANN classifiers, and evaluates them with a single
train/test split (ROC curves) and k-fold cross-validation.
"""

import random
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import clone, metrics, svm
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, recall_score, roc_auc_score
from sklearn.model_selection import KFold, train_test_split
from sklearn.neural_network import MLPClassifier

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
FINGERPRINTS = ("maccs", "morgan2", "morgan3")
SEED = 22
TEST_SIZE = 0.2
N_FOLDS = 5

random.seed(SEED)
np.random.seed(SEED)


def get_models():
    """Return fresh, unfitted models (parameters as in T007)."""
    return {
        "RF": RandomForestClassifier(n_estimators=100, criterion="entropy", random_state=SEED),
        "SVM": CalibratedClassifierCV(
            svm.SVC(kernel="rbf", C=1, gamma=0.1), ensemble=False
        ),
        "ANN": MLPClassifier(hidden_layer_sizes=(5, 3), random_state=SEED, max_iter=1000),
    }


def load_fingerprint_data(method):
    """Load X (bits), y (active label) and pIC50 from a fingerprint CSV."""
    df = pd.read_csv(DATA_DIR / f"DopamineD2_{method}.csv")
    bit_cols = [c for c in df.columns if c.startswith("bit_")]
    return df[bit_cols].to_numpy(), df["active"].to_numpy(), df["pIC50"].to_numpy()


def model_performance(model, test_x, test_y):
    """Accuracy, sensitivity, specificity and AUC on a test set."""
    prob = model.predict_proba(test_x)[:, 1]
    pred = model.predict(test_x)
    return {
        "accuracy": accuracy_score(test_y, pred),
        "sensitivity": recall_score(test_y, pred),
        "specificity": recall_score(test_y, pred, pos_label=0),
        "auc": roc_auc_score(test_y, prob),
    }


def single_split_evaluation(method, X, y):
    """Train all models on one stratified split; save a ROC plot."""
    train_x, test_x, train_y, test_y = train_test_split(
        X, y, test_size=TEST_SIZE, random_state=SEED, stratify=y
    )
    rows = []
    fig, ax = plt.subplots()
    for name, model in get_models().items():
        model.fit(train_x, train_y)
        rows.append({"fingerprint": method, "model": name, **model_performance(model, test_x, test_y)})
        fpr, tpr, _ = metrics.roc_curve(test_y, model.predict_proba(test_x)[:, 1])
        ax.plot(fpr, tpr, label=f"{name} (AUC = {rows[-1]['auc']:.2f})")
    ax.plot([0, 1], [0, 1], "r--")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(f"ROC curves, {method}")
    ax.legend(loc="lower right")
    fig.savefig(RESULTS_DIR / f"roc_{method}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return rows


def crossvalidation(method, X, y, n_folds=N_FOLDS):
    """Mean and std of each performance measure over k folds, per model."""
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=SEED)
    rows = []
    for name, model in get_models().items():
        t0 = time.time()
        folds = []
        for train_idx, test_idx in kf.split(X):
            fold_model = clone(model)  # fresh copy per fold
            fold_model.fit(X[train_idx], y[train_idx])
            folds.append(model_performance(fold_model, X[test_idx], y[test_idx]))
        fold_df = pd.DataFrame(folds)
        row = {"fingerprint": method, "model": name}
        for col in fold_df:
            row[f"{col}_mean"] = fold_df[col].mean()
            row[f"{col}_std"] = fold_df[col].std(ddof=0)
        row["time_s"] = time.time() - t0
        rows.append(row)
    return rows


if __name__ == "__main__":
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    split_rows, cv_rows = [], []
    for method in FINGERPRINTS:
        X, y, _ = load_fingerprint_data(method)
        print(f"{method}: X = {X.shape}, actives = {int(y.sum())}, inactives = {int(len(y) - y.sum())}")
        split_rows += single_split_evaluation(method, X, y)
        cv_rows += crossvalidation(method, X, y)

    split_df = pd.DataFrame(split_rows)
    cv_df = pd.DataFrame(cv_rows)
    split_df.to_csv(RESULTS_DIR / "ml_results_single_split.csv", index=False)
    cv_df.to_csv(RESULTS_DIR / "ml_results_crossvalidation.csv", index=False)

    pd.set_option("display.width", 200)
    print("\nSingle train/test split")
    print(split_df.round(2).to_string(index=False))
    print(f"\n{N_FOLDS}-fold cross-validation")
    print(cv_df.round(2).to_string(index=False))
