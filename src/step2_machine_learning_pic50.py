"""Regression of pIC50 for dopamine D2 compounds from fingerprints.

Reads the fingerprint CSVs written by step2_convert_to_fingerprint.py,
trains RF, SVR and ANN regressors, and evaluates them with a single train/test split
(predicted vs. measured plots) and k-fold cross-validation.
"""

import random
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import clone
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR

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
    """Return fresh, unfitted regressors."""
    return {
        "RF": RandomForestRegressor(n_estimators=100, random_state=SEED),
        "SVR": SVR(kernel="rbf", C=1.0, gamma="scale"),
        "ANN": MLPRegressor(hidden_layer_sizes=(64, 32), random_state=SEED, max_iter=2000),
    }


def load_fingerprint_data(method):
    """Load X (bits) and y (pIC50) from a fingerprint CSV."""
    df = pd.read_csv(DATA_DIR / f"DopamineD2_{method}.csv")
    bit_cols = [c for c in df.columns if c.startswith("bit_")]
    return df[bit_cols].to_numpy(dtype=float), df["pIC50"].to_numpy()


def regression_metrics(y_true, y_pred):
    """RMSE, MAE, R2 and Pearson r."""
    return {
        "rmse": np.sqrt(mean_squared_error(y_true, y_pred)),
        "mae": mean_absolute_error(y_true, y_pred),
        "r2": r2_score(y_true, y_pred),
        "pearson_r": np.corrcoef(y_true, y_pred)[0, 1],
    }


def single_split_evaluation(method, X, y):
    """Train all models on one split; save predicted-vs-measured plots."""
    train_x, test_x, train_y, test_y = train_test_split(
        X, y, test_size=TEST_SIZE, random_state=SEED
    )
    models = get_models()
    fig, axes = plt.subplots(1, len(models), figsize=(4 * len(models), 4), sharex=True, sharey=True)
    rows = []
    for ax, (name, model) in zip(axes, models.items()):
        model.fit(train_x, train_y)
        pred = model.predict(test_x)
        res = regression_metrics(test_y, pred)
        rows.append({"fingerprint": method, "model": name, **res})
        ax.scatter(test_y, pred, s=12, alpha=0.7)
        lims = [min(y.min(), pred.min()), max(y.max(), pred.max())]
        ax.plot(lims, lims, "r--")
        ax.set_title(f"{name} (R² = {res['r2']:.2f}, RMSE = {res['rmse']:.2f})")
        ax.set_xlabel("Measured pIC50")
    axes[0].set_ylabel("Predicted pIC50")
    fig.suptitle(method)
    fig.savefig(RESULTS_DIR / f"pred_vs_measured_{method}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return rows


def crossvalidation(method, X, y, n_folds=N_FOLDS):
    """Mean and std of each metric over k folds, per model."""
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=SEED)
    rows = []
    for name, model in get_models().items():
        t0 = time.time()
        folds = []
        for train_idx, test_idx in kf.split(X):
            fold_model = clone(model)  # fresh copy per fold
            fold_model.fit(X[train_idx], y[train_idx])
            folds.append(regression_metrics(y[test_idx], fold_model.predict(X[test_idx])))
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
        X, y = load_fingerprint_data(method)
        print(f"{method}: X = {X.shape}, pIC50 range = {y.min():.2f}-{y.max():.2f}")
        split_rows += single_split_evaluation(method, X, y)
        cv_rows += crossvalidation(method, X, y)

    split_df = pd.DataFrame(split_rows)
    cv_df = pd.DataFrame(cv_rows)
    split_df.to_csv(RESULTS_DIR / "reg_results_single_split.csv", index=False)
    cv_df.to_csv(RESULTS_DIR / "reg_results_crossvalidation.csv", index=False)

    pd.set_option("display.width", 200)
    print("\nSingle train/test split")
    print(split_df.round(2).to_string(index=False))
    print(f"\n{N_FOLDS}-fold cross-validation")
    print(cv_df.round(2).to_string(index=False))
