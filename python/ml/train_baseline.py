"""First industry-style ML baseline: scikit-learn GP frequency regression.

Counts are classified separately per C6 sector. Each sector's frequencies
are fitted by rank, allowing a *variable* total number of output modes.
All preprocessing is fitted on TRAIN only. VALIDATION is used for scoring;
TEST is never accessed here.

This is a baseline, not a robust mode tracker through crossings.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.preprocessing import StandardScaler

from .dataset import SECTORS, load_spectrum_records, validate_model_setup

CHANNEL = {0: ("Breathing", "breathing_weight"),
           1: ("CW", "cw_weight"), -1: ("CCW", "ccw_weight")}
WEIGHT_KEYS = ("dipole_weight",)


def sector_rows(record: dict, m: int) -> list[dict]:
    return sorted(
        (mode for mode in record["modes"] if mode["m_signed"] == m),
        key=lambda item: item["frequency"],
    )


def _fit_target(x: np.ndarray, y: np.ndarray):
    """Small-data fallback for slots represented by fewer than three samples."""
    scaler = StandardScaler().fit(x)
    if len(y) < 3 or np.std(y) < 1e-10:
        reg = DummyRegressor(strategy="mean").fit(scaler.transform(x), y)
    else:
        kernel = ConstantKernel(1.0, (1e-2, 1e2)) * Matern(
            length_scale=[1.0, 1.0],
            length_scale_bounds=(1e-2, 1e2),
            nu=1.5,
        ) + WhiteKernel(noise_level=1e-5, noise_level_bounds=(1e-8, 1e-2))
        reg = GaussianProcessRegressor(
            kernel=kernel, normalize_y=True, n_restarts_optimizer=0,
            random_state=42, alpha=1e-8,
        ).fit(scaler.transform(x), y)
    return {"scaler": scaler, "regressor": reg}


def _predict_target(model: dict, x: np.ndarray) -> float:
    return float(model["regressor"].predict(model["scaler"].transform(x))[0])


def train(records: list[dict]) -> dict:
    if len(records) < 5:
        raise ValueError("At least five training configurations are recommended.")
    validate_model_setup(records)
    x = np.array([[r["b"], r["xi"]] for r in records], dtype=float)
    all_sectors: dict[int, dict] = {}
    for m in SECTORS:
        targets = [sector_rows(r, m) for r in records]
        counts = np.array([len(rows) for rows in targets], dtype=int)
        if len(set(counts)) == 1:
            count_model = {"constant": int(counts[0])}
        else:
            classifier = RandomForestClassifier(
                n_estimators=200, min_samples_leaf=2, random_state=42
            ).fit(x, counts)
            count_model = {"classifier": classifier}
        slots: list[dict] = []
        for j in range(int(max(counts, default=0))):
            selected = [i for i, rows in enumerate(targets) if len(rows) > j]
            x_slot = x[selected]
            freq = np.array([targets[i][j]["frequency"] for i in selected])
            weights = {
                key: np.array([targets[i][j][key] for i in selected])
                for key in WEIGHT_KEYS
            }
            slots.append({
                "frequency": _fit_target(x_slot, freq),
                "weights": {
                    key: _fit_target(x_slot, np.log10(np.maximum(value, 1e-15)))
                    for key, value in weights.items()
                },
            })
        all_sectors[m] = {"count_model": count_model, "slots": slots}
    return {
        "schema_version": 1,
        "model": records[0]["model"],
        "omega_max": records[0]["omega_max"],
        "train_samples": len(records),
        "sectors": all_sectors,
    }


def predict(model: dict, b: float, xi: float, min_main_ratio: float = 0.001) -> list[dict]:
    x = np.array([[b, xi]], dtype=float)
    omega_max = float(model["omega_max"])
    predicted = []
    for m, entry in model["sectors"].items():
        count_model = entry["count_model"]
        n_modes = int(count_model["constant"]) if "constant" in count_model else int(
            count_model["classifier"].predict(x)[0]
        )
        for rank, slot in enumerate(entry["slots"][:n_modes]):
            omega = _predict_target(slot["frequency"], x)
            if not (0 < omega <= omega_max):
                continue
            weights = {key: float(10 ** np.clip(_predict_target(fit, x), -15, 10))
                       for key, fit in slot["weights"].items()}
            predicted.append({
                "frequency": omega, "m_signed": m, "m_mod6": m % 6,
                "rank_in_sector": rank, "label": "Multipolar",
                **weights,
            })
    # Assign only dominant observable modes by sector, as in Wolfram export.
    global_max = max((row["dipole_weight"] for row in predicted), default=0)
    for m, (name, key) in CHANNEL.items():
        relevant = [row for row in predicted if row["m_signed"] == m]
        if relevant:
            main = max(relevant, key=lambda row: row["dipole_weight"])
            if global_max > 0 and main["dipole_weight"] >= min_main_ratio * global_max:
                main["label"] = name
    return sorted(predicted, key=lambda row: row["frequency"])


def evaluate(model: dict, records: list[dict]) -> dict:
    """Compare variable-length lists *within* each symmetry sector."""
    frequency_errors = []
    n_true = n_pred = n_matched = 0
    exactly_right = 0
    for reference in records:
        predicted = predict(model, reference["b"], reference["xi"])
        predicted_by_m = {m: sorted((p for p in predicted if p["m_signed"] == m),
                                   key=lambda p: p["frequency"]) for m in SECTORS}
        exact = True
        for m in SECTORS:
            target = sector_rows(reference, m)
            actual = predicted_by_m[m]
            n_true += len(target)
            n_pred += len(actual)
            n = min(len(target), len(actual))
            n_matched += n
            frequency_errors.extend(abs(actual[j]["frequency"] - target[j]["frequency"])
                                    for j in range(n))
            exact &= len(target) == len(actual)
        exactly_right += bool(exact)
    return {
        "n_configurations": len(records),
        "n_true_modes": n_true,
        "n_predicted_modes": n_pred,
        "n_matched_modes": n_matched,
        "mae_frequency_matched_modes": float(np.mean(frequency_errors)) if frequency_errors else None,
        "max_frequency_error_matched_modes": float(np.max(frequency_errors)) if frequency_errors else None,
        "mode_count_precision": n_matched / n_pred if n_pred else 0.,
        "mode_count_recall": n_matched / n_true if n_true else 0.,
        "exact_sector_counts_fraction": exactly_right / len(records) if records else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="datasets/magnon_spectra/gamma")
    parser.add_argument("--output", default="models/magnon_surrogate/baseline.joblib")
    args = parser.parse_args()
    training = load_spectrum_records(args.data_dir, "train")
    validation = load_spectrum_records(args.data_dir, "validation")
    if not validation:
        raise SystemExit("No validation spectra: export the validation split first.")
    validate_model_setup(training + validation)
    if any(any(m["m_signed"] is None for m in row["modes"]) for row in training + validation):
        raise SystemExit("Some modes have unclassified m: inspect their C6 residuals before training.")
    model = train(training)
    report = evaluate(model, validation)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, path)
    path.with_suffix(".validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Saved model: {path}")
    print("Validation metrics (test split untouched):")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()