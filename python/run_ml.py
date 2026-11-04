"""CLI for the skyrmion-magnon scikit-learn surrogate.

Run from the project root:
    python python/run_ml.py summary
    python python/run_ml.py train
    python python/run_ml.py predict --b 0.5 --xi 0.88
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib

from ml.dataset import load_spectrum_records, print_summary, validate_model_setup
from ml.train_baseline import train, evaluate, predict

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = PROJECT_ROOT / "datasets" / "magnon_spectra" / "gamma"
DEFAULT_MODEL = PROJECT_ROOT / "models" / "magnon_surrogate" / "baseline.joblib"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    summary = commands.add_parser("summary", help="Show exported-spectra coverage")
    summary.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    fitting = commands.add_parser("train", help="Train on train; evaluate on validation")
    fitting.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    fitting.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    inference = commands.add_parser("predict", help="Predict a Gamma-point mode table")
    inference.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    inference.add_argument("--b", type=float, required=True)
    inference.add_argument("--xi", type=float, required=True)
    args = parser.parse_args()

    if args.command == "summary":
        print_summary(args.data_dir)
    elif args.command == "train":
        training = load_spectrum_records(args.data_dir, "train")
        validation = load_spectrum_records(args.data_dir, "validation")
        if not validation:
            parser.error("No validation JSON files: export validation spectra first")
        validate_model_setup(training + validation)
        if any(any(mode["m_signed"] is None for mode in sample["modes"])
               for sample in training + validation):
            parser.error("Unclassified m values: inspect C6 residuals before training")
        model = train(training)
        metrics = evaluate(model, validation)
        args.model.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(model, args.model)
        report = args.model.with_suffix(".validation.json")
        report.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
        print("Model:", args.model)
        print("Validation results:")
        print(json.dumps(metrics, indent=2))
        print("The test set was not used.")
    elif args.command == "predict":
        model = joblib.load(args.model)  # Only load trusted, locally trained models.
        modes = predict(model, args.b, args.xi)
        print(f"Gamma spectrum: b={args.b}, xi={args.xi}; {len(modes)} predicted modes")
        print(f"{'omega':>10} {'m':>4} {'label':>12} {'dipole':>12}")
        for mode in modes:
            print(f"{mode['frequency']:10.5f} {mode['m_signed']:4d} "
                  f"{mode['label']:>12} {mode['dipole_weight']:12.4e}")


if __name__ == "__main__":
    main()