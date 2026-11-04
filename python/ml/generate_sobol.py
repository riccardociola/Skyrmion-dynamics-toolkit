#!/usr/bin/env python3
"""Generate reproducible, independent Sobol designs for SkL magnon ML.

Each of the training, validation and test datasets uses its own scrambled
Sobol sequence. All sizes are powers of two, preserving Sobol balance.

Example (run from the project root; illustrative parameter bounds):

    python python/ml/generate_sobol.py \
        --b-min 0.3 --b-max 0.7 \
        --xi-min 0.6 --xi-max 1.1

Requires: numpy, scipy
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import qmc


# Independent scrambles; all counts are powers of two.
SPLITS = {
    "train": (6, 0),       # 2**6 = 64 points
    "validation": (4, 1),  # 2**4 = 16 points
    "test": (5, 2),        # 2**5 = 32 points
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate train/validation/test Sobol meshes for (b, xi)."
    )
    parser.add_argument("--b-min", type=float, required=True)
    parser.add_argument("--b-max", type=float, required=True)
    parser.add_argument("--xi-min", type=float, required=True)
    parser.add_argument("--xi-max", type=float, required=True)
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Base seed. Split seeds are seed, seed+1, seed+2 (default: 42).",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("datasets/magnon_spectra/sobol"),
        help="Output directory relative to the current working directory.",
    )
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            # Preserve high precision when saving the floating-point samples.
            formatted = {
                key: f"{value:.17g}" if isinstance(value, (float, np.floating)) else value
                for key, value in row.items()
            }
            writer.writerow(formatted)


def main() -> None:
    args = parse_args()
    if not (np.isfinite([args.b_min, args.b_max, args.xi_min, args.xi_max]).all()):
        raise SystemExit("All parameter bounds must be finite.")
    if args.b_min >= args.b_max or args.xi_min >= args.xi_max:
        raise SystemExit("Require b_min < b_max and xi_min < xi_max.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    lower = [args.b_min, args.xi_min]
    upper = [args.b_max, args.xi_max]
    all_rows: list[dict[str, object]] = []
    metadata = {
        "method": "independent scrambled Sobol sequences",
        "dimensions": ["b", "xi"],
        "bounds": {"b": [args.b_min, args.b_max], "xi": [args.xi_min, args.xi_max]},
        "base_seed": args.seed,
        "splits": {},
    }

    for split, (exponent, seed_offset) in SPLITS.items():
        split_seed = args.seed + seed_offset
        sampler = qmc.Sobol(d=2, scramble=True, seed=split_seed)
        samples = qmc.scale(sampler.random_base2(m=exponent), lower, upper)

        rows = [
            {"sample_id": f"{split}_{i:03d}", "b": float(b), "xi": float(xi)}
            for i, (b, xi) in enumerate(samples)
        ]
        name = f"{split}_{len(rows)}.csv"
        write_csv(args.output_dir / name, rows, ["sample_id", "b", "xi"])
        all_rows.extend({"split": split, **row} for row in rows)
        metadata["splits"][split] = {
            "count": len(rows), "scramble_seed": split_seed, "file": name
        }
        print(f"{split:10} {len(rows):3} points (seed={split_seed}) -> {name}")

    # These should be disjoint with independent random scrambles. Guard against
    # an accidental collision before passing the samples to expensive solvers.
    coordinates = [(row["b"], row["xi"]) for row in all_rows]
    if len(coordinates) != len(set(coordinates)):
        raise RuntimeError("Duplicate (b, xi) points found across splits.")

    write_csv(
        args.output_dir / "all_samples.csv",
        all_rows,
        ["split", "sample_id", "b", "xi"],
    )
    with (args.output_dir / "metadata.json").open("w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)
        f.write("\n")

    print(f"Saved {len(all_rows)} distinct points in: {args.output_dir.resolve()}")
    print("Also saved all_samples.csv and metadata.json")


if __name__ == "__main__":
    main()