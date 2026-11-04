"""Load variable-length Γ-point mode spectra exported by MLModeDataset.wl.

No padding or mode truncation is applied when reading data. The dataset
keeps the independent Sobol train/validation/test identities intact.
"""

from __future__ import annotations

import json
from pathlib import Path
import pandas as pd

SECTORS = (-2, -1, 0, 1, 2, 3)
SPLITS = ("train", "validation", "test")
REQUIRED_MODE_KEYS = {
    "mode_index", "frequency", "m_signed", "m_mod6", "c6_residual",
    "dipole_weight", "breathing_weight", "cw_weight", "ccw_weight", "label",
}


def load_spectrum_records(root: str | Path, split: str | None = None) -> list[dict]:
    """Read per-sample JSON records from ``datasets/magnon_spectra/gamma``.

    Parameters
    ----------
    root : directory
        Parent folder containing train/, validation/, test/.
    split : str | None
        One split or all three (if None).
    """
    folder = Path(root)
    if split is not None and split not in SPLITS:
        raise ValueError(f"Unknown split: {split}")
    splits = (split,) if split else SPLITS
    records: list[dict] = []
    seen: set[str] = set()
    for part in splits:
        for path in sorted((folder / part).glob("*.json")):
            with path.open("r", encoding="utf-8") as handle:
                record = json.load(handle)
            if record.get("schema_version") != 1:
                raise ValueError(f"Unsupported schema in {path}")
            if record.get("split") != part:
                raise ValueError(f"Split mismatch in {path}")
            ident = record.get("sample_id")
            if not isinstance(ident, str) or not ident.startswith(part + "_"):
                raise ValueError(f"Invalid sample_id in {path}")
            if ident in seen:
                raise ValueError(f"Duplicate sample_id: {ident}")
            seen.add(ident)
            if record.get("number_of_modes") != len(record.get("modes", [])):
                raise ValueError(f"Mode count mismatch: {path}")
            for mode in record["modes"]:
                if not REQUIRED_MODE_KEYS.issubset(mode):
                    raise ValueError(f"Missing mode fields in {path}")
                m = mode["m_signed"]
                if m is not None and m not in SECTORS:
                    raise ValueError(f"Invalid angular momentum {m} in {path}")
                if m is not None and mode["m_mod6"] != m % 6:
                    raise ValueError(f"Inconsistent m modulo 6 in {path}")
            records.append(record)
    return records


def samples_frame(records: list[dict]) -> pd.DataFrame:
    """One row per physical configuration; never uses individual modes as samples."""
    return pd.DataFrame.from_records(
        {
            "sample_id": r["sample_id"], "split": r["split"],
            "b": float(r["b"]), "xi": float(r["xi"]),
            "omega_max": float(r["omega_max"]),
            "number_of_modes": r["number_of_modes"],
            "number_unclassified": sum(m["m_signed"] is None for m in r["modes"]),
        }
        for r in records
    )


def modes_frame(records: list[dict]) -> pd.DataFrame:
    """One row per mode. Rank is frequency order *within* each C6 sector.

    Sector rank is an indexing baseline, not guaranteed to identify the same
    physical eigenbranch near avoided crossings or frequency-cutoff boundaries.
    """
    rows = []
    for record in records:
        sorted_modes = sorted(record["modes"], key=lambda m: m["frequency"])
        counts: dict[int, int] = {}
        for mode in sorted_modes:
            m = mode["m_signed"]
            rank = None
            if m is not None:
                rank = counts.get(m, 0)
                counts[m] = rank + 1
            rows.append({
                "sample_id": record["sample_id"],
                "split": record["split"],
                "b": float(record["b"]),
                "xi": float(record["xi"]),
                "omega_max": float(record["omega_max"]),
                "rank_in_sector": rank,
                **mode,
            })
    return pd.DataFrame.from_records(rows)


def validate_model_setup(records: list[dict]) -> None:
    """Reject training data mixing different Hamiltonians/cutoff settings."""
    if not records:
        raise ValueError("No mode records found.")
    reference = records[0]
    for record in records[1:]:
        if record["model"] != reference["model"]:
            raise ValueError(f"Model settings differ for {record['sample_id']}")
        if abs(record["omega_max"] - reference["omega_max"]) > 1e-10:
            raise ValueError("Different frequency cutoffs in one dataset")


def print_summary(root: str | Path) -> None:
    records = load_spectrum_records(root)
    if not records:
        print("No spectra exported yet.")
        return
    samples = samples_frame(records)
    modes = modes_frame(records)
    print("\nConfigurations by split:")
    print(samples.groupby("split", sort=False).agg(
        configurations=("sample_id", "count"),
        min_modes=("number_of_modes", "min"),
        max_modes=("number_of_modes", "max"),
        unclassified=("number_unclassified", "sum"),
    ).to_string())
    if not modes.empty:
        print("\nClassification counts:")
        print(modes.groupby(["split", "label"], dropna=False).size().to_string())


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="datasets/magnon_spectra/gamma")
    args = parser.parse_args()
    print_summary(args.data_dir)