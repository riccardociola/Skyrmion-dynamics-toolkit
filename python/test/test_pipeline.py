"""End-to-end synthetic schema and scikit-learn baseline smoke tests.

Synthetic toy modes validate software, not physical prediction accuracy.
"""

import json
import tempfile
import unittest
from pathlib import Path

from ml.dataset import load_spectrum_records, samples_frame, modes_frame
from ml.train_baseline import train, predict, evaluate


def sample(split, index):
    b = 0.32 + index * 0.03
    xi = 0.75 + index * 0.02
    modes = []
    for m, rank, omega, weight in [
        (0, 0, 0.8 + b, 0.9),
        (1, 0, 1.2 + xi, 0.8),
        (-1, 0, 1.4 + b, 0.7),
        (0, 1, 2.4 + b, 1e-7),
    ]:
        if m == 0 and rank == 1 and index % 2 == 0:
            continue
        modes.append({
            "mode_index": len(modes) + 1,
            "frequency": omega,
            "m_signed": m,
            "m_mod6": m % 6,
            "c6_residual": 1e-9,
            "dipole_weight": weight,
            "breathing_weight": weight if m == 0 else 0.,
            "cw_weight": weight if m == 1 else 0.,
            "ccw_weight": weight if m == -1 else 0.,
            "label": {0:"Breathing",1:"CW",-1:"CCW"}[m] if rank == 0 else "Multipolar",
        })
    return {
        "schema_version": 1,
        "sample_id": f"{split}_{index:03d}",
        "split": split,
        "b": b,
        "xi": xi,
        "omega_max": 3.,
        "model": {"dmi_type": "Bloch", "fourier_cutoff": 4,
                  "magnetic_representation": "SixFold", "cubic_anisotropy": False,
                  "rigid_lattice": False, "lambda": 500,
                  "Nx": 1/3, "Ny": 1/3, "Nz": 1/3},
        "number_of_modes": len(modes), "modes": modes,
    }


class PipelineTests(unittest.TestCase):
    def test_roundtrip_variable_modes_and_training(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            for split, n in (("train", 8), ("validation", 4), ("test", 2)):
                folder = root / split
                folder.mkdir()
                for i in range(n):
                    record = sample(split, i)
                    (folder / f"{record['sample_id']}.json").write_text(
                        json.dumps(record), encoding="utf-8")
            training = load_spectrum_records(root, "train")
            validation = load_spectrum_records(root, "validation")
            self.assertEqual(len(training), 8)
            self.assertEqual(len(samples_frame(training)), 8)
            self.assertIn("rank_in_sector", modes_frame(training).columns)
            self.assertEqual(len(set(r["number_of_modes"] for r in training)), 2)
            model = train(training)
            spectrum = predict(model, 0.4, 0.8)
            self.assertTrue(all(0 < mode["frequency"] < 3. for mode in spectrum))
            self.assertIn("m_signed", spectrum[0])
            report = evaluate(model, validation)
            self.assertEqual(report["n_configurations"], 4)
            self.assertIsNotNone(report["mae_frequency_matched_modes"])

    def test_reject_bad_modulo6(self):
        with tempfile.TemporaryDirectory() as td:
            folder = Path(td) / "train"
            folder.mkdir()
            record = sample("train", 0)
            record["modes"][0]["m_mod6"] = 5
            (folder / "train_000.json").write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaises(ValueError):
                load_spectrum_records(td, "train")


if __name__ == "__main__":
    unittest.main()