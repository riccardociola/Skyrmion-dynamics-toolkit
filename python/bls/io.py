"""Read full Cartesian Fourier magnon eigenmodes exported by Mathematica.

All Fourier arrays have the *same* node ordering from DynamicsNodes[N].
Modes are BdG-normalized in Mathematica; this module never rescales them.
"""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class MagneticModes:
    """Numerical data for one Cartesian Bloch wavevector.

    equilibrium_m: complex ndarray, shape (n_nodes, 3)
    eigenvectors: complex ndarray, shape (n_modes, n_nodes, 3)
    frequencies: dimensionless positive frequencies, shape (n_modes,)
    """

    nodes: np.ndarray
    reciprocal_vectors: np.ndarray
    wavevector: np.ndarray
    equilibrium_m: np.ndarray
    eigenvectors: np.ndarray
    frequencies: np.ndarray
    bdg_norms: np.ndarray
    legacy_norms: np.ndarray
    geometry: dict
    parameters: dict
    model: dict

    def node_index(self, n1: int, n2: int) -> int:
        """Find the Fourier-array index corresponding to G[n1,n2]."""
        matches = np.flatnonzero(
            (self.nodes[:, 0] == n1) & (self.nodes[:, 1] == n2)
        )
        if len(matches) != 1:
            raise KeyError(f"Reciprocal node {(n1, n2)} not in exported basis")
        return int(matches[0])

    def coefficient(self, mode_index: int, node: tuple[int, int] = (0, 0)) -> np.ndarray:
        """Complex Cartesian delta-M coefficient of one BdG-normalized mode."""
        return self.eigenvectors[mode_index, self.node_index(*node)]


def _complex_array(source: dict) -> np.ndarray:
    real = np.asarray(source["real"], dtype=float)
    imag = np.asarray(source["imag"], dtype=float)
    if real.shape != imag.shape:
        raise ValueError("Complex real/imaginary arrays have different shapes")
    return real + 1j * imag


def load_magnetic_modes(path: str | Path, *, norm_tolerance: float = 1e-4) -> MagneticModes:
    """Load and validate a single-k Fourier mode dataset.

    norm_tolerance checks |BdG norm| == 1, accepting either sign. The
    comparison with the legacy normFT is a diagnostic, not a renormalization.
    """
    with Path(path).open(encoding="utf-8") as stream:
        payload = json.load(stream)

    if payload.get("schema") != "skyrmion-magnetic-modes-v1":
        raise ValueError("Unknown magnon export schema")
    if payload.get("frame") != "simulation-cartesian":
        raise ValueError("Expected Cartesian, not affine, Fourier coefficients")
    if payload.get("fourier_convention") != "exp(+i (k+G).r)":
        raise ValueError("Unrecognized Fourier phase convention")

    nodes = np.asarray(payload["nodes"], dtype=int)
    reciprocal_vectors = np.asarray(payload["reciprocal_vectors_cartesian"], dtype=float)
    wavevector = np.asarray(payload["wavevector_cartesian"], dtype=float)
    equilibrium_m = _complex_array(payload["equilibrium"]["magnetization"])
    eigenvectors = _complex_array(payload["modes"]["eigenvectors"])
    frequencies = np.asarray(payload["modes"]["frequencies"], dtype=float)
    norms = np.asarray(payload["normalization"]["bdg_norms"], dtype=float)
    legacy = np.asarray(payload["normalization"]["legacy_normFT"], dtype=float)

    n_nodes = len(nodes)
    n_modes = len(frequencies)
    if nodes.shape != (n_nodes, 2):
        raise ValueError("Fourier indices must have shape (n_nodes, 2)")
    if reciprocal_vectors.shape != (n_nodes, 3) or wavevector.shape != (3,):
        raise ValueError("Reciprocal vectors/wavevector have incorrect dimensions")
    if equilibrium_m.shape != (n_nodes, 3):
        raise ValueError("Equilibrium Fourier magnetization has incorrect dimensions")
    if eigenvectors.shape != (n_modes, n_nodes, 3):
        raise ValueError("Eigenvector array must have shape (n_modes, n_nodes, 3)")
    if norms.shape != (n_modes,) or legacy.shape != (n_modes,):
        raise ValueError("Mode normalization arrays have incorrect dimensions")
    for value in (reciprocal_vectors, wavevector, equilibrium_m,
                  eigenvectors, frequencies, norms, legacy):
        if not np.all(np.isfinite(value)):
            raise ValueError("Export contains non-finite values")
    if n_modes and not np.allclose(np.abs(norms), 1., atol=norm_tolerance, rtol=0.):
        raise ValueError("Eigenvectors are not BdG-unit-normalized")

    return MagneticModes(
        nodes=nodes, reciprocal_vectors=reciprocal_vectors,
        wavevector=wavevector, equilibrium_m=equilibrium_m,
        eigenvectors=eigenvectors, frequencies=frequencies,
        bdg_norms=norms, legacy_norms=legacy,
        geometry=payload["equilibrium"]["geometry"],
        parameters=payload["equilibrium"]["parameters"],
        model=payload["equilibrium"]["model"],
    )