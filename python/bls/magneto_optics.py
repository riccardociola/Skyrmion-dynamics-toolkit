"""First-order, crossed-polarization BLS coupling (baseline only).

This is NOT yet a complete BLS spectrum. In particular, the second-order
G44 magneto-optic term and Stokes/anti-Stokes occupation factors are not
included. Optical momentum selection must be supplied by the caller through
the reciprocal node; it is not inferred from photon angles here.
"""

import numpy as np
from .io import MagneticModes


def first_order_amplitude(
    data: MagneticModes,
    mode_index: int,
    polarization_in: np.ndarray,
    polarization_out: np.ndarray,
    *,
    node: tuple[int, int] = (0, 0),
) -> complex:
    """Uncalibrated amplitude i(e_out* x e_in).deltaM_G for one mode.

    For linear polarization, the conjugation on e_out is immaterial.
    Complex polarizations are supported. The global magneto-optic constant
    is omitted, and the BdG-normalized mode amplitude is used unchanged.
    """
    e_in = np.asarray(polarization_in, dtype=complex)
    e_out = np.asarray(polarization_out, dtype=complex)
    if e_in.shape != (3,) or e_out.shape != (3,):
        raise ValueError("Polarization vectors must have three Cartesian components")
    delta_m = data.coefficient(mode_index, node)
    return 1j * np.dot(np.cross(np.conj(e_out), e_in), delta_m)


def first_order_intensities(
    data: MagneticModes,
    polarization_in: np.ndarray,
    polarization_out: np.ndarray,
    *,
    node: tuple[int, int] = (0, 0),
) -> np.ndarray:
    """Relative |amplitude|^2 for all modes; no thermal/optical prefactors."""
    return np.asarray([
        abs(first_order_amplitude(data, i, polarization_in,
                                  polarization_out, node=node)) ** 2
        for i in range(len(data.frequencies))
    ])