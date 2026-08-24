from dataclasses import dataclass
import numpy as np

@dataclass
class FourierBasis:
    """
    Represent the set of reciprocal-lattice Fourier modes used in the calculation.
    The basis contains all integer index pairs (n1, n2) within the cutoff range
        -cutoff <= n1, n2 <= cutoff.

    Each pair (n1, n2) labels one reciprocal-lattice vector G(n1, n2).
    For example, cutoff=5 gives an 11 x 11 grid, corresponding to 121
    reciprocal-space modes.

    This class only stores the reciprocal-space indices. The actual reciprocal
    vectors, Fourier coefficients, and symmetry reductions are handled elsewhere.
    """
    cutoff: int

    def __post_init__(self):
        self.indices = np.array(
            [
                (n1, n2)
                for n1 in range(-self.cutoff, self.cutoff + 1)
                for n2 in range(-self.cutoff, self.cutoff + 1)
            ],
            dtype=int,
        )

    @property
    def size(self):
        return len(self.indices)