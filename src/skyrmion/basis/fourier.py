from dataclasses import dataclass
import numpy as np

from skyrmion.basis.symmetry import independent_nodes


@dataclass
class FourierBasis:
    """
    Represent the reciprocal-space Fourier basis used in the calculation.

    The class stores:

    - all reciprocal-space index pairs (n1, n2) inside the square cutoff,
    - the subset of independent nodes after applying the reality condition
      and hexagonal cutoff,
    - the number of independent real parameters,
    - methods to reconstruct the full complex Fourier field.

    For a real vector field F(r), the Fourier components satisfy

        F(-G) = F(G)*

    so only one node from each pair G and -G needs to be stored independently.
    """

    cutoff: int

    def __post_init__(self):
        # Full square grid in affine reciprocal-space coordinates
        self.indices = np.array(
            [
                (n1, n2)
                for n1 in range(-self.cutoff, self.cutoff + 1)
                for n2 in range(-self.cutoff, self.cutoff + 1)
            ],
            dtype=int,
        )

        # Independent nodes inside the physical hexagonal cutoff
        self.independent_indices = independent_nodes(self.cutoff)

    @property
    def size(self):
        """
        Number of nodes in the full square reciprocal-space grid.
        """
        return len(self.indices)

    @property
    def independent_size(self):
        """
        Number of independent reciprocal-space nodes.
        """
        return len(self.independent_indices)

    @property
    def n_parameters(self):
        """
        Number of independent real parameters for a general complex vector field.

        Every non-zero independent Fourier node contains three complex
        components, corresponding to six real parameters.

        At G = 0 the coefficient must be purely real, so only three
        real parameters remain.
        """
        n_nonzero = self.independent_size - 1

        return 6 * n_nonzero + 3

    def unpack(self, x: np.ndarray) -> dict[tuple[int, int], np.ndarray]:
        """
        Convert a real parameter vector into the full complex Fourier field.

        The optimizer works with a one-dimensional real array x.

        For every non-zero independent node G, the six parameters are interpreted
        as

            Re(Fx), Im(Fx),
            Re(Fy), Im(Fy),
            Re(Fz), Im(Fz).

        The Fourier coefficient at -G is then reconstructed using

            F(-G) = F(G)*.

        At G = 0 only three real components are stored because the imaginary
        part must vanish.

        Returns
        -------
        dict
            Dictionary mapping (n1, n2) -> complex vector F_G.
        """

        if len(x) != self.n_parameters:
            raise ValueError(
                f"Expected {self.n_parameters} parameters, got {len(x)}."
            )

        field = {}
        index = 0

        for n1, n2 in self.independent_indices:
            node = (int(n1), int(n2))

            if node == (0, 0):
                coefficient = np.array(
                    [
                        x[index],
                        x[index + 1],
                        x[index + 2],
                    ],
                    dtype=complex,
                )

                field[node] = coefficient

                index += 3

            else:
                coefficient = np.array(
                    [
                        x[index] + 1j * x[index + 1],
                        x[index + 2] + 1j * x[index + 3],
                        x[index + 4] + 1j * x[index + 5],
                    ],
                    dtype=complex,
                )

                field[node] = coefficient
                field[(-node[0], -node[1])] = np.conjugate(coefficient)

                index += 6

        return field