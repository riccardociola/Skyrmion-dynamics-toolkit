import numpy as np


class AffineGeometry:
    """
    Describe the affine coordinate system used for the skyrmion lattice.

    Parameters
    ----------
    eta : float
        Relative length of the second in-plane basis vector.

    phi : float
        Parameter controlling the angle between the two in-plane
        basis vectors.

    alpha : float, optional
        Global rotation of the lattice around the z axis, in radians.

    Notes
    -----
    The affine basis vectors are stored as the rows of a 3x3 matrix.

    The class provides the metric tensor, dual basis, reciprocal-lattice
    vectors, scalar products, and transformations between affine and
    Cartesian coordinates.
    """

    def __init__(self, eta: float, phi: float, alpha: float = 0.0):
        self.eta = eta
        self.phi = phi
        self.alpha = alpha

    def rotation_matrix(self) -> np.ndarray:
        """
        Return the rotation matrix for a rotation around the z axis.
        """
        a = self.alpha

        return np.array([
            [np.cos(a), -np.sin(a), 0.0],
            [np.sin(a),  np.cos(a), 0.0],
            [0.0,        0.0,       1.0],
        ])

    def basis(self) -> np.ndarray:
        """
        Return the affine basis vectors as rows of a 3x3 matrix.

        The rows correspond to e1, e2, and e3.
        """
        theta = np.pi / 6 + self.phi / 2

        basis = np.array([
            [
                np.cos(theta),
                -np.sin(theta),
                0.0,
            ],
            [
                self.eta * np.cos(theta),
                self.eta * np.sin(theta),
                0.0,
            ],
            [
                0.0,
                0.0,
                1.0,
            ],
        ])

        return basis @ self.rotation_matrix()

    def metric(self) -> np.ndarray:
        """
        Return the metric tensor g_ab = e_a · e_b.
        """
        e = self.basis()

        return e @ e.T

    def metric_determinant(self) -> float:
        """
        Return the determinant of the metric tensor.
        """
        return np.linalg.det(self.metric())

    def unit_cell_volume(self) -> float:
        """
        Return the volume of the primitive unit cell.

        For the two-dimensional magnetic lattice this corresponds
        to the area of the in-plane primitive unit cell.
        """
        return np.sqrt(self.metric_determinant())

    def inverse_metric(self) -> np.ndarray:
        """
        Return the inverse metric tensor g^ab.
        """
        return np.linalg.inv(self.metric())

    def dual_basis(self) -> np.ndarray:
        """
        Return the dual affine basis.

        The dual basis E^a is defined by

            e_a · E^b = delta_a^b.
        """
        return np.linalg.inv(self.basis()).T

    def reciprocal_vector_affine(
        self,
        n1: int,
        n2: int,
    ) -> np.ndarray:
        """
        Return the covariant affine components of a reciprocal-lattice vector.

        In affine coordinates,

            G_a = 2*pi*(n1, n2, 0).
        """
        return 2 * np.pi * np.array(
            [n1, n2, 0.0],
            dtype=float,
        )

    def reciprocal_vector_cartesian(
        self,
        n1: int,
        n2: int,
    ) -> np.ndarray:
        """
        Return a reciprocal-lattice vector in Cartesian coordinates.
        """
        G_affine = self.reciprocal_vector_affine(n1, n2)

        return G_affine @ self.dual_basis()

    def dot_covariant(
        self,
        v1: np.ndarray,
        v2: np.ndarray,
    ):
        """
        Compute the scalar product of two vectors represented
        by covariant affine components.

        The contraction is

            v1_a g^ab v2_b.

        No complex conjugation is applied. This is important for
        Fourier-space expressions where the required conjugation
        is handled explicitly through the relation F(-G) = F(G)*.
        """
        return v1 @ self.inverse_metric() @ v2

    def dot_contravariant(
        self,
        v1: np.ndarray,
        v2: np.ndarray,
    ):
        """
        Compute the scalar product of two vectors represented
        by contravariant affine components.

        The contraction is

            v1^a g_ab v2^b.
        """
        return v1 @ self.metric() @ v2

    def covariant_to_cartesian(
        self,
        v: np.ndarray,
    ) -> np.ndarray:
        """
        Convert covariant affine components to Cartesian components.
        """
        return v @ self.dual_basis()

    def contravariant_to_cartesian(
        self,
        v: np.ndarray,
    ) -> np.ndarray:
        """
        Convert contravariant affine components to Cartesian components.
        """
        return v @ self.basis()

    def cartesian_to_covariant(
        self,
        v: np.ndarray,
    ) -> np.ndarray:
        """
        Convert a Cartesian vector to covariant affine components.
        """
        return self.basis() @ v

    def cartesian_to_contravariant(
        self,
        v: np.ndarray,
    ) -> np.ndarray:
        """
        Convert a Cartesian vector to contravariant affine components.
        """
        return self.dual_basis() @ v

    def reciprocal_indices(
        self, 
        cutoff: int,
    ) -> np.ndarray:
        """
        Return all reciprocal-lattice index pairs (n1, n2)
        inside the square affine grid defined by the cutoff.
        """
        return np.array(
            [
                (n1, n2)
                for n1 in range(-cutoff, cutoff + 1)
                for n2 in range(-cutoff, cutoff + 1)
            ],
            dtype=int,
        )