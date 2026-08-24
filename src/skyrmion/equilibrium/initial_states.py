import numpy as np


def triple_q_initial_guess(
    basis,
    geometry,
    m0: float = 0.5,
    amplitude: float = 0.3,
) -> np.ndarray:
    """
    Construct a simple Bloch-type triple-Q initial state.

    The real-space ansatz is

        m(r) = m0 * z_hat
             + A * sum_n [
                   z_hat cos(Q_n · r)
                   + (z_hat x Q_hat_n) sin(Q_n · r)
               ]

    where the three Q vectors belong to the first reciprocal-lattice shell.

    Parameters
    ----------
    basis:
        FourierBasis used in the calculation.

    geometry:
        AffineGeometry defining the reciprocal lattice.

    m0:
        Uniform magnetization along z.

    amplitude:
        Real-space amplitude A of each helix.

    Returns
    -------
    np.ndarray
        Real parameter vector compatible with FourierBasis.unpack().
    """

    x = np.zeros(basis.n_parameters)

    z_hat = np.array([0.0, 0.0, 1.0])

    # First reciprocal-lattice shell
    first_shell = {
        (1, 0),
        (0, 1),
        (1, 1),
    }

    index = 0

    for n1, n2 in basis.independent_indices:

        node = (int(n1), int(n2))

        # Uniform Fourier component
        if node == (0, 0):

            M0_cartesian = m0 * z_hat

            # Convert to covariant affine components
            M0 = geometry.cartesian_to_covariant(
                M0_cartesian
            )

            x[index:index + 3] = np.real(M0)

            index += 3

        # First-shell helices
        elif node in first_shell:

            Q_cartesian = geometry.reciprocal_vector_cartesian(
                n1,
                n2,
            )

            Q_hat = Q_cartesian / np.linalg.norm(Q_cartesian)

            transverse = np.cross(
                z_hat,
                Q_hat,
            )

            # Fourier coefficient:
            #
            # M_Q = A/2 [ z_hat - i (z_hat x Q_hat) ]
            #
            M_cartesian = (
                amplitude / 2
                * (
                    z_hat
                    - 1j * transverse
                )
            )

            # Our energy functions use covariant affine components
            M = geometry.cartesian_to_covariant(
                M_cartesian
            )

            x[index:index + 6] = [
                M[0].real,
                M[0].imag,
                M[1].real,
                M[1].imag,
                M[2].real,
                M[2].imag,
            ]

            index += 6

        # Higher Fourier shells start from zero
        else:

            index += 6

    return x