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
                    + 1j * transverse
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


def mathematica_skl_benchmark(basis):
    """
    Return a converged SkL configuration obtained from the Mathematica code.

    The Mathematica solution uses the restricted representation

        M(G) = (i*mu1, i*mu2, mu3).

    It is converted here to the general complex Fourier representation
    used by FourierBasis.
    """

    if basis.cutoff != 2:
        raise ValueError(
            "This benchmark was calculated with cutoff=2."
        )

    mu = {
        (-1, 1): (
            0.02891578169614092,
            0.02861567803910334,
            -0.03407175246218761,
        ),

        (0, 1): (
            -0.31019634996553976,
            -0.1550981749827669,
            0.26461130640141883,
        ),

        (0, 2): (
            0.005522794656405055,
            0.002761397328202453,
            -0.004297031869920622,
        ),

        (1, 1): (
            0.12602969098995563,
            -0.13661202177123362,
            -0.2159514674472051,
        ),

        (1, 2): (
            -0.028915781696137336,
            -0.0003001036570370284,
            0.03407175246218368,
        ),

        (2, 1): (
            0.0,
            -0.029823163992381357,
            -0.014605161742247315,
        ),

        (2, 2): (
            -0.01484053697016938,
            -0.0006133159675531446,
            0.01716173582622511,
        ),

        (1, 0): (
            0.12602969098995449,
            0.26264171276119223,
            0.21595146744720902,
        ),

        (2, 0): (
            0.014840536970169372,
            0.014227221002616081,
            0.017161735826225036,
        ),
    }

    mu0 = 0.4401195787361723

    x = np.zeros(basis.n_parameters)

    index = 0

    for n1, n2 in basis.independent_indices:

        node = (int(n1), int(n2))

        if node == (0, 0):

            # M(0) = (0, 0, mu3)
            x[index:index + 3] = [
                0.0,
                0.0,
                mu0,
            ]

            index += 3

        else:

            mu1, mu2, mu3 = mu[node]

            # General complex representation:
            #
            # M = (i*mu1, i*mu2, mu3)
            #
            x[index:index + 6] = [
                0.0, mu1,
                0.0, mu2,
                mu3, 0.0,
            ]

            index += 6

    return x