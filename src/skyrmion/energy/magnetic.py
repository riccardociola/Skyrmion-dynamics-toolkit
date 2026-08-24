import numpy as np


def exchange_energy(
    M: dict[tuple[int, int], np.ndarray],
    geometry,
    lattice_constant: float,
) -> float:
    """
    Compute the exchange energy in reciprocal space.

    Parameters
    ----------
    M:
        Dictionary mapping reciprocal-space nodes (n1, n2)
        to complex Fourier coefficients M_G.

    geometry:
        AffineGeometry object.

    lattice_constant:
        Skyrmion lattice constant a.

    Returns
    -------
    float
        Exchange energy.
    """
    energy = 0.0 + 0.0j

    for (n1, n2), M_G in M.items():

        G = geometry.reciprocal_vector_affine(n1, n2)

        M_minus_G = M[(-n1, -n2)]

        G2 = geometry.dot_covariant(G, G)

        MM = geometry.dot_covariant(
            M_G,
            M_minus_G,
        )

        energy += G2 * MM

    energy /= 2 * lattice_constant**2

    return float(np.real_if_close(energy))


def bloch_dmi_energy(
    M: dict[tuple[int, int], np.ndarray],
    geometry,
    lattice_constant: float,
) -> float:
    """
    Compute the Bloch-type Dzyaloshinskii-Moriya energy.
    """
    energy = 0.0 + 0.0j

    volume = geometry.unit_cell_volume()

    for (n1, n2), M_G in M.items():

        G_cart = geometry.reciprocal_vector_cartesian(n1, n2)

        M_G_cart = geometry.covariant_to_cartesian(M_G)
        M_minus_G_cart = geometry.covariant_to_cartesian(
            M[(-n1, -n2)]
        )

        term = np.dot(
            M_G_cart,
            np.cross(G_cart, M_minus_G_cart),
        )

        energy += term

    energy *= -1j / (
        lattice_constant * volume
    )

    return float(np.real_if_close(energy))


def neel_dmi_energy(
    M: dict[tuple[int, int], np.ndarray],
    geometry,
    lattice_constant: float,
) -> float:
    """
    Compute the Néel-type Dzyaloshinskii-Moriya energy.

    This implements the reciprocal-space expression

        E_DMI,N = (2 i / a) sum_G
                  M_z(-G) [ M(G) · G ]

    where the dot product is evaluated using the affine metric.

    Parameters
    ----------
    M:
        Dictionary mapping each reciprocal-space node (n1, n2)
        to the complex covariant Fourier coefficient M_G.

    geometry:
        AffineGeometry object used to compute scalar products
        in the affine basis.

    lattice_constant:
        Skyrmion-lattice constant a.

    Returns
    -------
    float
        Néel DMI energy.
    """
    energy = 0.0 + 0.0j

    for (n1, n2), M_G in M.items():

        G = geometry.reciprocal_vector_affine(n1, n2)
        M_minus_G = M[(-n1, -n2)]

        M_dot_G = geometry.dot_covariant(
            M_G,
            G,
        )

        energy += M_minus_G[2] * M_dot_G

    energy *= 2j / lattice_constant

    return float(np.real_if_close(energy))


def zeeman_energy(
    M: dict[tuple[int, int], np.ndarray],
    field: float,
) -> float:
    """
    Compute the Zeeman energy for a magnetic field along z.
    """
    M0 = M[(0, 0)]

    return float(
        np.real(-field * M0[2])
    )


def dipolar_energy_infinite(
    M: dict[tuple[int, int], np.ndarray],
    geometry,
    xi: float,
    demag_tensor: np.ndarray,
) -> float:
    """
    Compute the magnetostatic (dipolar) energy for an infinitely extended sample.

    Parameters
    ----------
    M:
        Dictionary mapping each reciprocal-space node (n1, n2)
        to the complex covariant Fourier coefficient M_G.

    geometry:
        AffineGeometry object.

    xi:
        Dimensionless magnetostatic coupling parameter.

    demag_tensor:
        3x3 demagnetization tensor given in Cartesian coordinates.

        For a mathematically infinite sample, use the zero matrix.

    Returns
    -------
    float
        Magnetostatic energy.
    """
    energy = 0.0 + 0.0j

    for (n1, n2), M_G in M.items():

        M_minus_G = M[(-n1, -n2)]

        # Uniform Fourier component
        if n1 == 0 and n2 == 0:

            M0_cart = geometry.covariant_to_cartesian(M_G)

            energy += (
                M0_cart
                @ demag_tensor
                @ M0_cart
            )

        # Finite-G Fourier components
        else:
            G = geometry.reciprocal_vector_affine(n1, n2)

            M_dot_G = geometry.dot_covariant(
                M_G,
                G,
            )

            M_minus_dot_G = geometry.dot_covariant(
                M_minus_G,
                G,
            )

            G2 = geometry.dot_covariant(
                G,
                G,
            )

            energy += (
                M_dot_G
                * M_minus_dot_G
                / G2
            )

    energy *= xi

    return float(np.real_if_close(energy))