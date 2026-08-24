import numpy as np


def quartic_constraint_energy(
    M: dict[tuple[int, int], np.ndarray],
    geometry,
    lam: float,
) -> float:
    """
    Compute the quartic constraint enforcing |M(r)| = 1.

    The constraint energy is

        E_constraint = lambda * [
            1
            - 2 * sum_G M(G) · M(-G)
            + sum_{G,P,S}
                (M(G) · M(P))
                (M(S) · M(-G-P-S))
        ]

    where all scalar products are evaluated using the affine metric.

    The Fourier field M is assumed to:
        - contain the full set of reciprocal-space nodes,
        - satisfy M(-G) = M(G)*,
        - contain only nodes allowed by the reciprocal-space cutoff.

    Parameters
    ----------
    M:
        Dictionary mapping reciprocal-space nodes (n1, n2)
        to complex covariant Fourier coefficients.

    geometry:
        AffineGeometry object.

    lam:
        Strength of the normalization constraint.

    Returns
    -------
    float
        Quartic constraint energy.
    """

    # --------------------------------
    # Quadratic contribution
    # --------------------------------

    quadratic = 0.0 + 0.0j

    for node, M_G in M.items():
        n1, n2 = node

        M_minus_G = M[(-n1, -n2)]

        quadratic += geometry.dot_covariant(
            M_G,
            M_minus_G,
        )

    # --------------------------------
    # Quartic contribution
    # --------------------------------

    quartic = 0.0 + 0.0j

    nodes = list(M.keys())

    for n1, n2 in nodes:
        M_n = M[(n1, n2)]

        for p1, p2 in nodes:
            M_p = M[(p1, p2)]

            M_n_dot_M_p = geometry.dot_covariant(
                M_n,
                M_p,
            )

            for s1, s2 in nodes:
                q = (
                    -n1 - p1 - s1,
                    -n2 - p2 - s2,
                )

                # The fourth Fourier component contributes only if
                # it lies inside the retained reciprocal-space basis.
                if q not in M:
                    continue

                M_s = M[(s1, s2)]
                M_q = M[q]

                M_s_dot_M_q = geometry.dot_covariant(
                    M_s,
                    M_q,
                )

                quartic += (
                    M_n_dot_M_p
                    * M_s_dot_M_q
                )

    energy = lam * (
        1
        - 2 * quadratic
        + quartic
    )

    return float(np.real_if_close(energy))