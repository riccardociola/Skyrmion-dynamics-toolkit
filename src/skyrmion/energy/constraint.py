import numpy as np


def quartic_constraint_energy(
    M: dict[tuple[int, int], np.ndarray],
    geometry,
    lam: float,
) -> float:
    """
    Compute the quartic constraint enforcing |M(r)| = 1.
    """

    g_inv = geometry.inverse_metric()

    def dot(v1, v2):
        return v1 @ g_inv @ v2

    # Quadratic term
    quadratic = 0.0 + 0.0j

    for (n1, n2), M_G in M.items():
        M_minus_G = M[(-n1, -n2)]

        quadratic += dot(
            M_G,
            M_minus_G,
        )

    # Quartic term
    quartic = 0.0 + 0.0j

    nodes = list(M.keys())

    for n1, n2 in nodes:
        M_n = M[(n1, n2)]

        for p1, p2 in nodes:
            M_p = M[(p1, p2)]

            M_n_dot_M_p = dot(M_n, M_p)

            for s1, s2 in nodes:

                q = (
                    -n1 - p1 - s1,
                    -n2 - p2 - s2,
                )

                if q not in M:
                    continue

                M_s = M[(s1, s2)]
                M_q = M[q]

                quartic += (
                    M_n_dot_M_p
                    * dot(M_s, M_q)
                )

    energy = lam * (
        1
        - 2 * quadratic
        + quartic
    )

    return float(np.real_if_close(energy))