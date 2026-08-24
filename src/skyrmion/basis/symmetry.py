import numpy as np


def is_in_hexagon(n1: int, n2: int, cutoff: int) -> bool:
    """
    Check whether a reciprocal-space node lies inside the hexagonal cutoff.

    The reciprocal lattice is represented using affine integer coordinates
    (n1, n2). In these coordinates, a naive finite cutoff produces a square:

        -cutoff <= n1 <= cutoff
        -cutoff <= n2 <= cutoff

    The physical reciprocal-space domain is instead hexagonal. The additional
    condition

        |n2 - n1| <= cutoff

    removes the two triangular corners of the square.

    This is necessary at finite Fourier cutoff to avoid artificial modes that
    do not belong to the intended hexagonal reciprocal-space domain.
    """
    return (
        abs(n1) <= cutoff
        and abs(n2) <= cutoff
        and abs(n2 - n1) <= cutoff
    )


def independent_nodes(cutoff: int) -> np.ndarray:
    """
    Return one representative from each pair of reciprocal-space nodes
    G and -G.

    For a real field, the Fourier coefficients satisfy

        F(-G) = F(G)*

    so the coefficients at G and -G are not independent.

    We choose as independent nodes:

        n2 > 0

    together with

        n2 == 0 and n1 >= 0.

    Only nodes inside the hexagonal cutoff are retained.
    """
    nodes = []

    for n1 in range(-cutoff, cutoff + 1):
        for n2 in range(-cutoff, cutoff + 1):

            if not is_in_hexagon(n1, n2, cutoff):
                continue

            if n2 > 0 or (n2 == 0 and n1 >= 0):
                nodes.append((n1, n2))

    return np.array(nodes, dtype=int)


def opposite_node(n1: int, n2: int) -> tuple[int, int]:
    """
    Return the reciprocal-space node corresponding to -G.
    """
    return -n1, -n2


def enforce_reality(coefficient: np.ndarray) -> np.ndarray:
    """
    Return the Fourier coefficient required at -G by the reality condition.

    For any real vector field F(r),

        F(-G) = F(G)*.

    Therefore the real part is even under G -> -G and the imaginary
    part is odd.
    """
    return np.conjugate(coefficient)