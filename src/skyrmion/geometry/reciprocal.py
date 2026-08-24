import numpy as np

def reciprocal_indices(cutoff: int) -> np.ndarray:
    """
    Return all reciprocal-lattice indices (n1, n2)
    with -cutoff <= (n1,n2) <= cutoff.
    """
    return np.array(
        [
            (n1, n2)
            for n1 in range(-cutoff, cutoff + 1)
            for n2 in range(-cutoff, cutoff + 1)
        ],
        dtype=int,
    )