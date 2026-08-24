import numpy as np

from skyrmion.geometry.affine import AffineGeometry


def test_affine_benchmarks():
    """
    Print some benchmark quantities for the undeformed hexagonal lattice.

    Parameters:
        eta = 1
        phi = 0
        alpha = 0

    For this case, the primitive unit-cell area should be sqrt(3)/2.
    """

    geometry = AffineGeometry(
        eta=1.0,
        phi=0.0,
        alpha=0.0,
    )

    print("\n--- Affine geometry benchmark ---")

    print("\nAffine basis:")
    print(geometry.basis())

    print("\nMetric tensor:")
    print(geometry.metric())

    print("\nInverse metric:")
    print(geometry.inverse_metric())

    print("\nDual basis:")
    print(geometry.dual_basis())

    print("\nMetric determinant:")
    print(geometry.metric_determinant())

    print("\nUnit-cell area:")
    print(geometry.unit_cell_volume())

    print("\nExpected unit-cell area:")
    print(np.sqrt(3) / 2)

    print("\nCheck dual basis e_a · E^b:")
    print(geometry.basis() @ geometry.dual_basis().T)

    print("\nReciprocal vector G(1, 0), affine components:")
    print(geometry.reciprocal_vector_affine(1, 0))

    print("\nReciprocal vector G(1, 0), Cartesian components:")
    print(geometry.reciprocal_vector_cartesian(1, 0))

    print("\nReciprocal vector G(0, 1), Cartesian components:")
    print(geometry.reciprocal_vector_cartesian(0, 1))

    print("\n-------------------------------")

    # Actual automated checks
    assert np.isclose(
        geometry.unit_cell_volume(),
        np.sqrt(3) / 2,
    )

    assert np.allclose(
        geometry.basis() @ geometry.dual_basis().T,
        np.eye(3),
    )

test_affine_benchmarks()