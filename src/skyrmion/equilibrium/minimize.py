from skyrmion.energy.magnetic import (
    exchange_energy,
    bloch_dmi_energy,
    zeeman_energy,
    dipolar_energy_infinite,
)

from skyrmion.energy.constraint import quartic_constraint_energy


def total_magnetic_energy(
    x,
    basis,
    geometry,
    lattice_constant,
    field,
    xi,
    demag_tensor,
    lam,
):
    """
    Total magnetic free energy used for equilibrium minimization.
    """

    M = basis.unpack(x)

    E_ex = exchange_energy(
        M,
        geometry,
        lattice_constant,
    )

    E_dmi = bloch_dmi_energy(
        M,
        geometry,
        lattice_constant,
    )

    E_z = zeeman_energy(
        M,
        field,
    )

    E_ms = dipolar_energy_infinite(
        M,
        geometry,
        xi,
        demag_tensor,
    )

    E_constraint = quartic_constraint_energy(
        M,
        geometry,
        lam,
    )

    return (
        E_ex
        + E_dmi
        + E_z
        + E_ms
        + E_constraint
    )