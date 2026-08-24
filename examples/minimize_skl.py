import numpy as np

from scipy.optimize import minimize

from skyrmion.geometry.affine import AffineGeometry
from skyrmion.basis.fourier import FourierBasis
from skyrmion.equilibrium.initial_states import triple_q_initial_guess

from skyrmion.energy.magnetic import (
    exchange_energy,
    bloch_dmi_energy,
    zeeman_energy,
    dipolar_energy_infinite,
)

from skyrmion.energy.constraint import quartic_constraint_energy

basis = FourierBasis(cutoff=1)

geometry = AffineGeometry(
    eta=1.0,
    phi=0.0,
    alpha=0.0,
)

x0 = triple_q_initial_guess(
    basis=basis,
    geometry=geometry,
    m0=0.5,
    amplitude=0.3,
)

M0 = basis.unpack(x0)


# Model parameters
a = 2 * np.pi
b = 0.5
xi = 0.65
lam = 100.0

N = np.diag([1 / 3, 1 / 3, 1 / 3])  # Demagnetization tensor for a spherical sample


print("\n--- Initial energy ---")

print("Exchange:", exchange_energy(M0, geometry, a))
print("Bloch DMI:", bloch_dmi_energy(M0, geometry, a))
print("Zeeman:", zeeman_energy(M0, b))
print("Dipolar:", dipolar_energy_infinite(M0, geometry, xi, N))
print("Constraint:", quartic_constraint_energy(M0, geometry, lam))

def total_energy(x):
    M = basis.unpack(x)

    return (
        exchange_energy(M, geometry, a)
        + bloch_dmi_energy(M, geometry, a)
        + zeeman_energy(M, b)
        + dipolar_energy_infinite(M, geometry, xi, N)
        + quartic_constraint_energy(M, geometry, lam)
    )

result = minimize(
    total_energy,
    x0,
    method="BFGS",
    options={
        "disp": True,
        "maxiter": 1000,
    },
)

print("\n--- Minimization result ---")
print("Success:", result.success)
print("Message:", result.message)
print("Initial energy:", total_energy(x0))
print("Final energy:", result.fun)

M_eq = basis.unpack(result.x)

print("\n--- Final energy contributions ---")

print("Exchange:", exchange_energy(M_eq, geometry, a))
print("Bloch DMI:", bloch_dmi_energy(M_eq, geometry, a))
print("Zeeman:", zeeman_energy(M_eq, b))
print("Dipolar:", dipolar_energy_infinite(M_eq, geometry, xi, N))
print("Constraint:", quartic_constraint_energy(M_eq, geometry, lam))

print("\n--- Final Fourier coefficients ---")

for node, coefficient in M_eq.items():
    print(node, coefficient)

print("\nGradient norm:")
print(np.linalg.norm(result.jac))