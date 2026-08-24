from skyrmion.basis.fourier import FourierBasis
from skyrmion.geometry.affine import AffineGeometry
from skyrmion.equilibrium.initial_states import triple_q_initial_guess


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