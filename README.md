# Skyrmion Dynamics Toolkit

This project builds on the theoretical research and numerical code I developed throughout my PhD, investigating the equilibrium properties and collective excitations of magnetic skyrmion lattices. It reorganizes that work into a modular **Wolfram Language (Mathematica)** toolkit for studying skyrmion-lattice energetics, magnon dynamics, and ultimately magnon–phonon coupling.

## Project structure

```text
wolfram/
├── Geometry.wl                   # Affine lattice geometry and coordinate transformations
├── FourierBasis.wl               # SixFold and General Fourier representations
├── EnergyMagnetic.wl             # Exchange, DMI, Zeeman, dipolar and constraint energies
├── EnergyAnisotropy.wl           # Cubic magnetocrystalline anisotropy (MCA)
├── BuildEnergy.wl                # Symbolic energy construction and caching
├── GroundState.wl                # Equilibrium minimization and parameter sweeps
├── DynamicsProjection.wl         # Projection onto physical transverse modes
├── DynamicsMagnetic.wl           # Isotropic magnon dynamics
├── DynamicsMagneticAnisotropy.wl # Magnon dynamics including MCA
├── ModeAnalysis.wl               # Eigenmode characterization
├── Dispersion.wl                 # Longitudinal and transverse magnon spectra
└── Plotting.wl                   # Equilibrium textures and eigenmode visualization

cache/energies/                 # Cached symbolic expressions
results/equilibria/            # Saved equilibrium configurations
notebooks/                     # Calculation and analysis notebooks
legacy/                        # Original research code
```

## Current status

- **Magnetic equilibrium:** Fourier-space free energies, affine lattice deformation, and numerical minimization for both restricted (`SixFold`) and unrestricted (`General`) magnetization representations.
- **Cubic anisotropy:** Field-direction-dependent MCA, in-plane lattice rotation, and rigid-versus-deformable lattice minimization; parameter sweeps can probe the transition toward an elongated skyrmion lattice (eSkL).
- **Magnon dynamics:** Isotropic Landau–Lifshitz dynamics with projection onto physical transverse modes, eigenmode analysis, and dispersion calculations.
- **MCA dynamics:** An initial implementation of the cubic-anisotropy Hessian and projected dynamical operator is available; numerical benchmarking and validation are ongoing.
- **Visualization:** Plots of equilibrium textures and magnon eigenfunctions, including deformed lattices.

## Planned development

1. **Validate anisotropic dynamics:** Test convergence and spectra across weak and strong MCA, including elongated skyrmion lattices and symmetry-broken modes.
2. **Introduce magnetoelastic coupling:** Add elastic degrees of freedom and derive coupled magnon–phonon dynamics, hybridization, and selection rules.
3. **Connect theory to experiments:** Develop observables relevant to Brillouin light scattering and ultrasound spectroscopy, including nonreciprocal propagation.
4. **Improve usability:** Expand automated tests, examples, documentation, and potentially a Python interface or port.

> **Status:** Active research code. The isotropic framework is established; the MCA extension is under validation, and magnetoelastic coupling is planned.