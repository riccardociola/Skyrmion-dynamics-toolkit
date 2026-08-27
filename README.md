# Skyrmion Dynamics Toolkit

Hybrid Mathematica/Python toolkit for constructing, minimizing and analyzing skyrmion-lattice models.

## Current functionality

### Geometry
- Affine lattice basis and dual basis
- Coordinate transformations between affine and Cartesian frames
- Reciprocal lattice vectors

### Fourier basis
- Hexagonal reciprocal-space cutoff
- Reality constraints
- Six-fold magnetic representation
- General magnetic and elastic representations
- Automatic generation of independent Fourier variables

### Magnetic energy
- Exchange interaction
- Bloch, Néel and antiskyrmion DMI
- Zeeman energy
- Magnetostatic energy
- Magnetization-length constraint

### Energy construction
- Automatic assembly of the magnetic energy
- Symbolic simplification
- Cached symbolic energies for reuse

### Ground-state calculation
- Newton minimization with arbitrary physical parameters
- User-defined initial ansatz
- Reuse of previous equilibria as initial states
- Automatic saving/loading of equilibrium configurations

### Visualization
- Fourier-to-real-space reconstruction
- Skyrmion configuration plots

## Next steps

- Standardize equilibrium result format
- Improve automatic continuation in magnetic field / model parameters
- Implement general-symmetry ground-state calculations
- Add magnetic anisotropy
- Add elasticity and magnetoelastic coupling
- Implement linearized dynamics
- Add observables such as BLS spectra and phonon velocities
- Build the Python interface and Wolfram bridge
