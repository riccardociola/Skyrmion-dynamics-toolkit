# Skyrmion Dynamics Toolkit

A modular Mathematica-based toolkit for calculating equilibrium configurations and magnon dynamics of skyrmion lattices.

The current implementation focuses on a hexagonal skyrmion lattice with magnetic interactions only. The code is designed to be extended later to magnetocrystalline anisotropy, magnetoelastic coupling, phonon dynamics, and experimental observables.

## Current functionality

### Equilibrium skyrmion lattice

- affine lattice geometry and reciprocal-space tools;
- Fourier representation of the magnetization;
- sixfold-symmetric and general magnetic Fourier bases;
- magnetic free-energy contributions:
  - exchange;
  - DMI;
  - Zeeman interaction;
  - dipolar interaction;
  - soft-spin normalization constraint;
- symbolic energy construction and caching;
- numerical minimization of the SkL equilibrium;
- saving, loading, and reusing equilibrium configurations;
- parameter sweeps with continuation.

### Magnon dynamics

- construction of the Fourier-space magnetic cross-product operator;
- projection onto the well-resolved transverse magnetic subspace;
- analytical construction of the magnetic Hessian;
- linearized Landau-Lifshitz dynamics;
- magnon eigenfrequencies and BdG-normalized eigenvectors;
- finite-wavevector support in the dynamical matrix.

### Visualization and mode analysis

- real-space SkL plots;
- real and imaginary parts of magnon eigenfunctions;
- dipole classification of zero-wavevector modes as:
  - Breathing;
  - CW;
  - CCW;
  - dipole-inactive;
- angular-momentum classification based on the sixfold rotational symmetry of the SkL.

## Main Wolfram modules

```text
wolfram/
├── Geometry.wl
├── FourierBasis.wl
├── EnergyMagnetic.wl
├── BuildEnergy.wl
├── GroundState.wl
├── DynamicsProjection.wl
├── DynamicsMagnetic.wl
├── Plotting.wl
└── ModeAnalysis.wl
```

## Typical workflow

```mathematica
emag = BuildMagneticEnergy[2, "Bloch"];

equilibrium =
    LoadMagneticEquilibrium[
        "Bloch_N2_b0.5_xi0.88_Nz0.333333.txt"
    ];

dynamics =
    PrepareMagneticDynamics[
        equilibrium
    ];

modes =
    SolveMagneticModes[
        dynamics,
        {0, 0, 0}
    ];

MagnonEigenfunctionPlot[
    dynamics,
    modes,
    2,
    "Part" -> "Both"
];
```

Mode analysis:

```mathematica
dipoleAnalysis =
    AnalyzeDipoleModes[
        dynamics,
        modes
    ];

angularAnalysis =
    AnalyzeAngularMomentumModes[
        dynamics,
        modes
    ];
```

## Roadmap

The next major development step is the inclusion of magnetoelastic coupling.

Planned extensions include:

- elastic degrees of freedom and acoustic phonons;
- coupled magnetic and displacement equations of motion;
- magnon-phonon hybridization;
- non-reciprocal phonon dynamics;
- magnetocrystalline anisotropy and distorted SkL configurations;
- finite-wavevector dispersion calculations and mode tracking;
- experimental observables, including:
  - BLS spectra;
  - microwave selection rules;
  - acoustic and ultrasound response;
  - mode-resolved coupling strengths.

## Project goal

The long-term goal is a reusable framework in which the physical model, equilibrium configuration, dynamics, and observables remain modular, so that new interactions or experimental probes can be added without rebuilding the full workflow.