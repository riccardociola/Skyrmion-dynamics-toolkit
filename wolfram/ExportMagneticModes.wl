(* ================================================================ *)
(* ExportMagneticModes.wl                                           *)
(* Export BdG-normalized magnetic modes for Python BLS calculations. *)
(*                                                                  *)
(* Dependencies (load first): Geometry.wl, FourierBasis.wl,        *)
(* GroundState.wl, DynamicsProjection.wl, DynamicsMagnetic.wl       *)
(*                                                                  *)
(* Convention: full Cartesian Fourier coefficients, ordered by      *)
(* dynamics["Nodes"], and exp[I (k+G).r] spatial dependence.        *)
(* Frequencies and Cartesian wavevectors are in the solver units.  *)
(* ================================================================ *)

ClearAll[
    ExportMagneticModes,
    ComputeAndExportMagneticModes,
    ExportMagneticModesHelp,
    blsExportComplexArray,
    blsExportScalar,
    blsExportArray
];

ExportMagneticModes::data =
    "Invalid dynamics/mode data: `1`.";
ExportMagneticModes::numeric =
    "The exported `1` contains unresolved or nonnumeric values.";
ExportMagneticModes::norm =
    "The exported modes are not all BdG-unit-normalized (maximum deviation from |norm|=1: `1`). No renormalization was applied.";
ExportMagneticModes::write =
    "Could not export the mode data to `1`.";

(* Keep all complex arrays as separate real and imaginary arrays in JSON. *)
blsExportComplexArray[array_] := <|
    "real" -> N[Re[array]],
    "imag" -> N[Im[array]]
|>;

blsExportScalar[rules_List, symbol_] := Module[{value},
    value = symbol /. rules;
    If[NumericQ[value], N[value], Null]
];

blsExportArray[array_] := N[array];

Options[ExportMagneticModes] = {
    "NormalizationTolerance" -> 10^-5,
    "CreateDirectories" -> True
};

(* Export already-solved modes. Does not recompute or renormalize them. *)
ExportMagneticModes[
    dynamics_Association,
    modes_Association,
    file_String,
    OptionsPattern[]
] := Module[
    {
        equilibrium, rules, nodes, numNodes,
        fullVectors, transverseVectors, numModes,
        k, frequencies, geometry, etaValue, phiValue,
        alphaValue, aValue, m0, reciprocalVectors,
        bdgNorms, legacyNorms, maxDeviation,
        model, parameters, output, parent, result
    },

    If[!AllTrue[
            {"Equilibrium", "Nodes", "TransverseData", "CrossOperatorFull"},
            KeyExistsQ[dynamics, #] &
        ] ||
       !AllTrue[
            {"Frequencies", "EigenvectorsFull", "EigenvectorsTransverse", "WaveVector"},
            KeyExistsQ[modes, #] &
        ],
        Message[ExportMagneticModes::data, "missing association keys"];
        Return[$Failed]
    ];

    equilibrium = dynamics["Equilibrium"];
    rules = EquilibriumRules[equilibrium];
    nodes = dynamics["Nodes"];
    numNodes = Length[nodes];
    fullVectors = modes["EigenvectorsFull"];
    transverseVectors = modes["EigenvectorsTransverse"];
    frequencies = modes["Frequencies"];
    numModes = Length[frequencies];
    k = modes["WaveVector"];

    If[numModes == 0 || !ListQ[nodes] || !MatrixQ[nodes, IntegerQ] ||
       Length[fullVectors] != numModes ||
       Length[transverseVectors] != numModes ||
       !AllTrue[fullVectors, VectorQ[#, NumericQ] && Length[#] == 3 numNodes &] ||
       !VectorQ[frequencies, NumericQ] ||
       !(VectorQ[k, NumericQ] && Length[k] == 3),
        Message[ExportMagneticModes::data, "inconsistent node/mode dimensions"];
        Return[$Failed]
    ];

    aValue = blsExportScalar[rules, a0];
    etaValue = blsExportScalar[rules, eta0];
    phiValue = blsExportScalar[rules, phi0];
    alphaValue = blsExportScalar[rules, alphaR0];
    If[MemberQ[{aValue, etaValue, phiValue, alphaValue}, Null],
        Message[ExportMagneticModes::numeric, "affine geometry"];
        Return[$Failed]
    ];

    (* M0 is stored internally as covariant affine Fourier components;
       full dynamical eigenvectors already use Cartesian components. *)
    m0 = Table[
        toCRT[
            EquilibriumMagnetization[equilibrium, Sequence @@ node],
            etaValue, phiValue, alphaValue
        ],
        {node, nodes}
    ];
    reciprocalVectors = Table[
        toCRT[G[Sequence @@ node], etaValue, phiValue, alphaValue] / aValue,
        {node, nodes}
    ];

    If[!MatrixQ[m0, NumericQ] || !MatrixQ[reciprocalVectors, NumericQ],
        Message[ExportMagneticModes::numeric, "magnetization or reciprocal vectors"];
        Return[$Failed]
    ];

    bdgNorms = Re /@ Table[
        MagneticBdGNorm[dynamics["TransverseData"], u],
        {u, transverseVectors}
    ];
    legacyNorms = Re /@ Table[
        I Conjugate[u] . dynamics["CrossOperatorFull"] . u,
        {u, fullVectors}
    ];
    If[!VectorQ[bdgNorms, NumericQ] || !VectorQ[legacyNorms, NumericQ],
        Message[ExportMagneticModes::numeric, "normalization diagnostics"];
        Return[$Failed]
    ];

    maxDeviation = If[numModes == 0, 0,
        Max[Abs[Abs[bdgNorms] - 1]]
    ];
    If[maxDeviation > OptionValue["NormalizationTolerance"],
        Message[ExportMagneticModes::norm, maxDeviation]
    ];

    geometry = <|
        "a0" -> aValue,
        "eta0" -> etaValue,
        "phi0" -> phiValue,
        "alphaR0" -> alphaValue
    |>;

    parameters = Association[
        (#[[1]] -> blsExportScalar[rules, #[[2]]]) & /@ {
            {"b", b}, {"xi", xi}, {"lambda", lambda},
            {"Nx", Nx}, {"Ny", Ny}, {"Nz", Nz},
            {"kappaCub", kappaCub}
        }
    ];

    model = <|
        "dmi_type" -> Lookup[equilibrium, "DMIType", "Unknown"],
        "fourier_cutoff" -> equilibrium["FourierCutoff"],
        "magnetic_representation" -> equilibrium["MagneticRepresentation"],
        "cubic_anisotropy" -> TrueQ[Lookup[equilibrium, "CubicAnisotropy", False]],
        "rigid_lattice" -> TrueQ[Lookup[equilibrium, "RigidLattice", False]]
    |>;

    output = <|
        "schema" -> "skyrmion-magnetic-modes-v1",
        "frame" -> "simulation-cartesian",
        "fourier_convention" -> "exp(+i (k+G).r)",
        "frequency_units" -> "dimensionless-solver",
        "wavevector_units" -> "inverse-dimensionless-length",
        "wavevector_cartesian" -> N[k],
        "nodes" -> nodes,
        "reciprocal_vectors_cartesian" -> N[reciprocalVectors],
        "equilibrium" -> <|
            "model" -> model,
            "parameters" -> parameters,
            "geometry" -> geometry,
            "magnetization" -> blsExportComplexArray[m0]
        |>,
        "modes" -> <|
            "frequencies" -> N[frequencies],
            "eigenvectors" -> blsExportComplexArray[
                Partition[#, 3] & /@ fullVectors
            ]
        |>,
        "normalization" -> <|
            "definition" -> "-i u_perp^dagger C_perp^(-1) u_perp",
            "legacy_normFT_definition" -> "real(i u_full^dagger C_full u_full)",
            "bdg_norms" -> N[bdgNorms],
            "legacy_normFT" -> N[legacyNorms],
            "max_abs_bdg_unit_norm_error" -> N[maxDeviation],
            "vectors_pre_normalized" -> True
        |>
    |>;

    parent = DirectoryName[ExpandFileName[file]];
    If[TrueQ[OptionValue["CreateDirectories"]] && !DirectoryQ[parent],
        CreateDirectory[parent, CreateIntermediateDirectories -> True]
    ];
    result = Quiet @ Check[
        Export[file, output, "RawJSON"],
        $Failed
    ];
    If[result === $Failed,
        Message[ExportMagneticModes::write, file];
        Return[$Failed]
    ];
    result
];

Options[ComputeAndExportMagneticModes] = Join[
    Options[SolveMagneticModes], Options[ExportMagneticModes]
];

(* Solve a single wavevector and immediately export the normalized modes. *)
ComputeAndExportMagneticModes[
    dynamics_Association,
    k_List,
    file_String,
    OptionsPattern[]
] := Module[{modes},
    modes = SolveMagneticModes[
        dynamics, k,
        "FrequencyRange" -> OptionValue["FrequencyRange"],
        "ImaginaryTolerance" -> OptionValue["ImaginaryTolerance"],
        "NormalizeModes" -> True
    ];
    If[modes === $Failed, Return[$Failed]];
    ExportMagneticModes[
        dynamics, modes, file,
        "NormalizationTolerance" -> OptionValue["NormalizationTolerance"],
        "CreateDirectories" -> OptionValue["CreateDirectories"]
    ]
];

ExportMagneticModesHelp[] := Column[{
    Style["ExportMagneticModes.wl", Bold, 16],
    "Export full Cartesian Fourier eigenvectors and equilibrium Fourier magnetization to JSON.",
    "All arrays use dynamics[\"Nodes\"] ordering. Modes must be BdG normalized.",
    "ExportMagneticModes[dynamics, modes, filename]",
    "ComputeAndExportMagneticModes[dynamics, k, filename, \"FrequencyRange\" -> {0, 3}]"
}];

Print["Loaded ExportMagneticModes.wl (isotropic magnon export)."];