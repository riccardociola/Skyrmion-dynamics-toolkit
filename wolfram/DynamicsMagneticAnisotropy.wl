(* ================================================================ *)
(* DynamicsMagneticAnisotropy.wl                                    *)
(*                                                                  *)
(* Linearized LL dynamics including cubic magnetocrystalline MCA.   *)
(* Works with either SixFold or General equilibrium Fourier data,   *)
(* and with rigid or deformed/elongated reciprocal lattices.        *)
(*                                                                  *)
(* Prerequisites, loaded in this order:                              *)
(* Geometry.wl, FourierBasis.wl, EnergyAnisotropy.wl,                *)
(* GroundState.wl, DynamicsProjection.wl, DynamicsMagnetic.wl        *)
(*                                                                  *)
(* Convention: the magnetic Hessian H here is HALF the second       *)
(* functional derivative of the free energy. In particular,        *)
(*                                                                  *)
(*   Fcub = kappaCub Sum_i <(vi . m)^4>                             *)
(*   Hcub_ab(r) = 6 kappaCub Sum_i (vi . m0(r))^2 vi_a vi_b.        *)
(*                                                                  *)
(* The projected LL matrix is D(k) = 2 I Cperp Hperp(k), where     *)
(* Cperp is the projection of (M0 x).                               *)
(*                                                                  *)
(* No existing isotropic solver functions are redefined.           *)
(* ================================================================ *)

ClearAll[
    BuildCubicAnisotropyHessian,
    BuildAnisotropicMagneticHessian,
    BuildProjectedAnisotropicMagneticHessian,
    PrepareMagneticAnisotropyDynamics,
    BuildAnisotropicLLDynamicalMatrix,
    SolveAnisotropicMagneticModes,
    DynamicsMagneticAnisotropyHelp
];

BuildCubicAnisotropyHessian::mca =
    "The supplied equilibrium does not have CubicAnisotropy -> True.";
BuildCubicAnisotropyHessian::kappa =
    "A numeric kappaCub parameter is missing from the saved equilibrium.";
BuildCubicAnisotropyHessian::field =
    "The equilibrium has no valid crystallographic field direction.";
BuildCubicAnisotropyHessian::dependency =
    "Load EnergyAnisotropy.wl and DynamicsProjection.wl before constructing the MCA Hessian.";
BuildCubicAnisotropyHessian::fourier =
    "Unresolved or nonnumeric equilibrium Fourier coefficients were encountered.";

PrepareMagneticAnisotropyDynamics::isotropic =
    "This solver requires an MCA equilibrium (CubicAnisotropy -> True). Use PrepareMagneticDynamics for isotropic equilibria.";
PrepareMagneticAnisotropyDynamics::dmi =
    "The current analytic quadratic Hessian supports Bloch DMI only (received `1`).";
PrepareMagneticAnisotropyDynamics::dependency =
    "Load DynamicsProjection.wl and DynamicsMagnetic.wl before preparing MCA dynamics.";

(* --------------------------------------------------------------- *)
(* Fourier-space cubic-anisotropy Hessian                           *)
(* --------------------------------------------------------------- *)

(*
    Let Gi, Gj be the retained magnetic fluctuation nodes. Define

        A_nu(G) = v_nu . M0_cart(G),
        S_nu(Q) = Sum_P A_nu(P) A_nu(Q-P).

    Then the full 3 x 3 block is

        Hcub[Gi,Gj] = 6 kappaCub Sum_nu
                        S_nu(Gi-Gj) Outer[Times,v_nu,v_nu].

    This depends only on Gi-Gj and NOT on the Bloch wavevector k.
    The projection to transverse modes is performed once per equilibrium.

    The Fourier cutoff refers to the equilibrium magnetic texture,
    irrespective of whether its Fourier parametrization is SixFold or
    General. It does not restrict the Bravais lattice to be hexagonal:
    eta0, phi0, and alphaR0 are always taken from the saved solution.
*)
BuildCubicAnisotropyHessian[
    equilibrium_Association,
    nodes_List
] :=
    Module[
        {
            rules, kappaValue, fieldDirection, cubicAxes,
            etaValue, phiValue, alphaValue,
            nNodes, magnetizationCartesian, projections, dyads,
            pairConvolution, nodeKey, key, value, zeroBlock,
            blocks, i, j
        },

        If[!TrueQ[Lookup[equilibrium, "CubicAnisotropy", False]],
            Message[BuildCubicAnisotropyHessian::mca];
            Return[$Failed]
        ];

        If[Length[DownValues[CubicCrystalAxes]] == 0 ||
           Length[DownValues[EquilibriumMagnetization]] == 0,
            Message[BuildCubicAnisotropyHessian::dependency];
            Return[$Failed]
        ];

        rules = EquilibriumRules[equilibrium];
        kappaValue = kappaCub /. rules;
        If[!NumericQ[kappaValue],
            Message[BuildCubicAnisotropyHessian::kappa];
            Return[$Failed]
        ];

        (* Use the saved input field, not a newly chosen direction. *)
        fieldDirection = Lookup[
            equilibrium,
            "FieldDirectionInput",
            Lookup[equilibrium, "FieldDirectionCrystal", Missing["Absent"]]
        ];
        If[!ListQ[fieldDirection] || Length[fieldDirection] != 3,
            Message[BuildCubicAnisotropyHessian::field];
            Return[$Failed]
        ];

        cubicAxes = CubicCrystalAxes[fieldDirection];
        If[cubicAxes === $Failed,
            Message[BuildCubicAnisotropyHessian::field];
            Return[$Failed]
        ];

        etaValue = eta0 /. rules;
        phiValue = phi0 /. rules;
        alphaValue = alphaR0 /. rules;
        nNodes = Length[nodes];

        (* Cached Cartesian Fourier amplitudes. *)
        magnetizationCartesian = Table[
            toCRT[
                EquilibriumMagnetization[
                    equilibrium,
                    Sequence @@ nodes[[i]]
                ],
                etaValue, phiValue, alphaValue
            ],
            {i, nNodes}
        ];

        If[!MatrixQ[N[magnetizationCartesian], NumericQ],
            Message[BuildCubicAnisotropyHessian::fourier];
            Return[$Failed]
        ];

        projections = Table[
            cubicAxes[[a]] . magnetizationCartesian[[i]],
            {i, nNodes}, {a, 1, 3}
        ];
        dyads = (Outer[Times, #, #] &) /@ cubicAxes;
        zeroBlock = ConstantArray[0, {3, 3}];

        (* Each integer pair {n1,n2} has a unique Gaussian-integer key. *)
        nodeKey[node_List] := node[[1]] + I node[[2]];

        (* Accumulate Hcub(Q) once, using all ordered pairs of nodes. *)
        pairConvolution = <||>;
        Do[
            key = nodeKey[nodes[[i]] + nodes[[j]]];
            value = 6 kappaValue Total[
                Table[
                    projections[[i, a]] projections[[j, a]] dyads[[a]],
                    {a, 1, 3}
                ]
            ];
            AssociateTo[
                pairConvolution,
                key -> (Lookup[pairConvolution, key, zeroBlock] + value)
            ],
            {i, nNodes}, {j, nNodes}
        ];

        blocks = Table[
            Lookup[
                pairConvolution,
                nodeKey[nodes[[i]] - nodes[[j]]],
                zeroBlock
            ],
            {i, nNodes}, {j, nNodes}
        ];

        Chop[ArrayFlatten[blocks]]
    ];

BuildCubicAnisotropyHessian[equilibrium_Association] :=
    BuildCubicAnisotropyHessian[
        equilibrium,
        DynamicsNodes[equilibrium["FourierCutoff"]]
    ];

(* --------------------------------------------------------------- *)
(* Full and projected magnetic Hessians                            *)
(* --------------------------------------------------------------- *)

BuildAnisotropicMagneticHessian[
    dynamicsData_Association,
    k_List
] :=
    Module[{quadratic},
        quadratic = BuildMagneticQuadraticHessian[
            dynamicsData["Equilibrium"], k, dynamicsData["Nodes"]
        ];
        If[quadratic === $Failed, Return[$Failed]];
        Chop[
            Normal[quadratic] +
            dynamicsData["ConstraintHessian"] +
            dynamicsData["AnisotropyHessian"]
        ]
    ];

BuildProjectedAnisotropicMagneticHessian[
    dynamicsData_Association,
    k_List
] :=
    Module[{quadratic, projectedQuadratic},
        quadratic = BuildMagneticQuadraticHessian[
            dynamicsData["Equilibrium"], k, dynamicsData["Nodes"]
        ];
        If[quadratic === $Failed, Return[$Failed]];
        projectedQuadratic = ProjectMagneticOperator[
            dynamicsData["TransverseData"], quadratic
        ];
        Chop[
            projectedQuadratic + dynamicsData["TextureHessianTransverse"]
        ]
    ];

(* --------------------------------------------------------------- *)
(* Prepare and project k-independent contributions once            *)
(* --------------------------------------------------------------- *)

Options[PrepareMagneticAnisotropyDynamics] = {
    "ProjectionCutoff" -> 0.97,
    "CheckInverse" -> False,
    "InverseTolerance" -> 10^-10
};

PrepareMagneticAnisotropyDynamics[
    equilibrium_Association,
    OptionsPattern[]
] :=
    Module[
        {
            nodes, crossOperator, transverseData,
            constraintHessian, anisotropyHessian,
            constraintTransverse, anisotropyTransverse
        },

        If[!TrueQ[Lookup[equilibrium, "CubicAnisotropy", False]],
            Message[PrepareMagneticAnisotropyDynamics::isotropic];
            Return[$Failed]
        ];
        If[Lookup[equilibrium, "DMIType", Missing[]] =!= "Bloch",
            Message[PrepareMagneticAnisotropyDynamics::dmi,
                Lookup[equilibrium, "DMIType", Missing[]]];
            Return[$Failed]
        ];
        If[Length[DownValues[BuildConstraintHessian]] == 0 ||
           Length[DownValues[BuildMagneticTransverseBasis]] == 0,
            Message[PrepareMagneticAnisotropyDynamics::dependency];
            Return[$Failed]
        ];

        nodes = DynamicsNodes[equilibrium["FourierCutoff"]];
        crossOperator = BuildM0CrossOperator[equilibrium, nodes];
        transverseData = BuildMagneticTransverseBasis[
            crossOperator,
            "ProjectionCutoff" -> OptionValue["ProjectionCutoff"],
            "CheckInverse" -> OptionValue["CheckInverse"],
            "InverseTolerance" -> OptionValue["InverseTolerance"]
        ];
        If[transverseData === $Failed, Return[$Failed]];

        constraintHessian = BuildConstraintHessian[equilibrium, nodes];
        anisotropyHessian = BuildCubicAnisotropyHessian[equilibrium, nodes];
        If[anisotropyHessian === $Failed, Return[$Failed]];

        (* The full fluctuation Hilbert space, not the restricted        *)
        (* equilibrium representation, is projected onto M0 x.        *)
        constraintTransverse = ProjectMagneticOperator[
            transverseData, constraintHessian
        ];
        anisotropyTransverse = ProjectMagneticOperator[
            transverseData, anisotropyHessian
        ];

        <|
            "Model" -> "CubicAnisotropy",
            "Equilibrium" -> equilibrium,
            "Nodes" -> nodes,
            "CrossOperatorFull" -> crossOperator,
            "TransverseData" -> transverseData,
            "ConstraintHessian" -> constraintHessian,
            "ConstraintHessianTransverse" -> constraintTransverse,
            "AnisotropyHessian" -> anisotropyHessian,
            "AnisotropyHessianTransverse" -> anisotropyTransverse,
            "TextureHessianTransverse" -> Chop[
                constraintTransverse + anisotropyTransverse
            ],
            "FullDimension" -> 3 Length[nodes],
            "TransverseDimension" -> transverseData["TransverseDimension"]
        |>
    ];

(* --------------------------------------------------------------- *)
(* LL matrix and eigenmodes                                        *)
(* --------------------------------------------------------------- *)

BuildAnisotropicLLDynamicalMatrix[
    dynamicsData_Association,
    k_List
] :=
    Module[{hessian, cross},
        hessian = BuildProjectedAnisotropicMagneticHessian[dynamicsData, k];
        If[hessian === $Failed, Return[$Failed]];
        cross = dynamicsData["TransverseData"]["CrossOperator"];
        <|
            "HessianTransverse" -> hessian,
            "DynamicalMatrix" -> Chop[2 I cross . hessian]
        |>
    ];

Options[SolveAnisotropicMagneticModes] = {
    "FrequencyRange" -> {0, Infinity},
    "ImaginaryTolerance" -> 10^-8,
    "NormalizeModes" -> True
};

(* Same selection, sorting, and BdG normalization as the             *)
(* isotropic SolveMagneticModes; only the Hessian differs.           *)
SolveAnisotropicMagneticModes[
    dynamicsData_Association,
    k_List : {0, 0, 0},
    OptionsPattern[]
] :=
    Module[
        {
            matrixData, dynamicalMatrix, eigenvalues, eigenvectors,
            frequencyRange, imaginaryTolerance, normalizeQ,
            modePairs, selectedPairs, frequencies, vectorsTransverse,
            norms, normalizedData, vectorsFull
        },

        matrixData = BuildAnisotropicLLDynamicalMatrix[dynamicsData, k];
        If[matrixData === $Failed, Return[$Failed]];
        dynamicalMatrix = matrixData["DynamicalMatrix"];
        {eigenvalues, eigenvectors} = Eigensystem[dynamicalMatrix];

        frequencyRange = OptionValue["FrequencyRange"];
        imaginaryTolerance = OptionValue["ImaginaryTolerance"];
        normalizeQ = OptionValue["NormalizeModes"];

        modePairs = Transpose[{eigenvalues, eigenvectors}];
        selectedPairs = Select[
            modePairs,
            frequencyRange[[1]] < Re[#[[1]]] <= frequencyRange[[2]] &&
            Abs[Im[#[[1]]]] <= imaginaryTolerance &
        ];
        selectedPairs = SortBy[selectedPairs, Re[#[[1]]] &];

        If[selectedPairs === {},
            Return[<|
                "WaveVector" -> k,
                "Frequencies" -> {},
                "EigenvectorsTransverse" -> {},
                "EigenvectorsFull" -> {},
                "BdGNormsBeforeNormalization" -> {},
                "DynamicalMatrix" -> dynamicalMatrix,
                "HessianTransverse" -> matrixData["HessianTransverse"],
                "NumberOfModes" -> 0
            |>]
        ];

        frequencies = Chop[Re[selectedPairs[[All, 1]]]];
        vectorsTransverse = selectedPairs[[All, 2]];

        If[TrueQ[normalizeQ],
            normalizedData = NormalizeMagneticMode[
                dynamicsData["TransverseData"], #
            ] & /@ vectorsTransverse;
            norms = normalizedData[[All, "OriginalNorm"]];
            vectorsTransverse = normalizedData[[All, "Vector"]],
            norms = MagneticBdGNorm[
                dynamicsData["TransverseData"], #
            ] & /@ vectorsTransverse
        ];

        vectorsFull = LiftMagneticVector[
            dynamicsData["TransverseData"], #
        ] & /@ vectorsTransverse;

        <|
            "WaveVector" -> k,
            "Frequencies" -> frequencies,
            "EigenvectorsTransverse" -> vectorsTransverse,
            "EigenvectorsFull" -> vectorsFull,
            "BdGNormsBeforeNormalization" -> norms,
            "DynamicalMatrix" -> dynamicalMatrix,
            "HessianTransverse" -> matrixData["HessianTransverse"],
            "NumberOfModes" -> Length[frequencies]
        |>
    ];

(* --------------------------------------------------------------- *)
(* Help                                                             *)
(* --------------------------------------------------------------- *)

DynamicsMagneticAnisotropyHelp[] :=
    Column[{
        Style["DynamicsMagneticAnisotropy.wl", Bold, 16],
        "Cubic-MCA Hessian and projected magnetic Landau-Lifshitz modes.",
        Grid[{
            {"Function", "Description"},
            {"BuildCubicAnisotropyHessian", "Construct the MCA convolution Hessian (independent of Bloch wavevector)."},
            {"PrepareMagneticAnisotropyDynamics", "Build and project cross operator, constraint, and cubic MCA Hessians once."},
            {"BuildAnisotropicMagneticHessian", "Assemble the full Cartesian Fourier-space Hessian."},
            {"BuildProjectedAnisotropicMagneticHessian", "Assemble the reduced Hessian for one Cartesian Bloch wavevector."},
            {"BuildAnisotropicLLDynamicalMatrix", "Return the projected LL matrix and Hessian."},
            {"SolveAnisotropicMagneticModes", "Find positive-frequency magnetic modes and normalize with the BdG metric."}
        }, Frame -> All, Alignment -> Left]
    }];

Print["Loaded DynamicsMagneticAnisotropy.wl (run DynamicsMagneticAnisotropyHelp[] for functions)."];