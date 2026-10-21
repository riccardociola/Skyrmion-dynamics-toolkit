(* ================================================================ *)
(* DynamicsMagnetic.wl                                              *)
(*                                                                  *)
(* Linearized Landau-Lifshitz dynamics for magnetic skyrmion       *)
(* lattices around a saved equilibrium configuration.               *)
(*                                                                  *)
(* The magnetic Hessian is constructed analytically and then        *)
(* projected onto the transverse subspace supplied by               *)
(* DynamicsProjection.wl.                                           *)
(* ================================================================ *)

(* Dependencies:                                                    *)
(*   Geometry.wl                                                    *)
(*   FourierBasis.wl                                                *)
(*   GroundState.wl                                                 *)
(*   DynamicsProjection.wl                                         *)

ClearAll[
    BlochDMIHessianBlock,
    MagneticQuadraticBlock,
    BuildMagneticQuadraticHessian,
    BuildConstraintHessian,
    BuildMagneticHessian,
    BuildProjectedMagneticHessian,
    PrepareMagneticDynamics,
    BuildLLDynamicalMatrix,
    MagneticBdGInnerProduct,
    MagneticBdGNorm,
    NormalizeMagneticMode,
    SolveMagneticModes,
    DynamicsMagneticHelp
];

(* --------------------------------------------------------------- *)
(* Bloch-DMI Hessian block                                         *)
(* --------------------------------------------------------------- *)

(*
    Return the 3 x 3 Cartesian Hessian block of the Bloch DMI
    for a Cartesian wavevector q.

    This is the analytical term appearing in the legacy dynamics
    notebook:

        -i epsilon_ijk q_k.
*)
BlochDMIHessianBlock[
    q_List
] :=
    -I TensorContract[
        TensorProduct[
            LeviCivitaTensor[3],
            q
        ],
        {{3, 4}}
    ];

(* --------------------------------------------------------------- *)
(* Quadratic magnetic Hessian block                                *)
(* --------------------------------------------------------------- *)

MagneticQuadraticBlock::dmi =
    "Dynamics for DMI type `1` is not implemented yet. Currently only \"Bloch\" is supported.";

(*
    Return the 3 x 3 Cartesian Hessian block associated with the
    quadratic magnetic terms at Cartesian wavevector q.

    Included terms:
        - quadratic part of the magnetization-length constraint;
        - exchange interaction;
        - DMI;
        - magnetostatic interaction.

    The Zeeman term does not contribute because its second
    functional derivative is zero.
*)
MagneticQuadraticBlock[
    q_List,
    lambdaValue_,
    xiValue_,
    nijCartesian_List,
    dmiType_String
] :=
    Module[
        {
            q2,
            dmiBlock,
            magnetostaticBlock
        },

        q2 =
            q . q;

        dmiBlock =
            Switch[
                dmiType,

                "Bloch",
                BlochDMIHessianBlock[q],

                _,
                Message[
                    MagneticQuadraticBlock::dmi,
                    dmiType
                ];
                Return[$Failed]
            ];

        (* At q = 0 use the sample demagnetization tensor. *)
        magnetostaticBlock =
            If[
                TrueQ[PossibleZeroQ[Chop[q2]]],

                xiValue nijCartesian,

                xiValue
                TensorProduct[q, q] / q2
            ];

        (
            -2 lambdaValue + q2/2
        ) IdentityMatrix[3]
        +
        dmiBlock
        +
        magnetostaticBlock
    ];

(* --------------------------------------------------------------- *)
(* Reciprocal-space diagonal Hessian                              *)
(* --------------------------------------------------------------- *)

(*
    Construct the reciprocal-node-diagonal part H_k of the magnetic
    Hessian.

    For each reciprocal node G, the 3 x 3 block is evaluated at

        q = k + G_cartesian/a0,

    where k is a Cartesian Bloch wavevector in the first Brillouin
    zone.

    The result is stored as a SparseArray because this contribution
    is block diagonal in reciprocal-node space.
*)
BuildMagneticQuadraticHessian[
    equilibrium_Association,
    k_List,
    nodes_List
] :=
    Module[
        {
            rules,
            aValue,
            etaValue,
            phiValue,
            alphaValue,
            lambdaValue,
            xiValue,
            nijCartesian,
            dmiType,
            numberOfNodes,
            blocks,
            rulesSparse
        },

        rules =
            EquilibriumRules[equilibrium];

        aValue =
            a0 /. rules;

        etaValue =
            eta0 /. rules;

        phiValue =
            phi0 /. rules;

        alphaValue =
            alphaR0 /. rules;

        lambdaValue =
            lambda /. rules;

        xiValue =
            xi /. rules;

        nijCartesian =
            DiagonalMatrix[
                {
                    Nx /. rules,
                    Ny /. rules,
                    Nz /. rules
                }
            ];

        dmiType =
            equilibrium["DMIType"];

        numberOfNodes =
            Length[nodes];

        blocks =
            Table[
                Module[
                    {q},

                    q =
                        k
                        +
                        1/aValue
                        toCRT[
                            G[
                                Sequence @@ nodes[[i]]
                            ],
                            etaValue,
                            phiValue,
                            alphaValue
                        ];

                    MagneticQuadraticBlock[
                        q,
                        lambdaValue,
                        xiValue,
                        nijCartesian,
                        dmiType
                    ]
                ],

                {i, numberOfNodes}
            ];

        If[
            MemberQ[blocks, $Failed, Infinity],
            Return[$Failed]
        ];

        (* Insert each 3 x 3 block directly into a sparse matrix. *)
        rulesSparse =
            Flatten[
                Table[
                    Table[
                        {
                            3 (i - 1) + a,
                            3 (i - 1) + b
                        }
                        ->
                        blocks[[i, a, b]],

                        {a, 1, 3},
                        {b, 1, 3}
                    ],
                    {i, numberOfNodes}
                ],
                2
            ];

        SparseArray[
            rulesSparse,
            {
                3 numberOfNodes,
                3 numberOfNodes
            }
        ]
    ];

(* Use the equilibrium cutoff when nodes are omitted. *)
BuildMagneticQuadraticHessian[
    equilibrium_Association,
    k_List
] :=
    BuildMagneticQuadraticHessian[
        equilibrium,
        k,
        DynamicsNodes[
            equilibrium["FourierCutoff"]
        ]
    ];

(* --------------------------------------------------------------- *)
(* Constraint / texture-dependent Hessian                         *)
(* --------------------------------------------------------------- *)

(*
    Construct the equilibrium-dependent Hessian contribution H_e
    generated by the quartic magnetization-length constraint.

    In Fourier space this term couples different reciprocal nodes
    through convolutions of the equilibrium magnetization M0.

    This contribution is independent of the Bloch wavevector k and
    therefore only needs to be calculated once for a given
    equilibrium configuration.
*)
BuildConstraintHessian[
    equilibrium_Association,
    nodes_List
] :=
    Module[
        {
            rules,
            cutoff,
            etaValue,
            phiValue,
            alphaValue,
            lambdaValue,
            numberOfNodes,
            mAffineTable,
            mCartesianTable,
            mAffine,
            mCartesian
        },

        rules =
            EquilibriumRules[equilibrium];

        cutoff =
            equilibrium["FourierCutoff"];

        etaValue =
            eta0 /. rules;

        phiValue =
            phi0 /. rules;

        alphaValue =
            alphaR0 /. rules;

        lambdaValue =
            lambda /. rules;

        numberOfNodes =
            Length[nodes];

        (*
            Cache all equilibrium Fourier components once. H_e uses
            the same M0 coefficients many times, so repeated rule
            substitution inside the triple convolution is avoided.
        *)
        mAffineTable =
            Table[
                EquilibriumMagnetization[
                    equilibrium,
                    n1,
                    n2
                ],
                {n1, -cutoff, cutoff},
                {n2, -cutoff, cutoff}
            ];

        mCartesianTable =
            Map[
                toCRT[
                    #,
                    etaValue,
                    phiValue,
                    alphaValue
                ] &,
                mAffineTable,
                {2}
            ];

        (* Return zero automatically outside the retained hexagon. *)
        mAffine[q_List] :=
            If[
                hexagonNodeQ[
                    q[[1]],
                    q[[2]],
                    cutoff
                ],
                mAffineTable[[
                    q[[1]] + cutoff + 1,
                    q[[2]] + cutoff + 1
                ]],
                {0, 0, 0}
            ];

        mCartesian[q_List] :=
            If[
                hexagonNodeQ[
                    q[[1]],
                    q[[2]],
                    cutoff
                ],
                mCartesianTable[[
                    q[[1]] + cutoff + 1,
                    q[[2]] + cutoff + 1
                ]],
                {0, 0, 0}
            ];

        ArrayFlatten[
            Table[
                Module[
                    {
                        q1,
                        q2,
                        scalarPart,
                        tensorPart
                    },

                    q1 =
                        nodes[[i]];

                    q2 =
                        nodes[[j]];

                    scalarPart =
                        Sum[
                            DotC[
                                mAffine[-q3],
                                mAffine[q1 - q2 + q3],
                                etaValue,
                                phiValue,
                                alphaValue
                            ],
                            {q3, nodes}
                        ];

                    tensorPart =
                        Sum[
                            TensorProduct[
                                mCartesian[-q3],
                                mCartesian[q1 - q2 + q3]
                            ],
                            {q3, nodes}
                        ];

                    2 lambdaValue
                    scalarPart
                    IdentityMatrix[3]
                    +
                    4 lambdaValue
                    tensorPart
                ],

                {i, numberOfNodes},
                {j, numberOfNodes}
            ]
        ]
    ];

(* Use the equilibrium cutoff when nodes are omitted. *)
BuildConstraintHessian[
    equilibrium_Association
] :=
    BuildConstraintHessian[
        equilibrium,
        DynamicsNodes[
            equilibrium["FourierCutoff"]
        ]
    ];

(* --------------------------------------------------------------- *)
(* Full magnetic Hessian                                          *)
(* --------------------------------------------------------------- *)

(*
    Construct the full magnetic Hessian for a prepared dynamics
    object.

    Additional magnetic contributions such as magnetocrystalline
    anisotropy can later be added here without modifying the LL
    equation-of-motion code.
*)
BuildMagneticHessian[
    dynamicsData_Association,
    k_List
] :=
    Module[
        {
            quadraticHessian
        },

        quadraticHessian =
            BuildMagneticQuadraticHessian[
                dynamicsData["Equilibrium"],
                k,
                dynamicsData["Nodes"]
            ];

        If[
            quadraticHessian === $Failed,
            Return[$Failed]
        ];

        Chop[
            Normal[quadraticHessian]
            +
            dynamicsData["ConstraintHessian"]
        ]
    ];

(* --------------------------------------------------------------- *)
(* Projected magnetic Hessian                                     *)
(* --------------------------------------------------------------- *)

(*
    Construct the magnetic Hessian directly in the transverse basis.

    The k-independent constraint contribution is projected once by
    PrepareMagneticDynamics. For every new k point only the sparse
    reciprocal-diagonal contribution H_k is rebuilt and projected.
*)
BuildProjectedMagneticHessian[
    dynamicsData_Association,
    k_List
] :=
    Module[
        {
            quadraticHessian,
            quadraticTransverse
        },

        quadraticHessian =
            BuildMagneticQuadraticHessian[
                dynamicsData["Equilibrium"],
                k,
                dynamicsData["Nodes"]
            ];

        If[
            quadraticHessian === $Failed,
            Return[$Failed]
        ];

        quadraticTransverse =
            ProjectMagneticOperator[
                dynamicsData["TransverseData"],
                quadraticHessian
            ];

        Chop[
            quadraticTransverse
            +
            dynamicsData["ConstraintHessianTransverse"]
        ]
    ];


(* --------------------------------------------------------------- *)
(* Prepare reusable dynamics data                                 *)
(* --------------------------------------------------------------- *)

Options[
    PrepareMagneticDynamics
] = {
    "ProjectionCutoff" -> 0.97,
    "CheckInverse" -> False,
    "InverseTolerance" -> 10^-10
};

(*
    This solver deliberately supports isotropic magnetic energies only.
    Reject MCA equilibria here, before building a Hessian that would
    omit the cubic-anisotropy contribution.
    Older saved isotropic equilibria without this metadata are accepted.
*)
PrepareMagneticDynamics::mca =
    "This solver supports only isotropic magnetic equilibria. The supplied equilibrium has CubicAnisotropy -> True; use the anisotropic dynamics solver instead.";

(*
    Precompute all quantities that are independent of the Bloch
    wavevector k:

        - dynamical reciprocal-space nodes;
        - Fourier-space M0 x operator;
        - transverse magnetic basis;
        - inverse of M0 x in the transverse basis;
        - texture-dependent constraint Hessian H_e.

    The returned Association can be reused for many k points.
*)
PrepareMagneticDynamics[
    equilibrium_Association,
    OptionsPattern[]
] :=
    Module[
        {
            nodes,
            crossOperator,
            transverseData,
            constraintHessian,
            constraintHessianTransverse
        },

        If[
            TrueQ[Lookup[equilibrium, "CubicAnisotropy", False]],
            Message[PrepareMagneticDynamics::mca];
            Return[$Failed]
        ];

        nodes =
            DynamicsNodes[
                equilibrium["FourierCutoff"]
            ];

        crossOperator =
            BuildM0CrossOperator[
                equilibrium,
                nodes
            ];

        transverseData =
            BuildMagneticTransverseBasis[
                crossOperator,
                "ProjectionCutoff" ->
                    OptionValue["ProjectionCutoff"],
                "CheckInverse" ->
                    OptionValue["CheckInverse"],
                "InverseTolerance" ->
                    OptionValue["InverseTolerance"]
            ];

        If[
            transverseData === $Failed,
            Return[$Failed]
        ];

        constraintHessian =
            BuildConstraintHessian[
                equilibrium,
                nodes
            ];

        (* Project the k-independent dense contribution only once. *)
        constraintHessianTransverse =
            ProjectMagneticOperator[
                transverseData,
                constraintHessian
            ];

        <|
            "Equilibrium" ->
                equilibrium,

            "Nodes" ->
                nodes,

            "CrossOperatorFull" ->
                crossOperator,

            "TransverseData" ->
                transverseData,

            "ConstraintHessian" ->
                constraintHessian,

            "ConstraintHessianTransverse" ->
                constraintHessianTransverse,

            "FullDimension" ->
                3 Length[nodes],

            "TransverseDimension" ->
                transverseData["TransverseDimension"]
        |>
    ];

(* --------------------------------------------------------------- *)
(* Linearized LL dynamical matrix                                 *)
(* --------------------------------------------------------------- *)

(*
    Construct the linearized Landau-Lifshitz matrix in the
    transverse magnetic basis:

        D(k) = 2 i C_perp H_perp(k),

    where C_perp is the projected M0 x operator and H_perp is the
    magnetic Hessian projected onto the same transverse subspace.
*)
BuildLLDynamicalMatrix[
    dynamicsData_Association,
    k_List
] :=
    Module[
        {
            hessianTransverse,
            crossTransverse
        },

        hessianTransverse =
            BuildProjectedMagneticHessian[
                dynamicsData,
                k
            ];

        If[
            hessianTransverse === $Failed,
            Return[$Failed]
        ];

        crossTransverse =
            dynamicsData[
                "TransverseData"
            ]["CrossOperator"];

        <|
            "HessianTransverse" ->
                hessianTransverse,

            "DynamicalMatrix" ->
                Chop[
                    2 I
                    crossTransverse
                    .
                    hessianTransverse
                ]
        |>
    ];

(* --------------------------------------------------------------- *)
(* Magnetic BdG inner product                                     *)
(* --------------------------------------------------------------- *)

(*
    BdG inner product for two vectors represented in the transverse
    magnetic basis:

        <u,v>_BdG = u^dagger [-i C_perp^-1] v.
*)
MagneticBdGInnerProduct[
    transverseData_Association,
    u_List,
    v_List
] :=
    Chop[
        Conjugate[u]
        .
        (
            -I
            transverseData["InverseCrossOperator"]
        )
        .
        v
    ];

(* BdG norm of one transverse magnetic mode. *)
MagneticBdGNorm[
    transverseData_Association,
    vector_List
] :=
    MagneticBdGInnerProduct[
        transverseData,
        vector,
        vector
    ];

(*
    Normalize one transverse magnetic eigenvector using the BdG norm.
    The function returns both the normalized vector and its original
    norm for diagnostics.
*)
NormalizeMagneticMode[
    transverseData_Association,
    vector_List
] :=
    Module[
        {norm},

        norm =
            MagneticBdGNorm[
                transverseData,
                vector
            ];

        <|
            "Vector" ->
                Chop[
                    vector / Sqrt[norm]
                ],

            "OriginalNorm" ->
                norm
        |>
    ];

(* --------------------------------------------------------------- *)
(* Solve the magnetic eigensystem                                 *)
(* --------------------------------------------------------------- *)

Options[
    SolveMagneticModes
] = {
    "FrequencyRange" -> {0, Infinity},
    "ImaginaryTolerance" -> 10^-8,
    "NormalizeModes" -> True
};

(*
    Diagonalize the projected LL dynamical matrix and return the
    positive-frequency modes in the requested frequency interval.

    This first implementation intentionally uses the full dense
    Eigensystem so that it can be benchmarked directly against the
    legacy notebook. Partial/iterative diagonalization can be added
    later without changing the Hessian construction.
*)
SolveMagneticModes[
    dynamicsData_Association,
    k_List : {0, 0, 0},
    OptionsPattern[]
] :=
    Module[
        {
            matrixData,
            dynamicalMatrix,
            eigenvalues,
            eigenvectors,
            frequencyRange,
            imaginaryTolerance,
            normalizeQ,
            modePairs,
            selectedPairs,
            frequencies,
            vectorsTransverse,
            norms,
            normalizedData,
            vectorsFull
        },

        matrixData =
            BuildLLDynamicalMatrix[
                dynamicsData,
                k
            ];

        If[
            matrixData === $Failed,
            Return[$Failed]
        ];

        dynamicalMatrix =
            matrixData["DynamicalMatrix"];

        {
            eigenvalues,
            eigenvectors
        } =
            Eigensystem[
                dynamicalMatrix
            ];

        frequencyRange =
            OptionValue["FrequencyRange"];

        imaginaryTolerance =
            OptionValue["ImaginaryTolerance"];

        normalizeQ =
            OptionValue["NormalizeModes"];

        modePairs =
            Transpose[
                {
                    eigenvalues,
                    eigenvectors
                }
            ];

        (* Keep stable modes with frequencies inside the requested range. *)
        selectedPairs =
            Select[
                modePairs,
                (
                    frequencyRange[[1]] < Re[#[[1]]] <= frequencyRange[[2]]
                    &&
                    Abs[Im[#[[1]]]] <= imaginaryTolerance
                ) &
            ];

        selectedPairs =
            SortBy[
                selectedPairs,
                Re[#[[1]]] &
            ];

        If[
            selectedPairs === {},

            Return[
                <|
                    "WaveVector" -> k,
                    "Frequencies" -> {},
                    "EigenvectorsTransverse" -> {},
                    "EigenvectorsFull" -> {},
                    "BdGNormsBeforeNormalization" -> {},
                    "DynamicalMatrix" -> dynamicalMatrix,
                    "HessianTransverse" -> matrixData["HessianTransverse"],
                    "NumberOfModes" -> 0
                |>
            ]
        ];

        frequencies =
            Chop[
                Re[
                    selectedPairs[[All, 1]]
                ]
            ];

        vectorsTransverse =
            selectedPairs[[All, 2]];

        If[
            TrueQ[normalizeQ],

            normalizedData =
                NormalizeMagneticMode[
                    dynamicsData["TransverseData"],
                    #
                ] & /@ vectorsTransverse;

            norms =
                normalizedData[[All, "OriginalNorm"]];

            vectorsTransverse =
                normalizedData[[All, "Vector"]],

            norms =
                MagneticBdGNorm[
                    dynamicsData["TransverseData"],
                    #
                ] & /@ vectorsTransverse
        ];

        vectorsFull =
            LiftMagneticVector[
                dynamicsData["TransverseData"],
                #
            ] & /@ vectorsTransverse;

        <|
            "WaveVector" ->
                k,

            "Frequencies" ->
                frequencies,

            "EigenvectorsTransverse" ->
                vectorsTransverse,

            "EigenvectorsFull" ->
                vectorsFull,

            "BdGNormsBeforeNormalization" ->
                norms,

            "DynamicalMatrix" ->
                dynamicalMatrix,

            "HessianTransverse" ->
                matrixData["HessianTransverse"],

            "NumberOfModes" ->
                Length[frequencies]
        |>
    ];

(* --------------------------------------------------------------- *)
(* Help                                                            *)
(* --------------------------------------------------------------- *)

DynamicsMagneticHelp[] :=
    Column[
        {
            Style[
                "DynamicsMagnetic.wl",
                Bold,
                16
            ],

            "Analytical magnetic Hessian construction and linearized Landau-Lifshitz dynamics around a saved equilibrium.",

            Grid[
                {
                    {"Function", "Description"},
                    {"BlochDMIHessianBlock", "Return the Cartesian Bloch-DMI Hessian block for one wavevector."},
                    {"MagneticQuadraticBlock", "Return the 3 x 3 quadratic magnetic Hessian block including exchange, DMI, dipolar and constraint terms."},
                    {"BuildMagneticQuadraticHessian", "Construct the reciprocal-node-diagonal contribution H_k."},
                    {"BuildConstraintHessian", "Construct the texture-dependent convolution contribution H_e from the quartic constraint."},
                    {"BuildMagneticHessian", "Assemble the full magnetic Hessian for a chosen Bloch wavevector."},
                    {"BuildProjectedMagneticHessian", "Construct the magnetic Hessian directly in the transverse basis, reusing the projected k-independent contribution."},
                    {"PrepareMagneticDynamics", "Precompute the k-independent projection and Hessian data for one equilibrium."},
                    {"BuildLLDynamicalMatrix", "Project the Hessian and construct the linearized LL dynamical matrix."},
                    {"MagneticBdGInnerProduct / MagneticBdGNorm", "Evaluate the magnetic BdG inner product and norm in the transverse basis."},
                    {"NormalizeMagneticMode", "Normalize one transverse magnetic eigenvector using the BdG norm."},
                    {"SolveMagneticModes", "Diagonalize the LL matrix and return positive-frequency normalized magnon modes."}
                },
                Frame -> All,
                Alignment -> Left
            ]
        }
    ];

Print[
    DynamicsMagneticHelp[]
];