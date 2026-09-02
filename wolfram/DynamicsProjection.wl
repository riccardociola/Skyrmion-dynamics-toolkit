(* ================================================================ *)
(* DynamicsProjection.wl                                            *)
(*                                                                  *)
(* Common Fourier-space projection tools for magnetic dynamics      *)
(* around a saved equilibrium configuration.                        *)
(*                                                                  *)
(* ================================================================ *)


ClearAll[
    DynamicsNodes,
    EquilibriumMagnetization,
    CrossProductMatrix,
    BuildM0CrossOperator,
    BuildMagneticTransverseBasis,
    ProjectMagneticVector,
    LiftMagneticVector,
    ProjectMagneticOperator,
    DynamicsProjectionHelp
];


(* --------------------------------------------------------------- *)
(* Dynamical reciprocal-space basis                               *)
(* --------------------------------------------------------------- *)

(*
    Return all reciprocal-lattice nodes inside the finite hexagonal cutoff.

    Unlike MagneticVariables, this is the full Fourier space used for 
    magnetic fluctuations. Each reciprocal node carries three magnetic components.
*)
DynamicsNodes[
    cutoff_Integer?NonNegative
] :=
    Select[
        Flatten[
            Table[
                {n1, n2},
                {n1, -cutoff, cutoff},
                {n2, -cutoff, cutoff}
            ],
            1
        ],

        hexagonNodeQ[
            #[[1]],
            #[[2]],
            cutoff
        ] &
    ];


(* --------------------------------------------------------------- *)
(* Equilibrium magnetization                                      *)
(* --------------------------------------------------------------- *)

(*
    Return the Fourier component M0[n1,n2] of a saved equilibrium.

    The cutoff and magnetic representation are read directly from
    the saved equilibrium metadata. MagneticFourierComponent is
    defined in FourierBasis.wl.
*)
EquilibriumMagnetization[
    equilibrium_Association,
    n1_Integer,
    n2_Integer
] :=
    Chop[
        MagneticFourierComponent[
            n1,
            n2,
            equilibrium["FourierCutoff"],
            equilibrium["MagneticRepresentation"]
        ]
        /. equilibrium["Solution"]
    ];


(* --------------------------------------------------------------- *)
(* Cartesian cross-product matrix                                 *)
(* --------------------------------------------------------------- *)

(*
    Return the 3 x 3 Cartesian matrix representing m x v.

        CrossProductMatrix[m].v == Cross[m,v]
*)
CrossProductMatrix[
    m_List
] :=
    {
        {0,      -m[[3]],  m[[2]]},
        {m[[3]],  0,      -m[[1]]},
        {-m[[2]], m[[1]],  0}
    };


(* --------------------------------------------------------------- *)
(* Fourier representation of M0 x                                *)
(* --------------------------------------------------------------- *)

(*
    Construct the full Fourier-space representation of M0(r) x.

    The block connecting reciprocal nodes G and G' contains M0[G-G'].
    The equilibrium magnetization is converted from affine covariant
    components to Cartesian components before the cross-product 
    matrix is constructed.
*)
BuildM0CrossOperator[
    equilibrium_Association,
    nodes_List
] :=
    Module[
        {
            rules,
            etaValue,
            phiValue,
            alphaValue,
            numberOfNodes
        },

        rules =
            EquilibriumRules[
                equilibrium
            ];

        etaValue =
            eta0 /. rules;

        phiValue =
            phi0 /. rules;

        alphaValue =
            alphaR0 /. rules;

        numberOfNodes =
            Length[
                nodes
            ];

        ArrayFlatten[
            Table[
                CrossProductMatrix[
                    toCRT[
                        EquilibriumMagnetization[
                            equilibrium,
                            Sequence @@
                                (
                                    nodes[[i]] - nodes[[j]]
                                )
                        ],
                        etaValue,
                        phiValue,
                        alphaValue
                    ]
                ],

                {i, numberOfNodes},
                {j, numberOfNodes}
            ]
        ]
    ];

(* Use the cutoff stored in the equilibrium when nodes are omitted. *)
BuildM0CrossOperator[
    equilibrium_Association
] :=
    BuildM0CrossOperator[
        equilibrium,
        DynamicsNodes[
            equilibrium["FourierCutoff"]
        ]
    ];


(* --------------------------------------------------------------- *)
(* Magnetic transverse subspace                                   *)
(* --------------------------------------------------------------- *)

Options[
    BuildMagneticTransverseBasis
] = {
    "ProjectionCutoff" -> 0.97,
    "CheckInverse" -> False,
    "InverseTolerance" -> 10^-10
};

BuildMagneticTransverseBasis::empty =
    "No transverse states satisfy the projection cutoff `1`.";

BuildMagneticTransverseBasis::inverse =
    "The spectral inverse residual is `1`, above the requested tolerance `2`. Falling back to a direct inverse in the reduced subspace.";

BuildMagneticTransverseBasis::inversefail =
    "The final inverse residual is `1`, still above the requested tolerance `2`. Check the equilibrium, projection cutoff, and numerical precision.";

(*
    Identify the well-resolved transverse magnetic subspace.

    For a real equilibrium texture, C = M0 x is anti-Hermitian.
    Therefore Hc = i C is Hermitian and has real eigenvalues kappa.
    Well-resolved transverse states have |kappa| close to 1, while
    longitudinal and cutoff-induced states lie closer to zero.

    The retained eigenvectors form the columns of Basis.

    In this basis:

        C_perp      = -i DiagonalMatrix[kappa]
        C_perp^-1   =  i DiagonalMatrix[1/kappa]

    The inverse is checked explicitly. If finite numerical errors
    make the spectral formula insufficiently accurate for the
    projected original operator, the code falls back to a direct
    inverse in the reduced transverse space.
*)
BuildMagneticTransverseBasis[
    crossOperator_?MatrixQ,
    OptionsPattern[]
] :=
    Module[
        {
            projectionCutoff,
            checkInverse,
            inverseTolerance,

            hermitianCrossOperator,
            eigenvalues,
            eigenvectors,

            selectedIndices,
            selectedEigenvalues,
            selectedVectors,

            basis,
            reducedCrossOperator,
            spectralInverse,
            inverseCrossOperator,

            spectralResidual,
            inverseResidual,
            inverseMethod,

            dimension
        },


        (* ------------------------------------------------------- *)
        (* Options                                                 *)
        (* ------------------------------------------------------- *)

        projectionCutoff =
            OptionValue["ProjectionCutoff"];

        checkInverse =
            OptionValue["CheckInverse"];

        inverseTolerance =
            OptionValue["InverseTolerance"];


        (* ------------------------------------------------------- *)
        (* Hermitian representation of M0 x                       *)
        (* ------------------------------------------------------- *)

        (*
            For a real equilibrium texture, M0 x is
            anti-Hermitian.

            Therefore

                i (M0 x)

            is Hermitian and can be diagonalized using
            Eigensystem with an orthonormal eigenbasis.
        *)
        hermitianCrossOperator =
            I crossOperator;


        {
            eigenvalues,
            eigenvectors
        } =
            Eigensystem[
                hermitianCrossOperator
            ];


        (* ------------------------------------------------------- *)
        (* Select well-resolved transverse states                 *)
        (* ------------------------------------------------------- *)

        (*
            Physical transverse states have eigenvalues
            approaching +/-1 in the complete Fourier basis.

            Finite Fourier truncation shifts some states away
            from these ideal values. Only sufficiently
            well-resolved states are retained.
        *)
        selectedIndices =
            Select[
                Range[
                    Length[eigenvalues]
                ],

                Abs[
                    eigenvalues[[#]]
                ] > projectionCutoff &
            ];


        selectedEigenvalues =
            eigenvalues[[
                selectedIndices
            ]];


        selectedVectors =
            eigenvectors[[
                selectedIndices
            ]];


        (* Eigenvectors returned by the Hermitian problem form
           an orthonormal basis. Store them as columns. *)
        basis =
            Transpose[
                selectedVectors
            ];


        (* ------------------------------------------------------- *)
        (* M0 x in the transverse basis                           *)
        (* ------------------------------------------------------- *)

        (*
            If

                i C v_n = kappa_n v_n,

            with C = M0 x, then

                C v_n = -i kappa_n v_n.

            Hence C is diagonal in this basis.
        *)
        reducedCrossOperator =
            DiagonalMatrix[
                -I selectedEigenvalues
            ];


        (* ------------------------------------------------------- *)
        (* Inverse of projected M0 x                              *)
        (* ------------------------------------------------------- *)

        (*
            Since the longitudinal sector has already been
            removed, all retained eigenvalues are nonzero.

            Therefore

                C^(-1) v_n = i/kappa_n v_n.
        *)
        spectralInverse =
            DiagonalMatrix[
                I / selectedEigenvalues
            ];


        inverseCrossOperator =
            spectralInverse;

        inverseMethod =
            "Spectral";


        (* ------------------------------------------------------- *)
        (* Optional validation of the inverse                     *)
        (* ------------------------------------------------------- *)

        spectralResidual =
            Missing["NotChecked"];

        inverseResidual =
            Missing["NotChecked"];


        If[
            TrueQ[checkInverse],

            dimension =
                Length[
                    selectedEigenvalues
                ];


            (* Check both left and right inverse relations. *)
            spectralResidual =
                Max[
                    Norm[
                        reducedCrossOperator .
                        spectralInverse
                        -
                        IdentityMatrix[dimension]
                    ],

                    Norm[
                        spectralInverse .
                        reducedCrossOperator
                        -
                        IdentityMatrix[dimension]
                    ]
                ];


            (*
                If the spectral inverse fails the requested
                tolerance, fall back to a direct inversion of
                the reduced operator.
            *)
            If[
                spectralResidual > inverseTolerance,

                inverseCrossOperator =
                    Inverse[
                        reducedCrossOperator
                    ];

                inverseMethod =
                    "Direct";
            ];


            inverseResidual =
                Max[
                    Norm[
                        reducedCrossOperator .
                        inverseCrossOperator
                        -
                        IdentityMatrix[dimension]
                    ],

                    Norm[
                        inverseCrossOperator .
                        reducedCrossOperator
                        -
                        IdentityMatrix[dimension]
                    ]
                ];
        ];


        (* ------------------------------------------------------- *)
        (* Return projection data                                 *)
        (* ------------------------------------------------------- *)

        <|
            "Basis" ->
                basis,

            "CrossOperator" ->
                reducedCrossOperator,

            "InverseCrossOperator" ->
                inverseCrossOperator,

            "AllCrossEigenvalues" ->
                eigenvalues,

            "SelectedEigenvalues" ->
                selectedEigenvalues,

            "ProjectionCutoff" ->
                projectionCutoff,

            "FullDimension" ->
                Length[crossOperator],

            "TransverseDimension" ->
                Length[selectedIndices],

            "InverseMethod" ->
                inverseMethod,

            "SpectralInverseResidual" ->
                spectralResidual,

            "InverseResidual" ->
                inverseResidual
        |>
    ];

(* --------------------------------------------------------------- *)
(* Transformations to and from the transverse basis               *)
(* --------------------------------------------------------------- *)

(* Project a full magnetic fluctuation vector onto the transverse basis. *)
ProjectMagneticVector[
    transverseData_Association,
    vector_List
] :=
    ConjugateTranspose[
        transverseData["Basis"]
    ] . vector;

(* Lift a transverse magnetic vector back to the full Fourier space. *)
LiftMagneticVector[
    transverseData_Association,
    vector_List
] :=
    transverseData["Basis"] . vector;

(* Project a magnetic operator onto the transverse magnetic subspace. *)
ProjectMagneticOperator[
    transverseData_Association,
    operator_?MatrixQ
] :=
    Chop[
        ConjugateTranspose[
            transverseData["Basis"]
        ] .
        operator .
        transverseData["Basis"]
    ];

(* --------------------------------------------------------------- *)
(* Help                                                            *)
(* --------------------------------------------------------------- *)

DynamicsProjectionHelp[] :=
    Column[
        {
            Style[
                "DynamicsProjection.wl",
                Bold,
                16
            ],

            "Common reciprocal-space and transverse-projection tools for magnetic and magnetoelastic dynamics.",

            Grid[
                {
                    {"Function", "Description"},
                    {"DynamicsNodes", "Return all reciprocal-lattice nodes inside the dynamical hexagonal cutoff."},
                    {"EquilibriumMagnetization", "Return a Fourier component of a saved equilibrium magnetization."},
                    {"CrossProductMatrix", "Return the Cartesian matrix representing m x v."},
                    {"BuildM0CrossOperator", "Construct the full Fourier-space M0 x operator."},
                    {"BuildMagneticTransverseBasis", "Construct the well-resolved transverse magnetic basis, projected M0 x operator and verified inverse."},
                    {"ProjectMagneticVector", "Project a full magnetic fluctuation vector onto the transverse basis."},
                    {"LiftMagneticVector", "Lift a transverse magnetic vector back to the full Fourier space."},
                    {"ProjectMagneticOperator", "Project a magnetic operator onto the transverse magnetic subspace."}
                },
                Frame -> All,
                Alignment -> Left
            ]
        }
    ];

Print[
    DynamicsProjectionHelp[]
];