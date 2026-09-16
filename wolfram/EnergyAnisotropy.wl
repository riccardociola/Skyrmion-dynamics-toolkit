(* ================================================================ *)
(* EnergyAnisotropy.wl                                             *)
(*                                                                  *)
(* Magnetocrystalline-anisotropy contributions to the magnetic     *)
(* free energy.                                                     *)
(*                                                                  *)
(* Coordinate convention                                           *)
(* ---------------------                                           *)
(* The simulation Cartesian z axis is always chosen parallel to    *)
(* the external magnetic field and to the skyrmion tubes.          *)
(*                                                                  *)
(* The user specifies the field direction in cubic crystal         *)
(* coordinates, e.g.                                               *)
(*                                                                  *)
(*     {0,0,1}   -> H || [001]                                     *)
(*     {1,1,0}   -> H || [110]                                     *)
(*     {1,1,1}   -> H || [111]                                     *)
(*                                                                  *)
(* or any other nonzero 3-vector.                                  *)
(*                                                                  *)
(* Cubically equivalent field directions [100], [010], and [001]   *)
(* are represented internally as [001]. The same is done for their *)
(* opposite directions, since the present cubic anisotropy is even *)
(* under reversal of a cubic axis.                                 *)
(*                                                                  *)
(* For a general field direction, [100] is used as the fixed       *)
(* reference axis defining alphaR = 0. Its projection onto the     *)
(* plane perpendicular to H defines the simulation x direction.    *)
(* The y direction completes a right-handed Cartesian frame.       *)
(*                                                                  *)
(* alphaR is NOT fixed by this file. It enters through the affine   *)
(* transformation and can therefore be treated as a minimization   *)
(* variable once cubic anisotropy is included.                     *)
(*                                                                  *)
(* Cubic anisotropy convention                                     *)
(* ----------------------------                                    *)
(*                                                                  *)
(*     Ecub = kappaCub (m_[100]^4 + m_[010]^4 + m_[001]^4)         *)
(*                                                                  *)
(* where the three magnetization components are taken along the    *)
(* cubic crystallographic axes.                                    *)
(*                                                                  *)
(* Dependencies:                                                    *)
(*     Geometry.wl                                                 *)
(*     FourierBasis.wl                                             *)
(* ================================================================ *)

ClearAll[
    CanonicalCubicFieldDirection,
    CubicLabFrame,
    CubicCrystalAxes,
    AnisotropyFourierNodes,
    CubicAnisotropyEnergy,
    EnergyAnisotropyHelp
];

(* --------------------------------------------------------------- *)
(* Messages                                                        *)
(* --------------------------------------------------------------- *)

CanonicalCubicFieldDirection::dim =
    "Field direction must be a three-component vector. Received `1`.";

CanonicalCubicFieldDirection::zero =
    "Field direction must be nonzero.";

CubicLabFrame::frame =
    "Could not construct a stable Cartesian frame for field direction `1`.";

CubicAnisotropyEnergy::rep =
    "Unknown magnetic representation `1`. Expected \"SixFold\" or \"General\".";

(* --------------------------------------------------------------- *)
(* Field-direction convention                                     *)
(* --------------------------------------------------------------- *)

ClearAll[cubicParallelQ];

cubicParallelQ[
    v1_List,
    v2_List,
    tolerance_: 10^-12
] :=
    Module[
        {cross},

        cross =
            Cross[
                v1,
                v2
            ];

        If[
            VectorQ[cross, NumericQ],

            Norm[N[cross]] < tolerance,

            TrueQ[
                Simplify[
                    cross == {0, 0, 0}
                ]
            ]
        ]
    ];

(*
    Normalize the field direction in cubic crystal coordinates.

    Any field parallel to one of the three cubic coordinate axes is
    represented internally as {0,0,1}. This makes H || [100],
    H || [010], and H || [001] explicitly equivalent in the code.
*)
CanonicalCubicFieldDirection[
    fieldDirection_List
] :=
    Module[
        {
            normalizedDirection,
            norm,
            cubicAxes
        },

        If[
            Length[fieldDirection] != 3,

            Message[
                CanonicalCubicFieldDirection::dim,
                fieldDirection
            ];

            Return[$Failed]
        ];

        norm =
            Sqrt[
                fieldDirection . fieldDirection
            ];

        If[
            TrueQ[Simplify[norm == 0]]
            ||
            (
                NumericQ[N[norm]]
                &&
                Chop[N[norm]] == 0
            ),

            Message[
                CanonicalCubicFieldDirection::zero
            ];

            Return[$Failed]
        ];

        normalizedDirection =
            Simplify[
                fieldDirection / norm
            ];

        cubicAxes =
            {
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1}
            };

        If[
            AnyTrue[
                cubicAxes,
                cubicParallelQ[
                    normalizedDirection,
                    #
                ] &
            ],

            {0, 0, 1},

            normalizedDirection
        ]
    ];

(* --------------------------------------------------------------- *)
(* Crystal orientation relative to the simulation frame           *)
(* --------------------------------------------------------------- *)

(*
    Construct the simulation Cartesian frame in crystal coordinates.

    zLab:
        parallel to the external magnetic field.

    xLab:
        normalized projection of the cubic [100] axis onto the plane
        perpendicular to the field.

    yLab:
        chosen so that {xLab,yLab,zLab} is right handed.

    alphaR is deliberately absent here. It rotates the AFFINE SkL
    basis relative to this fixed Cartesian simulation frame and is
    handled later by toAff[...,eta,phi,alphaR].
*)
CubicLabFrame[
    fieldDirection_List
] :=
    Module[
        {
            zLabCrystal,
            referenceAxisCrystal,
            xProjection,
            xProjectionNorm,
            xLabCrystal,
            yLabCrystal,
            crystalToLabMatrix,
            labToCrystalMatrix
        },

        zLabCrystal =
            CanonicalCubicFieldDirection[
                fieldDirection
            ];

        If[
            zLabCrystal === $Failed,
            Return[$Failed]
        ];

        referenceAxisCrystal =
            {1, 0, 0};

        xProjection =
            Simplify[
                referenceAxisCrystal
                -
                (
                    referenceAxisCrystal . zLabCrystal
                )
                zLabCrystal
            ];

        xProjectionNorm =
            Simplify[
                Sqrt[
                    xProjection . xProjection
                ]
            ];

        If[
            TrueQ[Simplify[xProjectionNorm == 0]]
            ||
            (
                NumericQ[N[xProjectionNorm]]
                &&
                Chop[N[xProjectionNorm]] == 0
            ),

            Message[
                CubicLabFrame::frame,
                fieldDirection
            ];

            Return[$Failed]
        ];

        xLabCrystal =
            Simplify[
                xProjection / xProjectionNorm
            ];

        yLabCrystal =
            Simplify[
                Cross[
                    zLabCrystal,
                    xLabCrystal
                ]
            ];

        crystalToLabMatrix =
            Simplify[
                {
                    xLabCrystal,
                    yLabCrystal,
                    zLabCrystal
                }
            ];

        labToCrystalMatrix =
            Simplify[
                Transpose[
                    crystalToLabMatrix
                ]
            ];

        <|
            "InputFieldDirectionCrystal" ->
                fieldDirection,

            "FieldDirectionCrystal" ->
                zLabCrystal,

            "ReferenceAxisCrystal" ->
                referenceAxisCrystal,

            "XLabCrystal" ->
                xLabCrystal,

            "YLabCrystal" ->
                yLabCrystal,

            "ZLabCrystal" ->
                zLabCrystal,

            "CrystalToLabMatrix" ->
                crystalToLabMatrix,

            "LabToCrystalMatrix" ->
                labToCrystalMatrix
        |>
    ];

(*
    Return the cubic [100], [010], and [001] axes expressed in the
    simulation Cartesian frame.
*)
CubicCrystalAxes[
    fieldDirection_List
] :=
    Module[
        {
            frame,
            crystalToLabMatrix,
            crystalAxes
        },

        frame =
            CubicLabFrame[
                fieldDirection
            ];

        If[
            frame === $Failed,
            Return[$Failed]
        ];

        crystalToLabMatrix =
            frame[
                "CrystalToLabMatrix"
            ];

        crystalAxes =
            {
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1}
            };

        Simplify[
            (
                crystalToLabMatrix . #
            )
            &
            /@
            crystalAxes
        ]
    ];

(* --------------------------------------------------------------- *)
(* Reciprocal-space support                                       *)
(* --------------------------------------------------------------- *)

(*
    Reciprocal-lattice nodes inside the same hexagonal Fourier cutoff
    used by FourierBasis.wl.

    This imposes only the cutoff geometry. It does NOT impose sixfold
    symmetry on the magnetic Fourier amplitudes.
*)
AnisotropyFourierNodes[
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
        ]
        &
    ];

(* --------------------------------------------------------------- *)
(* Cubic magnetocrystalline anisotropy                            *)
(* --------------------------------------------------------------- *)

(*
    Cubic anisotropy energy

        Ecub = kappaCub Sum_i < (m . vi)^4 >

    where vi are the cubic [100], [010], [001] axes expressed in the
    simulation Cartesian frame.

    Each crystal axis is converted to covariant affine components via

        toAff[vi,eta,phi,alphaR]

    so alphaR explicitly controls the in-plane orientation of the SkL
    relative to the fixed crystallographic axes.

    The Fourier-space quartic convolution is evaluated through a pair
    convolution to avoid the direct cubic loop over reciprocal nodes.

    representation:
        "SixFold" -> use the sixfold magnetic parametrization;
        "General" -> use the fully general complex magnetic Fourier
                     parametrization.

    Automatic:
        use $MagneticRepresentation when configured;
        otherwise default to "SixFold".

    This function never changes the representation automatically when
    cubic anisotropy is active.
*)
CubicAnisotropyEnergy[
    kappaCub_,
    fieldDirection_List,
    cutoff_Integer?NonNegative,
    eta_,
    phi_,
    alphaR_,
    representation_: Automatic
] :=
    Module[
        {
            resolvedRepresentation,
            nodes,
            crystalAxesCartesian,
            crystalAxesAffine,
            nodeKey,
            axisEnergy
        },

        resolvedRepresentation =
            Replace[
                representation,

                Automatic :>
                    If[
                        ValueQ[$MagneticRepresentation],
                        $MagneticRepresentation,
                        "SixFold"
                    ]
            ];

        If[
            ! MemberQ[
                {"SixFold", "General"},
                resolvedRepresentation
            ],

            Message[
                CubicAnisotropyEnergy::rep,
                resolvedRepresentation
            ];

            Return[$Failed]
        ];

        nodes =
            AnisotropyFourierNodes[
                cutoff
            ];

        crystalAxesCartesian =
            CubicCrystalAxes[
                fieldDirection
            ];

        If[
            crystalAxesCartesian === $Failed,
            Return[$Failed]
        ];

        crystalAxesAffine =
            (
                toAff[
                    #,
                    eta,
                    phi,
                    alphaR
                ]
            )
            &
            /@
            crystalAxesCartesian;

        (*
            Compact unique key for node {n1,n2}. Integer reciprocal
            indices make n1 + I n2 unambiguous.
        *)
        nodeKey[
            node_List
        ] :=
            node[[1]] + I node[[2]];

        axisEnergy[
            axisAffine_List
        ] :=
            Module[
                {
                    projectionAssociation,
                    pairRules,
                    pairConvolutionAssociation
                },

                projectionAssociation =
                    Association[
                        Table[
                            nodeKey[node]
                            ->
                            DotC[
                                MagneticFourierComponent[
                                    node[[1]],
                                    node[[2]],
                                    cutoff,
                                    resolvedRepresentation
                                ],

                                axisAffine,
                                eta,
                                phi,
                                alphaR
                            ],

                            {node, nodes}
                        ]
                    ];

                pairRules =
                    Flatten[
                        Table[
                            nodeKey[
                                node1 + node2
                            ]
                            ->
                            (
                                Lookup[
                                    projectionAssociation,
                                    nodeKey[node1],
                                    0
                                ]

                                Lookup[
                                    projectionAssociation,
                                    nodeKey[node2],
                                    0
                                ]
                            ),

                            {node1, nodes},
                            {node2, nodes}
                        ],
                        1
                    ];

                pairConvolutionAssociation =
                    Merge[
                        pairRules,
                        Total
                    ];

                Total[
                    KeyValueMap[
                        Function[
                            {qKey, value},

                            value
                            Lookup[
                                pairConvolutionAssociation,
                                -qKey,
                                0
                            ]
                        ],

                        pairConvolutionAssociation
                    ]
                ]
            ];

        kappaCub
        Total[
            axisEnergy
                /@
            crystalAxesAffine
        ]
    ];

(* --------------------------------------------------------------- *)
(* Help                                                            *)
(* --------------------------------------------------------------- *)

EnergyAnisotropyHelp[] :=
    Column[
        {
            Style[
                "EnergyAnisotropy.wl",
                Bold,
                16
            ],

            "Cubic magnetocrystalline-anisotropy tools for the skyrmion-lattice model.",

            Grid[
                {
                    {"Function", "Description"},
                    {"CanonicalCubicFieldDirection", "Normalize the field direction in crystal coordinates and map the equivalent [100], [010], and [001] directions to the internal [001] convention."},
                    {"CubicLabFrame", "Construct the fixed simulation Cartesian frame from the field direction using [100] as the reference direction for alphaR = 0."},
                    {"CubicCrystalAxes", "Return the cubic [100], [010], and [001] axes expressed in the simulation Cartesian frame."},
                    {"AnisotropyFourierNodes", "Return all reciprocal-lattice nodes inside the hexagonal Fourier cutoff without imposing sixfold symmetry on the Fourier amplitudes."},
                    {"CubicAnisotropyEnergy", "Calculate the quartic cubic anisotropy energy for either the SixFold or General magnetic Fourier representation."}
                },
                Frame -> All,
                Alignment -> Left
            ],

            "Example field directions:",

            Grid[
                {
                    {"H || [001]", "{0,0,1}"},
                    {"H || [110]", "{1,1,0}"},
                    {"H || [111]", "{1,1,1}"},
                    {"Arbitrary direction", "{h,k,l}"}
                },
                Frame -> All,
                Alignment -> Left
            ]
        }
    ];

Print[
    EnergyAnisotropyHelp[]
];