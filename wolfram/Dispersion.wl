(* ================================================================ *)
(* Dispersion.wl                                                    *)
(*                                                                  *)
(* Magnon dispersion calculations and mode tracking.               *)
(*                                                                  *)
(* Implemented dispersions:                                        *)
(*     - longitudinal: k = {0,0,qz};                               *)
(*     - transverse: Gamma-M-K-Gamma in the first magnetic BZ.     *)
(*                                                                  *)
(* All eigensystems are calculated independently in parallel.      *)
(* Mode identity is fixed at Gamma and then propagated locally by  *)
(* comparing BdG overlaps only between neighboring k points.       *)
(*                                                                  *)
(* Dependencies:                                                    *)
(*     Geometry.wl                                                 *)
(*     GroundState.wl                                              *)
(*     DynamicsMagnetic.wl                                         *)
(*     ModeAnalysis.wl                                             *)
(* ================================================================ *)

ClearAll[
    DispersionModeColor,
    MagneticBdGOverlapMatrix,
    MatchNeighboringModes,
    PropagateModeIdentification,

    FirstBrillouinZoneData,
    GammaMKPathData,

    LongitudinalMagnonSpectrum,
    TransverseMagnonSpectrum,

    DispersionHelp
];

(* --------------------------------------------------------------- *)
(* Messages                                                        *)
(* --------------------------------------------------------------- *)

LongitudinalMagnonSpectrum::qmax =
    "qzMax must be a positive number. Received `1`.";

LongitudinalMagnonSpectrum::omax =
    "omegaMax must be a positive number. Received `1`.";

LongitudinalMagnonSpectrum::npoints =
    "nPoints must be an integer greater than or equal to 3. Received `1`.";

LongitudinalMagnonSpectrum::odd =
    "nPoints = `1` is even. It has been increased to `2` so that qz = 0 is included exactly in the uniform mesh.";

LongitudinalMagnonSpectrum::nomodes =
    "No positive-frequency Gamma-point modes were found below omegaMax = `1`.";

LongitudinalMagnonSpectrum::parallel =
    "At least one longitudinal eigensystem calculation failed.";

FirstBrillouinZoneData::geometry =
    "The first Brillouin zone could not be reconstructed from the six nearest reciprocal-lattice vectors.";

TransverseMagnonSpectrum::omax =
    "omegaMax must be a positive number. Received `1`.";

TransverseMagnonSpectrum::npoints =
    "nPoints must be an integer greater than or equal to 4. Received `1`.";

TransverseMagnonSpectrum::nomodes =
    "No positive-frequency Gamma-point modes were found below omegaMax = `1`.";

TransverseMagnonSpectrum::parallel =
    "At least one transverse eigensystem calculation failed.";

(* --------------------------------------------------------------- *)
(* Publication color cycle                                        *)
(* --------------------------------------------------------------- *)

(*
    Six-color publication-friendly cycle based on the Okabe-Ito
    colorblind-safe palette.

    Colors repeat every six modes:
        mode 1 -> color 1
        ...
        mode 6 -> color 6
        mode 7 -> color 1
        ...
*)
$DispersionColors = {
    RGBColor["#0072B2"],   (* blue *)
    RGBColor["#D55E00"],   (* vermillion *)
    RGBColor["#009E73"],   (* bluish green *)
    RGBColor["#CC79A7"],   (* reddish purple *)
    RGBColor["#E69F00"],   (* orange *)
    RGBColor["#56B4E9"]    (* sky blue *)
};

DispersionModeColor[
    modeIndex_Integer?Positive
] :=
    $DispersionColors[[
        Mod[
            modeIndex - 1,
            Length[$DispersionColors]
        ]
        + 1
    ]];

(* --------------------------------------------------------------- *)
(* Local BdG mode tracking                                        *)
(* --------------------------------------------------------------- *)

(*
    Calculate the squared magnetic BdG overlap between two sets of
    transverse eigenvectors.

    Rows correspond to the previously identified branches.
    Columns correspond to the raw modes at the neighboring k point.
*)
MagneticBdGOverlapMatrix[
    dynamicsData_Association,
    previousVectors_List,
    currentVectors_List
] :=
    Module[
        {
            transverseData
        },

        transverseData =
            dynamicsData["TransverseData"];


        Table[
            Abs[
                MagneticBdGInnerProduct[
                    transverseData,
                    previousVectors[[previousIndex]],
                    currentVectors[[currentIndex]]
                ]
            ]^2,

            {
                previousIndex,
                Length[previousVectors]
            },

            {
                currentIndex,
                Length[currentVectors]
            }
        ]
    ];

(*
    Match the previously identified branches to the raw modes at the
    next k point.

    Matching is one-to-one. All possible BdG overlaps are sorted from
    largest to smallest and the strongest still-available pair is
    accepted first.

    This avoids assigning the same current eigenmode to two branches,
    while remaining negligible in cost compared with diagonalization.
*)
MatchNeighboringModes[
    dynamicsData_Association,
    previousVectors_List,
    currentVectors_List
] :=
    Module[
        {
            numberOfPrevious,
            numberOfCurrent,

            overlapMatrix,
            rankedPairs,

            assignedCurrentIndices,
            assignedOverlaps,
            currentUsed,

            overlapValue,
            previousIndex,
            currentIndex
        },


        numberOfPrevious =
            Length[
                previousVectors
            ];

        numberOfCurrent =
            Length[
                currentVectors
            ];


        If[
            numberOfPrevious == 0,

            Return[
                <|
                    "OverlapMatrix" -> {},
                    "CurrentIndices" -> {},
                    "Overlaps" -> {}
                |>
            ]
        ];


        If[
            numberOfCurrent == 0,

            Return[
                <|
                    "OverlapMatrix" -> {},
                    "CurrentIndices" ->
                        ConstantArray[
                            Missing["NotAssigned"],
                            numberOfPrevious
                        ],
                    "Overlaps" ->
                        ConstantArray[
                            Missing["NotAssigned"],
                            numberOfPrevious
                        ]
                |>
            ]
        ];


        overlapMatrix =
            MagneticBdGOverlapMatrix[
                dynamicsData,
                previousVectors,
                currentVectors
            ];


        rankedPairs =
            Reverse[
                SortBy[
                    Flatten[
                        Table[
                            {
                                overlapMatrix[[
                                    previousIndex,
                                    currentIndex
                                ]],
                                previousIndex,
                                currentIndex
                            },

                            {
                                previousIndex,
                                numberOfPrevious
                            },

                            {
                                currentIndex,
                                numberOfCurrent
                            }
                        ],
                        1
                    ],
                    First
                ]
            ];


        assignedCurrentIndices =
            ConstantArray[
                Missing["NotAssigned"],
                numberOfPrevious
            ];

        assignedOverlaps =
            ConstantArray[
                Missing["NotAssigned"],
                numberOfPrevious
            ];

        currentUsed =
            ConstantArray[
                False,
                numberOfCurrent
            ];


        Do[
            overlapValue =
                pair[[1]];

            previousIndex =
                pair[[2]];

            currentIndex =
                pair[[3]];


            If[
                MissingQ[
                    assignedCurrentIndices[[
                        previousIndex
                    ]]
                ]
                &&
                ! TrueQ[
                    currentUsed[[
                        currentIndex
                    ]]
                ],

                assignedCurrentIndices[[
                    previousIndex
                ]] =
                    currentIndex;

                assignedOverlaps[[
                    previousIndex
                ]] =
                    overlapValue;

                currentUsed[[
                    currentIndex
                ]] =
                    True;
            ],

            {
                pair,
                rankedPairs
            }
        ];


        <|
            "OverlapMatrix" ->
                overlapMatrix,

            "CurrentIndices" ->
                assignedCurrentIndices,

            "Overlaps" ->
                assignedOverlaps
        |>
    ];

(*
    Propagate Gamma-point mode identity through an ORDERED list of
    neighboring eigensystems.

    rawModes[[1]] must be the Gamma-point eigensystem.

    referenceModeIndices specifies which Gamma modes are followed.
    These are normally the modes whose Gamma-point frequency lies
    below the requested plotting limit omegaMax.

    At each new k point the comparison is made with the eigenvectors
    selected at the immediately preceding point, not with Gamma.

    If a branch is temporarily unavailable, its most recently known
    eigenvector is retained only as a recovery reference for the next
    point. The branch frequency at the missing point remains Missing.
*)
PropagateModeIdentification[
    dynamicsData_Association,
    rawModes_List,
    referenceModeIndices_List
] :=
    Module[
        {
            gammaModes,
            numberOfBranches,

            previousVectors,
            trackedData,

            currentModes,
            currentVectors,
            match,

            currentIndices,
            currentFrequencies,
            nextPreviousVectors,

            branchIndex,
            currentIndex
        },


        If[
            rawModes === {},

            Return[{}]
        ];


        gammaModes =
            First[
                rawModes
            ];

        numberOfBranches =
            Length[
                referenceModeIndices
            ];


        previousVectors =
            gammaModes[
                "EigenvectorsTransverse"
            ][[
                referenceModeIndices
            ]];


        trackedData =
            {
                <|
                    "CurrentIndices" ->
                        referenceModeIndices,

                    "Frequencies" ->
                        gammaModes[
                            "Frequencies"
                        ][[
                            referenceModeIndices
                        ]],

                    "Overlaps" ->
                        ConstantArray[
                            1.,
                            numberOfBranches
                        ],

                    "OverlapMatrix" ->
                        IdentityMatrix[
                            numberOfBranches
                        ]
                |>
            };


        Do[
            currentModes =
                rawModes[[
                    pointIndex
                ]];

            currentVectors =
                currentModes[
                    "EigenvectorsTransverse"
                ];


            match =
                MatchNeighboringModes[
                    dynamicsData,
                    previousVectors,
                    currentVectors
                ];


            currentIndices =
                match[
                    "CurrentIndices"
                ];


            currentFrequencies =
                Table[
                    currentIndex =
                        currentIndices[[
                            branchIndex
                        ]];


                    If[
                        MissingQ[
                            currentIndex
                        ],

                        Missing[
                            "NotTracked"
                        ],

                        currentModes[
                            "Frequencies"
                        ][[
                            currentIndex
                        ]]
                    ],

                    {
                        branchIndex,
                        numberOfBranches
                    }
                ];


            nextPreviousVectors =
                Table[
                    currentIndex =
                        currentIndices[[
                            branchIndex
                        ]];


                    If[
                        MissingQ[
                            currentIndex
                        ],

                        previousVectors[[
                            branchIndex
                        ]],

                        currentVectors[[
                            currentIndex
                        ]]
                    ],

                    {
                        branchIndex,
                        numberOfBranches
                    }
                ];


            AppendTo[
                trackedData,

                <|
                    "CurrentIndices" ->
                        currentIndices,

                    "Frequencies" ->
                        currentFrequencies,

                    "Overlaps" ->
                        match[
                            "Overlaps"
                        ],

                    "OverlapMatrix" ->
                        match[
                            "OverlapMatrix"
                        ]
                |>
            ];


            previousVectors =
                nextPreviousVectors,

            {
                pointIndex,
                2,
                Length[
                    rawModes
                ]
            }
        ];


        trackedData
    ];

(* --------------------------------------------------------------- *)
(* Longitudinal magnon spectrum                                   *)
(* --------------------------------------------------------------- *)

Options[LongitudinalMagnonSpectrum] = {
    ImageSize -> 800,
    "ShowLegend" -> True,
    "SymmetryTolerance" -> 10^-6
};

(*
    Calculate the longitudinal magnon spectrum for k = {0,0,qz}.

    All mesh-point eigensystems, including Gamma, are evaluated
    independently with ParallelTable.

    The Gamma point fixes:
        - mode number;
        - color;
        - angular-momentum label m.

    Mode identity is then propagated locally:
        Gamma -> +qz
        Gamma -> -qz

    so every state is compared with the nearest previously identified
    wavevector state rather than directly with Gamma.

    The eigensystem is solved internally up to 1.15 omegaMax.
    The plot itself is restricted to omegaMax.
*)
LongitudinalMagnonSpectrum[
    dynamicsData_Association,
    qzMax_: 1,
    omegaMax_: 3,
    nPoints_: 401,
    OptionsPattern[]
] :=
    Module[
        {
            actualNPoints,
            internalOmegaMax,

            imageSize,
            showLegend,
            symmetryTolerance,

            qzValues,
            zeroIndex,

            parallelDynamicsData,
            rawQData,
            rawModes,

            gammaModes,
            gammaFrequencies,
            referenceModeIndices,

            angularAnalysis,
            angularMomentumIndices,

            positiveTracking,
            negativeTrackingOutward,
            trackedData,

            branches,
            branchFrequencies,

            modeColors,
            modeLabels,
            plot
        },


        (* Validate inputs. *)
        If[
            ! NumericQ[qzMax]
            ||
            qzMax <= 0,

            Message[
                LongitudinalMagnonSpectrum::qmax,
                qzMax
            ];

            Return[$Failed]
        ];


        If[
            ! NumericQ[omegaMax]
            ||
            omegaMax <= 0,

            Message[
                LongitudinalMagnonSpectrum::omax,
                omegaMax
            ];

            Return[$Failed]
        ];


        If[
            ! IntegerQ[nPoints]
            ||
            nPoints < 3,

            Message[
                LongitudinalMagnonSpectrum::npoints,
                nPoints
            ];

            Return[$Failed]
        ];


        actualNPoints =
            If[
                OddQ[nPoints],

                nPoints,

                Message[
                    LongitudinalMagnonSpectrum::odd,
                    nPoints,
                    nPoints + 1
                ];

                nPoints + 1
            ];


        internalOmegaMax =
            1.15
            omegaMax;


        imageSize =
            OptionValue[
                ImageSize
            ];

        showLegend =
            OptionValue[
                "ShowLegend"
            ];

        symmetryTolerance =
            OptionValue[
                "SymmetryTolerance"
            ];


        (* Uniform symmetric mesh with exact Gamma point. *)
        qzValues =
            N[
                Subdivide[
                    -qzMax,
                    qzMax,
                    actualNPoints - 1
                ]
            ];


        zeroIndex =
            Quotient[
                actualNPoints + 1,
                2
            ];

        qzValues[[
            zeroIndex
        ]] =
            0;


        (* Data required by SolveMagneticModes on parallel kernels. *)
        parallelDynamicsData =
            KeyTake[
                dynamicsData,
                {
                    "Equilibrium",
                    "Nodes",
                    "TransverseData",
                    "ConstraintHessianTransverse"
                }
            ];


        DistributeDefinitions[
            parallelDynamicsData,
            SolveMagneticModes
        ];


        (* Every qz point is solved independently in parallel. *)
        rawQData =
            ParallelTable[
                <|
                    "Qz" ->
                        qz,

                    "Modes" ->
                        SolveMagneticModes[
                            parallelDynamicsData,
                            {
                                0,
                                0,
                                qz
                            },

                            "FrequencyRange" ->
                                {
                                    0,
                                    internalOmegaMax
                                }
                        ]
                |>,

                {
                    qz,
                    qzValues
                }
            ];


        If[
            MemberQ[
                Lookup[
                    rawQData,
                    "Modes"
                ],
                $Failed
            ],

            Message[
                LongitudinalMagnonSpectrum::parallel
            ];

            Return[$Failed]
        ];


        rawModes =
            Lookup[
                rawQData,
                "Modes"
            ];


        (* Gamma-point reference modes. *)
        gammaModes =
            rawModes[[
                zeroIndex
            ]];

        gammaFrequencies =
            gammaModes[
                "Frequencies"
            ];


        referenceModeIndices =
            Select[
                Range[
                    Length[
                        gammaFrequencies
                    ]
                ],

                gammaFrequencies[[#]]
                    <=
                omegaMax
                &
            ];


        If[
            referenceModeIndices === {},

            Message[
                LongitudinalMagnonSpectrum::nomodes,
                omegaMax
            ];

            Return[$Failed]
        ];


        (* Gamma-point angular-momentum labels. *)
        angularAnalysis =
            AnalyzeAngularMomentumModes[
                dynamicsData,
                gammaModes,

                "SymmetryTolerance" ->
                    symmetryTolerance
            ];


        If[
            angularAnalysis === $Failed,
            Return[$Failed]
        ];


        angularMomentumIndices =
            Lookup[
                angularAnalysis[
                    "Modes"
                ],
                "AngularMomentumClassification"
            ];


        (*
            Propagate independently from Gamma toward positive and
            negative qz. Both tracking passes therefore begin from
            the same Gamma-point reference.
        *)
        positiveTracking =
            PropagateModeIdentification[
                dynamicsData,

                rawModes[[
                    zeroIndex ;;
                    Length[rawModes]
                ]],

                referenceModeIndices
            ];


        negativeTrackingOutward =
            PropagateModeIdentification[
                dynamicsData,

                Reverse[
                    rawModes[[
                        1 ;;
                        zeroIndex
                    ]]
                ],

                referenceModeIndices
            ];


        (*
            Restore ascending qz order:
                -qzMax ... 0 ... +qzMax
        *)
        trackedData =
            Join[
                Most[
                    Reverse[
                        negativeTrackingOutward
                    ]
                ],

                positiveTracking
            ];


        (* One plotted branch for every Gamma mode below omegaMax. *)
        branches =
            Table[
                branchFrequencies =
                    (
                        #[
                            "Frequencies"
                        ][[
                            branchIndex
                        ]]
                        &
                        /@
                        trackedData
                    );


                Cases[
                    Transpose[
                        {
                            qzValues,
                            branchFrequencies
                        }
                    ],

                    {
                        x_?NumericQ,
                        y_?NumericQ
                    }
                ],

                {
                    branchIndex,
                    Length[
                        referenceModeIndices
                    ]
                }
            ];


        modeColors =
            DispersionModeColor
                /@
            referenceModeIndices;


        modeLabels =
            Table[
                With[
                    {
                        referenceIndex =
                            referenceModeIndices[[
                                branchIndex
                            ]]
                    },

                    Row[
                        {
                            "Mode ",
                            referenceIndex,
                            "   (m = ",
                            If[
                                MissingQ[
                                    angularMomentumIndices[[
                                        referenceIndex
                                    ]]
                                ],

                                "?",

                                angularMomentumIndices[[
                                    referenceIndex
                                ]]
                            ],
                            ")"
                        }
                    ]
                ],

                {
                    branchIndex,
                    Length[
                        referenceModeIndices
                    ]
                }
            ];


        plot =
            ListLinePlot[
                branches,

                Joined ->
                    True,

                PlotStyle ->
                    (
                        Directive[
                            #,
                            Thick
                        ]
                        &
                        /@
                        modeColors
                    ),

                PlotRange ->
                    {
                        {
                            -qzMax,
                            qzMax
                        },

                        {
                            0,
                            omegaMax
                        }
                    },

                Frame ->
                    True,

                Axes ->
                    False,

                FrameLabel ->
                    {
                        Style[
                            Subscript[
                                        "q",
                                        "z"
                                    ],
                            15
                        ],

                        Style[
                            Row[
                                {
                                    "\[Omega]/",
                                    Subscript[
                                        "\[Omega]",
                                        "c2"
                                    ]
                                }
                            ],
                            15
                        ]
                    },

                LabelStyle ->
                    13,

                PlotLegends ->
                    If[
                        TrueQ[
                            showLegend
                        ],

                        Placed[
                            LineLegend[
                                modeColors,
                                modeLabels,
                                LegendLayout ->
                                    "Column"
                            ],
                            Right
                        ],

                        None
                    ],

                ImageSize ->
                    imageSize
            ];


        <|
            "QzMax" ->
                qzMax,

            "OmegaMax" ->
                omegaMax,

            "InternalOmegaMax" ->
                internalOmegaMax,

            "RequestedNumberOfPoints" ->
                nPoints,

            "NumberOfPoints" ->
                actualNPoints,

            "QzValues" ->
                qzValues,

            "GammaModes" ->
                gammaModes,

            "ReferenceModeIndices" ->
                referenceModeIndices,

            "AngularMomentumAnalysis" ->
                angularAnalysis,

            "AngularMomentumIndices" ->
                angularMomentumIndices,

            "RawData" ->
                rawQData,

            "TrackedData" ->
                trackedData,

            "Branches" ->
                branches,

            "ModeColors" ->
                modeColors,

            "ModeLabels" ->
                modeLabels,

            "Plot" ->
                plot
        |>
    ];

(* --------------------------------------------------------------- *)
(* First hexagonal Brillouin zone                                 *)
(* --------------------------------------------------------------- *)

(*
    Reconstruct the first magnetic Brillouin zone directly from the
    reciprocal lattice of the saved equilibrium.

    All returned wavevectors are Cartesian and use the same units as
    the Bloch wavevector k expected by SolveMagneticModes.

    The construction is the reciprocal-space Wigner-Seitz cell:

        1. generate nearby reciprocal-lattice vectors;
        2. retain the six shortest nonzero vectors;
        3. order them by polar angle;
        4. intersect the perpendicular bisectors of adjacent vectors.

    The selected M point is the midpoint of one BZ edge. The chosen
    edge is the one whose outward normal is closest to the positive
    Cartesian x direction.

    K is chosen as the counter-clockwise endpoint of that same edge.
    For a sixfold-symmetric SkL all equivalent choices are related by
    C6 symmetry.
*)
FirstBrillouinZoneData[
    dynamicsData_Association
] :=
    Module[
        {
            equilibrium,
            solutionRules,

            aValue,
            etaValue,
            phiValue,
            alphaValue,

            candidateNodes,
            candidateVectors,
            nearestVectors,
            orderedNeighbors,

            neighborAngles,
            vertices,
            orderedVertices,

            selectedNeighborIndex,
            selectedNeighbor,
            mPoint,

            edgeTolerance,
            adjacentVertices,
            kPoint,

            gammaPoint,
            reciprocalBasis
        },


        equilibrium =
            dynamicsData["Equilibrium"];

        solutionRules =
            EquilibriumRules[
                equilibrium
            ];

        aValue =
            a0 /. solutionRules;

        etaValue =
            eta0 /. solutionRules;

        phiValue =
            phi0 /. solutionRules;

        alphaValue =
            alphaR0 /. solutionRules;


        (* Primitive reciprocal basis in Cartesian coordinates. *)
        reciprocalBasis =
            {
                toCRT[
                    G[1, 0],
                    etaValue,
                    phiValue,
                    alphaValue
                ] / aValue,

                toCRT[
                    G[0, 1],
                    etaValue,
                    phiValue,
                    alphaValue
                ] / aValue
            };


        (*
            Generate a small neighborhood of reciprocal lattice nodes.
            A range of +/-2 is more than sufficient to identify the
            six nearest neighbors of a hexagonal reciprocal lattice.
        *)
        candidateNodes =
            DeleteCases[
                Flatten[
                    Table[
                        {n1, n2},
                        {n1, -2, 2},
                        {n2, -2, 2}
                    ],
                    1
                ],
                {0, 0}
            ];


        candidateVectors =
            (
                toCRT[
                    G[
                        Sequence @@ #
                    ],
                    etaValue,
                    phiValue,
                    alphaValue
                ] / aValue
            ) & /@ candidateNodes;


        (* Six reciprocal vectors defining the hexagonal Wigner-Seitz cell. *)
        nearestVectors =
            Take[
                SortBy[
                    candidateVectors,
                    Norm
                ],
                6
            ];


        (* Order the edge normals counter-clockwise. *)
        orderedNeighbors =
            SortBy[
                nearestVectors,
                ArcTan[
                    #[[1]],
                    #[[2]]
                ] &
            ];


        neighborAngles =
            ArcTan[
                #[[1]],
                #[[2]]
            ] & /@ orderedNeighbors;


        (*
            Each BZ corner is the intersection of the two Bragg planes

                r . g_i     = |g_i|^2 / 2
                r . g_(i+1) = |g_(i+1)|^2 / 2.
        *)
        vertices =
            Table[
                Module[
                    {
                        g1,
                        g2,
                        vertex2D
                    },

                    g1 =
                        orderedNeighbors[[i]];

                    g2 =
                        orderedNeighbors[[
                            Mod[
                                i,
                                Length[orderedNeighbors]
                            ]
                            + 1
                        ]];


                    vertex2D =
                        Quiet[
                            Check[
                                LinearSolve[
                                    {
                                        g1[[1 ;; 2]],
                                        g2[[1 ;; 2]]
                                    },

                                    {
                                        g1 . g1 / 2,
                                        g2 . g2 / 2
                                    }
                                ],

                                $Failed
                            ]
                        ];


                    If[
                        vertex2D === $Failed,

                        $Failed,

                        {
                            vertex2D[[1]],
                            vertex2D[[2]],
                            0
                        }
                    ]
                ],

                {
                    i,
                    Length[orderedNeighbors]
                }
            ];


        If[
            MemberQ[
                vertices,
                $Failed
            ],

            Message[
                FirstBrillouinZoneData::geometry
            ];

            Return[$Failed]
        ];


        orderedVertices =
            SortBy[
                vertices,
                ArcTan[
                    #[[1]],
                    #[[2]]
                ] &
            ];


        (*
            Select the BZ edge whose normal points closest to +x.
            M is exactly half of that reciprocal neighbor vector.
        *)
        selectedNeighborIndex =
            First[
                Ordering[
                    Abs[
                        neighborAngles
                    ],
                    1
                ]
            ];


        selectedNeighbor =
            orderedNeighbors[[
                selectedNeighborIndex
            ]];


        mPoint =
            selectedNeighbor / 2;


        (*
            Locate the two BZ vertices lying on the selected edge.
            The counter-clockwise endpoint is used as K.
        *)
        edgeTolerance =
            10^-8
            Max[
                1,
                selectedNeighbor . selectedNeighbor
            ];


        adjacentVertices =
            Select[
                orderedVertices,

                Abs[
                    selectedNeighbor[[1 ;; 2]]
                    .
                    #[[1 ;; 2]]
                    -
                    selectedNeighbor . selectedNeighbor / 2
                ]
                <
                edgeTolerance
                &
            ];


        If[
            Length[
                adjacentVertices
            ]
            !=
            2,

            Message[
                FirstBrillouinZoneData::geometry
            ];

            Return[$Failed]
        ];


        kPoint =
            First[
                MaximalBy[
                    adjacentVertices,

                    (
                        selectedNeighbor[[1]]
                        #[[2]]
                        -
                        selectedNeighbor[[2]]
                        #[[1]]
                    ) &
                ]
            ];


        gammaPoint =
            {0, 0, 0};


        <|
            "Gamma" ->
                gammaPoint,

            "M" ->
                mPoint,

            "K" ->
                kPoint,

            "Vertices" ->
                orderedVertices,

            "ReciprocalBasis" ->
                reciprocalBasis,

            "NearestReciprocalVectors" ->
                orderedNeighbors,

            "SelectedEdgeNormal" ->
                selectedNeighbor
        |>
    ];

(* --------------------------------------------------------------- *)
(* Gamma-M-K-Gamma path                                           *)
(* --------------------------------------------------------------- *)

(*
    Construct a mesh along the standard high-symmetry path

        Gamma -> M -> K -> Gamma.

    nPoints is the total number of mesh points, including both Gamma
    endpoints and the M and K points.

    The available intervals are distributed approximately in
    proportion to the geometric length of each path segment, while
    guaranteeing at least one interval per segment.

    PathCoordinate is the cumulative Cartesian distance along the
    high-symmetry path and is used as the horizontal plotting axis.
*)
GammaMKPathData[
    bzData_Association,
    nPoints_Integer
] :=
    Module[
        {
            gammaPoint,
            mPoint,
            kPoint,

            pathVertices,
            segmentLengths,
            totalIntervals,
            remainingIntervals,

            rawAdditionalIntervals,
            additionalIntervals,
            leftoverIntervals,
            fractionalParts,
            intervalOrder,
            intervalsPerSegment,

            segmentStartCoordinates,
            pathEntries,
            segmentEntries,

            highSymmetryCoordinates
        },


        gammaPoint =
            bzData["Gamma"];

        mPoint =
            bzData["M"];

        kPoint =
            bzData["K"];


        pathVertices =
            {
                gammaPoint,
                mPoint,
                kPoint,
                gammaPoint
            };


        segmentLengths =
            Norm[
                pathVertices[[# + 1]]
                -
                pathVertices[[#]]
            ] & /@
                Range[3];


        totalIntervals =
            nPoints - 1;


        (*
            Reserve one interval for each segment first, then allocate
            the remaining intervals according to segment length.
        *)
        remainingIntervals =
            totalIntervals - 3;


        rawAdditionalIntervals =
            If[
                remainingIntervals > 0,

                remainingIntervals
                segmentLengths
                /
                Total[
                    segmentLengths
                ],

                ConstantArray[
                    0,
                    3
                ]
            ];


        additionalIntervals =
            Floor[
                rawAdditionalIntervals
            ];


        leftoverIntervals =
            remainingIntervals
            -
            Total[
                additionalIntervals
            ];


        fractionalParts =
            FractionalPart[
                rawAdditionalIntervals
            ];


        If[
            leftoverIntervals > 0,

            intervalOrder =
                Ordering[
                    fractionalParts,
                    -leftoverIntervals
                ];

            Do[
                additionalIntervals[[
                    intervalIndex
                ]] += 1,

                {
                    intervalIndex,
                    intervalOrder
                }
            ];
        ];


        intervalsPerSegment =
            1
            +
            additionalIntervals;


        segmentStartCoordinates =
            {
                0,
                segmentLengths[[1]],
                segmentLengths[[1]]
                    +
                    segmentLengths[[2]]
            };


        pathEntries =
            {};


        Do[
            segmentEntries =
                Table[
                    <|
                        "WaveVector" ->
                            (
                                pathVertices[[segmentIndex]]
                                +
                                t
                                (
                                    pathVertices[[segmentIndex + 1]]
                                    -
                                    pathVertices[[segmentIndex]]
                                )
                            ),

                        "PathCoordinate" ->
                            (
                                segmentStartCoordinates[[
                                    segmentIndex
                                ]]
                                +
                                t
                                segmentLengths[[
                                    segmentIndex
                                ]]
                            ),

                        "Segment" ->
                            segmentIndex
                    |>,

                    {
                        t,
                        N[
                            Subdivide[
                                0,
                                1,
                                intervalsPerSegment[[
                                    segmentIndex
                                ]]
                            ]
                        ]
                    }
                ];


            (* Avoid duplicating M and K at adjacent segments. *)
            If[
                segmentIndex > 1,
                segmentEntries =
                    Rest[
                        segmentEntries
                    ];
            ];


            pathEntries =
                Join[
                    pathEntries,
                    segmentEntries
                ],

            {
                segmentIndex,
                1,
                3
            }
        ];


        highSymmetryCoordinates =
            {
                0,
                segmentLengths[[1]],
                segmentLengths[[1]]
                    +
                    segmentLengths[[2]],
                Total[
                    segmentLengths
                ]
            };


        <|
            "Entries" ->
                pathEntries,

            "WaveVectors" ->
                Lookup[
                    pathEntries,
                    "WaveVector"
                ],

            "PathCoordinates" ->
                Lookup[
                    pathEntries,
                    "PathCoordinate"
                ],

            "IntervalsPerSegment" ->
                intervalsPerSegment,

            "SegmentLengths" ->
                segmentLengths,

            "HighSymmetryCoordinates" ->
                highSymmetryCoordinates,

            "HighSymmetryLabels" ->
                {
                    "\[CapitalGamma]",
                    "M",
                    "K",
                    "\[CapitalGamma]"
                }
        |>
    ];


(* --------------------------------------------------------------- *)
(* Transverse Gamma-M-K-Gamma magnon spectrum                     *)
(* --------------------------------------------------------------- *)

Options[TransverseMagnonSpectrum] = {
    ImageSize -> 800,
    "ShowLegend" -> True,
    "SymmetryTolerance" -> 10^-6
};

(*
    Calculate the in-plane magnon dispersion along

        Gamma -> M -> K -> Gamma

    in the first magnetic Brillouin zone.

    All path-point eigensystems are evaluated independently with
    ParallelTable.

    The first Gamma point fixes:
        - mode number;
        - color;
        - angular-momentum label m.

    Mode identity is then propagated sequentially along the path by
    comparing every mode only with the modes selected at the previous
    neighboring path point.

    The eigensystem is solved internally up to 1.15 omegaMax.
    The plot itself is restricted to omegaMax.
*)
TransverseMagnonSpectrum[
    dynamicsData_Association,
    omegaMax_: 3,
    nPoints_: 401,
    OptionsPattern[]
] :=
    Module[
        {
            internalOmegaMax,

            imageSize,
            showLegend,
            symmetryTolerance,

            bzData,
            pathData,
            pathEntries,
            pathCoordinates,
            highSymmetryCoordinates,
            highSymmetryLabels,

            parallelDynamicsData,
            rawPathData,
            rawModes,

            gammaModes,
            gammaFrequencies,
            referenceModeIndices,

            angularAnalysis,
            angularMomentumIndices,

            trackedData,

            branches,
            branchFrequencies,

            modeColors,
            modeLabels,

            highSymmetryTicks,
            plot
        },


        (* Validate inputs. *)
        If[
            ! NumericQ[
                omegaMax
            ]
            ||
            omegaMax <= 0,

            Message[
                TransverseMagnonSpectrum::omax,
                omegaMax
            ];

            Return[$Failed]
        ];


        If[
            ! IntegerQ[
                nPoints
            ]
            ||
            nPoints < 4,

            Message[
                TransverseMagnonSpectrum::npoints,
                nPoints
            ];

            Return[$Failed]
        ];


        internalOmegaMax =
            1.15
            omegaMax;


        imageSize =
            OptionValue[
                ImageSize
            ];

        showLegend =
            OptionValue[
                "ShowLegend"
            ];

        symmetryTolerance =
            OptionValue[
                "SymmetryTolerance"
            ];


        (* Reconstruct BZ and build Gamma-M-K-Gamma path. *)
        bzData =
            FirstBrillouinZoneData[
                dynamicsData
            ];


        If[
            bzData === $Failed,
            Return[$Failed]
        ];


        pathData =
            GammaMKPathData[
                bzData,
                nPoints
            ];


        pathEntries =
            pathData[
                "Entries"
            ];

        pathCoordinates =
            pathData[
                "PathCoordinates"
            ];

        highSymmetryCoordinates =
            pathData[
                "HighSymmetryCoordinates"
            ];

        highSymmetryLabels =
            pathData[
                "HighSymmetryLabels"
            ];


        (* Data required by SolveMagneticModes on parallel kernels. *)
        parallelDynamicsData =
            KeyTake[
                dynamicsData,
                {
                    "Equilibrium",
                    "Nodes",
                    "TransverseData",
                    "ConstraintHessianTransverse"
                }
            ];


        DistributeDefinitions[
            parallelDynamicsData,
            SolveMagneticModes
        ];


        (* Every path point, including both Gamma endpoints, is solved in parallel. *)
        rawPathData =
            ParallelTable[
                <|
                    "PathCoordinate" ->
                        pathEntry[
                            "PathCoordinate"
                        ],

                    "WaveVector" ->
                        pathEntry[
                            "WaveVector"
                        ],

                    "Segment" ->
                        pathEntry[
                            "Segment"
                        ],

                    "Modes" ->
                        SolveMagneticModes[
                            parallelDynamicsData,

                            pathEntry[
                                "WaveVector"
                            ],

                            "FrequencyRange" ->
                                {
                                    0,
                                    internalOmegaMax
                                }
                        ]
                |>,

                {
                    pathEntry,
                    pathEntries
                }
            ];


        If[
            MemberQ[
                Lookup[
                    rawPathData,
                    "Modes"
                ],
                $Failed
            ],

            Message[
                TransverseMagnonSpectrum::parallel
            ];

            Return[$Failed]
        ];


        rawModes =
            Lookup[
                rawPathData,
                "Modes"
            ];


        (* The first path point is Gamma and defines the mode identities. *)
        gammaModes =
            First[
                rawModes
            ];

        gammaFrequencies =
            gammaModes[
                "Frequencies"
            ];


        referenceModeIndices =
            Select[
                Range[
                    Length[
                        gammaFrequencies
                    ]
                ],

                gammaFrequencies[[#]]
                    <=
                omegaMax
                &
            ];


        If[
            referenceModeIndices === {},

            Message[
                TransverseMagnonSpectrum::nomodes,
                omegaMax
            ];

            Return[$Failed]
        ];


        angularAnalysis =
            AnalyzeAngularMomentumModes[
                dynamicsData,
                gammaModes,

                "SymmetryTolerance" ->
                    symmetryTolerance
            ];


        If[
            angularAnalysis === $Failed,
            Return[$Failed]
        ];


        angularMomentumIndices =
            Lookup[
                angularAnalysis[
                    "Modes"
                ],
                "AngularMomentumClassification"
            ];


        (*
            Sequential identification pass over already-computed data.
            This is inexpensive: only BdG overlap matrices are built.
        *)
        trackedData =
            PropagateModeIdentification[
                dynamicsData,
                rawModes,
                referenceModeIndices
            ];


        branches =
            Table[
                branchFrequencies =
                    (
                        #[
                            "Frequencies"
                        ][[
                            branchIndex
                        ]]
                        &
                        /@
                        trackedData
                    );


                Cases[
                    Transpose[
                        {
                            pathCoordinates,
                            branchFrequencies
                        }
                    ],

                    {
                        x_?NumericQ,
                        y_?NumericQ
                    }
                ],

                {
                    branchIndex,
                    Length[
                        referenceModeIndices
                    ]
                }
            ];


        modeColors =
            DispersionModeColor
                /@
            referenceModeIndices;


        modeLabels =
            Table[
                With[
                    {
                        referenceIndex =
                            referenceModeIndices[[
                                branchIndex
                            ]]
                    },

                    Row[
                        {
                            "Mode ",
                            referenceIndex,
                            "   (m = ",
                            If[
                                MissingQ[
                                    angularMomentumIndices[[
                                        referenceIndex
                                    ]]
                                ],

                                "?",

                                angularMomentumIndices[[
                                    referenceIndex
                                ]]
                            ],
                            ")"
                        }
                    ]
                ],

                {
                    branchIndex,
                    Length[
                        referenceModeIndices
                    ]
                }
            ];


        highSymmetryTicks =
            MapThread[
                {
                    #1,
                    Style[
                        #2,
                        14
                    ]
                }
                &,

                {
                    highSymmetryCoordinates,
                    highSymmetryLabels
                }
            ];


        plot =
            ListLinePlot[
                branches,

                Joined ->
                    True,

                PlotStyle ->
                    (
                        Directive[
                            #,
                            Thick
                        ]
                        &
                        /@
                        modeColors
                    ),

                PlotRange ->
                    {
                        {
                            First[
                                highSymmetryCoordinates
                            ],

                            Last[
                                highSymmetryCoordinates
                            ]
                        },

                        {
                            0,
                            omegaMax
                        }
                    },

                Frame ->
                    True,

                Axes ->
                    False,

                FrameTicks ->
                    {
                        {
                            Automatic,
                            None
                        },

                        {
                            highSymmetryTicks,
                            None
                        }
                    },

                GridLines ->
                    {
                        highSymmetryCoordinates,
                        None
                    },

                GridLinesStyle ->
                    Directive[
                        GrayLevel[0.8],
                        Dashed
                    ],

                FrameLabel ->
                    {
                        None,

                        Style[
                            Row[
                                {
                                    "\[Omega]/",
                                    Subscript[
                                        "\[Omega]",
                                        "c2"
                                    ]
                                }
                            ],
                            15
                        ]
                    },

                LabelStyle ->
                    13,

                PlotLegends ->
                    If[
                        TrueQ[
                            showLegend
                        ],

                        Placed[
                            LineLegend[
                                modeColors,
                                modeLabels,
                                LegendLayout ->
                                    "Column"
                            ],
                            Right
                        ],

                        None
                    ],

                ImageSize ->
                    imageSize
            ];


        <|
            "OmegaMax" ->
                omegaMax,

            "InternalOmegaMax" ->
                internalOmegaMax,

            "NumberOfPoints" ->
                nPoints,

            "BrillouinZone" ->
                bzData,

            "Path" ->
                pathData,

            "GammaModes" ->
                gammaModes,

            "ReferenceModeIndices" ->
                referenceModeIndices,

            "AngularMomentumAnalysis" ->
                angularAnalysis,

            "AngularMomentumIndices" ->
                angularMomentumIndices,

            "RawData" ->
                rawPathData,

            "TrackedData" ->
                trackedData,

            "Branches" ->
                branches,

            "ModeColors" ->
                modeColors,

            "ModeLabels" ->
                modeLabels,

            "Plot" ->
                plot
        |>
    ];

(* --------------------------------------------------------------- *)
(* Help                                                            *)
(* --------------------------------------------------------------- *)

DispersionHelp[] :=
    Column[
        {
            Style[
                "Dispersion.wl",
                Bold,
                16
            ],

            "Magnon dispersion calculations with parallel eigensystems and nearest-neighbor BdG mode tracking.",

            Grid[
                {
                    {"Function", "Description"},
                    {"DispersionModeColor", "Return the fixed six-color cyclic publication palette. Mode 7 reuses the color of mode 1."},
                    {"MagneticBdGOverlapMatrix", "Calculate squared magnetic BdG overlaps between two sets of transverse eigenvectors."},
                    {"MatchNeighboringModes", "Perform a one-to-one greedy matching between modes at neighboring wavevectors using the largest BdG overlaps."},
                    {"PropagateModeIdentification", "Propagate Gamma-point branch identities through an ordered sequence of neighboring eigensystems."},
                    {"FirstBrillouinZoneData", "Reconstruct the first hexagonal magnetic Brillouin zone from the reciprocal lattice and return Gamma, M, K and the BZ vertices."},
                    {"GammaMKPathData", "Construct the wavevector mesh and cumulative path coordinate along Gamma-M-K-Gamma."},
                    {"LongitudinalMagnonSpectrum", "Calculate the longitudinal spectrum in parallel and propagate Gamma-point identities independently toward positive and negative qz."},
                    {"TransverseMagnonSpectrum", "Calculate the Gamma-M-K-Gamma spectrum in parallel and propagate Gamma-point identities sequentially along the path."}
                },
                Frame -> All,
                Alignment -> Left
            ]
        }
    ];

Print[
    DispersionHelp[]
];