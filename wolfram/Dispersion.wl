(* ================================================================ *)
(* Dispersion.wl                                                    *)
(*                                                                  *)
(* Magnon dispersion calculations and mode tracking.               *)
(*                                                                  *)
(* Current implementation:                                         *)
(*     longitudinal dispersion for k = {0,0,qz}.                   *)
(*                                                                  *)
(* Mode identity is fixed at the Gamma point. Finite-qz modes are   *)
(* assigned directly to the Gamma-point modes using the largest     *)
(* magnetic BdG overlap.                                           *)
(*                                                                  *)
(* Dependencies:                                                    *)
(*     DynamicsMagnetic.wl                                         *)
(*     ModeAnalysis.wl                                             *)
(* ================================================================ *)

ClearAll[
    DispersionModeColor,
    AssignModesToGamma,
    LongitudinalMagnonSpectrum,
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
    "At least one finite-qz calculation failed. Inspect the returned raw data for the affected wavevectors.";

(* --------------------------------------------------------------- *)
(* Publication color cycle                                        *)
(* --------------------------------------------------------------- *)

(*
    Six-color publication-friendly cycle based on the Okabe-Ito
    colorblind-safe palette.

    Colors repeat every six modes, so mode 7 uses the same color as
    mode 1, mode 8 the same as mode 2, and so on.
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
(* BdG assignment to Gamma-point modes                            *)
(* --------------------------------------------------------------- *)

(*
    Assign all finite-wavevector eigenmodes to the reference modes
    calculated at the Gamma point.

    For each finite-q mode, calculate its squared magnetic BdG overlap
    with every Gamma-point mode,

        |<Gamma_j | mode_i(q)>_BdG|^2,

    and assign mode_i(q) to the Gamma mode with the largest overlap.

    More than one finite-q mode can in principle be assigned to the
    same Gamma reference. For the plotted branch, only the state with
    the largest overlap with that reference is retained.

    The full assignment information is returned for diagnostics.
*)
AssignModesToGamma[
    dynamicsData_Association,
    gammaModes_Association,
    finiteModes_Association
] :=
    Module[
        {
            transverseData,
            gammaVectors,
            finiteVectors,
            finiteFrequencies,

            numberOfGammaModes,
            numberOfFiniteModes,

            overlapMatrix,
            assignments,
            selectedByReference,

            bestReference,
            candidateAssignments
        },


        transverseData =
            dynamicsData["TransverseData"];

        gammaVectors =
            gammaModes["EigenvectorsTransverse"];

        finiteVectors =
            finiteModes["EigenvectorsTransverse"];

        finiteFrequencies =
            finiteModes["Frequencies"];

        numberOfGammaModes =
            Length[
                gammaVectors
            ];

        numberOfFiniteModes =
            Length[
                finiteVectors
            ];


        (* No finite-frequency modes in the requested window. *)
        If[
            numberOfFiniteModes == 0,

            Return[
                <|
                    "OverlapMatrix" ->
                        {},

                    "Assignments" ->
                        {},

                    "SelectedByReference" ->
                        ConstantArray[
                            Missing["NotAssigned"],
                            numberOfGammaModes
                        ]
                |>
            ]
        ];


        (*
            Rows correspond to finite-q modes.
            Columns correspond to Gamma-point reference modes.
        *)
        overlapMatrix =
            Table[
                Abs[
                    MagneticBdGInnerProduct[
                        transverseData,
                        gammaVectors[[referenceIndex]],
                        finiteVectors[[finiteIndex]]
                    ]
                ]^2,

                {
                    finiteIndex,
                    numberOfFiniteModes
                },

                {
                    referenceIndex,
                    numberOfGammaModes
                }
            ];


        (* Assign each finite-q mode to its strongest Gamma overlap. *)
        assignments =
            Table[
                bestReference =
                    First[
                        Ordering[
                            overlapMatrix[[finiteIndex]],
                            -1
                        ]
                    ];

                <|
                    "FiniteModeIndex" ->
                        finiteIndex,

                    "Frequency" ->
                        finiteFrequencies[[finiteIndex]],

                    "ReferenceModeIndex" ->
                        bestReference,

                    "Overlap" ->
                        overlapMatrix[[
                            finiteIndex,
                            bestReference
                        ]]
                |>,

                {
                    finiteIndex,
                    numberOfFiniteModes
                }
            ];


        (*
            Construct at most one plotted state for every Gamma
            reference mode.

            If several finite-q states were assigned to the same
            reference, retain the one with the largest overlap.
        *)
        selectedByReference =
            Table[
                candidateAssignments =
                    Select[
                        assignments,
                        #["ReferenceModeIndex"]
                            ==
                        referenceIndex
                        &
                    ];


                If[
                    candidateAssignments === {},

                    Missing["NotAssigned"],

                    First[
                        MaximalBy[
                            candidateAssignments,
                            #["Overlap"] &
                        ]
                    ]
                ],

                {
                    referenceIndex,
                    numberOfGammaModes
                }
            ];


        <|
            "OverlapMatrix" ->
                overlapMatrix,

            "Assignments" ->
                assignments,

            "SelectedByReference" ->
                selectedByReference
        |>
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
    Calculate the longitudinal magnon spectrum for

        k = {0,0,qz}.

    Arguments:
        dynamicsData
            Prepared magnetic dynamics from PrepareMagneticDynamics.

        qzMax
            Maximum absolute longitudinal wavevector.
            Default: 1.

        omegaMax
            Maximum frequency included in the spectrum.
            Default: 3.

        nPoints
            Number of points in the uniform mesh.
            Default: 201.

    The mesh is always symmetric around qz = 0 and always contains
    qz = 0. If an even nPoints is supplied, it is increased by one.

    Workflow:
        1. Solve the Gamma-point eigensystem.
        2. Determine the C6 angular-momentum index m of every
           Gamma-point mode.
        3. Solve all qz != 0 points in parallel using ParallelTable.
        4. Assign finite-qz states directly to the Gamma modes using
           the largest magnetic BdG overlap.
        5. Plot every Gamma-defined branch with a fixed cyclic color.

    The function returns an Association containing both numerical data
    and the final plot.
*)
LongitudinalMagnonSpectrum[
    dynamicsData_Association,
    qzMax_: 1,
    omegaMax_: 3,
    nPoints_: 201,
    OptionsPattern[]
] :=
    Module[
        {
            actualNPoints,
            qzValues,
            zeroIndex,
            nonzeroQzValues,

            imageSize,
            showLegend,
            symmetryTolerance,

            gammaModes,
            angularAnalysis,
            angularMomentumIndices,
            numberOfReferenceModes,

            parallelDynamicsData,
            rawFiniteQData,
            assignedFiniteQData,

            branches,
            selectedState,
            finitePoints,

            modeColors,
            modeLabels,
            plot
        },


        (* ------------------------------------------------------- *)
        (* Validate inputs                                         *)
        (* ------------------------------------------------------- *)

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


        imageSize =
            OptionValue[ImageSize];

        showLegend =
            OptionValue["ShowLegend"];

        symmetryTolerance =
            OptionValue["SymmetryTolerance"];


        (* ------------------------------------------------------- *)
        (* Uniform symmetric qz mesh                              *)
        (* ------------------------------------------------------- *)

        qzValues =
            N[
                Subdivide[
                    -qzMax,
                    qzMax,
                    actualNPoints - 1
                ]
            ];


        (* Explicitly set the central mesh point to exact zero. *)
        zeroIndex =
            First[
                Ordering[
                    Abs[qzValues],
                    1
                ]
            ];

        qzValues[[
            zeroIndex
        ]] =
            0;


        nonzeroQzValues =
            Delete[
                qzValues,
                zeroIndex
            ];


        (* ------------------------------------------------------- *)
        (* Gamma-point reference modes                            *)
        (* ------------------------------------------------------- *)

        gammaModes =
            SolveMagneticModes[
                dynamicsData,
                {0, 0, 0},
                "FrequencyRange" ->
                    {
                        0,
                        omegaMax
                    }
            ];


        numberOfReferenceModes =
            gammaModes["NumberOfModes"];


        If[
            numberOfReferenceModes == 0,

            Message[
                LongitudinalMagnonSpectrum::nomodes,
                omegaMax
            ];

            Return[$Failed]
        ];


        (* Determine m for every Gamma-point reference mode. *)
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
                angularAnalysis["Modes"],
                "AngularMomentumClassification"
            ];


        (* ------------------------------------------------------- *)
        (* Parallel finite-qz eigensystems                        *)
        (* ------------------------------------------------------- *)

        (*
            Only the fields required by SolveMagneticModes are sent
            to the parallel kernels. This avoids distributing the
            large full constraint Hessian and full cross operator.
        *)
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


        (*
            Distribute the reusable numerical dynamics object once.
            SolveMagneticModes and its dependent Global` definitions
            are then available to the parallel kernels.
        *)
        DistributeDefinitions[
            parallelDynamicsData,
            SolveMagneticModes
        ];


        rawFiniteQData =
            ParallelTable[
                Module[
                    {
                        finiteModes
                    },

                    finiteModes =
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
                                    omegaMax
                                }
                        ];


                    <|
                        "Qz" ->
                            qz,

                        "Modes" ->
                            finiteModes
                    |>
                ],

                {
                    qz,
                    nonzeroQzValues
                }
            ];


        If[
            MemberQ[
                Lookup[
                    rawFiniteQData,
                    "Modes"
                ],
                $Failed
            ],

            Message[
                LongitudinalMagnonSpectrum::parallel
            ];
        ];


        (* ------------------------------------------------------- *)
        (* BdG assignment to Gamma-point references               *)
        (* ------------------------------------------------------- *)

        assignedFiniteQData =
            Map[
                Function[
                    qData,

                    If[
                        qData["Modes"] === $Failed,

                        <|
                            "Qz" ->
                                qData["Qz"],

                            "Modes" ->
                                $Failed,

                            "Assignment" ->
                                $Failed
                        |>,

                        <|
                            "Qz" ->
                                qData["Qz"],

                            "Modes" ->
                                qData["Modes"],

                            "Assignment" ->
                                AssignModesToGamma[
                                    dynamicsData,
                                    gammaModes,
                                    qData["Modes"]
                                ]
                        |>
                    ]
                ],

                rawFiniteQData
            ];


        (* ------------------------------------------------------- *)
        (* Construct one branch per Gamma-point mode              *)
        (* ------------------------------------------------------- *)

        branches =
            Table[

                finitePoints =
                    DeleteMissing[
                        Table[
                            If[
                                qData["Assignment"] === $Failed,

                                Missing[
                                    "CalculationFailed"
                                ],

                                selectedState =
                                    qData[
                                        "Assignment"
                                    ][
                                        "SelectedByReference"
                                    ][[
                                        referenceIndex
                                    ]];


                                If[
                                    MissingQ[
                                        selectedState
                                    ],

                                    Missing[
                                        "NotAssigned"
                                    ],

                                    {
                                        qData["Qz"],
                                        selectedState["Frequency"]
                                    }
                                ]
                            ],

                            {
                                qData,
                                assignedFiniteQData
                            }
                        ]
                    ];


                SortBy[
                    Join[
                        finitePoints,

                        {
                            {
                                0,
                                gammaModes[
                                    "Frequencies"
                                ][[
                                    referenceIndex
                                ]]
                            }
                        }
                    ],

                    First
                ],

                {
                    referenceIndex,
                    numberOfReferenceModes
                }
            ];


        (* ------------------------------------------------------- *)
        (* Fixed mode colors and labels                           *)
        (* ------------------------------------------------------- *)

        modeColors =
            Table[
                DispersionModeColor[
                    referenceIndex
                ],

                {
                    referenceIndex,
                    numberOfReferenceModes
                }
            ];


        modeLabels =
            Table[
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
                ],

                {
                    referenceIndex,
                    numberOfReferenceModes
                }
            ];


        (* ------------------------------------------------------- *)
        (* Spectrum plot                                           *)
        (* ------------------------------------------------------- *)

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
                            "q_z",
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


        (* ------------------------------------------------------- *)
        (* Return data and plot                                   *)
        (* ------------------------------------------------------- *)

        <|
            "QzMax" ->
                qzMax,

            "OmegaMax" ->
                omegaMax,

            "RequestedNumberOfPoints" ->
                nPoints,

            "NumberOfPoints" ->
                actualNPoints,

            "QzValues" ->
                qzValues,

            "GammaModes" ->
                gammaModes,

            "AngularMomentumAnalysis" ->
                angularAnalysis,

            "AngularMomentumIndices" ->
                angularMomentumIndices,

            "RawFiniteQData" ->
                rawFiniteQData,

            "AssignedFiniteQData" ->
                assignedFiniteQData,

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

            "Longitudinal magnon dispersion calculations with Gamma-point mode identification and BdG-overlap tracking.",

            Grid[
                {
                    {"Function", "Description"},
                    {"DispersionModeColor", "Return the fixed six-color cyclic publication palette for a mode index. Mode 7 reuses the color of mode 1."},
                    {"AssignModesToGamma", "Assign finite-wavevector magnon modes to Gamma-point reference modes using the largest magnetic BdG overlap."},
                    {"LongitudinalMagnonSpectrum", "Calculate the longitudinal spectrum for k = {0,0,qz}, classify Gamma-point modes by angular momentum, evaluate all nonzero qz points in parallel, and return tracked branches plus the final plot."}
                },
                Frame -> All,
                Alignment -> Left
            ]
        }
    ];

Print[
    DispersionHelp[]
];