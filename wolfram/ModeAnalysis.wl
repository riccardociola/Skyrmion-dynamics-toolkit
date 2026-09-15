(* ================================================================ *)
(* ModeAnalysis.wl                                                  *)
(*                                                                  *)
(* Analysis and classification of magnetic magnon eigenmodes.      *)
(*                                                                  *)
(* Implemented analyses at k = 0:                                 *)
(*     - dipole-moment classification;                              *)
(*     - C6 angular-momentum classification from rotation symmetry. *)
(*                                                                  *)
(* Further mode classifiers can be added without changing the       *)
(* existing dipole or angular-momentum interfaces.                  *)
(*                                                                  *)
(* Dependencies:                                                    *)
(*     Geometry.wl                                                  *)
(*     FourierBasis.wl                                              *)
(*     DynamicsProjection.wl                                       *)
(*     DynamicsMagnetic.wl                                         *)
(* ================================================================ *)

ClearAll[
    ModeDipoleMoment,
    ClassifyDipoleMode,
    AnalyzeDipoleModes,

    ModeAngularMomentum,
    AnalyzeAngularMomentumModes,

    ModeAnalysisHelp
];

(* --------------------------------------------------------------- *)
(* Messages                                                        *)
(* --------------------------------------------------------------- *)

ModeDipoleMoment::nonzeroK =
    "Dipole classification is defined here only for k = 0. The supplied mode data use k = `1`.";

ModeDipoleMoment::mode =
    "Mode index `1` is invalid. The supplied eigensystem contains `2` modes.";

ModeDipoleMoment::zeronode =
    "The reciprocal-space node {0,0} was not found in the dynamical basis.";

ModeDipoleMoment::convention =
    "Unknown time convention `1`. Use \"Exp[-Iomega t]\" or \"Exp[+Iomega t]\".";

AnalyzeDipoleModes::convention =
    "Unknown time convention `1`. Use \"Exp[-Iomega t]\" or \"Exp[+Iomega t]\".";

ModeAngularMomentum::nonzeroK =
    "Angular-momentum classification is defined here only for k = 0. The supplied mode data use k = `1`.";

ModeAngularMomentum::mode =
    "Mode index `1` is invalid. The supplied eigensystem contains `2` modes.";

ModeAngularMomentum::samples =
    "At least one valid real-space sample point is required.";

ModeAngularMomentum::zeronorm =
    "The sampled eigenfunction has vanishing norm, so its C6 symmetry cannot be classified.";

(* --------------------------------------------------------------- *)
(* Dipole moment of one mode                                      *)
(* --------------------------------------------------------------- *)

Options[ModeDipoleMoment] = {
    "WaveVectorTolerance" -> 10^-10,
    "TimeConvention" -> "Exp[-Iomega t]"
};

(*
    Extract the uniform magnetic dipole component of one normalized
    magnon eigenmode.

    This analysis is restricted to k = 0.

    At k = 0,

        integral_uc delta m(r) d^2r

    is proportional to the Fourier coefficient at reciprocal node
    G = 0. Therefore no real-space integration is required.

    The eigenvectors returned by SolveMagneticModes are already
    BdG-normalized and stored in Cartesian components.

    Circular amplitudes are defined as

        dPlus  = (dx + i dy)/Sqrt[2]
        dMinus = (dx - i dy)/Sqrt[2].

    For the default time convention

        delta m(t) = Re[delta m Exp[-i omega t]],

    dMinus corresponds to CCW rotation and dPlus to CW rotation.
*)
ModeDipoleMoment[
    dynamicsData_Association,
    modes_Association,
    modeIndex_Integer,
    OptionsPattern[]
] :=
    Module[
        {
            waveVectorTolerance,
            timeConvention,

            k,
            nodes,
            numberOfModes,
            zeroNodePosition,

            modeVector,
            modeCoefficients,
            dipoleVector,

            dx,
            dy,
            dz,

            circularPlus,
            circularMinus,

            cwAmplitude,
            ccwAmplitude
        },


        (* Options. *)
        waveVectorTolerance =
            OptionValue["WaveVectorTolerance"];

        timeConvention =
            OptionValue["TimeConvention"];


        (* Restrict the analysis to the Gamma point. *)
        k =
            modes["WaveVector"];

        If[
            Norm[N[k]] > waveVectorTolerance,

            Message[
                ModeDipoleMoment::nonzeroK,
                k
            ];

            Return[$Failed]
        ];


        (* Validate the requested mode index. *)
        numberOfModes =
            Length[
                modes["Frequencies"]
            ];

        If[
            modeIndex < 1
            ||
            modeIndex > numberOfModes,

            Message[
                ModeDipoleMoment::mode,
                modeIndex,
                numberOfModes
            ];

            Return[$Failed]
        ];


        (* Locate the G = 0 reciprocal-space node. *)
        nodes =
            dynamicsData["Nodes"];

        zeroNodePosition =
            FirstPosition[
                nodes,
                {0, 0},
                Missing["NotFound"]
            ];

        If[
            MissingQ[
                zeroNodePosition
            ],

            Message[
                ModeDipoleMoment::zeronode
            ];

            Return[$Failed]
        ];

        zeroNodePosition =
            First[
                zeroNodePosition
            ];


        (* Extract the Cartesian G = 0 Fourier coefficient. *)
        modeVector =
            modes["EigenvectorsFull"][[modeIndex]];

        modeCoefficients =
            Partition[
                modeVector,
                3
            ];

        dipoleVector =
            modeCoefficients[[
                zeroNodePosition
            ]];


        {
            dx,
            dy,
            dz
        } =
            dipoleVector;


        (* Circular in-plane components. *)
        circularPlus =
            (
                dx
                +
                I dy
            ) / Sqrt[2];

        circularMinus =
            (
                dx
                -
                I dy
            ) / Sqrt[2];


        (* Map circular amplitudes onto CW / CCW labels. *)
        {
            cwAmplitude,
            ccwAmplitude
        } =
            Switch[
                timeConvention,

                "Exp[-Iomega t]",
                {
                    circularPlus,
                    circularMinus
                },

                "Exp[+Iomega t]",
                {
                    circularMinus,
                    circularPlus
                },

                _,

                Message[
                    ModeDipoleMoment::convention,
                    timeConvention
                ];

                Return[$Failed]
            ];


        (* Return the dipole analysis of this mode. *)
        <|
            "DipoleVector" ->
                dipoleVector,

            "Dx" ->
                dx,

            "Dy" ->
                dy,

            "Dz" ->
                dz,

            "CircularPlusAmplitude" ->
                circularPlus,

            "CircularMinusAmplitude" ->
                circularMinus,

            "CWAmplitude" ->
                cwAmplitude,

            "CCWAmplitude" ->
                ccwAmplitude,

            "BreathingWeight" ->
                Abs[dz]^2,

            "CWWeight" ->
                Abs[cwAmplitude]^2,

            "CCWWeight" ->
                Abs[ccwAmplitude]^2,

            "TotalDipoleWeight" ->
                Total[
                    Abs[dipoleVector]^2
                ],

            "TimeConvention" ->
                timeConvention
        |>
    ];

(* --------------------------------------------------------------- *)
(* Classification of one dipole analysis                          *)
(* --------------------------------------------------------------- *)

(*
    Classify one mode from its dipole weights.

    activityThreshold is an absolute threshold for the total dipole
    weight. Modes below it are labeled "Dipole Inactive".

    For active modes, the largest of

        BreathingWeight,
        CWWeight,
        CCWWeight

    determines the dipole classification.
*)
ClassifyDipoleMode[
    dipoleAnalysis_Association,
    activityThreshold_?NumericQ
] :=
    Module[
        {
            totalWeight,
            channelWeights,
            dominantChannel
        },

        totalWeight =
            dipoleAnalysis[
                "TotalDipoleWeight"
            ];


        If[
            totalWeight <= activityThreshold,

            Return[
                "Dipole Inactive"
            ]
        ];


        channelWeights =
            {
                dipoleAnalysis["BreathingWeight"],
                dipoleAnalysis["CWWeight"],
                dipoleAnalysis["CCWWeight"]
            };


        dominantChannel =
            First[
                Ordering[
                    channelWeights,
                    -1
                ]
            ];


        Switch[
            dominantChannel,

            1,
            "Breathing",

            2,
            "CW",

            3,
            "CCW"
        ]
    ];

(* --------------------------------------------------------------- *)
(* Analyze all k = 0 modes                                        *)
(* --------------------------------------------------------------- *)

Options[AnalyzeDipoleModes] = {
    "RelativeDipoleThreshold" -> 10^-6,
    "AbsoluteDipoleThreshold" -> 0,
    "WaveVectorTolerance" -> 10^-10,
    "TimeConvention" -> "Exp[-Iomega t]"
};

(*
    Analyze and classify all modes in a k = 0 eigensystem.

    The activity threshold is chosen as

        max[
            AbsoluteDipoleThreshold,
            RelativeDipoleThreshold * maximum total dipole weight
        ].

    This avoids treating tiny numerical residuals in otherwise dark
    modes as physically meaningful dipole moments.

    The returned structure is designed to be extendable. Future
    angular-momentum information can be added alongside the existing
    "Dipole" field for every mode.
*)
AnalyzeDipoleModes[
    dynamicsData_Association,
    modes_Association,
    OptionsPattern[]
] :=
    Module[
        {
            relativeThreshold,
            absoluteThreshold,
            waveVectorTolerance,
            timeConvention,

            numberOfModes,
            dipoleData,

            maximumDipoleWeight,
            activityThreshold,

            modeResults
        },


        (* Options. *)
        relativeThreshold =
            OptionValue[
                "RelativeDipoleThreshold"
            ];

        absoluteThreshold =
            OptionValue[
                "AbsoluteDipoleThreshold"
            ];

        waveVectorTolerance =
            OptionValue[
                "WaveVectorTolerance"
            ];

        timeConvention =
            OptionValue[
                "TimeConvention"
            ];


        If[
            ! MemberQ[
                {
                    "Exp[-Iomega t]",
                    "Exp[+Iomega t]"
                },
                timeConvention
            ],

            Message[
                AnalyzeDipoleModes::convention,
                timeConvention
            ];

            Return[$Failed]
        ];


        numberOfModes =
            Length[
                modes["Frequencies"]
            ];


        (* Extract dipole information for every mode. *)
        dipoleData =
            Table[
                ModeDipoleMoment[
                    dynamicsData,
                    modes,
                    modeIndex,

                    "WaveVectorTolerance" ->
                        waveVectorTolerance,

                    "TimeConvention" ->
                        timeConvention
                ],

                {
                    modeIndex,
                    1,
                    numberOfModes
                }
            ];


        If[
            MemberQ[
                dipoleData,
                $Failed
            ],

            Return[$Failed]
        ];


        (* Determine the threshold for dipole activity. *)
        maximumDipoleWeight =
            Max[
                Lookup[
                    dipoleData,
                    "TotalDipoleWeight"
                ]
            ];


        activityThreshold =
            Max[
                absoluteThreshold,
                relativeThreshold
                maximumDipoleWeight
            ];


        (* Build the mode-by-mode analysis. *)
        modeResults =
            Table[
                <|
                    "ModeIndex" ->
                        modeIndex,

                    "Frequency" ->
                        modes["Frequencies"][[modeIndex]],

                    "Dipole" ->
                        dipoleData[[modeIndex]],

                    "DipoleClassification" ->
                        ClassifyDipoleMode[
                            dipoleData[[modeIndex]],
                            activityThreshold
                        ]
                |>,

                {
                    modeIndex,
                    1,
                    numberOfModes
                }
            ];


        (* Return the complete dipole analysis. *)
        <|
            "WaveVector" ->
                modes["WaveVector"],

            "TimeConvention" ->
                timeConvention,

            "RelativeDipoleThreshold" ->
                relativeThreshold,

            "AbsoluteDipoleThreshold" ->
                absoluteThreshold,

            "MaximumDipoleWeight" ->
                maximumDipoleWeight,

            "DipoleActivityThreshold" ->
                activityThreshold,

            "Modes" ->
                modeResults
        |>
    ];



(* --------------------------------------------------------------- *)
(* C6 angular-momentum classification                              *)
(* --------------------------------------------------------------- *)

Options[ModeAngularMomentum] = {
    "WaveVectorTolerance" -> 10^-10,

    (* Maximum accepted normalized mismatch for a C6 eigenstate. *)
    "SymmetryTolerance" -> 10^-6,

    (*
        Real-space points used to test the rotation eigenvalue.
        Automatic uses several generic points around the SkL center.
        User-supplied points must be Cartesian {x,y} coordinates
        expressed in units of a0.
    *)
    "SamplePoints" -> Automatic,

    (*
        Rotation center in Cartesian coordinates measured in units
        of a0. Automatic uses the center of the magnetic unit cell.
    *)
    "RotationCenter" -> Automatic
};

(*
    Determine the angular-momentum quantum number m modulo 6 of one
    magnetic eigenmode at k = 0.

    For a hexagonal skyrmion lattice the mode is tested against

        R6 deltaM_m(R6^-1 r)
            = Exp[i 2 Pi m/6] deltaM_m(r),

    with m = -2,-1,0,1,2,3.

    R6 acts on both:
        1. real-space coordinates, through R6^-1 r;
        2. the three Cartesian spin components, through R6.

    The rotation is performed about the center of the magnetic unit
    cell. Explicitly, for a point r,

        rPreimage =
            center + R6^-1 . (r - center).

    This is the automatic equivalent of the double coordinate shift
    used in the legacy notebook.

    The symmetry is evaluated at several real-space points. For each
    candidate m the normalized RMS mismatch is

        sqrt[
            Sum_p ||R6 deltaM(R6^-1 r_p)
                    - Exp[i 2 Pi m/6] deltaM(r_p)||^2
            /
            Sum_p ||deltaM(r_p)||^2
        ].

    The m with the smallest mismatch is returned. The classification
    is accepted only if that residual is below "SymmetryTolerance".

    Notes:
        - m is defined modulo 6. We use the representatives
          -2,-1,0,1,2,3 throughout the classification.
        - The equivalent non-negative modulo-6 value is also returned
          for backward compatibility.
        - Near an exact degeneracy, a numerical eigensolver can return
          linear combinations of symmetry eigenstates. Such a mode can
          therefore remain "Unclassified" even when the degenerate
          subspace itself has a well-defined C6 representation.
*)
ModeAngularMomentum[
    dynamicsData_Association,
    modes_Association,
    modeIndex_Integer,
    OptionsPattern[]
] :=
    Module[
        {
            waveVectorTolerance,
            symmetryTolerance,
            samplePointsOption,
            rotationCenterOption,

            k,
            equilibrium,
            nodes,
            numberOfModes,
            solutionRules,

            etaValue,
            phiValue,
            alphaValue,

            center,
            samplePoints,

            modeVector,
            modeCoefficients,
            reciprocalWaveVectors,

            rotationAngle,
            rotation2DInverse,
            rotation3D,

            modeField,
            originalFields,
            rotatedFields,

            denominator,
            mValues,
            residuals,
            bestPosition,
            bestM,
            bestResidual,

            acceptedM,
            signedM,
            classification
        },


        (* Options. *)
        waveVectorTolerance =
            OptionValue[
                "WaveVectorTolerance"
            ];

        symmetryTolerance =
            OptionValue[
                "SymmetryTolerance"
            ];

        samplePointsOption =
            OptionValue[
                "SamplePoints"
            ];

        rotationCenterOption =
            OptionValue[
                "RotationCenter"
            ];


        (* Restrict the symmetry classification to the Gamma point. *)
        k =
            modes["WaveVector"];

        If[
            Norm[N[k]] > waveVectorTolerance,

            Message[
                ModeAngularMomentum::nonzeroK,
                k
            ];

            Return[$Failed]
        ];


        (* Validate the requested mode index. *)
        numberOfModes =
            Length[
                modes["Frequencies"]
            ];

        If[
            modeIndex < 1
            ||
            modeIndex > numberOfModes,

            Message[
                ModeAngularMomentum::mode,
                modeIndex,
                numberOfModes
            ];

            Return[$Failed]
        ];


        equilibrium =
            dynamicsData["Equilibrium"];

        nodes =
            dynamicsData["Nodes"];

        solutionRules =
            EquilibriumRules[
                equilibrium
            ];

        etaValue =
            eta0 /. solutionRules;

        phiValue =
            phi0 /. solutionRules;

        alphaValue =
            alphaR0 /. solutionRules;


        (* Center of the C6 rotation. *)
        center =
            If[
                rotationCenterOption === Automatic,

                {
                    x0[
                        etaValue,
                        phiValue,
                        alphaValue
                    ],

                    y0[
                        etaValue,
                        phiValue,
                        alphaValue
                    ]
                },

                rotationCenterOption
            ];


        (*
            Default sample points are deliberately generic and not
            related by C6 symmetry. They are defined relative to the
            magnetic unit-cell center in dimensionless Cartesian
            coordinates.
        *)
        samplePoints =
            If[
                samplePointsOption === Automatic,

                (
                    center + #
                ) & /@
                    {
                        {0.11, 0.07},
                        {0.23, -0.13},
                        {-0.17, 0.19},
                        {-0.29, -0.08},
                        {0.06, -0.31},
                        {0.31, 0.22},
                        {-0.08, 0.27},
                        {0.18, 0.32},
                        {-0.26, -0.24},
                        {0.34, -0.03},
                        {-0.14, 0.35},
                        {0.27, -0.29}
                    },

                samplePointsOption
            ];


        If[
            ! ListQ[samplePoints]
            ||
            Length[samplePoints] == 0,

            Message[
                ModeAngularMomentum::samples
            ];

            Return[$Failed]
        ];


        (* Full Cartesian Fourier eigenvector of the selected mode. *)
        modeVector =
            modes["EigenvectorsFull"][[modeIndex]];

        modeCoefficients =
            Partition[
                modeVector,
                3
            ];


        (*
            At k = 0 the real-space phase contains only the Cartesian
            reciprocal-lattice vector G. Coordinates are in units of
            a0, so no additional factor of a0 is required.
        *)
        reciprocalWaveVectors =
            (
                toCRT[
                    G[
                        Sequence @@ #
                    ],
                    etaValue,
                    phiValue,
                    alphaValue
                ]
            ) & /@ nodes;


        (* Numerical real-space eigenfunction in Cartesian components. *)
        modeField[
            xValue_?NumericQ,
            yValue_?NumericQ
        ] :=
            Total[
                MapThread[
                    #1
                    Exp[
                        I #2 .
                        {xValue, yValue, 0}
                    ] &,

                    {
                        modeCoefficients,
                        reciprocalWaveVectors
                    }
                ]
            ];


        (* C6 rotation: counter-clockwise by 2 Pi / 6. *)
        rotationAngle =
            2 Pi / 6;

        rotation2DInverse =
            RotationMatrix[
                -rotationAngle
            ];

        rotation3D =
            RotationMatrix[
                rotationAngle,
                {0, 0, 1}
            ];


        (*
            Evaluate the original field and the fully rotated field

                R6 deltaM[R6^-1 (r-center) + center]

            at the same set of real-space points.
        *)
        originalFields =
            (
                modeField[
                    #[[1]],
                    #[[2]]
                ]
            ) & /@ samplePoints;


        rotatedFields =
            (
                Module[
                    {
                        relativePosition,
                        preimagePosition
                    },

                    relativePosition =
                        # - center;

                    preimagePosition =
                        center
                        +
                        rotation2DInverse .
                        relativePosition;

                    rotation3D .
                        modeField[
                            preimagePosition[[1]],
                            preimagePosition[[2]]
                        ]
                ]
            ) & /@ samplePoints;


        (* Normalization used for the symmetry mismatch. *)
        denominator =
            Total[
                Norm[#]^2
                & /@ originalFields
            ];


        If[
            Chop[denominator] == 0,

            Message[
                ModeAngularMomentum::zeronorm
            ];

            Return[$Failed]
        ];


        (* Test all six one-dimensional irreducible C6 characters. *)
        mValues =
            Range[
                -2,
                3
            ];


        residuals =
            Table[
                Sqrt[
                    Total[
                        MapThread[
                            Norm[
                                #1
                                -
                                Exp[
                                    I rotationAngle m
                                ]
                                #2
                            ]^2 &,

                            {
                                rotatedFields,
                                originalFields
                            }
                        ]
                    ]
                    /
                    denominator
                ],

                {
                    m,
                    mValues
                }
            ];


        (* Candidate with the smallest normalized mismatch. *)
        bestPosition =
            First[
                Ordering[
                    residuals,
                    1
                ]
            ];

        bestM =
            mValues[[
                bestPosition
            ]];

        bestResidual =
            residuals[[
                bestPosition
            ]];


        (* Accept the label only when the C6 relation is accurate. *)
        acceptedM =
            If[
                bestResidual <= symmetryTolerance,
                bestM,
                Missing["Unclassified"]
            ];


        signedM =
            acceptedM;


        classification =
            If[
                MissingQ[
                    acceptedM
                ],

                "Unclassified",

                "m = "
                <>
                ToString[
                    acceptedM
                ]
                <>
                " (mod 6)"
            ];


        <|
            "AngularMomentumSigned" ->
                signedM,

            "AngularMomentumModulo6" ->
                If[
                    MissingQ[
                        acceptedM
                    ],
                    Missing["Unclassified"],
                    Mod[
                        acceptedM,
                        6
                    ]
                ],

            "BestCandidateM" ->
                bestM,

            "BestResidual" ->
                bestResidual,

            "SymmetryTolerance" ->
                symmetryTolerance,

            "ResidualsByM" ->
                AssociationThread[
                    mValues,
                    residuals
                ],

            "RotationEigenvalue" ->
                If[
                    MissingQ[
                        acceptedM
                    ],

                    Missing["Unclassified"],

                    Exp[
                        I rotationAngle acceptedM
                    ]
                ],

            "RotationCenter" ->
                center,

            "SamplePoints" ->
                samplePoints,

            "Classification" ->
                classification
        |>
    ];

(* --------------------------------------------------------------- *)
(* Analyze angular momentum of all k = 0 modes                    *)
(* --------------------------------------------------------------- *)

Options[AnalyzeAngularMomentumModes] = {
    "WaveVectorTolerance" -> 10^-10,
    "SymmetryTolerance" -> 10^-6,
    "SamplePoints" -> Automatic,
    "RotationCenter" -> Automatic
};

(*
    Apply ModeAngularMomentum to every mode of a k = 0 eigensystem.

    The result is kept separate from AnalyzeDipoleModes so that the
    two classification layers can be validated independently. They
    can later be merged into a higher-level AnalyzeModes function.
*)
AnalyzeAngularMomentumModes[
    dynamicsData_Association,
    modes_Association,
    OptionsPattern[]
] :=
    Module[
        {
            numberOfModes,
            waveVectorTolerance,
            symmetryTolerance,
            samplePoints,
            rotationCenter,
            angularData
        },

        numberOfModes =
            Length[
                modes["Frequencies"]
            ];

        waveVectorTolerance =
            OptionValue[
                "WaveVectorTolerance"
            ];

        symmetryTolerance =
            OptionValue[
                "SymmetryTolerance"
            ];

        samplePoints =
            OptionValue[
                "SamplePoints"
            ];

        rotationCenter =
            OptionValue[
                "RotationCenter"
            ];


        angularData =
            Table[
                ModeAngularMomentum[
                    dynamicsData,
                    modes,
                    modeIndex,

                    "WaveVectorTolerance" ->
                        waveVectorTolerance,

                    "SymmetryTolerance" ->
                        symmetryTolerance,

                    "SamplePoints" ->
                        samplePoints,

                    "RotationCenter" ->
                        rotationCenter
                ],

                {
                    modeIndex,
                    1,
                    numberOfModes
                }
            ];


        If[
            MemberQ[
                angularData,
                $Failed
            ],

            Return[$Failed]
        ];


        <|
            "WaveVector" ->
                modes["WaveVector"],

            "Symmetry" ->
                "C6",

            "SymmetryTolerance" ->
                symmetryTolerance,

            "Modes" ->
                Table[
                    <|
                        "ModeIndex" ->
                            modeIndex,

                        "Frequency" ->
                            modes["Frequencies"][[modeIndex]],

                        "AngularMomentum" ->
                            angularData[[modeIndex]],

                        "AngularMomentumClassification" ->
                            angularData[[modeIndex]][
                                "AngularMomentumSigned"
                            ]
                    |>,

                    {
                        modeIndex,
                        1,
                        numberOfModes
                    }
                ]
        |>
    ];

(* --------------------------------------------------------------- *)
(* Help                                                            *)
(* --------------------------------------------------------------- *)

ModeAnalysisHelp[] :=
    Column[
        {
            Style[
                "ModeAnalysis.wl",
                Bold,
                16
            ],

            "Analysis and classification of magnetic magnon eigenmodes at k = 0 using dipole moments and C6 rotation symmetry.",

            Grid[
                {
                    {"Function", "Description"},
                    {"ModeDipoleMoment", "Extract the G = 0 Fourier component of one k = 0 magnon mode and compute its breathing, CW and CCW dipole weights."},
                    {"ClassifyDipoleMode", "Classify one dipole analysis as Breathing, CW, CCW or Dipole Inactive."},
                    {"AnalyzeDipoleModes", "Analyze and classify all modes of a k = 0 eigensystem using relative and/or absolute dipole-activity thresholds."},
                    {"ModeAngularMomentum", "Determine the C6 angular-momentum label m modulo 6 of one k = 0 mode from its real-space rotation symmetry."},
                    {"AnalyzeAngularMomentumModes", "Apply the C6 angular-momentum classification to every mode of a k = 0 eigensystem."}
                },
                Frame -> All,
                Alignment -> Left
            ]
        }
    ];

Print[
    ModeAnalysisHelp[]
];