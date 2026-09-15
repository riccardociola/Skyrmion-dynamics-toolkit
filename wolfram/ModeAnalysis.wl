(* ================================================================ *)
(* ModeAnalysis.wl                                                  *)
(*                                                                  *)
(* Analysis and classification of magnetic magnon eigenmodes.      *)
(*                                                                  *)
(* First implementation:                                           *)
(*     dipole-moment classification at k = 0.                      *)
(*                                                                  *)
(* Future extensions can add angular-momentum analysis without      *)
(* changing the dipole-analysis interface defined here.             *)
(*                                                                  *)
(* Dependencies:                                                    *)
(*     DynamicsProjection.wl                                       *)
(*     DynamicsMagnetic.wl                                         *)
(* ================================================================ *)

ClearAll[
    ModeDipoleMoment,
    ClassifyDipoleMode,
    AnalyzeDipoleModes,
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

(* --------------------------------------------------------------- *)
(* Dipole moment of one mode                                      *)
(* --------------------------------------------------------------- *)

Options[ModeDipoleMoment] = {
    "WaveVectorTolerance" -> 10^-10,
    "TimeConvention" -> "Exp[-I omega t]"
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
    weight. Modes below it are labeled "DipoleInactive".

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
                "DipoleInactive"
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

            "Analysis and classification of magnetic magnon eigenmodes. The current version implements dipole-moment classification at k = 0 and is structured for later angular-momentum analysis.",

            Grid[
                {
                    {"Function", "Description"},
                    {"ModeDipoleMoment", "Extract the G = 0 Fourier component of one k = 0 magnon mode and compute its breathing, CW and CCW dipole weights."},
                    {"ClassifyDipoleMode", "Classify one dipole analysis as Breathing, CW, CCW or DipoleInactive."},
                    {"AnalyzeDipoleModes", "Analyze and classify all modes of a k = 0 eigensystem using relative and/or absolute dipole-activity thresholds."}
                },
                Frame -> All,
                Alignment ->Left
            ]
        }
    ];

Print[
    ModeAnalysisHelp[]
];