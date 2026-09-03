(* ================================================================ *)
(* Plotting.wl                                                      *)
(*                                                                  *)
(* Real-space visualization of equilibrium skyrmion textures and   *)
(* magnetic magnon eigenfunctions.                                 *)
(*                                                                  *)
(* Dependencies:                                                    *)
(*     Geometry.wl                                                  *)
(*     FourierBasis.wl                                              *)
(*     GroundState.wl                                               *)
(*     DynamicsProjection.wl                                       *)
(*     DynamicsMagnetic.wl                                         *)
(* ================================================================ *)

ClearAll[
    MagnetizationAffine,
    MagnetizationCartesian,

    neighborVectors,
    Tpoint,
    WZC,

    SkyrmionConfigurationPlot,

    MagnonEigenfunctionCartesian,
    MagnonEigenfunctionPlot,

    PlottingHelp
];

(* --------------------------------------------------------------- *)
(* Equilibrium magnetization                                      *)
(* --------------------------------------------------------------- *)

(*
    Reconstruct the equilibrium magnetization at affine coordinates
    (xi1,xi2).

    The Fourier cutoff and magnetic representation are read directly
    from the saved equilibrium Association.
*)
MagnetizationAffine[
    xi1Coord_,
    xi2Coord_,
    equilibrium_Association
] :=
    Module[
        {
            cutoff,
            representation,
            solutionRules
        },

        cutoff =
            equilibrium["FourierCutoff"];

        representation =
            equilibrium["MagneticRepresentation"];

        solutionRules =
            EquilibriumRules[
                equilibrium
            ];

        Sum[
            MagneticFourierComponent[
                n1,
                n2,
                cutoff,
                representation
            ]
            Exp[
                I G[n1, n2] .
                {xi1Coord, xi2Coord, 0}
            ],

            {n1, -cutoff, cutoff},
            {n2, -cutoff, cutoff}
        ]
        /. solutionRules
    ];

(*
    Reconstruct the equilibrium magnetization in Cartesian components
    at Cartesian coordinates (x,y), expressed in units of a0.
*)
MagnetizationCartesian[
    x_,
    y_,
    equilibrium_Association
] :=
    Module[
        {
            solutionRules,
            etaValue,
            phiValue,
            alphaValue,
            xi1Value,
            xi2Value,
            mAffine
        },

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

        xi1Value =
            xi1[
                x,
                y,
                etaValue,
                phiValue,
                alphaValue
            ];

        xi2Value =
            xi2[
                x,
                y,
                etaValue,
                phiValue,
                alphaValue
            ];

        mAffine =
            MagnetizationAffine[
                xi1Value,
                xi2Value,
                equilibrium
            ];

        Chop[
            Re[
                toCRT[
                    mAffine,
                    etaValue,
                    phiValue,
                    alphaValue
                ]
            ]
        ]
    ];

(* --------------------------------------------------------------- *)
(* Wigner-Seitz cell                                              *)
(* --------------------------------------------------------------- *)

(* Shortest lattice vectors surrounding the origin. *)
neighborVectors[
    eta_,
    phi_,
    alphaR_
] :=
    SortBy[
        {
            e1[eta, phi, alphaR],
            e1[eta, phi, alphaR] + e2[eta, phi, alphaR],
            e1[eta, phi, alphaR] - e2[eta, phi, alphaR],
            e2[eta, phi, alphaR],

            -e1[eta, phi, alphaR],
            -(e1[eta, phi, alphaR] + e2[eta, phi, alphaR]),
            -(e1[eta, phi, alphaR] - e2[eta, phi, alphaR]),
            -e2[eta, phi, alphaR]
        }[[All, 1 ;; 2]],

        Norm
    ];

(* Vertices used to construct the Wigner-Seitz cell. *)
Tpoint[
    eta_,
    phi_,
    alphaR_
] :=
    SortBy[
        SortBy[
            Table[
                {
                    x0[eta, phi, alphaR],
                    y0[eta, phi, alphaR]
                }
                +
                neighborVectors[
                    eta,
                    phi,
                    alphaR
                ][[idx]],

                {idx, 1, 6}
            ],
            First
        ],
        Last
    ];

(* Wigner-Seitz cell of the magnetic lattice. *)
WZC[
    eta_,
    phi_,
    alphaR_
] :=
    RegionIntersection[

        Region[
            Triangle[
                {
                    Tpoint[eta, phi, alphaR][[1]],
                    Tpoint[eta, phi, alphaR][[4]],
                    Tpoint[eta, phi, alphaR][[5]]
                }
            ]
        ],

        Region[
            Triangle[
                {
                    Tpoint[eta, phi, alphaR][[2]],
                    Tpoint[eta, phi, alphaR][[3]],
                    Tpoint[eta, phi, alphaR][[6]]
                }
            ]
        ]
    ];

(* --------------------------------------------------------------- *)
(* Equilibrium configuration plot                                 *)
(* --------------------------------------------------------------- *)

Options[SkyrmionConfigurationPlot] = {
    ImageSize -> 600,
    "ShowFrameAndLabels" -> True,
    "ShowTitle" -> True,
    "Range" -> 1.0,
    "PlotPoints" -> 80,
    "VectorPoints" -> 19
};

(*
    Plot a saved equilibrium configuration.

    Color:
        out-of-plane magnetization mz.

    Arrows:
        in-plane Cartesian magnetization (mx,my).

    Dashed boundary:
        Wigner-Seitz cell.
*)
SkyrmionConfigurationPlot[
    equilibrium_Association,
    OptionsPattern[]
] :=
    Module[
        {
            solutionRules,
            dmiType,
            bValue,
            xiValue,
            etaValue,
            phiValue,
            alphaValue,
            centerX,
            centerY,
            imageSize,
            showFrameAndLabels,
            showTitle,
            range,
            plotPoints,
            vectorPoints,
            frameLabels,
            plotTitle,
            densityPlot,
            vectorPlot,
            cellBoundary
        },

        solutionRules =
            EquilibriumRules[
                equilibrium
            ];

        dmiType =
            equilibrium["DMIType"];

        bValue =
            b /. equilibrium["Parameters"];

        xiValue =
            xi /. equilibrium["Parameters"];

        etaValue =
            eta0 /. solutionRules;

        phiValue =
            phi0 /. solutionRules;

        alphaValue =
            alphaR0 /. solutionRules;

        (* Plot options. *)
        imageSize =
            OptionValue[ImageSize];

        showFrameAndLabels =
            OptionValue["ShowFrameAndLabels"];

        showTitle =
            OptionValue["ShowTitle"];

        range =
            OptionValue["Range"];

        plotPoints =
            OptionValue["PlotPoints"];

        vectorPoints =
            OptionValue["VectorPoints"];

        (* Center the plot around the primitive cell. *)
        centerX =
            x0[
                etaValue,
                phiValue,
                alphaValue
            ];

        centerY =
            y0[
                etaValue,
                phiValue,
                alphaValue
            ];

        frameLabels =
            If[
                showFrameAndLabels,
                {"x/a0", "y/a0"},
                None
            ];

        plotTitle =
            If[
                showTitle,
                Style[
                    Row[
                        {
                            dmiType,
                            " DMI   |   ",
                            "b = ",
                            NumberForm[N[bValue], {5, 3}],
                            "   |   ",
                            "\[Xi] = ",
                            NumberForm[N[xiValue], {5, 3}]
                        }
                    ],
                    15,
                    Bold
                ],
                None
            ];

        (* Out-of-plane magnetization. *)
        densityPlot =
            DensityPlot[
                MagnetizationCartesian[x, y, equilibrium][[3]],
                {
                    x,
                    centerX - range,
                    centerX + range
                },
                {
                    y,
                    centerY - range,
                    centerY + range
                },

                PlotPoints ->
                    plotPoints,

                ColorFunction ->
                    "LightTemperatureMap",

                ColorFunctionScaling ->
                    True,

                AspectRatio ->
                    1,

                Frame ->
                    showFrameAndLabels,

                FrameLabel ->
                    frameLabels,

                LabelStyle ->
                    16,

                PlotRange ->
                    All,

                PlotRangePadding ->
                    None,

                PlotLegends ->
                    BarLegend[
                        Automatic,
                        LegendLabel -> Subscript["m", "z"]
                    ]
            ];


        (* In-plane magnetization. *)
        vectorPlot =
            VectorPlot[
                MagnetizationCartesian[x, y, equilibrium][[1 ;; 2]],
                {
                    x,
                    centerX - range,
                    centerX + range
                },
                {
                    y,
                    centerY - range,
                    centerY + range
                },

                VectorScaling ->
                    Automatic,

                VectorAspectRatio ->
                    0.33,

                VectorColorFunction ->
                    None,

                VectorStyle ->
                    {GrayLevel[0.2], Thick},

                VectorPoints ->
                    vectorPoints,

                Frame ->
                    showFrameAndLabels,

                FrameLabel ->
                    frameLabels,

                LabelStyle ->
                    16,

                PlotRangePadding ->
                    None
            ];


        (* Wigner-Seitz-cell boundary. *)
        cellBoundary =
            RegionPlot[
                WZC[
                    etaValue,
                    phiValue,
                    alphaValue
                ],

                PlotStyle ->
                    None,

                BoundaryStyle ->
                    {Thick, Dashed, Black}
            ];


        Show[
            densityPlot,
            vectorPlot,
            cellBoundary,

            ImageSize ->
                imageSize,

            Frame ->
                showFrameAndLabels,

            PlotLabel ->
                plotTitle
        ]
    ];

(* --------------------------------------------------------------- *)
(* Magnon eigenfunction reconstruction                            *)
(* --------------------------------------------------------------- *)

MagnonEigenfunctionPlot::mode =
    "Mode index `1` is invalid. The supplied eigensystem contains `2` positive-frequency modes.";

MagnonEigenfunctionPlot::part =
    "Unknown plot part `1`. Use \"Real\", \"Imaginary\", or \"Both\".";

(*
    Reconstruct one normalized magnon eigenfunction in Cartesian
    components.

    The full eigenvectors returned by SolveMagneticModes are ordered
    according to dynamicsData["Nodes"] and already contain Cartesian
    magnetic components.

    Wavevector convention:
        k is the Cartesian Bloch wavevector stored in modes["WaveVector"].

        For reciprocal node G, the total physical wavevector is

            q = k + G_cart/a0.

    Plotting coordinates (x,y) are expressed in units of a0.
    Therefore the real-space phase is

        exp[i (a0 k + G_cart) . (x,y,0)].

    At k = 0 this reduces to the expression used in the legacy
    eigenfunction plots.
*)
MagnonEigenfunctionCartesian[
    x_,
    y_,
    dynamicsData_Association,
    modes_Association,
    modeIndex_Integer
] :=
    Module[
        {
            equilibrium,
            nodes,
            solutionRules,
            aValue,
            etaValue,
            phiValue,
            alphaValue,
            k,
            modeVector,
            modeCoefficients,
            phaseWaveVectors
        },

        equilibrium =
            dynamicsData["Equilibrium"];

        nodes =
            dynamicsData["Nodes"];

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

        k =
            modes["WaveVector"];

        modeVector =
            modes["EigenvectorsFull"][[modeIndex]];

        (* One Cartesian 3-vector for every reciprocal-space node. *)
        modeCoefficients =
            Partition[
                modeVector,
                3
            ];

        (*
            Dimensionless wavevectors appropriate for coordinates
            measured in units of a0.
        *)
        phaseWaveVectors =
            (
                aValue k
                +
                toCRT[G[Sequence @@ #],
                    etaValue,
                    phiValue,
                    alphaValue
                ]
            ) & /@ nodes;

        Chop[
            Total[
                MapThread[
                    #1
                    Exp[
                        I #2 .
                        {x, y, 0}
                    ] &,

                    {
                        modeCoefficients,
                        phaseWaveVectors
                    }
                ]
            ]
        ]
    ];

(* --------------------------------------------------------------- *)
(* Magnon eigenfunction plot                                      *)
(* --------------------------------------------------------------- *)

Options[MagnonEigenfunctionPlot] = {

    ImageSize -> 600,

    (* "Real", "Imaginary", or "Both". *)
    "Part" -> "Real",

    "ShowFrameAndLabels" -> True,

    "ShowTitle" -> True,

    "Range" -> 0.8,

    "PlotPoints" -> 80,

    "VectorPoints" -> 25,

    (*
        Maximum absolute value used for the symmetric dm_z color
        scale. Automatic estimates it numerically over the plot area.
    *)
    "ColorMaximum" -> Automatic,

    (* Number of sample points per axis used for Automatic scaling. *)
    "ColorSamplingPoints" -> 41
};

(*
    Plot one normalized magnon eigenfunction.

    Color:
        real or imaginary part of delta m_z.

    Arrows:
        corresponding in-plane components
        (delta m_x, delta m_y).

    "Part" -> "Both" returns real and imaginary panels side-by-side
    using the same symmetric color scale.
*)
MagnonEigenfunctionPlot[
    dynamicsData_Association,
    modes_Association,
    modeIndex_Integer,
    OptionsPattern[]
] :=
    Module[
        {
            equilibrium,
            nodes,
            solutionRules,
            dmiType,
            bValue,
            xiValue,
            etaValue,
            phiValue,
            alphaValue,
            aValue,
            centerX,
            centerY,
            k,
            modeVector,
            frequencyValue,
            modeCoefficients,
            phaseWaveVectors,
            imageSize,
            requestedPart,
            showFrameAndLabels,
            showTitle,
            range,
            plotPoints,
            vectorPoints,
            colorMaximum,
            colorSamplingPoints,
            frameLabels,
            xSamples,
            ySamples,
            modeField,
            makeSinglePlot
        },

        (* Validate mode index before accessing the eigenvector. *)
        If[
            modeIndex < 1
            ||
            modeIndex > modes["NumberOfModes"],

            Message[
                MagnonEigenfunctionPlot::mode,
                modeIndex,
                modes["NumberOfModes"]
            ];

            Return[$Failed]
        ];


        requestedPart =
            OptionValue["Part"];

        If[
            ! MemberQ[
                {"Real", "Imaginary", "Both"},
                requestedPart
            ],

            Message[
                MagnonEigenfunctionPlot::part,
                requestedPart
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


        (* Physical metadata used in the title. *)
        dmiType =
            equilibrium["DMIType"];

        bValue =
            b /. equilibrium["Parameters"];

        xiValue =
            xi /. equilibrium["Parameters"];


        (* Equilibrium lattice parameters. *)
        aValue =
            a0 /. solutionRules;

        etaValue =
            eta0 /. solutionRules;

        phiValue =
            phi0 /. solutionRules;

        alphaValue =
            alphaR0 /. solutionRules;

        centerX =
            x0[
                etaValue,
                phiValue,
                alphaValue
            ];

        centerY =
            y0[
                etaValue,
                phiValue,
                alphaValue
            ];

        (* Mode data. *)
        k =
            modes["WaveVector"];

        modeVector =
            modes["EigenvectorsFull"][[modeIndex]];

        frequencyValue =
            modes["Frequencies"][[modeIndex]];

        modeCoefficients =
            Partition[
                modeVector,
                3
            ];

        (*
            Because x and y are plotted in units of a0, multiply the
            Cartesian Bloch wavevector k by a0 before combining it
            with the dimensionless Cartesian reciprocal vector.
        *)
        phaseWaveVectors =
            (
                aValue k
                +
                toCRT[G[Sequence @@ #],
                    etaValue,
                    phiValue,
                    alphaValue
                ]
            ) & /@ nodes;

        (* Plot options. *)
        imageSize =
            OptionValue[ImageSize];

        showFrameAndLabels =
            OptionValue["ShowFrameAndLabels"];

        showTitle =
            OptionValue["ShowTitle"];

        range =
            OptionValue["Range"];

        plotPoints =
            OptionValue["PlotPoints"];

        vectorPoints =
            OptionValue["VectorPoints"];

        colorMaximum =
            OptionValue["ColorMaximum"];

        colorSamplingPoints =
            OptionValue["ColorSamplingPoints"];

        frameLabels =
            If[
                showFrameAndLabels,
                {"x/a0", "y/a0"},
                None
            ];

        (*
            Local numerical reconstruction used repeatedly by
            DensityPlot and VectorPlot. All mode-dependent quantities
            have already been extracted, avoiding repeated Association
            lookups and Fourier-vector construction.
        *)
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
                        phaseWaveVectors
                    }
                ]
            ];


        (* Estimate one common symmetric color scale for Re and Im. *)
        If[
            colorMaximum === Automatic,

            xSamples =
                Subdivide[
                    centerX - range,
                    centerX + range,
                    colorSamplingPoints - 1
                ];

            ySamples =
                Subdivide[
                    centerY - range,
                    centerY + range,
                    colorSamplingPoints - 1
                ];

            colorMaximum =
                Max[
                    Flatten[
                        Table[
                            Abs[
                                modeField[xValue, yValue][[3]]
                            ],
                            {xValue, xSamples},
                            {yValue, ySamples}
                        ]
                    ]
                ];
        ];


        (* Avoid a singular color scale for an exactly vanishing mz. *)
        If[
            ! NumericQ[colorMaximum]
            ||
            colorMaximum <= 0,

            colorMaximum =
                1;
        ];

        (*
            Construct either the real or imaginary eigenfunction
            panel. Both use the same colorMaximum.
        *)
        makeSinglePlot[
            part_String
        ] :=
            Module[
                {
                    partFunction,
                    partLabel,
                    plotTitle,
                    colorFunction,
                    densityPlot,
                    vectorPlot,
                    cellBoundary
                },

                {
                    partFunction,
                    partLabel
                } =
                    Switch[
                        part,

                        "Real",
                        {Re, "Real part"},

                        "Imaginary",
                        {Im, "Imaginary part"}
                    ];

                plotTitle =
                    If[
                        showTitle,
                        Style[
                            Row[
                                {
                                    partLabel,
                                    "   |   Mode ",
                                    modeIndex,
                                    "   |   ",
                                    "\[Omega]/",
                                    Subscript["\[Omega]", "c2"],
                                    " = ",
                                    NumberForm[
                                        N[frequencyValue],
                                        {6, 4}
                                    ],
                                    "   |   ",
                                    dmiType,
                                    " DMI   |   ",
                                    "b = ",
                                    NumberForm[
                                        N[bValue],
                                        {5, 3}
                                    ],
                                    "   |   ",
                                    "\[Xi] = ",
                                    NumberForm[
                                        N[xiValue],
                                        {5, 3}
                                    ]
                                }
                            ],
                            14,
                            Bold
                        ],
                        None
                    ];

                colorFunction =
                    Function[
                        value,
                        ColorData["LightTemperatureMap"][
                            Clip[
                                (value / colorMaximum + 1) / 2,
                                {0, 1}
                            ]
                        ]
                    ];

                densityPlot =
                    DensityPlot[
                        partFunction[
                            modeField[x, y][[3]]
                        ],

                        {
                            x,
                            centerX - range,
                            centerX + range
                        },

                        {
                            y,
                            centerY - range,
                            centerY + range
                        },

                        PlotPoints ->
                            plotPoints,

                        PlotRange ->
                            All,

                        ColorFunction ->
                            colorFunction,

                        ColorFunctionScaling ->
                            False,

                        AspectRatio ->
                            1,

                        Frame ->
                            showFrameAndLabels,

                        FrameLabel ->
                            frameLabels,

                        LabelStyle ->
                            16,

                        PlotRangePadding ->
                            None,

                        PlotLegends ->
                            BarLegend[
                                {
                                    colorFunction,
                                    {
                                        -colorMaximum,
                                        colorMaximum
                                    }
                                },

                                LegendLabel -> Subscript["\[Delta]m", "z"]
                            ]
                    ];


                vectorPlot =
                    VectorPlot[
                        partFunction[
                            modeField[
                                x,
                                y
                            ][[1 ;; 2]]
                        ],

                        {
                            x,
                            centerX - range,
                            centerX + range
                        },

                        {
                            y,
                            centerY - range,
                            centerY + range
                        },

                        VectorScaling ->
                            Automatic,

                        VectorAspectRatio ->
                            0.33,

                        VectorColorFunction ->
                            None,

                        VectorStyle ->
                            {
                                GrayLevel[0.25],
                                Thick
                            },

                        VectorPoints ->
                            vectorPoints,

                        Frame ->
                            showFrameAndLabels,

                        FrameLabel ->
                            frameLabels,

                        LabelStyle ->
                            16,

                        PlotRangePadding ->
                            None
                    ];


                cellBoundary =
                    RegionPlot[
                        WZC[
                            etaValue,
                            phiValue,
                            alphaValue
                        ],

                        PlotStyle ->
                            None,

                        BoundaryStyle ->
                            {
                                Thick,
                                Dashed,
                                Black
                            }
                    ];


                Show[
                    densityPlot,
                    vectorPlot,
                    cellBoundary,

                    Frame ->
                        showFrameAndLabels,

                    PlotLabel ->
                        plotTitle,

                    ImageSize ->
                        imageSize
                ]
            ];


        Switch[
            requestedPart,

            "Real",
            makeSinglePlot["Real"],

            "Imaginary",
            makeSinglePlot["Imaginary"],

            "Both",
            GraphicsRow[
                {
                    makeSinglePlot["Real"],
                    makeSinglePlot["Imaginary"]
                },

                ImageSize ->
                    2.5 imageSize
            ]
        ]
    ];

(* --------------------------------------------------------------- *)
(* Help                                                            *)
(* --------------------------------------------------------------- *)

PlottingHelp[] :=
    Column[
        {
            Style[
                "Plotting.wl",
                Bold,
                16
            ],

            "Real-space visualization of saved skyrmion equilibria and normalized magnetic magnon eigenfunctions.",

            Grid[
                {
                    {"Function", "Description"},
                    {"MagnetizationAffine", "Reconstruct the saved equilibrium magnetization in affine coordinates."},
                    {"MagnetizationCartesian", "Reconstruct the saved equilibrium magnetization in Cartesian components."},
                    {"WZC", "Construct the real-space Wigner-Seitz cell."},
                    {"SkyrmionConfigurationPlot", "Plot a saved equilibrium configuration with mz color map and in-plane magnetization arrows."},
                    {"MagnonEigenfunctionCartesian", "Reconstruct one normalized magnon eigenfunction in Cartesian components, including finite Bloch wavevector k."},
                    {"MagnonEigenfunctionPlot", "Plot the real, imaginary, or both parts of a magnon eigenfunction with common symmetric color scaling."},
                    {"neighborVectors / Tpoint", "Internal geometric helpers used to construct the Wigner-Seitz cell."}
                },
                Frame -> All,
                Alignment -> Left
            ]
        }
    ];

Print[
    PlottingHelp[]
];