(* ================================================================ *)
(* Plotting.wl                                                      *)
(*                                                                  *)
(* Real-space reconstruction and visualization of skyrmion fields. *)
(* ================================================================ *)


ClearAll[
    MagnetizationAffine,
    MagnetizationCartesian,

    neighborVectors,
    Tpoint,
    WZC,

    SkyrmionConfigurationPlot
];


(* --------------------------------------------------------------- *)
(* Real-space magnetization                                       *)
(* --------------------------------------------------------------- *)

(*
    Reconstruct the magnetic field from its Fourier components
    at affine coordinates (xi1, xi2).
*)
MagnetizationAffine[
    xi1Coord_,
    xi2Coord_,
    solutionRules_List
] :=
    Sum[
        M[n1, n2]
        Exp[
            I G[n1, n2] . {xi1Coord, xi2Coord, 0}
        ],

        {n1, -$FourierCutoff, $FourierCutoff},
        {n2, -$FourierCutoff, $FourierCutoff}
    ] /. solutionRules;


(*
    Reconstruct the magnetic field in Cartesian components
    at Cartesian coordinates (x,y).
*)
MagnetizationCartesian[
    x_,
    y_,
    solutionRules_List
] :=
    Module[
        {
            etaValue,
            phiValue,
            alphaValue,
            xi1Value,
            xi2Value,
            mAffine
        },

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
                solutionRules
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
(* Skyrmion configuration plot                                    *)
(* --------------------------------------------------------------- *)

Options[SkyrmionConfigurationPlot] = {

    ImageSize -> 500,

    "ShowFrameAndLabels" -> True,

    "Range" -> 1.0,

    "PlotPoints" -> 80,

    "VectorPoints" -> 19
};


(*
    Plot the magnetic configuration.

    Color:
        out-of-plane magnetization mz.

    Arrows:
        in-plane Cartesian magnetization (mx,my).

    Dashed boundary:
        Wigner-Seitz cell.
*)
SkyrmionConfigurationPlot[
    solutionRules_List,
    OptionsPattern[]
] :=
    Module[
        {
            etaValue,
            phiValue,
            alphaValue,

            centerX,
            centerY,
            range,

            imageSize,
            showFrameAndLabels,
            plotPoints,
            vectorPoints,

            frameLabels,
            densityPlot,
            vectorPlot,
            cellBoundary
        },


        (* Equilibrium lattice parameters. *)
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


        (* Out-of-plane magnetization. *)
        densityPlot =
            DensityPlot[
                Evaluate[
                    MagnetizationCartesian[
                        x,
                        y,
                        solutionRules
                    ][[3]]
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

                PlotPoints -> plotPoints,

                ColorFunction -> "LightTemperatureMap",
                ColorFunctionScaling -> True,

                AspectRatio -> 1,

                Frame -> showFrameAndLabels,
                FrameLabel -> frameLabels,

                LabelStyle -> 16,

                PlotRange -> All,
                PlotRangePadding -> None,

                PlotLegends ->
                    BarLegend[
                        Automatic,
                        LegendLabel -> "mz"
                    ]
            ];


        (* In-plane magnetization. *)
        vectorPlot =
            VectorPlot[
                Evaluate[
                    MagnetizationCartesian[
                        x,
                        y,
                        solutionRules
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

                VectorScale ->
                    {0.03, Automatic, Automatic},

                VectorColorFunction -> None,

                VectorStyle ->
                    {GrayLevel[0.2], Thick},

                VectorPoints -> vectorPoints,

                Frame -> showFrameAndLabels,
                FrameLabel -> frameLabels,

                LabelStyle -> 16,

                PlotRangePadding -> None
            ];


        (* Wigner-Seitz-cell boundary. *)
        cellBoundary =
            RegionPlot[
                WZC[
                    etaValue,
                    phiValue,
                    alphaValue
                ],

                PlotStyle -> None,

                BoundaryStyle ->
                    {Thick, Dashed, Black}
            ];


        Show[
            densityPlot,
            vectorPlot,
            cellBoundary,

            ImageSize -> imageSize,

            Frame -> showFrameAndLabels
        ]
    ];