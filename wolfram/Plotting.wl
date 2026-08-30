(* ================================================================ *)
(* Plotting.wl                                                      *)
(*                                                                  *)
(* Real-space reconstruction and visualization of saved            *)
(* skyrmion-lattice equilibrium configurations.                     *)
(* ================================================================ *)

(*
    Dependencies:
        Geometry.wl
        FourierBasis.wl
        GroundState.wl

    The main plotting routine expects a saved equilibrium Association
    produced by SaveMagneticEquilibrium and loaded using
    LoadMagneticEquilibrium.

    Typical workflow:

        equilibrium =
            LoadMagneticEquilibrium["configuration_name.txt"];

        SkyrmionConfigurationPlot[equilibrium];
*)


ClearAll[
    MagnetizationAffine,
    MagnetizationCartesian,

    neighborVectors,
    Tpoint,
    WZC,

    SkyrmionConfigurationPlot,
    PlottingHelp
];


(* --------------------------------------------------------------- *)
(* Real-space magnetization                                       *)
(* --------------------------------------------------------------- *)

(*
    Reconstruct the magnetic field from its Fourier components
    at affine coordinates (xi1, xi2).

    solutionRules must contain both the optimized Fourier
    coefficients and the fixed physical/lattice parameters.

    The returned components are covariant affine components.
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

    The Cartesian position is first transformed to affine
    coordinates, the Fourier series is reconstructed there,
    and the resulting field is transformed back to Cartesian
    components.

    Small numerical imaginary parts are removed using Re and Chop.
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


        (* Cartesian -> affine coordinates. *)
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


        (* Reconstruct Fourier field in affine components. *)
        mAffine =
            MagnetizationAffine[
                xi1Value,
                xi2Value,
                solutionRules
            ];


        (* Affine -> Cartesian field components. *)
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

(*
    Return the shortest lattice vectors surrounding the origin.
    This is an internal helper used to construct the Wigner-Seitz
    cell.
*)
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


(*
    Construct the vertices used to define the Wigner-Seitz cell.
    Internal helper for WZC.
*)
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


(*
    Construct the real-space Wigner-Seitz cell of the magnetic
    lattice.
*)
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

    (* Overall size of the final figure. *)
    ImageSize -> 600,

    (* Show frame, coordinate labels and physical-parameter title. *)
    "ShowFrameAndLabels" -> True,

    (* Cartesian plotting range around the unit-cell center. *)
    "Range" -> 1.0,

    (* Sampling density of the mz density plot. *)
    "PlotPoints" -> 80,

    (* Sampling density of the in-plane vector field. *)
    "VectorPoints" -> 19
};


(*
    Plot a saved magnetic equilibrium configuration.

    Input:
        equilibrium

            Association produced by SaveMagneticEquilibrium and
            loaded using LoadMagneticEquilibrium.

    Visualization:
        - color map: out-of-plane magnetization mz;
        - arrows: in-plane Cartesian magnetization (mx,my);
        - dashed line: Wigner-Seitz-cell boundary.

    When "ShowFrameAndLabels" -> True, the plot also displays:
        - Cartesian coordinate labels;
        - DMI type;
        - magnetic field b;
        - dipolar strength xi.
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
            range,

            imageSize,
            showFrameAndLabels,
            plotPoints,
            vectorPoints,

            frameLabels,
            plotTitle,

            densityPlot,
            vectorPlot,
            cellBoundary
        },


        (* ------------------------------------------------------- *)
        (* Extract saved equilibrium data                         *)
        (* ------------------------------------------------------- *)

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


        (* ------------------------------------------------------- *)
        (* Equilibrium lattice parameters                         *)
        (* ------------------------------------------------------- *)

        etaValue =
            eta0 /. solutionRules;

        phiValue =
            phi0 /. solutionRules;

        alphaValue =
            alphaR0 /. solutionRules;


        (* ------------------------------------------------------- *)
        (* Plot options                                            *)
        (* ------------------------------------------------------- *)

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


        (* ------------------------------------------------------- *)
        (* Plot center                                             *)
        (* ------------------------------------------------------- *)

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


        (* ------------------------------------------------------- *)
        (* Frame labels and title                                  *)
        (* ------------------------------------------------------- *)

        frameLabels =
            If[
                showFrameAndLabels,
                {"x/a0", "y/a0"},
                None
            ];


        plotTitle =
            If[
                showFrameAndLabels,

                Style[
                    Row[
                        {
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
                    16,
                    Bold
                ],

                None
            ];


        (* ------------------------------------------------------- *)
        (* Out-of-plane magnetization                             *)
        (* ------------------------------------------------------- *)

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
                        LegendLabel -> "mz"
                    ]
            ];


        (* ------------------------------------------------------- *)
        (* In-plane magnetization                                 *)
        (* ------------------------------------------------------- *)

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

                VectorColorFunction ->
                    None,

                VectorStyle ->
                    {
                        GrayLevel[0.2],
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


        (* ------------------------------------------------------- *)
        (* Wigner-Seitz-cell boundary                             *)
        (* ------------------------------------------------------- *)

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


        (* ------------------------------------------------------- *)
        (* Final combined plot                                     *)
        (* ------------------------------------------------------- *)

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
(* Help                                                            *)
(* --------------------------------------------------------------- *)

PlottingHelp[] :=
    Column[
        {
            Style["Plotting.wl", Bold, 16],
            "Real-space reconstruction and visualization of saved magnetic skyrmion-lattice equilibria.",
            Grid[
                {
                    {"Function", "Description"},
                    {"MagnetizationAffine", "Reconstruct the magnetization in covariant affine components."},
                    {"MagnetizationCartesian", "Reconstruct the magnetization in Cartesian components at (x,y)."},
                    {"WZC", "Construct the real-space Wigner-Seitz cell."},
                    {"SkyrmionConfigurationPlot", "Plot a loaded saved equilibrium, including mz, in-plane magnetization and the Wigner-Seitz boundary."},
                    {"neighborVectors / Tpoint", "Internal geometric helpers used to construct the Wigner-Seitz cell."}
                },
                Frame -> All,
                Alignment -> Left
            ]
        }
    ];

Print[PlottingHelp[]];