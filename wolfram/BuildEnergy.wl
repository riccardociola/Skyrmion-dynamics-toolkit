(* ================================================================ *)
(* BuildEnergy.wl                                                   *)
(*                                                                  *)
(* Construction, per-term simplification, and caching of symbolic   *)
(* magnetic free energies.                                         *)
(*                                                                  *)
(* Dependencies:                                                    *)
(*   Geometry.wl, FourierBasis.wl, EnergyMagnetic.wl                *)
(*   EnergyAnisotropy.wl (only when cubic MCA is enabled).           *)
(* ================================================================ *)

ClearAll[
    BuildMagneticEnergy,
    MagneticEnergyCacheFile,
    SimplifyMagneticEnergyContribution,
    BuildEnergyHelp
];

(* --------------------------------------------------------------- *)
(* Options                                                         *)
(* --------------------------------------------------------------- *)

Options[BuildMagneticEnergy] = {
    "MagneticRepresentation" -> "SixFold",
    "CubicAnisotropy" -> False,
    "FieldDirection" -> {0, 0, 1},
    "SimplificationTimeLimit" -> 200,
    "ForceRebuild" -> False
};

Options[MagneticEnergyCacheFile] = Options[BuildMagneticEnergy];

BuildMagneticEnergy::dmi =
    "Unknown DMI type `1`. Use \"Bloch\", \"Neel\", or \"Anti\".";

BuildMagneticEnergy::rep =
    "Unknown magnetic representation `1`. Use \"SixFold\" or \"General\".";

BuildMagneticEnergy::mca =
    "EnergyAnisotropy.wl must be loaded to enable cubic anisotropy.";

BuildMagneticEnergy::fourier =
    "FourierBasis.wl does not provide the MagneticFourierComponent function required by CubicAnisotropyEnergy.";

BuildMagneticEnergy::field =
    "Invalid crystal field direction `1`. It must be a nonzero three-component vector.";

BuildMagneticEnergy::limit =
    "SimplificationTimeLimit must be positive or Infinity. Received `1`.";

BuildMagneticEnergy::cache =
    "Unable to write the energy cache file `1`.";

(* --------------------------------------------------------------- *)
(* Cache location and model signature                              *)
(* --------------------------------------------------------------- *)

(* Evaluated when this source file is loaded. *)
$BuildEnergyDirectory =
    DirectoryName[
        ExpandFileName[$InputFileName]
    ];

$ProjectRoot =
    DirectoryName[$BuildEnergyDirectory];

(*
    Cache version 2 deliberately does not reuse old, SixFold-only
    caches. Cache identity includes the Fourier representation and
    whether cubic MCA is present. If MCA is present, it also includes
    the normalized, canonicalized field direction.

    Symmetry-equivalent [100], [010], [001] orientations are thus
    assigned the same model cache. Strength kappaCub stays symbolic,
    so its numerical value does not need to enter the cache name.
*)
MagneticEnergyCacheFile[
    cutoff_Integer?NonNegative,
    dmiType_String,
    OptionsPattern[]
] :=
    Module[
        {
            representation,
            useMCA,
            fieldDirection,
            canonicalDirection,
            orientationTag,
            anisotropyTag
        },

        representation =
            OptionValue["MagneticRepresentation"];

        useMCA =
            TrueQ[OptionValue["CubicAnisotropy"]];

        fieldDirection =
            OptionValue["FieldDirection"];

        anisotropyTag =
            If[useMCA, "cub", "nocub"];

        orientationTag =
            If[
                useMCA,

                canonicalDirection =
                    CanonicalCubicFieldDirection[fieldDirection];

                If[canonicalDirection === $Failed,
                    Return[$Failed]
                ];

                "_H" <>
                    StringTake[
                        IntegerString[
                            Hash[canonicalDirection, "SHA256"],
                            16,
                            64
                        ],
                        12
                    ],

                ""
            ];

        FileNameJoin[
            {
                $ProjectRoot,
                "cache",
                "energies",
                StringJoin[
                    "magnetic_v2_",
                    ToLowerCase[dmiType],
                    "_N", ToString[cutoff],
                    "_", ToLowerCase[representation],
                    "_", anisotropyTag,
                    orientationTag,
                    ".txt"
                ]
            }
        ]
    ];

(* --------------------------------------------------------------- *)
(* Simplify one contribution independently                         *)
(* --------------------------------------------------------------- *)

(*
    Apply a separate time limit to symbolic Simplify of each term.
    On timeout, return that term in its unsimplified form instead of
    discarding it. The time limit does not include initial symbolic
    construction of the term.
*)
SimplifyMagneticEnergyContribution[
    label_String,
    contribution_,
    assumptions_,
    timeLimit_
] :=
    Module[
        {simplified, timedOut = False},

        Print["Simplifying ", label, " (limit: ", timeLimit, " s)..."];

        simplified =
            TimeConstrained[
                Simplify[
                    contribution,
                    Assumptions -> assumptions
                ],
                timeLimit,
                (
                    timedOut = True;
                    contribution
                )
            ];

        If[timedOut,
            Print[
                "  Time limit reached for ", label,
                "; retaining unsimplified contribution."
            ]
        ];

        simplified
    ];

(* --------------------------------------------------------------- *)
(* Magnetic energy construction                                    *)
(* --------------------------------------------------------------- *)

(*
    Usage:
      BuildMagneticEnergy[2, "Bloch"]

      BuildMagneticEnergy[
          2, "Bloch",
          "MagneticRepresentation" -> "General",
          "CubicAnisotropy" -> True,
          "FieldDirection" -> {1, 1, 0}
      ]

    The cubic anisotropy coefficient kappaCub is kept symbolic.
    Its value is assigned later through parameter replacement rules.

    Either Fourier representation may be used with or without MCA.
    alphaR0 is not fixed: it remains available to the minimizer.
*)
BuildMagneticEnergy[
    cutoff_Integer?NonNegative,
    dmiType_String : "Bloch",
    OptionsPattern[]
] :=
    Module[
        {
            representation,
            useMCA,
            fieldDirection,
            timeLimit,
            forceRebuild,
            canonicalDirection,

            cacheFile,
            dmiEnergy,
            nijCartesian,
            assumptions,
            terms,
            energy
        },

        representation =
            OptionValue["MagneticRepresentation"];

        useMCA =
            TrueQ[OptionValue["CubicAnisotropy"]];

        fieldDirection =
            OptionValue["FieldDirection"];

        timeLimit =
            OptionValue["SimplificationTimeLimit"];

        forceRebuild =
            TrueQ[OptionValue["ForceRebuild"]];

        If[!MemberQ[{"SixFold", "General"}, representation],
            Message[BuildMagneticEnergy::rep, representation];
            Return[$Failed]
        ];

        If[!MemberQ[{"Bloch", "Neel", "Anti"}, dmiType],
            Message[BuildMagneticEnergy::dmi, dmiType];
            Return[$Failed]
        ];

        If[
            !(timeLimit === Infinity ||
                (NumericQ[timeLimit] && TrueQ[timeLimit > 0])),
            Message[BuildMagneticEnergy::limit, timeLimit];
            Return[$Failed]
        ];

        If[useMCA,
            If[Length[DownValues[CubicAnisotropyEnergy]] == 0 ||
                Length[DownValues[CanonicalCubicFieldDirection]] == 0,
                Message[BuildMagneticEnergy::mca];
                Return[$Failed]
            ];

            (* Do not silently build symbolic energy with an
               unevaluated Fourier component. *)
            If[Length[DownValues[MagneticFourierComponent]] == 0,
                Message[BuildMagneticEnergy::fourier];
                Return[$Failed]
            ];

            If[!ListQ[fieldDirection] || Length[fieldDirection] != 3,
                Message[BuildMagneticEnergy::field, fieldDirection];
                Return[$Failed]
            ];

            canonicalDirection =
                CanonicalCubicFieldDirection[fieldDirection];

            If[canonicalDirection === $Failed,
                Message[BuildMagneticEnergy::field, fieldDirection];
                Return[$Failed]
            ];
        ];

        (* Configure the same active Fourier field used by the
           existing EnergyMagnetic.wl contributions. *)
        If[ConfigureFourierBasis[cutoff, representation] === $Failed,
            Return[$Failed]
        ];

        cacheFile =
            MagneticEnergyCacheFile[
                cutoff,
                dmiType,
                "MagneticRepresentation" -> representation,
                "CubicAnisotropy" -> useMCA,
                "FieldDirection" -> fieldDirection
            ];

        If[cacheFile === $Failed, Return[$Failed]];

        If[FileExistsQ[cacheFile] && !forceRebuild,
            Print["Loading cached magnetic energy: ", cacheFile];
            Return[Get[cacheFile]]
        ];

        Print[
            "Constructing energy: N = ", cutoff,
            ", DMI = ", dmiType,
            ", representation = ", representation,
            ", cubic MCA = ", useMCA
        ];

        assumptions = {
            a0 > 0,
            eta0 > 0,
            -Pi <= phi0 < Pi,
            0 <= alphaR0 <= 2 Pi,
            xi >= 0,
            lambda > 0,
            0 <= Nx <= 1,
            0 <= Ny <= 1,
            0 <= Nz <= 1,
            Nx + Ny + Nz == 1,
            Element[b, Reals],
            Element[kappaCub, Reals]
        };

        (* Fixed Cartesian demagnetization tensor in the lab frame. *)
        nijCartesian = DiagonalMatrix[{Nx, Ny, Nz}];

        dmiEnergy =
            Switch[
                dmiType,
                "Bloch", BlochDMIEnergy[a0, eta0, phi0, alphaR0],
                "Neel",  NeelDMIEnergy[a0, eta0, phi0, alphaR0],
                "Anti",  AntiDMIEnergy[a0, eta0, phi0, alphaR0]
            ];

        (* Simplify each physical contribution separately. Never
           perform a second, global Simplify on the entire sum. *)
        terms = {
            SimplifyMagneticEnergyContribution[
                "exchange",
                ExchangeEnergy[a0, eta0, phi0, alphaR0],
                assumptions, timeLimit
            ],

            SimplifyMagneticEnergyContribution[
                "DMI",
                dmiEnergy,
                assumptions, timeLimit
            ],

            SimplifyMagneticEnergyContribution[
                "Zeeman",
                ZeemanEnergy[b],
                assumptions, timeLimit
            ],

            SimplifyMagneticEnergyContribution[
                "magnetostatic",
                MagnetostaticInfiniteEnergy[
                    xi, nijCartesian, eta0, phi0, alphaR0
                ],
                assumptions, timeLimit
            ],

            SimplifyMagneticEnergyContribution[
                "constraint",
                ConstraintEnergy[lambda, eta0, phi0, alphaR0],
                assumptions, timeLimit
            ]
        };

        If[useMCA,
            AppendTo[
                terms,
                SimplifyMagneticEnergyContribution[
                    "cubic anisotropy",
                    CubicAnisotropyEnergy[
                        kappaCub,
                        canonicalDirection,
                        cutoff,
                        eta0, phi0, alphaR0,
                        representation
                    ],
                    assumptions, timeLimit
                ]
            ]
        ];

        energy = Total[terms];

        If[!DirectoryQ[DirectoryName[cacheFile]],
            CreateDirectory[
                DirectoryName[cacheFile],
                CreateIntermediateDirectories -> True
            ]
        ];

        Print["Saving symbolic energy to: ", cacheFile];
        If[
            Check[Put[energy, cacheFile]; True, False] === False,
            Message[BuildMagneticEnergy::cache, cacheFile]
        ];

        energy
    ];

(* --------------------------------------------------------------- *)
(* Help                                                            *)
(* --------------------------------------------------------------- *)

BuildEnergyHelp[] :=
    Column[
        {
            Style["BuildEnergy.wl", Bold, 16],
            "Assembly, per-contribution simplification and caching of symbolic magnetic energies.",

            Grid[
                {
                    {"Function", "Description"},
                    {"BuildMagneticEnergy", "Build/cache the energy for SixFold or General magnetization, optionally including cubic MCA."},
                    {"MagneticEnergyCacheFile", "Return the cache path for the specified DMI, cutoff, Fourier representation and MCA field orientation."},
                    {"SimplifyMagneticEnergyContribution", "Simplify one term with its own time limit; retain the original term on timeout."}
                },
                Frame -> All,
                Alignment -> Left
            ],

            "Main options: MagneticRepresentation, CubicAnisotropy, FieldDirection, SimplificationTimeLimit, ForceRebuild."
        }
    ];

Print[BuildEnergyHelp[]];