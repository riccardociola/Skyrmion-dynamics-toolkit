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
    MagneticFieldDirectionTag,
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
    Cache names:
        magnetic_bloch_N2_sixfold.txt
        magnetic_bloch_N2_general_Kcub_H001.txt
        magnetic_bloch_N2_general_Kcub_H110.txt

    Isotropic models carry no anisotropy or field-direction suffix.
    With cubic anisotropy enabled, H identifies the normalized field
    direction in crystal coordinates. The cubic-axis-equivalent
    [100], [010], [001] directions are all represented as H001.

    An integer or rational Miller direction is reduced to coprime
    integer components. Simple single-digit directions produce short
    names like H110; larger or signed indices are encoded without
    ambiguity. Non-rational directions receive a stable short hash of
    the normalized direction, avoiding collisions from rounding.

    kappaCub remains symbolic, so its magnitude is not part of the
    filename. Both magnetic representations are supported independently
    of whether cubic anisotropy is enabled.
*)
MagneticFieldDirectionTag[
    fieldDirection_List
] :=
    Module[
        {
            direction,
            exactDirection,
            integerDirection,
            commonDenominator,
            commonDivisor,
            coordinates,
            hash
        },

        direction =
            CanonicalCubicFieldDirection[fieldDirection];

        If[direction === $Failed,
            Return[$Failed]
        ];

        (* All cubic-axis-equivalent directions have one cache name. *)
        If[direction === {0, 0, 1},
            Return["H001"]
        ];

        (* Avoid introducing approximate rational equivalences. *)
        exactDirection =
            VectorQ[
                fieldDirection,
                (IntegerQ[#] || Head[#] === Rational) &
            ];

        If[exactDirection,
            commonDenominator =
                Apply[LCM, Denominator /@ fieldDirection];

            integerDirection =
                commonDenominator fieldDirection;

            commonDivisor =
                Apply[GCD, Abs[integerDirection]];

            integerDirection =
                integerDirection / commonDivisor;

            coordinates =
                If[
                    Max[Abs[integerDirection]] <= 9,
                    StringJoin[
                        (If[# < 0, "m" <> ToString[-#], ToString[#]] &)
                            /@ integerDirection
                    ],
                    StringRiffle[
                        (If[# < 0, "m" <> ToString[-#], ToString[#]] &)
                            /@ integerDirection,
                        "_"
                    ]
                ];

            Return["H" <> coordinates];
        ];

        (* Arbitrary real/symbolic directions: retain a unique tag. *)
        hash =
            StringTake[
                IntegerString[
                    Hash[direction, "SHA256"],
                    16,
                    64
                ],
                12
            ];

        "Hdir_" <> hash
    ];

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
            suffix,
            directionTag
        },

        representation =
            OptionValue["MagneticRepresentation"];

        useMCA =
            TrueQ[OptionValue["CubicAnisotropy"]];

        fieldDirection =
            OptionValue["FieldDirection"];

        suffix = "";

        If[useMCA,
            If[!ListQ[fieldDirection] || Length[fieldDirection] != 3,
                Message[BuildMagneticEnergy::field, fieldDirection];
                Return[$Failed]
            ];

            directionTag =
                MagneticFieldDirectionTag[fieldDirection];

            If[directionTag === $Failed,
                Return[$Failed]
            ];

            suffix = "_Kcub_" <> directionTag;
        ];

        FileNameJoin[
            {
                $ProjectRoot,
                "cache",
                "energies",
                StringJoin[
                    "magnetic_",
                    ToLowerCase[dmiType],
                    "_N", ToString[cutoff],
                    "_", ToLowerCase[representation],
                    suffix,
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
                    {"MagneticEnergyCacheFile", "Return a descriptive cache path; append _Kcub_H... only when cubic MCA is enabled."},
                    {"MagneticFieldDirectionTag", "Build a stable crystallographic direction label for the MCA cache."},
                    {"SimplifyMagneticEnergyContribution", "Simplify one term with its own time limit; retain the original term on timeout."}
                },
                Frame -> All,
                Alignment -> Left
            ],

            "Main options: MagneticRepresentation, CubicAnisotropy, FieldDirection, SimplificationTimeLimit, ForceRebuild."
        }
    ];

Print[BuildEnergyHelp[]];