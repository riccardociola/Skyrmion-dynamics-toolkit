(* ================================================================ *)
(* BuildEnergy.wl                                                   *)
(*                                                                  *)
(* Construction and caching of symbolic energy functionals.         *)
(* ================================================================ *)


ClearAll[
    BuildMagneticEnergy,
    MagneticEnergyCacheFile
];


(* --------------------------------------------------------------- *)
(* Cache filename                                                  *)
(* --------------------------------------------------------------- *)

(* Directory containing this source file, evaluated when the file is loaded. *)
$BuildEnergyDirectory =
    DirectoryName[
        ExpandFileName[$InputFileName]
    ];

(* Root directory of the project. *)
$ProjectRoot =
    DirectoryName[$BuildEnergyDirectory];

(* Path of the cached symbolic magnetic energy. *)
MagneticEnergyCacheFile[
    cutoff_Integer,
    dmiType_String
] :=
    FileNameJoin[
        {
            $ProjectRoot,
            "cache",
            "energies",
            StringJoin[
                "magnetic_",
                ToLowerCase[dmiType],
                "_N",
                ToString[cutoff],
                ".txt"
            ]
        }
    ];


(* --------------------------------------------------------------- *)
(* Magnetic energy construction                                   *)
(* --------------------------------------------------------------- *)

BuildMagneticEnergy[
    cutoff_Integer,
    dmiType_String : "Bloch"
] :=
    Module[
        {
            cacheFile,
            dmiEnergy,
            energy
        },

        ConfigureFourierBasis[
            cutoff,
            "SixFold"
        ];

        cacheFile =
            MagneticEnergyCacheFile[
                cutoff,
                dmiType
            ];


        (* Load already constructed symbolic energy if available. *)
        If[
            FileExistsQ[cacheFile],

            Print[
                "Loading cached magnetic energy: ",
                cacheFile
            ];

            Return[
                Get[cacheFile]
            ];
        ];


        Print[
            "Constructing magnetic energy for cutoff N = ",
            cutoff,
            "..."
        ];


        (* Select the appropriate DMI contribution. *)
        dmiEnergy =
            Switch[
                dmiType,

                "Bloch",
                BlochDMIEnergy[
                    a0,
                    eta0,
                    phi0,
                    alphaR0
                ],

                "Neel",
                NeelDMIEnergy[
                    a0,
                    eta0,
                    phi0,
                    alphaR0
                ],

                "Anti",
                AntiDMIEnergy[
                    a0,
                    eta0,
                    phi0,
                    alphaR0
                ],

                _,
                Message[
                    BuildMagneticEnergy::dmi,
                    dmiType
                ];
                Return[$Failed]
            ];

        (* Symbolic Cartesian demagnetization tensor. *)
        NijCartesian =
            DiagonalMatrix[
                {Nx, Ny, Nz}
            ];

        energy =
            ExchangeEnergy[
                a0,
                eta0,
                phi0,
                alphaR0
            ]
            +
            dmiEnergy
            +
            ZeemanEnergy[b]
            +
            MagnetostaticInfiniteEnergy[
                xi,
                NijCartesian,
                eta0,
                phi0,
                alphaR0
            ]
            +
            ConstraintEnergy[
                lambda,
                eta0,
                phi0,
                alphaR0
            ];

        Print["Simplifying symbolic energy..."];

        energy =
            TimeConstrained[
                Simplify[
                    energy,
                    Assumptions -> {
                        a0 > 0,
                        eta0 > 0,
                        -Pi <= phi0 < Pi,
                        0 <= alphaR0 <= 2 Pi,
                        xi >= 0,
                        lambda > 0,
                        0 <= Nx <= 1,
                        0 <= Ny <= 1,
                        0 <= Nz <= 1,
                        Nx + Ny + Nz == 1
                    }
                ],
                300,
                energy
            ];

        Print[
            "Saving symbolic energy to: ",
            cacheFile
        ];

        If[
            !DirectoryQ[DirectoryName[cacheFile]],
            CreateDirectory[
                DirectoryName[cacheFile],
                CreateIntermediateDirectories -> True
            ]
        ];

        Put[
            energy,
            cacheFile
        ];

        energy
    ];


BuildMagneticEnergy::dmi =
    "Unknown DMI type `1`. Use \"Bloch\", \"Neel\" or \"Anti\".";