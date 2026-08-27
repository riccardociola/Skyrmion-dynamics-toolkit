(* ================================================================ *)
(* GroundState.wl                                                   *)
(*                                                                  *)
(* Ground-state / metastable-equilibrium minimization routines.      *)
(* ================================================================ *)


ClearAll[
    FindMagneticEquilibrium,
    SaveMagneticEquilibrium,
    LoadMagneticEquilibrium
];


(* --------------------------------------------------------------- *)
(* Magnetic equilibrium minimization                              *)
(* --------------------------------------------------------------- *)

FindMagneticEquilibrium[
    genericEnergy_,
    parameterRules_List,
    initialValues_List
] :=
    Module[
        {specificEnergy},

        (* Insert physical parameters into the cached symbolic energy. *)
        specificEnergy =
            Chop[genericEnergy /. parameterRules];

        (* Minimize from the supplied initial ansatz. *)
        FindMinimum[
            specificEnergy,
            initialValues,
            Method -> "Newton",
            AccuracyGoal -> 12,
            PrecisionGoal -> 10,
            WorkingPrecision -> 40,
            MaxIterations -> Infinity,
            StepMonitor :>
                Print[
                    {a0, eta0, phi0, alphaR0} /. parameterRules
                ]
        ]
    ];


(* --------------------------------------------------------------- *)
(* Save equilibrium configuration                                  *)
(* --------------------------------------------------------------- *)


EquilibriumFileName[
    dmiType_String,
    parameterRules_List
] :=
    Module[
        {
            bValue,
            xiValue,
            nzValue
        },

        bValue =
            b /. parameterRules;

        xiValue =
            xi /. parameterRules;

        nzValue =
            Nz /. parameterRules;

        StringJoin[
            dmiType,
            "_N", ToString[$FourierCutoff],

            "_b", ToString[N[bValue]],

            "_xi", ToString[N[xiValue]],

            "_Nz", ToString[N[nzValue]]
        ]
    ];


SaveMagneticEquilibrium[
    result_,
    parameterRules_List,
    dmiType_String
] :=
    Module[
        {
            resultDirectory,
            fileName,
            filePath,
            data
        },

        resultDirectory =
            FileNameJoin[
                {
                    $ProjectRoot,
                    "results",
                    "equilibria"
                }
            ];

        If[
            !DirectoryQ[resultDirectory],
            CreateDirectory[
                resultDirectory,
                CreateIntermediateDirectories -> True
            ]
        ];

        fileName =
            EquilibriumFileName[
                dmiType,
                parameterRules
            ];

        filePath =
            FileNameJoin[
                {
                    resultDirectory,
                    fileName <> ".txt"
                }
            ];

        data =
            <|
                "Energy" -> result[[1]],
                "Solution" -> result[[2]],
                "Parameters" -> parameterRules,
                "DMIType" -> dmiType,
                "FourierCutoff" -> $FourierCutoff,
                "MagneticRepresentation" -> $MagneticRepresentation
            |>;

        Put[
            data,
            filePath
        ];

        Print[
            "Equilibrium configuration saved to: ",
            filePath
        ];

        filePath
    ];


(* --------------------------------------------------------------- *)
(* Load saved equilibrium                                          *)
(* --------------------------------------------------------------- *)

LoadMagneticEquilibrium[
    fileName_String
] :=
    Module[
        {filePath},

        filePath =
            FileNameJoin[
                {
                    $ProjectRoot,
                    "results",
                    "equilibria",
                    fileName
                }
            ];

        If[
            !FileExistsQ[filePath],
            Print["Equilibrium file not found: ", filePath];
            Return[$Failed];
        ];

        Get[filePath]
    ];