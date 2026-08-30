(* ================================================================ *)
(* GroundState.wl                                                   *)
(*                                                                  *)
(* Ground-state / metastable-equilibrium minimization, storage,    *)
(* loading, and continuation routines.                              *)
(* ================================================================ *)


ClearAll[
    FindMagneticEquilibrium,
    EquilibriumFileName,
    EquilibriumFilePath,
    EquilibriumExistsQ,
    SaveMagneticEquilibrium,
    LoadMagneticEquilibrium,
    EquilibriumRules,
    EquilibriumAnsatz,
    ParameterSweep
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
            Chop[
                genericEnergy /. parameterRules
            ];

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
(* Equilibrium file naming                                        *)
(* --------------------------------------------------------------- *)

(*
    Construct a descriptive file name from the model and
    the most relevant physical parameters.
*)
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

            "_N",
            ToString[$FourierCutoff],

            "_b",
            ToString[N[bValue]],

            "_xi",
            ToString[N[xiValue]],

            "_Nz",
            ToString[N[nzValue]]
        ]
    ];


(* --------------------------------------------------------------- *)
(* Equilibrium file path                                          *)
(* --------------------------------------------------------------- *)

EquilibriumFilePath[
    dmiType_String,
    parameterRules_List
] :=
    FileNameJoin[
        {
            $ProjectRoot,
            "results",
            "equilibria",

            EquilibriumFileName[
                dmiType,
                parameterRules
            ] <> ".txt"
        }
    ];


(* --------------------------------------------------------------- *)
(* Check whether equilibrium already exists                       *)
(* --------------------------------------------------------------- *)

EquilibriumExistsQ[
    dmiType_String,
    parameterRules_List
] :=
    FileExistsQ[
        EquilibriumFilePath[
            dmiType,
            parameterRules
        ]
    ];


(* --------------------------------------------------------------- *)
(* Save equilibrium configuration                                  *)
(* --------------------------------------------------------------- *)

SaveMagneticEquilibrium[
    result_,
    parameterRules_List,
    dmiType_String
] :=
    Module[
        {
            filePath,
            resultDirectory,
            data
        },

        filePath =
            EquilibriumFilePath[
                dmiType,
                parameterRules
            ];

        resultDirectory =
            DirectoryName[
                filePath
            ];

        If[
            !DirectoryQ[resultDirectory],

            CreateDirectory[
                resultDirectory,
                CreateIntermediateDirectories -> True
            ]
        ];

        data =
            <|
                "Energy" ->
                    result[[1]],

                "Solution" ->
                    result[[2]],

                "Parameters" ->
                    parameterRules,

                "DMIType" ->
                    dmiType,

                "FourierCutoff" ->
                    $FourierCutoff,

                "MagneticRepresentation" ->
                    $MagneticRepresentation
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
(* Load saved equilibrium by file name                             *)
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

            Print[
                "Equilibrium file not found: ",
                filePath
            ];

            Return[$Failed];
        ];

        Get[
            filePath
        ]
    ];


(* --------------------------------------------------------------- *)
(* Load saved equilibrium from its physical parameters             *)
(* --------------------------------------------------------------- *)

LoadMagneticEquilibrium[
    dmiType_String,
    parameterRules_List
] :=
    Module[
        {filePath},

        filePath =
            EquilibriumFilePath[
                dmiType,
                parameterRules
            ];

        If[
            !FileExistsQ[filePath],

            Print[
                "Equilibrium file not found: ",
                filePath
            ];

            Return[$Failed];
        ];

        Get[
            filePath
        ]
    ];


(* --------------------------------------------------------------- *)
(* Complete equilibrium rules                                     *)
(* --------------------------------------------------------------- *)

(*
    Combine optimized variables and fixed physical parameters.
    Useful for plotting and real-space reconstruction.
*)
EquilibriumRules[
    equilibrium_Association
] :=
    Join[
        equilibrium["Solution"],
        equilibrium["Parameters"]
    ];


(* --------------------------------------------------------------- *)
(* Convert equilibrium into a minimization ansatz                  *)
(* --------------------------------------------------------------- *)

(*
    Convert replacement rules

        variable -> value

    into FindMinimum initial values

        {variable, value}.
*)
EquilibriumAnsatz[
    equilibrium_Association
] :=
    Module[
        {},

        If[
            equilibrium["FourierCutoff"] =!= $FourierCutoff,

            Print[
                "Saved Fourier cutoff N = ",
                equilibrium["FourierCutoff"],
                " does not match current cutoff N = ",
                $FourierCutoff,
                "."
            ];

            Return[$Failed];
        ];

        If[
            equilibrium["MagneticRepresentation"] =!=
                $MagneticRepresentation,

            Print[
                "Saved magnetic representation ",
                equilibrium["MagneticRepresentation"],
                " does not match current representation ",
                $MagneticRepresentation,
                "."
            ];

            Return[$Failed];
        ];

        List @@@ equilibrium["Solution"]
    ];


(* --------------------------------------------------------------- *)
(* Generic parameter sweep / continuation                          *)
(* --------------------------------------------------------------- *)

(*
    Sweep a physical parameter while following the same
    metastable equilibrium branch.

    At each step:

        1. construct the current parameter set;
        2. check whether the equilibrium already exists;
        3. load it if available;
        4. otherwise minimize and optionally save it;
        5. use the resulting equilibrium as the ansatz for the
           next parameter value.
*)
ParameterSweep[
    genericEnergy_,
    dmiType_String,
    baseParameterRules_List,
    sweepParameter_Symbol,
    parameterValues_List,
    initialValues_List,
    save_: True
] :=
    Module[
        {
            currentAnsatz,
            currentParameters,
            result,
            equilibrium
        },

        currentAnsatz =
            initialValues;

        Table[

            currentParameters =
                Join[
                    DeleteCases[
                        baseParameterRules,
                        Rule[sweepParameter, _]
                    ],

                    {
                        sweepParameter -> parameterValue
                    }
                ];


            Print[
                sweepParameter,
                " = ",
                parameterValue
            ];


            If[
                EquilibriumExistsQ[
                    dmiType,
                    currentParameters
                ],

                (* Existing equilibrium. *)
                Print[
                    "Existing equilibrium found. Loading..."
                ];

                equilibrium =
                    LoadMagneticEquilibrium[
                        dmiType,
                        currentParameters
                    ];

                currentAnsatz =
                    EquilibriumAnsatz[
                        equilibrium
                    ];


                <|
                    "SweepParameter" ->
                        sweepParameter,

                    "ParameterValue" ->
                        parameterValue,

                    "Energy" ->
                        equilibrium["Energy"],

                    "Solution" ->
                        equilibrium["Solution"],

                    "Parameters" ->
                        equilibrium["Parameters"],

                    "Source" ->
                        "Loaded"
                |>,


                (* No saved equilibrium: perform minimization. *)
                Print[
                    "No existing equilibrium found. Minimizing..."
                ];

                result =
                    FindMagneticEquilibrium[
                        genericEnergy,
                        currentParameters,
                        currentAnsatz
                    ];

                currentAnsatz =
                    List @@@ result[[2]];

                If[
                    save,

                    SaveMagneticEquilibrium[
                        result,
                        currentParameters,
                        dmiType
                    ]
                ];


                <|
                    "SweepParameter" ->
                        sweepParameter,

                    "ParameterValue" ->
                        parameterValue,

                    "Energy" ->
                        result[[1]],

                    "Solution" ->
                        result[[2]],

                    "Parameters" ->
                        currentParameters,

                    "Source" ->
                        "Calculated"
                |>
            ],

            {parameterValue, parameterValues}
        ]
    ];