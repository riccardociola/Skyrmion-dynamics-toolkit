(* ================================================================ *)
(* GroundState.wl                                                   *)
(*                                                                  *)
(* Magnetic equilibrium minimization, storage and continuation.    *)
(* Supports SixFold / General Fourier bases and optional cubic MCA. *)
(*                                                                  *)
(* Dependencies: FourierBasis.wl, BuildEnergy.wl,                    *)
(*               EnergyAnisotropy.wl for cubic MCA.                 *)
(* ================================================================ *)

ClearAll[
    GroundStateModelSettings,
    GroundStateInitialValues,
    GroundStateParameterLabel,
    FindMagneticEquilibrium,
    EquilibriumFileName,
    EquilibriumFilePath,
    EquilibriumExistsQ,
    SaveMagneticEquilibrium,
    LoadMagneticEquilibrium,
    EquilibriumRules,
    EquilibriumAnsatz,
    ParameterSweep,
    GroundStateHelp
];

(* Model options match BuildMagneticEnergy. Automatic representation   *)
(* and cutoff use the active Fourier configuration.                  *)
$GroundStateModelOptions = {
    "FourierCutoff" -> Automatic,
    "MagneticRepresentation" -> Automatic,
    "CubicAnisotropy" -> Automatic,
    "FieldDirection" -> {0, 0, 1}
};

Options[EquilibriumFileName] = $GroundStateModelOptions;
Options[EquilibriumFilePath] = $GroundStateModelOptions;
Options[EquilibriumExistsQ] = $GroundStateModelOptions;
Options[SaveMagneticEquilibrium] = $GroundStateModelOptions;
Options[LoadMagneticEquilibrium] = $GroundStateModelOptions;

Options[FindMagneticEquilibrium] = Join[
    $GroundStateModelOptions,
    {"Save" -> False, "DMIType" -> Automatic, "PrintSteps" -> False}
];

Options[EquilibriumAnsatz] = Join[
    $GroundStateModelOptions,
    {"AllowRepresentationChange" -> False}
];

Options[ParameterSweep] = $GroundStateModelOptions;

GroundStateModelSettings::cutoff =
    "Invalid Fourier cutoff `1`. Build an energy first or specify a nonnegative FourierCutoff.";
GroundStateModelSettings::rep =
    "Invalid magnetic representation `1`. Expected \"SixFold\" or \"General\".";
GroundStateModelSettings::field =
    "Invalid cubic field direction `1`. Check that EnergyAnisotropy.wl is loaded.";
GroundStateModelSettings::builder =
    "BuildEnergy.wl must be loaded to construct model-specific equilibrium filenames.";
FindMagneticEquilibrium::kappa =
    "Cubic MCA was requested but parameterRules does not contain kappaCub -> value.";
FindMagneticEquilibrium::energy =
    "CubicAnisotropy option does not agree with the symbolic energy: presence of kappaCub = `1`, option = `2`.";
FindMagneticEquilibrium::vars =
    "Initial values must specify a0, eta0 and phi0; missing `1`.";
FindMagneticEquilibrium::ansatz =
    "Unable to prepare the minimization ansatz.";
FindMagneticEquilibrium::save =
    "To save automatically, specify \"DMIType\" -> \"Bloch\", \"Neel\", or \"Anti\".";
EquilibriumAnsatz::cutoff =
    "Saved cutoff `1` differs from the current cutoff `2`.";
EquilibriumAnsatz::rep =
    "Saved representation `1` differs from the requested representation `2`. Set \"AllowRepresentationChange\" -> True to expand a SixFold seed into a General one.";
EquilibriumAnsatz::conversion =
    "Automatic conversion from `1` to `2` is not supported.";
SaveMagneticEquilibrium::result =
    "Expected a successful FindMinimum result {energy, solutionRules}.";
SaveMagneticEquilibrium::rotation =
    "Cubic MCA equilibrium must contain alphaR0 among the optimized solution variables.";
ParameterSweep::ansatz =
    "Could not construct a compatible continuation ansatz for `1` = `2`.";

(* --------------------------------------------------------------- *)
(* Model metadata                                                  *)
(* --------------------------------------------------------------- *)

(* Resolve the model once and pass the same choices to naming,      *)
(* saving, loading, and minimization.                                *)
GroundStateModelSettings[
    parameterRules_List,
    cutoffInput_: Automatic,
    representationInput_: Automatic,
    cubicInput_: Automatic,
    fieldDirectionInput_: {0, 0, 1}
] := Module[
    {cutoff, representation, useMCA, direction},

    cutoff = Replace[
        cutoffInput,
        Automatic :> If[ValueQ[$FourierCutoff], $FourierCutoff, $Failed]
    ];
    representation = Replace[
        representationInput,
        Automatic :> If[ValueQ[$MagneticRepresentation],
            $MagneticRepresentation, "SixFold"]
    ];
    useMCA = Replace[
        cubicInput,
        Automatic :> !FreeQ[parameterRules, HoldPattern[kappaCub -> _]]
    ];

    If[!IntegerQ[cutoff] || cutoff < 0,
        Message[GroundStateModelSettings::cutoff, cutoff];
        Return[$Failed]
    ];
    If[!MemberQ[{"SixFold", "General"}, representation],
        Message[GroundStateModelSettings::rep, representation];
        Return[$Failed]
    ];

    direction = Missing["NotApplicable"];
    If[TrueQ[useMCA],
        If[!ListQ[fieldDirectionInput] || Length[fieldDirectionInput] != 3 ||
            Length[DownValues[CanonicalCubicFieldDirection]] == 0,
            Message[GroundStateModelSettings::field, fieldDirectionInput];
            Return[$Failed]
        ];
        direction = CanonicalCubicFieldDirection[fieldDirectionInput];
        If[direction === $Failed,
            Message[GroundStateModelSettings::field, fieldDirectionInput];
            Return[$Failed]
        ];
    ];

    <|
        "FourierCutoff" -> cutoff,
        "MagneticRepresentation" -> representation,
        "CubicAnisotropy" -> TrueQ[useMCA],
        "FieldDirectionCrystal" -> direction
    |>
];

(* --------------------------------------------------------------- *)
(* Minimization ansatz                                             *)
(* --------------------------------------------------------------- *)

(* Include ALL independent Fourier variables of the chosen basis.  *)
(* Missing coefficients are seeded with zero (useful when promoting *)
(* a SixFold equilibrium to the General representation).            *)
(* For isotropic energies alphaR0 is fixed by a parameter rule;     *)
(* for MCA alphaR0 is optimized, with an initial seed of zero.       *)
GroundStateInitialValues[
    initialValues_List,
    settings_Association,
    parameterRules_List
] := Module[
    {seedRules, seed, fixedValues, useMCA, variables, missingGeometry,
     alphaSeed, prepared, parameters},

    seedRules = Cases[initialValues, {_, _}];
    seed = Association[Rule @@@ seedRules];
    fixedValues = Association[parameterRules];
    useMCA = settings["CubicAnisotropy"];

    missingGeometry = Select[{a0, eta0, phi0}, !KeyExistsQ[seed, #] &];
    If[missingGeometry =!= {},
        Message[FindMagneticEquilibrium::vars, missingGeometry];
        Return[$Failed]
    ];

    alphaSeed = If[useMCA,
        Lookup[seed, alphaR0, Lookup[fixedValues, alphaR0, 0]],
        Lookup[fixedValues, alphaR0, Lookup[seed, alphaR0, 0]]
    ];

    parameters = DeleteCases[parameterRules, HoldPattern[alphaR0 -> _]];
    If[!useMCA,
        parameters = Append[parameters, alphaR0 -> alphaSeed]
    ];

    variables = Join[
        MagneticVariables[
            settings["FourierCutoff"], settings["MagneticRepresentation"]
        ],
        {a0, eta0, phi0},
        If[useMCA, {alphaR0}, {}]
    ];

    prepared = Table[
        {variable, If[variable === alphaR0, alphaSeed,
            Lookup[seed, variable, 0]]},
        {variable, variables}
    ];

    <|"InitialValues" -> prepared, "Parameters" -> parameters|>
];

(* --------------------------------------------------------------- *)
(* Magnetic equilibrium minimization                              *)
(* --------------------------------------------------------------- *)

(* genericEnergy comes from BuildMagneticEnergy with matching       *)
(* representation, MCA flag and field direction. kappaCub is        *)
(* substituted numerically here, not in BuildEnergy.wl.            *)
FindMagneticEquilibrium[
    genericEnergy_,
    parameterRules_List,
    initialValues_List,
    OptionsPattern[]
] := Module[
    {useMCA, hasCubicTerm, settings, prepared, effectiveParameters,
     startingValues, specificEnergy, result, dmiType, printSteps},

    hasCubicTerm = !FreeQ[genericEnergy, kappaCub];
    useMCA = Replace[OptionValue["CubicAnisotropy"],
        Automatic -> hasCubicTerm];

    If[TrueQ[useMCA] =!= hasCubicTerm,
        Message[FindMagneticEquilibrium::energy, hasCubicTerm, useMCA];
        Return[$Failed]
    ];
    If[TrueQ[useMCA] &&
       FreeQ[parameterRules, HoldPattern[kappaCub -> _]],
        Message[FindMagneticEquilibrium::kappa];
        Return[$Failed]
    ];

    settings = GroundStateModelSettings[
        parameterRules,
        OptionValue["FourierCutoff"],
        OptionValue["MagneticRepresentation"],
        useMCA,
        OptionValue["FieldDirection"]
    ];
    If[settings === $Failed, Return[$Failed]];

    prepared = GroundStateInitialValues[initialValues, settings, parameterRules];
    If[prepared === $Failed,
        Message[FindMagneticEquilibrium::ansatz];
        Return[$Failed]
    ];

    effectiveParameters = prepared["Parameters"];
    startingValues = prepared["InitialValues"];
    specificEnergy = Chop[genericEnergy /. effectiveParameters];
    printSteps = TrueQ[OptionValue["PrintSteps"]];

    result = FindMinimum[
        specificEnergy,
        startingValues,
        Method -> "Newton",
        AccuracyGoal -> 12,
        PrecisionGoal -> 10,
        WorkingPrecision -> 40,
        MaxIterations -> Infinity,
        StepMonitor :> If[printSteps,
            Print[{a0, eta0, phi0, alphaR0} /. effectiveParameters]
        ]
    ];

    If[TrueQ[OptionValue["Save"]] && MatchQ[result, {_, _List}],
        dmiType = OptionValue["DMIType"];
        If[dmiType === Automatic,
            Message[FindMagneticEquilibrium::save],
            SaveMagneticEquilibrium[
                result, effectiveParameters, dmiType,
                "FourierCutoff" -> settings["FourierCutoff"],
                "MagneticRepresentation" -> settings["MagneticRepresentation"],
                "CubicAnisotropy" -> useMCA,
                "FieldDirection" -> OptionValue["FieldDirection"]
            ]
        ]
    ];
    result
];

(* --------------------------------------------------------------- *)
(* Equilibrium file naming                                        *)
(* --------------------------------------------------------------- *)

(* Human-readable model prefix matches MagneticEnergyCacheFile.     *)
(* A short hash of all fixed physical parameters prevents different *)
(* lambda or demagnetization settings sharing an equilibrium file.  *)
GroundStateParameterLabel[parameterRules_List, parameter_, label_String] :=
    Module[{value, printable},
        value = parameter /. parameterRules;
        If[value === parameter, Return[""]];
        printable = ToString[InputForm[N[value, 7]]];
        printable = StringReplace[printable, {
            RegularExpression["`[0-9.]*"] -> "",
            "*^" -> "e", " " -> "", "-" -> "m", "/" -> "_"
        }];
        "_" <> label <> printable
    ];

EquilibriumFileName[
    dmiType_String,
    parameterRules_List,
    OptionsPattern[]
] := Module[
    {settings, modelFile, modelPrefix, fixed, signature, digest,
     cutoff, rep, useMCA, direction},

    settings = GroundStateModelSettings[
        parameterRules,
        OptionValue["FourierCutoff"],
        OptionValue["MagneticRepresentation"],
        OptionValue["CubicAnisotropy"],
        OptionValue["FieldDirection"]
    ];
    If[settings === $Failed, Return[$Failed]];

    If[Length[DownValues[MagneticEnergyCacheFile]] == 0,
        Message[GroundStateModelSettings::builder];
        Return[$Failed]
    ];

    cutoff = settings["FourierCutoff"];
    rep = settings["MagneticRepresentation"];
    useMCA = settings["CubicAnisotropy"];
    direction = OptionValue["FieldDirection"];
    modelFile = MagneticEnergyCacheFile[
        cutoff, dmiType,
        "MagneticRepresentation" -> rep,
        "CubicAnisotropy" -> useMCA,
        "FieldDirection" -> direction
    ];
    If[modelFile === $Failed, Return[$Failed]];
    modelPrefix = FileBaseName[modelFile];

    fixed = SortBy[parameterRules, ToString[InputForm[First[#]]] &];
    signature = ToString[InputForm[fixed]];
    digest = StringTake[IntegerString[Hash[signature, "SHA256"], 16, 64], 10];

    StringJoin[
        modelPrefix,
        GroundStateParameterLabel[parameterRules, b, "b"],
        GroundStateParameterLabel[parameterRules, xi, "xi"],
        GroundStateParameterLabel[parameterRules, Nz, "Nz"],
        If[useMCA, GroundStateParameterLabel[parameterRules, kappaCub, "K"], ""],
        "_id", digest
    ]
];

EquilibriumFilePath[
    dmiType_String,
    parameterRules_List,
    OptionsPattern[]
] := Module[{name},
    name = EquilibriumFileName[
        dmiType, parameterRules,
        "FourierCutoff" -> OptionValue["FourierCutoff"],
        "MagneticRepresentation" -> OptionValue["MagneticRepresentation"],
        "CubicAnisotropy" -> OptionValue["CubicAnisotropy"],
        "FieldDirection" -> OptionValue["FieldDirection"]
    ];
    If[name === $Failed, Return[$Failed]];
    FileNameJoin[{$ProjectRoot, "results", "equilibria", name <> ".txt"}]
];

EquilibriumExistsQ[
    dmiType_String,
    parameterRules_List,
    OptionsPattern[]
] := Module[{path},
    path = EquilibriumFilePath[
        dmiType, parameterRules,
        "FourierCutoff" -> OptionValue["FourierCutoff"],
        "MagneticRepresentation" -> OptionValue["MagneticRepresentation"],
        "CubicAnisotropy" -> OptionValue["CubicAnisotropy"],
        "FieldDirection" -> OptionValue["FieldDirection"]
    ];
    path =!= $Failed && FileExistsQ[path]
];

(* --------------------------------------------------------------- *)
(* Save / load equilibria                                          *)
(* --------------------------------------------------------------- *)

SaveMagneticEquilibrium[
    result_,
    parameterRules_List,
    dmiType_String,
    OptionsPattern[]
] := Module[{settings, path, data},
    If[!MatchQ[result, {_, _List}],
        Message[SaveMagneticEquilibrium::result];
        Return[$Failed]
    ];

    settings = GroundStateModelSettings[
        parameterRules,
        OptionValue["FourierCutoff"],
        OptionValue["MagneticRepresentation"],
        OptionValue["CubicAnisotropy"],
        OptionValue["FieldDirection"]
    ];
    If[settings === $Failed, Return[$Failed]];

    effectiveParameters = If[settings["CubicAnisotropy"],
        DeleteCases[parameterRules, HoldPattern[alphaR0 -> _]],
        parameterRules
    ];
    If[settings["CubicAnisotropy"] &&
       FreeQ[effectiveParameters, HoldPattern[kappaCub -> _]],
        Message[FindMagneticEquilibrium::kappa];
        Return[$Failed]
    ];
    If[!settings["CubicAnisotropy"] &&
       FreeQ[parameterRules, HoldPattern[alphaR0 -> _]],
        effectiveParameters = Append[effectiveParameters, alphaR0 -> 0]
    ];
    If[settings["CubicAnisotropy"] &&
       FreeQ[result[[2]], HoldPattern[alphaR0 -> _]],
        Message[SaveMagneticEquilibrium::rotation];
        Return[$Failed]
    ];

    path = EquilibriumFilePath[
        dmiType, effectiveParameters,
        "FourierCutoff" -> settings["FourierCutoff"],
        "MagneticRepresentation" -> settings["MagneticRepresentation"],
        "CubicAnisotropy" -> settings["CubicAnisotropy"],
        "FieldDirection" -> OptionValue["FieldDirection"]
    ];
    If[path === $Failed, Return[$Failed]];
    If[!DirectoryQ[DirectoryName[path]],
        CreateDirectory[DirectoryName[path],
            CreateIntermediateDirectories -> True]
    ];

    data = <|
        "Energy" -> result[[1]],
        "Solution" -> result[[2]],
        "Parameters" -> effectiveParameters,
        "DMIType" -> dmiType,
        "FourierCutoff" -> settings["FourierCutoff"],
        "MagneticRepresentation" -> settings["MagneticRepresentation"],
        "CubicAnisotropy" -> settings["CubicAnisotropy"],
        "FieldDirectionCrystal" -> settings["FieldDirectionCrystal"]
    |>;
    Put[data, path];
    Print["Equilibrium configuration saved to: ", path];
    path
];

(* Loading by filename remains compatible with older .txt files.   *)
LoadMagneticEquilibrium[fileName_String] := Module[{path},
    path = If[FileExistsQ[fileName], fileName,
        FileNameJoin[{$ProjectRoot, "results", "equilibria", fileName}]
    ];
    If[!FileExistsQ[path],
        Print["Equilibrium file not found: ", path];
        Return[$Failed]
    ];
    Get[path]
];

LoadMagneticEquilibrium[
    dmiType_String,
    parameterRules_List,
    OptionsPattern[]
] := Module[{path},
    path = EquilibriumFilePath[
        dmiType, parameterRules,
        "FourierCutoff" -> OptionValue["FourierCutoff"],
        "MagneticRepresentation" -> OptionValue["MagneticRepresentation"],
        "CubicAnisotropy" -> OptionValue["CubicAnisotropy"],
        "FieldDirection" -> OptionValue["FieldDirection"]
    ];
    If[path === $Failed, Return[$Failed]];
    LoadMagneticEquilibrium[path]
];

EquilibriumRules[equilibrium_Association] :=
    Join[equilibrium["Solution"], equilibrium["Parameters"]];

(* --------------------------------------------------------------- *)
(* Restart / continuation ansatz                                   *)
(* --------------------------------------------------------------- *)

(* An existing SixFold equilibrium can initialize a General one:   *)
(* shared Fourier coefficients are reused; extra ones start at 0.  *)
(* alphaR0 becomes an optimization variable if cubic MCA is on.     *)
EquilibriumAnsatz[
    equilibrium_Association,
    OptionsPattern[]
] := Module[{cutoff, rep, oldRep, oldCutoff, useMCA, settings,
             seedRules, prepared},

    oldCutoff = equilibrium["FourierCutoff"];
    oldRep = equilibrium["MagneticRepresentation"];
    cutoff = Replace[OptionValue["FourierCutoff"],
        Automatic :> If[ValueQ[$FourierCutoff], $FourierCutoff, oldCutoff]];
    rep = Replace[OptionValue["MagneticRepresentation"],
        Automatic :> If[ValueQ[$MagneticRepresentation], $MagneticRepresentation, oldRep]];

    If[oldCutoff =!= cutoff,
        Message[EquilibriumAnsatz::cutoff, oldCutoff, cutoff];
        Return[$Failed]
    ];
    If[oldRep =!= rep,
        If[!TrueQ[OptionValue["AllowRepresentationChange"]],
            Message[EquilibriumAnsatz::rep, oldRep, rep];
            Return[$Failed]
        ];
        If[!(oldRep === "SixFold" && rep === "General"),
            Message[EquilibriumAnsatz::conversion, oldRep, rep];
            Return[$Failed]
        ]
    ];

    useMCA = Replace[OptionValue["CubicAnisotropy"],
        Automatic :> Lookup[equilibrium, "CubicAnisotropy", False]];
    settings = GroundStateModelSettings[
        equilibrium["Parameters"], cutoff, rep, useMCA,
        OptionValue["FieldDirection"]
    ];
    If[settings === $Failed, Return[$Failed]];

    seedRules = List @@@ equilibrium["Solution"];
    prepared = GroundStateInitialValues[
        seedRules, settings, equilibrium["Parameters"]
    ];
    If[prepared === $Failed, Return[$Failed]];
    prepared["InitialValues"]
];

(* --------------------------------------------------------------- *)
(* Generic parameter sweep / continuation                          *)
(* --------------------------------------------------------------- *)

ParameterSweep[
    genericEnergy_,
    dmiType_String,
    baseParameterRules_List,
    sweepParameter_Symbol,
    parameterValues_List,
    initialValues_List,
    save_: True,
    OptionsPattern[]
] := Module[
    {settings, useMCA, ansatz, currentParameters, currentPrepared,
     effectiveParameters, path, equilibrium, result, source, rows = {}},

    useMCA = Replace[OptionValue["CubicAnisotropy"],
        Automatic -> !FreeQ[genericEnergy, kappaCub]];
    settings = GroundStateModelSettings[
        baseParameterRules,
        OptionValue["FourierCutoff"],
        OptionValue["MagneticRepresentation"],
        useMCA,
        OptionValue["FieldDirection"]
    ];
    If[settings === $Failed, Return[$Failed]];

    ansatz = initialValues;
    Do[
        currentParameters = Join[
            DeleteCases[baseParameterRules, HoldPattern[Rule[sweepParameter, _]]],
            {sweepParameter -> value}
        ];
        currentPrepared = GroundStateInitialValues[
            ansatz, settings, currentParameters
        ];
        If[currentPrepared === $Failed,
            Message[ParameterSweep::ansatz, sweepParameter, value];
            Return[$Failed]
        ];
        effectiveParameters = currentPrepared["Parameters"];

        path = EquilibriumFilePath[
            dmiType, effectiveParameters,
            "FourierCutoff" -> settings["FourierCutoff"],
            "MagneticRepresentation" -> settings["MagneticRepresentation"],
            "CubicAnisotropy" -> useMCA,
            "FieldDirection" -> OptionValue["FieldDirection"]
        ];
        If[path === $Failed, Return[$Failed]];

        If[FileExistsQ[path],
            Print["Loading existing equilibrium: ", path];
            equilibrium = LoadMagneticEquilibrium[path];
            If[equilibrium === $Failed, Return[$Failed]];
            result = {equilibrium["Energy"], equilibrium["Solution"]};
            source = "Loaded";
            ,
            Print[sweepParameter, " = ", value, ": minimizing..."];
            result = FindMagneticEquilibrium[
                genericEnergy, effectiveParameters,
                currentPrepared["InitialValues"],
                "FourierCutoff" -> settings["FourierCutoff"],
                "MagneticRepresentation" -> settings["MagneticRepresentation"],
                "CubicAnisotropy" -> useMCA,
                "FieldDirection" -> OptionValue["FieldDirection"]
            ];
            If[!MatchQ[result, {_, _List}], Return[$Failed]];
            source = "Calculated";
            If[TrueQ[save],
                SaveMagneticEquilibrium[
                    result, effectiveParameters, dmiType,
                    "FourierCutoff" -> settings["FourierCutoff"],
                    "MagneticRepresentation" -> settings["MagneticRepresentation"],
                    "CubicAnisotropy" -> useMCA,
                    "FieldDirection" -> OptionValue["FieldDirection"]
                ]
            ]
        ];

        ansatz = List @@@ result[[2]];
        AppendTo[rows, <|
            "SweepParameter" -> sweepParameter,
            "ParameterValue" -> value,
            "Energy" -> result[[1]],
            "Solution" -> result[[2]],
            "Parameters" -> effectiveParameters,
            "Source" -> source
        |>];,
        {value, parameterValues}
    ];
    rows
];

(* --------------------------------------------------------------- *)
(* Help                                                            *)
(* --------------------------------------------------------------- *)

GroundStateHelp[] := Column[{
    Style["GroundState.wl", Bold, 16],
    "Magnetic equilibrium minimization, saving, loading and continuation; supports cubic MCA and either Fourier representation.",
    Grid[{
        {"Function", "Description"},
        {"FindMagneticEquilibrium", "Minimize a magnetic energy; optimize alphaR0 when MCA is active."},
        {"GroundStateInitialValues", "Complete the Fourier ansatz and handle the alphaR0 variable."},
        {"EquilibriumFileName / EquilibriumFilePath", "Distinct names for model, MCA orientation and fixed parameters."},
        {"SaveMagneticEquilibrium / LoadMagneticEquilibrium", "Store and recover equilibria with model metadata."},
        {"EquilibriumRules / EquilibriumAnsatz", "Read solutions or prepare a compatible restart ansatz."},
        {"ParameterSweep", "Continue an equilibrium branch while varying a physical parameter."}
    }, Frame -> All, Alignment -> Left],
    "Model options: FourierCutoff, MagneticRepresentation, CubicAnisotropy, FieldDirection."
}];

Print[GroundStateHelp[]];