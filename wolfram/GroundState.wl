(* ================================================================ *)
(* GroundState.wl                                                   *)
(* Magnetic equilibrium minimization, storage and continuation.    *)
(* Supports MCA, SixFold/General bases, and optional rigid lattice. *)
(* Dependencies: FourierBasis.wl, BuildEnergy.wl, EnergyAnisotropy.wl *)
(* ================================================================ *)

ClearAll[
    GroundStateFieldDirectionTag, GroundStateModelSettings,
    GroundStateInitialValues, GroundStateCompletedSolution,
    GroundStateParameterLabel, GroundStateNzString,
    GroundStateParametersMatchQ, GroundStateValuesEqualQ,
    GroundStateLocateSavedEquilibrium,
    GroundStateMatchesRequestQ, FindMagneticEquilibrium,
    EquilibriumFileName, EquilibriumFilePath, EquilibriumExistsQ,
    SaveMagneticEquilibrium, LoadMagneticEquilibrium,
    EquilibriumRules, EquilibriumAnsatz, ParameterSweep, GroundStateHelp
];

$GroundStateModelOptions = {
    "FourierCutoff" -> Automatic,
    "MagneticRepresentation" -> Automatic,
    "CubicAnisotropy" -> Automatic,
    "FieldDirection" -> {0, 0, 1},
    "RigidLattice" -> False
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
    {"FourierCutoff" -> Automatic,
     "MagneticRepresentation" -> Automatic,
     "CubicAnisotropy" -> Automatic,
     "FieldDirection" -> {0, 0, 1},
     "RigidLattice" -> Automatic},
    {"AllowRepresentationChange" -> False}
];
Options[ParameterSweep] = $GroundStateModelOptions;

GroundStateModelSettings::cutoff = "Invalid Fourier cutoff `1`.";
GroundStateModelSettings::rep = "Invalid magnetic representation `1`.";
GroundStateModelSettings::field = "Invalid field direction `1`.";
GroundStateModelSettings::rigid = "RigidLattice must be True or False (received `1`).";
EquilibriumFileName::params = "Equilibrium filenames require b, xi, Nz and lambda; missing `1`.";
FindMagneticEquilibrium::kappa = "Cubic MCA requested but kappaCub -> value is missing.";
FindMagneticEquilibrium::energy = "CubicAnisotropy = `2` disagrees with the energy (kappaCub present = `1`).";
FindMagneticEquilibrium::vars = "Initial geometry values are missing: `1`.";
FindMagneticEquilibrium::ansatz = "Could not prepare minimization variables.";
FindMagneticEquilibrium::save = "Set DMIType to Bloch, Neel or Anti to save automatically.";
EquilibriumAnsatz::cutoff = "Saved Fourier cutoff `1` differs from requested cutoff `2`.";
EquilibriumAnsatz::rep = "Saved representation `1` differs from requested `2`. Set AllowRepresentationChange -> True to expand a SixFold seed into General.";
EquilibriumAnsatz::conversion = "Unsupported representation change `1` -> `2`.";
SaveMagneticEquilibrium::result = "Expected {energy, solutionRules} with valid equilibrium geometry.";
SaveMagneticEquilibrium::rotation = "Cubic MCA requires alphaR0 in Solution.";
SaveMagneticEquilibrium::collision = "File `1` already contains different model metadata or physical parameters. Refusing to overwrite.";
LoadMagneticEquilibrium::metadata = "Saved equilibrium in `1` does not match the requested parameters/model.";
ParameterSweep::ansatz = "Could not prepare the continuation ansatz for `1` = `2`.";

(* --------------------------------------------------------------- *)
(* Model metadata and crystallographic field label                  *)
(* --------------------------------------------------------------- *)

(* Use the user's raw indices: H100 and H001 are different filenames. *)
(* Minus signs and decimal points become filename-safe text.         *)
GroundStateFieldDirectionTag[direction_List] /; Length[direction] == 3 :=
    Module[{parts},
        parts = ToString[InputForm[#]] & /@ direction;
        parts = StringReplace[parts, {
            "-" -> "m", "." -> "p", " " -> "", "/" -> "d",
            "[" -> "", "]" -> "", "(" -> "", ")" -> "",
            "*" -> "x", "^" -> "pow", "`" -> ""
        }];
        "H" <> If[AllTrue[parts, StringMatchQ[#, DigitCharacter] &],
            StringJoin[parts], StringRiffle[parts, "_"]]
    ];
GroundStateFieldDirectionTag[_] := $Failed;

(* The orientation is canonicalized ONLY for the physical crystal frame. *)
(* Filename and title refer to FieldDirectionInput, not its canonical copy. *)
GroundStateModelSettings[
    parameterRules_List, cutoffInput_: Automatic,
    representationInput_: Automatic, cubicInput_: Automatic,
    fieldDirectionInput_: {0, 0, 1}, rigidInput_: False
] := Module[{cutoff, representation, useMCA, direction, directionTag},
    cutoff = Replace[cutoffInput,
        Automatic :> If[ValueQ[$FourierCutoff], $FourierCutoff, $Failed]];
    representation = Replace[representationInput,
        Automatic :> If[ValueQ[$MagneticRepresentation], $MagneticRepresentation, "SixFold"]];
    useMCA = Replace[cubicInput,
        Automatic :> !FreeQ[parameterRules, HoldPattern[kappaCub -> _]]];

    If[!IntegerQ[cutoff] || cutoff < 0,
        Message[GroundStateModelSettings::cutoff, cutoff]; Return[$Failed]];
    If[!MemberQ[{"SixFold", "General"}, representation],
        Message[GroundStateModelSettings::rep, representation]; Return[$Failed]];
    If[!MemberQ[{True, False}, rigidInput],
        Message[GroundStateModelSettings::rigid, rigidInput]; Return[$Failed]];

    direction = Missing["NotApplicable"];
    directionTag = Missing["NotApplicable"];
    If[TrueQ[useMCA],
        If[!ListQ[fieldDirectionInput] || Length[fieldDirectionInput] != 3 ||
           Length[DownValues[CanonicalCubicFieldDirection]] == 0,
            Message[GroundStateModelSettings::field, fieldDirectionInput];
            Return[$Failed]
        ];
        direction = CanonicalCubicFieldDirection[fieldDirectionInput];
        directionTag = GroundStateFieldDirectionTag[fieldDirectionInput];
        If[direction === $Failed || directionTag === $Failed,
            Message[GroundStateModelSettings::field, fieldDirectionInput];
            Return[$Failed]
        ];
    ];

    <|
        "FourierCutoff" -> cutoff,
        "MagneticRepresentation" -> representation,
        "CubicAnisotropy" -> TrueQ[useMCA],
        "RigidLattice" -> TrueQ[rigidInput],
        "FieldDirectionInput" -> If[TrueQ[useMCA], fieldDirectionInput, Missing["NotApplicable"]],
        "FieldDirectionCrystal" -> direction,
        "FieldDirectionTag" -> directionTag
    |>
];

(* --------------------------------------------------------------- *)
(* Minimization ansatz                                              *)
(* --------------------------------------------------------------- *)

(* eta0 and phi0 are optimized normally, but fixed to 1 and 0 for  *)
(* RigidLattice. They ALWAYS belong to the saved Solution, never   *)
(* Parameters. Isotropic configurations ALWAYS fix alphaR0 -> 0.  *)
GroundStateInitialValues[
    initialValues_List, settings_Association, parameterRules_List
] := Module[{seed, useMCA, rigid, variables, missingGeometry,
             effectiveParameters, alphaSeed, prepared},
    seed = Association[Rule @@@ Cases[initialValues, {_, _}]];
    useMCA = settings["CubicAnisotropy"];
    rigid = settings["RigidLattice"];

    missingGeometry = Select[
        If[rigid, {a0}, {a0, eta0, phi0}], !KeyExistsQ[seed, #] &];
    If[missingGeometry =!= {},
        Message[FindMagneticEquilibrium::vars, missingGeometry]; Return[$Failed]];

    (* Remove geometry and rotation rules supplied as fixed parameters. *)
    effectiveParameters = Select[parameterRules,
        !MemberQ[{eta0, phi0, alphaR0}, First[#]] &];
    If[!useMCA, AppendTo[effectiveParameters, alphaR0 -> 0]];

    alphaSeed = If[useMCA,
        Lookup[seed, alphaR0, Lookup[Association[parameterRules], alphaR0, 0]], 0];

    variables = Join[
        MagneticVariables[settings["FourierCutoff"], settings["MagneticRepresentation"]],
        {a0}, If[rigid, {}, {eta0, phi0}],
        If[useMCA, {alphaR0}, {}]
    ];
    prepared = Table[
        {var, If[var === alphaR0, alphaSeed, Lookup[seed, var, 0]]},
        {var, variables}
    ];
    <|
        "InitialValues" -> prepared,
        "Parameters" -> effectiveParameters,
        "FixedGeometryRules" -> If[rigid, {eta0 -> 1, phi0 -> 0}, {}]
    |>
];

(* Complete the solution in a standard format for both lattice types. *)
GroundStateCompletedSolution[solution_List, rigid_, useMCA_] :=
    Module[{rules},
        rules = Select[solution,
            !MemberQ[If[rigid, {eta0, phi0}, {}], First[#]] &];
        If[rigid, rules = Join[rules, {eta0 -> 1, phi0 -> 0}]];
        If[!useMCA,
            rules = Select[rules, First[#] =!= alphaR0 &]];
        If[FreeQ[rules, HoldPattern[eta0 -> _]] ||
           FreeQ[rules, HoldPattern[phi0 -> _]], Return[$Failed]];
        If[useMCA && FreeQ[rules, HoldPattern[alphaR0 -> _]], Return[$Failed]];
        rules
    ];

(* --------------------------------------------------------------- *)
(* Minimization                                                    *)
(* --------------------------------------------------------------- *)

FindMagneticEquilibrium[
    genericEnergy_, parameterRules_List, initialValues_List, OptionsPattern[]
] := Module[{hasCubicTerm, useMCA, settings, prepared, parameters,
             startingValues, specificEnergy, result, completedSolution,
             dmiType, printSteps},
    hasCubicTerm = !FreeQ[genericEnergy, kappaCub];
    useMCA = Replace[OptionValue["CubicAnisotropy"], Automatic -> hasCubicTerm];
    If[TrueQ[useMCA] =!= hasCubicTerm,
        Message[FindMagneticEquilibrium::energy, hasCubicTerm, useMCA];
        Return[$Failed]];
    If[TrueQ[useMCA] && FreeQ[parameterRules, HoldPattern[kappaCub -> _]],
        Message[FindMagneticEquilibrium::kappa]; Return[$Failed]];

    settings = GroundStateModelSettings[
        parameterRules, OptionValue["FourierCutoff"],
        OptionValue["MagneticRepresentation"], useMCA,
        OptionValue["FieldDirection"], OptionValue["RigidLattice"]];
    If[settings === $Failed, Return[$Failed]];
    prepared = GroundStateInitialValues[initialValues, settings, parameterRules];
    If[prepared === $Failed,
        Message[FindMagneticEquilibrium::ansatz]; Return[$Failed]];

    parameters = prepared["Parameters"];
    startingValues = prepared["InitialValues"];
    specificEnergy = Chop[genericEnergy /. parameters /. prepared["FixedGeometryRules"]];
    printSteps = TrueQ[OptionValue["PrintSteps"]];
    result = FindMinimum[
        specificEnergy, startingValues,
        Method -> "Newton", AccuracyGoal -> 12, PrecisionGoal -> 10,
        WorkingPrecision -> 40, MaxIterations -> Infinity,
        StepMonitor :> If[printSteps,
            Print[{a0, eta0, phi0, alphaR0} /. prepared["FixedGeometryRules"] /. parameters]]
    ];
    If[!MatchQ[result, {_, _List}], Return[result]];

    completedSolution = GroundStateCompletedSolution[
        result[[2]], settings["RigidLattice"], useMCA];
    If[completedSolution === $Failed,
        Message[SaveMagneticEquilibrium::result]; Return[$Failed]];
    result = {result[[1]], completedSolution};

    If[TrueQ[OptionValue["Save"]],
        dmiType = OptionValue["DMIType"];
        If[dmiType === Automatic,
            Message[FindMagneticEquilibrium::save],
            SaveMagneticEquilibrium[result, parameters, dmiType,
                "FourierCutoff" -> settings["FourierCutoff"],
                "MagneticRepresentation" -> settings["MagneticRepresentation"],
                "CubicAnisotropy" -> useMCA,
                "FieldDirection" -> OptionValue["FieldDirection"],
                "RigidLattice" -> settings["RigidLattice"]]
        ]
    ];
    result
];

(* --------------------------------------------------------------- *)
(* Equilibrium file names                                           *)
(* --------------------------------------------------------------- *)

(* Three DECIMAL places at most: 1/3 -> 0.333, 0.4 -> 0.4. *)
(* Formatting through integer arithmetic avoids precision markers and *)
(* machine-precision tails in the resulting filename.             *)
GroundStateNzString[value_?NumericQ] := Module[{scaled, magnitude, fraction, sign},
    scaled = Round[1000 N[value, 16]];
    sign = If[scaled < 0, "-", ""];
    magnitude = Abs[scaled];
    fraction = StringReplace[
        IntegerString[Mod[magnitude, 1000], 10, 3],
        RegularExpression["0+$"] -> ""
    ];
    sign <> ToString[Quotient[magnitude, 1000]] <>
        If[fraction === "", "", "." <> fraction]
];

(* Preserve meaningful decimals for non-integral lambda.              *)
(* A numerically integral value, e.g. 500. or 500`30, becomes "500". *)
GroundStateParameterLabel[parameterRules_List, parameter_, label_String] :=
    Module[{value, printable},
        value = parameter /. parameterRules;
        If[value === parameter, Return[""]];
        printable = Which[
            parameter === Nz && NumericQ[value],
                GroundStateNzString[value],
            parameter === lambda && NumericQ[value] &&
                TrueQ[Abs[N[value - Round[value], 16]] < 10^-10],
                ToString[InputForm[Round[value]]],
            True,
                ToString[InputForm[N[value, 7]]]
        ];
        printable = StringReplace[printable, {
            RegularExpression["`[0-9.]*"] -> "",
            RegularExpression["\\.$"] -> "",
            "*^" -> "e", " " -> "", "-" -> "m", "/" -> "_"
        }];
        "_" <> label <> printable
    ];

(* Equilibrium filenames do not reuse energy-cache naming/canonicalization. *)
EquilibriumFileName[
    dmiType_String, parameterRules_List, OptionsPattern[]
] := Module[{settings, name, missing},
    settings = GroundStateModelSettings[
        parameterRules, OptionValue["FourierCutoff"],
        OptionValue["MagneticRepresentation"], OptionValue["CubicAnisotropy"],
        OptionValue["FieldDirection"], OptionValue["RigidLattice"]];
    If[settings === $Failed, Return[$Failed]];
    missing = Select[{b, xi, Nz, lambda},
        !MemberQ[First /@ parameterRules, #] &];
    If[missing =!= {},
        Message[EquilibriumFileName::params, missing]; Return[$Failed]];
    If[settings["CubicAnisotropy"] &&
       FreeQ[parameterRules, HoldPattern[kappaCub -> _]],
        Message[FindMagneticEquilibrium::kappa]; Return[$Failed]];
    name = StringJoin[
        dmiType,
        "_N", ToString[settings["FourierCutoff"]],
        GroundStateParameterLabel[parameterRules, b, "b"],
        GroundStateParameterLabel[parameterRules, xi, "xi"],
        GroundStateParameterLabel[parameterRules, Nz, "Nz"],
        GroundStateParameterLabel[parameterRules, lambda, "lambda"],
        If[settings["CubicAnisotropy"],
            GroundStateParameterLabel[parameterRules, kappaCub, "Kcub"] <>
                "_" <> settings["FieldDirectionTag"], ""],
        "_", ToLowerCase[settings["MagneticRepresentation"]],
        If[settings["RigidLattice"], "_rigidLattice", ""]
    ];
    name
];

EquilibriumFilePath[dmiType_String, parameterRules_List, OptionsPattern[]] :=
    Module[{name},
        name = EquilibriumFileName[dmiType, parameterRules,
            "FourierCutoff" -> OptionValue["FourierCutoff"],
            "MagneticRepresentation" -> OptionValue["MagneticRepresentation"],
            "CubicAnisotropy" -> OptionValue["CubicAnisotropy"],
            "FieldDirection" -> OptionValue["FieldDirection"],
            "RigidLattice" -> OptionValue["RigidLattice"]];
        If[name === $Failed, Return[$Failed]];
        FileNameJoin[{$ProjectRoot, "results", "equilibria", name <> ".txt"}]
    ];

(* Compare full physical parameters (also checks Nx and Ny, not in filename). *)
(* Compare physical values rather than their printed representations. *)
(* This permits restarting a sweep with numerically equivalent input. *)
GroundStateValuesEqualQ[a_, c_] :=
    If[NumericQ[a] && NumericQ[c],
        TrueQ[
            Abs[N[a - c, 16]] <=
                10^-10 Max[1, Abs[N[a, 16]], Abs[N[c, 16]]]
        ],
        TrueQ[Simplify[a == c]]
    ];

GroundStateParametersMatchQ[first_List, second_List] := Module[{a, c},
    a = Association[Select[first,
        !MemberQ[{alphaR0, eta0, phi0}, First[#]] &]];
    c = Association[Select[second,
        !MemberQ[{alphaR0, eta0, phi0}, First[#]] &]];
    Sort[Keys[a]] === Sort[Keys[c]] &&
        AllTrue[Keys[a], GroundStateValuesEqualQ[a[#], c[#]] &]
];

GroundStateMatchesRequestQ[
    data_Association, dmiType_String, parameters_List, settings_Association
] :=
    Lookup[data, "DMIType", Missing[]] === dmiType &&
    Lookup[data, "FourierCutoff", Missing[]] === settings["FourierCutoff"] &&
    Lookup[data, "MagneticRepresentation", Missing[]] === settings["MagneticRepresentation"] &&
    TrueQ[Lookup[data, "CubicAnisotropy", False]] === settings["CubicAnisotropy"] &&
    TrueQ[Lookup[data, "RigidLattice", False]] === settings["RigidLattice"] &&
    (!settings["CubicAnisotropy"] ||
        Lookup[data, "FieldDirectionInput", Missing[]] === settings["FieldDirectionInput"]) &&
    GroundStateParametersMatchQ[Lookup[data, "Parameters", {}], parameters];

(* Locate a saved equilibrium under the new filename or an older name. *)
(* Legacy filenames are accepted ONLY when their complete model metadata *)
(* and physical parameters match the requested configuration.          *)
GroundStateLocateSavedEquilibrium[
    dmiType_String, parameterRules_List, settings_Association
] := Module[{path, directory, pattern, candidates, match},
    path = EquilibriumFilePath[dmiType, parameterRules,
        "FourierCutoff" -> settings["FourierCutoff"],
        "MagneticRepresentation" -> settings["MagneticRepresentation"],
        "CubicAnisotropy" -> settings["CubicAnisotropy"],
        "FieldDirection" -> settings["FieldDirectionInput"],
        "RigidLattice" -> settings["RigidLattice"]];
    If[path === $Failed, Return[$Failed]];
    If[FileExistsQ[path], Return[path]];

    directory = DirectoryName[path];
    If[!DirectoryQ[directory], Return[$Failed]];
    pattern = dmiType <> "_N" <> ToString[settings["FourierCutoff"]] <>
        "_*_" <> ToLowerCase[settings["MagneticRepresentation"]] <>
        If[settings["RigidLattice"], "_rigidLattice", ""] <> ".txt";
    candidates = FileNames[pattern, directory];
    match = SelectFirst[candidates,
        Function[file,
            Module[{data = Quiet[Check[Get[file], $Failed]]},
                AssociationQ[data] &&
                    GroundStateMatchesRequestQ[data, dmiType, parameterRules, settings]
            ]
        ],
        $Failed
    ];
    match
];

EquilibriumExistsQ[dmiType_String, parameterRules_List, OptionsPattern[]] :=
    Module[{path, data, settings},
        settings = GroundStateModelSettings[parameterRules,
            OptionValue["FourierCutoff"], OptionValue["MagneticRepresentation"],
            OptionValue["CubicAnisotropy"], OptionValue["FieldDirection"],
            OptionValue["RigidLattice"]];
        If[settings === $Failed, Return[False]];
        path = GroundStateLocateSavedEquilibrium[dmiType, parameterRules, settings];
        If[path === $Failed, Return[False]];
        data = Quiet[Check[Get[path], $Failed]];
        AssociationQ[data] &&
            GroundStateMatchesRequestQ[data, dmiType, parameterRules, settings]
    ];

(* --------------------------------------------------------------- *)
(* Save / load                                                      *)
(* --------------------------------------------------------------- *)

SaveMagneticEquilibrium[
    result_, parameterRules_List, dmiType_String, OptionsPattern[]
] := Module[{settings, solution, parameters, path, data, existing},
    If[!MatchQ[result, {_, _List}],
        Message[SaveMagneticEquilibrium::result]; Return[$Failed]];
    settings = GroundStateModelSettings[parameterRules,
        OptionValue["FourierCutoff"], OptionValue["MagneticRepresentation"],
        OptionValue["CubicAnisotropy"], OptionValue["FieldDirection"],
        OptionValue["RigidLattice"]];
    If[settings === $Failed, Return[$Failed]];
    parameters = Select[parameterRules,
        !MemberQ[{eta0, phi0, alphaR0}, First[#]] &];
    If[settings["CubicAnisotropy"],
        If[FreeQ[parameters, HoldPattern[kappaCub -> _]],
            Message[FindMagneticEquilibrium::kappa]; Return[$Failed]],
        AppendTo[parameters, alphaR0 -> 0]
    ];
    solution = GroundStateCompletedSolution[
        result[[2]], settings["RigidLattice"], settings["CubicAnisotropy"]];
    If[solution === $Failed,
        Message[SaveMagneticEquilibrium::result]; Return[$Failed]];

    path = EquilibriumFilePath[dmiType, parameters,
        "FourierCutoff" -> settings["FourierCutoff"],
        "MagneticRepresentation" -> settings["MagneticRepresentation"],
        "CubicAnisotropy" -> settings["CubicAnisotropy"],
        "FieldDirection" -> OptionValue["FieldDirection"],
        "RigidLattice" -> settings["RigidLattice"]];
    If[path === $Failed, Return[$Failed]];
    If[FileExistsQ[path],
        existing = Quiet[Check[Get[path], $Failed]];
        If[!AssociationQ[existing] ||
           !GroundStateMatchesRequestQ[existing, dmiType, parameters, settings],
            Message[SaveMagneticEquilibrium::collision, path]; Return[$Failed]]
    ];
    If[!DirectoryQ[DirectoryName[path]],
        CreateDirectory[DirectoryName[path], CreateIntermediateDirectories -> True]];

    data = <|
        "Energy" -> result[[1]],
        "Solution" -> solution,
        "Parameters" -> parameters,
        "DMIType" -> dmiType,
        "FourierCutoff" -> settings["FourierCutoff"],
        "MagneticRepresentation" -> settings["MagneticRepresentation"],
        "CubicAnisotropy" -> settings["CubicAnisotropy"],
        "RigidLattice" -> settings["RigidLattice"],
        "FieldDirectionInput" -> settings["FieldDirectionInput"],
        "FieldDirectionCrystal" -> settings["FieldDirectionCrystal"],
        "FieldDirectionTag" -> settings["FieldDirectionTag"]
    |>;
    Put[data, path];
    Print["Equilibrium configuration saved to: ", path];
    path
];

(* Direct filenames remain usable for older equilibria. *)
LoadMagneticEquilibrium[fileName_String] := Module[{path},
    path = If[FileExistsQ[fileName], fileName,
        FileNameJoin[{$ProjectRoot, "results", "equilibria", fileName}]];
    If[!FileExistsQ[path],
        Print["Equilibrium file not found: ", path]; Return[$Failed]];
    Get[path]
];

LoadMagneticEquilibrium[
    dmiType_String, parameterRules_List, OptionsPattern[]
] := Module[{settings, path, data},
    settings = GroundStateModelSettings[parameterRules,
        OptionValue["FourierCutoff"], OptionValue["MagneticRepresentation"],
        OptionValue["CubicAnisotropy"], OptionValue["FieldDirection"],
        OptionValue["RigidLattice"]];
    If[settings === $Failed, Return[$Failed]];
    path = GroundStateLocateSavedEquilibrium[dmiType, parameterRules, settings];
    If[path === $Failed,
        Print["No matching saved equilibrium found for ", dmiType,
            " and parameters ", parameterRules];
        Return[$Failed]
    ];
    data = LoadMagneticEquilibrium[path];
    If[!AssociationQ[data] ||
       !GroundStateMatchesRequestQ[data, dmiType, parameterRules, settings],
        Message[LoadMagneticEquilibrium::metadata, path]; Return[$Failed]];
    data
];

EquilibriumRules[equilibrium_Association] :=
    Join[equilibrium["Solution"], equilibrium["Parameters"]];

(* --------------------------------------------------------------- *)
(* Restart / continuation                                           *)
(* --------------------------------------------------------------- *)

EquilibriumAnsatz[equilibrium_Association, OptionsPattern[]] :=
    Module[{oldCutoff, oldRep, cutoff, rep, useMCA, rigid, settings, prepared},
        oldCutoff = equilibrium["FourierCutoff"];
        oldRep = equilibrium["MagneticRepresentation"];
        (* A saved equilibrium is self-describing: do not use unrelated globals. *)
        cutoff = Replace[OptionValue["FourierCutoff"], Automatic -> oldCutoff];
        rep = Replace[OptionValue["MagneticRepresentation"], Automatic -> oldRep];
        If[oldCutoff =!= cutoff,
            Message[EquilibriumAnsatz::cutoff, oldCutoff, cutoff]; Return[$Failed]];
        If[oldRep =!= rep,
            If[!TrueQ[OptionValue["AllowRepresentationChange"]],
                Message[EquilibriumAnsatz::rep, oldRep, rep]; Return[$Failed]];
            If[!(oldRep === "SixFold" && rep === "General"),
                Message[EquilibriumAnsatz::conversion, oldRep, rep]; Return[$Failed]]
        ];
        useMCA = Replace[OptionValue["CubicAnisotropy"],
            Automatic :> Lookup[equilibrium, "CubicAnisotropy", False]];
        rigid = Replace[OptionValue["RigidLattice"],
            Automatic :> Lookup[equilibrium, "RigidLattice", False]];
        settings = GroundStateModelSettings[
            equilibrium["Parameters"], cutoff, rep, useMCA,
            OptionValue["FieldDirection"], rigid];
        If[settings === $Failed, Return[$Failed]];
        prepared = GroundStateInitialValues[
            List @@@ equilibrium["Solution"], settings, equilibrium["Parameters"]];
        If[prepared === $Failed, Return[$Failed]];
        prepared["InitialValues"]
    ];

ParameterSweep[
    genericEnergy_, dmiType_String, baseParameterRules_List,
    sweepParameter_Symbol, parameterValues_List, initialValues_List,
    save_: True, OptionsPattern[]
] := Module[{useMCA, settings, ansatz, currentParameters, prepared,
             parameters, path, equilibrium, result, source, rows = {}},
    useMCA = Replace[OptionValue["CubicAnisotropy"],
        Automatic -> !FreeQ[genericEnergy, kappaCub]];
    settings = GroundStateModelSettings[baseParameterRules,
        OptionValue["FourierCutoff"], OptionValue["MagneticRepresentation"],
        useMCA, OptionValue["FieldDirection"], OptionValue["RigidLattice"]];
    If[settings === $Failed, Return[$Failed]];
    ansatz = initialValues;
    Do[
        currentParameters = Join[
            DeleteCases[baseParameterRules, HoldPattern[Rule[sweepParameter, _]]],
            {sweepParameter -> value}];
        prepared = GroundStateInitialValues[ansatz, settings, currentParameters];
        If[prepared === $Failed,
            Message[ParameterSweep::ansatz, sweepParameter, value]; Return[$Failed]];
        parameters = prepared["Parameters"];
        path = GroundStateLocateSavedEquilibrium[dmiType, parameters, settings];
        If[path =!= $Failed,
            equilibrium = Quiet[Check[Get[path], $Failed]];
            If[!AssociationQ[equilibrium] ||
                !GroundStateMatchesRequestQ[equilibrium, dmiType, parameters, settings],
                Message[LoadMagneticEquilibrium::metadata, path];
                Return[$Failed]
            ];
            Print["Loading saved equilibrium: ", FileNameTake[path]];
            result = {equilibrium["Energy"], equilibrium["Solution"]};
            source = "Loaded";
            ,
            result = FindMagneticEquilibrium[
                genericEnergy, parameters, prepared["InitialValues"],
                "FourierCutoff" -> settings["FourierCutoff"],
                "MagneticRepresentation" -> settings["MagneticRepresentation"],
                "CubicAnisotropy" -> useMCA,
                "FieldDirection" -> OptionValue["FieldDirection"],
                "RigidLattice" -> settings["RigidLattice"]];
            If[!MatchQ[result, {_, _List}], Return[$Failed]];
            source = "Calculated";
            If[TrueQ[save],
                If[SaveMagneticEquilibrium[result, parameters, dmiType,
                    "FourierCutoff" -> settings["FourierCutoff"],
                    "MagneticRepresentation" -> settings["MagneticRepresentation"],
                    "CubicAnisotropy" -> useMCA,
                    "FieldDirection" -> OptionValue["FieldDirection"],
                    "RigidLattice" -> settings["RigidLattice"]] === $Failed,
                    Return[$Failed]]]
        ];
        ansatz = List @@@ result[[2]];
        AppendTo[rows, <|
            "SweepParameter" -> sweepParameter, "ParameterValue" -> value,
            "Energy" -> result[[1]], "Solution" -> result[[2]],
            "Parameters" -> parameters, "Source" -> source|>];
        , {value, parameterValues}
    ];
    rows
];

(* --------------------------------------------------------------- *)
(* Help                                                            *)
(* --------------------------------------------------------------- *)

GroundStateHelp[] := Column[{
    Style["GroundState.wl", Bold, 16],
    "Equilibrium minimization, storage and sweeps with optional MCA and rigid lattice.",
    Grid[{
        {"Function", "Description"},
        {"FindMagneticEquilibrium", "Minimize the energy; optionally fix eta0=1, phi0=0 via RigidLattice."},
        {"GroundStateInitialValues", "Generate a complete Fourier ansatz, with alphaR0 fixed for isotropic models."},
        {"EquilibriumFileName / EquilibriumFilePath", "Descriptive names: Nz to 3 decimals, integer-valued lambda without decimal, raw H indices, representation and rigidLattice."},
        {"SaveMagneticEquilibrium / LoadMagneticEquilibrium", "Preserve geometry in Solution and validate model metadata."},
        {"EquilibriumRules / EquilibriumAnsatz", "Read equilibria or create restart ansatz, including SixFold-to-General expansion."},
        {"ParameterSweep", "Continue an equilibrium branch; recognize older save filenames by validating metadata."}
    }, Frame -> All, Alignment -> Left],
    "Options: FourierCutoff, MagneticRepresentation, CubicAnisotropy, FieldDirection, RigidLattice."
}];

Print[GroundStateHelp[]];