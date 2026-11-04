(* ================================================================== *)
(* MLSeedSweep.wl                                                    *)
(*                                                                  *)
(* Batch equilibrium calculations from Python-generated Sobol CSVs. *)
(* Saves only to results/equilibria/ml_seed/, without changing the   *)
(* normal GroundState.wl save/load paths.                            *)
(*                                                                  *)
(* Requires: Geometry.wl, FourierBasis.wl, GroundState.wl,           *)
(*           plus the energy modules already used by the project.   *)
(* ================================================================== *)

ClearAll[
    ReadMLSeedCSV,
    MLSeedEquilibriumPath,
    MLSeedParameterSweep,
    MLSeedMatchesRequestQ,
    MLSeedNumericField,
    MLSeedSweepHelp
];

ReadMLSeedCSV::file = "Cannot read the Sobol CSV file `1`.";
ReadMLSeedCSV::header = "CSV `1` requires sample_id, b and xi columns.";
ReadMLSeedCSV::row = "Invalid or inconsistent Sobol CSV row `1` in `2`.";
ReadMLSeedCSV::duplicate = "Duplicate sample IDs found in CSV `1`.";
MLSeedParameterSweep::model = "Requested CubicAnisotropy setting does not match the symbolic energy.";
MLSeedParameterSweep::root = "Set $ProjectRoot to your project directory before running the ML sweep.";
MLSeedParameterSweep::collision = "Existing ML equilibrium `1` does not match the requested sample or model. No overwrite performed.";
MLSeedParameterSweep::solve = "Equilibrium minimization failed for sample `1`.";
MLSeedParameterSweep::save = "Could not save equilibrium for sample `1` to `2`.";
MLSeedParameterSweep::limit = "MaxSamples must be All or a positive integer (received `1`).";

(* Only numeric CSV fields are parsed as Wolfram expressions. *)
MLSeedNumericField[value_] := Module[{s, converted},
    s = StringTrim[If[StringQ[value], value, ToString[value]]];
    If[!StringMatchQ[s, RegularExpression[
         "^[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?$"
        ]], Return[$Failed]];
    converted = Quiet[Check[
        ToExpression[StringReplace[s, {"e" -> "*^", "E" -> "*^"}]],
        $Failed
    ]];
    If[!NumericQ[converted] || !TrueQ[Abs[N[converted]] < Infinity],
        $Failed, N[converted]]
];

(* Accept split-specific CSVs and the combined all_samples.csv. *)
ReadMLSeedCSV[csvPath_String] := Module[
    {rows, header, required, parsed, i},
    If[!FileExistsQ[csvPath],
        Message[ReadMLSeedCSV::file, csvPath]; Return[$Failed]];
    rows = Quiet[Check[Import[csvPath, "CSV"], $Failed]];
    If[rows === $Failed || !ListQ[rows] || Length[rows] < 2,
        Message[ReadMLSeedCSV::file, csvPath]; Return[$Failed]];
    header = StringTrim /@ (If[StringQ[#], #, ToString[#]] & /@ First[rows]);
    required = {"sample_id", "b", "xi"};
    If[!TrueQ[And @@ (MemberQ[header, #] & /@ required)],
        Message[ReadMLSeedCSV::header, csvPath]; Return[$Failed]];

    parsed = Table[
        Module[{row = rows[[i]], rowFields, rowID, rowSplit, rowB, rowXi},
            If[Length[row] =!= Length[header],
                Message[ReadMLSeedCSV::row, i, csvPath]; Return[$Failed]];
            rowFields = AssociationThread[header, row];
            rowID = StringTrim[If[StringQ[rowFields["sample_id"]],
                rowFields["sample_id"], ToString[rowFields["sample_id"]]]];
            rowSplit = First[StringSplit[rowID, "_"]];
            rowB = MLSeedNumericField[rowFields["b"]];
            rowXi = MLSeedNumericField[rowFields["xi"]];
            If[!StringMatchQ[rowID, RegularExpression[
                  "^(train|validation|test)_[0-9]+$"
                ]] || rowB === $Failed || rowXi === $Failed ||
                (KeyExistsQ[rowFields, "split"] &&
                    rowFields["split"] =!= rowSplit),
                Message[ReadMLSeedCSV::row, i, csvPath]; Return[$Failed]];
            <|"SampleID" -> rowID, "Split" -> rowSplit,
              "b" -> rowB, "xi" -> rowXi|>
        ],
        {i, 2, Length[rows]}
    ];
    If[MemberQ[parsed, $Failed], Return[$Failed]];
    If[Length[DeleteDuplicates[Lookup[parsed, "SampleID"]]] =!= Length[parsed],
        Message[ReadMLSeedCSV::duplicate, csvPath]; Return[$Failed]];
    parsed
];

(* Name includes split + sample ID to guarantee unique files despite
   compact formatting of b and xi in EquilibriumFileName. *)
MLSeedEquilibriumPath[
    sample_Association, dmiType_String,
    parameterRules_List, settings_Association
] := Module[{baseName},
    baseName = EquilibriumFileName[dmiType, parameterRules,
        "FourierCutoff" -> settings["FourierCutoff"],
        "MagneticRepresentation" -> settings["MagneticRepresentation"],
        "CubicAnisotropy" -> settings["CubicAnisotropy"],
        "FieldDirection" -> If[settings["CubicAnisotropy"],
            settings["FieldDirectionInput"], {0, 0, 1}],
        "RigidLattice" -> settings["RigidLattice"]];
    If[baseName === $Failed, Return[$Failed]];
    FileNameJoin[{$ProjectRoot, "results", "equilibria", "ml_seed",
        sample["SampleID"] <> "_" <> baseName <> ".txt"}]
];

(* Data remain compatible with EquilibriumRules, plotting and dynamics.
   MLSeedSampleID and MLSeedSplit distinguish the ML datasets. *)
MLSeedMatchesRequestQ[
    data_Association, sample_Association, dmiType_String,
    parameters_List, settings_Association
] :=
    Lookup[data, "MLSeedSampleID", Missing["Absent"]] ===
        sample["SampleID"] &&
    Lookup[data, "MLSeedSplit", Missing["Absent"]] ===
        sample["Split"] &&
    GroundStateMatchesRequestQ[data, dmiType, parameters, settings] &&
    ListQ[Lookup[data, "Solution", Missing["Absent"]]] &&
    !FreeQ[data["Solution"], HoldPattern[eta0 -> _]] &&
    !FreeQ[data["Solution"], HoldPattern[phi0 -> _]];

Options[MLSeedParameterSweep] = Join[
    $GroundStateModelOptions,
    {"UseNearestSeed" -> True, "PrintSteps" -> False,
     "MaxSamples" -> All}
];

(* Run one CSV's Sobol configurations and save/reuse each result.

   Signature:
     MLSeedParameterSweep[energy, "Bloch", fixedRules,
         "datasets/magnon_spectra/sobol/train_64.csv",
         initialValues, "FourierCutoff" -> 4, ...]

   fixedRules normally include lambda, Nx, Ny, Nz (not b, xi).
   Missing Fourier components in an ansatz are initialized by
   GroundStateInitialValues. The nearest completed sample can serve
   as a restart seed; set UseNearestSeed -> False to disable.
*)
MLSeedParameterSweep[
    genericEnergy_, dmiType_String, fixedParameterRules_List,
    csvPath_String, initialValues_List, OptionsPattern[]
] := Module[
    {samples, settings, hasMCA, useMCA, completed = {}, rows = {},
     sample, parameters, prepared, path, saved, result, source,
     seed, nearest, i, data, directory, tempPath, bScale, xiScale,
     maxSamples, writeOK},

    If[!ValueQ[$ProjectRoot] || !DirectoryQ[$ProjectRoot],
        Message[MLSeedParameterSweep::root]; Return[$Failed]];
    samples = ReadMLSeedCSV[csvPath];
    If[samples === $Failed, Return[$Failed]];
    maxSamples = OptionValue["MaxSamples"];
    If[maxSamples =!= All &&
       (!IntegerQ[maxSamples] || maxSamples < 1),
        Message[MLSeedParameterSweep::limit, maxSamples]; Return[$Failed]];
    If[IntegerQ[maxSamples], samples = Take[samples, UpTo[maxSamples]]];
    bScale = Max[10^-12, Max[Lookup[samples, "b"]] - Min[Lookup[samples, "b"]]];
    xiScale = Max[10^-12, Max[Lookup[samples, "xi"]] - Min[Lookup[samples, "xi"]]];

    hasMCA = !FreeQ[genericEnergy, kappaCub];
    useMCA = Replace[OptionValue["CubicAnisotropy"], Automatic -> hasMCA];
    If[TrueQ[useMCA] =!= hasMCA,
        Message[MLSeedParameterSweep::model]; Return[$Failed]];
    settings = GroundStateModelSettings[
        fixedParameterRules,
        OptionValue["FourierCutoff"],
        OptionValue["MagneticRepresentation"], useMCA,
        OptionValue["FieldDirection"], OptionValue["RigidLattice"]];
    If[settings === $Failed, Return[$Failed]];

    directory = FileNameJoin[{$ProjectRoot, "results", "equilibria", "ml_seed"}];
    If[!DirectoryQ[directory],
        CreateDirectory[directory, CreateIntermediateDirectories -> True]];

    Do[
        sample = samples[[i]];
        parameters = Join[
            Select[fixedParameterRules,
                !MemberQ[{b, xi, eta0, phi0, alphaR0}, First[#]] &],
            {b -> sample["b"], xi -> sample["xi"]}
        ];

        (* CSV order is not spatial order. Nearest finished configuration
           is typically a better warm start than the previous CSV row. *)
        seed = initialValues;
        If[TrueQ[OptionValue["UseNearestSeed"]] && completed =!= {},
            nearest = First[Ordering[
                (((#["b"] - sample["b"])/bScale)^2 +
                 ((#["xi"] - sample["xi"])/xiScale)^2) & /@ completed,
                1
            ]];
            seed = List @@@ completed[[nearest]]["Solution"]
        ];
        prepared = GroundStateInitialValues[seed, settings, parameters];
        If[prepared === $Failed, Return[$Failed]];
        parameters = prepared["Parameters"];
        path = MLSeedEquilibriumPath[sample, dmiType, parameters, settings];
        If[path === $Failed, Return[$Failed]];

        Print["[", i, "/", Length[samples], "] ", sample["SampleID"],
            "  (b = ", sample["b"], ", xi = ", sample["xi"], ")"];

        If[FileExistsQ[path],
            saved = Quiet[Check[Get[path], $Failed]];
            If[!AssociationQ[saved] ||
               !MLSeedMatchesRequestQ[saved, sample, dmiType,
                   parameters, settings],
                Message[MLSeedParameterSweep::collision, path];
                Return[$Failed]];
            Print["  Loaded existing configuration."];
            result = {saved["Energy"], saved["Solution"]};
            source = "Loaded",

            result = FindMagneticEquilibrium[
                genericEnergy, parameters, prepared["InitialValues"],
                "Save" -> False,
                "FourierCutoff" -> settings["FourierCutoff"],
                "MagneticRepresentation" -> settings["MagneticRepresentation"],
                "CubicAnisotropy" -> settings["CubicAnisotropy"],
                "FieldDirection" -> OptionValue["FieldDirection"],
                "RigidLattice" -> settings["RigidLattice"],
                "PrintSteps" -> OptionValue["PrintSteps"]
            ];
            If[!MatchQ[result, {_, _List}],
                Message[MLSeedParameterSweep::solve, sample["SampleID"]];
                Return[$Failed]];

            data = <|
                "Energy" -> result[[1]],
                "Solution" -> result[[2]],
                "Parameters" -> parameters,
                "DMIType" -> dmiType,
                "FourierCutoff" -> settings["FourierCutoff"],
                "MagneticRepresentation" -> settings["MagneticRepresentation"],
                "CubicAnisotropy" -> settings["CubicAnisotropy"],
                "RigidLattice" -> settings["RigidLattice"],
                "FieldDirectionInput" -> settings["FieldDirectionInput"],
                "FieldDirectionCrystal" -> settings["FieldDirectionCrystal"],
                "FieldDirectionTag" -> settings["FieldDirectionTag"],
                "MLSeedSampleID" -> sample["SampleID"],
                "MLSeedSplit" -> sample["Split"],
                "MLSeedSourceCSV" -> FileNameTake[csvPath]
            |>;
            (* Avoid leaving a half-written result if computation is interrupted. *)
            tempPath = path <> ".tmp-" <> CreateUUID[];
            writeOK = Quiet[Check[
                Put[data, tempPath];
                FileExistsQ[tempPath] && RenameFile[tempPath, path] =!= $Failed,
                False
            ]];
            If[!TrueQ[writeOK] || !FileExistsQ[path],
                If[FileExistsQ[tempPath], DeleteFile[tempPath]];
                Message[MLSeedParameterSweep::save, sample["SampleID"], path];
                Return[$Failed]
            ];
            Print["  Saved: ", FileNameTake[path]];
            source = "Calculated"
        ];

        AppendTo[completed, <|"b" -> sample["b"], "xi" -> sample["xi"],
            "Solution" -> result[[2]]|>];
        AppendTo[rows, <|
            "SampleID" -> sample["SampleID"], "Split" -> sample["Split"],
            "b" -> sample["b"], "xi" -> sample["xi"],
            "Energy" -> result[[1]], "Solution" -> result[[2]],
            "Path" -> path, "Source" -> source
        |>],
        {i, Length[samples]}
    ];
    rows
];

MLSeedSweepHelp[] := Column[{
    Style["MLSeedSweep.wl", Bold, 16],
    "Run or resume train / validation / test Sobol equilibrium calculations.",
    "Input CSVs: train_64.csv, validation_16.csv, test_32.csv, or all_samples.csv.",
    "Output: results/equilibria/ml_seed/<sample_id>_<configuration>.txt",
    "The saved objects match GroundState equilibrium associations, with MLSeedSampleID and MLSeedSplit metadata.",
    "Options: FourierCutoff, MagneticRepresentation, CubicAnisotropy, FieldDirection, RigidLattice, UseNearestSeed, MaxSamples, PrintSteps."
}];

Print["Loaded MLSeedSweep.wl. Call MLSeedSweepHelp[] for usage."];