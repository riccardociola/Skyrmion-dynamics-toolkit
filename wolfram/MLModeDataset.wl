(* ================================================================ *)
(* MLModeDataset.wl                                                *)
(* Export Gamma-point spectra for Sobol equilibrium samples.        *)
(*                                                                  *)
(* Load first: Geometry, FourierBasis, GroundState,                *)
(* DynamicsProjection, DynamicsMagnetic, ModeAnalysis, MLSeedSweep. *)
(* Isotropic Bloch skyrmion lattices only.                           *)
(* ================================================================ *)

ClearAll[
    MLModeSampleFile, MLModeLabels, MLModeJSONRecord,
    MLModeMatchesRequestQ, MLModeDatasetSweep, MLModeDatasetHelp
];

MLModeDatasetSweep::csv = "Cannot parse Sobol CSV `1`. Load MLSeedSweep.wl first.";
MLModeDatasetSweep::root = "Set $ProjectRoot to the project root directory.";
MLModeDatasetSweep::ambiguous = "Expected one equilibrium for `1`, found `2`: `3`.";
MLModeDatasetSweep::badmodel = "Equilibrium `1` is not an isotropic Bloch configuration.";
MLModeDatasetSweep::conflict = "Existing mode export `1` does not match this request; delete it or use Overwrite -> True.";
MLModeDatasetSweep::calc = "Failed to compute or analyze modes for `1`.";
MLModeDatasetSweep::write = "Failed to export `1`.";

Options[MLModeDatasetSweep] = {
    "FrequencyMax" -> 3.,
    "ProjectionCutoff" -> 0.95,
    "SymmetryTolerance" -> 10^-6,
    "MainBranchThreshold" -> 0.001,
    "Overwrite" -> False,
    "SkipMissing" -> True,
    "MaxSamples" -> All
};

(* A saved Sobol equilibrium carries a unique sample ID. *)
MLModeSampleFile[sample_Association] := Module[{paths, folder},
    folder = FileNameJoin[{$ProjectRoot, "results", "equilibria", "ml_seed"}];
    paths = FileNames[sample["SampleID"] <> "_*.txt", folder];
    paths
];

(* The three principal dipole-active modes are the strongest members
   of m = 0, +1, -1 in their respective coupling channels. We only
   attach these labels if their weight is significant relative to the
   maximum dipolar coupling. All other classified modes are Multipolar.
   No mode is discarded or filtered based on its dipole weight. *)
MLModeLabels[modeRows_List, threshold_?NumericQ] := Module[
    {channels, peakWeights, globalMax, labels, indices, best, weight},
    channels = {{0, "Breathing", "breathing_weight"},
                {1, "CW", "cw_weight"}, {-1, "CCW", "ccw_weight"}};
    labels = Table[If[IntegerQ[row["m_signed"]], "Multipolar",
                      "Unclassified"], {row, modeRows}];
    peakWeights = Flatten[Table[
        Lookup[modeRows, key, {}], {key,
            {"breathing_weight", "cw_weight", "ccw_weight"}}]];
    globalMax = If[peakWeights === {}, 0., Max[peakWeights]];
    Do[
        indices = Select[Range[Length[modeRows]],
            modeRows[[#, "m_signed"]] === channel[[1]] &];
        If[indices =!= {},
            best = First[MaximalBy[indices,
                modeRows[[#, channel[[3]]]] &]];
            weight = modeRows[[best, channel[[3]]]];
            If[globalMax > 0 && weight >= threshold globalMax,
                labels[[best]] = channel[[2]]]
        ],
        {channel, channels}
    ];
    MapThread[Append[#1, "label" -> #2] &, {modeRows, labels}]
];

MLModeJSONRecord[
    sample_Association, equilibrium_Association,
    modes_Association, angular_Association, dipole_Association,
    omegaMax_?NumericQ, projectionCutoff_?NumericQ,
    symmetryTolerance_?NumericQ, threshold_?NumericQ, equilibriumPath_String
] := Module[{frequencies, rows, a, d, m, residual, model},
    frequencies = modes["Frequencies"];
    rows = Table[
        a = angular["Modes"][[j, "AngularMomentum"]];
        d = dipole["Modes"][[j, "Dipole"]];
        m = a["AngularMomentumSigned"];
        residual = a["BestResidual"];
        <|
            "mode_index" -> j,
            "frequency" -> N[frequencies[[j]]],
            "m_signed" -> If[IntegerQ[m], m, Null],
            "m_mod6" -> If[IntegerQ[m], Mod[m, 6], Null],
            "c6_residual" -> If[NumericQ[residual], N[residual], Null],
            "dipole_weight" -> N[d["TotalDipoleWeight"]],
            "breathing_weight" -> N[d["BreathingWeight"]],
            "cw_weight" -> N[d["CWWeight"]],
            "ccw_weight" -> N[d["CCWWeight"]]
        |>,
        {j, Length[frequencies]}
    ];
    rows = MLModeLabels[rows, threshold];
    model = <|
        "dmi_type" -> equilibrium["DMIType"],
        "fourier_cutoff" -> equilibrium["FourierCutoff"],
        "magnetic_representation" -> equilibrium["MagneticRepresentation"],
        "rigid_lattice" -> TrueQ[Lookup[equilibrium, "RigidLattice", False]],
        "lambda" -> N[lambda /. equilibrium["Parameters"]],
        "Nx" -> N[Nx /. equilibrium["Parameters"]],
        "Ny" -> N[Ny /. equilibrium["Parameters"]],
        "Nz" -> N[Nz /. equilibrium["Parameters"]],
        "cubic_anisotropy" -> False
    |>;
    <|
        "schema_version" -> 1,
        "sample_id" -> sample["SampleID"],
        "split" -> sample["Split"],
        "b" -> N[sample["b"]],
        "xi" -> N[sample["xi"]],
        "omega_max" -> N[omegaMax],
        "symmetry_tolerance" -> N[symmetryTolerance],
    "projection_cutoff" -> N[projectionCutoff],
        "main_branch_threshold" -> N[threshold],
        "equilibrium_file" -> FileNameTake[equilibriumPath],
        "model" -> model,
        "number_of_modes" -> Length[rows],
        "modes" -> rows
    |>
];

MLModeMatchesRequestQ[data_Association, sample_Association,
    omegaMax_, projectionCutoff_, symmetryTolerance_, threshold_, equilibriumPath_] :=
    Lookup[data, "schema_version", None] === 1 &&
    Lookup[data, "sample_id", None] === sample["SampleID"] &&
    Lookup[data, "split", None] === sample["Split"] &&
    StringQ[Lookup[data, "equilibrium_file", None]] &&
    data["equilibrium_file"] === FileNameTake[equilibriumPath] &&
    AllTrue[{{"b", sample["b"]}, {"xi", sample["xi"]},
             {"omega_max", omegaMax}, {"projection_cutoff", projectionCutoff},
             {"symmetry_tolerance", symmetryTolerance},
             {"main_branch_threshold", threshold}},
        KeyExistsQ[data, #[[1]]] &&
        NumericQ[data[#[[1]]]] &&
        Abs[data[#[[1]]] - N[#[[2]]]] < 10^-10 &];

(* Export a CSV split, allowing incomplete equilibrium sweeps. *)
MLModeDatasetSweep[csvPath_String, OptionsPattern[]] := Module[
    {samples, maxSamples, omegaMax, projectionCutoff, symmetryTolerance,
     threshold, rows = {}, sample, eqFiles, eqPath, eq, outputFolder,
     outputPath, existing, dynamics, modes, angular, dipole, record,
     tmp, success, source},

    If[!ValueQ[$ProjectRoot] || !StringQ[$ProjectRoot],
        Message[MLModeDatasetSweep::root]; Return[$Failed]];
    samples = ReadMLSeedCSV[csvPath];
    If[samples === $Failed,
        Message[MLModeDatasetSweep::csv, csvPath]; Return[$Failed]];
    maxSamples = OptionValue["MaxSamples"];
    If[maxSamples =!= All && !(IntegerQ[maxSamples] && maxSamples > 0),
        Return[$Failed]];
    If[IntegerQ[maxSamples], samples = Take[samples, UpTo[maxSamples]]];
    omegaMax = OptionValue["FrequencyMax"];
    projectionCutoff = OptionValue["ProjectionCutoff"];
    symmetryTolerance = OptionValue["SymmetryTolerance"];
    threshold = OptionValue["MainBranchThreshold"];

    Do[
        sample = samples[[i]];
        eqFiles = MLModeSampleFile[sample];
        If[eqFiles === {},
            Print["[", i, "/", Length[samples], "] Missing equilibrium: ",
                sample["SampleID"]];
            AppendTo[rows, <|"SampleID" -> sample["SampleID"],
                "Status" -> "Missing equilibrium"|>];
            If[!TrueQ[OptionValue["SkipMissing"]], Return[$Failed]];
            Continue[]
        ];
        If[Length[eqFiles] != 1,
            Message[MLModeDatasetSweep::ambiguous,
                sample["SampleID"], Length[eqFiles], eqFiles];
            Return[$Failed]
        ];
        eqPath = First[eqFiles];
        eq = Quiet[Check[Get[eqPath], $Failed]];
        If[!AssociationQ[eq] ||
            Lookup[eq, "MLSeedSampleID", None] =!= sample["SampleID"] ||
            Lookup[eq, "MLSeedSplit", None] =!= sample["Split"] ||
            Abs[N[b /. eq["Parameters"]] - sample["b"]] > 10^-10 ||
            Abs[N[xi /. eq["Parameters"]] - sample["xi"]] > 10^-10 ||
            TrueQ[Lookup[eq, "CubicAnisotropy", False]] ||
            Lookup[eq, "DMIType", ""] =!= "Bloch",
            Message[MLModeDatasetSweep::badmodel, eqPath]; Return[$Failed]
        ];
        outputFolder = FileNameJoin[{$ProjectRoot, "datasets",
            "magnon_spectra", "gamma", sample["Split"]}];
        If[!DirectoryQ[outputFolder],
            CreateDirectory[outputFolder, CreateIntermediateDirectories -> True]];
        outputPath = FileNameJoin[{outputFolder, sample["SampleID"] <> ".json"}];

        If[FileExistsQ[outputPath] && !TrueQ[OptionValue["Overwrite"]] &&
            AbsoluteTime[FileDate[outputPath, "Modification"]] >=
            AbsoluteTime[FileDate[eqPath, "Modification"]],
            existing = Quiet[Check[Import[outputPath, "RawJSON"], $Failed]];
            If[!AssociationQ[existing] ||
                !MLModeMatchesRequestQ[existing, sample, omegaMax,
                    projectionCutoff, symmetryTolerance, threshold, eqPath],
                Message[MLModeDatasetSweep::conflict, outputPath];
                Return[$Failed]
            ];
            Print["[", i, "/", Length[samples], "] Reusing ", sample["SampleID"]];
            AppendTo[rows, <|"SampleID" -> sample["SampleID"],
                "Status" -> "Loaded", "Path" -> outputPath|>];
            Continue[]
        ];

        Print["[", i, "/", Length[samples], "] Computing Gamma modes: ",
            sample["SampleID"]];
        dynamics = PrepareMagneticDynamics[eq,
            "ProjectionCutoff" -> projectionCutoff];
        If[dynamics === $Failed,
            Message[MLModeDatasetSweep::calc, sample["SampleID"]];
            Return[$Failed]
        ];
        modes = SolveMagneticModes[dynamics, {0, 0, 0},
            "FrequencyRange" -> {0, omegaMax}];
        If[modes === $Failed,
            Message[MLModeDatasetSweep::calc, sample["SampleID"]];
            Return[$Failed]
        ];
        If[modes["NumberOfModes"] == 0,
            angular = <|"Modes" -> {}|>;
            dipole = <|"Modes" -> {}|>,
            angular = AnalyzeAngularMomentumModes[dynamics, modes,
                "SymmetryTolerance" -> symmetryTolerance];
            dipole = AnalyzeDipoleModes[dynamics, modes]
        ];
        If[angular === $Failed || dipole === $Failed,
            Message[MLModeDatasetSweep::calc, sample["SampleID"]];
            Return[$Failed]
        ];
        record = MLModeJSONRecord[sample, eq, modes, angular, dipole,
            omegaMax, projectionCutoff, symmetryTolerance, threshold, eqPath];
        tmp = outputPath <> ".tmp-" <> CreateUUID[];
        success = Quiet[Check[
            Export[tmp, record, "RawJSON"];
            If[FileExistsQ[outputPath], DeleteFile[outputPath]];
            FileExistsQ[tmp] && RenameFile[tmp, outputPath] =!= $Failed,
            False
        ]];
        If[!TrueQ[success] || !FileExistsQ[outputPath],
            If[FileExistsQ[tmp], DeleteFile[tmp]];
            Message[MLModeDatasetSweep::write, outputPath];
            Return[$Failed]
        ];
        AppendTo[rows, <|"SampleID" -> sample["SampleID"],
            "Status" -> "Calculated", "Path" -> outputPath,
            "NumberOfModes" -> modes["NumberOfModes"]|>];
        , {i, Length[samples]}
    ];
    rows
];

MLModeDatasetHelp[] := Column[{
    Style["MLModeDataset.wl", Bold, 16],
    "MLModeDatasetSweep[csvPath, options] exports Gamma spectra from ML Sobol equilibria.",
    "Required: $ProjectRoot, loaded dependencies and MLSeedSweep.wl.",
    "Output: datasets/magnon_spectra/gamma/{train,validation,test}/sample_id.json.",
    "Options: FrequencyMax, ProjectionCutoff, SymmetryTolerance, MainBranchThreshold, Overwrite, SkipMissing, MaxSamples."
}];

Print["Loaded MLModeDataset.wl (run MLModeDatasetHelp[] for usage)."];