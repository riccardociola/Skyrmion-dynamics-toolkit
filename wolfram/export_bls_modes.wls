(* Export isotropic magnon modes to JSON with fine-grained progress logging.
   Run from the project root using Python's run_bls.py, or directly:
     wolframscript -local -file wolfram/export_bls_modes.wls INPUT OUTPUT
*)

ClearAll[blsProgress, blsTimed];
args = Rest[$ScriptCommandLine];
If[Length[args] != 2,
    Print["ERROR: Expected INPUT EQUILIBRIUM and OUTPUT JSON. Script args: ",
        InputForm[$ScriptCommandLine]];
    Exit[2]
];

inputPath = ExpandFileName[args[[1]]];
outputPath = ExpandFileName[args[[2]]];
scriptPath = ExpandFileName[First[$ScriptCommandLine]];
$ProjectRoot = DirectoryName[DirectoryName[scriptPath]];
progressPath = outputPath <> ".progress.txt";

If[!DirectoryQ[DirectoryName[outputPath]],
    CreateDirectory[DirectoryName[outputPath], CreateIntermediateDirectories -> True]
];

blsProgress[message_String] := Module[{line},
    line = DateString[{"ISODate", " ", "Time"}] <> " | " <> message;
    Export[progressPath, line, "Text"];
    Print[line];
];
SetAttributes[blsTimed, HoldRest];
blsTimed[description_String, computation_] := Module[{start, result},
    blsProgress[description <> " - STARTED"];
    start = AbsoluteTime[];
    result = computation;
    blsProgress[description <> " - DONE (" <>
        ToString[Round[AbsoluteTime[] - start, 0.1], InputForm] <> " s)"];
    result
];

blsProgress["Wolfram kernel running. Project root: " <> $ProjectRoot];
If[!FileExistsQ[inputPath],
    blsProgress["ERROR: Equilibrium not found: " <> inputPath];
    Exit[1]
];

modules = {"Geometry.wl", "FourierBasis.wl", "GroundState.wl",
    "DynamicsProjection.wl", "DynamicsMagnetic.wl", "ExportMagneticModes.wl"};

Do[
    modulePath = FileNameJoin[{$ProjectRoot, "wolfram", module}];
    If[!FileExistsQ[modulePath],
        blsProgress["ERROR: Missing module: " <> modulePath];
        Exit[1]
    ];
    blsTimed["Loading " <> module, Get[modulePath]];
, {module, modules}];

equilibrium = blsTimed["Reading equilibrium file", LoadMagneticEquilibrium[inputPath]];
If[!AssociationQ[equilibrium],
    blsProgress["ERROR: Equilibrium is not an Association"];
    Exit[1]
];
If[TrueQ[Lookup[equilibrium, "CubicAnisotropy", False]],
    blsProgress["ERROR: This exporter only handles isotropic dynamics"];
    Exit[1]
];

nCutoff = equilibrium["FourierCutoff"];
blsProgress["Fourier cutoff N=" <> ToString[nCutoff] <> "; representation=" <>
    ToString[equilibrium["MagneticRepresentation"]]];

(* Exactly the same operations, defaults, and output keys as
   PrepareMagneticDynamics, split into individual timing checkpoints. *)
nodes = blsTimed["Generating reciprocal nodes", DynamicsNodes[nCutoff]];
blsProgress["Reciprocal nodes=" <> ToString[Length[nodes]] <>
    "; full matrix dimension=" <> ToString[3 Length[nodes]]];

crossOperator = blsTimed["Building M0 x operator", BuildM0CrossOperator[equilibrium, nodes]];
If[!MatrixQ[crossOperator],
    blsProgress["ERROR: M0 x is not a matrix"];
    Exit[1]
];

transverseData = blsTimed[
    "Diagonalizing M0 x and building transverse projection",
    BuildMagneticTransverseBasis[crossOperator]
];
If[!AssociationQ[transverseData],
    blsProgress["ERROR: Transverse projection failed"];
    Exit[1]
];
blsProgress["Transverse dimension=" <> ToString[transverseData["TransverseDimension"]]];

constraintHessian = blsTimed[
    "Building dense constraint Hessian (Fourier convolution)",
    BuildConstraintHessian[equilibrium, nodes]
];
If[!MatrixQ[constraintHessian],
    blsProgress["ERROR: Constraint Hessian is not a matrix"];
    Exit[1]
];

constraintHessianTransverse = blsTimed[
    "Projecting constraint Hessian",
    ProjectMagneticOperator[transverseData, constraintHessian]
];
If[!MatrixQ[constraintHessianTransverse],
    blsProgress["ERROR: Projected Hessian is not a matrix"];
    Exit[1]
];

dynamics = <|
    "Equilibrium" -> equilibrium,
    "Nodes" -> nodes,
    "CrossOperatorFull" -> crossOperator,
    "TransverseData" -> transverseData,
    "ConstraintHessian" -> constraintHessian,
    "ConstraintHessianTransverse" -> constraintHessianTransverse,
    "FullDimension" -> 3 Length[nodes],
    "TransverseDimension" -> transverseData["TransverseDimension"]
|>;

modes = blsTimed[
    "Constructing LL matrix and solving modes at Gamma",
    SolveMagneticModes[dynamics, {0, 0, 0}, "FrequencyRange" -> {0, 3}]
];
If[!AssociationQ[modes],
    blsProgress["ERROR: Mode solution failed"];
    Exit[1]
];
blsProgress["Selected positive-frequency modes=" <> ToString[Length[modes["Frequencies"]]]];

partialPath = outputPath <> ".partial.json";
If[FileExistsQ[partialPath], DeleteFile[partialPath]];
result = blsTimed["Exporting complex Fourier arrays to JSON",
    ExportMagneticModes[dynamics, modes, partialPath]
];
If[result === $Failed || !FileExistsQ[partialPath],
    blsProgress["ERROR: Export failed"];
    Exit[1]
];
If[FileExistsQ[outputPath], DeleteFile[outputPath]];
RenameFile[partialPath, outputPath];
blsProgress["FINISHED: " <> outputPath];
Exit[0];
