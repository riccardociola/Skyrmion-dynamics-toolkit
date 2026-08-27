(* ================================================================ *)
(* GroundState.wl                                                   *)
(*                                                                  *)
(* Ground-state / metastable-equilibrium minimization routines.      *)
(* ================================================================ *)


ClearAll[
    FindMagneticEquilibrium
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