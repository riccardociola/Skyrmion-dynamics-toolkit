(* ================================================================ *)
(* FourierBasis.wl                                                  *)
(*                                                                  *)
(* Fourier-space representations of magnetic and elastic fields.    *)
(* Handles reciprocal-space cutoff, reality conditions and field    *)
(* symmetry.                                                        *)
(* ================================================================ *)


ClearAll[
    hexagonNodeQ,
    independentNodeQ,
    independentNodes,
    muR, 
    muI,
    uR, 
    uI,
    MSixFold,
    MGeneral,
    UGeneral,
    M,
    U,
    MagneticVariables,
    ElasticVariables,
    ConfigureFourierBasis,
    MagneticFourierComponent,
    FourierBasisHelp
];


(* --------------------------------------------------------------- *)
(* Reciprocal-space domain                                        *)
(* --------------------------------------------------------------- *)

(* True if (n1,n2) lies inside the finite hexagonal cutoff. *)
hexagonNodeQ[n1_Integer, n2_Integer, cutoff_Integer] :=
    Abs[n1] <= cutoff &&
    Abs[n2] <= cutoff &&
    Abs[n2 - n1] <= cutoff;


(* Select one representative from each pair G and -G. *)
independentNodeQ[n1_Integer, n2_Integer] :=
    n2 > 0 ||
    (n2 == 0 && n1 >= 0);


(* List of independent reciprocal-space nodes inside the hexagon. *)
independentNodes[cutoff_Integer] :=
    Select[
        Join[
            Flatten[
                Table[
                    {n1, n2},
                    {n1, -cutoff, cutoff},
                    {n2, 1, cutoff}
                ],
                1
            ],
            Table[
                {n1, 0},
                {n1, 0, cutoff}
            ]
        ],
        hexagonNodeQ[#[[1]], #[[2]], cutoff] &
    ];


(* --------------------------------------------------------------- *)
(* Sixfold magnetic representation                                *)
(* --------------------------------------------------------------- *)

(*
    High-symmetry magnetic Fourier field:

        M_G = {i muI1_G, i muI2_G, muR3_G}

    Reality condition M_-G = Conjugate[M_G] implies

        muI1_-G = -muI1_G
        muI2_-G = -muI2_G
        muR3_-G =  muR3_G

    At G = 0 only muR3 is allowed.
*)

MSixFold[n1_Integer, n2_Integer, cutoff_Integer] :=
    Which[

        (* Outside reciprocal-space cutoff. *)
        !hexagonNodeQ[n1, n2, cutoff],
        {0, 0, 0},

        (* Uniform Fourier component. *)
        n1 == 0 && n2 == 0,
        {
            0,
            0,
            muR[3, 0, 0]
        },

        (* Independent half of reciprocal space. *)
        independentNodeQ[n1, n2],
        {
            I muI[1, n1, n2],
            I muI[2, n1, n2],
            muR[3, n1, n2]
        },

        (* Dependent half obtained from the reality condition. *)
        True,
        {
            -I muI[1, -n1, -n2],
            -I muI[2, -n1, -n2],
            muR[3, -n1, -n2]
        }
    ];


(* --------------------------------------------------------------- *)
(* General magnetic representation                                *)
(* --------------------------------------------------------------- *)

(*
    General complex magnetic Fourier field:

        M_G =
        {
            muR1_G + i muI1_G,
            muR2_G + i muI2_G,
            muR3_G + i muI3_G
        }

    Reality condition:
        M_-G = Conjugate[M_G]

    The zero-mode therefore has only real components.
*)

MGeneral[n1_Integer, n2_Integer, cutoff_Integer] :=
    Which[

        !hexagonNodeQ[n1, n2, cutoff],
        {0, 0, 0},

        n1 == 0 && n2 == 0,
        {
            muR[1, 0, 0],
            muR[2, 0, 0],
            muR[3, 0, 0]
        },

        independentNodeQ[n1, n2],
        {
            muR[1, n1, n2] + I muI[1, n1, n2],
            muR[2, n1, n2] + I muI[2, n1, n2],
            muR[3, n1, n2] + I muI[3, n1, n2]
        },

        True,
        {
            muR[1, -n1, -n2] - I muI[1, -n1, -n2],
            muR[2, -n1, -n2] - I muI[2, -n1, -n2],
            muR[3, -n1, -n2] - I muI[3, -n1, -n2]
        }
    ];


(* --------------------------------------------------------------- *)
(* General elastic representation                                 *)
(* --------------------------------------------------------------- *)

(*
    Complex Fourier representation of the displacement field U.

    Reality condition:
        U_-G = Conjugate[U_G]
*)

UGeneral[n1_Integer, n2_Integer, cutoff_Integer] :=
    Which[

        !hexagonNodeQ[n1, n2, cutoff],
        {0, 0, 0},

        n1 == 0 && n2 == 0,
        {
            uR[1, 0, 0],
            uR[2, 0, 0],
            uR[3, 0, 0]
        },

        independentNodeQ[n1, n2],
        {
            uR[1, n1, n2] + I uI[1, n1, n2],
            uR[2, n1, n2] + I uI[2, n1, n2],
            uR[3, n1, n2] + I uI[3, n1, n2]
        },

        True,
        {
            uR[1, -n1, -n2] - I uI[1, -n1, -n2],
            uR[2, -n1, -n2] - I uI[2, -n1, -n2],
            uR[3, -n1, -n2] - I uI[3, -n1, -n2]
        }
    ];


(* --------------------------------------------------------------- *)
(* Independent optimization variables                             *)
(* --------------------------------------------------------------- *)

(* Independent magnetic variables for the sixfold representation. *)
MagneticVariables[cutoff_Integer, "SixFold"] :=
    Flatten[
        Map[
            Function[node,
                If[
                    node == {0, 0},

                    {
                        muR[3, 0, 0]
                    },

                    {
                        muI[1, node[[1]], node[[2]]],
                        muI[2, node[[1]], node[[2]]],
                        muR[3, node[[1]], node[[2]]]
                    }
                ]
            ],
            independentNodes[cutoff]
        ]
    ];


(* Independent magnetic variables for the general representation. *)
MagneticVariables[cutoff_Integer, "General"] :=
    Flatten[
        Map[
            Function[node,
                If[
                    node == {0, 0},

                    {
                        muR[1, 0, 0],
                        muR[2, 0, 0],
                        muR[3, 0, 0]
                    },

                    {
                        muR[1, node[[1]], node[[2]]],
                        muI[1, node[[1]], node[[2]]],

                        muR[2, node[[1]], node[[2]]],
                        muI[2, node[[1]], node[[2]]],

                        muR[3, node[[1]], node[[2]]],
                        muI[3, node[[1]], node[[2]]]
                    }
                ]
            ],
            independentNodes[cutoff]
        ]
    ];


(* Independent elastic-field variables. *)
ElasticVariables[cutoff_Integer] :=
    Flatten[
        Map[
            Function[node,
                If[
                    node == {0, 0},

                    {
                        uR[1, 0, 0],
                        uR[2, 0, 0],
                        uR[3, 0, 0]
                    },

                    {
                        uR[1, node[[1]], node[[2]]],
                        uI[1, node[[1]], node[[2]]],

                        uR[2, node[[1]], node[[2]]],
                        uI[2, node[[1]], node[[2]]],

                        uR[3, node[[1]], node[[2]]],
                        uI[3, node[[1]], node[[2]]]
                    }
                ]
            ],
            independentNodes[cutoff]
        ]
    ];


(* --------------------------------------------------------------- *)
(* Active Fourier basis                                           *)
(* --------------------------------------------------------------- *)

ConfigureFourierBasis::rep =
    "Unknown magnetic representation `1`. Use \"SixFold\" or \"General\".";


ConfigureFourierBasis[
    cutoff_Integer?NonNegative,
    representation_String : "SixFold"
] :=
    Module[{},
        If[
            !MemberQ[{"SixFold", "General"}, representation],
            Message[ConfigureFourierBasis::rep, representation];
            Return[$Failed];
        ];

        $FourierCutoff = cutoff;
        $MagneticRepresentation = representation;
    ];


MagneticFourierComponent[
    n1_Integer,
    n2_Integer,
    cutoff_Integer?NonNegative,
    representation_String
] :=
    Switch[
        representation,

        "SixFold",
        MSixFold[
            n1,
            n2,
            cutoff
        ],

        "General",
        MGeneral[
            n1,
            n2,
            cutoff
        ],

        _,
        Message[
            MagneticFourierComponent::representation,
            representation
        ];

        $Failed
    ];


MagneticFourierComponent::representation =
    "Unknown magnetic representation `1`.";


(* Active magnetic Fourier field. *)
M[n1_Integer, n2_Integer] :=
    MagneticFourierComponent[n1, n2, $FourierCutoff, $MagneticRepresentation];


(* Elastic displacement field. *)
U[n1_Integer, n2_Integer] :=
    UGeneral[n1, n2, $FourierCutoff];


(* --------------------------------------------------------------- *)
(* Help                                                            *)
(* --------------------------------------------------------------- *)

FourierBasisHelp[] :=
    Column[
        {
            Style["FourierBasis.wl", Bold, 16],
            "Fourier-space field representations, reciprocal-space cutoff and reality constraints.",
            Grid[
                {
                    {"Function", "Description"},
                    {"ConfigureFourierBasis", "Set the reciprocal cutoff and magnetic representation (SixFold or General)."},
                    {"M", "Return the active magnetic Fourier component M[n1,n2]."},
                    {"U", "Return the elastic displacement Fourier component U[n1,n2]."},
                    {"MagneticVariables", "Return the independent magnetic optimization variables."},
                    {"ElasticVariables", "Return the independent elastic optimization variables."},
                    {"MagneticFourierComponent", "Return a magnetic Fourier component for an explicit cutoff and representation."},
                    {"independentNodes", "Return the independent reciprocal-space nodes inside the hexagonal cutoff."},
                    {"hexagonNodeQ", "Test whether a reciprocal node lies inside the finite hexagonal cutoff."}
                },
                Frame -> All,
                Alignment -> Left
            ]
        }
    ];

Print[FourierBasisHelp[]];