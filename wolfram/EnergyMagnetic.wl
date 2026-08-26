(* ================================================================ *)
(* EnergyMagnetic.wl                                                *)
(*                                                                  *)
(* Magnetic contributions to the free-energy functional.            *)
(* The active Fourier representation and cutoff are provided by     *)
(* FourierBasis.wl                                                  *)
(* ================================================================ *)


ClearAll[
    ExchangeEnergy,
    BlochDMIEnergy,
    NeelDMIEnergy,
    AntiDMIEnergy,
    ZeemanEnergy,
    MagnetostaticInfiniteEnergy,
    MagnetostaticFilmEnergy,
    MagnetostaticPeriodicEnergy,
    ConstraintEnergy
];


(* --------------------------------------------------------------- *)
(* Exchange interaction                                           *)
(* --------------------------------------------------------------- *)

(* Exchange energy  *)
ExchangeEnergy[a_, eta_, phi_, alphaR_] :=
    1/(2 a^2) Sum[
        DotC[
            G[n1, n2],
            G[n1, n2],
            eta, phi, alphaR
        ]
        DotC[
            M[n1, n2],
            M[-n1, -n2],
            eta, phi, alphaR
        ],
        {n1, -$FourierCutoff, $FourierCutoff},
        {n2, -$FourierCutoff, $FourierCutoff}
    ];


(* --------------------------------------------------------------- *)
(* Dzyaloshinskii-Moriya interaction                              *)
(* --------------------------------------------------------------- *)

(* Bulk/Bloch DMI  *)
BlochDMIEnergy[a_, eta_, phi_, alphaR_] :=
    -(I/(a Vpuc[eta, phi, alphaR])) Sum[
        TensorContract[
            TensorProduct[
                LeviCivitaTensor[3],
                M[n1, n2],
                G[n1, n2],
                M[-n1, -n2]
            ],
            {{1, 4}, {2, 5}, {3, 6}}
        ],
        {n1, -$FourierCutoff, $FourierCutoff},
        {n2, -$FourierCutoff, $FourierCutoff}
    ];


(* Interfacial/Neel DMI  *)
NeelDMIEnergy[a_, eta_, phi_, alphaR_] :=
    (2 I)/a Sum[
        M[-n1, -n2][[3]]
        DotC[
            M[n1, n2],
            G[n1, n2],
            eta, phi, alphaR
        ],
        {n1, -$FourierCutoff, $FourierCutoff},
        {n2, -$FourierCutoff, $FourierCutoff}
    ];


(* DMI supporting antiskyrmions  *)
AntiDMIEnergy[a_, eta_, phi_, alphaR_] :=
    I/(a Vpuc[eta, phi, alphaR]) Sum[
        TensorContract[
            TensorProduct[
                LeviCivitaTensor[3],

                (
                    (
                        #[[1]] toCRT[
                            G[n1, n2],
                            eta, phi, alphaR
                        ][[1]]
                        -
                        #[[2]] toCRT[
                            G[n1, n2],
                            eta, phi, alphaR
                        ][[2]]
                    ) & /@
                    ee[eta, phi, alphaR]
                ),

                M[-n1, -n2],
                M[n1, n2]
            ],
            {{1, 4}, {2, 5}, {3, 6}}
        ],
        {n1, -$FourierCutoff, $FourierCutoff},
        {n2, -$FourierCutoff, $FourierCutoff}
    ];


(* --------------------------------------------------------------- *)
(* Zeeman interaction                                             *)
(* --------------------------------------------------------------- *)

(* Zeeman energy for a magnetic field along the z direction  *)
ZeemanEnergy[b_] :=
    -b M[0, 0][[3]];


(* --------------------------------------------------------------- *)
(* Magnetostatic interaction                                      *)
(* --------------------------------------------------------------- *)

(*
    Magnetostatic energy for an infinitely extended sample.

    xi:
        dimensionless magnetostatic coupling.

    NijCartesian:
        demagnetization tensor in Cartesian components.
        A zero tensor corresponds to a boundless sample.
*)
MagnetostaticInfiniteEnergy[
    xi_,
    NijCartesian_List,
    eta_,
    phi_,
    alphaR_
] :=
    Module[
        {NijCovariant},

        (* Convert Cartesian demagnetization tensor to covariant affine components. *)
        NijCovariant =
            ee[eta, phi, alphaR] .
            NijCartesian .
            Transpose[ee[eta, phi, alphaR]];

        xi Sum[

            If[
                n1 == 0 && n2 == 0,

                TensorContract[
                    TensorProduct[
                        NijCovariant,
                        dualMetric[eta, phi, alphaR],
                        dualMetric[eta, phi, alphaR],
                        M[0, 0],
                        M[0, 0]
                    ],
                    {{1, 3}, {2, 5}, {4, 7}, {6, 8}}
                ],

                DotC[
                    M[n1, n2],
                    G[n1, n2],
                    eta, phi, alphaR
                ]
                DotC[
                    M[-n1, -n2],
                    G[n1, n2],
                    eta, phi, alphaR
                ]
                /
                DotC[
                    G[n1, n2],
                    G[n1, n2],
                    eta, phi, alphaR
                ]
            ],

            {n1, -$FourierCutoff, $FourierCutoff},
            {n2, -$FourierCutoff, $FourierCutoff}
        ]
    ];


(*
    Magnetostatic energy of a film with finite thickness L 
*)
MagnetostaticFilmEnergy[
    xi_,
    a_,
    L_,
    eta_,
    phi_,
    alphaR_
] :=
    Module[
        {chi},

        xi Sum[

            If[
                n1 == 0 && n2 == 0,

                M[0, 0][[3]]^2,

                chi =
                    (L/a) Sqrt[
                        DotC[
                            G[n1, n2],
                            G[n1, n2],
                            eta, phi, alphaR
                        ]
                    ];

                DotC[
                    M[n1, n2],
                    G[n1, n2],
                    eta, phi, alphaR
                ]
                DotC[
                    M[-n1, -n2],
                    G[n1, n2],
                    eta, phi, alphaR
                ]
                /
                DotC[
                    G[n1, n2],
                    G[n1, n2],
                    eta, phi, alphaR
                ]
                *
                (Exp[-chi] - 1 + chi)/chi

                +

                M[n1, n2][[3]]
                M[-n1, -n2][[3]]
                (1 - Exp[-chi])/chi
            ],

            {n1, -$FourierCutoff, $FourierCutoff},
            {n2, -$FourierCutoff, $FourierCutoff}
        ]
    ];


(*
    Magnetostatic energy of an infinite periodic multilayer.

    L:
        thickness of the ferromagnetic layer.

    h:
        thickness of the spacer.
*)
MagnetostaticPeriodicEnergy[
    xi_,
    a_,
    L_,
    h_,
    eta_,
    phi_,
    alphaR_
] :=
    Module[
        {chi, sigma},

        xi Sum[

            If[
                n1 == 0 && n2 == 0,

                M[0, 0][[3]]^2 (1 - 1/(1 + h/L)),

                chi =
                    (L/a) Sqrt[
                        DotC[
                            G[n1, n2],
                            G[n1, n2],
                            eta, phi, alphaR
                        ]
                    ];

                sigma =
                    (h/a) Sqrt[
                        DotC[
                            G[n1, n2],
                            G[n1, n2],
                            eta, phi, alphaR
                        ]
                    ];

                DotC[
                    M[n1, n2],
                    G[n1, n2],
                    eta, phi, alphaR
                ]
                DotC[
                    M[-n1, -n2],
                    G[n1, n2],
                    eta, phi, alphaR
                ]
                /
                DotC[
                    G[n1, n2],
                    G[n1, n2],
                    eta, phi, alphaR
                ]
                *
                (
                    (Exp[-chi] - 1 + chi)/chi
                    +
                    (4/chi)
                    Sinh[chi/2]^2/
                    (Exp[chi + sigma] - 1)
                )

                +

                M[n1, n2][[3]]
                M[-n1, -n2][[3]]
                *
                (
                    (1 - Exp[-chi])/chi
                    -
                    (4/chi)
                    Sinh[chi/2]^2/
                    (Exp[chi + sigma] - 1)
                )
            ],

            {n1, -$FourierCutoff, $FourierCutoff},
            {n2, -$FourierCutoff, $FourierCutoff}
        ]
    ];


(* --------------------------------------------------------------- *)
(* Magnetization-length constraint                                *)
(* --------------------------------------------------------------- *)

(*
    Quartic penalty enforcing |m(r)| = 1.
*)
ConstraintEnergy[lambda_, eta_, phi_, alphaR_] :=
    lambda (
        1
        -
        2 Sum[
            DotC[
                M[n1, n2],
                M[-n1, -n2],
                eta, phi, alphaR
            ],
            {n1, -$FourierCutoff, $FourierCutoff},
            {n2, -$FourierCutoff, $FourierCutoff}
        ]
        +
        ParallelSum[
            DotC[
                M[n1, n2],
                M[p1, p2],
                eta, phi, alphaR
            ]
            DotC[
                M[s1, s2],
                M[
                    -n1 - p1 - s1,
                    -n2 - p2 - s2
                ],
                eta, phi, alphaR
            ],

            {n1, -$FourierCutoff, $FourierCutoff},
            {n2, -$FourierCutoff, $FourierCutoff},

            {p1, -$FourierCutoff, $FourierCutoff},
            {p2, -$FourierCutoff, $FourierCutoff},

            {s1, -$FourierCutoff, $FourierCutoff},
            {s2, -$FourierCutoff, $FourierCutoff}
        ]
    );