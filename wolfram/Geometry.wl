(* ================================================================ *)
(* Geometry.wl                                                      *)
(*                                                                  *)
(* Affine-coordinate tools for skyrmion-lattice calculations.       *)
(*                                                                  *)
(*                                                                  *)
(* Notebook / mathematical notation       Source code           *)
(* ─────────────────────────────────────────────────────        *)
(* η                                      eta                   *)
(* φ                                      phi                   *)
(* αR                                     alphaR                *)
(* gαβ                                    metric                *)
(* Gαβ                                    dualMetric            *)
(* ξ1                                     xi1                   *)
(* μ1                                     mu1                   *)
(* λ                                      lambda                *)
(* ================================================================ *)

ClearAll[
    rotAlpha,
    ee,
    e1, e2, e3,
    metric,
    metricDet,
    Vpuc,
    dualMetric,
    EE,
    E1, E2, E3,
    xi1, xi2,
    X, Y,
    x0, y0,
    DotC,
    DotCn,
    G,
    Gc,
    FT,
    toCRTc,
    toCRT,
    toAff,
    toAffc,
    GeometryHelp
];


(* --------------------------------------------------------------- *)
(* Assumptions                                                     *)
(* --------------------------------------------------------------- *)

$Assumptions =
    eta > 0 &&
    -Pi <= phi < Pi &&
    0 <= alphaR <= 2 Pi;


(* --------------------------------------------------------------- *)
(* Affine basis                                                   *)
(* --------------------------------------------------------------- *)

(* Rotation around the Cartesian z axis. *)
rotAlpha[alphaR_] := {
    {Cos[alphaR], -Sin[alphaR], 0},
    {Sin[alphaR],  Cos[alphaR], 0},
    {0, 0, 1}
};


(* Affine basis vectors {e1,e2,e3}, stored as rows. *)
ee[eta_, phi_, alphaR_] := {
    {
        Cos[Pi/6 + phi/2],
        -Sin[Pi/6 + phi/2],
        0
    },
    {
        eta Cos[Pi/6 + phi/2],
        eta Sin[Pi/6 + phi/2],
        0
    },
    {
        0,
        0,
        1
    }
} . rotAlpha[alphaR];


(* Individual affine basis vectors. *)
e1[eta_, phi_, alphaR_] := ee[eta, phi, alphaR][[1]];
e2[eta_, phi_, alphaR_] := ee[eta, phi, alphaR][[2]];
e3 := {0, 0, 1};


(* --------------------------------------------------------------- *)
(* Metric and dual basis                                          *)
(* --------------------------------------------------------------- *)

(* Metric tensor g_ab = e_a . e_b. *)
metric[eta_, phi_, alphaR_] :=
    Outer[
        Dot,
        ee[eta, phi, alphaR],
        ee[eta, phi, alphaR],
        1
    ] // FullSimplify;


(* Determinant of the metric tensor. *)
metricDet[eta_, phi_, alphaR_] :=
    Det[metric[eta, phi, alphaR]] // FullSimplify;


(* Primitive unit-cell volume (area for the 2D magnetic lattice). *)
Vpuc[eta_, phi_, alphaR_] :=
    Sqrt[metricDet[eta, phi, alphaR]] // FullSimplify;


(* Dual metric tensor g^ab. *)
dualMetric[eta_, phi_, alphaR_] :=
    Inverse[metric[eta, phi, alphaR]];


(* Dual basis vectors, stored as rows. *)
EE[eta_, phi_, alphaR_] :=
    Transpose[
        Inverse[
            ee[eta, phi, alphaR]
        ]
    ];


(* Individual dual basis vectors. *)
E1[eta_, phi_, alphaR_] := EE[eta, phi, alphaR][[1]];
E2[eta_, phi_, alphaR_] := EE[eta, phi, alphaR][[2]];
E3 := {0, 0, 1};


(* --------------------------------------------------------------- *)
(* Coordinate transformations                                     *)
(* --------------------------------------------------------------- *)

(* Affine coordinates (xi1, xi2) corresponding to Cartesian (x,y). *)
xi1[x_, y_, eta_, phi_, alphaR_] :=
    {x, y, 0} . E1[eta, phi, alphaR];

xi2[x_, y_, eta_, phi_, alphaR_] :=
    {x, y, 0} . E2[eta, phi, alphaR];


(* Cartesian x coordinate from affine coordinates (xi1, xi2). *)
X[xi1_, xi2_, eta_, phi_, alphaR_] :=
    (
        xi1 e1[eta, phi, alphaR] +
        xi2 e2[eta, phi, alphaR]
    )[[1]];


(* Cartesian y coordinate from affine coordinates (xi1, xi2). *)
Y[xi1_, xi2_, eta_, phi_, alphaR_] :=
    (
        xi1 e1[eta, phi, alphaR] +
        xi2 e2[eta, phi, alphaR]
    )[[2]];


(* Center of the parallelogram primitive unit cell. *)
x0[eta_, phi_, alphaR_] :=
    (
        (e1[eta, phi, alphaR] + e2[eta, phi, alphaR])/2
    )[[1]];

y0[eta_, phi_, alphaR_] :=
    (
        (e1[eta, phi, alphaR] + e2[eta, phi, alphaR])/2
    )[[2]];


(* --------------------------------------------------------------- *)
(* Scalar products                                                *)
(* --------------------------------------------------------------- *)

(* Scalar product of two vectors given in covariant components. *)
DotC[v1_List, v2_List, eta_, phi_, alphaR_] :=
    v1 . dualMetric[eta, phi, alphaR] . v2;


(* Scalar product of two vectors given in contravariant components. *)
DotCn[v1_List, v2_List, eta_, phi_, alphaR_] :=
    v1 . metric[eta, phi, alphaR] . v2;


(* --------------------------------------------------------------- *)
(* Reciprocal lattice                                             *)
(* --------------------------------------------------------------- *)

(* Covariant affine components of reciprocal vector G(n1,n2). *)
G[n1_, n2_] :=
    2 Pi {n1, n2, 0};


(* Contravariant affine components of reciprocal vector G(n1,n2). *)
Gc[n1_Integer, n2_Integer, eta_, phi_, alphaR_] :=
    dualMetric[eta, phi, alphaR] . G[n1, n2];


(* --------------------------------------------------------------- *)
(* Fourier transform                                              *)
(* --------------------------------------------------------------- *)

(* Fourier transform over one primitive affine unit cell. *)
FT[f_, {n1_, n2_}] :=
    NIntegrate[
        f[xi1, xi2] Exp[
            -I G[n1, n2] . {xi1, xi2, 0}
        ],
        {xi1, 0, 1},
        {xi2, 0, 1}
    ];


(* --------------------------------------------------------------- *)
(* Basis conversions                                              *)
(* --------------------------------------------------------------- *)

(* Contravariant affine components -> Cartesian components. *)
toCRTc[vc_List, eta_, phi_, alphaR_] :=
    Inner[
        Times,
        vc,
        ee[eta, phi, alphaR]
    ];


(* Covariant affine components -> Cartesian components. *)
toCRT[v_List, eta_, phi_, alphaR_] :=
    Inner[
        Times,
        v,
        EE[eta, phi, alphaR]
    ];


(* Cartesian components -> covariant affine components. *)
toAff[v_List, eta_, phi_, alphaR_] :=
    (# . v) & /@ ee[eta, phi, alphaR];


(* Cartesian components -> contravariant affine components. *)
toAffc[v_List, eta_, phi_, alphaR_] :=
    (# . v) & /@ EE[eta, phi, alphaR];


(* --------------------------------------------------------------- *)
(* Help                                                            *)
(* --------------------------------------------------------------- *)

GeometryHelp[] :=
    Column[
        {
            Style["Geometry.wl", Bold, 16],
            "Affine geometry, reciprocal-lattice tools and coordinate transformations.",
            Grid[
                {
                    {"Function", "Description"},
                    {"ee / EE", "Return the affine basis and dual basis vectors."},
                    {"metric / dualMetric", "Return the affine metric and inverse metric tensors."},
                    {"Vpuc", "Return the primitive unit-cell area."},
                    {"xi1 / xi2", "Convert Cartesian coordinates to affine coordinates."},
                    {"X / Y", "Convert affine coordinates to Cartesian coordinates."},
                    {"DotC / DotCn", "Scalar products for covariant / contravariant affine components."},
                    {"G / Gc", "Return reciprocal-lattice vectors in affine components."},
                    {"toCRT / toCRTc", "Convert affine vector components to Cartesian components."},
                    {"toAff / toAffc", "Convert Cartesian vector components to affine components."},
                    {"FT", "Numerically compute a Fourier coefficient over one primitive cell."}
                },
                Frame -> All,
                Alignment -> Left
            ]
        }
    ];

Print[GeometryHelp[]];