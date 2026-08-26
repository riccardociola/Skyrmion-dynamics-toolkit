(* ================================================================ *)
(* Geometry.wl                                                       *)
(*                                                                  *)
(* Affine-coordinate tools for skyrmion-lattice calculations.       *)
(* The basis vectors are treated as rows of the affine basis matrix. *)
(* ================================================================ *)


ClearAll[
    rotα,
    ee,
    e1, e2, e3,
    gαβ,
    g,
    Vpuc,
    Gαβ,
    EE,
    E1, E2, E3,
    ξ1, ξ2,
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
    toAffc
];


(* --------------------------------------------------------------- *)
(* Assumptions                                                     *)
(* --------------------------------------------------------------- *)

$Assumptions =
    η > 0 &&
    -Pi <= φ < Pi &&
    0 <= αR <= 2 Pi;


(* --------------------------------------------------------------- *)
(* Affine basis                                                   *)
(* --------------------------------------------------------------- *)

(* Rotation around the Cartesian z axis. *)
rotα[αR_] := {
    {Cos[αR], -Sin[αR], 0},
    {Sin[αR],  Cos[αR], 0},
    {0,        0,       1}
};


(* Affine basis vectors {e1,e2,e3}, stored as rows. *)
ee[η_, φ_, αR_] := {
    {
        Cos[Pi/6 + φ/2],
        -Sin[Pi/6 + φ/2],
        0
    },
    {
        η Cos[Pi/6 + φ/2],
        η Sin[Pi/6 + φ/2],
        0
    },
    {
        0,
        0,
        1
    }
} . rotα[αR];


(* Individual affine basis vectors. *)
e1[η_, φ_, αR_] := ee[η, φ, αR][[1]];
e2[η_, φ_, αR_] := ee[η, φ, αR][[2]];
e3 := {0, 0, 1};


(* --------------------------------------------------------------- *)
(* Metric and dual basis                                          *)
(* --------------------------------------------------------------- *)

(* Metric tensor g_ab = e_a . e_b. *)
gαβ[η_, φ_, αR_] :=
    Outer[
        Dot,
        ee[η, φ, αR],
        ee[η, φ, αR],
        1
    ] // FullSimplify;


(* Determinant of the metric tensor. *)
g[η_, φ_, αR_] :=
    Det[gαβ[η, φ, αR]] // FullSimplify;


(* Primitive unit-cell volume (area for the 2D magnetic lattice). *)
Vpuc[η_, φ_, αR_] :=
    Sqrt[g[η, φ, αR]] // FullSimplify;


(* Dual metric tensor g^ab. *)
Gαβ[η_, φ_, αR_] :=
    Inverse[gαβ[η, φ, αR]];


(* Dual basis vectors, stored as rows. *)
EE[η_, φ_, αR_] :=
    Transpose[
        Inverse[
            ee[η, φ, αR]
        ]
    ];


(* Individual dual basis vectors. *)
E1[η_, φ_, αR_] := EE[η, φ, αR][[1]];
E2[η_, φ_, αR_] := EE[η, φ, αR][[2]];
E3 := {0, 0, 1};


(* --------------------------------------------------------------- *)
(* Coordinate transformations                                     *)
(* --------------------------------------------------------------- *)

(* Affine coordinates ξ1, ξ2 corresponding to Cartesian (x,y). *)
ξ1[x_, y_, η_, φ_, αR_] :=
    {x, y, 0} . E1[η, φ, αR];

ξ2[x_, y_, η_, φ_, αR_] :=
    {x, y, 0} . E2[η, φ, αR];


(* Cartesian x coordinate from affine coordinates (ξ1,ξ2). *)
X[ξ1_, ξ2_, η_, φ_, αR_] :=
    (
        ξ1 e1[η, φ, αR] +
        ξ2 e2[η, φ, αR]
    )[[1]];


(* Cartesian y coordinate from affine coordinates (ξ1,ξ2). *)
Y[ξ1_, ξ2_, η_, φ_, αR_] :=
    (
        ξ1 e1[η, φ, αR] +
        ξ2 e2[η, φ, αR]
    )[[2]];


(* Center of the parallelogram primitive unit cell. *)
x0[η_, φ_, αR_] :=
    (
        (e1[η, φ, αR] + e2[η, φ, αR])/2
    )[[1]];

y0[η_, φ_, αR_] :=
    (
        (e1[η, φ, αR] + e2[η, φ, αR])/2
    )[[2]];


(* --------------------------------------------------------------- *)
(* Scalar products                                                *)
(* --------------------------------------------------------------- *)

(* Scalar product of two vectors given in covariant components. *)
DotC[v1_List, v2_List, η_, φ_, αR_] :=
    TensorContract[
        Gαβ[η, φ, αR] \[TensorProduct] v1 \[TensorProduct] v2,
        {{1, 3}, {2, 4}}
    ];


(* Scalar product of two vectors given in contravariant components. *)
DotCn[v1_List, v2_List, η_, φ_, αR_] :=
    TensorContract[
        gαβ[η, φ, αR] \[TensorProduct] v1 \[TensorProduct] v2,
        {{1, 3}, {2, 4}}
    ];


(* --------------------------------------------------------------- *)
(* Reciprocal lattice                                             *)
(* --------------------------------------------------------------- *)

(* Covariant affine components of reciprocal vector G(n1,n2). *)
G[n1_, n2_] :=
    2 Pi {n1, n2, 0};


(* Contravariant affine components of reciprocal vector G(n1,n2). *)
Gc[n1_Integer, n2_Integer, η_, φ_, αR_] :=
    TensorContract[
        Gαβ[η, φ, αR] \[TensorProduct] G[n1, n2],
        {{2, 3}}
    ];


(* --------------------------------------------------------------- *)
(* Fourier transform                                              *)
(* --------------------------------------------------------------- *)

(* Fourier transform over one primitive affine unit cell. *)
FT[f_, {n1_, n2_}] :=
    NIntegrate[
        f[ξ1, ξ2] Exp[
            -I G[n1, n2] . {ξ1, ξ2, 0}
        ],
        {ξ1, 0, 1},
        {ξ2, 0, 1}
    ];


(* --------------------------------------------------------------- *)
(* Basis conversions                                              *)
(* --------------------------------------------------------------- *)

(* Contravariant affine components -> Cartesian components. *)
toCRTc[vc_List, η_, φ_, αR_] :=
    Inner[
        Times,
        vc,
        ee[η, φ, αR]
    ];


(* Covariant affine components -> Cartesian components. *)
toCRT[v_List, η_, φ_, αR_] :=
    Inner[
        Times,
        v,
        EE[η, φ, αR]
    ];


(* Cartesian components -> covariant affine components. *)
toAff[v_List, η_, φ_, αR_] :=
    (# . v) & /@ ee[η, φ, αR];


(* Cartesian components -> contravariant affine components. *)
toAffc[v_List, η_, φ_, αR_] :=
    (# . v) & /@ EE[η, φ, αR];