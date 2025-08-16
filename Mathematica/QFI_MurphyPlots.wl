(* ::Package:: *)

(* QFI_MurphyPlots.wl
   Helpers for quotient matrices, spectral tilt, occupancy LDP, Murphy horizons,
   and action-betweenness. Mathematica (Wolfram Language) package. *)

BeginPackage["QFIHelpers`"]

QuotientMatrix::usage =
  "QuotientMatrix[kernelSamples, partition, measureWeights] -> P.
   kernelSamples: list of functions or arrays K[x, Aj] or conditional samples.
   partition: list of cell indicator functions or an index map cell[x].
   measureWeights: weights/probabilities for sampling x within a cell.";

TiltSpectralRadius::usage =
  "TiltSpectralRadius[P_, t_, S_] returns spectral radius rho(diag(exp(t*1_S)) . P).";

OccupancyRate::usage =
  "OccupancyRate[P_, S_, theta_] computes I_S(theta) = sup_t { t theta - log rho(t) }.";

CriticalHorizon::usage =
  "CriticalHorizon[P_, S_, theta_, alpha_:1.0, beta_:0.05] returns L_c ≈ log(1/beta)/(alpha I_S(theta)).";

MurphyPlot::usage =
  "MurphyPlot[P_, S_List, thetas_List, alpha_:1.0, beta_:0.05] plots L_c across theta for each S in S_List.";

ActionBetweenness::usage =
  "ActionBetweenness[P_, alpha_:1.0, weights_:Automatic] computes AB_alpha(k) over nodes k.";

Begin["`Private`"]

ClearAll[ToCellIndex]
ToCellIndex[partition_] := Module[{},
  Which[
    VectorQ[partition, IntegerQ], (* already indices *)
      Function[x, partition[[x]]],
    True, (* assume list of indicator functions *)
      Function[x, First@FirstPosition[partition /. f_ :> Boole[f[x] == 1], 1]]
  ]
];

(* Simple aggregator when we have a sample of x's per cell and access to K[x, Aj]. *)
ClearAll[QuotientMatrix]
QuotientMatrix[kernelSamples_, partition_, measureWeights_:Automatic] := Module[
  {cells, m, idx, reps, P},
  idx = ToCellIndex[partition];
  cells = DeleteDuplicates @ (idx /@ kernelSamples[[All, 1]]);
  m = Max[cells];
  reps = Table[Select[kernelSamples, idx[#[[1]]] == i &], {i, 1, m}];
  P = ConstantArray[0., {m, m}];
  Do[
    If[Length[reps[[i]]] > 0,
      With[{Xi = reps[[i, All, 1]], Ki = reps[[i, All, 2]]},
        (* Average K(x, Aj) over x ~ µ(· | Ai) *)
        Do[
          P[[i, j]] = Mean[Ki /. x_ :> x; (* ensure K[x, Aj] yields numeric *)
                           Table[Ki[[k]][j], {k, Length[Ki]}]
                         ],
          {j, 1, m}
        ];
      ];
    ],
    {i, 1, m}
  ];
  (* Row-stochastic normalization safeguard *)
  P = Normalize[#, Total] & /@ P;
  P
];

ClearAll[TiltSpectralRadius]
TiltSpectralRadius[P_?MatrixQ, t_?NumericQ, S_List] := Module[
  {m = Length[P], D, T, vals},
  D = DiagonalMatrix[Table[If[MemberQ[S, i], Exp[t], 1.0], {i, 1, Length[P]}]];
  T = D . P;
  Max[Abs[Eigenvalues[T]]]
];

ClearAll[OccupancyRate]
OccupancyRate[P_?MatrixQ, S_List, theta_?NumericQ] := Module[
  {f, tstar, val},
  f[t_] := t*theta - Log[TiltSpectralRadius[P, t, S]];
  (* Unimodal convex dual; bracket search *)
  tstar = t /. Quiet@FindMaximum[f[t], {t, 0}, Method -> "Newton"][[2, 1]];
  val = f[tstar];
  If[NumericQ[val], val, $Failed]
];

ClearAll[CriticalHorizon]
CriticalHorizon[P_?MatrixQ, S_List, theta_?NumericQ, alpha_:1.0, beta_:0.05] := Module[
  {I = OccupancyRate[P, S, theta]},
  If[I === $Failed || I <= 0, Infinity, Log[1./beta]/(alpha*I)]
];

ClearAll[MurphyPlot]
MurphyPlot[P_?MatrixQ, S_List, thetas_List, alpha_:1.0, beta_:0.05] := Module[
  {data, curves},
  data = Table[
    {θ, CriticalHorizon[P, Sset, θ, alpha, beta]},
    {Sset, S}, {θ, thetas}
  ];
  curves = Table[
    ListLinePlot[data[[k]], PlotLegends -> Placed[Automatic, {0.8, 0.8}],
      PlotRange -> All, AxesLabel -> {"θ", "L_c"}, PlotLabel -> Row[{"S = ", S[[k]]}]],
    {k, Length[S]}
  ];
  Show[curves]
];

(* Action-Betweenness via delta of optimal tilted costs *)
ClearAll[ActionBetweenness]
ActionBetweenness[P_?MatrixQ, alpha_:1.0, weights_:Automatic] := Module[
  {m = Length[P], W, costAvoid, costAll, ab},
  W = If[weights === Automatic, ConstantArray[1., {m, m}], weights];
  costAll[i_, j_] := Module[{S = {j}, θ = 1.0}, Log[TiltSpectralRadius[P, alpha, S]]];
  costAvoid[i_, j_, k_] := Module[{Pk = P},
    Pk[[All, k]] = 0.; Pk[[k, All]] = 0.; (* avoid k by zeroing transitions *)
    Log[TiltSpectralRadius[Pk, alpha, {j}]]
  ];
  ab = Table[
    Total@Flatten@Table[
      If[i == k || j == k || i == j, 0,
        W[[i, j]]*Max[0, costAll[i, j] - costAvoid[i, j, k]]
      ],
      {i, 1, m}, {j, 1, m}
    ],
    {k, 1, m}
  ];
  ab
];

End[]; (* `Private` *)
EndPackage[];