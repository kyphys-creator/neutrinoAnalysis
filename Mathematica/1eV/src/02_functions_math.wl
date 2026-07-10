(* ========================================================================== *)
(*  02_functions_math  (1 eV threshold)
*)
(*  General helpers: plot styles, box smearing, linear interpolation, trapezoid integration, data->flux builders, sampling integrals and Ev grids. From Functions.nb, ExpRate.nb and 1_Main.nb (final definitions kept where the legacy chain redefined a function).
*)
(* ========================================================================== *)

(* ---- plot styles ---- *)
clr[x_, y_] := ColorData[x, "ColorList"][[y]];

Inword[word_][a_][x_, y_][color_] := Inset[Style[word, a, color, FontFamily -> "Times"], ImageScaled[{x, y}],
    {Center, Center}];

PlotfS1[f_][x1_, x2_, xrange_, yrange_][B_, L_, T_, fsize_][scaling_, ImSize_, Styling_] :=
   Plot[f, {x, x1, x2}, PlotRange -> {xrange, yrange}, ImageSize -> ImSize,
    ScalingFunctions -> Switch[scaling, "N", None, "L", "Log", "LL", {"Log", "Log"}], Frame -> True,
    LabelStyle -> Directive[FontFamily -> "Times", Black, fsize],
    FrameLabel -> {Style[B, fsize, SingleLetterItalics -> False], Style[L, fsize], Style[T, fsize + 2]}, Styling];

BinPlotfS1[xlist_, ylist_, xrange_, yrange_][B_, L_, T_, fsize_][scaling_, ImSize_, Styling_] :=
  Module[{leng = Length[ylist], flist}, flist = ConstantArray[0, 2*leng];
    For[i = 1, i <= leng, i++, flist[[2*i - 1]] = {xlist[[i]], ylist[[i]]};
      flist[[2*i]] = {xlist[[i + 1]], ylist[[i]]}; ]; ListLinePlot[flist, PlotRange -> {xrange, yrange},
     ImageSize -> ImSize, ScalingFunctions -> Switch[scaling, "N", None, "L", "Log", "LL", {"Log", "Log"}],
     Frame -> True, LabelStyle -> Directive[FontFamily -> "Times", Black, fsize], PlotLabel -> Style[T, fsize + 2],
     FrameLabel -> {Style[B, fsize, SingleLetterItalics -> False], Style[L, fsize]}, Styling]]

(* ---- energy resolution: box smearing ---- *)
Gbox[x_, mu_, sig_] := (1/(2*sig))*UnitBox[(x - mu)/(2*sig)];

IntGboxmu[ER_, e1_, e2_, sig_] :=
  (1/(2*sig))*(((e2 - ER + sig)*UnitStep[e2 - ER + sig] - (e2 - ER - sig)*UnitStep[e2 - ER - sig]) -
    ((e1 - ER + sig)*UnitStep[e1 - ER + sig] - (e1 - ER - sig)*UnitStep[e1 - ER - sig]))

(* ---- linear interpolation / integration ---- *)
LIntpl1D[x1_, y1_] := Module[{l1, Leng}, Leng = Length[x1]; l1 = ConstantArray[0, Leng];
    For[i = 1, i <= Leng, i++, l1[[i]] = {{x1[[i]]}, y1[[i]]}; ]; ListInterpolation[l1, InterpolationOrder -> 1]]

LInttrapSm[f_][x1_, x2_][n_] := Module[{dx, tab}, dx = (x2 - x1)/n; tab = Table[N[f[x]], {x, x1, x2, dx}];
    (Total[tab] - (1/2)*(f[x1] + f[x2]))*dx]

LInttrapPk[f_][x1_, x2_, xpeak_][n_] := Module[{xlist1, xlist2, tab},
   xlist1 = Table[xpeak*(1 - (x1/xpeak)^(i/n)) + x1, {i, 0, n}]; xlist2 = Table[xpeak*(x2/xpeak)^(i/n), {i, 0, n}];
    Sum[(xlist1[[1 + j]] - xlist1[[j]])*(1/2)*(f[xlist1[[1 + j]]] + f[xlist1[[j]]]), {j, 1, n}] +
     Sum[(xlist2[[1 + j]] - xlist2[[j]])*(1/2)*(f[xlist2[[1 + j]]] + f[xlist2[[j]]]), {j, 1, n}]]

(* ---- flux-like interpolations from data ---- *)
dFlux[Elist_, dPhilist_] := LIntpl1D[Elist, dPhilist];

IntFlux[Elist_, dPhilist_][dE_] := Module[{Emin, Emax, dPhi, intphi, elis},
   Emin = Elist[[1]]; Emax = Elist[[-1]]; dPhi = LIntpl1D[Elist, dPhilist];
    intphi = Table[NIntegrate[dPhi[y], {y, Eval, Emax}], {Eval, Emin, Emax, dE}];
    elis = Table[Eval, {Eval, Emin, Emax, dE}]; LIntpl1D[elis, intphi]]

(* ---- sampling integrals over Ev grids ---- *)
SampleInt[f_][Ev1_, dEv_, n_] := Table[NIntegrate[f[Ev], {Ev, Ev1 + m*dEv, Ev1 + dEv + m*dEv}],
   {m, 0, n - 1, 1}]

SampleIntByEdges[f_][(edges_List)?VectorQ] := Module[{e = edges},
   If[Length[e] < 2, Return[{}]]; Table[NIntegrate[f[Ev], {Ev, e[[i]], e[[i + 1]]}], {i, 1, Length[e] - 1}]]
SampleIntByEdges::edges = "`1`";

makeEvEdgesByCount[ev1_, ev2_, (nLow_Integer)?Positive, (nHigh_Integer)?Positive, threshold_:2] :=
  Module[{lowEdges, highEdges}, Which[ev2 <= threshold, Subdivide[ev1, ev2, Max[nLow, 1]], ev1 >= threshold,
    Subdivide[ev1, ev2, Max[nHigh, 1]], True, lowEdges = Subdivide[ev1, threshold, Max[nLow, 1]];
     highEdges = Subdivide[threshold, ev2, Max[nHigh, 1]]; Join[lowEdges, Rest[highEdges]]]]

(* ---- legacy blocks from Functions.nb (not used by this pipeline) ---- *)
\[Xi] = 0;
f0[\[Epsilon]_] := 1/(Exp[\[Epsilon] - \[Xi]] + 1);

H[T_, \[Eta]_, \[Beta]_] := If[T > 5*MeV, \[Eta]*(T/(5*MeV))^\[Beta]*(Sqrt[8*Pi^3*(gstar/90)]/Mpl)*T^2,
    (Sqrt[8*Pi^3*(gstar/90)]/Mpl)*T^2];

eta = 1;
beta = 0;
SetCos[cosmo_] := Switch[cosmo, "Std", eta = 1; beta = 0; , "ST1",
    eta = 7.4*10^5; beta = -0.8; , "ST2", eta = 0.03; beta = 0; , "K", eta = 1; beta = 1; ];

\[CapitalGamma]int[\[Epsilon]_, T_] := d\[Alpha]*GF^2*\[Epsilon]*T^5;
VT[\[Epsilon]_, T_] := (-B)*\[Epsilon]*T^5;
VD[L_, T_] := 0.34*GF*T^3*L;
sinsq2\[Theta]m[L_][\[Epsilon]_, T_] := sins2\[Theta]/(sins2\[Theta] + (Sqrt[1 - sins2\[Theta]] - ((2*\[Epsilon]*T)/ms^2)*(VT[\[Epsilon], T] + VD[L, T]))^2);
lm[L_][\[Epsilon]_, T_] := (((ms^2/(2*\[Epsilon]*T))*sins2\[Theta])^2 + ((ms^2/(2*\[Epsilon]*T))*Sqrt[1 - sins2\[Theta]] - (VT[\[Epsilon], T] + VD[L, T]))^
     2)^(-2^(-1))
