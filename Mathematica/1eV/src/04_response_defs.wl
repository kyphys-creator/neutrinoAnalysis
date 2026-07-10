(* ========================================================================== *)
(*  04_response_defs  (1 eV threshold)
*)
(*  Response-function definitions: expected rate Exprate (+resolution variants), differential responses CurlR/CurlRbox and the CRloop builders. The per-bin response functions themselves are built in 06_response_functions. From Functions.nb, ExpRate.nb.
*)
(* ========================================================================== *)

(* ---- expected rate for one recoil energy ---- *)
Exprate[ER_, emax_, dPhi_] := Module[{emin}, emin = (ER + Sqrt[ER*(ER + 2*Mn)])/2;
    NT*NIntegrate[diffCross[Mn][ER, Ev]*dPhi[Ev], {Ev, emin, emax}]]

ExprateResolution[Ep_, emax_, dPhi_, sigmaE_] := Module[{ERmin, ERmax},
   (ERmin = Ep - sigmaE/2; ERmax = Ep + sigmaE/2; (ERmin = Max[0, ERmin]);
     (epsilon[Ep]*NIntegrate[(1/sigmaE)*UnitBox[(ER - Ep)/sigmaE]*Exprate[ER, emax, dPhi], {ER, ERmin, ERmax}]))]

ExprateResolutionBox[Ep_, emax_, dPhi_, sigmaE_] :=
  epsilon[Ep]*NIntegrate[Exprate[ER, emax, dPhi]*Gbox[ER, Ep, sigmaE], {ER, eV/10^4, 200*eV}]

ExprateNoG[ER_, emax_, intf_] := Module[{emin}, emin = (ER + Sqrt[ER*(ER + 2*Mn)])/2;
    NIntegrate[NT*DvdiffCross[Mn][ER, Ev]*intf[Ev], {Ev, emin, emax}]]

(* ---- differential response dcurlyR/dE' per recoil bin ---- *)
CurlRbox[Ev_, e1_, e2_, sig_] := Module[{ERmax}, ERmax = (2*Ev^2)/(Mn + 2*Ev);
    NIntegrate[NT*DvdiffCross[Mn][ER, Ev]*IntGboxmu[ER, e1, e2, sig], {ER, eV/10^4, ERmax}]]

CurlR[Ev_, e1_, e2_, sig_] := Module[{ERmax}, ERmax = (2*Ev^2)/(Mn + 2*Ev);
    If[Max[1*eV, e1] < Min[e2, ERmax], NIntegrate[NT*DvdiffCross[Mn][ER, Ev], {ER, Max[1*eV, e1], Min[e2, ERmax]}], 0]]

CRloopBox[ev1_, ev2_, dev_, sigER_] := Module[{tabe, tabR, funclist, evList, datERpair},
   evList = Range[ev1, ev2, dev]; datERpair = Partition[datERbin, 2, 1];
    funclist = Table[tabR = ParallelTable[CurlRbox[ev, datERpair[[j,1]]*eV, datERpair[[j,2]]*eV, sigER], {ev, evList}];
       tabe = evList; LIntpl1D[tabe, tabR], {j, Length[datERpair]}]; funclist]

CRloop[ev1_, ev2_, dev_, sigER_] := Module[{tabe, tabR, funclist, evList, datERpair},
   evList = Range[ev1, ev2, dev]; (datERpair = Partition[datERbin, 2, 1]);
    (funclist = Table[tabR = ParallelTable[CurlR[ev, datERpair[[j,1]]*eV, datERpair[[j,2]]*eV, sigER],
          {ev, evList}]; tabe = evList; LIntpl1D[tabe, tabR], {j, Length[datERpair]}]); funclist]

CRloopBoxData[ev1_, ev2_, dev_, sigER_] := Module[{datlist, evList, datERpairs},
    evList = Range[ev1, ev2, dev]; datERpairs = Transpose[{Most[datERbin], Rest[datERbin]}]*eV;
     datlist = ParallelTable[Table[CurlRbox[ev, datERpair[[1]], datERpair[[2]], sigER], {ev, evList}],
       {datERpair, datERpairs}]; datlist]

CRloopData[ev1_, ev2_, dev_, sigER_] := Module[{datlist, evList, datERpairs},
    evList = Range[ev1, ev2, dev]; datERpairs = Transpose[{Most[datERbin], Rest[datERbin]}]*eV;
     datlist = ParallelTable[Table[CurlR[ev, datERpair[[1]], datERpair[[2]], sigER], {ev, evList}],
       {datERpair, datERpairs}]; datlist]

(* ---- response with Ev kept as a parameter -------------------------------
   The Ev dependence of CurlRbox/CurlR enters only through the prefactor
   1/Ev^3 and the kinematic upper limit ERmax(Ev) = 2 Ev^2/(Mn + 2 Ev), so
   the recoil integral can be done ONCE per bin:

       curlyR_j(Ev) = G_j(ERmax(Ev)) / Ev^3 ,
       G_j(X) = Integrate[kern(ER) Box_j(ER), {ER, X0, X}] .

   G_j is obtained by a single NDSolve in the log variable ER = X0 Exp[s]
   (the log variable resolves the low-ER bins; validated against kink-aware
   high-precision NIntegrate to <~1e-4 relative). A fully symbolic Integrate
   exists but its Ei/Ci closed form suffers catastrophic branch-cut and
   cancellation problems in the qR << 1 regime relevant here, so the
   one-shot numeric antiderivative is used instead. ---- *)

CRkernel[ER_] := NT*((GF^2*Mn)/(4*Pi))*qw^2*Mn*ER*F[Sqrt[2*Mn*ER]]^2;
ERmaxOf[Ev_] := (2*Ev^2)/(Mn + 2*Ev);

CurlRboxG[e1_, e2_, sig_] := Module[{g, X0 = eV/10^4, X1, smax},
  X1 = ERmaxOf[7*MeV]*1.05; smax = Log[X1/X0];
  g = NDSolveValue[{gg'[s] == X0*Exp[s]*CRkernel[X0*Exp[s]]*IntGboxmu[X0*Exp[s], e1, e2, sig],
        gg[0] == 0}, gg, {s, 0, smax},
        AccuracyGoal -> 20, PrecisionGoal -> 12, MaxSteps -> 10^6,
        MaxStepSize -> smax/4000];
  Function[Ev, g[Log[Min[ERmaxOf[Ev], X1]/X0]]/Ev^3]];

CurlRG[e1_, e2_] := Module[{g, a, X1, smax},
  a = Max[1*eV, e1]; X1 = ERmaxOf[7*MeV]*1.05;
  If[a >= Min[e2, X1], Return[Function[Ev, 0]]];
  smax = Log[X1/a];
  g = NDSolveValue[{gg'[s] == a*Exp[s]*CRkernel[a*Exp[s]], gg[0] == 0},
        gg, {s, 0, smax},
        AccuracyGoal -> 20, PrecisionGoal -> 12, MaxSteps -> 10^6,
        MaxStepSize -> smax/4000];
  Function[Ev, With[{X = Min[e2, ERmaxOf[Ev], X1]},
    If[X <= a, 0, g[Log[X/a]]/Ev^3]]]];
