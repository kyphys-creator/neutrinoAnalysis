(* ========================================================================== *)
(*  05_curlyR_plot  (5 eV threshold)
*)
(*  Plot the bin-integrated response functions curlyR_j of 04_response_defs
    over the continuous neutrino energy, via CurlRboxG (recoil integral done
    once per bin, Ev kept as a parameter). y axis in cm^2/(MeV kg).
*)
(* ========================================================================== *)

datERpairs = Partition[datERbin, 2, 1];

curlyRfns = Table[CurlRboxG[p[[1]]*eV, p[[2]]*eV, 1*eV], {p, datERpairs}];

(* natural units -> cm^2/(MeV kg) *)
toPhys = MeV*(10^3*grams)/cm^2;

curlyRPlot = Plot[
  Evaluate[Table[curlyRfns[[j]][x*MeV]*toPhys, {j, Length[datERpairs]}]],
  {x, 0.41, 7}, PlotRange -> Full, PlotPoints -> 500, MaxRecursion -> 5,
  Frame -> True,
  FrameLabel -> {"\!\(\*SubscriptBox[\(E\), \(\[Nu]\)]\) [MeV]",
    "\[ScriptCapitalR] [\!\(\*SuperscriptBox[\(cm\), \(2\)]\)/(MeV kg)]"}];
Export[FileNameJoin[{outputDir, "curlyR.pdf"}], curlyRPlot];
