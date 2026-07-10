(* ========================================================================== *)
(*  03b_flux_plot  (1 eV threshold)
*)
(*  Plots of the differential and integrated reactor fluxes built in 03_functions_response. From 5_Calculation.nb / 1_Main.nb (Flux data).
*)
(* ========================================================================== *)

dfluxPlot = Plot[{dfluxsolid[x*MeV]*cm^2*sec*MeV, dfluxdash[x*MeV]*cm^2*sec*MeV}, {x, 0.18, 7},
  ScalingFunctions -> {"Log", None}, PlotStyle -> {{Black}, {Black, Dashed}}, Frame -> True,
  FrameLabel -> {"\!\(\*SubscriptBox[\(E\), \(v\)]\)[MeV]",
    "d\[CurlyPhi]/\!\(\*SubscriptBox[\(dE\), \(v\)]\)[1/\!\(\*SuperscriptBox[\(cm\), \(2\)]\)/s/MeV]"}];
Export[FileNameJoin[{outputDir, "flux_differential.pdf"}], dfluxPlot];

intfluxPlot = Plot[{intfluxsolid[x*MeV]*cm^2*sec, intfluxdash[x*MeV]*cm^2*sec}, {x, 0.18, 7},
  ScalingFunctions -> {"Log", None}, PlotStyle -> {{Black}, {Black, Dashed}}, Frame -> True,
  FrameLabel -> {"\!\(\*SubscriptBox[\(E\), \(v\)]\)[MeV]", "\[CapitalPhi][1/\!\(\*SuperscriptBox[\(cm\), \(2\)]\)/s]"}];
Export[FileNameJoin[{outputDir, "flux_integrated.pdf"}], intfluxPlot];
