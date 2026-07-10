(* ========================================================================== *)
(*  07_response_function_plot  (1 eV threshold)
*)
(*  Load the response functions saved by 06_response_functions (.wdx) and plot them over the continuous neutrino energy. From 1_Main.nb.
*)
(* ========================================================================== *)

CRlistbox = Import[FileNameJoin[{outputDir, "CRlistbox.wdx"}]];

CRplotlist = Table[CRlistbox[[i]][x*MeV], {i, 1, Length[CRlistbox]}];

responsePlot = Plot[CRplotlist, {x, 0.18, 7}, PlotRange -> Full];
Export[FileNameJoin[{outputDir, "response_functions.pdf"}], responsePlot];

responsePlotScaled = Plot[CRplotlist*MeV*10^3*(grams/cm^2), {x, 0.18, 7},
  PlotRange -> {{0.18, 3}, {0, 5/10^17}}, MaxRecursion -> 60];
Export[FileNameJoin[{outputDir, "response_functions_originalUnit.pdf"}], responsePlotScaled];
