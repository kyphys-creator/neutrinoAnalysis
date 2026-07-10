(* ========================================================================== *)
(*  08_response_matrix  (5 eV threshold)
*)
(*  Choose the number of neutrino-energy intervals, integrate the response functions of 06_response_functions (.wdx) over each interval and export the response matrices as csv to ../output. From 1_Main.nb; all interval grids start at the 5 eV neutrino-energy minimum 0.41 MeV.
*)
(* ========================================================================== *)

CRlistbox = Import[FileNameJoin[{outputDir, "CRlistbox.wdx"}]];
CRlist = Import[FileNameJoin[{outputDir, "CRlist.wdx"}]];

(* number of intervals *)
m = Length[ERbin] - 1;
n = 360;

(* ---- n intervals up to 2 MeV ---- *)
CRmat = {};
t4R = 0;
For[i = 1, i <= Length[CRlistbox], i++,
  edges = makeEvEdgesByCount[0.41*MeV, 2*MeV, n, 20, 2*MeV];
  t4R = SampleIntByEdges[CRlistbox[[i]]][edges];
  AppendTo[CRmat, t4R]];
Export[FileNameJoin[{outputDir, "CRmat360.csv"}], CRmat];
CRmatOrg = CRmat*10^3*(grams/cm^2);
Export[FileNameJoin[{outputDir, "CRmat360_originalUnit.csv"}], CRmatOrg];

(* ---- same interval count extended to 7 MeV ---- *)
CRmat = {};
t4R = 0;
For[i = 1, i <= Length[CRlistbox], i++,
  edges = makeEvEdgesByCount[0.41*MeV, 7*MeV, n, 20, 2*MeV];
  t4R = SampleIntByEdges[CRlistbox[[i]]][edges];
  AppendTo[CRmat, t4R]];
Export[FileNameJoin[{outputDir, "CRmat640.csv"}], CRmat];
CRmatOrg = CRmat*10^3*(grams/cm^2);
Export[FileNameJoin[{outputDir, "CRmat640_originalUnit.csv"}], CRmatOrg];

(* ---- uniform sampling variant (CRlist, no box smearing) ---- *)
CRmat = {};
t4R = 0;
For[i = 1, i <= Length[CRlist], i++,
  t4R = SampleInt[CRlist[[i]]][0.41*MeV, (6.59/n)*MeV, n];
  AppendTo[CRmat, t4R]];
Export[FileNameJoin[{outputDir, "CRmat80.csv"}], CRmat];
