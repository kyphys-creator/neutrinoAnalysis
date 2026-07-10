(* ========================================================================== *)
(*  06_response_functions  (1 eV threshold)
*)
(*  Integrate the differential response of 04_response_defs over each experimental recoil bin, giving one response function of the neutrino energy per bin (CRlistbox: box-smeared, CRlist: sharp). The interpolating-function data is stored as .wdx in ../output.
*)
(* ========================================================================== *)

CRlistbox = CRloopBox[0.18*MeV, 7*MeV, 0.01*MeV, 1*eV];
CRlist = CRloop[0.18*MeV, 7*MeV, 0.01*MeV, 1*eV];

Export[FileNameJoin[{outputDir, "CRlistbox.wdx"}], CRlistbox];
Export[FileNameJoin[{outputDir, "CRlist.wdx"}], CRlist];
