(* ========================================================================== *)
(*  05b_curlyR_export  (5 eV threshold)
*)
(*  Evaluate the bin-integrated response functions curlyR_j (CurlRboxG)
    on 2000 neutrino-energy points and export as csv to ../output:
    column 1 = Ev [MeV], one column per recoil bin, values in
    cm^2/(MeV kg).
*)
(* ========================================================================== *)

datERpairs = Partition[datERbin, 2, 1];

curlyRfns = Table[CurlRboxG[p[[1]]*eV, p[[2]]*eV, 1*eV], {p, datERpairs}];

toPhys = MeV*(10^3*grams)/cm^2;

npts = 2000;
evgrid = N[Subdivide[0.41, 7.0, npts - 1]];

curlyRheader = Prepend[
   Table[ToString[Round[datERpairs[[j,1]]]] <> "-" <>
         ToString[Round[datERpairs[[j,2]]]] <> " eV",
     {j, Length[datERpairs]}], "Ev_MeV"];

curlyRtable = Table[
   Prepend[Table[curlyRfns[[j]][ev*MeV]*toPhys, {j, Length[curlyRfns]}], ev],
   {ev, evgrid}];

Export[FileNameJoin[{outputDir, "curlyR_table.csv"}],
  Prepend[curlyRtable, curlyRheader]];
