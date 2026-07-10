(* ========================================================================== *)
(*  01_setup  (1 eV threshold)
*)
(*  Inputs: directories, energy scale, units, target/flux constants, data bins and measured dN/dE spectra (../input). From Units.nb, Constants_NuFlux.nb, 2_3_Input_Data.nb and 1_Main.nb (Setting).
*)
(* ========================================================================== *)

(* ---- directories ---- *)
Direc = DirectoryName[$InputFileName];
inputDir = FileNameJoin[{ParentDirectory[Direc], "input"}];
outputDir = FileNameJoin[{ParentDirectory[Direc], "output"}];
Quiet[CreateDirectory[outputDir]];

(* ---- energy scale ---- *)
GeV = 10^15;

(* ---- units ---- *)
meV = GeV/10^12;
eV = GeV/10^9;
keV = GeV/10^6;
MeV = GeV/10^3;

grams = 5.62*10^23*GeV;
Mpl = 1.22*10^19*GeV;
Kel = (8.62*GeV)/10^14;

cm = ((1.98*GeV)/10^14)^(-1);
fm = cm/10^13;
sec = ((6.58*GeV)/10^25)^(-1);
yr = 365*24*3600*sec;

GF = 1.166/(10^5*GeV^2);
h = 0.7;
\[Rho]DM = 0.26*(1.05/10^5)*h^2*GeV*cm^(-3);
n\[Nu]\[Alpha] = 112*cm^(-3);

(* ---- target & flux constants (Ge) ---- *)
A = 72;
Z = 32;
Nn = A - Z;
a = 0.7*fm;
RA = 1.2*A^(1/3)*fm;
sin2W = 0.23873;
qw = Nn - (1 - 4*sin2W)*Z;
rho0 = (3*A)/(4*3.14*RA^3);
Mn = (0.93149410372*A - 0.0725)*GeV;
NT = (1/(71.78724990521931*grams))*6.022*10^23;
(* NOTE: legacy Constants_NuFlux.nb also carried scratch re-assignments
   (e.g. RA = 1.2 A^(1/3) without fm) below the canonical block; they are
   intentionally NOT reproduced here. *)

(* ---- experimental recoil bins [eV] and neutrino spectra ---- *)
datERbin = {1., 3., 5., 7., 9., 11., 13., 15., 17., 19., 21., 23., 25., 27., 29., 31., 33., 35., 37., 39.,
    41., 43., 45., 47., 49., 51., 56., 61., 66., 71., 81., 120.};
datERbin2 = {1., 3., 5., 7., 9., 11., 13., 15., 17., 19., 21., 23., 25., 27., 29., 31., 33., 35., 37., 39.,
    41., 43., 45., 47., 49., 51., 56., 61., 66., 71., 81., 116.};
ERbin = datERbin*1*eV;
ERbin2 = datERbin2*1*eV;

datdNdEsolid = Import[FileNameJoin[{inputDir, "dNdEsolid.csv"}], "Data"];
datdNdEdash = Import[FileNameJoin[{inputDir, "dNdEdash.csv"}], "Data"];

(* sanity: recoil energies for E_nu = 0.3 / 0.5 MeV (from 1_Main "Setting") *)
(2*(0.3*MeV)^2)/(Mn + 2*(0.3*MeV))/eV
(2*(0.5*MeV)^2)/(Mn + 2*(0.5*MeV))/eV
