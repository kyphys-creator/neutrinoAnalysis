(* ========================================================================== *)
(*  03_functions_response  (5 eV threshold)
*)
(*  Building blocks for the response: nuclear form factor, CEvNS differential cross sections, threshold efficiency, and the reactor-flux interpolations built from the imported spectra. From Functions.nb and 5_Calculation.nb.
*)
(* ========================================================================== *)

(* ---- nuclear form factor & cross sections ---- *)
F[q_] := (((4*Pi*rho0)/(A*q^3))*(Sin[q*RA] - q*RA*Cos[q*RA]))/(1 + (a*q)^2);

diffCross[Mn_][ER_, Ev_] := ((GF^2*Mn)/(4*Pi))*qw^2*(1 - (Mn*ER)/(2*Ev^2))*F[Sqrt[2*Mn*ER]]^2;
DvdiffCross[Mn_][ER_, Ev_] := ((GF^2*Mn)/(4*Pi))*qw^2*((Mn*ER)/Ev^3)*F[Sqrt[2*Mn*ER]]^2;

(* ---- detector threshold efficiency (legacy keeps 1 eV here even in the
   5 eV set; the recoil bins start at 5 eV so the cut is inert) ---- *)
epsilon[Ep_] := If[Ep >= 1*eV, 1, 0]

(* ---- reactor-antineutrino flux from the imported spectra ---- *)
datEvsolid = Transpose[datdNdEsolid][[1]]*MeV;
datdPhisolid = (Transpose[datdNdEsolid][[2]]/MeV)*(((2.65*10^22*(MeV/sec))/(205.3*MeV))*(1/(4*Pi))*
     (1/(7200*cm)^2 + 1/(10200*cm)^2));

dfluxsolid = dFlux[datEvsolid, datdPhisolid];
intfluxsolid = IntFlux[datEvsolid, datdPhisolid][0.001*MeV];

datEvdash = Transpose[datdNdEdash][[1]]*MeV;
datdPhidash = (Transpose[datdNdEdash][[2]]/MeV)*(((2.65*10^22*(MeV/sec))/(205.3*MeV))*(1/(4*Pi))*
     (1/(7200*cm)^2 + 1/(10200*cm)^2));

dfluxdash = dFlux[datEvdash, datdPhidash];
intfluxdash = IntFlux[datEvdash, datdPhidash][0.001*MeV];
