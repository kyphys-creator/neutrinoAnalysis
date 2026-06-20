"""Extract one panel's plotting arrays for the final CGB-reconstruction figure.

Run from inside the threshold folder so the data + module resolve, e.g.

    cd 1eV && PYTHONPATH=. python3 FinalFigures/_extract_panel.py 1eV \
        /abs/path/1eV/FinalFigures/panel_1eV.npz
    cd 5eV && PYTHONPATH=. python3 /abs/path/1eV/FinalFigures/_extract_panel.py 5eV \
        /abs/path/1eV/FinalFigures/panel_5eV.npz

Scenario b (Exp Bkg) anchors the best-fit step and the Delta-chi^2 < 1e-3
degeneracy band; the 2sigma bands come from the saved standard band JSONs
for scenario a (No Bkg) and b (Exp Bkg). Everything is stored in physical
units (cm^-2 s^-1), which are GeV-invariant.
"""
import sys
import glob
import json
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')

from neutrino_analysis_band import NeutrinoAnalysis, load_band

THR_LABEL = sys.argv[1]
OUT       = sys.argv[2]
T_EXP     = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0   # exposure [kg yr]
SCEN2     = sys.argv[4] if len(sys.argv) > 4 else 'b'          # 2nd background scenario
TLABEL    = f'T{int(T_EXP)}' if float(T_EXP).is_integer() else f'T{T_EXP}'
LEVEL     = 0.954                      # 2 sigma
SCEN_LABELS = {'a': 'No Bkg', 'b': 'Exp Bkg', 'flat': 'Flat Bkg',
               'b2': '10x Exp Bkg', 'b100': '100x Exp Bkg'}
SCEN2_LABEL = SCEN_LABELS.get(SCEN2, SCEN2)
DEG_THR   = 1e-3                       # degeneracy band: Delta chi^2 / c < 1e-3

# Reconstructed-flux energy grid per recoil threshold. A higher E'_thr raises
# the minimum probeable neutrino energy: 1 eV -> 0.18-2.0 MeV, 5 eV -> 0.41-2.0.
ERANGE = {'1eV': (0.18, 2.0), '5eV': (0.41, 2.0)}
DEG_IDX   = [0, 1, 2, 3, 4, 5, 8, 13, 16, 20, 25, 30, 40, 50, 60, 70, 80, 90,
             100, 110, 120, 130, 140, 150, 160, 167]


def degeneracy_interval(a, idx, v0, chi2min, res_x):
    """[lower, upper] in RAW units where the profiled Delta chi^2/c < DEG_THR.
    The profiled chi^2 is convex in the fixed value, so a single crossing per
    side: bracket outward until above threshold, then bisect."""
    data = a.data_vector

    def dchi(v):
        r = a.optimize_with_fixed_parameter(data, idx, float(v), x0=res_x.copy())
        return r.fun / a.c - chi2min

    def edge(direction, factor=1.4, nbracket=40, nbisect=40):
        v_in = v0                                  # inside (dchi < thr)
        v_out = None
        v = v0
        for _ in range(nbracket):
            v = v * factor if direction > 0 else v / factor
            if v <= 0:
                return 0.0
            if dchi(v) > DEG_THR:
                v_out = v
                break
            v_in = v
        if v_out is None:
            return v_in                            # no crossing within range
        lo, hi = (v_in, v_out) if direction > 0 else (v_out, v_in)
        for _ in range(nbisect):
            m = 0.5 * (lo + hi)
            inside = dchi(m) < DEG_THR
            if direction > 0:
                lo, hi = (m, hi) if inside else (lo, m)
            else:
                lo, hi = (lo, m) if inside else (m, hi)
        return 0.5 * (lo + hi)

    return edge(-1), edge(+1)


def band_rows(pattern, level, eb):
    """Sorted (energy, lo, hi) physical-unit rows for one level."""
    rows = []
    for f in sorted(glob.glob(pattern)):
        b = load_band(f)
        for lv, (lo, hi) in b['band_physical'].items():
            if abs(lv - level) < 1e-6:
                rows.append((eb[b['index']], lo, hi))
    rows.sort()
    return np.array(rows) if rows else np.empty((0, 3))


def main():
    a = NeutrinoAnalysis(background_scenario=SCEN2, intervals='180',
                         GeV=0.32e16, solver='osqp', T=T_EXP)
    res = a.optimize(a.data_vector)
    unit = a.cm ** 2 * a.sec
    emin, emax = ERANGE.get(THR_LABEL, (0.18, 2.0))
    eb = np.linspace(emin, emax, a.n)
    chi2min = res.fun / a.c

    # true input flux: with NC (solid) and without NC (dashed)
    x_nc, Phi_nc, x_no, Phi_no = a._calculate_integrated_flux()

    # best-fit step (scenario b)
    bestfit_phys = res.x * unit

    # degeneracy band Delta chi^2/c < 1e-3 (scenario b)
    deg = []
    for idx in DEG_IDX:
        lo, hi = degeneracy_interval(a, idx, float(res.x[idx]), chi2min, res.x)
        deg.append((eb[idx], lo * unit, hi * unit))
        print(f"[{THR_LABEL}] deg idx={idx:3d} E={eb[idx]:.3f} "
              f"[{lo*unit:.3e}, {hi*unit:.3e}]", flush=True)
    deg = np.array(deg)

    # 2 sigma bands from saved standard JSONs (strict idxNNN glob skips _bkgpen).
    # band_a = No Bkg baseline (scenario a); band_b = the 2nd scenario (SCEN2).
    band_a = band_rows(f'{TLABEL}/scenario_bkg_a/bands/band_bkga_idx[0-9][0-9][0-9].json', LEVEL, eb)
    band_b = band_rows(f'{TLABEL}/scenario_bkg_{SCEN2}/bands/band_bkg{SCEN2}_idx[0-9][0-9][0-9].json', LEVEL, eb)

    # sanity: result.x*unit should match a saved best_fit_physical
    chk = sorted(glob.glob(f'{TLABEL}/scenario_bkg_{SCEN2}/bands/band_bkg{SCEN2}_idx040.json'))
    if chk:
        bf = load_band(chk[0])['best_fit_physical']
        print(f"[{THR_LABEL}] unit check: bestfit[40]={bestfit_phys[40]:.4e} "
              f"vs json={bf:.4e} (ratio {bestfit_phys[40]/bf:.4f})", flush=True)

    np.savez(OUT, thr=THR_LABEL, eb=eb, unit=unit,
             x_nc=x_nc, Phi_nc=Phi_nc, x_no=x_no, Phi_no=Phi_no,
             bestfit_phys=bestfit_phys, deg=deg, band_a=band_a, band_b=band_b,
             level=LEVEL, deg_thr=DEG_THR, scen2_label=SCEN2_LABEL)
    print(f"[{THR_LABEL}] WROTE {OUT}  (deg={len(deg)}, "
          f"band_a={len(band_a)}, band_b={len(band_b)})", flush=True)


if __name__ == '__main__':
    main()
