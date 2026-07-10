"""Band-width vs parameter-count scan, with the background-penalty option on.

For the 1 eV case, sweep the number of flux parameters n (= CRmat columns,
selected via ``intervals``) with ``bkg_penalty=True`` and measure the
confidence-band width at a few fixed physical energies. The data bin count
m = 29 is fixed, so lowering n reduces the underdetermination.

Run from the 1eV/ directory so the data files resolve:

    cd 1eV && PYTHONPATH=. python3 param_scan_bkgpen/run_param_scan.py

Writes param_scan_bkgpen/param_scan_results.json (one record per
scenario x n x target-energy x confidence-level).
"""
import os
import json
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')

from neutrino_analysis_band import NeutrinoAnalysis

OUTDIR    = 'param_scan_bkgpen'
SCENARIOS = ['a', 'b']
NPARS     = [180]
TARGET_E  = [0.40, 0.60, 0.90]          # MeV; band measured at nearest index
LEVELS    = (0.678, 0.954)
MC = dict(num_pseudo_data=40, n_pseudo_edge=300, step=1.5, rel_tol=0.04,
          seed=42, n_jobs=1, verbose=False)

os.makedirs(OUTDIR, exist_ok=True)
records = []
for scen in SCENARIOS:
    for n in NPARS:
        a = NeutrinoAnalysis(background_scenario=scen, intervals=str(n),
                             GeV=0.32e16, solver='osqp', T=3, bkg_penalty=True)
        a.optimize(a.data_vector)
        eb = np.linspace(0.18, 2, a.n)
        for E in TARGET_E:
            idx = int(np.argmin(np.abs(eb - E)))
            band = a.find_confidence_band(idx, levels=LEVELS, **MC)
            for lv in LEVELS:
                lo, hi = band['band_physical'][lv]
                records.append(dict(
                    scenario=scen, n=int(n), target_E=float(E),
                    idx=idx, E=float(eb[idx]), level=float(lv),
                    lower=float(lo), upper=float(hi), width=float(hi - lo),
                    best_fit=float(band['best_fit_physical']),
                    bkg_penalty=True))
            print(f"scen={scen} n={n:3d} E~{E:.2f} idx={idx:3d} "
                  f"(E={eb[idx]:.3f}) done", flush=True)

out = os.path.join(OUTDIR, 'param_scan_results_nNum.json')
with open(out, 'w') as f:
    json.dump(records, f, indent=1)
print(f"WROTE {out}  records={len(records)}", flush=True)
