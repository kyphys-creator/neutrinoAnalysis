"""Worker: generate one penalty-mode confidence band (temporary helper)."""
import sys
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
from neutrino_analysis_band import NeutrinoAnalysis

scen, idx = sys.argv[1], int(sys.argv[2])
a = NeutrinoAnalysis(background_scenario=scen, intervals='180',
                     GeV=0.32e16, solver='osqp', T=3, bkg_penalty=True)
a.optimize(a.data_vector)
a.find_and_save_band(
    idx,
    levels=(0.678, 0.954),
    num_pseudo_data=50,
    n_pseudo_edge=500,
    step=1.5,
    rel_tol=0.03,
    seed=42,
    n_jobs=1,
    verbose=False,
)
print(f"done {scen} idx{idx:03d}", flush=True)
