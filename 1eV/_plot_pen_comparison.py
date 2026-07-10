"""Penalty vs no-penalty band comparison figure (temporary helper).

One axis per confidence level: scenarios a (no detector bkg, B=h only) and
b (exponential detector bkg), each with and without the background penalty.
Scenario = colour, penalty = line style, so the penalty effect is readable.
"""
import glob
import json
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

from neutrino_analysis_band import NeutrinoAnalysis, load_band

NORM = 1e12
LEVELS = {0.954: '2sigma', 0.678: '1sigma'}
SCEN = {'a': ('C0', r'bkg a ($b_i=0$)'), 'b': ('C3', 'bkg b (exp.)')}

a0 = NeutrinoAnalysis(background_scenario='a', intervals='180',
                      GeV=0.32e16, solver='osqp', T=3)
eb = np.linspace(0.18, 2, a0.n)
x, Phi, x2, Phid = a0._calculate_integrated_flux()


def rows_for(scen, pen, level):
    suf = '_bkgpen' if pen else ''
    files = sorted(glob.glob(
        f'T3/scenario_bkg_{scen}/bands/band_bkg{scen}_idx*{suf}.json'))
    if not pen:
        files = [f for f in files if '_bkgpen' not in f]
    out = []
    for f in files:
        b = load_band(f)
        for lv in b['band_physical']:
            if abs(lv - level) < 1e-6:
                lo, hi = b['band_physical'][lv]
                out.append((eb[b['index']], lo, hi))
    out.sort()
    return (np.array([r[0] for r in out]),
            np.array([r[1] for r in out]) / NORM,
            np.array([r[2] for r in out]) / NORM)


for level, tag in LEVELS.items():
    plt.figure(figsize=(8, 6))
    plt.plot(x, Phi / NORM, color='black', lw=2, label='With NC')
    plt.plot(x2, Phid / NORM, color='black', lw=2, ls='dotted', label='Without NC')

    for scen, (color, slabel) in SCEN.items():
        for pen in (False, True):
            xs, lo, hi = rows_for(scen, pen, level)
            if len(xs) == 0:
                print(f'[warn] no bands for {scen} pen={pen}')
                continue
            ok = np.isfinite(lo) & np.isfinite(hi)
            if pen:
                plt.plot(xs[ok], lo[ok], color=color, lw=2.2, ls='--')
                plt.plot(xs[ok], hi[ok], color=color, lw=2.2, ls='--',
                         label=f'{slabel}, with penalty')
            else:
                plt.fill_between(xs[ok], lo[ok], hi[ok], color=color,
                                 alpha=0.18, zorder=1)
                plt.plot(xs[ok], lo[ok], color=color, lw=1.6)
                plt.plot(xs[ok], hi[ok], color=color, lw=1.6,
                         label=f'{slabel}, no penalty')

    plt.xscale('log')
    plt.xlim(0.15, 2.2)
    plt.ylim(0, 3.0)
    plt.xlabel(r'$E_\nu$ [MeV]', fontsize=22)
    plt.ylabel(r'$\Phi(>E_\nu)$ [$10^{12}$ cm$^{-2}$ s$^{-1}$]', fontsize=22)
    nsig = '2' if level > 0.9 else '1'
    plt.title(rf'{nsig}$\sigma$ CL bands, T=3', fontsize=18)
    plt.gca().xaxis.set_minor_locator(
        LogLocator(base=10.0, subs=np.arange(1.0, 10) * 0.1, numticks=20))
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.legend(loc='upper right', fontsize=12, frameon=False)
    plt.tight_layout()
    fn = f'band_comparison_T3_{tag}_pen_vs_nopen.pdf'
    plt.savefig(fn)
    print(f'saved {fn}')

# quantitative summary: width ratios and a-b separation per index
print('\nlevel scen  median-width-ratio(pen/nopen)')
for level in LEVELS:
    for scen in SCEN:
        xs0, lo0, hi0 = rows_for(scen, False, level)
        xs1, lo1, hi1 = rows_for(scen, True, level)
        common = sorted(set(np.round(xs0, 6)) & set(np.round(xs1, 6)))
        r = []
        for cx in common:
            i0 = np.argmin(abs(xs0 - cx)); i1 = np.argmin(abs(xs1 - cx))
            w0 = hi0[i0] - lo0[i0]; w1 = hi1[i1] - lo1[i1]
            if np.isfinite(w0) and np.isfinite(w1) and w0 > 0:
                r.append(w1 / w0)
        print(f'{level:.3f}  {scen}   x{np.median(r):.3f}  (n={len(r)})')
