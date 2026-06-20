"""Build the CGB-reconstruction figure, one per recoil-threshold, each saved in
its own threshold folder. The plot style mimics
``NeutrinoAnalysis.plot_band_comparison`` (physrev style sheet, lw=3 black
theory curves, alpha=0.22 fills with thin same-colour boundary lines, large
axis fonts, log x-axis, upper-right frameless legend).

Reads:
  <repo>/1eV/FinalFigures/panel_1eV.npz  ->  <repo>/1eV/FinalFigures/cgb_reconstruction_1eV.*
  <repo>/5eV/FinalFigures/panel_5eV.npz  ->  <repo>/5eV/FinalFigures/cgb_reconstruction_5eV.*

Per panel:
  - black solid / dashed : true input flux with / without the NC component
  - green x points       : best-fit {Phi_j} (scenario b, Exp Bkg)  [not connected]
  - purple region        : degeneracy band Delta chi^2 < 1e-3 (scenario b)
  - yellow band          : 2 sigma pointwise CL, No Bkg  (scenario a)
  - light-blue band      : 2 sigma pointwise CL, Exp Bkg (scenario b)

Exposure 3 kg yr -> energy window 0.1-3 MeV.
"""
import os
import numpy as np
import pandas as pd
from scipy import integrate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))          # <repo>/1eV/FinalFigures
REPO = os.path.dirname(os.path.dirname(HERE))              # <repo>

# Same Physical Review style sheet the module applies on import.
_STYLE = os.path.join(REPO, '1eV', 'physrev.mplstyle')
if os.path.exists(_STYLE):
    try:
        plt.style.use(_STYLE)
    except Exception:
        pass

NORM = 1e12
XLIM = (0.1, 3.0)          # 3 kg yr energy window
YLIM = (0.0, 3.0)

def _panel(thr, exposure, suffix=''):
    """One panel spec. ``suffix`` distinguishes exposures in the file names
    (empty for the default 3 kg yr)."""
    folder = os.path.join(REPO, thr, 'FinalFigures')
    npz = f"panel_{thr}{suffix}.npz"
    out = f"cgb_reconstruction_{thr}{suffix}"
    title = (rf"$E^{{\prime}}_{{\rm thr}} = {thr[:-2]}~{{\rm eV}}$"
             rf"   (${exposure:g}~{{\rm kg\,yr}}$)")
    return dict(thr=thr, npz=os.path.join(folder, npz), outdir=folder,
                outname=out, title=title)


PANELS = [
    _panel('1eV', 3),
    _panel('5eV', 3),
    _panel('1eV', 30,  '_30kgyr'),
    _panel('5eV', 300, '_300kgyr'),
    dict(thr='5eV',
         npz=os.path.join(REPO, '5eV', 'FinalFigures', 'panel_5eV_flat.npz'),
         outdir=os.path.join(REPO, '5eV', 'FinalFigures'),
         outname='cgb_reconstruction_5eV_flat',
         title=r"$E^{\prime}_{\rm thr} = 5~{\rm eV}$   ($3~{\rm kg\,yr}$, Flat Bkg)",
         single=True,    # show the Flat Bkg scenario alone (no No Bkg band)
         unfolding_dir=os.path.join(REPO, '5eV', 'FinalFigures')),
]

C_DEG   = '#7030a0'      # purple  degeneracy band
C_NOBKG = 'gold'         # yellow  No Bkg 2 sigma
C_EXPB  = 'lightblue'    # light-blue Exp Bkg 2 sigma
C_BEST  = 'green'
C_UNF   = '#d62728'      # red  regularized-unfolding comparison band


def unfolding_band(dirpath):
    """Regularized-unfolding comparison band: integrate the upper/lower fission
    spectra (scenario 2 CSVs) above each energy, converting to flux via the
    distance/normalisation constant. Returns (E, upper, lower) in cm^-2 s^-1,
    i.e. Phi(E)-Phi(2 MeV) -- directly comparable to the figure's y-axis."""
    up = pd.read_csv(os.path.join(dirpath, 'upperBand_scenario2.csv'))
    lo = pd.read_csv(os.path.join(dirpath, 'lowerBand_scenario2.csv'))
    deff = 1.0 / np.sqrt(1 / (72e2) ** 2 + 1 / (102e2) ** 2)
    const = (205.3 / 2.65e22) * (4 * np.pi * deff ** 2)
    upMeV, upF = np.array(up['MeV']), np.array(up['fissionMeV'])
    loMeV, loF = np.array(lo['MeV']), np.array(lo['fissionMeV'])
    MeVlist = np.append(upMeV, [2.0])
    U = np.interp(MeVlist, upMeV, upF).astype(float)
    L = np.interp(MeVlist, loMeV, loF).astype(float)
    for i in range(len(MeVlist)):
        U[i] = integrate.trapezoid(U[i:], MeVlist[i:]) / const
        L[i] = integrate.trapezoid(L[i:], MeVlist[i:]) / const
    return MeVlist, U, L


def _phi_ylabel(norm):
    e = int(round(np.log10(norm)))
    if e == 0:
        return r"$\Phi$ [cm$^{-2}$sec$^{-1}$]"
    return rf"$\Phi$ [$10^{{{e}}}$ cm$^{{-2}}$sec$^{{-1}}$]"


def band(ax, rows, color, z):
    """fill_between with alpha=0.22 + thin same-colour boundary lines, exactly
    like plot_band_comparison's 'fill' style."""
    if rows is None or len(rows) == 0:
        return
    E, lo, hi = rows[:, 0], rows[:, 1] / NORM, rows[:, 2] / NORM
    ok = np.isfinite(lo) & np.isfinite(hi)
    ax.fill_between(E[ok], lo[ok], hi[ok], color=color, alpha=0.22, zorder=z)
    ax.plot(E[ok], lo[ok], color=color, lw=1.0, zorder=z + 1)
    ax.plot(E[ok], hi[ok], color=color, lw=1.0, zorder=z + 1)


def make_panel(panel):
    d = np.load(panel['npz'], allow_pickle=True)
    fig, ax = plt.subplots(figsize=(8, 6))

    single = bool(panel.get('single', False))
    # widest behind, degeneracy core last
    band(ax, d['band_b'], C_EXPB, 2)         # 2nd-scenario 2 sigma
    if not single:
        band(ax, d['band_a'], C_NOBKG, 4)    # No Bkg 2 sigma
    band(ax, d['deg'],    C_DEG, 6)          # degeneracy

    # label for the 2nd-scenario band (Exp Bkg by default, e.g. Flat Bkg)
    scen2 = str(d['scen2_label']) if 'scen2_label' in d.files else 'Exp Bkg'

    # optional regularized-unfolding comparison band (upper/lower edges + fill)
    udir = panel.get('unfolding_dir')
    if udir:
        mv, U, L = unfolding_band(udir)
        ax.fill_between(mv, L / NORM, U / NORM, color=C_UNF, alpha=0.12, zorder=7)
        ax.plot(mv, U / NORM, color=C_UNF, lw=1.6, zorder=7.1)
        ax.plot(mv, L / NORM, color=C_UNF, lw=1.6, zorder=7.1)

    # best-fit as points (not connected)
    eb = d['eb']; bf = d['bestfit_phys'] / NORM
    ax.scatter(eb, bf, s=30, marker='x', color=C_BEST, zorder=8)

    # true input flux: black lw=3 solid / dashed
    ax.plot(d['x_nc'], d['Phi_nc'] / NORM, color='black', lw=3, zorder=9)
    ax.plot(d['x_no'], d['Phi_no'] / NORM, color='black', lw=3, ls='dashed', zorder=9)

    ax.set_xscale('log')
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    plt.rcParams['ytick.labelsize'] = 20
    plt.rcParams['xtick.labelsize'] = 20
    ax.set_xlabel(r"$E_\nu$ [MeV]", fontsize=30)
    ax.set_ylabel(_phi_ylabel(NORM), fontsize=30)
    ax.xaxis.set_minor_locator(
        LogLocator(base=10.0, subs=np.arange(1.0, 10) * 0.1, numticks=20))
    ax.tick_params(axis='both', which='major', labelsize=23)
    ax.tick_params(axis='both', which='minor', labelsize=23)

    handles = [
        Line2D([0], [0], color='black', lw=3, label='With NC'),
        Line2D([0], [0], color='black', lw=3, ls='dashed', label='Without NC'),
        Line2D([0], [0], color=C_BEST, marker='x', lw=0, label=r'Best-fit $\{\Phi_j\}$'),
        Patch(facecolor=C_DEG, alpha=0.22, edgecolor=C_DEG, label=r'$\Delta\chi^2<10^{-3}$'),
    ]
    if not single:
        handles.append(Patch(facecolor=C_NOBKG, alpha=0.22, edgecolor=C_NOBKG,
                             label=r'$2\sigma$ No Bkg'))
    handles.append(Patch(facecolor=C_EXPB, alpha=0.22, edgecolor=C_EXPB,
                         label=rf'$2\sigma$ {scen2}'))
    if panel.get('unfolding_dir'):
        handles.append(Line2D([0], [0], color=C_UNF, lw=1.6,
                              label='Reg. unfolding'))
    ax.set_title(panel['title'], fontsize=20)
    ax.legend(handles=handles, loc='upper right', fontsize=15, frameon=False)

    for ext in ('pdf', 'png'):
        out = os.path.join(panel['outdir'], f"{panel['outname']}.{ext}")
        fig.savefig(out, bbox_inches='tight')
        print('saved', out)
    plt.close(fig)


def main():
    for panel in PANELS:
        if not os.path.exists(panel['npz']):
            print('skip (missing)', panel['npz'])
            continue
        make_panel(panel)


if __name__ == '__main__':
    main()
