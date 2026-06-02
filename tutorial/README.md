# Tutorial

A self-contained, hands-on tour of the `neutrino_analysis_band` API. The
module and the data it needs (`CRmat/`, `Ratebin/`, `Danny's files/`) are
bundled here, so run the notebooks **with this `tutorial/` folder as the
working directory** — no setup beyond the dependencies.

## Dependencies

```bash
pip install numpy scipy matplotlib pandas numba
pip install cvxpy osqp        # for solver='osqp'
pip install highspy           # optional; fastest staircase solver (SciPy is the bundled fallback)
```

## Lessons (run in order)

| # | Notebook | Topic |
|---|---|---|
| 00 | `00_setup.ipynb` | Setup, orientation, the constructor arguments |
| 01 | `01_optimize_and_plot.ipynb` | `optimize`, reading the result, `plot_flux_comparison` |
| 02 | `02_solvers_and_staircase.ipynb` | `scipy` vs `osqp`, vertex selection (staircase flux) |
| 03 | `03_chi2_scan.ipynb` | `scan_fixed_parameter` (manual Δχ² grid) |
| 04 | `04_confidence_band.ipynb` | `find_confidence_band` (root finding) |
| 05 | `05_save_and_overlay.ipynb` | `find_and_save_band`, `plot_flux_with_bands` |
| 06 | `06_scenario_comparison.ipynb` | `plot_band_comparison` across scenarios |
| 07 | `07_varying_T.ipynb` | the `T` argument and `T<T>/` output foldering |
| 08 | `08_pitfalls_and_troubleshooting.ipynb` | GeV precision, `inf` edges, no-staircase, noise |

The Monte-Carlo sizes in the lessons are kept small so each notebook runs
quickly. For real results use the larger values noted in lesson 08 (and in the
top-level `README.md`).

`_build_tutorial.py` is the generator script for these notebooks; it is not
needed to follow the tutorial.
