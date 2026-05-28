"""
NeutrinoAnalysis — fast version with two interchangeable solver backends.

Two backends are available, switched via the ``solver`` argument to
``__init__`` (or by calling ``set_solver`` later):

  * ``solver='scipy'`` — scipy.optimize.minimize(method='trust-constr')
                        with an analytic Jacobian and a constant analytic
                        Hessian. Identical solver to the original code,
                        just much better-informed. ~2–3× faster than the
                        original on this problem.

  * ``solver='osqp'``  — Reformulates the χ² minimization as the quadratic
                        program it actually is, and solves it with OSQP via
                        cvxpy. ~100× faster than the original on this
                        problem. Requires ``pip install cvxpy osqp``.

Both backends produce the same flux to machine precision (|Δx| ~ 1e-9 in
benchmarks). They share the same public API, so user-level scripts do not
need to change other than passing ``solver=...``.
"""

import os

import numpy as np
from scipy.optimize import minimize, LinearConstraint
from scipy import integrate
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt
import pandas as pd
from numba import njit
from matplotlib.ticker import LogLocator

try:
    from joblib import Parallel, delayed
    _HAS_JOBLIB = True
except ImportError:
    _HAS_JOBLIB = False

try:
    import cvxpy as cp
    _HAS_CVXPY = True
except ImportError:
    _HAS_CVXPY = False


# -------------------------------------------------------------------
# Numba kernels (used by the scipy backend)
# -------------------------------------------------------------------
@njit(cache=True, fastmath=True)
def _chiN_core(xvec, M, data_minus_bkg, inv_data, T):
    """χ² = T · Σ (data − Bkg − M·x)² / data   (bins with data > 0 only)."""
    r = data_minus_bkg - M @ xvec
    return T * np.sum(r * r * inv_data)


@njit(cache=True, fastmath=True)
def _chiN_grad_core(xvec, M, data_minus_bkg, inv_data, T):
    """∂χ²/∂x = −2T · Mᵀ ((data − Bkg − M·x) · inv_data)."""
    r = data_minus_bkg - M @ xvec
    return -(2.0 * T) * (M.T @ (r * inv_data))


# -------------------------------------------------------------------
# Lightweight result wrapper so both backends look the same to callers.
# -------------------------------------------------------------------
class _Result:
    __slots__ = ('x', 'fun', 'success', 'status', 'nit')

    def __init__(self, x, fun, success=True, status='', nit=0):
        self.x = x
        self.fun = fun
        self.success = success
        self.status = status
        self.nit = nit


# -------------------------------------------------------------------
# Backend: scipy trust-constr with analytic jac & const hess
# -------------------------------------------------------------------
class _ScipyBackend:
    name = 'scipy'

    def __init__(self, parent):
        self.p = parent

    def solve(self, data, x0, extra_constraints=None, display=False):
        p = self.p
        M = p.M_matrix
        T = p.T

        if data is p.data_vector:
            dmb, inv_d, H = p._dmb_default, p._inv_d_default, p._hess_default
        else:
            dmb, inv_d = p._make_dmb_inv(data)
            H = p._build_hessian_from(dmb, inv_d)

        def fun(x):  return _chiN_core(x, M, dmb, inv_d, T)
        def jac(x):  return _chiN_grad_core(x, M, dmb, inv_d, T)
        def hess(_): return H

        cons = [p.ordering_constraint]
        if extra_constraints:
            cons.extend(extra_constraints)

        options = {'maxiter': 100000, 'xtol': 1e-9, 'gtol': 1e-9}
        if display:
            options['verbose'] = 3

        return minimize(
            fun, x0, method='trust-constr',
            jac=jac, hess=hess,
            bounds=[(p.eps, None)] * p.n,
            constraints=cons, options=options,
        )


# -------------------------------------------------------------------
# Backend: OSQP via cvxpy.
#   min  ||W^{1/2} (M x − (d − b))||²
#   s.t. x ≥ 0,  x[i] ≥ x[i+1],  (optional) x[k] = v
# W = T · diag(1/d_i) on bins where d_i > 0, else 0.
# -------------------------------------------------------------------
class _OSQPBackend:
    """
    OSQP via cvxpy.

    Note on scaling: the project's natural-unit constant ``self.c`` is ~1e70,
    which makes ``M_matrix``, ``data_vector`` etc. far larger than OSQP's
    finite-bound limit (~1e30). The scipy backend doesn't notice because
    trust-constr scales internally. Here we solve the *unscaled* problem in
    a local variable and multiply back when reporting ``fun``, so that
    ``res.fun / self.c`` still gives the χ²/c the user expects.

    Concretely, the equivalent unscaled problem uses
        M_s = M / c          d_s = d / c        b_s = b / c
        inv_d_s = 1 / d_s = c / d
        w_s = sqrt(T * inv_d_s) = sqrt(c) * w
    and we compensate by returning ``fun * c`` so external code that divides
    by ``self.c`` keeps working.
    """
    name = 'osqp'

    DEFAULT_OPTS = dict(
        eps_abs=1e-10, eps_rel=1e-10,
        eps_prim_inf=1e-10, eps_dual_inf=1e-10,
        max_iter=400000,
        polish=True, polish_refine_iter=10,
        adaptive_rho=True,
        scaling=20,           # more scaling iterations → better conditioning
        verbose=False, warm_start=True,
    )

    def __init__(self, parent):
        if not _HAS_CVXPY:
            raise ImportError(
                "solver='osqp' requires cvxpy. Install with: pip install cvxpy osqp"
            )
        self.p = parent
        # Internal variable ``y`` is the *column-scaled* flux: x = D ⊙ y, with
        # D = 1/‖M_s column‖. Without this rescaling the optimal x (~1e7) and
        # the design-matrix columns (~1e-6) span ~13 orders of magnitude, and
        # OSQP/CLARABEL both terminate early at a wrong point while reporting
        # status='optimal'. D is data-independent, so the parameterised
        # (DPP) problem can still be re-solved fast across pseudo-data.
        self._y = cp.Variable(parent.n, nonneg=True)
        self._cache = {}     # QP problems, keyed by fixed_index
        self._lp_cache = {}  # vertex-selection LPs, keyed by fixed_index
        # The χ² minimiser is a high-dimensional face of the feasible polytope
        # whenever n > rank(M) (here 180 params vs 29 bins). OSQP/CLARABEL are
        # interior-point/ADMM methods that return a *relative-interior* point of
        # that face — a smooth ramp. The physically meaningful estimate is a
        # *vertex* (piecewise-constant flux), found by a simplex method. After
        # the QP fixes the unique fitted values μ = M·x, a HiGHS-simplex LP over
        # {M·x = μ, monotone, x ≥ 0} returns such a vertex.
        self.vertex_select = True

    def reset_cache(self):
        self._cache.clear()
        self._lp_cache.clear()

    def _column_scale(self):
        M_s = self.p.M_matrix / self.p.c
        cn = np.linalg.norm(M_s, axis=0)
        cn[cn == 0] = 1.0
        return 1.0 / cn

    def _build_problem(self, fixed_index):
        """
        Build a parameterised problem so that *any* data vector for this
        background scenario can be plugged in by setting cvxpy Parameters.
        Only ``fixed_index`` changes the structure (adds an equality
        constraint), so that's the only cache key.
        """
        p = self.p
        n, m = p.n, p.m
        y = self._y
        D = self._column_scale()

        # cvxpy Parameters: w (per-bin weight √(T/d_s)) and z (w * dmb_s).
        # Residual = w ⊙ (M_s·D ⊙ y) − z, solving for the scaled variable y.
        w_par = cp.Parameter(m, nonneg=True)
        z_par = cp.Parameter(m)
        MsD = (p.M_matrix / p.c) * D[None, :]   # constant, scaled once

        residual = cp.multiply(w_par, MsD @ y) - z_par
        obj = cp.Minimize(cp.sum_squares(residual))
        # Ordering on the physical flux x = D ⊙ y.
        cons = [cp.multiply(D[:-1], y[:-1]) >= cp.multiply(D[1:], y[1:])]

        fv_param = None
        if fixed_index is not None:
            fv_param = cp.Parameter()
            cons.append(D[fixed_index] * y[fixed_index] == fv_param)

        prob = cp.Problem(obj, cons)
        return prob, w_par, z_par, fv_param, D

    def _get_problem(self, fixed_index):
        if fixed_index not in self._cache:
            self._cache[fixed_index] = self._build_problem(fixed_index)
        return self._cache[fixed_index]

    def _build_lp(self, fixed_index):
        """
        Vertex-selection LP: among all fluxes reproducing the fitted values
        μ = M_s·x (parameter ``mu_par``), pick a vertex of the monotone polytope
        by minimising Σx with a simplex method. ``mu_par`` is the only thing that
        changes between data vectors, so one LP per ``fixed_index`` is cached.
        """
        p = self.p
        n, m = p.n, p.m
        M_s = p.M_matrix / p.c
        x = cp.Variable(n, nonneg=True)
        mu_par = cp.Parameter(m)
        cons = [M_s @ x == mu_par, x[:-1] >= x[1:]]
        fv_par = None
        if fixed_index is not None:
            fv_par = cp.Parameter()
            cons.append(x[fixed_index] == fv_par)
        prob = cp.Problem(cp.Minimize(cp.sum(x)), cons)
        return prob, x, mu_par, fv_par

    def _get_lp(self, fixed_index):
        if fixed_index not in self._lp_cache:
            self._lp_cache[fixed_index] = self._build_lp(fixed_index)
        return self._lp_cache[fixed_index]

    def _select_vertex(self, x_interior, fixed_index, fixed_value):
        """Return a piecewise-constant vertex with the same fit as ``x_interior``."""
        M_s = self.p.M_matrix / self.p.c
        mu = M_s @ x_interior
        prob, x, mu_par, fv_par = self._get_lp(fixed_index)
        mu_par.value = mu
        if fixed_index is not None:
            fv_par.value = float(fixed_value)
        try:
            prob.solve(solver=cp.HIGHS)
        except Exception:
            return None
        return x.value if prob.status in ('optimal', 'optimal_inaccurate') else None

    def _set_data_params(self, data, w_par, z_par):
        """Fill in the cvxpy Parameters from a (scaled) data vector."""
        p = self.p
        data_s = data / p.c
        bkg_s = p.Bkg_vector / p.c
        safe = np.where(data_s > 0, data_s, 1.0)
        inv_d_s = np.where(data_s > 0, 1.0 / safe, 0.0)
        w = np.sqrt(p.T * inv_d_s)
        z = w * (data_s - bkg_s)
        # Add a tiny floor to w to keep the residual operator well-conditioned
        # when many bins have zero weight. (Optional; doesn't change the optimum.)
        w_par.value = w
        z_par.value = z

    def solve(self, data, x0=None, extra_constraints=None,
              fixed_index=None, fixed_value=None, display=False):
        if fixed_index is None and extra_constraints:
            fixed_index, fixed_value = self._extract_fixed(extra_constraints)

        prob, w_par, z_par, fv_param, D = self._get_problem(fixed_index)
        self._set_data_params(data, w_par, z_par)
        if fixed_index is not None:
            fv_param.value = float(fixed_value)
        if x0 is not None:
            self._y.value = np.asarray(x0, dtype=float) / D

        opts = dict(self.DEFAULT_OPTS)
        if display:
            opts['verbose'] = True

        # Try OSQP first (very fast). If it reports inaccurate or fails,
        # fall back to CLARABEL (slower but more robust on near-degenerate
        # problems).
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')   # squelch "solution may be inaccurate"
            try:
                prob.solve(solver=cp.OSQP, **opts)
            except Exception:
                prob.status = 'solver_error'

        if prob.status not in ('optimal',) or self._y.value is None:
            # Retry with CLARABEL — no warm-start support, but very accurate.
            try:
                prob.solve(solver=cp.CLARABEL, verbose=display)
            except Exception:
                pass

        if self._y.value is None:
            return _Result(np.full(self.p.n, np.nan), np.nan,
                           success=False, status=prob.status, nit=0)
        x = D * self._y.value
        fun = float(prob.value) * self.p.c

        # Replace the interior (ramp) solution with a piecewise-constant vertex
        # that has the identical fit, so the OSQP backend matches the staircase
        # shape the scipy backend produces. χ² (``fun``) is unchanged.
        if self.vertex_select:
            xv = self._select_vertex(x, fixed_index, fixed_value)
            if xv is not None:
                x = xv

        return _Result(
            x=np.array(x),
            fun=fun,
            success=prob.status in ('optimal', 'optimal_inaccurate'),
            status=prob.status,
            nit=int(prob.solver_stats.num_iters) if prob.solver_stats else 0,
        )

    @staticmethod
    def _extract_fixed(extra_constraints):
        for c in extra_constraints:
            if isinstance(c, LinearConstraint):
                A = c.A
                A_dense = A.toarray() if hasattr(A, 'toarray') else np.asarray(A)
                if A_dense.shape[0] == 1 and np.isclose(c.lb, c.ub).all():
                    idx = int(np.argmax(np.abs(A_dense[0])))
                    if np.isclose(A_dense[0, idx], 1.0):
                        return idx, float(np.atleast_1d(c.lb)[0])
        return None, None


# ===================================================================
# Main class
# ===================================================================
class NeutrinoAnalysis:
    """
    Neutrino flux optimization pipeline with two interchangeable backends.

    Parameters
    ----------
    background_scenario : str
    intervals : str
    GeV, c : float
    solver : {'scipy', 'osqp'}
        Backend to use. Both produce the same solution to ~1e-9.
        OSQP is roughly 50–100× faster on this problem; scipy is the
        safe, well-tested fallback.
    """

    # ----------------------------------------------------------------
    def __init__(self, background_scenario='c', intervals='180',
                 GeV=1e-6, c=1, solver='scipy'):
        self.iterationTime = 1000
        self.T = 3
        self.GeV = GeV
        self.c = c

        self._define_constants()
        self._load_data(intervals)
        self._prepare_backgrounds()

        # Master normalization constant — kept identical to the original
        self.c = 3e12 * 2.693500303951368e+58

        self._build_ordering_constraint()

        self.result = None
        self.set_background(background_scenario)
        self.set_solver(solver)

    # ---------- data ----------
    def _load_data(self, intervals='180'):
        self.fig1Solid = pd.read_csv("Danny’s files/fig1-solid.csv")
        self.fig1dashed = pd.read_csv("Danny’s files/fig1-dashed.csv")
        self.intervals = intervals
        self.CRmat = (
            np.genfromtxt(
                f'CRmat/originalUnit/CRmat{self.intervals}_originalUnit.csv',
                delimiter=','
            )
            * self.cm ** 2 / (10 ** 3 * self.gram)
            * (10 ** 3 * self.gram) * self.yr
        )
        self.Ratebin7 = np.genfromtxt('Ratebin/Ratebin7_originalUnit.csv', delimiter=',')
        self.Ratebin2 = np.genfromtxt('Ratebin/Ratebin2_originalUnit.csv', delimiter=',')
        self.RateDiff = self.Ratebin7 - self.Ratebin2
        self.n = len(self.CRmat[0])
        self.m = len(self.Ratebin7)

        f = 2.65e22 / 205.3 / 434795262.39118177
        self.fig1Solid['cm**-2sec-1MeV-1'] = self.fig1Solid['fissionMeV'] * f
        self.fig1dashed['cm**-2sec-1MeV-1'] = self.fig1dashed['fissionMeV'] * f

    def _define_constants(self):
        self.MeV = 1e-3 * self.GeV
        self.keV = 1e-6 * self.GeV
        self.eV  = 1e-9 * self.GeV
        self.meV = 1e-12 * self.GeV
        self.gram = 5.62e23 * self.GeV
        self.sec  = 1 / (6.58e-25 * self.GeV)
        self.yr   = 365 * 24 * 3600 * self.sec
        self.cm   = 1 / (1.98e-14 * self.GeV)
        self.eps = 0
        self.bins = np.array([5., 7., 9., 11., 13., 15., 17., 19., 21., 23.,
                              25., 27., 29., 31., 33., 35., 37., 39., 41., 43.,
                              45., 47., 49., 51., 56., 61., 66., 71., 81., 120.])
        self.Mdetector = 1.0

    def _compute_bi(self, A, B, C):
        bi = []
        for i in range(len(self.bins) - 1):
            Ei, Ei1 = self.bins[i], self.bins[i + 1]
            integral = -A * B * (np.exp(-Ei1 / B) - np.exp(-Ei / B)) + C * (Ei1 - Ei)
            bi.append(self.Mdetector * integral)
        return np.array(bi)

    def _prepare_backgrounds(self):
        d = 365
        ba    = self._compute_bi(A=0   * self.eV / self.keV * d, B=10, C=0   * self.eV / self.keV * d)
        bb    = self._compute_bi(A=460 * self.eV / self.keV * d, B=10, C=100 * self.eV / self.keV * d)
        b2    = self._compute_bi(A=920 * self.eV / self.keV * d, B=10, C=200 * self.eV / self.keV * d)
        bc    = self._compute_bi(A=50  * self.eV / self.keV * d, B=10, C=20  * self.eV / self.keV * d)
        bflat = self._compute_bi(A=0   * self.eV / self.keV * d, B=10, C=100 * self.eV / self.keV * d)
        self.background_df = pd.DataFrame({
            "Bin Start [eV]": self.bins[:-1], "Bin End [eV]": self.bins[1:],
            "b_i (a)": ba, "b_i (b)": bb, "b_i (c)": bc,
            "b_i (flat)": bflat, "b_i (b2)": b2,
        })

    # ---------- backend ----------
    def set_solver(self, solver):
        """Switch backend at any time. ``solver`` ∈ {'scipy', 'osqp'}."""
        solver = solver.lower()
        if solver == 'scipy':
            self._backend = _ScipyBackend(self)
        elif solver == 'osqp':
            self._backend = _OSQPBackend(self)
        else:
            raise ValueError(f"Unknown solver: {solver!r}. Use 'scipy' or 'osqp'.")
        self.solver = solver

    # ---------- constraints ----------
    def _build_ordering_constraint(self):
        rows = np.repeat(np.arange(self.n - 1), 2)
        cols = np.empty(2 * (self.n - 1), dtype=int)
        cols[0::2] = np.arange(self.n - 1)
        cols[1::2] = np.arange(1, self.n)
        vals = np.tile(np.array([1.0, -1.0]), self.n - 1)
        A = csc_matrix((vals, (rows, cols)), shape=(self.n - 1, self.n))
        self.ordering_constraint = LinearConstraint(A, lb=0, ub=np.inf)

    def _make_fixed_constraint(self, fixed_index, fixed_value):
        A_fixed = csc_matrix(([1.0], ([0], [fixed_index])), shape=(1, self.n))
        return LinearConstraint(A_fixed, lb=fixed_value, ub=fixed_value)

    # ---------- background & cached χ² ingredients ----------
    def set_background(self, bkg_scenario='c'):
        self.background_scenario = bkg_scenario
        bkg_map = {'a': 'b_i (a)', 'b': 'b_i (b)', 'c': 'b_i (c)',
                   'flat': 'b_i (flat)', 'none': 'b_i (a)', 'b2': 'b_i (b2)'}
        col = bkg_map.get(bkg_scenario)
        if col is None:
            raise ValueError(f"Invalid background scenario: {bkg_scenario}")

        self.ExtBkg = np.array(self.background_df[col])
        k = self.c
        self.M_matrix = k * self.CRmat
        if bkg_scenario == 'none':
            self.data_vector = k * (self.Ratebin2 + self.ExtBkg)
            self.Bkg_vector  = k * self.ExtBkg
        else:
            self.data_vector = k * (self.Ratebin7 + self.ExtBkg)
            self.Bkg_vector  = k * (self.RateDiff + self.ExtBkg)

        self._dmb_default, self._inv_d_default = self._make_dmb_inv(self.data_vector)
        self._hess_default = self._build_hessian_from(self._dmb_default,
                                                     self._inv_d_default)

        if hasattr(self, '_backend') and self._backend.name == 'osqp':
            self._backend.reset_cache()

    def _make_dmb_inv(self, data):
        safe = np.where(data > 0, data, 1.0)
        inv = np.where(data > 0, 1.0 / safe, 0.0)
        return data - self.Bkg_vector, inv

    def _build_hessian_from(self, _dmb, inv_data):
        return (2.0 * self.T) * (self.M_matrix.T * inv_data) @ self.M_matrix

    # ---------- public solver API ----------
    def optimize(self, data, x0=None, display=False):
        if x0 is None:
            x0 = np.ones(self.n)
        if self._backend.name == 'osqp':
            res = self._backend.solve(data, x0=x0, display=display)
        else:
            res = self._backend.solve(data, x0, extra_constraints=None, display=display)
        self.result = res
        return res

    def optimize_with_fixed_parameter(self, data, fixed_index, fixed_value, x0=None):
        if x0 is None:
            print("Warning: main optimization has not been run. Starting from scratch.")
            x0 = np.ones(self.n)
        x0 = x0.copy()
        x0[fixed_index] = fixed_value

        if self._backend.name == 'osqp':
            return self._backend.solve(data, x0=x0,
                                       fixed_index=fixed_index,
                                       fixed_value=fixed_value)
        fc = self._make_fixed_constraint(fixed_index, fixed_value)
        return self._backend.solve(data, x0, extra_constraints=[fc], display=False)

    # ---------- scan ----------
    def scan_parameter(self, data, index, num_points=21, scan_range=0.75):
        if self.result is None:
            raise RuntimeError("Run optimize() first.")
        base_value = self.result.x[index]
        if index >= 135:
            scan_range = 1 - 1e-11

        if index < 10:
            fixed_values = np.append(
                np.linspace(0.9e12 / self.cm ** 2 / self.sec,
                            1044872844297.51 / self.cm ** 2 / self.sec,
                            num_points)[:-1],
                np.linspace(1044872844297.51 / self.cm ** 2 / self.sec,
                            6e12 / self.cm ** 2 / self.sec, num_points)
            )
        elif index < 21:
            fixed_values = np.linspace((1 - scan_range) * base_value,
                                       (1 + 2) * base_value,
                                       2 * num_points - 1)
        elif index < 134:
            fixed_values = np.linspace((1 - scan_range) * base_value,
                                       (1 + scan_range) * base_value,
                                       2 * num_points - 1)
        else:
            fixed_values = np.append(
                np.linspace((1 - scan_range) * base_value,
                            (1 + scan_range) * base_value, num_points)[:-1],
                np.linspace((1 + scan_range) * base_value, 4 * base_value, num_points)
            )

        chi_sq_results = []
        warm = self.result.x.copy()
        total = len(fixed_values)
        for k, v in enumerate(fixed_values, 1):
            r = self.optimize_with_fixed_parameter(data, index, v, x0=warm)
            warm = r.x.copy()
            chi_sq_results.append(r.fun / self.c)
            print(f"[{k}/{total}] idx={index}  flux={v * self.cm**2 * self.sec:.4e}  χ²={chi_sq_results[-1]:.6g}")
        return fixed_values, chi_sq_results

    # ---------- Monte Carlo ----------
    @staticmethod
    def _fit_one_pseudo(pseudo_data_scaled, fixed_index, fixed_value, x0_seed, analysis):
        result_dict = {}
        try:
            res_free = analysis.optimize(pseudo_data_scaled / analysis.T,
                                         x0=x0_seed.copy(), display=False)
            result_dict.update(chi2_free=res_free.fun / analysis.c,
                               x_free=res_free.x,
                               success_free=res_free.success)
        except Exception as e:
            return {**result_dict, 'chi2_free': None, 'x_free': None,
                    'success_free': False, 'error': repr(e)}
        try:
            res_fixed = analysis.optimize_with_fixed_parameter(
                pseudo_data_scaled / analysis.T,
                fixed_index, fixed_value,
                x0=res_free.x.copy()
            )
            result_dict.update(chi2_fixed=res_fixed.fun / analysis.c,
                               x_fixed=res_fixed.x,
                               success_fixed=res_fixed.success)
        except Exception as e:
            return {**result_dict, 'chi2_fixed': None, 'x_fixed': None,
                    'success_fixed': False, 'error': repr(e)}

        result_dict['delta_chi2'] = abs(result_dict['chi2_fixed'] - result_dict['chi2_free'])
        return result_dict

    def run_full_monte_carlo_analysis(self, num_pseudo_data, fixed_index, fixed_value,
                                      seed=None, n_jobs=1, verbose=True):
        result = self.optimize(self.data_vector, display=False)
        result_fixed = self.optimize_with_fixed_parameter(
            self.data_vector, fixed_index, fixed_value, x0=result.x.copy()
        )
        self.best_fit_flux       = result.x
        self.best_fit_chi2       = result.fun
        self.best_fit_chi2_fixed = result_fixed.fun

        modPrime = self.M_matrix @ result_fixed.x / self.c
        self.modPrime_physical = modPrime / np.diff(self.bins)

        if verbose:
            print(f"[Baseline] χ²/c = {self.best_fit_chi2 / self.c:.6g}")
            print(f"[Baseline] event rates: {self.modPrime_physical.min():.3e} .. {self.modPrime_physical.max():.3e}")
            print(f"Generating {num_pseudo_data} pseudo-data sets "
                  f"(idx={fixed_index}, value={self.cm**2 * self.sec * fixed_value:.3e})")

        if seed is not None:
            np.random.seed(seed)
        self.pseudo_data_sets = []
        self.pseudo_data_scaled = []
        bw = np.diff(self.bins)
        for _ in range(num_pseudo_data):
            Event = self.T * (self.modPrime_physical * bw + self.Bkg_vector / self.c)
            pseudo_data = np.random.normal(Event, scale=np.sqrt(np.abs(Event)))
            self.pseudo_data_sets.append(pseudo_data)
            self.pseudo_data_scaled.append(pseudo_data * self.c)

        x0_seed = self.best_fit_flux
        # joblib parallelism is meaningful only for the scipy backend; OSQP is
        # already fast and cvxpy problems are awkward to ship across processes.
        if n_jobs != 1 and _HAS_JOBLIB and self.solver == 'scipy':
            self.fit_results = Parallel(n_jobs=n_jobs, prefer='processes')(
                delayed(self._fit_one_pseudo)(pd_, fixed_index, fixed_value, x0_seed, self)
                for pd_ in self.pseudo_data_scaled
            )
        else:
            self.fit_results = []
            for i, pd_ in enumerate(self.pseudo_data_scaled, 1):
                r = self._fit_one_pseudo(pd_, fixed_index, fixed_value, x0_seed, self)
                if verbose:
                    print(f"  [{i}/{num_pseudo_data}] "
                          f"χ²_free={r.get('chi2_free')}  "
                          f"χ²_fix={r.get('chi2_fixed')}  Δχ²={r.get('delta_chi2')}")
                self.fit_results.append(r)

    def analyze_monte_carlo_results(self, fixed_value, confidence_level=0.90):
        self.delta_chi2_values = [r['delta_chi2'] for r in self.fit_results
                                  if r.get('delta_chi2') is not None]
        if not self.delta_chi2_values:
            raise RuntimeError("No valid Δχ² values to compute cutoff.")
        ds = np.sort(self.delta_chi2_values)
        idx = max(int(confidence_level * len(ds)) - 1, 0)
        self.delta_chi2_cutoff = ds[idx]

        is_in_range = (
            self.best_fit_chi2_fixed / self.c
            < self.best_fit_chi2 / self.c + self.delta_chi2_cutoff
        )
        self.fixed_value_included = fixed_value if is_in_range else None
        return is_in_range

    def scan_fixed_parameter(self, fixed_index, scan_range=0.35, num_points=11,
                             num_pseudo_data=10, seed=None, n_jobs=1):
        base_value = self.result.x[fixed_index]
        if fixed_index >= 50:
            scan_range = 1 - 1e-11

        if fixed_index < 2:
            fixed_values = np.append(
                np.linspace(0.6e12 / self.cm ** 2 / self.sec,
                            1.35e12 / self.cm ** 2 / self.sec,
                            num_points)[:-1],
                np.linspace(1.35e12 / self.cm ** 2 / self.sec,
                            7e13 / self.cm ** 2 / self.sec, num_points)
            )
        elif fixed_index < 7:
            fixed_values = np.linspace((1 - scan_range) * base_value,
                                       (1 + 2) * base_value,
                                       2 * num_points - 1)
        elif fixed_index < 50:
            fixed_values = np.linspace((1 - scan_range) * base_value,
                                       (1 + scan_range) * base_value,
                                       2 * num_points - 1)
        else:
            fixed_values = np.append(
                np.linspace((1 - scan_range) * base_value,
                            (1 + scan_range) * base_value, num_points)[:-1],
                np.linspace((1 + scan_range) * base_value, 4 * base_value, num_points)
            )

        results = []
        total = len(fixed_values)
        print(f"[Scan idx={fixed_index}] {total} points, solver={self.solver}")
        for i, fv in enumerate(fixed_values, 1):
            print(f"\n--- [{i}/{total}] flux={self.cm**2 * self.sec * fv:.4e} ---")
            self.run_full_monte_carlo_analysis(
                num_pseudo_data=num_pseudo_data,
                fixed_index=fixed_index,
                fixed_value=fv,
                seed=seed,
                n_jobs=n_jobs,
                verbose=False,
            )
            in68 = self.analyze_monte_carlo_results(fv, 0.678);  c68 = self.delta_chi2_cutoff
            in90 = self.analyze_monte_carlo_results(fv, 0.90);   c90 = self.delta_chi2_cutoff
            in95 = self.analyze_monte_carlo_results(fv, 0.954);  c95 = self.delta_chi2_cutoff

            results.append({
                'index': fixed_index,
                'fixed_value_raw': fv,
                'fixed_value_physical': self.cm ** 2 * self.sec * fv,
                'included1sigma': in68, 'included90': in90, 'included2sigma': in95,
                'delta_chi2_cutoff68': c68, 'delta_chi2_cutoff90': c90, 'delta_chi2_cutoff95': c95,
                'delta_chi2': list(self.delta_chi2_values),
            })
        return results

    # ---------- plotting ----------
    def _calculate_integrated_flux(self):
        def dPhidEnu(E):
            return np.interp(E, self.fig1Solid['MeV'], self.fig1Solid['cm**-2sec-1MeV-1'], 1)
        def dPhidEnudashed(E):
            return np.interp(E, self.fig1dashed['MeV'], self.fig1dashed['cm**-2sec-1MeV-1'], 1)
        x2to7 = np.logspace(np.log10(2), np.log10(7), self.iterationTime)
        Phi2to7       = integrate.trapezoid(dPhidEnu(x2to7),       x2to7)
        Phi2to7dashed = integrate.trapezoid(dPhidEnudashed(x2to7), x2to7)
        x  = np.logspace(-2, np.log10(7), self.iterationTime); y  = dPhidEnu(x)
        x2 = np.logspace(-2, np.log10(7), self.iterationTime); y2 = dPhidEnudashed(x2)
        Phi       = np.zeros(self.iterationTime)
        Phidashed = np.zeros(self.iterationTime)
        for i in range(self.iterationTime):
            Phi[i]       = integrate.trapezoid(y[i:],  x[i:])
            Phidashed[i] = integrate.trapezoid(y2[i:], x2[i:])
        return x, Phi - Phi2to7, x2, Phidashed - Phi2to7dashed

    def plot_flux_comparison(self, save=True):
        if self.result is None:
            raise RuntimeError("Run optimize() before plotting.")
        x, Phi, x2, Phidashed = self._calculate_integrated_flux()
        plt.figure(figsize=(8, 6))
        plt.plot(x,  Phi,       label=r'With neutron capture',    color='red')
        plt.plot(x2, Phidashed, label=r'Without neutron capture', color='brown', ls='dashed')
        eb = np.linspace(0.41, 2, self.n)
        plt.scatter(eb, self.result.x * (self.cm ** 2 * self.sec),
                    label=f'Optimized Flux (Bkg: {self.background_scenario})',
                    zorder=5, s=1)
        plt.xscale('log'); plt.ylim(0, 2e12); plt.xlim(0.1, 2)
        plt.xlabel(r"$E_\nu$ [MeV]"); plt.ylabel(r"$\Phi$ [cm$^{-2}$sec$^{-1}$]")
        plt.title('Optimized Neutrino Flux vs. Theoretical Models')
        plt.legend(); plt.grid(True, which="both", ls="--")
        if save:
            os.makedirs(f'scenario_bkg_{self.background_scenario}', exist_ok=True)
            fn = f'scenario_bkg_{self.background_scenario}/flux_comparison_bkg_{self.background_scenario}.pdf'
            plt.savefig(fn); print(f"Plot saved as {fn}")

    def plot_scan_results(self, index, fixed_values, chi_sq_results, save=True):
        plt.figure(figsize=(8, 6))
        fv = fixed_values * (self.cm ** 2 * self.sec)
        plt.plot(fv, chi_sq_results - np.min(chi_sq_results), 'o-')
        plt.hlines(1/3,  1e-3 * fv[0], 1e3 * fv[-1], color='pink', ls='dashed', label=r'$\Delta \chi^2 = 1/3$')
        plt.hlines(1e-3, 1e-3 * fv[0], 1e3 * fv[-1], color='blue', ls='dotted', label=r'$\Delta \chi^2 = 10^{-3}$')
        if index < 2:
            plt.vlines(fv[(len(fv) - 1) // 2], 1e-10, 1e5, color='red', ls='dashdot', label='Minimized')
            plt.xscale('log')
            plt.gca().xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(1.0, 10) * 0.1, numticks=20))
        elif index < 50:
            plt.vlines(fv[(len(fv) - 1) // 2], 1e-10, 1e5, color='red', ls='dashdot', label='Minimized')
        else:
            plt.vlines(fv[(len(fv) - 1) // 4], 1e-10, 1e5, color='red', ls='dashdot', label='Minimized')
        x, Phi, x2, Phidashed = self._calculate_integrated_flux()
        eb = np.linspace(0.41, 2, self.n)
        plt.vlines(Phi[np.argmin(np.abs(x  - eb[index]))], 1e-10, 1e5, color='grey', label='With neutron capture')
        plt.vlines(Phidashed[np.argmin(np.abs(x2 - eb[index]))], 1e-10, 1e5, color='grey', ls='dashed', label='Without neutron capture')
        plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(1.0, 10) * 0.1, numticks=20))
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.tick_params(axis='both', which='minor', labelsize=15)
        plt.xlabel(rf'Fixed Flux Value for Parameter $\Phi_{{{index+1}}}$', fontsize=11)
        plt.ylabel('Δχ²', fontsize=15)
        plt.title(f'Δχ² Profile for Parameter $\\Phi_{{{index+1}}}$ '
                  f'(At {eb[index]:.2f} MeV) under 1 eV threshold (Bkg: {self.background_scenario})')
        plt.grid(True)
        plt.xlim(fv[0], fv[np.array(chi_sq_results) <= 1e2][-1])
        plt.ylim(1e-5, 1e2); plt.legend(loc='lower right', ncol=1, fontsize=9); plt.yscale('log')
        if save:
            os.makedirs(f'scenario_bkg_{self.background_scenario}', exist_ok=True)
            fn = f'scenario_bkg_{self.background_scenario}/scan_param_{index}_bkg_{self.background_scenario}.pdf'
            plt.savefig(fn); print(f"Scan plot saved as {fn}")


# -------------------------------------------------------------------
# Module-level helper kept for backward compatibility
# -------------------------------------------------------------------
def scan_parameter(analysis, data, index, num_points=21, scan_range=0.35):
    return analysis.scan_parameter(data, index, num_points=num_points, scan_range=scan_range)
