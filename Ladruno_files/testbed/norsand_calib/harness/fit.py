"""Constrained fit of (chi, h, N, N_bar, rho, rho_bar) (task item 4).

Refusal rules as constraints (owner-approved, sheet S.39 / §4.2; O2 params.validate):
  N_bar <= N,   rho/rho_bar >= beta = (1-N)/(1-N_bar),   WW: rho, rho_bar in (1/2, 1].
They are built into a sequential reparametrisation, so scipy least_squares only sees the box u in [0, 1]^k and
every u decodes to an admissible set. Decoding order chi, h, N, N_bar, rho_bar, rho; each free parameter is
lo + u (hi - lo) on an interval that depends on the ones already decoded (and on the fixed ones):
  chi      [CHI_LO, CHI_HI] (linear)               h   [H_LO, H_HI] (log)
  N        [max(N_LO, N_bar if fixed, 1 - (1-N_bar) rho/rho_bar if N_bar, rho, rho_bar all fixed), N_HI]
  N_bar    [0, min(N, 1 - (1-N) rho_bar/rho if rho and rho_bar both fixed)]
  rho_bar  [RHO_LO, min(1, rho/beta if rho fixed)]
  rho      [max(RHO_LO, beta rho_bar), 1]
(RHO_LO = 1/2 + 1e-3: rho = 1/2 itself is refused.) An empty interval is an infeasible fixed set: ValueError.

Optimiser: least_squares(method='trf', bounds [0, 1]) from a scrambled Sobol multi-start in u (+ optional given
starts), in parallel processes. Identifiability at the optimum: the residual Jacobian in log-scaled physical
coordinates (J_j = d r / d ln|theta_j|, central differences), its singular values and condition number, the
correlation matrix of (J^T J)^-1, and profile scans (cost with one parameter fixed on a grid, the others re-fit).
"""
from __future__ import annotations

import math
import os
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import qmc

from .model import FIT_NAMES

CHI_LO, CHI_HI = -8.0, -0.3
H_LO, H_HI = 5.0, 5000.0
N_LO, N_HI = 0.0, 0.8
RHO_LO = 0.5 + 1e-3
# u is ordered in the DECODING order (rho_bar before rho), which is also the order of FitSpec.free
DECODE_ORDER = ("chi", "h", "N", "N_bar", "rho_bar", "rho")


@dataclass
class FitSpec:
    free: tuple = FIT_NAMES
    fixed: dict = field(default_factory=dict)
    bounds: dict = field(default_factory=lambda: dict(chi=(CHI_LO, CHI_HI), h=(H_LO, H_HI), N=(N_LO, N_HI)))

    def __post_init__(self):
        unknown = [n for n in self.free if n not in DECODE_ORDER]
        if unknown:
            raise ValueError(f"unknown parameters {unknown}")
        self.free = tuple(n for n in DECODE_ORDER if n in self.free)
        miss = [n for n in FIT_NAMES if n not in self.free and n not in self.fixed]
        if miss:
            raise ValueError(f"parameters neither free nor fixed: {miss}")

    # ---- the sequential map -------------------------------------------------------------------
    def _intervals_and_values(self, u=None, theta=None):
        """Walk the decoding order. With u: decode -> theta. With theta: encode -> u. Returns (theta, u)."""
        fx = self.fixed
        th, uo = {}, {}
        it = iter(u) if u is not None else None

        def take(name, lo, hi, log=False):
            if name not in self.free:
                th[name] = float(fx[name])
                return
            if not (hi > lo):
                raise ValueError(f"empty interval for {name}: [{lo}, {hi}] given the fixed set {fx}")
            a, b = (math.log(lo), math.log(hi)) if log else (lo, hi)
            if it is not None:
                uu = float(next(it))
                val = a + uu * (b - a)
                th[name] = math.exp(val) if log else val
            else:
                val = math.log(theta[name]) if log else float(theta[name])
                uu = (val - a) / (b - a)
                th[name] = float(theta[name])
            uo[name] = uu

        take("chi", *self.bounds["chi"])
        take("h", *self.bounds["h"], log=True)
        nlo, nhi = self.bounds["N"]
        if "N" in self.free:
            if "N_bar" in fx:
                nlo = max(nlo, fx["N_bar"])
                if "rho" in fx and "rho_bar" in fx:
                    nlo = max(nlo, 1.0 - (1.0 - fx["N_bar"]) * fx["rho"] / fx["rho_bar"])
            elif "rho" in fx and "rho_bar" in fx:
                # N_bar free in [0, 1 - (1-N) rho_bar/rho]: non-empty only for N > 1 - rho/rho_bar
                nlo = max(nlo, 1.0 - fx["rho"] / fx["rho_bar"] + 1e-6)
        take("N", nlo, nhi)
        N = th["N"]
        nbhi = N
        if "rho" in fx and "rho_bar" in fx:
            nbhi = min(nbhi, 1.0 - (1.0 - N) * fx["rho_bar"] / fx["rho"])
        take("N_bar", 0.0, nbhi)
        beta = (1.0 - N) / (1.0 - th["N_bar"])
        rbhi = 1.0 if "rho" not in fx else min(1.0, fx["rho"] / beta)
        take("rho_bar", RHO_LO, rbhi)
        take("rho", max(RHO_LO, beta * th["rho_bar"]), 1.0)
        return th, np.array([uo[n] for n in self.free])

    def decode(self, u) -> dict:
        return self._intervals_and_values(u=u)[0]

    def encode(self, theta: dict) -> np.ndarray:
        return self._intervals_and_values(theta=theta)[1]


def check_admissible(theta: dict, zeta: str = "WW"):
    """The owner's refusal rules, stated once more independently of the map (asserted after decoding)."""
    N, Nb, r, rb = theta["N"], theta["N_bar"], theta["rho"], theta["rho_bar"]
    beta = (1.0 - N) / (1.0 - Nb)
    ok = Nb <= N + 1e-15 and r / rb >= beta * (1 - 1e-15)
    if zeta == "WW":
        ok = ok and 0.5 < r <= 1.0 and 0.5 < rb <= 1.0
    return ok and theta["chi"] < 0.0 and theta["h"] > 0.0


def soft_flags(theta: dict) -> list:
    """Admissible-but-warned conditions (sheet §11.2): reported with every fit result, never refused."""
    out = []
    if theta["rho"] > theta["rho_bar"]:
        out.append(f"rho {theta['rho']:.4g} > rho_bar {theta['rho_bar']:.4g} (violates AB06's psi_c <= phi_c reading)")
    return out


# ---- one least-squares run (top level so it pickles) --------------------------------------------
_CTX = {}


def _init_worker(objective, spec, ls_kw, progress=None):
    _CTX.update(objective=objective, spec=spec, ls_kw=ls_kw, progress=progress)


PROGRESS_EVERY = 25


def _progress(r):
    """Every PROGRESS_EVERY residual evaluations, append "time start evals cost" to the progress file (if any)."""
    path = _CTX.get("progress")
    _CTX["k_eval"] = _CTX.get("k_eval", 0) + 1
    if path and _CTX["k_eval"] % PROGRESS_EVERY == 0:
        with open(path, "a") as f:
            f.write(f"{time.strftime('%H:%M:%S')} start {_CTX.get('start')} evals {_CTX['k_eval']} "
                    f"cost {0.5 * float(np.dot(r, r)):.6e}\n")


FAIL_RESIDUAL = 1.0e3


def _fun(u):
    """Residuals at u. An exception inside the model (it should not happen: O2 never raises for numerical
    trouble and the map only yields admissible sets) returns FAIL_RESIDUAL everywhere and is counted."""
    spec, obj = _CTX["spec"], _CTX["objective"]
    th = spec.decode(np.clip(u, 0.0, 1.0))
    try:
        r = obj.residuals(th)
        _CTX["n_res"] = r.size
        _progress(r)
        return r
    except Exception as e:  # noqa: BLE001
        _CTX["n_fail"] = _CTX.get("n_fail", 0) + 1
        _CTX["last_fail"] = f"{type(e).__name__}: {e}"
        if "n_res" not in _CTX:
            raise
        return np.full(_CTX["n_res"], FAIL_RESIDUAL)


def _one_start(args):
    k, u0 = args
    t0 = time.perf_counter()
    obj = _CTX["objective"]
    n0 = obj.n_evals
    _CTX["n_fail"], _CTX["last_fail"], _CTX["start"], _CTX["k_eval"] = 0, "", k, 0
    res = least_squares(_fun, u0, bounds=(0.0, 1.0), **_CTX["ls_kw"])
    th = _CTX["spec"].decode(res.x)
    return dict(start=k, u0=list(map(float, u0)), u=list(map(float, res.x)), theta=th, cost=float(res.cost),
                nfev=int(res.nfev), njev=int(res.njev or 0), status=int(res.status), message=res.message,
                residual_evals=obj.n_evals - n0, seconds=time.perf_counter() - t0,
                optimality=float(res.optimality), soft_flags=soft_flags(th),
                model_exceptions=_CTX.get("n_fail", 0), last_exception=_CTX.get("last_fail", ""))


# max_nfev counts the residual calls of the iteration only (TRF does not count the 2-point Jacobian's k calls),
# so one start costs up to max_nfev * (k + 1) residual evaluations.
LS_DEFAULT = dict(method="trf", x_scale=1.0, diff_step=1e-4, ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=200)


def multistart(objective, spec: FitSpec, n_starts: int = 8, seed: int = 0, workers: int = 1, extra_starts=(),
               ls_kw: dict | None = None, margin: float = 0.05, progress: str | None = None):
    """Scrambled Sobol starts in [margin, 1 - margin]^k plus extra_starts (theta dicts). Returns the runs sorted
    by cost."""
    ls = dict(LS_DEFAULT)
    ls.update(ls_kw or {})
    k = len(spec.free)
    starts = []
    if n_starts > 0:
        sob = qmc.Sobol(d=k, scramble=True, seed=seed)
        U = sob.random(n_starts)
        starts += list(margin + (1.0 - 2.0 * margin) * U)
    for th in extra_starts:
        starts.append(np.clip(spec.encode(th), 0.0, 1.0))
    jobs = list(enumerate(starts))
    if workers <= 1:
        _init_worker(objective, spec, ls, progress)
        runs = [_one_start(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                                 initargs=(objective, spec, ls, progress)) as ex:
            runs = list(ex.map(_one_start, jobs))
    runs.sort(key=lambda r: r["cost"])
    return runs


# ---- identifiability ------------------------------------------------------------------------------
def jacobian_log(objective, theta: dict, names, rel: float = 1e-4):
    """J[:, j] = d r / d ln|theta_j| by central differences (step rel in ln|theta_j|)."""
    cols = []
    for n in names:
        tp, tm = dict(theta), dict(theta)
        tp[n] = theta[n] * math.exp(rel)
        tm[n] = theta[n] * math.exp(-rel)
        cols.append((objective.residuals(tp) - objective.residuals(tm)) / (2.0 * rel))
    return np.array(cols).T


def identifiability(objective, theta: dict, names, rel: float = 1e-4):
    J = jacobian_log(objective, theta, names, rel)
    s = np.linalg.svd(J, compute_uv=False)
    JTJ = J.T @ J
    cov = np.linalg.pinv(JTJ)
    sd = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    corr = cov / np.outer(sd, sd) if np.all(sd > 0) else np.full_like(cov, np.nan)
    _, _, Vt = np.linalg.svd(J, full_matrices=False)
    return dict(names=list(names), singular_values=s.tolist(), condition=float(s[0] / s[-1]) if s[-1] > 0 else math.inf,
                rel_sd_unit_residual=dict(zip(names, sd.tolist())), correlation=corr.tolist(),
                weakest_direction=dict(zip(names, Vt[-1].tolist())), n_residuals=int(J.shape[0]))


def profile(objective, spec: FitSpec, theta_best: dict, name: str, factors=(0.8, 0.9, 0.95, 1.05, 1.1, 1.2),
            ls_kw: dict | None = None, workers: int = 1, progress: str | None = None):
    """Profile scan: fix `name` at theta_best[name] * f, re-fit the other free parameters from theta_best.
    With a small max_nfev in ls_kw the profile cost is an UPPER bound of the true profile (stated with the result)."""
    out = []
    jobs = []
    for f in factors:
        val = theta_best[name] * f
        fixed = dict(spec.fixed)
        fixed[name] = val
        free = tuple(n for n in spec.free if n != name)
        sub = FitSpec(free=free, fixed=fixed, bounds=spec.bounds)
        try:
            u0 = np.clip(sub.encode(theta_best), 0.0, 1.0)
        except ValueError as e:
            out.append(dict(factor=f, value=val, cost=None, note=str(e)))
            continue
        jobs.append((f, val, sub, u0, progress))
    ls = dict(LS_DEFAULT)
    ls.update(ls_kw or {})
    if workers <= 1:
        res = [_profile_point(objective, j, ls) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            res = list(ex.map(_profile_point, [objective] * len(jobs), jobs, [ls] * len(jobs)))
    out += res
    out.sort(key=lambda r: r["factor"])
    return out


def _profile_point(objective, job, ls):
    f, val, sub, u0, progress = job
    _init_worker(objective, sub, ls, progress)
    r = _one_start((f"profile x{f}", u0))
    return dict(factor=f, value=val, cost=r["cost"], theta=r["theta"], nfev=r["nfev"], seconds=r["seconds"])
