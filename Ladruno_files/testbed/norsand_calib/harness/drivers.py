"""Element drivers on O2 (the C++ contract) and their O1 counterparts (task item 1).

Axes (principal, fixed; no shear anywhere): index 0 = axial (strain controlled), 1 and 2 lateral.
Compression negative inside; the Curve accessors report the lab convention (compression positive).

  PS   drained plane strain compression: eps_11 = 0 (index 1, the plane-strain axis, sigma_2 intermediate),
       sigma_22 = sigma_3' held at its initial value (index 2), eps_00 prescribed (compression).
  TC   drained triaxial compression: sigma_11 = sigma_22 held, eps_00 prescribed (compression).
  TE   drained triaxial extension, lateral stress held: eps_00 prescribed (extension), sigma_11 = sigma_22 held
       (axial unloading; the axial stress becomes the minor principal stress).
  PSU  constant-volume plane strain (TIMs rung C6): d_eps = (da, 0, -da).
  TCU  constant-volume triaxial: d_eps = (da, -da/2, -da/2).

O2 mixed control (drained kinds): per increment, the unknown lateral strain increment x (one scalar: the same on
every stress-controlled index T, which all drained kinds allow by symmetry) solves sigma_T(x) = sigma_T,0 by
Newton with O2's step and its consistent tangent (tangent() of the trial state: the chained tangent when O2
substepped, sheet §9.6), from the tangent predictor of the last converged state; when Newton stalls (O2's
response is piecewise smooth in x: the substep ladder can switch between trials) it falls back to a bracket +
Brent solve (_solve_lateral). Converged when |sigma_T - sigma_T,0| <= MIX_TOL_REL * p_init. Every trial restarts
from the committed state. Failure statuses: 'refused:<O2 reason>', 'mixed_nobracket', 'mixed_jump'.

O1 counterpart: o1_rate.integrate_increment with smask (zero stress rate on T, Mandel order xx,yy,zz,...),
i.e. the same mixed problem solved inside the rate equations (Radau, rtol 1e-10).
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from .model import oracle_module

MIX_TOL_REL = 1.0e-11

KINDS = {
    "PS": dict(drained=True, sign=-1.0, T=(2,), S_ratio=(1.0, 0.0, None)),
    "TC": dict(drained=True, sign=-1.0, T=(1, 2), S_ratio=(1.0, None, None)),
    "TE": dict(drained=True, sign=+1.0, T=(1, 2), S_ratio=(1.0, None, None)),
    "PSU": dict(drained=False, sign=-1.0, T=(), S_ratio=(1.0, 0.0, -1.0)),
    "TCU": dict(drained=False, sign=-1.0, T=(), S_ratio=(1.0, -0.5, -0.5)),
}


@dataclass
class Curve:
    kind: str
    eps: np.ndarray            # (m+1, 3) principal total strain, compression negative
    sig: np.ndarray            # (m+1, 3) principal stress, compression negative
    pi_i: np.ndarray           # (m+1,)
    v: np.ndarray              # (m+1,)
    D: np.ndarray              # (m+1,) dissipation of the increment (0 at the start)
    status: str = "ok"
    n_target: int = 0
    stats: dict = field(default_factory=dict)

    @property
    def complete(self) -> bool:
        return self.status == "ok" and len(self.eps) == self.n_target + 1

    @property
    def eps_a_pct(self):
        """Axial strain, % , compression positive (TE: negative)."""
        return -100.0 * self.eps[:, 0]

    @property
    def eps_v_pct(self):
        """Volumetric strain, %, compression positive (dilation negative, Tatsuoka's figures)."""
        return -100.0 * self.eps.sum(axis=1)

    @property
    def sr(self):
        """sigma1'/sigma3' = major/minor principal stress magnitude."""
        s = -self.sig
        return s.max(axis=1) / s.min(axis=1)

    @property
    def p(self):
        return -self.sig.mean(axis=1)

    @property
    def q(self):
        s = self.sig
        return np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2 + (s[:, 2] - s[:, 0]) ** 2))


def _diag(st):
    return np.array([st.sigma[0, 0], st.sigma[1, 1], st.sigma[2, 2]])


def _tangent_TT(P, st, T, S):
    from o2_algo import api as A
    C = A.tangent(P, st)
    Cd = np.array([[C[i, i, j, j] for j in range(3)] for i in range(3)])
    return Cd[np.ix_(T, T)], Cd[np.ix_(T, S)]


def run_o2(P, st0, kind: str, eps_a_total: float, n: int) -> Curve:
    """O2 element test: n equal increments of the axial strain magnitude eps_a_total (> 0; the sign comes from
    the kind). Stops at the first failure (status 'refused:<reason>' / 'mixed_nobracket' / 'mixed_jump')."""
    from o2_algo import api as A
    spec = KINDS[kind]
    da = spec["sign"] * abs(eps_a_total) / n
    T = list(spec["T"])
    S = [i for i in range(3) if i not in T]
    p_init = abs(float(np.trace(st0.sigma)) / 3.0)
    tol = MIX_TOL_REL * p_init
    target = _diag(st0)[T] if T else None
    st = st0
    eps = [np.zeros(3)]
    sig = [_diag(st0)]
    pi = [st0.pi_i]
    v = [st0.v]
    Dl = [0.0]
    status = "ok"
    n_steps = n_newton = n_sub = n_fallback = 0
    max_jump = 0.0
    sub_at, brk_at = [], []
    t0 = time.perf_counter()
    for _ in range(n):
        d = np.zeros(3)
        ratio = spec["S_ratio"]
        for i in S:
            d[i] = ratio[i] * da
        if not T:
            new = A.step(P, st, np.diag(d))
            n_steps += 1
        else:
            new, x, info = _solve_lateral(P, st, d, T, S, target, tol, p_init)
            n_steps += info["steps"]
            n_newton += info["steps"]
            n_fallback += info["fallback"]
            max_jump = max(max_jump, info["jump"])
            d[T] = x
            if info["fail"]:
                status = info["fail"]
                break
            if new.flags["refused"]:
                status = "refused:" + new.flags["reason"]
                break
        if new.flags.get("substeps", 1) > 1:
            n_sub += 1
            sub_at.append(len(eps) - 1)
        if T and info["fallback"]:
            brk_at.append(len(eps) - 1)
        st = new
        eps.append(eps[-1] + d)
        sig.append(_diag(st))
        pi.append(st.pi_i)
        v.append(st.v)
        Dl.append(st.D)
    return Curve(kind, np.array(eps), np.array(sig), np.array(pi), np.array(v), np.array(Dl), status, n,
                 dict(o2_steps=n_steps, newton=n_newton, substepped_incr=n_sub, bracketed_incr=n_fallback,
                      max_mixed_residual_rel=max_jump, substepped_at=sub_at, bracketed_at=brk_at,
                      seconds=time.perf_counter() - t0, final_state=st))


def _solve_lateral(P, st, d, T, S, target, tol, p_init):
    """The drained mixed-control solve of one increment. Every drained kind has ONE unknown: the lateral strain
    increment x, the same on every stress-controlled index (PS: one index; TC/TE: the two lateral indices with
    equal targets, so the axisymmetric solution has equal lateral strains). f(x) = mean_T(sigma_T(x)) - target.
      1. Newton with the slope J = sum_b C[t, t, b, b] (O2's consistent tangent of the trial state), from the
         tangent predictor; at most NEWTON_MAX iterations, each must reduce |f|.
      2. Otherwise (O2's response is only piecewise smooth in x: the substep ladder can switch between trials):
         bracket the sign change by expanding steps from the best x, then scipy brentq to x round-off.
         If the bracket holds a jump instead of a root (no |f| <= tol), the end with the smaller |f| is taken when
         |f| <= JUMP_ACCEPT_REL * p_init (counted: stats max_mixed_residual_rel), else the increment fails.
    Returns (state, x, info)."""
    from o2_algo import api as A
    from scipy.optimize import brentq
    t = T[0]
    info = dict(steps=0, fallback=0, jump=0.0, fail="")
    cache = {}

    def ev(x):
        if x in cache:
            return cache[x]
        dd = d.copy()
        dd[T] = x
        new = A.step(P, st, np.diag(dd))
        info["steps"] += 1
        f = None if new.flags["refused"] else float(np.mean(_diag(new)[T]) - target[0])
        cache[x] = (f, new)
        return f, new

    def ev_strict(z):
        f, new = ev(z)
        if f is None:
            raise _Refused(new)
        return f

    def slope(state):
        C = A.tangent(P, state)
        return float(sum(C[t, t, b, b] for b in T))

    CTT, CTS = _tangent_TT(P, st, T, S)
    J0 = float(CTT.sum(axis=1)[0])
    x = (target[0] - _diag(st)[t] - float((CTS @ d[S])[0])) / J0
    try:
        f, new = ev(x)
        if f is None:                     # the predictor itself refused: retry from x = 0 (no lateral strain)
            x = 0.0
            f = ev_strict(x)
            new = cache[x][1]
        best = (abs(f), x, f, new)
        for _ in range(NEWTON_MAX):
            if abs(f) <= tol:
                return new, x, info
            J = slope(new)
            if not (J > 0.0):
                break
            xn = x - f / J
            fn, nn = ev(xn)
            if fn is None or abs(fn) >= abs(f):
                break
            x, f, new = xn, fn, nn
            best = (abs(f), x, f, new)
        if abs(f) <= tol:
            return new, x, info
        # 2. bracket + Brent
        info["fallback"] = 1
        _, x, f, new = best
        J = slope(new)
        h = max(abs(f) / J if J > 0.0 else 0.0, 1e-9)
        dirn = -1.0 if f > 0.0 else 1.0
        a, fa = x, f
        b = fb = None
        for _ in range(60):
            xb = a + dirn * h
            fb_, nb = ev(xb)
            if fb_ is None:               # too far: a refused trial; shrink the expansion
                h *= 0.25
                continue
            if fb_ == 0.0 or (fb_ > 0.0) != (fa > 0.0):
                b, fb = xb, fb_
                break
            if abs(fb_) < abs(fa):
                a, fa = xb, fb_
            h *= 2.0
        if b is None:
            info["fail"] = "mixed_nobracket"
            return new, x, info
        xr = brentq(ev_strict, min(a, b), max(a, b), xtol=1e-18, rtol=4.0 * np.finfo(float).eps,
                    maxiter=200)
        fr, nr = ev(xr)
        if abs(fr) <= tol:
            return nr, xr, info
        cand = min(((abs(v[0]), z, v[1]) for z, v in cache.items() if v[0] is not None), key=lambda c: c[0])
        info["jump"] = cand[0] / p_init
        if cand[0] <= JUMP_ACCEPT_REL * p_init:
            return cand[2], cand[1], info
        info["fail"] = "mixed_jump"
        return cand[2], cand[1], info
    except _Refused as e:
        info["fail"] = "refused:" + e.state.flags["reason"]
        return e.state, x, info


class _Refused(Exception):
    def __init__(self, state):
        super().__init__(state.flags["reason"])
        self.state = state


NEWTON_MAX = 12
# A bracket that holds a jump (no root) is accepted when the residual left is <= this fraction of p_init. The jump
# is O2's: its response to the increment is discontinuous where the substep ladder changes level (measured
# 2.8e-4 and 1.3e-3 p_init at sigma3' 4.9 kPa from the apex start, Esmeralda 2026-10-02; with the default
# ramp_end start the diagnosis curves (the truth and the stalled point of the first smoke) held none). Accepting it costs a stress error of that size at one
# increment (the order of the step's own discretisation error) instead of a refused curve.
JUMP_ACCEPT_REL = 1.0e-3


def run_o1(P, st0, kind: str, eps_a_total: float, n: int, rtol: float = 1e-10) -> Curve:
    """O1 counterpart of run_o2 (the truth; the endpoint does not depend on n)."""
    from o1_rate.integrator import integrate_increment
    spec = KINDS[kind]
    da = spec["sign"] * abs(eps_a_total) / n
    T = list(spec["T"])
    S = [i for i in range(3) if i not in T]
    smask = [i in T for i in range(3)] + [False, False, False]
    st = st0
    eps = [np.zeros(3)]
    sig = [_diag(st0)]
    pi = [st0.pi_i]
    v = [st0.v]
    Dl = [0.0]
    status = "ok"
    t0 = time.perf_counter()
    for _ in range(n):
        d = np.zeros(3)
        for i in S:
            d[i] = spec["S_ratio"][i] * da
        st = integrate_increment(P, st, np.diag(d), smask=smask if T else None, rtol=rtol)
        if st.flags["status"] != "ok":
            status = st.flags["status"]
            break
        e_tot = np.diag(st.flags["eps_total"]).copy()
        eps.append(e_tot)
        sig.append(_diag(st))
        pi.append(st.pi_i)
        v.append(st.v)
        Dl.append(st.D)
    return Curve(kind, np.array(eps), np.array(sig), np.array(pi), np.array(v), np.array(Dl), status, n,
                 dict(seconds=time.perf_counter() - t0, final_state=st))


def simulate(setup, theta: dict, kind: str, sigma3: float, e_init: float, eps_a_total: float, n: int,
             oracle: str = "O2") -> Curve:
    """One element test from the isotropic state sigma3 (kPa > 0), void ratio e_init (v0 = 1 + e)."""
    from .model import initial
    P, st0 = initial(oracle, setup, theta, sigma3, e_init)
    if oracle == "O2":
        return run_o2(P, st0, kind, eps_a_total, n)
    return run_o1(P, st0, kind, eps_a_total, n)
