"""Shared-interface entry points of O1 (initial_state, run_path, tangent, triaxial)."""
from __future__ import annotations

import math

import numpy as np

from .integrator import State, integrate_increment
from .model import (continuum_tangent, elastic_strain_from_stress, energy, floor_project, invariants,
                    plastic, zeta_fun, SQ32, SQ6, R_TOL_REL)
from .params import Params


def pi_of_eta(P: Params, p: float, eta: float) -> float:
    """Inverse of (S.12) at (p, eta): the image pressure of the yield surface through stress ratio eta."""
    if P.N == 0.0:
        return p * math.exp(eta / P.M - 1.0)
    if not eta < P.M / P.N:
        raise ValueError(f"eta = {eta} >= M/N: no yield surface through the state")
    return p * ((1.0 - P.N) / (1.0 - eta * P.N / P.M)) ** ((1.0 - P.N) / P.N)


def initial_state(params: Params, sigma0, v0: float, pi_i0: float | None = None,
                  pi_rule: str = "surface") -> State:
    """State at sigma0 (3x3, compression negative), specific volume v0.

    eps^e is the inverse of the energy (co-axial; (S.5) Newton under BA06, (S.5h'') closed form under
    HAR), then the p' floor Pi_f (S.48) when p_min > 0 (counted: n_f_init, eps_f_v; sheet 9.7
    'initialState').  Every rule below uses the floored stress.
    pi_i0 given: used as is (the deck's -pi0 overrides any rule).
    pi_i0 None, pi_rule = "surface" (pre-round-3 O1 default): the yield surface through the stress,
      eta = zeta(theta) q / |p| ((S.12) inverse, BA06 2.8).
    pi_i0 None, pi_rule = "S53": the unified rule (S.53), owner decision (c)/(d): the surface through
      (p_init, eta*), eta* = max(eta_init, c2 M) with c2 = c2 (smooth cap), c1 (planar, chi_cap), 0 (no
      cap; then identical to "surface"); eta_init = 0 on the axis (3.2); refused if eta* >= M/N; the
      B > 0 guard of section 7 is checked at (pi_i0, psi_i0).
    """
    P = params
    sigma0 = np.asarray(sigma0, dtype=float)
    ee = elastic_strain_from_stress(sigma0, P)
    ee, f_init, dvf_init, _ = floor_project(ee, P)
    el = energy(ee, P, tangent=False)
    if pi_rule not in ("surface", "S53"):
        raise ValueError(f"pi_rule must be 'surface' or 'S53', got {pi_rule!r}")
    if pi_i0 is None:
        p, s, R = invariants(el.sig)
        q = SQ32 * R
        if R > 0:
            y = float(np.sum(s * (s @ s))) / R ** 3
        else:
            y = -1.0 / SQ6
        theta = math.acos(min(1.0, max(-1.0, SQ6 * y))) / 3.0
        z, _ = zeta_fun(theta, y, P.rho, P.zeta)
        eta = -z * q / p
        if pi_rule == "S53":
            if R < R_TOL_REL * abs(p):
                eta = 0.0                                   # the axis (3.2)
            c2 = {"none": 0.0, "planar": P.c1, "smooth": P.c2}[P.cap]
            eta = max(eta, c2 * P.M)
        pi = pi_of_eta(P, p, eta)
    else:
        pi = float(pi_i0)
    flags = dict(v0=float(v0), eps_total=np.zeros((3, 3)), W=0.0, W_abs=0.0, Dp_total=0.0,
                 status="initial", plastic=False, H_zero=[])
    if P.pmin > 0.0:
        flags.update(p_min=P.pmin, n_f_init=int(f_init), eps_f_v=dvf_init, eps_f_v_init=dvf_init,
                     W_f=P.pmin * dvf_init, at_floor=bool(el.p > -P.pmin * (1.0 + 1e-10)), floor=False,
                     n_floor_increments=0)
    pq = plastic(el.sig, pi, v0, P)
    if pi_i0 is None and pi_rule == "S53" and P.N > 0 and not pq.B > 0.0:
        raise ValueError(f"S53 initial state refused: B = {pq.B} <= 0 at (pi_i0, psi_i0) (section 7 guard)")
    flags.update(F_rel=pq.F / (P.M * abs(pq.p)), psi=pq.psi, pistar=pq.pistar, H=pq.H)
    return State(sigma=el.sig, eps_e=ee, pi_i=pi, v=float(v0), D=0.0, eps_p_v=0.0,
                 eps_p_s=0.0, flags=flags)


def run_path(params: Params, state0: State, deps, rtol: float = 1e-10,
             kin: str = "small") -> list[State]:
    """Apply n strain increments deps (n,3,3) (full symmetric tensors, tensor shear), each a
    linear ramp over unit pseudo-time integrated by Radau.  Returns the state after each."""
    deps = np.asarray(deps, dtype=float)
    out = []
    st = state0
    for k in range(deps.shape[0]):
        d = 0.5 * (deps[k] + deps[k].T)
        st = integrate_increment(params, st, d, rtol=rtol, kin=kin)
        st.flags["step"] = k + 1
        out.append(st)
        if st.flags["status"] != "ok":
            break
    return out


def tangent(params: Params, state: State, plastic_branch: bool | None = None,
            floor_branch: bool | None = None):
    """CONTINUUM tangent (3,3,3,3): a^ep of (S.42) (loading branch) if the last increment
    ended in the plastic mode (state.flags['plastic']), else a^e (S.3)/(S.33).  If it ended on the
    floor (state.flags['floor'], sheet 9.7 / S.55) the floor mechanism is added: the exact one-sided
    linearisation, with 1:C = 0 (no bulk stiffness at a floored state; owner decision (c))."""
    if plastic_branch is None:
        plastic_branch = bool(state.flags.get("plastic", False))
    if floor_branch is None:
        floor_branch = bool(state.flags.get("floor", False))
    C, _ = continuum_tangent(state.eps_e, state.pi_i, state.v, params, plastic_branch, floor_branch)
    return C


def triaxial(params: Params, state0: State, kind: str, axial_strain_total: float,
             n_incr: int, rtol: float = 1e-10) -> list[State]:
    """Single-point triaxial, axial direction x (index 0).

    axial_strain_total is the SIGNED total axial tensor strain (compression negative, e.g.
    -0.2 for 20 % compression; positive = extension).
      drained   : eps_xx prescribed, sigma_yy' = sigma_zz' = 0 (lateral stress constant),
                  shear strains zero (mixed control inside the rate equations).
      undrained : isochoric, eps_yy = eps_zz = -eps_xx/2, all strain controlled.
    """
    if kind not in ("drained", "undrained"):
        raise ValueError("kind must be 'drained' or 'undrained'")
    da = axial_strain_total / n_incr
    out = []
    st = state0
    for k in range(n_incr):
        if kind == "undrained":
            d = np.diag([da, -0.5 * da, -0.5 * da])
            st = integrate_increment(params, st, d, rtol=rtol)
        else:
            d = np.diag([da, 0.0, 0.0])
            smask = [False, True, True, False, False, False]
            st = integrate_increment(params, st, d, smask=smask, rtol=rtol)
        st.flags["step"] = k + 1
        out.append(st)
        if st.flags["status"] != "ok":
            break
    return out
