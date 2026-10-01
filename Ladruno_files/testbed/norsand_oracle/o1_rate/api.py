"""Shared-interface entry points of O1 (initial_state, run_path, tangent, triaxial)."""
from __future__ import annotations

import math

import numpy as np

from .integrator import State, integrate_increment
from .model import (continuum_tangent, elastic_strain_from_stress, energy, invariants, plastic,
                    zeta_fun, SQ32, SQ6)
from .params import Params


def initial_state(params: Params, sigma0, v0: float, pi_i0: float | None = None) -> State:
    """State at sigma0 (3x3, compression negative), specific volume v0.

    pi_i0 None: pi_i is put ON the yield surface through sigma0, by the inverse of (S.12)
    (BA06 2.8) at eta = -zeta(theta) q / p.  eps^e is the inverse of the energy (co-axial).
    """
    P = params
    sigma0 = np.asarray(sigma0, dtype=float)
    ee = elastic_strain_from_stress(sigma0, P)
    el = energy(ee, P, tangent=False)
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
        if P.N == 0.0:
            pi = p * math.exp(eta / P.M - 1.0)
        else:
            if not eta < P.M / P.N:
                raise ValueError(f"eta = {eta} >= M/N: no yield surface through sigma0")
            pi = p * ((1.0 - P.N) / (1.0 - eta * P.N / P.M)) ** ((1.0 - P.N) / P.N)
    else:
        pi = float(pi_i0)
    flags = dict(v0=float(v0), eps_total=np.zeros((3, 3)), W=0.0, W_abs=0.0, Dp_total=0.0,
                 status="initial", plastic=False, H_zero=[])
    pq = plastic(el.sig, pi, v0, P)
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


def tangent(params: Params, state: State, plastic_branch: bool | None = None):
    """CONTINUUM tangent (3,3,3,3): a^ep of (S.42) (loading branch) if the last increment
    ended in the plastic mode (state.flags['plastic']), else a^e (S.3)/(S.33)."""
    if plastic_branch is None:
        plastic_branch = bool(state.flags.get("plastic", False))
    C, _ = continuum_tangent(state.eps_e, state.pi_i, state.v, params, plastic_branch)
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
