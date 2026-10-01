"""Shared O1/O2 interface: State, initial_state, run_path, tangent, drivers.

Tensor convention: 3x3 symmetric numpy arrays, compression negative, TENSOR strains
(no engineering shear anywhere inside; sheet §1.2).

Small-strain mode (default): eps^e_tr = eps^e_n + d_eps (any symmetric d_eps, non-coaxial
allowed), v = v0 (1 + tr eps), d v / d eps~_b = v0 (S.31).
Finite-strain mode (State.finite = True, used only by k2_path): the increments are principal
LOG-strain increments of a DIAGONAL, fixed-direction deformation (the (S.43) protocol), for
which b^e,tr = f b^e_n f^T reduces to eps~_a = eps^e_n,a + ln f_a; v = v0 J = v0 exp(tr eps),
d v / d eps~_b = v (§1.4). The return map is byte-identical (§1.4); only vfac and the
tangent assembly (S.34) differ.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import kernel as K
from .params import Params


@dataclass
class State:
    sigma: np.ndarray                     # 3x3, compression negative
    eps_e: np.ndarray                     # 3x3 elastic strain
    pi_i: float                           # image pressure (< 0)
    v: float                              # specific volume
    D: float = 0.0                        # plastic dissipation of the last step (S.38)
    eps_p_v: float = 0.0                  # accumulated plastic volumetric strain
    eps_p_s: float = 0.0                  # accumulated plastic deviatoric strain (sum of dlam sqrt(2/3) Omega)
    flags: dict = field(default_factory=lambda: dict(plastic=False, vertex=False, cap_active=False,
                                                     refused=False, local_iters=0, pi_iters=0, reason=""))
    eps_p: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))   # accumulated plastic strain tensor
    v0: float = 0.0
    finite: bool = False
    # diagnostics of the last step (converged values; 0 on elastic steps where meaningless)
    dlam: float = 0.0
    eta: float = 0.0
    psi_i: float = 0.0
    pi_star: float = 0.0
    F_p: float = 0.0
    Omega: float = 0.0
    cache: dict = field(default_factory=dict, repr=False)   # spectral data for tangent()

    def copy(self) -> "State":
        return State(self.sigma.copy(), self.eps_e.copy(), self.pi_i, self.v, self.D, self.eps_p_v,
                     self.eps_p_s, dict(self.flags), self.eps_p.copy(), self.v0, self.finite, self.dlam,
                     self.eta, self.psi_i, self.pi_star, self.F_p, self.Omega, dict(self.cache))


def _eigh(A: np.ndarray):
    w, V = np.linalg.eigh(0.5 * (A + A.T))
    return w, V


def _from_principal(vals: np.ndarray, V: np.ndarray) -> np.ndarray:
    return (V * vals) @ V.T


def initial_state(params: Params, sigma0: np.ndarray, v0: float, pi_i0: float | None,
                  finite: bool = False) -> State:
    """State at sigma0 with eps^p = 0. pi_i0 None: pi_i placed so that F(sigma0, pi_i) = 0
    (inverse of (S.12); hydrostatic sigma0 => pi_i = p (1-N)^((1-N)/N), the apex through p)."""
    params.validate()
    sig = np.asarray(sigma0, float)
    w, V = _eigh(sig)
    eps_e_p = K.invert_elastic(params, w)
    inv = K.invariants(w)
    if pi_i0 is None:
        if inv.vertex:
            eta = 0.0
        else:
            z, _, _ = K.zeta_y(inv.theta, params.rho, params.zeta)
            eta = -z * inv.q / inv.p
        pi_i0 = K.pi_of_eta(params, inv.p, eta)
    st = State(_from_principal(w, V), _from_principal(eps_e_p, V), float(pi_i0), float(v0),
               v0=float(v0), finite=finite)
    psi, _ = K.csl(params, v0, pi_i0)
    st.psi_i = psi
    try:
        st.eta = K.eta_of(params, inv.p, pi_i0)
    except ValueError:
        pass
    st.cache = dict(eps_tr=eps_e_p, nvec=V, sig=w, atilde=K.elastic(params, eps_e_p).ae)
    return st


def step(params: Params, st: State, deps: np.ndarray) -> State:
    """One backward-Euler increment (sheet §9.1) with bounded substepping. Never raises for
    numerical trouble.

    CONTRACT (substepping): the increment is first attempted whole. If the return map refuses
    it (any refusal reason), it is retried as 2, 4, ..., 2^MAX_SUBSTEP_HALVINGS equal
    sub-increments, each a full backward-Euler step chained on the previous one (own spectral
    decomposition, own trial check, v composed exactly). The first level at which every
    sub-increment is accepted is the answer; flags['substeps'] records that count (1 = no
    substepping). If the finest level is still refused, the state returned is the frozen
    trial-elastic state of the WHOLE increment, flags['refused'] = True and reason =
    "<finest reason> (substeps exhausted at 2^MAX_SUBSTEP_HALVINGS)".
    On a substepped increment the consistent tangent (cache -> tangent()) is the closed-form
    CTO of the LAST sub-increment: a consistent linearisation of that sub-step, not of the whole
    increment (the chain rule through the earlier sub-steps is not assembled). D is the sum of
    the sub-increment dissipations (each >= 0), plastic = any sub-increment plastic,
    local_iters / pi_iters are summed, vertex / cap_active are the last sub-increment's."""
    deps = np.asarray(deps, float)
    first = None
    for j in range(K.MAX_SUBSTEP_HALVINGS + 1):
        m = 2 ** j
        cur = st
        D, iters, piters, plastic, ok = 0.0, 0, 0, False, True
        for _ in range(m):
            cur = _step_once(params, cur, deps / m)
            if cur.flags["refused"]:
                ok = False
                break
            D += cur.D
            iters += cur.flags["local_iters"]
            piters += cur.flags["pi_iters"]
            plastic = plastic or cur.flags["plastic"]
        if ok:
            cur.D = D
            cur.flags.update(local_iters=iters, pi_iters=piters, plastic=plastic, substeps=m)
            return cur
        if first is None:
            first = cur
        last_reason = cur.flags["reason"]
    first.flags["reason"] = f"{last_reason} (substeps exhausted at 2^{K.MAX_SUBSTEP_HALVINGS})"
    first.flags["substeps"] = 2 ** K.MAX_SUBSTEP_HALVINGS
    return first


def _step_once(params: Params, st: State, deps: np.ndarray) -> State:
    """One backward-Euler increment without substepping."""
    tr = float(np.trace(deps))
    if st.finite:
        v = st.v * math.exp(tr)
        vfac = v
    else:
        v = st.v + st.v0 * tr
        vfac = st.v0
    eps_tr = st.eps_e + 0.5 * (deps + deps.T)
    w, V = _eigh(eps_tr)
    res = K.return_map(params, w, st.pi_i, v, vfac)
    new = State(_from_principal(res.sig, V), _from_principal(res.eps_e, V), res.pi, v,
                v0=st.v0, finite=st.finite, eps_p=st.eps_p.copy(),
                eps_p_v=st.eps_p_v, eps_p_s=st.eps_p_s)
    new.flags = dict(plastic=res.plastic, vertex=res.vertex, cap_active=res.cap_active,
                     refused=res.refused, local_iters=res.local_iters, pi_iters=res.pi_iters,
                     reason=res.reason, substeps=1)
    if res.refused:
        # freeze: return the trial-elastic stress with the refusal flag; the driver stops.
        new.cache = dict(eps_tr=w, nvec=V, sig=res.sig, atilde=res.atilde)
        return new
    new.D = res.D
    new.dlam = res.dlam
    new.eta = res.eta
    new.psi_i = res.psi
    new.pi_star = res.ps
    new.F_p = res.F_p
    new.Omega = res.Om
    if res.plastic:
        dep = res.dlam * _from_principal(res.q_a, V)
        new.eps_p = st.eps_p + dep
        new.eps_p_v = st.eps_p_v + res.dlam * float(res.q_a.sum())
        new.eps_p_s = st.eps_p_s + res.dlam * K.SQ23 * res.Om
    new.cache = dict(eps_tr=w, nvec=V, sig=res.sig, atilde=res.atilde)
    return new


def run_path(params: Params, state0: State, deps: np.ndarray) -> list[State]:
    """n increments -> n states. After a refusal the remaining entries are copies of the
    refused state (flags.refused True); nothing is raised."""
    params.validate()
    deps = np.asarray(deps, float)
    out = []
    st = state0
    for i in range(deps.shape[0]):
        if st.flags.get("refused", False):
            out.append(st.copy())
            continue
        st = step(params, st, deps[i])
        out.append(st)
    return out


def tangent(params: Params, state: State) -> np.ndarray:
    """Small-strain consistent tangent (S.33) of the last step (elastic a^e if the last step
    was elastic, or for a fresh initial state)."""
    c = state.cache
    return K.tangent_small(c["atilde"], c["sig"], c["eps_tr"], c["nvec"])


def tangent_finite(params: Params, state: State) -> np.ndarray:
    """Finite-strain spatial tangent a^ep = c~ + tau(+)1 (S.34) of the last step; the state
    must have been produced in finite mode (eps~ = log stretches, v = v0 J)."""
    c = state.cache
    return K.tangent_finite(c["atilde"], c["sig"], c["eps_tr"], c["nvec"])


# --------------------------------------------------------------------------------------
# drivers (sheet §13 K1.6-K1.8, §14)
# --------------------------------------------------------------------------------------
TRIAX_LAT_TOL_REL = 1.0e-11   # drained: |sigma_lat - sigma_lat0| <= tol*|p0| per increment
TRIAX_LAT_MAX_ITERS = 40


def triaxial(params: Params, state0: State, kind: str, axial_strain_total: float, n_incr: int) -> list[State]:
    """Axisymmetric triaxial about axis 3 (eps_11 = eps_22 = lateral), n_incr equal axial-strain
    increments. axial_strain_total < 0 is compression (TXC, theta = pi/3), > 0 extension.
    kind='drained': sigma_11 = sigma_22 held at their initial value (Newton on the lateral
    strain with the consistent tangent); kind='undrained': isochoric, eps_lat = -eps_ax/2."""
    params.validate()
    states = []
    st = state0
    da = axial_strain_total / n_incr
    sig_lat0 = state0.sigma[0, 0]
    for _ in range(n_incr):
        if kind == "undrained":
            dl = -0.5 * da
            st = step(params, st, np.diag([dl, dl, da]))
        elif kind == "drained":
            dl = 0.0
            if st.cache:
                Ct = tangent(params, st)
                k = Ct[0, 0, 0, 0] + Ct[0, 0, 1, 1]
                dl = -(Ct[0, 0, 2, 2] * da) / k          # elastic/tangent predictor
            ok = False
            for _it in range(TRIAX_LAT_MAX_ITERS):
                trial = step(params, st, np.diag([dl, dl, da]))
                if trial.flags["refused"]:
                    break
                res = trial.sigma[0, 0] - sig_lat0
                if abs(res) <= TRIAX_LAT_TOL_REL * abs(params.p0):
                    ok = True
                    break
                Ct = tangent(params, trial)
                k = Ct[0, 0, 0, 0] + Ct[0, 0, 1, 1]
                dl -= res / k
            st = trial
            if not ok and not st.flags["refused"]:
                st.flags["refused"] = True
                st.flags["reason"] = "triaxial_lateral_noconv"
        else:
            raise ValueError("kind must be 'drained' or 'undrained'")
        states.append(st)
        if st.flags["refused"]:
            break
    return states


def k2_path(params: Params, state0: State, n_max: int, n1: int = 10, lam1: float = 1.0e-3,
            lam2: float = 4.0e-4, acoustic_kwargs: dict | None = None):
    """Sheet §14 (S.43) protocol in finite strain: f1 = diag(1+lam2, 1-lam1, 1) for n1 steps,
    then f2 = diag(1, 1-lam2, 1+lam1) until localization or n_max. Returns a dict with the
    states, the per-step min det of the finite-strain acoustic tensor (S.44), its
    normalisation by the step-n1 value, and the first-localization step under both criteria:
      n_first  : first step n with min_det <= 0
      n_interp : linear-interpolated zero crossing of the normalised curve (fractional)
    None when no crossing within n_max."""
    from .acoustic import acoustic_min_det
    params.validate()
    st = state0
    if not st.finite:
        st = st.copy()
        st.finite = True
    f1 = np.log(np.array([1.0 + lam2, 1.0 - lam1, 1.0]))
    f2 = np.log(np.array([1.0, 1.0 - lam2, 1.0 + lam1]))
    states, dets, nvecs = [], [], []
    kw = acoustic_kwargs or {}
    n_first = None
    for n in range(1, n_max + 1):
        d = f1 if n <= n1 else f2
        st = step(params, st, np.diag(d))
        states.append(st)
        if st.flags["refused"]:
            break
        a4 = tangent_finite(params, st)
        md, nv = acoustic_min_det(params, st, a4, **kw)
        dets.append(md)
        nvecs.append(nv)
        if md <= 0.0 and n_first is None:
            n_first = n
            break
    dets = np.array(dets)
    norm = dets[n1 - 1] if len(dets) >= n1 else 1.0
    n_interp = None
    if n_first is not None and n_first >= 2:
        d0, d1 = dets[n_first - 2], dets[n_first - 1]
        n_interp = (n_first - 1) + d0 / (d0 - d1)
    return dict(states=states, min_det=dets, min_det_normalised=dets / norm, n_vec=nvecs,
                n_first=n_first, n_interp=n_interp)
