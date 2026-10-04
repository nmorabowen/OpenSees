"""Shared O1/O2 interface: State, initial_state, run_path, tangent, drivers.

Tensor convention: 3x3 symmetric numpy arrays, compression negative, TENSOR strains
(no engineering shear anywhere inside; sheet §1.2).

Specific volume (sheet §1.2, G2 owner decision 2026-10-01, exponential update in BOTH modes):
v_{n+1} = v_n exp(tr d_eps) (= v0 exp(tr eps)), with the trial-strain derivative
d v_{n+1} / d eps~_b = v_{n+1} (S.31): the return map is called with vfac = v = v_{n+1}, the
converged specific volume of the step (never v0, never v_n). v0 is carried as a committed datum
only (initial_state; enters no derivative). Supersedes the G0/G1 linear rule v = v0 (1 + tr eps),
vfac = v0.
Small-strain mode (default): eps^e_tr = eps^e_n + d_eps (any symmetric d_eps, non-coaxial
allowed).
Finite-strain mode (State.finite = True, used only by k2_path): the increments are principal
LOG-strain increments of a DIAGONAL, fixed-direction deformation (the (S.43) protocol), for
which b^e,tr = f b^e_n f^T reduces to eps~_a = eps^e_n,a + ln f_a; v = v0 J = v0 exp(tr eps) is
then the same law (§1.4). The return map and its v-handling are identical in the two modes; only
the tangent assembly (S.34) differs.
"""
from __future__ import annotations

import dataclasses
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
    # --- §9.7 floor counters (committed, cumulative; sheet §9.7 "Determinism, counting, interfaces") ---
    eps_f_v: float = 0.0                  # sum of d eps^f_v over the history (init, trial and post events)
    W_f: float = 0.0                      # = p_min eps_f_v, the energy-source bound (S.52)
    n_f_tr: int = 0                       # number of trial floor events
    n_f_post: int = 0                     # number of post floor events
    n_f_init: int = 0                     # 1 if initial_state projected

    def copy(self) -> "State":
        return dataclasses.replace(self, sigma=self.sigma.copy(), eps_e=self.eps_e.copy(), flags=dict(self.flags),
                                   eps_p=self.eps_p.copy(), cache=dict(self.cache))

    @property
    def floor(self) -> dict:
        """The shell's `floor` response: (at_floor, n_f_tr, n_f_post, eps_f_v, W_f) (+ n_f_init)."""
        return dict(at_floor=bool(self.flags.get("at_floor", False)), n_f_tr=self.n_f_tr, n_f_post=self.n_f_post,
                    n_f_init=self.n_f_init, eps_f_v=self.eps_f_v, W_f=self.W_f)


def _eigh(A: np.ndarray):
    w, V = np.linalg.eigh(0.5 * (A + A.T))
    return w, V


def _from_principal(vals: np.ndarray, V: np.ndarray) -> np.ndarray:
    return (V * vals) @ V.T


def initial_state(params: Params, sigma0: np.ndarray, v0: float, pi_i0: float | None,
                  finite: bool = False, pi0_rule: str = "unified") -> State:
    """State at sigma0 with eps^p = 0.

    Floor (sheet §9.7): eps^e := Pi_f(invert(sigma0)) (counted, n_f_init = 1); at a floored point sigma0 is replaced
    by sigma(eps^e_f) (p = -p_min; under HAR q also scaled). sigma0 must have p < 0 (the inverse map needs it; under
    HAR also eps* > 0, which the closed form (S.5h'') guarantees for p < 0).

    pi_i0 None: the unified rule (S.53) (sheet §5.4, owner decision (d)): the surface through
    (p_init, eta*) with eta* = max(eta_init, c2 M), eta_init = zeta(theta) q/|p| of the FLOORED initial stress (0 on
    the axis) and c2 the cap's upper bound (c2 := 0 for cap = 'none', c1 = c2 for 'planar'); refused (ValueError) if
    eta* >= M/N (no surface through the state) or if the §7 guard B <= 0 fails at (pi_i0, psi_i0).
    pi0_rule = 'legacy': eta* = eta_init (the pre-round-3 rule, the apex through p_init for an isotropic start);
    identical to 'unified' when cap = 'none'. An explicit pi_i0 overrides both."""
    params.validate()
    sig = np.asarray(sigma0, float)
    w, V = _eigh(sig)
    eps_e_p = K.invert_elastic(params, w)
    eps_pre = eps_e_p
    fl = K.floor_project(params, eps_e_p)
    n_f_init, eps_f_v = 0, 0.0
    if fl.active:
        eps_e_p = fl.eps_f
        w = K.elastic(params, eps_e_p).sig
        n_f_init, eps_f_v = 1, fl.dfv
    inv = K.invariants(w)
    if pi_i0 is None:
        if pi0_rule not in ("unified", "legacy"):
            raise ValueError(f"pi0_rule must be 'unified' or 'legacy', got {pi0_rule!r}")
        if inv.vertex:
            eta = 0.0
        else:
            z, _, _ = K.zeta_y(inv.theta, params.rho, params.zeta)
            eta = -z * inv.q / inv.p
        if pi0_rule == "unified":
            c2 = 0.0 if params.cap == "none" else params.c2
            eta = max(eta, c2 * params.M)                                     # (S.53)
        pi_i0 = K.pi_of_eta(params, inv.p, eta)                               # raises if eta >= M/N
        if pi0_rule == "unified":
            psi0, _ = K.csl(params, v0, pi_i0)
            try:
                K.pistar(params, inv.p, K.flow(params, inv, pi_i0).Om, psi0)  # the §7 guard B > 0 at (pi_i0, psi_i0)
            except K.EvalError as e:
                raise ValueError(f"initial state refused: {e} at (pi_i0={pi_i0:.6g}, psi_i0={psi0:.4g}) (sheet §5.4, §7)")
    st = State(_from_principal(w, V), _from_principal(eps_e_p, V), float(pi_i0), float(v0),
               v0=float(v0), finite=finite, eps_f_v=eps_f_v, W_f=params.p_min * eps_f_v, n_f_init=n_f_init)
    psi, _ = K.csl(params, v0, pi_i0)
    st.psi_i = psi
    try:
        st.eta = K.eta_of(params, inv.p, pi_i0)
    except ValueError:
        pass
    st.flags["at_floor"] = bool(params.p_min > 0.0 and inv.p > -params.p_min * (1.0 + K.AT_FLOOR_TOL))
    st.flags.update(floor_tr=0, floor_post=0, fpattern="")
    ae = K.elastic(params, eps_e_p).ae
    st.cache = dict(eps_tr=eps_e_p, nvec=V, sig=w, atilde=(ae @ fl.Phi if fl.active else ae),
                    floor_events=([(eps_pre, eps_e_p, fl.dfv)] if fl.active else []))
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

    CONTRACT (tangent, sheet §9.6, owner decision 2026-10-01): on a substepped increment
    (substeps > 1) tangent() returns the CHAINED consistent tangent, the exact derivative of the
    increment's final stress with respect to the TOTAL strain increment, propagated through every
    accepted sub-increment by the recursion (S.46) and assembled by (S.47). It is built in the
    small-strain mode only (finite mode keeps the last sub-increment's (S.34) tangent, see README).
    A refused level discards its sensitivities; the finer level restarts from S_0 = 0. On a
    non-substepped increment (substeps = 1) the tangent is the closed-form CTO (S.31)-(S.33),
    bit-identical to the pre-§9.6 oracle (the chain reduces to it, see step_fractions).
    D is the sum of the sub-increment dissipations (each >= 0), plastic = any sub-increment
    plastic, local_iters / pi_iters are summed, vertex / cap_active are the last sub-increment's."""
    deps = np.asarray(deps, float)
    first = None
    for j in range(K.MAX_SUBSTEP_HALVINGS + 1):
        m = 2 ** j
        cur, ok, C = _run_fractions(params, st, deps, [1.0 / m] * m, chain=(m > 1 and not st.finite))
        if ok:
            cur.flags["substeps"] = m
            if C is not None:
                cur.cache["C_chain"] = C
            return cur
        if first is None:
            first = cur
        last_reason = cur.flags["reason"]
    first.flags["reason"] = f"{last_reason} (substeps exhausted at 2^{K.MAX_SUBSTEP_HALVINGS})"
    first.flags["substeps"] = 2 ** K.MAX_SUBSTEP_HALVINGS
    return first


def step_fractions(params: Params, st: State, deps: np.ndarray, fractions, chain: bool = True) -> State:
    """The increment deps taken as the given sub-increments alpha_k * deps (any alpha_k > 0 with
    sum = 1, e.g. a recursive-halving shape (1/2, 1/4, 1/4)), with NO ladder: a refusal of any
    sub-increment returns that refused state. With chain=True the §9.6 chained tangent is
    assembled for every m, including m = 1 (so the reduction "C_chain(m = 1) = (S.33)" can be
    measured; the ladder itself never chains m = 1). tangent() on the returned state gives the
    chain; cache['atilde'] etc. still describe the last sub-increment."""
    deps = np.asarray(deps, float)
    fr = [float(a) for a in fractions]
    if abs(sum(fr) - 1.0) > 1e-12 or any(a <= 0.0 for a in fr):
        raise ValueError("fractions must be positive and sum to 1")
    cur, ok, C = _run_fractions(params, st, deps, fr, chain=(chain and not st.finite))
    cur.flags["substeps"] = len(fr)
    if ok and C is not None:
        cur.cache["C_chain"] = C
    return cur


def _run_fractions(params: Params, st: State, deps: np.ndarray, fractions, chain: bool):
    """Sub-increments fractions[k] * deps chained on each other; returns (state, ok, C_chain).
    ok = False: `state` is the refused sub-increment's state. The sensitivities (S.46) are
    carried only when chain is True (every accepted sub-increment advances them; an elastic
    sub-increment through the elastic line)."""
    cur = st
    D, iters, piters, plastic, pattern = 0.0, 0, 0, False, ""
    ftr, fpost, fpat, fevents = 0, 0, [], []
    if chain:
        S_eps, S_pi, cum = K.chain_start()
    for a in fractions:
        cur = _step_once(params, cur, a * deps)
        if cur.flags["refused"]:
            return cur, False, None
        D += cur.D
        iters += cur.flags["local_iters"]
        piters += cur.flags["pi_iters"]
        plastic = plastic or cur.flags["plastic"]
        pattern += "P" if cur.flags["plastic"] else "E"      # branch pattern, e.g. "EPPP" (diagnostic)
        ftr += cur.flags["floor_tr"]                         # §9.7: a substepped increment sums its sub-increments
        fpost += cur.flags["floor_post"]
        fpat.append(cur.flags["fpattern"])                   # e.g. "FPf,-P-" (trial floored / plastic / post floored)
        fevents += cur.cache.get("floor_events", [])
        if chain:
            # v_{k+1} = cur.v: the sub-increment's own converged specific volume (S^v closed form, §9.6 G2)
            S_eps, S_pi, cum = K.chain_propagate(S_eps, S_pi, cum, a, cur.v, cur.cache["res"], cur.cache["nvec"])
    cur.D = D
    cur.flags.update(local_iters=iters, pi_iters=piters, plastic=plastic, pattern=pattern,
                     floor_tr=ftr, floor_post=fpost, fpattern=",".join(fpat))
    cur.cache["floor_events"] = fevents
    C = K.chain_assemble(params, cur.cache["res"], cur.cache["nvec"], S_eps) if chain else None
    return cur, True, C


def _step_once(params: Params, st: State, deps: np.ndarray) -> State:
    """One backward-Euler increment without substepping (the floor of §9.7 inside return_map: steps 1f / 5f)."""
    tr = float(np.trace(deps))
    # (S.26) v_{n+1} = v_n exp(tr d_eps), the same law in small-strain and finite mode (§1.2, §1.4; G2
    # owner decision 2026-10-01); d v_{n+1}/d eps~_b = v_{n+1}, so vfac = v (S.31). Was v + v0 tr, vfac = v0.
    v = st.v * math.exp(tr)
    vfac = v
    eps_tr = st.eps_e + 0.5 * (deps + deps.T)
    w, V = _eigh(eps_tr)
    res = K.return_map(params, w, st.pi_i, v, vfac)
    new = State(_from_principal(res.sig, V), _from_principal(res.eps_e, V), res.pi, v,
                v0=st.v0, finite=st.finite, eps_p=st.eps_p.copy(),
                eps_p_v=st.eps_p_v, eps_p_s=st.eps_p_s,
                eps_f_v=st.eps_f_v, W_f=st.W_f, n_f_tr=st.n_f_tr, n_f_post=st.n_f_post, n_f_init=st.n_f_init)
    new.flags = dict(plastic=res.plastic, vertex=res.vertex, cap_active=res.cap_active,
                     refused=res.refused, local_iters=res.local_iters, pi_iters=res.pi_iters,
                     reason=res.reason, substeps=1, floor_tr=0, floor_post=0, fpattern="",
                     at_floor=bool(st.flags.get("at_floor", False)))
    if res.refused:
        # freeze: return the trial-elastic stress with the refusal flag; the driver stops. Counts nothing (§9.7).
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
    # §9.7 counters: per (sub-)increment floor_tr / floor_post in {0, 1}; cumulative n_f, eps^f_v, W_f; at_floor
    new.eps_f_v = st.eps_f_v + res.dfv_tr + res.dfv_post
    new.W_f = params.p_min * new.eps_f_v
    new.n_f_tr = st.n_f_tr + int(res.floor_tr)
    new.n_f_post = st.n_f_post + int(res.floor_post)
    p_c = float(res.sig.sum()) / 3.0
    new.flags.update(floor_tr=int(res.floor_tr), floor_post=int(res.floor_post),
                     fpattern=("F" if res.floor_tr else "-") + ("P" if res.plastic else "E") + ("f" if res.floor_post else "-"),
                     at_floor=bool(params.p_min > 0.0 and p_c > -params.p_min * (1.0 + K.AT_FLOOR_TOL)))
    new.cache = dict(eps_tr=w, nvec=V, sig=res.sig, atilde=res.atilde, res=res, floor_events=list(res.floor_events))
    return new


def floor_energy(params: Params, state: State) -> float:
    """E_f (S.52) of the last increment (sum over its floor events, sub-increments included; the init event for a
    fresh initial_state): the on-demand response (owner decision (c), 2026-10-03), computed from the closed-form
    Psi; the bound p_min d eps^f_v is used for an out-of-domain HAR trial. Always <= the increment's W_f share."""
    return float(sum(K.floor_energy(params, pre, post, dfv) for pre, post, dfv in state.cache.get("floor_events", [])))


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
    """Small-strain consistent tangent of the last increment: the chained tangent (S.47) when
    the increment was substepped (cache['C_chain'], sheet §9.6), else the closed-form CTO (S.33)
    (elastic a^e if the last step was elastic, or for a fresh initial state)."""
    c = state.cache
    if "C_chain" in c:
        return c["C_chain"]
    return K.tangent_small(c["atilde"], c["sig"], c["eps_tr"], c["nvec"])


def tangent_last_substep(params: Params, state: State) -> np.ndarray:
    """The (S.33) CTO of the LAST sub-increment alone (the pre-§9.6 behaviour on substepped
    increments); kept for the measurement of how far it is from the chained tangent."""
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
                if abs(res) <= TRIAX_LAT_TOL_REL * params.p_ref:
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
