"""WP-230 (TIMs F25): the low-pressure CLAMP as a state machine, a DEFINITION, driven around the
unmodified WP-134 oracle (the WP-151 testbed copy `sanisand_r1`) -- the twin of WP-152's separation
driver (`../sanisand_tension_cutoff/tc_oracle.py`), with ONE thing changed: a point that loses its
confinement KEEPS ITS DEVIATOR at p = p_floor instead of going to the isotropic state at p_min.

State machine (Ladruno_implementation/230_sanisand_lowp_clamp.md).  Triggers, gates, the closing law
and the re-contact condition are WP-152's, unchanged:

  NORMAL   integrate the increment.  status p_floor (the exact trajectory reaches p -> 0 inside the
           increment: SAS-ME's code 6, trigger E1) at a committed model p0 <= p0_max under a
           NON-COMPRESSING increment -> CLAMPED at the END of the increment:
               sigma = p_floor I + s_k,   s_k = dev(sigma_committed)  (the last committed deviator),
               fabric kept, e follows the strain, tr(eps_entry) := tr(eps).
           The gates hold exactly as in WP-152 (refused_highp, refused_compressing).
  CLAMPED  g = tr(eps) - tr(eps_entry) (compression positive).
               sigma = p_c(g) I + s_k,   p_c(g) = p_floor + K(p_contact) max(g, 0).
           The deviator is FROZEN: it does not follow the strain (a clamped point is not integrated;
           it holds the shear it had when it lost its confinement), so it is bounded by what the
           committed state carried.  g >= g_c = (p_contact - p_floor) / K(p_contact) -> NORMAL at the
           end of this increment with sigma = p_re I + s_k, p_re = p_c(g) >= p_contact.

alpha and alpha_in of a clamped point (`alpha_rule`).  Neither adds a stress to a stress ratio:
  "leave"    (the default) alpha and alpha_in stay as committed while clamped: no plastic flow
             happens, so no internal variable moves.  The clamped state is then OUTSIDE the yield
             cone, and f / (sqrt(2/3) m p) -- the distance in cone radii -- is recorded per step.
             At RE-CONTACT the integrator needs a start with f <= 0, so alpha is moved there by the
             closest-point projection below, if (and only if) the exit state is outside the cone.
  "project"  at entry and at every clamped step alpha is moved to the closest point of the yield
             cone around the clamped stress:
                 x = s_k - p alpha,   n = x / |x|,   alpha := s_k / p - sqrt(2/3) m n     (f = 0)
             and alpha_in := alpha.  The stress is then ON the cone throughout, and the price moves
             to rho_alpha = |alpha| / alpha^b, which is recorded.
Under both rules alpha_in := alpha at re-contact (a new loading process starts there; WP-152 does
the same with alpha = alpha_in = 0).

`mode="separation"` reproduces tc_oracle.drive() record for record (checked by `--selfcheck`).
Strain increments are compression positive, Voigt engineering shear (the oracle's)."""
from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "sanisand_reseat_r1"))
sys.path.insert(0, os.path.join(HERE, "..", "sanisand_tension_cutoff"))
import r1common as C  # noqa: E402
from sanisand_r1 import State, integrate  # noqa: E402
from sanisand_r1.model import (SQ23, dev, elastic_moduli, norm, rho_alpha_of, e_critical,  # noqa: E402
                               g_lode, t2v)

I3 = np.eye(3)
Z3 = np.zeros((3, 3))


def _e_after(e, dv, P, O):
    ref = P.e_init if O.void_ratio_law == "initial" else e
    return e - (1.0 + ref) * dv


def model_p(st, O):
    return float(np.trace(st.sigma)) / 3.0 + O.p_residual


def f_ratio(sig, alpha, P):
    """f / (sqrt(2/3) m p): the distance of the stress from the yield cone, in cone radii
    (<= 0 inside or on the cone)."""
    p = float(np.trace(sig)) / 3.0
    return (norm(dev(sig) - p * alpha) - SQ23 * P.m * p) / (SQ23 * P.m * p)


def bound_ratio(sig, e, P, O):
    """Stress ratio over the BOUNDING stress ratio at the stress's own Lode angle:
    sqrt(3/2)|s|/p / (g(theta_s) Mc exp(-nb psi)).  Above 1 the stress is outside the bounding
    surface whatever alpha is: the part of the clamp's price no treatment of alpha removes."""
    p = float(np.trace(sig)) / 3.0
    s = dev(sig)
    ns = norm(s)
    if ns <= 1.0e-14 * p:
        return 0.0
    u = s / ns
    c3 = max(-1.0, min(1.0, math.sqrt(6.0) * float(np.trace(u @ u @ u))))
    psi = e - e_critical(p + O.p_residual, P)
    return math.sqrt(1.5) * (ns / p) / (g_lode(c3, P.c) * P.Mc * math.exp(-P.nb * psi))


def project_alpha(sig, alpha, P):
    """Closest point of the yield cone around `sig` to `alpha` (unchanged if f <= 0)."""
    p = float(np.trace(sig)) / 3.0
    s = dev(sig)
    x = s - p * alpha
    nx = norm(x)
    if nx <= SQ23 * P.m * p:
        return alpha, False
    return s / p - SQ23 * P.m * (x / nx), True


def drive(state0, deps_list, O, p_sep, p_contact, P=None, p0_max=None, mode="clamp", p_floor=0.1,
          alpha_rule="leave"):
    """One record per increment: mode after it ('N' normal / 'S' separated / 'C' clamped), event,
    sigma (Voigt, compression positive), alpha, alpha_in, model p, g, f_ratio, rho_alpha, status."""
    assert mode in ("clamp", "separation") and alpha_rule in ("leave", "project")
    P = C.P if P is None else P
    p0_max = p_contact if p0_max is None else p0_max
    p_base = p_floor if mode == "clamp" else O.p_min
    st, md, g, s_k = state0, "N", 0.0, None
    out = []
    for k, de in enumerate(deps_list):
        de = [float(x) for x in de]
        dv = de[0] + de[1] + de[2]
        ev, status = None, None
        if md == "N":
            r = integrate(st, de, P, O, record=False)
            status = r.status
            if status == "ok":
                st = r.state
            elif status == "p_floor" and model_p(st, O) > p0_max:
                out.append(dict(k=k, mode=md, event="refused_highp", status=status, p0=model_p(st, O)))
                break
            elif status == "p_floor" and dv > 1.0e-10 * float(np.linalg.norm(de)):
                out.append(dict(k=k, mode=md, event="refused_compressing", status=status,
                                p0=model_p(st, O)))
                break
            elif status == "p_floor":                                   # E1 (tension)
                e_end = _e_after(st.e, dv, P, O)
                if mode == "separation":
                    st = State((p_base - O.p_residual) * I3, Z3.copy(), st.z.copy(), e_end, Z3.copy())
                    md, g, ev = "S", 0.0, "enter_tension"
                else:
                    s_k = dev(st.sigma)
                    sig = (p_base - O.p_residual) * I3 + s_k
                    al, ain = st.alpha.copy(), st.alpha_in.copy()
                    if alpha_rule == "project":
                        al, _ = project_alpha(sig, al, P)
                        ain = al.copy()
                    st = State(sig, al, st.z.copy(), e_end, ain)
                    md, g, ev = "C", 0.0, "enter_tension"
            else:
                out.append(dict(k=k, mode=md, event="refused", status=status))
                break
        else:
            e_end = _e_after(st.e, dv, P, O)
            g += dv
            _, K = elastic_moduli(p_contact, e_end, P, O)
            gc = (p_contact - p_base) / K
            p_c = p_base + K * max(g, 0.0)                              # continuous closing
            if md == "S":
                st = State((p_c - O.p_residual) * I3, Z3.copy(), st.z.copy(), e_end, Z3.copy())
            else:
                sig = (p_c - O.p_residual) * I3 + s_k
                al, ain = st.alpha, st.alpha_in
                if alpha_rule == "project":
                    al, _ = project_alpha(sig, al, P)
                    ain = al.copy()
                st = State(sig, al.copy(), st.z.copy(), e_end, ain.copy())
            if g >= gc:
                if md == "C":                                           # hand the integrator f <= 0
                    al, _ = project_alpha(st.sigma, st.alpha, P)
                    st = State(st.sigma, al, st.z, st.e, al.copy())
                md, ev = "N", "recontact"
        psi = st.e - e_critical(model_p(st, O), P)
        out.append(dict(k=k, mode=md, event=ev, status=status, sigma=t2v(st.sigma).tolist(),
                        alpha=t2v(st.alpha).tolist(), alpha_in=t2v(st.alpha_in).tolist(),
                        p=model_p(st, O), g=g, e=st.e, f_ratio=f_ratio(st.sigma, st.alpha, P),
                        rho_alpha=rho_alpha_of(st.alpha, psi, P),
                        bound=bound_ratio(st.sigma, st.e, P, O),
                        q=math.sqrt(1.5) * norm(dev(st.sigma))))
    return out


def work(records, deps_list, sigma0):
    """Net work sum sigma_mid : deps (trapezoid; Voigt engineering shear, so shear counts once)."""
    s_prev = np.array(sigma0, dtype=float)
    W = 0.0
    for r, de in zip(records, deps_list):
        if "sigma" not in r:
            break
        s = np.array(r["sigma"])
        W += float(np.dot(0.5 * (s + s_prev), np.array(de, dtype=float)))
        s_prev = s
    return W


# ---------------------------------------------------------------- the element paths
# WP-152's three (isotropic, uniaxial extension, isotropic cycles) carry no shear at the entry, or
# little; two more DO, since the kept deviator is what this mode is about.  Normal strains only, so
# tests/wp152_cutoff_tools.py's single-brick driver can run them.
def _dev(d):
    return [d, -0.5 * d, -0.5 * d, 0, 0, 0]          # isochoric triaxial shear, x the axis


def path_dev_open(n_sh=30, n_out=80, n_hold=20, n_back=120, d=1.0e-5):
    """Triaxial shear to a mobilised stress ratio, isotropic opening past the entry, more shear
    WHILE clamped (the frozen deviator must not move), then isotropic closing past re-contact."""
    return [_dev(d)] * n_sh + [[-d, -d, -d, 0, 0, 0]] * n_out + [_dev(d)] * n_hold + [[d, d, d, 0, 0, 0]] * n_back


def path_te_dev(n_out=60, n_sh=20, n_back=120, d=1.0e-5):
    """Uniaxial extension past the entry, triaxial shear while clamped, uniaxial recompression."""
    return [[-d, 0, 0, 0, 0, 0]] * n_out + [_dev(d)] * n_sh + [[d, 0, 0, 0, 0, 0]] * n_back


def paths():
    import tc_oracle as T
    return dict(iso=T.path_iso(), te=T.path_te(), cyc=T.path_cycles(),
                dev_open=path_dev_open(), te_dev=path_te_dev())


def flip_state(worktree):
    import json
    flip = json.load(open(os.path.join(worktree, "tests", "data", "wp152_flip_state.json")))
    s = flip["sigma"]
    sig = np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])
    return flip, State(sig, Z3.copy(), Z3.copy(), flip["e"], Z3.copy())


def fixture(worktree):
    """tests/data/wp230_oracle_paths.json: the five element paths from WP-152's post-flip state, R1
    (the full set), p_sep 0.5, p_contact 1.0, p_floor 0.1, under "clamp/leave" and "clamp/project"."""
    import json
    flip, st0 = flip_state(worktree)
    O = C.variants()["T1B1S"]
    p_sep, p_contact, p_floor = 0.5, 1.0, 0.1
    out = dict(note="WP-230 oracle expectations: clamp_oracle.py --fixture (unmodified WP-134 oracle, R1 "
                    "T1B1S, p_sep 0.5 p_contact 1.0 p_floor 0.1, p_min 0.0101, p_r 0)",
               flip=flip, p_sep=p_sep, p_contact=p_contact, p_floor=p_floor, paths={})
    for name, path in paths().items():
        for rule in ("leave", "project"):
            rec = drive(st0, path, O, p_sep, p_contact, mode="clamp", p_floor=p_floor, alpha_rule=rule)
            cl = [r for r in rec if r.get("mode") == "C"]
            out["paths"][f"{name}|{rule}"] = dict(
                increments=path, records=rec, events=[(r["k"], r["event"]) for r in rec if r.get("event")],
                work=work(rec, path, t2v(st0.sigma)),
                max_f_ratio=max([r["f_ratio"] for r in cl], default=0.0),
                max_rho_alpha=max([r["rho_alpha"] for r in cl], default=0.0),
                max_bound=max([r["bound"] for r in cl], default=0.0))
            o = out["paths"][f"{name}|{rule}"]
            print(f"{name:10s} {rule:8s} events {o['events']} last {rec[-1].get('mode')} "
                  f"work {o['work']:.4g} max f/(cone) {o['max_f_ratio']:.4g} max rho_alpha {o['max_rho_alpha']:.4g} "
                  f"max bound {o['max_bound']:.4g}")
    json.dump(out, open(os.path.join(worktree, "tests", "data", "wp230_oracle_paths.json"), "w"))


def selfcheck(worktree):
    """mode='separation' must be tc_oracle.drive(), record for record."""
    import tc_oracle as T
    _, st0 = flip_state(worktree)
    O = C.variants()["T1B1S"]
    for name, path in (("iso", T.path_iso()), ("te", T.path_te()), ("cyc", T.path_cycles())):
        a = T.drive(st0, path, O, 0.5, 1.0)
        b = drive(st0, path, O, 0.5, 1.0, mode="separation")
        assert len(a) == len(b), name
        for x, y in zip(a, b):
            for key in x:
                assert x[key] == y[key], (name, x["k"], key, x[key], y[key])
        print(f"separation mode == tc_oracle.drive on {name}: {len(a)} records identical")


if __name__ == "__main__":
    wt = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
    if len(sys.argv) > 1 and sys.argv[1] == "--fixture":
        fixture(wt)
    elif len(sys.argv) > 1 and sys.argv[1] == "--selfcheck":
        selfcheck(wt)
    else:
        selfcheck(wt)
        fixture(wt)
