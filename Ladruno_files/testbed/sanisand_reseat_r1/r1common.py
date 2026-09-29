"""R1 prototype: shared setup (paths, variants, state loaders).

Everything runs on the scratchpad COPY `sanisand_r1` of the WP-134 oracle
(Ladruno_scripts/sanisand_reference @ 6cef73cc8), whose R1 toggles are all OFF by
default (bit-identical to the oracle: check_identity.py).
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT3 = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
W = os.environ.get("R1_WORKTREE", _ROOT3 if os.path.basename(os.path.dirname(HERE)) == "testbed"
                   else r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees\sharp-chandrasekhar-d6ff72")
ESM = os.environ.get("R1_ESMERALDA_ANALYSIS",
                     r"C:\Users\nmora\AppData\Local\Temp\claude\C--Users-nmora-Github-OpenSees-Compile-OpenSees--claude-worktrees-tims-implementation-review-3733c6\e034e494-142a-4e05-862e-9d267344b7ad\scratchpad\esmeralda_wp138\analysis")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, HERE)

from sanisand_r1 import CAMPAIGN, Options, State, integrate, quantities  # noqa: E402
from sanisand_r1 import ring as rring  # noqa: E402
from sanisand_r1.model import SQ23, v2t, t2v, dev, norm, ddot, I3  # noqa: E402

P = CAMPAIGN
UW_MODEL = rring.ring_variants()["uw_model"]      # the WP-134 oracle ("DM04" below)
CONE = SQ23 * P.m                                  # yield-cone radius in alpha space, 0.00408


def variants(cone=None):
    """name -> Options.  DM04 = the oracle.  Naming: B = bounded h (denominator
    floor eps, in cone radii; `add` = the PM4Sand-like additive form), T = threshold
    re-seat (delta, in cone radii), S = softening cap (kappa).  `cone` = the yield
    cone radius sqrt(2/3) m of the parameter set in use (default: the campaign's)."""
    b = UW_MODEL
    c = CONE if cone is None else cone
    v = {
        "DM04": b,
        # bounded h only, DM04 exactly for a >= eps
        "B1": b.with_(h_reg="max", h_eps=c),
        "B0.25": b.with_(h_reg="max", h_eps=0.25 * c),
        # PM4Sand-like additive floor (changes h everywhere by eps/a)
        "Badd1": b.with_(h_reg="add", h_eps=c),
        # bounded h + softening cap (H >= kappa X)
        "B1S": b.with_(h_reg="max", h_eps=c, h_soft_kappa=0.5),
        # threshold re-seat + bounded h
        "T0.5B1": b.with_(h_reg="max", h_eps=c, reseat_delta=0.5 * c),
        "T1B1": b.with_(h_reg="max", h_eps=c, reseat_delta=c),
        "T2B1": b.with_(h_reg="max", h_eps=c, reseat_delta=2 * c),
        "T2B0.25": b.with_(h_reg="max", h_eps=0.25 * c, reseat_delta=2 * c),
        # all three
        "T2B1S": b.with_(h_reg="max", h_eps=c, h_soft_kappa=0.5, reseat_delta=2 * c),
        "T1B1S": b.with_(h_reg="max", h_eps=c, h_soft_kappa=0.5, reseat_delta=c),
        # WP-150 memo (draft #892) R1 exactly: floor only where b:n <= 0 (c_A = 1),
        # and R1 + R1b (hysteretic re-seat a_rev = one cone radius)
        "R150": b.with_(h_reg="max_soft", h_eps=c),
        "R150+T1": b.with_(h_reg="max_soft", h_eps=c, reseat_delta=c),
        # c_A sensitivity at c_rev = 1 (with and without the cap)
        "T1B0.5": b.with_(h_reg="max", h_eps=0.5 * c, reseat_delta=c),
        "T1B2": b.with_(h_reg="max", h_eps=2 * c, reseat_delta=c),
        "T1B0.5S": b.with_(h_reg="max", h_eps=0.5 * c, h_soft_kappa=0.5, reseat_delta=c),
        "T1B2S": b.with_(h_reg="max", h_eps=2 * c, h_soft_kappa=0.5, reseat_delta=c),
    }
    return v


# ---------------------------------------------------------------------------
# states
# ---------------------------------------------------------------------------
CSV8 = os.path.join(W, "Ladruno_implementation", "_tims_2d_model_requests_2026-09-25",
                    "ring_points_b8.csv")


def ring_rows(mesh="b8"):
    return rrows(os.path.join(os.path.dirname(CSV8), f"ring_points_{mesh}.csv"))


def rrows(path):
    return rring.load_ring_csv(path)


def ring_state(mesh, element, gp):
    for r in ring_rows(mesh):
        if r["element"] == element and r["gp"] == gp:
            return rring.row_state(r), r
    raise KeyError((mesh, element, gp))


def reproducer():
    return rring.reproducer_state(), list(rring.REPRODUCER_DEPS)


def _p_from_psi(e, psi):
    ec = e - psi
    arg = (P.e0 - ec) / P.lambda_c
    return P.P_atm * max(arg, 0.0) ** (1.0 / P.xi)


REFUSER_CSV = os.path.join(HERE, "data", "refuser_states.csv")
_REF_CACHE = {}


def _refuser_rows():
    if not _REF_CACHE and os.path.exists(REFUSER_CSV):
        import csv as _csv
        for r in _csv.DictReader(open(REFUSER_CSV, newline="")):
            _REF_CACHE[(r["leg"], int(r["k"]))] = r
    return _REF_CACHE


def ckpt_state(leg, k, which="field_last_converged.npz"):
    """Committed wall state of GP index k: from data/refuser_states.csv when it holds
    it (the committed testbed), else from the Esmeralda checkpoint (ckpt_state_npz)."""
    rows = _refuser_rows()
    if which == "field_last_converged.npz" and (leg, k) in rows:
        r = rows[(leg, k)]
        g = lambda n: [float(r[f"{n}_{i}"]) for i in range(6)]
        st = State.from_voigt(g("sigma"), g("alpha"), g("z"), float(r["e"]), g("alpha_in"))
        info = dict(leg=leg, k=k, gx=float(r["gx"]), gy=float(r["gy"]), step=int(r["step"]),
                    s_over_B=float(r["s_over_B"]), p=float(r["p"]), psi=float(r["psi"]),
                    f_read=float("nan"))
        return st, info
    return ckpt_state_npz(leg, k, which)


def ckpt_state_npz(leg, k, which="field_last_converged.npz"):
    """Committed state of GP index k from the Esmeralda field checkpoint (the
    WP-138 deck's save_field: sig = OpenSees plane-strain stress [xx yy xy],
    compression NEGATIVE; st = getState() = [eps_e(6) alpha(6) z(6) alpha_in(6)
    e dgamma]; sigma_zz recovered from psi and e with p_r = 0, as the deck does).
    Returns (State, info dict incl. the last converged strain increment)."""
    z = np.load(os.path.join(ESM, "ck", leg, "ckpt", which))
    sig, st, psi = z["sig"][k], z["st"][k], float(z["psi"][k])
    e = float(st[24])
    p = _p_from_psi(e, psi)
    s6 = np.zeros(6)
    s6[0], s6[1], s6[3] = -sig[0], -sig[1], -sig[2]
    s6[2] = 3.0 * p - s6[0] - s6[1]
    state = State.from_voigt(s6, st[6:12], st[12:18], e, st[18:24])
    info = dict(leg=leg, k=k, tag=float(z["tag"][k]), gx=float(z["gx"][k]),
                gy=float(z["gy"][k]), step=int(z["step"]), s_over_B=float(z["s_over_B"]),
                p=p, psi=psi, f_read=float(z["f"][k]),
                raw=dict(sigma=s6.tolist(), alpha=st[6:12].tolist(), z=st[12:18].tolist(),
                         alpha_in=st[18:24].tolist(), e=e))
    return state, info


def last_increment(leg, k, prev, last="field_last_converged.npz"):
    """The last converged strain increment at GP k (compression positive, Voigt
    engineering shear, plane strain): -(eps_last - eps_prev).  From
    data/refuser_states.csv when it holds it."""
    rows = _refuser_rows()
    if last == "field_last_converged.npz" and (leg, k) in rows and rows[(leg, k)].get("deps_0"):
        return [float(rows[(leg, k)][f"deps_{i}"]) for i in range(6)]
    a = np.load(os.path.join(ESM, "ck", leg, "ckpt", prev))["eps"][k]
    b = np.load(os.path.join(ESM, "ck", leg, "ckpt", last))["eps"][k]
    d = -(b - a)
    return [d[0], d[1], 0.0, d[2], 0.0, 0.0]


# the loadingNonPosH refusers at the walls (orchestrator's analysis/tables/floor_refusers.csv)
REFUSERS = [
    # leg, k, element, gp, previous checkpoint for the loading direction
    ("E_B", 7516, 1880, 1, "field_step00375.npz"),
    ("E_B", 7512, 1879, 1, "field_step00375.npz"),
    ("E_D", 7844, 1962, 1, None),
    ("E_D", 8228, 2058, 1, None),
    ("E_B16", 31279, 7820, 4, None),
]


def state_report(st, O=UW_MODEL):
    q = quantities(st.sigma, st.alpha, st.z, st.e, st.alpha_in, P, O)
    eta = math.sqrt(1.5) * norm(q.s) / q.p
    return dict(p=q.p, eta=eta, f=q.f, a=q.a, bn=q.bn, rho_alpha=q.rho_alpha,
                rho_b=q.rho_b, psi=q.psi, cos3t=q.cos3t, X=q.X, b0=q.b0,
                a_over_cone=q.a / CONE, alpha_norm=norm(st.alpha),
                ain_norm=norm(st.alpha_in))


def fmt(x, n=3):
    if isinstance(x, float):
        if x != x:
            return "nan"
        if abs(x) >= 1e4 or (abs(x) < 1e-3 and x != 0.0):
            return f"{x:.{n-1}e}"
        return f"{x:.{n}g}"
    return str(x)
