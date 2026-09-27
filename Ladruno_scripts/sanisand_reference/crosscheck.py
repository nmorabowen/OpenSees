"""Reference vs the fork's C++ (ladrunoSANISANDReplay) on benign states, with
every difference attributed to a named UW addition by toggling it alone.

Benign states: ON the yield surface, alpha well inside the bounding surface,
alpha_in = 0 (so (alpha - alpha_in):n > 0), z = 0, at p0 = 20 / 50 / 100 kPa,
three stress directions.  Probes: one strain increment of size delta in four
directions (triaxial loading, triaxial unloading, simple shear, plane-strain
isotropic compression), delta = 1e-5 / 1e-4 / 1e-3.
"""
from __future__ import annotations

import math
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from .integrator import Control, integrate
from .model import CAMPAIGN, Options, State, on_yield_state, t2v, v2t

UW_TOGGLES = {
    # name: the Options field and its PAPER value (switching it back to paper
    # inside the all-UW preset isolates that one addition's effect)
    "U1 d_factor": ("d_factor", False),
    "U3 p_min": ("p_min", 0.0),
    "U4 G(e_init)": ("g_void_ratio", "current"),
    "U5 e-law(e_init)": ("void_ratio_law", "current"),
    "U6 alpha_in rule": ("alpha_in_rule", "paper"),
    "U7 h cap": ("h_cap", None),
    "U8+U9 frozen K,G": ("elastic_moduli", "continuous"),
    "U9 frozen K,G in plastic part": ("elastic_moduli", "frozen"),
}

DIRECTIONS = {
    "TC": np.diag([2.0, -1.0, -1.0]),
    "TCshear": np.array([[1.0, 0.6, 0.0], [0.6, -0.5, 0.0], [0.0, 0.0, -0.5]]),
    "TE": np.diag([-2.0, 1.0, 1.0]),
}
ETAS = {"TC": 0.5, "TCshear": 0.8, "TE": 0.4}   # times Mc (TE: times c*Mc)


def benign_states(P=CAMPAIGN, p0s=(20.0, 50.0, 100.0), e=0.72):
    out = []
    for p0 in p0s:
        for name, dirn in DIRECTIONS.items():
            eta = ETAS[name] * P.Mc * (P.c if name == "TE" else 1.0)
            out.append((f"p{int(p0)}_{name}", on_yield_state(p0, eta, e, P, dirn)))
    return out


def probes(delta):
    return {
        "txLoad": [delta, 0.0, 0.0, 0.0, 0.0, 0.0],
        "txUnload": [-delta, 0.0, 0.0, 0.0, 0.0, 0.0],
        "shear": [0.0, 0.0, 0.0, delta, 0.0, 0.0],
        "isoComp": [delta, delta, 0.0, 0.0, 0.0, 0.0],
    }


def job_of(state, deps, proto):
    return dict(proto=proto, sigma=t2v(state.sigma).tolist(),
                alpha=t2v(state.alpha).tolist(), alpha_in=t2v(state.alpha_in).tolist(),
                z=t2v(state.z).tolist(), e=state.e, deps=list(deps))


def _ref_one(args):
    state, deps, P, O, rtol = args
    r = integrate(state, Control.strain(deps), P, O, rtol=rtol, record=False)
    return dict(status=r.status, sigma=t2v(r.state.sigma).tolist(),
                alpha=t2v(r.state.alpha).tolist(), z=t2v(r.state.z).tolist(),
                e=r.state.e, f_end=r.f_end, max_rho_b=r.max_rho_b,
                segments=[(s["mode"], s["event"]) for s in r.segments],
                reseats=len(r.reseats), notes=r.notes)


def run_reference(cases, P, variants, rtol=1.0e-10, workers=None):
    """cases: list of (state, deps). variants: dict name -> Options.
    Returns {variant: [result per case]}."""
    tasks, keys = [], []
    for vn, O in variants.items():
        for i, (st, de) in enumerate(cases):
            tasks.append((st, de, P, O, rtol))
            keys.append((vn, i))
    out = {vn: [None] * len(cases) for vn in variants}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for (vn, i), res in zip(keys, ex.map(_ref_one, tasks, chunksize=2)):
            out[vn][i] = res
    return out


def rel_incr_diff(start, a, b, key):
    """|| (a - start) - (b - start) || / || b - start ||  on the tensor `key`
    (Frobenius, tensor components), with the absolute difference too."""
    s = v2t(start[key]) if isinstance(start[key], list) else start[key]
    da = v2t(a[key]) - s
    db = v2t(b[key]) - s
    num = float(np.linalg.norm(da - db))
    den = float(np.linalg.norm(db))
    return (num / den if den > 0 else float("nan")), num


def uw_variants(pmin=0.0101):
    """paper; uw (the RK45 comparator); uw_me (the ModifiedEuler comparator);
    and uw_me with ONE addition switched back to the paper ("uw_me-<name>")."""
    uw = Options.uw(p_min=pmin)
    uw_me = Options.uw_me(p_min=pmin)
    v = {"paper": Options(), "uw": uw, "uw_me": uw_me}
    for name, (fld, paper_val) in UW_TOGGLES.items():
        v["uw_me-" + name] = uw_me.with_(**{fld: paper_val})
    return v
