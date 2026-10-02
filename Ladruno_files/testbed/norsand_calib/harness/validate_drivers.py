"""Driver validation (task item 1): every O2 element driver against O1 on one path, first order as in G1
(plan §5.1 "O2 -> O1 with first-order convergence as d_eps -> 0"; G1 test_g1_convergence_tangents.py criteria).

Per driver kind (PS, TC, TE, PSU, TCU): O1 truth (Radau rtol 1e-10, N_TRUTH increments; the endpoint does not
depend on the chunking) and O2 at n = NS increments to the same axial strain. Errors at the endpoint:
  e_sig = ||sigma_O2 - sigma_O1|| / ||sigma_O1||,  e_pi = |pi_O2 - pi_O1| / |pi_O1|,
  e_eps = ||eps_O2 - eps_O1|| / ||eps_O1|| (the mixed-control unknowns, drained kinds).
Gate (G1's): each error strictly decreases along NS and its least-squares log-log order is in [0.8, 1.3].
Extra (undrained kinds, K1.7 closed form, sheet §13): at |eps_a| = 200 % the O2 endpoint has
p -> p_cs = -p_a ((e0_csl - e)/lambda_c)^(1/xi) and zeta(theta) q/|p| -> M.

Parameters: Toyoura PLACEHOLDER CSL (sand.py), BA06 via the energy plug ('per_test'), WW, smooth cap,
theta_V below; start sigma3' 49 kPa, e 0.716 (Tatsuoka Fig. 16a state, data/tatsuoka1986/source/compare.md),
pi_i0 = 0.8 p (rule 'ratio': an elastic start, so O1 has no vertex start); the 'on_surface' (apex) and
'ramp_end' starts are reported separately (census).
Run:  python -m harness.validate_drivers [--quick]     (from norsand_calib/)
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

from .drivers import simulate
from .model import Setup
from .sand import TOYOURA_PLACEHOLDER as TOY

THETA_V = dict(chi=-3.0, h=150.0, N=0.3, N_bar=0.2, rho=0.712, rho_bar=0.75)
SIG3, E_INIT = 49.0, 0.716
# axial strain magnitude per kind. TE stops at 1.5 % (pre-peak; its peak is at ~2.2 % here): run to 5 % its
# endpoint sits on the post-peak asymptote, where the endpoint error decays faster than first order (order 1.64
# for sigma and pi measured 2026-10-02, Esmeralda) and so does not test the step.
EPS_A = dict(PS=0.05, TC=0.05, TE=0.015, PSU=0.05, TCU=0.05)
NS = (25, 50, 100, 200, 400)
N_TRUTH = 20
KINDS = ("PS", "TC", "TE", "PSU", "TCU")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out")


def order(ns, errs):
    return float(np.polyfit(np.log(1.0 / np.asarray(ns, float)), np.log(np.asarray(errs, float)), 1)[0])


def run(quick=False):
    setup = Setup(TOY, pi0_rule="ratio", pi0_ratio=0.8)
    ns = NS[:3] if quick else NS
    report = dict(date=time.strftime("%Y-%m-%d"), theta=THETA_V, sigma3=SIG3, e=E_INIT, eps_a=EPS_A, ns=ns,
                  sand=TOY.name, kinds={})
    all_ok = True
    for kind in KINDS:
        t0 = time.perf_counter()
        tr = simulate(setup, THETA_V, kind, SIG3, E_INIT, EPS_A[kind], N_TRUTH, oracle="O1")
        t_o1 = time.perf_counter() - t0
        rec = dict(o1_status=tr.status, o1_seconds=t_o1)
        if not tr.complete:
            rec["ok"] = False
            all_ok = False
            report["kinds"][kind] = rec
            print(f"{kind}: O1 truth incomplete ({tr.status})")
            continue
        sT, piT, eT = tr.sig[-1], tr.pi_i[-1], tr.eps[-1]
        es, ep, ee, secs, sts = [], [], [], [], []
        for n in ns:
            c = simulate(setup, THETA_V, kind, SIG3, E_INIT, EPS_A[kind], n, oracle="O2")
            sts.append(c.status)
            secs.append(c.stats["seconds"])
            if not c.complete:
                es.append(math.nan); ep.append(math.nan); ee.append(math.nan)
                continue
            es.append(float(np.linalg.norm(c.sig[-1] - sT) / np.linalg.norm(sT)))
            ep.append(float(abs(c.pi_i[-1] - piT) / abs(piT)))
            ee.append(float(np.linalg.norm(c.eps[-1] - eT) / np.linalg.norm(eT)))
        checks = {"e_sig": es, "e_pi": ep} | ({"e_eps": ee} if kind in ("PS", "TC", "TE") else {})
        ok = all(s == "ok" for s in sts)
        orders = {}
        for nm, errs in checks.items():
            if any(not (e > 0.0) for e in errs):
                ok = False
                orders[nm] = None
                continue
            dec = all(errs[i + 1] < errs[i] for i in range(len(errs) - 1))
            p = order(ns, errs)
            orders[nm] = p
            ok = ok and dec and 0.8 <= p <= 1.3
        rec.update(o2_status=sts, o2_seconds=secs, errors=checks, orders=orders, ok=ok,
                   o1_end=dict(sig=sT.tolist(), pi_i=float(piT), eps=eT.tolist()))
        all_ok = all_ok and ok
        report["kinds"][kind] = rec
        print(f"{kind}: O1 {tr.status} ({t_o1:.1f} s); O2 n={list(ns)} status {set(sts)}")
        for nm, errs in checks.items():
            print(f"    {nm}: " + "  ".join(f"{e:.3e}" for e in errs) + f"   order {orders[nm]}")
        print(f"    O2 seconds: " + " ".join(f"{s:.2f}" for s in secs) + f"   -> {'PASS' if ok else 'FAIL'}")
        sys.stdout.flush()

    # K1.7 closed form for the undrained drivers (loose-ish start so p_cs is moderate)
    from o2_algo import kernel as K
    e_cs = 0.90
    for kind in ("PSU", "TCU"):
        c = simulate(setup, THETA_V, kind, SIG3, e_cs, 2.0, 2000, oracle="O2")
        st = c.stats["final_state"]
        p = float(np.trace(st.sigma)) / 3.0
        p_cs = -TOY.p_a * ((TOY.e0 - e_cs) / TOY.lambda_c) ** (1.0 / TOY.xi)
        w = np.linalg.eigvalsh(st.sigma)
        inv = K.invariants(w)
        z, _, _ = K.zeta_y(inv.theta, THETA_V["rho"], "WW")
        eta = z * inv.q / abs(p)
        rec = dict(status=c.status, p_end=p, p_cs=p_cs, rel_p=abs(p - p_cs) / abs(p_cs), eta_zeta=eta,
                   M=TOY.M, rel_eta=abs(eta - TOY.M) / TOY.M, theta_deg=math.degrees(inv.theta),
                   seconds=c.stats["seconds"])
        report["kinds"].setdefault(kind, {})["k17"] = rec
        print(f"{kind} K1.7 at |eps_a| 200 %: p {p:.4f} vs p_cs {p_cs:.4f} (rel {rec['rel_p']:.2e}); "
              f"zeta q/|p| {eta:.5f} vs M {TOY.M} (rel {rec['rel_eta']:.2e}); theta {rec['theta_deg']:.2f} deg; "
              f"{c.status}, {c.stats['seconds']:.1f} s")

    # the pi_i0 rules (model.Setup), reported: which drivers complete from each start, and where O2 substeps
    census = {}
    for rule in ("on_surface", "ramp_end"):
        s2 = Setup(TOY, pi0_rule=rule)
        for kind in KINDS:
            c = simulate(s2, THETA_V, kind, SIG3, E_INIT, EPS_A[kind], 100, oracle="O2")
            census[f"{rule}:{kind}"] = dict(status=c.status, steps_done=len(c.eps) - 1,
                                            substepped_at=c.stats.get("substepped_at"),
                                            bracketed_at=c.stats.get("bracketed_at"))
            print(f"{rule} start {kind}: {c.status} after {len(c.eps) - 1} of 100 increments; "
                  f"O2 substepped at increments {c.stats.get('substepped_at')}")
    report["start_census"] = census
    report["all_ok"] = all_ok
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "validate_drivers.json"), "w") as f:
        json.dump(report, f, indent=1, default=float)
    print("DRIVERS", "PASS" if all_ok else "FAIL")
    return all_ok


if __name__ == "__main__":
    ok = run(quick="--quick" in sys.argv)
    sys.exit(0 if ok else 1)
