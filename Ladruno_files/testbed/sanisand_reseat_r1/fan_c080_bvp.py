"""Route 2 at footing scale: the c = 0.80 control's OWN wall.

C080_EB_off (Esmeralda, build bd93c558d, E_B settings with the Lode parameter c = 0.80 and
R1 OFF) stops at s/B 0.0480 (FLOOR; 82 loadingNonPosH + 287 maxSubsteps refusals).
fan_c080.py measured c-sensitivity on the c = 0.71 wall states; this runs the same wall fan
(32 directions x {3e-6, 3e-5}) on the c = 0.80 run's own committed states at its last
converged step (308), with c = 0.80, for DM04 and R1.  Questions: are its refusers the same
Zeno re-seat set (singular for DM04 in the exact oracle), which Lode side are they on, and
does R1 integrate them?

The states are the loadingNonPosH refusers the run's log names (element, gp):
  the final wall (steps 307-308): 1833/4, 1834/4, 1832/4, 1832/2 (x -0.8..-0.9 m, y -1.5..-2.0 m)
  earlier (steps 248-306):        1817/1, 1817/2, 1973/1, 1973/2, 1829/3 (x +-1.0..1.3 m, y -0.2..-0.3 m)
They are stored in data/refuser_states_c080.csv (compression POSITIVE, as refuser_states.csv),
so the fan reruns without the checkpoint.  To rebuild the CSV, set R1_ESMERALDA_ANALYSIS to a
folder holding ck/C080_EB_off/ckpt/field_last_converged.npz (from
~/ladruno_r1/deck/runs/C080_EB_off/ckpt on Esmeralda) and pass --extract.

Usage: py -3.11 fan_c080_bvp.py [--extract] [variants...]   (default: DM04 B1 T1B1 T1B1S)
Output: out/fan_c080_bvp.json"""
from __future__ import annotations

import collections
import csv
import dataclasses
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import r1common as C

_ARGV = sys.argv[1:]
sys.argv = sys.argv[:1]          # a3_fan parses sys.argv at import time
from a3_fan import deps_of, fib_sphere  # noqa: E402

LEG = "C080_EB_off"
C080 = 0.80
# (element, gp, x, y) as the run's log prints them; k = 4 (element - 1) + (gp - 1)
REFUSERS = [(1833, 4, -0.90, -1.73), (1834, 4, -0.90, -1.54), (1832, 4, -0.90, -1.91),
            (1832, 2, -0.79, -2.02),
            (1817, 1, -1.27, -0.34), (1817, 2, -1.16, -0.34), (1973, 1, +1.16, -0.34),
            (1973, 2, +1.27, -0.34), (1829, 3, -0.98, -0.23)]
CSV_PATH = os.path.join(C.HERE, "data", "refuser_states_c080.csv")
NAMES = ("sigma", "alpha", "alpha_in", "z")
P080 = dataclasses.replace(C.P, c=C080)


def kof(el, gp):
    return 4 * (el - 1) + (gp - 1)


def extract():
    rows = []
    for (el, gp, x, y) in REFUSERS:
        k = kof(el, gp)
        st, info = C.ckpt_state_npz(LEG, k)
        assert abs(info["gx"] - x) < 0.011 and abs(info["gy"] - y) < 0.011, (el, gp, info["gx"], info["gy"])
        r = dict(leg=LEG, k=k, element=el, gp=gp, step=info["step"], s_over_B=info["s_over_B"],
                 gx=info["gx"], gy=info["gy"], p=info["p"], psi=info["psi"], e=info["raw"]["e"])
        for n in NAMES:
            for i, v in enumerate(info["raw"][n]):
                r[f"{n}_{i}"] = repr(float(v))
        rows.append(r)
    os.makedirs(os.path.dirname(CSV_PATH), exist_ok=True)
    with open(CSV_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {CSV_PATH} ({len(rows)} states)")


def load_states():
    out = {}
    for r in csv.DictReader(open(CSV_PATH, newline="")):
        g = lambda n: [float(r[f"{n}_{i}"]) for i in range(6)]
        st = C.State.from_voigt(g("sigma"), g("alpha"), g("z"), float(r["e"]), g("alpha_in"))
        out[int(r["k"])] = (st, dict(element=int(r["element"]), gp=int(r["gp"]), gx=float(r["gx"]),
                                     gy=float(r["gy"]), step=int(r["step"]),
                                     s_over_B=float(r["s_over_B"])))
    return out


def report(st):
    q = C.quantities(st.sigma, st.alpha, st.z, st.e, st.alpha_in, P080, C.UW_MODEL)
    eta = math.sqrt(1.5) * C.norm(q.s) / q.p
    return dict(p=q.p, eta=eta, f=q.f, a_over_cone=q.a / C.CONE, bn=q.bn, rho_alpha=q.rho_alpha,
                rho_b=q.rho_b, psi=q.psi, cos3t=q.cos3t, X=q.X)


_STATES = None


def job(a):
    global _STATES
    if _STATES is None:
        _STATES = load_states()
    k, vname, idir, delta, deps = a
    st, _ = _STATES[k]
    O = C.variants()[vname]
    t0 = time.time()
    try:
        r = C.integrate(st, deps, P080, O, record=False)
        out = dict(status=r.status, reseats=len(r.reseats), max_rho_alpha=r.max_rho_alpha,
                   min_hx=getattr(r, "min_h_over_x", float("nan")))
    except Exception as ex:  # a crash is a finding
        out = dict(status="exception:" + type(ex).__name__, err=str(ex)[:200])
    out.update(k=k, variant=vname, idir=idir, delta=delta, secs=time.time() - t0)
    return out


def trace(k=7331):
    """The re-seat sequence of the first DM04 failure (smallest magnitude) at GP k, as
    zeno_trace.py does at E_B 1880/1: re-seat times, b:n, |dalpha/dt|, cos3theta.
    Needs out/fan_c080_bvp.json.  Output: out/fan_c080_bvp_trace.json"""
    import numpy as np
    from sanisand_r1.integrator import Control, _Increment
    fan = json.load(open(os.path.join(C.OUT, "fan_c080_bvp.json")))["trials"]
    bad = sorted((t for t in fan if t["variant"] == "DM04" and t["k"] == k and t["status"] != "ok"),
                 key=lambda t: (t["delta"], t["idir"]))[0]
    st, info = load_states()[k]
    de = deps_of(fib_sphere(32)[bad["idir"]], bad["delta"])
    out = []
    for vname in ("DM04", "B1", "T1B1S"):
        O = C.variants()[vname]
        r = C.integrate(st, de, P080, O, record=True)
        T = np.array(r.path["t"]); Y = np.array(r.path["y"])
        inc = _Increment(st, Control.strain(de), P080, O, 1e-10, 1e-2, "Radau")
        rows = []
        for j0, rs in enumerate(r.reseats):
            j = int(np.argmin(np.abs(T - rs["t"])))
            y = Y[j]
            q = C.quantities(C.v2t(y[0:6]), C.v2t(y[6:12]), C.v2t(y[12:18]), float(y[18]),
                             C.v2t(y[6:12]), P080, O)      # just re-seated: alpha_in = alpha
            inc.alpha_in = C.v2t(y[6:12]).copy()
            dy, _ = inc.rates(y, "plastic", q)
            rows.append(dict(k=j0, t=float(rs["t"]), bn=float(q.bn),
                             dalpha=float(C.norm(C.v2t(dy[6:12]))), cos3t=float(q.cos3t)))
        out.append(dict(variant=vname, status=r.status, t_end=float(r.t_end),
                        n_reseats=len(r.reseats), reseats=rows))
        print(f"{vname:>6} {r.status:>14} t_end {r.t_end:.6f} re-seats {len(r.reseats)}")
        for w in rows[:10]:
            print(f"        t {w['t']:.9f}  b:n {w['bn']:.2e}  |dalpha/dt| {w['dalpha']:.2e}  cos3t {w['cos3t']:.3f}")
    json.dump(dict(state=info, idir=bad["idir"], delta=bad["delta"], runs=out),
              open(os.path.join(C.OUT, "fan_c080_bvp_trace.json"), "w"), indent=1)


def main(argv):
    if "--trace" in argv:
        trace()
        return
    if "--extract" in argv:
        extract()
        argv = [a for a in argv if a != "--extract"]
    states = load_states()
    V = argv or ["DM04", "B1", "T1B1", "T1B1S"]
    print(f"{LEG}, c = {C080}: committed states at step {next(iter(states.values()))[1]['step']}, "
          f"s/B {next(iter(states.values()))[1]['s_over_B']:.5f}")
    print(f"{'ele/gp':>8} {'x':>6} {'y':>6} {'p':>7} {'eta':>6} {'a/cone':>8} {'b:n':>9} "
          f"{'rho_a':>6} {'cos3t':>7} {'X':>9}")
    diag = {}
    for k, (st, info) in states.items():
        d = report(st)
        diag[k] = dict(info, **d)
        print(f"{info['element']:>5}/{info['gp']} {info['gx']:>6.2f} {info['gy']:>6.2f} {d['p']:>7.2f} "
              f"{d['eta']:>6.3f} {d['a_over_cone']:>8.3f} {d['bn']:>9.2e} {d['rho_alpha']:>6.3f} "
              f"{d['cos3t']:>7.3f} {d['X']:>9.2e}")
    dirs = fib_sphere(32)
    tasks = [(k, v, i, d, deps_of(u, d)) for k in states for v in V
             for d in (3e-6, 3e-5) for i, u in enumerate(dirs)]
    t0 = time.time()
    with Pool(6, maxtasksperchild=50) as pool:
        res = pool.map(job, tasks, chunksize=4)
    json.dump(dict(leg=LEG, c=C080, states={str(k): v for k, v in diag.items()}, trials=res),
              open(os.path.join(C.OUT, "fan_c080_bvp.json"), "w"))
    print(f"{len(res)} trials in {time.time() - t0:.0f} s")
    for v in V:
        rr = [r for r in res if r["variant"] == v]
        print(f"{v:>6}: failed {sum(1 for r in rr if r['status'] != 'ok')}/{len(rr)}",
              dict(collections.Counter(r["status"] for r in rr)))
        for k, (st, info) in states.items():
            s = [r for r in rr if r["k"] == k]
            print(f"        {info['element']}/{info['gp']}: failed "
                  f"{sum(1 for r in s if r['status'] != 'ok')}/{len(s)}, "
                  f"max re-seats {max(r.get('reseats', 0) for r in s)}")


if __name__ == "__main__":
    main(_ARGV)
