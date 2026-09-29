"""Why floor + hysteresis WITHOUT the cap still walls at footing scale.

R1_EB_fh (Esmeralda, build bd93c558d, E_B settings, -sasHFloor 1 -sasReseatHyst 1,
no cap) stops at s/B 0.0525 with 93 loadingNonPosH refusals at 11 Gauss points; fhc
(the same plus -sasSoftCap 0.5) runs on past s/B 0.054 with none.  On those 11 points'
committed states at fh's last converged step (317), c = 0.71 (the campaign set):
  - the state: p, eta, cos3theta(n), (alpha-alpha_in):n / rho_c, b:n, rho_alpha, and H/X
    with the FLOORED h (1 + Kp/X, Kp = 2/3 p b0/max(a, c_A rho_c) b:n): < 0 means the
    floor alone cannot keep H > 0 there -- a genuine b:n < 0 (post-peak) state, not a
    leftover Zeno sequence (which needs b:n -> 0+);
  - the wall fan (32 directions x {3e-6, 3e-5}) for DM04, fh (T1B1), fhc (T1B1S, kappa 0.5)
    and fhc with kappa 0.25 / 0.75 (sensitivity of the end stress to kappa).
States: data/refuser_states_fh.csv (compression POSITIVE); rebuild with --extract and
R1_ESMERALDA_ANALYSIS -> a folder with ck/R1_EB_fh/ckpt/field_last_converged.npz.
Output: out/fan_fh_wall.json"""
from __future__ import annotations

import collections
import csv
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

LEG = "R1_EB_fh"
# (element, gp, x, y, first refusal step) from the run's log
REFUSERS = [(1879, 4, -0.15, -2.10, 262), (1882, 1, -0.15, -1.65, 262), (1871, 2, -0.23, -1.46, 263),
            (1908, 2, +0.34, -1.27, 263), (1893, 3, +0.15, -1.73, 270), (1905, 2, +0.34, -1.84, 270),
            (1880, 4, -0.15, -1.91, 279), (1879, 2, -0.04, -2.21, 282), (705, 4, +0.23, -2.29, 295),
            (1817, 1, -1.27, -0.34, 301), (1891, 3, +0.15, -2.10, 308)]
CSV_PATH = os.path.join(C.HERE, "data", "refuser_states_fh.csv")
NAMES = ("sigma", "alpha", "alpha_in", "z")


def variants():
    v = C.variants()
    b, c = C.UW_MODEL, C.CONE
    return {"DM04": v["DM04"], "fh": v["T1B1"], "fhc": v["T1B1S"],
            "fhc_k0.25": b.with_(h_reg="max", h_eps=c, h_soft_kappa=0.25, reseat_delta=c),
            "fhc_k0.75": b.with_(h_reg="max", h_eps=c, h_soft_kappa=0.75, reseat_delta=c)}


def kof(el, gp):
    return 4 * (el - 1) + (gp - 1)


def extract():
    rows = []
    for (el, gp, x, y, s0) in REFUSERS:
        st, info = C.ckpt_state_npz(LEG, kof(el, gp))
        assert abs(info["gx"] - x) < 0.011 and abs(info["gy"] - y) < 0.011, (el, gp, info["gx"], info["gy"])
        r = dict(leg=LEG, k=kof(el, gp), element=el, gp=gp, first_refusal_step=s0, step=info["step"],
                 s_over_B=info["s_over_B"], gx=info["gx"], gy=info["gy"], p=info["p"], psi=info["psi"],
                 e=info["raw"]["e"])
        for n in NAMES:
            for i, v in enumerate(info["raw"][n]):
                r[f"{n}_{i}"] = repr(float(v))
        rows.append(r)
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
                                     gy=float(r["gy"]), step=int(r["step"]), s_over_B=float(r["s_over_B"]),
                                     first=int(r["first_refusal_step"])))
    return out


def report(st):
    q = C.quantities(st.sigma, st.alpha, st.z, st.e, st.alpha_in, C.P, C.UW_MODEL)
    eta = math.sqrt(1.5) * C.norm(q.s) / q.p
    h_floor = q.b0 / max(q.a, C.CONE)                      # c_A = 1
    kp = (2.0 / 3.0) * q.p * h_floor * q.bn
    return dict(p=q.p, eta=eta, cos3t=q.cos3t, a_over_cone=q.a / C.CONE, bn=q.bn,
                rho_alpha=q.rho_alpha, X=q.X, HX_floor=1.0 + kp / q.X if q.X > 0 else float("nan"),
                psi=q.psi)


_STATES = None


def job(a):
    global _STATES
    if _STATES is None:
        _STATES = load_states()
    k, vname, idir, delta, deps = a
    st, _ = _STATES[k]
    t0 = time.time()
    try:
        r = C.integrate(st, deps, C.P, variants()[vname], record=False)
        out = dict(status=r.status, reseats=len(r.reseats), min_hx=getattr(r, "min_h_over_x", float("nan")),
                   sig=C.t2v(r.state.sigma).tolist())
    except Exception as ex:  # a crash is a finding
        out = dict(status="exception:" + type(ex).__name__, err=str(ex)[:200])
    out.update(k=k, variant=vname, idir=idir, delta=delta, secs=time.time() - t0)
    return out


def main(argv):
    if "--extract" in argv:
        extract()
        argv = [a for a in argv if a != "--extract"]
    states = load_states()
    V = argv or list(variants().keys())
    s0 = next(iter(states.values()))[1]
    print(f"{LEG}, c = {C.P.c}: committed states at step {s0['step']}, s/B {s0['s_over_B']:.5f}")
    print(f"{'ele/gp':>8} {'1st':>4} {'x':>6} {'y':>6} {'p':>7} {'eta':>6} {'cos3t':>6} {'a/rc':>7} "
          f"{'b:n':>9} {'rho_a':>6} {'H/X|floor':>9} {'psi':>7}")
    diag = {}
    for k, (st, info) in states.items():
        d = report(st)
        diag[k] = dict(info, **d)
        print(f"{info['element']:>5}/{info['gp']} {info['first']:>4} {info['gx']:>6.2f} {info['gy']:>6.2f} "
              f"{d['p']:>7.2f} {d['eta']:>6.3f} {d['cos3t']:>6.3f} {d['a_over_cone']:>7.3f} {d['bn']:>9.2e} "
              f"{d['rho_alpha']:>6.3f} {d['HX_floor']:>9.3f} {d['psi']:>7.3f}")
    dirs = fib_sphere(32)
    tasks = [(k, v, i, d, deps_of(u, d)) for k in states for v in V
             for d in (3e-6, 3e-5) for i, u in enumerate(dirs)]
    t0 = time.time()
    with Pool(6, maxtasksperchild=50) as pool:
        res = pool.map(job, tasks, chunksize=4)
    json.dump(dict(leg=LEG, c=C.P.c, states={str(k): v for k, v in diag.items()}, trials=res),
              open(os.path.join(C.OUT, "fan_fh_wall.json"), "w"))
    print(f"{len(res)} trials in {time.time() - t0:.0f} s")
    for v in V:
        rr = [r for r in res if r["variant"] == v]
        print(f"{v:>10}: failed {sum(1 for r in rr if r['status'] != 'ok')}/{len(rr)}",
              dict(collections.Counter(r["status"] for r in rr)))
    # kappa sensitivity of the end stress, on trials every capped variant integrates
    by = collections.defaultdict(dict)
    for r in res:
        if r["status"] == "ok" and r["variant"].startswith("fhc"):
            by[(r["k"], r["idir"], r["delta"])][r["variant"]] = r["sig"]
    rel = collections.defaultdict(list)
    for key, d in by.items():
        if len(d) == 3:
            ref = d["fhc"]
            n0 = math.sqrt(sum(x * x for x in ref))
            for kv in ("fhc_k0.25", "fhc_k0.75"):
                rel[kv].append(math.sqrt(sum((a - b) ** 2 for a, b in zip(d[kv], ref))) / n0)
    for kv, xs in rel.items():
        xs.sort()
        print(f"end stress |{kv} - fhc| / |fhc|: median {xs[len(xs) // 2]:.2e}, max {xs[-1]:.2e} over {len(xs)} trials")


if __name__ == "__main__":
    main(_ARGV)
