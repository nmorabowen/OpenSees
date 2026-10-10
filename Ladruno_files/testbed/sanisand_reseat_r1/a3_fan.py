"""(a) singular states, part 1: a FAN of trial strain increments from the committed
states of the real loadingNonPosH refusers (E_B 1880/1, 1879/1; E_D 1962/1,
2058/1; E_B16 7820/4 -- orchestrator's floor_refusers.csv), at the Newton-iterate
scale.  For each (state, direction, variant): integrate one increment exactly
(Radau, rtol 1e-10) and record the status, re-seats, min H/X, max rho_alpha and
the end stress.  Output: out/a3_fan.json + out/a3_fan.md.

Directions: a Fibonacci sphere in the plane-strain strain space, Mandel-scaled
(eps_xx, eps_yy, gamma_xy/sqrt2), magnitude DELTA (Voigt engineering shear on
input to the oracle)."""
from __future__ import annotations

import json
import math
import sys
import time
from multiprocessing import Pool

import numpy as np

import r1common as C

NDIR = int(sys.argv[1]) if len(sys.argv) > 1 else 32
ONLY = sys.argv[2:]          # variant names (default: every variant in r1common.variants())
DELTAS = [3e-6, 3e-5]
WORKERS = 6


def fib_sphere(n):
    pts = []
    ga = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        y = 1.0 - 2.0 * (i + 0.5) / n
        r = math.sqrt(max(0.0, 1.0 - y * y))
        th = ga * i
        pts.append((r * math.cos(th), y, r * math.sin(th)))
    return pts


def deps_of(u, delta):
    ex, ey, gs = u                       # Mandel: gamma/sqrt2
    return [delta * ex, delta * ey, 0.0, delta * gs * math.sqrt(2.0), 0.0, 0.0]


def job(args):
    sid, vname, idir, delta, deps = args
    leg, k = sid
    st, _ = C.ckpt_state(leg, k)
    O = C.variants()[vname]
    t0 = time.time()
    try:
        r = C.integrate(st, deps, C.P, O, record=False)
        out = dict(status=r.status, t_end=r.t_end, reseats=len(r.reseats),
                   min_hx=getattr(r, "min_h_over_x", float("nan")),
                   n_cap=getattr(r, "n_cap", [0, 0, 0]),
                   max_rho_alpha=r.max_rho_alpha, f_end=r.f_end,
                   dsig=(C.t2v(r.state.sigma) - C.t2v(st.sigma)).tolist(),
                   dalpha=(C.t2v(r.state.alpha) - C.t2v(st.alpha)).tolist(),
                   segments=len(r.segments), notes=r.notes[-2:])
    except Exception as ex:  # a crash is a finding, keep going
        out = dict(status="exception:" + type(ex).__name__, err=str(ex)[:200])
    out.update(leg=leg, k=k, variant=vname, idir=idir, delta=delta, deps=deps,
               secs=time.time() - t0)
    return out


def main():
    dirs = fib_sphere(NDIR)
    V = ONLY or list(C.variants().keys())
    tasks = []
    for (leg, k, el, gp, prev) in C.REFUSERS:
        for delta in DELTAS:
            for i, u in enumerate(dirs):
                for v in V:
                    tasks.append(((leg, k), v, i, delta, deps_of(u, delta)))
    print(len(tasks), "tasks", flush=True)
    t0 = time.time()
    res = []
    with Pool(WORKERS) as pool:
        for n, r in enumerate(pool.imap_unordered(job, tasks, chunksize=4)):
            res.append(r)
            if (n + 1) % 200 == 0:
                print(f"{n+1}/{len(tasks)}  {time.time()-t0:.0f}s", flush=True)
    json.dump(res, open(C.OUT + ("/a3_fan_" + "+".join(ONLY) + ".json" if ONLY else "/a3_fan.json"), "w"))
    print("done", time.time() - t0)


if __name__ == "__main__":
    main()
