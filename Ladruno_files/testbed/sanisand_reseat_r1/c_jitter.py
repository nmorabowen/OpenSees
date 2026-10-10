"""(c) Near-neutral jittering loading from the real near-peak refuser states.

c1  CHAINS: N committed increments of norm DELTA from each refuser state,
    d_k = unit(d0 + sigma * xi_k) (xi_k random unit, seeded), sigma in
    {0.3, 1, 3}; plus a REVERSAL chain d_k = s_k d0 + 0.3 xi_k with random
    signs s_k (small back-and-forth around the state).  Per chain: increments
    completed before the first non-ok status, re-seats, min H/X, max rho_alpha.
c2  OBJECTIVITY: E_B 1880/1 and E_B16 7820/4, a fixed noise sequence, jitter
    amplitude lam in {0, 0.01, 0.03, 0.1, 0.3, 1}: d_k = d0 + lam xi_k (not
    normalised).  The end-state stress after N increments vs lam = 0.
c3  CONTINUITY: single trial increments d(theta) = DELTA (cos theta d0 + sin
    theta d1), 181 angles; the stress increment as a function of the direction.
d0 = the GP's last converged strain increment (checkpoint difference).
Output: out/c1.json, out/c2.json, out/c3.json"""
from __future__ import annotations

import json
import math
import sys
import time
from multiprocessing import Pool

import numpy as np

import r1common as C

DELTA = 3e-6
N_CHAIN = 60
WORKERS = 6
VAR = ["DM04", "B1", "R150", "R150+T1", "T1B1", "T1B1S", "T2B1S"]
PREV = {"E_B": "field_step00375.npz", "E_D": "field_step00385.npz",
        "E_B16": "field_step00070.npz"}


def mandel_unit(v):
    """plane-strain Voigt (xx, yy, 0, gamma, 0, 0) -> unit in Mandel metric."""
    m = np.array([v[0], v[1], v[3] / math.sqrt(2.0)])
    return m / np.linalg.norm(m)


def voigt(m, scale):
    return [scale * m[0], scale * m[1], 0.0, scale * m[2] * math.sqrt(2.0), 0.0, 0.0]


def d0_of(leg, k):
    return mandel_unit(C.last_increment(leg, k, PREV[leg]))


def rand_units(seed, n):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 3))
    return x / np.linalg.norm(x, axis=1)[:, None]


def run_chain(st, incs, O):
    n_ok, reseats, min_hx, max_ra = 0, 0, float("inf"), 0.0
    status = "ok"
    qs = []
    for de in incs:
        try:
            r = C.integrate(st, de, C.P, O, record=False)
        except Exception as ex:
            status = "exception:" + type(ex).__name__
            break
        reseats += len(r.reseats)
        min_hx = min(min_hx, r.min_h_over_x)
        max_ra = max(max_ra, r.max_rho_alpha)
        if r.status != "ok":
            status = r.status
            break
        st = r.state
        n_ok += 1
        s = st.sigma
        p = np.trace(s) / 3.0
        qs.append([float(p), float(math.sqrt(1.5) * C.norm(s - p * np.eye(3)))])
    return dict(status=status, n_ok=n_ok, reseats=reseats, min_hx=min_hx,
                max_rho_alpha=max_ra, sigma_end=C.t2v(st.sigma).tolist(), pq=qs)


def job_c1(a):
    leg, k, kind, sig, v = a
    st, _ = C.ckpt_state(leg, k)
    d0 = d0_of(leg, k)
    xi = rand_units(1000 + k, N_CHAIN)
    if kind == "jitter":
        incs = [voigt((d0 + sig * x) / np.linalg.norm(d0 + sig * x), DELTA) for x in xi]
    else:  # reversal chain
        rng = np.random.default_rng(2000 + k)
        sg = rng.choice([-1.0, 1.0], size=N_CHAIN)
        incs = [voigt((s * d0 + 0.3 * x) / np.linalg.norm(s * d0 + 0.3 * x), DELTA)
                for s, x in zip(sg, xi)]
    t0 = time.time()
    out = run_chain(st, incs, C.variants()[v])
    out.update(leg=leg, k=k, kind=kind, sigma=sig, variant=v, secs=time.time() - t0)
    return ("c1", out)


def job_c2(a):
    leg, k, lam, v = a
    st, _ = C.ckpt_state(leg, k)
    d0 = d0_of(leg, k)
    xi = rand_units(3000 + k, N_CHAIN)
    incs = [voigt(d0 + lam * x, DELTA) for x in xi]
    t0 = time.time()
    out = run_chain(st, incs, C.variants()[v])
    out.update(leg=leg, k=k, lam=lam, variant=v, secs=time.time() - t0)
    return ("c2", out)


def job_c3(a):
    leg, k, th, v = a
    st, _ = C.ckpt_state(leg, k)
    d0 = d0_of(leg, k)
    e = np.array([0.0, 0.0, 1.0]) if abs(d0[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    d1 = e - (e @ d0) * d0
    d1 /= np.linalg.norm(d1)
    de = voigt(math.cos(th) * d0 + math.sin(th) * d1, DELTA)
    try:
        r = C.integrate(st, de, C.P, C.variants()[v], record=False)
        out = dict(status=r.status, reseats=len(r.reseats), min_hx=r.min_h_over_x,
                   dsig=(C.t2v(r.state.sigma) - C.t2v(st.sigma)).tolist(),
                   dalpha=(C.t2v(r.state.alpha) - C.t2v(st.alpha)).tolist())
    except Exception as ex:
        out = dict(status="exception:" + type(ex).__name__)
    out.update(leg=leg, k=k, theta=th, variant=v)
    return ("c3", out)


def dispatch(a):
    return {"c1": job_c1, "c2": job_c2, "c3": job_c3}[a[0]](a[1])


def main():
    which = sys.argv[1:] or ["c1", "c2", "c3"]
    tasks = []
    if "c1" in which:
        for (leg, k, el, gp, prev) in C.REFUSERS:
            for kind, sig in (("jitter", 0.3), ("jitter", 1.0), ("jitter", 3.0), ("reversal", 0.3)):
                for v in VAR:
                    tasks.append(("c1", (leg, k, kind, sig, v)))
    if "c2" in which:
        for (leg, k) in (("E_B", 7516), ("E_B16", 31279)):
            for lam in (0.0, 0.01, 0.03, 0.1, 0.3, 1.0):
                for v in VAR:
                    tasks.append(("c2", (leg, k, lam, v)))
    if "c3" in which:
        for (leg, k) in (("E_B", 7516), ("E_B16", 31279)):
            for i in range(181):
                th = 2 * math.pi * i / 180
                for v in ("DM04", "B1", "R150", "T1B1", "T1B1S"):
                    tasks.append(("c3", (leg, k, th, v)))
    print(len(tasks), "tasks", flush=True)
    res = {"c1": [], "c2": [], "c3": []}
    t0 = time.time()
    with Pool(WORKERS) as pool:
        for n, (kind, r) in enumerate(pool.imap_unordered(dispatch, tasks, chunksize=2)):
            res[kind].append(r)
            if (n + 1) % 100 == 0:
                print(f"{n+1}/{len(tasks)} {time.time()-t0:.0f}s", flush=True)
    for kind in which:
        json.dump(res[kind], open(f"{C.OUT}/{kind}.json", "w"))
    print("done", time.time() - t0)


if __name__ == "__main__":
    main()
