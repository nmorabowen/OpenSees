"""Does the extension-side non-convexity (c = 0.71 < 7/9) feed the Zeno wall?
The wall fan (a3_fan.py: 5 refuser states x 32 directions x {3e-6, 3e-5}) with the
Lode parameter set to c = 0.80 (convex) and to the campaign's 0.71, DM04 only
(the paper model; the states are the c = 0.71 committed wall states, so this
measures SENSITIVITY, not a c = 0.80 history).  Output: out/fan_c080.json"""
from __future__ import annotations

import dataclasses
import json
import time
from multiprocessing import Pool

import r1common as C
from a3_fan import deps_of, fib_sphere


def job(a):
    (leg, k), cval, idir, delta, deps = a
    P = dataclasses.replace(C.P, c=cval)
    st, _ = C.ckpt_state(leg, k)
    t0 = time.time()
    try:
        r = C.integrate(st, deps, P, C.variants()["DM04"], record=False)
        q = C.quantities(st.sigma, st.alpha, st.z, st.e, st.alpha_in, P, C.UW_MODEL)
        out = dict(status=r.status, reseats=len(r.reseats), rho_alpha0=q.rho_alpha, rho_b0=q.rho_b)
    except Exception as ex:
        out = dict(status="exception:" + type(ex).__name__)
    out.update(leg=leg, k=k, c=cval, idir=idir, delta=delta, secs=time.time() - t0)
    return out


if __name__ == "__main__":
    dirs = fib_sphere(32)
    tasks = [((leg, k), cval, i, d, deps_of(u, d))
             for (leg, k, el, gp, prev) in C.REFUSERS for cval in (0.71, 0.80)
             for d in (3e-6, 3e-5) for i, u in enumerate(dirs)]
    with Pool(6, maxtasksperchild=50) as pool:
        res = pool.map(job, tasks, chunksize=4)
    json.dump(res, open(C.OUT + "/fan_c080.json", "w"))
    import collections
    for cval in (0.71, 0.80):
        rr = [r for r in res if r["c"] == cval]
        print(f"c = {cval}: failed {sum(1 for r in rr if r['status'] != 'ok')}/{len(rr)}",
              dict(collections.Counter(r["status"] for r in rr)))
        for (leg, k, el, gp, prev) in C.REFUSERS:
            s = [r for r in rr if r["leg"] == leg and r["k"] == k]
            print(f"   {leg} {el}/{gp}: failed {sum(1 for r in s if r['status'] != 'ok')}/64, "
                  f"rho_alpha0 {s[0].get('rho_alpha0', float('nan')):.3f}, rho_b0 {s[0].get('rho_b0', float('nan')):.3f}")
