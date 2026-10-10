"""(a) singular states, part 2: the task's named states.
  a1  b8 1950/3 (WP-134's 0/0 point under shear) and 1950/2, the WP-127 probes
      isoComp (d, d, 0) and shear gamma_xy = d at d = 1e-5 / 1e-4 / 1e-3;
  a2  WP-128's smallest reproducer: sigma = 0.0101 I, alpha = alpha_in = z = 0,
      one plane-strain d eps_yy = 1e-4 (and 1e-5, 3e-4), plus p_s = 1 kPa.
Every variant; exact integration.  Output: out/a12.json"""
from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool

import r1common as C

WORKERS = 3
ONLY = sys.argv[1:]          # variant names (default: all)


def job(args):
    case, vname, st_key, deps = args
    if st_key[0] == "ring":
        st, _ = C.ring_state("b8", st_key[1], st_key[2])
    else:
        from sanisand_r1 import ring as rring
        st = rring.reproducer_state(p_s=st_key[1])
    O = C.variants()[vname]
    t0 = time.time()
    try:
        r = C.integrate(st, deps, C.P, O, record=False)
        s = r.summary()
        out = dict(status=r.status, t_end=r.t_end, reseats=len(r.reseats),
                   min_hx=r.min_h_over_x, n_cap=r.n_cap, max_rho_alpha=r.max_rho_alpha,
                   rho_alpha_end=s["rho_alpha_end"], eta_end=s["eta_end"], p_end=s["p_end"],
                   f_end=r.f_end, sigma=s["sigma"], alpha=s["alpha"], notes=r.notes[-2:])
    except Exception as ex:
        out = dict(status="exception:" + type(ex).__name__, err=str(ex)[:300])
    out.update(case=case, variant=vname, secs=time.time() - t0)
    return out


def main():
    V = ONLY or list(C.variants().keys())
    tasks = []
    for el, gp in [(1950, 3), (1950, 2)]:
        for d in (1e-5, 1e-4, 1e-3):
            for pn, de in (("isoComp", [d, d, 0, 0, 0, 0]), ("shear", [0, 0, 0, d, 0, 0])):
                for v in V:
                    tasks.append((f"b8 {el}/{gp} {pn} {d:g}", v, ("ring", el, gp), de))
    for ps, d in ((0.0101, 1e-5), (0.0101, 1e-4), (0.0101, 3e-4), (1.0, 1e-4)):
        for v in V:
            tasks.append((f"reproducer p_s {ps:g} deps_yy {d:g}", v, ("repro", ps), [0, d, 0, 0, 0, 0]))
    print(len(tasks), "tasks", flush=True)
    t0 = time.time()
    res = []
    with Pool(WORKERS) as pool:
        for n, r in enumerate(pool.imap_unordered(job, tasks)):
            res.append(r)
            if (n + 1) % 20 == 0:
                print(f"{n+1}/{len(tasks)} {time.time()-t0:.0f}s", flush=True)
    json.dump(res, open(C.OUT + ("/a12_extra.json" if ONLY else "/a12.json"), "w"))
    print("done", time.time() - t0)


if __name__ == "__main__":
    main()
