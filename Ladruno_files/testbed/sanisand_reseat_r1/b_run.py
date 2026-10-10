"""(b) runner: every (test, variant) -> out/b/<test>__<variant>.json (resumable)."""
from __future__ import annotations

import json
import os
import sys
import time
from multiprocessing import Pool

import r1common as C
import b_tests as B

OUTB = os.path.join(C.OUT, "b")
os.makedirs(OUTB, exist_ok=True)

CYC = {
    "CTXu_e0.6944_csr0.20": lambda O: B.cyclic_stress("CTXu", 100.0, 0.6944, 0.20 * 200.0, O, n_max=20, span=0.03),
    "CTXu_e0.80_csr0.10": lambda O: B.cyclic_stress("CTXu", 100.0, 0.80, 0.10 * 200.0, O, n_max=20, span=0.03),
    "CSSu_e0.80_csr0.08": lambda O: B.cyclic_stress("CSSu", 100.0, 0.80, 0.08 * 100.0, O, n_max=20, span=0.03),
    "CSSu_e0.6944_csr0.25": lambda O: B.cyclic_stress("CSSu", 100.0, 0.6944, 0.25 * 100.0, O, n_max=15, span=0.03),
    "CSSd_g0.001": lambda O: B.cyclic_strain("CSSd", 100.0, 0.6944, 0.001, O, n_cyc=10),
    "CSSd_g0.005": lambda O: B.cyclic_strain("CSSd", 100.0, 0.6944, 0.005, O, n_cyc=10),
    "CTXd_e0.0005": lambda O: B.cyclic_strain("CTXd", 100.0, 0.6944, 0.0005, O, n_cyc=6),
}
ALL = dict(B.TESTS)
ALL.update(CYC)
VARIANTS = ["DM04", "B1", "Badd1", "T1B1", "T2B1", "T2B0.25", "T2B1S"]


def job(a):
    test, v = a
    path = os.path.join(OUTB, f"{test}__{v}.json")
    if os.path.exists(path):
        return test, v, "cached", 0.0
    t0 = time.time()
    try:
        r = ALL[test](C.variants()[v])
    except Exception as ex:
        r = dict(status="exception:" + type(ex).__name__ + ": " + str(ex)[:300])
    r["secs"] = time.time() - t0
    json.dump(r, open(path, "w"), default=float)
    return test, v, r.get("status"), r["secs"]


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "--variants":
        VARIANTS = args[1].split(",")
        args = args[2:]
    only = args or list(ALL)
    tasks = [(t, v) for t in only for v in VARIANTS]
    # longest first
    tasks.sort(key=lambda a: 0 if a[0] in CYC else 1)
    t0 = time.time()
    with Pool(6) as pool:
        for n, (t, v, s, secs) in enumerate(pool.imap_unordered(job, tasks)):
            print(f"{n+1}/{len(tasks)} {t} {v}: {s} ({secs:.0f}s)  [{time.time()-t0:.0f}s]", flush=True)
