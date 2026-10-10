"""DM04-only pilot: cycles to liquefaction for candidate amplitudes (to pick the
(b) cyclic amplitudes).  Output: out/cyc_pilot.json"""
import json
import sys
import time
from multiprocessing import Pool

import r1common as C
import b_tests as B


def job(a):
    kind, p0, e0, ratio = a
    O = C.variants()["DM04"]
    t0 = time.time()
    amp = ratio * (2 * p0 if kind == "CTXu" else p0)   # CSR = q/(2 p0) or tau/sigma_v0
    try:
        r = B.cyclic_stress(kind, p0, e0, amp, O, n_max=25)
        out = dict(status=r["status"], n_liq=r["n_liq"], halves=r["halves"], da=r["da"],
                   reseats=r["reseats"], p_end=r["p_end"])
    except Exception as ex:
        out = dict(status="exception:" + str(ex)[:200])
    out.update(kind=kind, p0=p0, e0=e0, csr=ratio, secs=time.time() - t0)
    print(out, flush=True)
    return out


if __name__ == "__main__":
    tasks = [("CTXu", 100.0, 0.80, r) for r in (0.10, 0.15, 0.20)] + \
            [("CTXu", 100.0, 0.6944, r) for r in (0.20, 0.30)] + \
            [("CSSu", 100.0, 0.80, r) for r in (0.08, 0.12)] + \
            [("CSSu", 100.0, 0.6944, r) for r in (0.15, 0.25)]
    with Pool(3) as pool:
        res = pool.map(job, tasks)
    json.dump(res, open(C.OUT + "/cyc_pilot.json", "w"), indent=1)
