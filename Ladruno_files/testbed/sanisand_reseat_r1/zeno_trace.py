"""The Zeno accumulation of alpha_in re-seats, traced (for the memo's figure 1).

E_B 1880/1 (committed wall state), one trial increment of norm 3e-6 along
fib_sphere(12)[4] (a chatter direction of the smoke fan).  For DM04, B1 (floor
only), T1B1 and T1B1S: the re-seat times t_k, and at each re-seat b:n, H/X and
|dalpha/dt| (per unit pseudo-time of the increment).  Output: out/zeno_trace.json"""
from __future__ import annotations

import json

import numpy as np

import r1common as C
from a3_fan import deps_of, fib_sphere
from sanisand_r1.integrator import Control, _Increment


def trace(vname):
    st, _ = C.ckpt_state("E_B", 7516)
    de = deps_of(fib_sphere(12)[4], 3e-6)
    O = C.variants()[vname]
    r = C.integrate(st, de, C.P, O, record=True)
    T = np.array(r.path["t"]); Y = np.array(r.path["y"])
    inc = _Increment(st, Control.strain(de), C.P, O, 1e-10, 1e-2, "Radau")
    rows = []
    for k, rs in enumerate(r.reseats):
        t = rs["t"]
        j = int(np.argmin(np.abs(T - t)))
        y = Y[j]
        ain_before = None
        q = C.quantities(C.v2t(y[0:6]), C.v2t(y[6:12]), C.v2t(y[12:18]), float(y[18]),
                         C.v2t(y[6:12]), C.P, O)          # just re-seated: alpha_in = alpha
        inc.alpha_in = C.v2t(y[6:12]).copy()
        dy, info = inc.rates(y, "plastic", q)
        rows.append(dict(k=k, t=float(t), bn=float(q.bn), a=float(q.a),
                         dalpha=float(C.norm(C.v2t(dy[6:12]))), N=float(info["N"]),
                         cos3t=float(q.cos3t), rho_b=float(q.rho_b)))
    return dict(variant=vname, status=r.status, t_end=float(r.t_end), n_reseats=len(r.reseats),
                reseats=rows, min_hx=float(r.min_h_over_x))


if __name__ == "__main__":
    out = [trace(v) for v in ("DM04", "B1", "T1B1", "T1B1S")]
    for o in out:
        print(o["variant"], o["status"], o["t_end"], o["n_reseats"],
              [f"{r['t']:.9f}/{r['bn']:.2e}/{r['dalpha']:.2e}" for r in o["reseats"][:8]])
    json.dump(out, open(C.OUT + "/zeno_trace.json", "w"), indent=1)
