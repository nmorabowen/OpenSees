"""#893 review M1: R1's effect on the SMALL-STRAIN response (exact oracle, campaign set).

The h floor bounds h right after a re-seat, where DM04 has h = inf (a purely elastic
start with the cone dragged along).  So R1 starts plastic flow at once, over roughly
one cone radius of back-stress travel, and the small-strain stiffness drops a little.
From the isotropic state p0 = 100 kPa, e0 = 0.6944 (alpha = alpha_in = 0, fabric 0):
  CVSS  constant-volume simple shear, one exact integration to each gamma
  TCu   undrained (isochoric) triaxial compression, to each eps_a
  cyc   10 strain-controlled constant-volume cycles, gamma = +-1e-5 (legs integrated
        exactly), tau compared at the same strain
Variants: DM04, R1 (c_A 1, c_rev 1, kappa 0.5) and R1 with c_rev 2 (T2B1S).
Output: out/m1_small_strain.json"""
from __future__ import annotations

import json
import math
import os

import numpy as np

import r1common as C
from sanisand_r1 import driver
from sanisand_r1.integrator import Control, integrate, path_table

E0 = 0.6944
V = {"DM04": C.variants()["DM04"], "R1": C.variants()["T1B1S"], "R1_crev2": C.variants()["T2B1S"]}


def iso():
    return driver.isotropic_state(100.0, E0)


def one(ctl, O):
    r = integrate(iso(), ctl, C.P, O, rtol=1e-12)
    s = C.t2v(r.state.sigma)
    return r.status, s, len(r.reseats)


def cvss_cycles(amp, O, n_cyc=10):
    st = iso()
    legs = [amp] + [(-2 * amp if k % 2 == 0 else 2 * amp) for k in range(2 * n_cyc)]
    g_acc, G, T, reseats, held = 0.0, [], [], 0, 0
    for d in legs:
        r = integrate(st, Control.strain([0, 0, 0, d, 0, 0]), C.P, O, rtol=1e-10)
        reseats += len(r.reseats)
        held += getattr(r, "n_held", 0)
        for row in path_table(r, C.P, O)[1:]:
            G.append(g_acc + row["eps"][3])
            T.append(row["sigma"][3])
        g_acc += d
        st = r.state
        if r.status != "ok":
            return dict(status=r.status, g=G, tau=T, reseats=reseats)
    return dict(status="ok", g=G, tau=T, reseats=reseats)


def main():
    out = {}
    print("CVSS: tau(R1)/tau(DM04) - 1 at the same gamma")
    for g in (1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4):
        row = {}
        for name, O in V.items():
            stt, s, nr = one(Control.strain([0, 0, 0, g, 0, 0]), O)
            row[name] = dict(status=stt, tau=s[3], reseats=nr)
        d1 = row["R1"]["tau"] / row["DM04"]["tau"] - 1.0
        d2 = row["R1_crev2"]["tau"] / row["DM04"]["tau"] - 1.0
        out[f"CVSS/{g:g}"] = row
        print(f"   gamma {g:.0e}: DM04 tau {row['DM04']['tau']:.5g} kPa; R1 {d1:+.3%}; R1 c_rev 2 {d2:+.3%}")
    print("TCu: q(R1)/q(DM04) - 1 at the same eps_a")
    for ea in (1e-6, 3e-6, 1e-5, 3e-5, 1e-4):
        row = {}
        for name, O in V.items():
            stt, s, nr = one(Control.strain([ea, -ea / 2, -ea / 2, 0, 0, 0]), O)
            row[name] = dict(status=stt, q=s[0] - s[1], reseats=nr)
        d1 = row["R1"]["q"] / row["DM04"]["q"] - 1.0
        out[f"TCu/{ea:g}"] = row
        print(f"   eps_a {ea:.0e}: DM04 q {row['DM04']['q']:.5g} kPa; R1 {d1:+.3%}")
    print("cyclic CVSS, gamma = +-1e-5, 10 cycles: max |tau_R1 - tau_DM04| / max|tau_DM04|, same strain")
    cyc = {name: cvss_cycles(1e-5, O) for name, O in V.items()}
    ref = cyc["DM04"]
    tmax = max(abs(t) for t in ref["tau"])
    for name in ("R1", "R1_crev2"):
        c = cyc[name]
        # the legs are identical in strain; compare on a common gamma grid per leg by
        # interpolation along the path (both paths are monotone within a leg)
        diffs = []
        gr, tr = np.array(ref["g"]), np.array(ref["tau"])
        gc, tc = np.array(c["g"]), np.array(c["tau"])
        # split into legs at the direction changes of gamma
        def legs_of(g):
            idx = [0] + [i for i in range(1, len(g) - 1) if (g[i] - g[i - 1]) * (g[i + 1] - g[i]) < 0] + [len(g) - 1]
            return list(zip(idx[:-1], idx[1:]))
        for (a0, a1), (b0, b1) in zip(legs_of(gr), legs_of(gc)):
            ga, ta = gr[a0:a1 + 1], tr[a0:a1 + 1]
            gb, tb = gc[b0:b1 + 1], tc[b0:b1 + 1]
            if ga[-1] < ga[0]:
                ga, ta, gb, tb = ga[::-1], ta[::-1], gb[::-1], tb[::-1]
            grid = np.linspace(max(ga[0], gb[0]), min(ga[-1], gb[-1]), 50)
            diffs.append(np.max(np.abs(np.interp(grid, gb, tb) - np.interp(grid, ga, ta))))
        print(f"   {name}: status {c['status']}, re-seats {c['reseats']} (DM04 {ref['reseats']}); "
              f"max |dtau| / tau_max = {max(diffs) / tmax:.2e}")
        out[f"cyc/{name}"] = dict(status=c["status"], reseats=c["reseats"], rel=max(diffs) / tmax)
    json.dump(out, open(os.path.join(C.OUT, "m1_small_strain.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
