"""CTXu gate, calibration-independent control: DM04's own TOYOURA set (A0 0.704,
m 0.01, c 0.712).  Stage 1: DM04 pilot over (e0, CSR) at p0 = 100 kPa; stage 2:
for the pilot cases with 3 <= N <= 20, DM04 / B1 / T1B1 / T1B1S x sigma_zz
perturbation {-1e-9, 0, +1e-9}.  The R1 floor and threshold scale with Toyoura's
own cone radius sqrt(2/3) * 0.01.  Output: out/cyc_toyoura.json"""
from __future__ import annotations

import dataclasses
import json
import time
from multiprocessing import Pool

import numpy as np

import r1common as C
import b_tests as B
from sanisand_r1.model import TOYOURA, SQ23

_ISO0 = B._iso          # pristine: never wrap a wrapper (perturbations must not compound)


def run(a):
    stage, vname, e0, csr, pert = a
    P = dataclasses.replace(TOYOURA, e_init=e0)
    C.P = P
    B.C.P = P
    def iso(p0, e_=e0):
        st = _ISO0(p0, e_)
        st.sigma[2, 2] *= (1 + pert)
        return st
    B._iso = iso
    O = C.variants(cone=SQ23 * P.m)[vname]
    t0 = time.time()
    try:
        r = B.cyclic_stress("CTXu", 100.0, e0, csr * 200.0, O, n_max=20, span=0.03)
        S = np.array(r["hist"]["sigma"])
        out = dict(status=r["status"], n_liq=r["n_liq"], da=float(r["da"]), halves=r["halves"],
                   max_asym=float(np.abs(S[:, 1] - S[:, 2]).max()), reseats=r["reseats"],
                   p_end=r["p_end"])
    except Exception as ex:
        out = dict(status="exception:" + str(ex)[:200])
    out.update(stage=stage, variant=vname, e0=e0, csr=csr, pert=pert, secs=time.time() - t0)
    print(out, flush=True)
    return out


if __name__ == "__main__":
    pilot = [("pilot", "DM04", e0, csr, 0.0) for e0, csrs in ((0.808, (0.10, 0.15, 0.20)),
                                                            (0.735, (0.20, 0.30)))
             for csr in csrs]
    with Pool(3, maxtasksperchild=1) as pool:
        res = pool.map(run, pilot, chunksize=1)
    picks = [(r["e0"], r["csr"]) for r in res
             if r.get("n_liq") is not None and 3 <= r["n_liq"] <= 20][:2]
    print("picked", picks, flush=True)
    gate = [("gate", v, e0, csr, pz) for (e0, csr) in picks
            for v in ("DM04", "B1", "T1B1", "T1B1S") for pz in (-1e-9, 0.0, 1e-9)]
    with Pool(3, maxtasksperchild=1) as pool:
        res += pool.map(run, gate, chunksize=1)
    json.dump(res, open(C.OUT + "/cyc_toyoura.json", "w"), indent=1)
