"""BLOCKING-GATE study for the CTXu e0.6944 CSR 0.2 outlier.
(1) every variant x perturbation of sigma_zz (relative) {-1e-8,-1e-9,0,+1e-9,+1e-8}:
    does each variant show the same two branches as DM04?
(2) control with c = 0.80 (> 7/9, a convex Lode function): DM04 and T1B1S x
    {-1e-9, 0, +1e-9}: does the path stay axisymmetric and the branches merge?
Records N at 5% DA (or None), DA after the run, max |sigma_yy - sigma_zz|, and the
first eps_a where |sigma_yy - sigma_zz| exceeds 1e-3 kPa (symmetry breaking).
Output: out/cyc_gate.json"""
import json, time, dataclasses
from multiprocessing import Pool
import numpy as np
import r1common as C, b_tests as B
_ISO0 = B._iso          # pristine: never wrap a wrapper (perturbations must not compound)
_P0 = C.P
def job(a):
    tag, vname, pert, cval = a
    C.P = dataclasses.replace(_P0, c=cval) if cval is not None else _P0   # explicit, every task
    B.C.P = C.P
    def iso(p0, e0=B.E0):
        st = _ISO0(p0, e0); st.sigma[2, 2] *= (1 + pert); return st
    B._iso = iso
    t0 = time.time()
    try:
        r = B.cyclic_stress("CTXu", 100.0, 0.6944, 40.0, C.variants()[vname], n_max=20, span=0.03)
        h = r["hist"]; S = np.array(h["sigma"]); E = np.array(h["eps"])[:, 0]
        asy = np.abs(S[:, 1] - S[:, 2]); k = int(np.argmax(asy > 1e-3)) if (asy > 1e-3).any() else None
        out = dict(status=r["status"], n_liq=r["n_liq"], da=float(r["da"]), halves=r["halves"],
                   max_asym=float(asy.max()), eps_break=(float(E[k]) if k is not None else None),
                   half_break=(int(h["half"][k]) if k is not None else None), reseats=r["reseats"])
    except Exception as ex:
        out = dict(status="exception:" + str(ex)[:200])
    out.update(tag=tag, variant=vname, pert=pert, c=cval, secs=time.time() - t0)
    print(out, flush=True); return out
if __name__ == "__main__":
    tasks = []
    for v in ["DM04", "B1", "T1B1", "T1B1S", "T2B1"]:
        for pz in (-1e-8, -1e-9, 0.0, 1e-9, 1e-8):
            tasks.append(("c0.71", v, pz, None))
    for v in ["DM04", "T1B1S"]:
        for pz in (-1e-9, 0.0, 1e-9):
            tasks.append(("c0.80", v, pz, 0.80))
    with Pool(4, maxtasksperchild=1) as pool: res = pool.map(job, tasks, chunksize=1)
    json.dump(res, open(C.OUT + "/cyc_gate.json", "w"), indent=1)
