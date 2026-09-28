"""Is CTXu (e 0.6944, CSR 0.2) well-conditioned under DM04 ITSELF?  DM04 at three
rtol values and with a 1e-9 relative perturbation of sigma_zz at the start: if the
cycles-to-5%-DA spread is as large as the variant spread, the test cannot judge
the fix (the axisymmetric extension path is unstable).  Output: out/cyc_sensitivity.json"""
import json, sys, time
from multiprocessing import Pool
import numpy as np
import r1common as C, b_tests as B
def job(a):
    tag, rtol, pert = a
    B.RTOL_CYC = rtol
    orig = B._iso
    if pert:
        def iso(p0, e0=B.E0):
            st = orig(p0, e0); st.sigma[2, 2] *= (1 + pert); return st
        B._iso = iso
    t0 = time.time()
    r = B.cyclic_stress("CTXu", 100.0, 0.6944, 40.0, C.variants()["DM04"], n_max=20, span=0.03)
    h = r["hist"]; S = np.array(h["sigma"]); asym = float(np.max(np.abs(S[:,1]-S[:,2])))
    out = dict(tag=tag, rtol=rtol, pert=pert, status=r["status"], n_liq=r["n_liq"], da=r["da"], halves=r["halves"],
               max_asym_yy_zz=asym, secs=time.time()-t0)
    print(out, flush=True); return out
if __name__ == "__main__":
    tasks = [("rtol1e-8",1e-8,0.0),("rtol1e-9",1e-9,0.0),("rtol1e-10",1e-10,0.0),("pert+1e-9",1e-9,1e-9),("pert-1e-9",1e-9,-1e-9)]
    with Pool(5, maxtasksperchild=1) as pool: res = pool.map(job, tasks, chunksize=1)
    json.dump(res, open(C.OUT + "/cyc_sensitivity.json","w"), indent=1)
