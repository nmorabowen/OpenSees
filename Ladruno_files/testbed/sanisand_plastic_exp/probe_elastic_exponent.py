"""WP-229 probe (NOT an option, nothing committed to the model): the 17-test gate with the ELASTIC exponent
free as well, to reconcile the gate with the F26 intake, whose "all 17 within +-10 %" is its cell XC
(nh 0.206, ng 0.756). The elastic law is monkey-patched here: G = G0 p_at F(e) (p/p_at)^ng.
Same sampling as gate_nh_oracle.py (500 equal increments).   python probe_elastic_exponent.py
"""
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "Ladruno_scripts"))
sys.path.insert(0, HERE)
CELLS = {"DM04": (0.5, 0.5), "XC": (0.206, 0.756), "XB exponents, DM04 constants": (0.251, 0.619),
         "nh 0.25 alone": (0.25, 0.5)}


def one(job):
    name, nh, ng, label, drained, e0, p0 = job
    import numpy as np
    import sanisand_reference.model as M
    import sanisand_reference.integrator as I
    from sanisand_reference import TOYOURA
    from sanisand_reference.driver import isotropic_state
    from sanisand_reference.ring import ring_variants
    # The patch is set for EVERY job, ng = 0.5 included: pool workers are reused, and a patch left by an
    # earlier job would otherwise leak into the next one (it did, in the first version of this probe).
    def em(p, e, P, O):
        pG = max(max(p, O.p_min) if O.p_min > 0.0 else p, 0.0)
        eG = P.e_init if O.g_void_ratio == "initial" else e
        G = P.G0 * P.P_atm * (2.97 - eG) ** 2 / (1.0 + eG) * (pG / P.P_atm) ** ng
        return G, 2.0 * (1.0 + P.nu) / (3.0 * (1.0 - 2.0 * P.nu)) * G
    M.elastic_moduli = em
    I.elastic_moduli = em
    P = replace(TOYOURA, e_init=e0)
    O = ring_variants()["uw_model"].with_(plastic_exp=nh, start_outside="plastic")
    da = 0.25 / 500
    st = isotropic_state(p0, e0)
    ctl = (I.Control((True, False, False, True, True, True), (da, 0, 0, 0, 0, 0)) if drained
           else I.Control.strain([da, -0.5 * da, -0.5 * da, 0, 0, 0]))
    q, status = [0.0], "ok"
    for k in range(500):
        r = I.integrate(st, ctl, P, O, rtol=1e-8, record=False)
        if r.status != "ok":
            status = r.status
            break
        st = r.state
        pk = float(np.trace(st.sigma)) / 3.0
        q.append(math.sqrt(1.5) * M.norm(st.sigma - pk * np.eye(3)))
    return name, label, status, max(q), q[-1]


if __name__ == "__main__":
    from gate_nh_oracle import TESTS
    jobs = [(n, a, b) + t for n, (a, b) in CELLS.items() for t in TESTS]
    with ProcessPoolExecutor(12, max_tasks_per_child=1) as ex:
        R = {(r[0], r[1]): r for r in ex.map(one, jobs)}
    L = ["| test | " + " | ".join(f"{n} (nh {a}, ng {b})" for n, (a, b) in CELLS.items() if n != "DM04") + " |",
         "|---|" + "---|" * (len(CELLS) - 1)]
    for t in TESTS:
        ref = R[("DM04", t[0])]
        row = f"| {t[0]} |"
        for n in CELLS:
            if n == "DM04":
                continue
            r = R[(n, t[0])]
            row += f" {r[3] / ref[3]:.3f} / {r[4] / ref[4]:.3f}{'' if r[2] == 'ok' else ' (' + r[2] + ')'} |"
        L.append(row)
    for n in CELLS:
        if n == "DM04":
            continue
        bad = [t[0] for t in TESTS if R[(n, t[0])][2] != "ok" or abs(R[(n, t[0])][3] / R[("DM04", t[0])][3] - 1) > 0.1
               or abs(R[(n, t[0])][4] / R[("DM04", t[0])][4] - 1) > 0.1]
        L.append(f"- {n}: {len(bad)} of 17 outside +-10 % {bad}")
    txt = "q_peak / q at 25 %, ratio to DM04 (both exponents 0.5)\n\n" + "\n".join(L)
    open(os.path.join(HERE, "out_probe_elastic_exponent.md"), "w").write(txt + "\n")
    print(txt)
