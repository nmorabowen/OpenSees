"""WP-129: generate the WP-134 oracle fixture for the SAS-ME tests.

Runs WP-134's independent reference integrator (`sanisand_reference`, draft PR
#872: DM04 rate equations, SciPy Radau rtol 1e-10 with event detection) with
the `uw_model` preset -- the UW constitutive additions U1-U5, the PAPER's
alpha_in rule, continuous moduli: the target a corrected C++ integrator should
reproduce -- on the cases the SAS-ME tests replay, and writes
tests/data/wp129_uw_model_reference.json.

Needs numpy + scipy (CPython 3.11 on the dev box; the 3.12 test runner has no
scipy, which is why the tests read a fixture instead of calling the oracle):

    SANISAND_REF=<dir holding the sanisand_reference package> \
        py -3.11 Ladruno_files/testbed/wp129_sasme/gen_uw_model_reference.py
"""
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.environ.get("SANISAND_REF", os.path.join(ROOT, "Ladruno_scripts")))

from sanisand_reference import CAMPAIGN, State, integrate  # noqa: E402
from sanisand_reference.ring import ring_variants  # noqa: E402

O = ring_variants()["uw_model"]
NU = 0.312885
E_INIT = 0.6944


def k0_state(p0):
    K0 = NU / (1.0 - NU)
    sv = 3.0 * p0 / (1.0 + 2.0 * K0)
    sig = [K0 * sv, sv, K0 * sv, 0.0, 0.0, 0.0]
    p = sum(sig[:3]) / 3.0
    al = [sig[0] / p - 1, sig[1] / p - 1, sig[2] / p - 1, 0.0, 0.0, 0.0]
    return dict(sigma=sig, alpha=al, alpha_in=list(al), z=[0.0] * 6, e=E_INIT)


def run(case):
    st = State.from_voigt(case["sigma"], case["alpha"], case["z"], case["e"], case["alpha_in"])
    r = integrate(st, case["dstrain"], CAMPAIGN, O)
    s = r.summary()
    return {k: s[k] for k in ("status", "t_end", "f_end", "max_rho_alpha", "rho_alpha_end",
                              "max_rho_b", "rho_b_end", "eta_end", "p_end", "sigma", "alpha",
                              "z", "alpha_in", "n_reseats")}


def probes(delta):
    return {"isoComp": [delta, delta, 0, 0, 0, 0], "isoExt": [-delta, -delta, 0, 0, 0, 0],
            "shear+": [0, 0, 0, delta, 0, 0], "shear-": [0, 0, 0, -delta, 0, 0]}


def main():
    cases = []
    cases.append(dict(kind="reproducer", sigma=[0.0101] * 3 + [0.0] * 3, alpha=[0.0] * 6,
                      alpha_in=[0.0] * 6, z=[0.0] * 6, e=0.697787979641054,
                      dstrain=[0, 1e-4, 0, 0, 0, 0]))
    dirs = {"active": [0.3, 1.0, 0, 0, 0, 0], "passive": [1.0, -0.3, 0, 0, 0, 0],
            "shear": [0, 0, 0, 1.0, 0, 0]}
    for p0 in (20.0, 50.0, 100.0):
        for dn, d in dirs.items():
            for delta in (1e-5, 1e-4):
                c = k0_state(p0)
                c.update(kind="benign", p0=p0, dir=dn, delta=delta,
                         dstrain=[delta * x for x in d])
                cases.append(c)
    att = os.path.join(ROOT, "Ladruno_implementation", "_tims_2d_model_requests_2026-09-25")
    for mesh in ("b8", "b16"):
        with open(os.path.join(att, f"ring_points_{mesh}.csv"), newline="") as fh:
            for row in csv.DictReader(fh):
                el, gp = int(row["element"]), int(row["gp"])
                if mesh == "b8" and el == 1950 and gp in (2, 3):
                    continue        # inadmissible: refused, nothing to compare
                base = dict(kind="ring", mesh=mesh, element=el, gp=gp,
                            sigma=[float(row[f"sigma_{i}"]) for i in range(6)],
                            alpha=[float(row[f"alpha_{i}"]) for i in range(6)],
                            alpha_in=[float(row[f"alpha_in_{i}"]) for i in range(6)],
                            z=[float(row[f"z_{i}"]) for i in range(6)], e=float(row["e"]))
                for delta in (1e-6, 1e-5):
                    for pn, de in probes(delta).items():
                        c = dict(base)
                        c.update(probe=pn, delta=delta, dstrain=de)
                        cases.append(c)
    for i, c in enumerate(cases):
        c["ref"] = run(c)
        if i % 50 == 0:
            print(i, len(cases), c["kind"], c["ref"]["status"], flush=True)
    out = os.path.join(ROOT, "tests", "data", "wp129_uw_model_reference.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        json.dump(dict(oracle="WP-134 sanisand_reference, preset uw_model, Radau rtol 1e-10",
                       cases=cases), fh)
    print("wrote", out, len(cases))


if __name__ == "__main__":
    main()
