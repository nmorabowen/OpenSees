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


def deadend_chain(p0=20.0, n=120):
    """Review numerics item 2 (reviewer's p5.py): alpha = r just inside the
    bounding surface (rho_alpha 0.999) at p0, then PROPORTIONAL elastic
    compression; psi rises and the bounding surface shrinks around the fixed
    alpha. The oracle's chain, each increment from its own previous state."""
    u = [2 / math.sqrt(6), -1 / math.sqrt(6), -1 / math.sqrt(6), 0, 0, 0]

    def state(k):
        a = [k * x for x in u]
        return dict(sigma=[p0 * (1 + a[i]) if i < 3 else 0.0 for i in range(6)], alpha=a,
                    alpha_in=list(a), z=[0.0] * 6, e=E_INIT)

    def rho(st):
        r = integrate(State.from_voigt(st["sigma"], st["alpha"], st["z"], st["e"], st["alpha_in"]),
                      [0.0] * 6, CAMPAIGN, O)
        return r.summary()["rho_alpha_end"]
    lo, hi = 0.1, 2.0
    for _ in range(50):
        k = 0.5 * (lo + hi)
        if rho(state(k)) > 1.0:
            hi = k
        else:
            lo = k
    st = state(lo * 0.999)
    # proportional elastic direction: C^-1 sigma (K, G scale together with sqrt p)
    G = 264.32 * 101.0 * (2.97 - E_INIT) ** 2 / (1 + E_INIT) * math.sqrt(p0 / 101.0)
    K = 2.0 / 3.0 * (1 + NU) / (1 - 2 * NU) * G
    S = st["sigma"]
    pS = sum(S[:3]) / 3.0
    d = [(S[i] - pS) / (2 * G) + pS / (3 * K) if i < 3 else 0.0 for i in range(6)]
    m = max(abs(x) for x in d)
    d = [x / m for x in d]
    out = []
    cur = st
    for j in range(n):
        de = [1e-4 * (1 + j // 10) * x for x in d]
        c = dict(cur)
        c.update(kind="deadend", step=j, dstrain=de)
        c["ref"] = run(c)
        out.append(c)
        if c["ref"]["status"] != "ok":
            break
        r = c["ref"]
        cur = dict(sigma=r["sigma"], alpha=r["alpha"], alpha_in=r["alpha_in"], z=r["z"],
                   e=cur["e"] - (1 + E_INIT) * sum(de[:3]))
    return out


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
    # review of #871 (numerics): elastic predictor cases (iso Delta p, loading and
    # unloading, from alpha = r at three cone offsets), the convergence set
    # (K0 x 8 directions x 1e-4/1e-3), and the psi-driven "dead-end" chain.
    cone = math.sqrt(2.0 / 3.0) * 0.005
    u = [2 / math.sqrt(6), -1 / math.sqrt(6), -1 / math.sqrt(6), 0, 0, 0]
    for k in (0.5, 2.0, 1.2):
        for p0 in (20.0, 100.0):
            for eps in (5e-5, 1.15e-4, 3e-4):
                for sgn in (1, -1):
                    r0 = [k * cone * x for x in u]
                    sig = [p0 * (1 + r0[0]), p0 * (1 + r0[1]), p0 * (1 + r0[2]), 0, 0, 0]
                    cases.append(dict(kind="elastic", sigma=sig, alpha=r0, alpha_in=list(r0),
                                      z=[0.0] * 6, e=E_INIT, dstrain=[sgn * eps] * 3 + [0, 0, 0],
                                      k=k, p0=p0, eps=sgn * eps))
    cdirs = {"isoC": [1, 1, 1, 0, 0, 0], "triaxC": [1, -0.3, -0.3, 0, 0, 0],
             "volshear": [1, 1, 1, 0.5, 0, 0], "act": [0.3, 1, 0, 0, 0, 0],
             "pas": [1, -0.3, 0, 0, 0, 0], "ext": [-0.2, -1, -0.2, 0, 0, 0],
             "shear": [0, 0, 0, 1, 0, 0], "vcomp_sh": [0.5, 1, 0.5, 0.3, 0, 0]}
    for p0 in (5.0, 20.0, 100.0):
        for dn, d in cdirs.items():
            for delta in (1e-4, 1e-3):
                c = k0_state(p0)
                c.update(kind="conv", p0=p0, dir=dn, delta=delta, dstrain=[delta * x for x in d])
                cases.append(c)
    for i, c in enumerate(cases):
        c["ref"] = run(c)
        if i % 50 == 0:
            print(i, len(cases), c["kind"], c["ref"]["status"], flush=True)
    cases += deadend_chain()
    out = os.path.join(ROOT, "tests", "data", "wp129_uw_model_reference.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        json.dump(dict(oracle="WP-134 sanisand_reference, preset uw_model, Radau rtol 1e-10",
                       cases=cases), fh)
    print("wrote", out, len(cases))


if __name__ == "__main__":
    main()
