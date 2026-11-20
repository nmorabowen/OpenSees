"""WP-229 (TIMs F26) gate, oracle side: DM04's 17 Toyoura triaxial simulations of Verdugo &
Ishihara (1996) with the plastic-modulus pressure exponent nh in b0 = G0 h0 (1 - ch e) (p/p_at)^-nh.

The matrix is WP-150's gate 0 (`gate0_toyoura_oracle.py`, branch wp/150-sanisand-regularization-memo):
DM04 Figs. 5-9, triaxial compression from isotropic states to 25 % axial strain. Constants: DM04 Table 1
(`sanisand_reference.TOYOURA`, p_at 100 kPa), `uw_model` options (the oracle the C++ SAS-ME targets).

Gate: for nh in 0.2 / 0.25 / 0.3, q_peak and q at 25 % strain within +-10 % of the nh = 0.5 result.
nh = 0.5 takes the original sqrt branch, so it IS the unpatched oracle.

  python gate_nh_oracle.py            -> out_gate_nh_oracle.json, out_gate_nh_oracle.md
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "Ladruno_scripts"))

NH = (0.5, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
NSTEP = 500
TESTS = ([(f"F5_u_e0.735_p{p}", False, 0.735, float(p)) for p in (100, 1000, 2000, 3000)]
         + [(f"F6_u_e0.833_p{p}", False, 0.833, float(p)) for p in (100, 1000, 2000, 3000)]
         + [(f"F7_u_e0.907_p{p}", False, 0.907, float(p)) for p in (100, 1000, 2000)]
         + [(f"F8_d_p500_e{e}", True, e, 500.0) for e in (0.960, 0.886, 0.810)]
         + [(f"F9_d_p100_e{e}", True, e, 100.0) for e in (0.996, 0.917, 0.831)])


def one(job):
    nh, label, drained, e0, p0 = job
    from sanisand_reference import TOYOURA
    import math
    import numpy as np
    from sanisand_reference.driver import isotropic_state
    from sanisand_reference.integrator import Control, integrate
    from sanisand_reference.model import norm
    from sanisand_reference.ring import ring_variants
    P = replace(TOYOURA, e_init=e0)
    # start_outside="plastic": a restart from the previous end state carries that segment's consistency
    # drift (about 1e-8 of the cone radius, i.e. > the absolute 1e-7 kPa start tolerance at MPa pressures);
    # it is on the surface and is integrated as such.
    O = ring_variants()["uw_model"].with_(plastic_exp=nh, start_outside="plastic")
    t0 = time.time()
    # NSTEP equal increments, each its own integrate() call from the previous end state, so the curve is
    # sampled at the same strains as the C++ side (gate_nh_cxx.py) and the two can be compared point by point.
    da = 0.25 / NSTEP
    st = isotropic_state(p0, e0)
    if drained:
        ctl = Control((True, False, False, True, True, True), (da, 0.0, 0.0, 0.0, 0.0, 0.0))
    else:
        ctl = Control.strain([da, -0.5 * da, -0.5 * da, 0.0, 0.0, 0.0])
    eps, q, p, status, fmax = [0.0], [0.0], [p0], "ok", 0.0
    for k in range(NSTEP):
        r = integrate(st, ctl, P, O, rtol=1e-8, record=False)
        if r.status != "ok":
            status = f"{r.status} at step {k}"
            break
        st = r.state
        pk = float(np.trace(st.sigma)) / 3.0
        fmax = max(fmax, abs(r.f_end) / (math.sqrt(2.0 / 3.0) * P.m * pk))
        eps.append((k + 1) * da)
        p.append(pk)
        q.append(math.sqrt(1.5) * norm(st.sigma - pk * np.eye(3)))
    return dict(nh=nh, label=label, status=status, q_peak=max(q), q_end=q[-1],
                eps_peak=eps[q.index(max(q))], p_end=p[-1], e_end=st.e, eps_a=eps, q=q, p=p,
                max_f_rel=fmax,
                wall=time.time() - t0)


def main():
    jobs = [(nh,) + t for nh in NH for t in TESTS]
    with ProcessPoolExecutor(max_workers=min(12, os.cpu_count() or 4)) as ex:
        recs = list(ex.map(one, jobs))
    by = {(r["nh"], r["label"]): r for r in recs}
    json.dump(recs, open(os.path.join(HERE, "out_gate_nh_oracle.json"), "w"))
    lines = ["| test | nh 0.5: q_peak / q_25% | " + " | ".join(
        f"nh {nh}: q_peak ratio / q_25% ratio" for nh in NH[1:]) + " |",
        "|---|---|" + "---|" * (len(NH) - 1)]
    worst, bad = 0.0, []
    for label, *_ in TESTS:
        ref = by[(0.5, label)]
        row = f"| {label} | {ref['q_peak']:.1f} / {ref['q_end']:.1f} ({ref['status']}) |"
        for nh in NH[1:]:
            r = by[(nh, label)]
            rp, re_ = r["q_peak"] / ref["q_peak"], r["q_end"] / ref["q_end"]
            ok = r["status"] == "ok" and ref["status"] == "ok"
            worst = max(worst, abs(rp - 1), abs(re_ - 1))
            if not ok or abs(rp - 1) > 0.10 or abs(re_ - 1) > 0.10:
                bad.append((nh, label, rp, re_, r["status"]))
            row += f" {rp:.4f} / {re_:.4f}{'' if ok else ' (' + r['status'] + ')'} |"
        lines.append(row)
    lines.append("")
    lines.append(f"largest |f| / (sqrt(2/3) m p) at an increment end, all runs: "
                 f"{max(r['max_f_rel'] for r in recs):.2e}")
    lines.append(f"worst |ratio - 1| = {worst:.4f}; outside +-10 % or not ok: {len(bad)} of {17 * (len(NH) - 1)}")
    for b in bad:
        lines.append(f"- nh {b[0]} {b[1]}: q_peak x{b[2]:.4f}, q_25% x{b[3]:.4f}, status {b[4]}")
    txt = "\n".join(lines)
    open(os.path.join(HERE, "out_gate_nh_oracle.md"), "w").write(txt + "\n")
    print(txt)




# ---------------------------------------------------------------------------------------------------
# tests/data/wp229_oracle_nh.json: short STRAIN-controlled undrained triaxial chains (no mixed control,
# so the C++ side is a plain chain of replays), Toyoura, uw_model, at nh 0.5 and 0.25, from p0 10 / 100 /
# 1000 kPa.  40 increments of 0.05 % axial strain.  tests/test_ladruno_sanisand_plastic_exp.py compares
# the C++ against it.          python gate_nh_oracle.py --fixture
def fixture():
    import math
    import numpy as np
    from sanisand_reference import TOYOURA
    from sanisand_reference.driver import isotropic_state
    from sanisand_reference.integrator import Control, integrate
    from sanisand_reference.model import norm, t2v
    from sanisand_reference.ring import ring_variants
    da, n = 5.0e-4, 40
    out = dict(note="WP-229 oracle expectations: gate_nh_oracle.py --fixture (WP-134 oracle, uw_model, "
                    "Toyoura DM04 Table 1, p_at 100, p_min 0.0101; undrained triaxial, strain controlled)",
               params=TOYOURA.as_opensees(), da=da, n=n, paths={})
    for e0, p0 in ((0.833, 10.0), (0.833, 100.0), (0.735, 1000.0)):
        for nh in (0.5, 0.25):
            P = replace(TOYOURA, e_init=e0)
            O = ring_variants()["uw_model"].with_(plastic_exp=nh, start_outside="plastic")
            st = isotropic_state(p0, e0)
            ctl = Control.strain([da, -0.5 * da, -0.5 * da, 0.0, 0.0, 0.0])
            q, p = [], []
            for k in range(n):
                r = integrate(st, ctl, P, O, rtol=1e-9, record=False)
                assert r.status == "ok", (e0, p0, nh, k, r.status)
                st = r.state
                pk = float(np.trace(st.sigma)) / 3.0
                p.append(pk)
                q.append(math.sqrt(1.5) * norm(st.sigma - pk * np.eye(3)))
            out["paths"][f"e{e0}_p{p0:g}_nh{nh}"] = dict(e0=e0, p0=p0, nh=nh, p=p, q=q,
                                                          sigma_end=t2v(st.sigma).tolist())
            print(f"e {e0} p0 {p0:g} nh {nh}: q_end {q[-1]:.4f} p_end {p[-1]:.4f}")
    json.dump(out, open(os.path.join(ROOT, "tests", "data", "wp229_oracle_nh.json"), "w"), indent=0)


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "--fixture":
    fixture()
elif __name__ == "__main__":
    main()
