"""WP-150 GATE 0 (a), oracle side: DM04's own Toyoura simulations of Verdugo & Ishihara (1996), re-run EXACTLY.

Dafalias & Manzari (2004) JEM 130(6):622, Table 1 constants (verified against the PDF, p. 5 = journal p. 626):
G0 125, nu 0.05, M 1.25, c 0.712, lambda_c 0.019, e0 0.934, xi 0.7, m 0.01, h0 7.05, ch 0.968, nb 1.1, A0 0.704,
nd 3.5, zmax 4, cz 600. p_at is not tabulated; 100 kPa is used (the OpenSees example; an explicit choice).
Test matrix = the paper's Figs. 5-9 (journal pp. 629-630), all triaxial compression from isotropic states, loaded
to 25 % axial strain (the unloading branches are not reproduced):
  Fig. 5  undrained, e 0.735, p0 100 / 1000 / 2000 / 3000
  Fig. 6  undrained, e 0.833, p0 100 / 1000 / 2000 / 3000
  Fig. 7  undrained, e 0.907, p0 100 / 1000 / 2000      (panel b is a duplicate of Fig. 6b in the journal)
  Fig. 8  drained, p0 500, e0 0.960 / 0.886 / 0.810
  Fig. 9  drained, p0 100, e0 0.996 / 0.917 / 0.831
Two option sets: `paper` (the DM04 equations as published: G with the current e, de = -(1+e) dev) and `uw_model`
(the UW additions U1-U5 the fork's C++ carries: G with e_init, de = -(1+e_init) dev, the low-p D sigmoid, p_min).
Writes out_gate0_oracle.json: {optset: {test: {eps_a, p, q, e, ev}}}.
"""
import json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "Ladruno_scripts"))
from dataclasses import replace  # noqa: E402
from sanisand_reference import TOYOURA, Options  # noqa: E402
from sanisand_reference.driver import triaxial  # noqa: E402
from sanisand_reference.ring import ring_variants  # noqa: E402

TESTS = ([(f"F5_u_e0.735_p{p}", False, 0.735, float(p)) for p in (100, 1000, 2000, 3000)]
         + [(f"F6_u_e0.833_p{p}", False, 0.833, float(p)) for p in (100, 1000, 2000, 3000)]
         + [(f"F7_u_e0.907_p{p}", False, 0.907, float(p)) for p in (100, 1000, 2000)]
         + [(f"F8_d_p500_e{e}", True, e, 500.0) for e in (0.960, 0.886, 0.810)]
         + [(f"F9_d_p100_e{e}", True, e, 100.0) for e in (0.996, 0.917, 0.831)])


def main():
    out = {}
    for opt_name in ("paper", "uw_model"):
        out[opt_name] = {}
        for label, drained, e0, p0 in TESTS:
            # e_init matters only under the UW options (U4/U5): it is the test's own initial e
            P = replace(TOYOURA, e_init=e0)
            O = Options() if opt_name == "paper" else ring_variants()["uw_model"]
            t0 = time.time()
            res, tab = triaxial(P, p0, e0, 0.25, drained=drained, O=O, n_out=400, rtol=1e-8)
            rec = dict(status=res.status, eps_a=[t["eps"][0] for t in tab], p=[t["p"] for t in tab],
                       q=[t["q"] for t in tab], e=[t["e"] for t in tab], ev=[sum(t["eps"][:3]) for t in tab])
            out[opt_name][label] = rec
            print(f"{opt_name:8s} {label:18s} {res.status:14s} q_peak {max(rec['q']):7.1f} q_end {rec['q'][-1]:7.1f} "
                  f"p_end {rec['p'][-1]:7.1f} e_end {rec['e'][-1]:.4f} ({time.time()-t0:.1f}s)", flush=True)
    json.dump(out, open(os.path.join(HERE, "out_gate0_oracle.json"), "w"))


if __name__ == "__main__":
    main()
