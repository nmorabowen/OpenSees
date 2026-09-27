"""WP-134 validation (1)+(2): elastic sanity, critical-state limit, and DM04's
qualitative Toyoura behaviour (drained / undrained triaxial compression).

    python Ladruno_files/testbed/sanisand_reference/run_paper.py
Writes out/paper.md and out/paper_paths.json (sampled paths for plotting)."""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "Ladruno_scripts"))

import numpy as np  # noqa: E402

from sanisand_reference import TOYOURA, CAMPAIGN, Options, State, e_critical, elastic_moduli, integrate  # noqa: E402
from sanisand_reference.driver import triaxial  # noqa: E402

OUT = os.path.join(HERE, "out")
TESTS = [  # (label, drained, e0, p0)
    ("dense drained", True, 0.735, 100.0),
    ("loose drained", True, 0.95, 100.0),
    ("e=0.833 undrained p0=100", False, 0.833, 100.0),
    ("e=0.833 undrained p0=1000", False, 0.833, 1000.0),
    ("e=0.833 undrained p0=2000", False, 0.833, 2000.0),
    ("e=0.833 undrained p0=3000", False, 0.833, 3000.0),
    ("e=0.907 undrained p0=100", False, 0.907, 100.0),
    ("e=0.735 undrained p0=100", False, 0.735, 100.0),
]


def main():
    os.makedirs(OUT, exist_ok=True)
    L = ["## Elastic sanity", ""]
    O = Options(g_void_ratio="initial")
    p0, ev = 20.0, 3.0e-3
    st = State(p0 * np.eye(3), np.zeros((3, 3)), np.zeros((3, 3)), 0.7, np.zeros((3, 3)))
    r = integrate(st, [ev / 3] * 3 + [0, 0, 0], CAMPAIGN, O)
    _, K0 = elastic_moduli(p0, 0.7, CAMPAIGN, O)
    k = K0 / math.sqrt(p0)
    exact = (math.sqrt(p0) + 0.5 * k * ev) ** 2
    L.append(f"Isotropic compression, campaign set, p0 = 20 kPa, eps_v = 3e-3 (G at e_init so K = k sqrt(p)): "
             f"reference p = {r.end['p']:.10g} kPa, closed form {(exact):.10g} kPa, rel. error "
             f"{abs(r.end['p'] - exact) / exact:.1e}; one elastic segment, alpha untouched.")
    L += ["", "## DM04 Toyoura set, triaxial compression to 30 % axial strain (rtol 1e-8)", "",
          "| test | ψ0 | status | q_peak | p_min | p_end | q_end | η_end/M | ψ_end | ε_v end | max ρ_α | max \\|f\\| plastic | wall s |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    paths = {}
    for label, drained, e0, pp in TESTS:
        t0 = time.time()
        res, tab = triaxial(TOYOURA, pp, e0, 0.3, drained=drained, rtol=1e-8)
        wall = time.time() - t0
        last = tab[-1]
        psi0 = e0 - e_critical(pp, TOYOURA)
        L.append(f"| {label} | {psi0:+.3f} | {res.status} | {max(t['q'] for t in tab):.1f} | "
                 f"{min(t['p'] for t in tab):.1f} | {last['p']:.1f} | {last['q']:.1f} | "
                 f"{last['eta'] / TOYOURA.Mc:.4f} | {last['psi']:+.4f} | {sum(last['eps'][:3]):+.4f} | "
                 f"{res.max_rho_alpha:.3f} | {res.max_abs_f_plastic:.1e} | {wall:.1f} |")
        paths[label] = [dict(eps_a=t["eps"][0], p=t["p"], q=t["q"], e=t["e"],
                             eps_v=sum(t["eps"][:3])) for t in tab]
    txt = "\n".join(L) + "\n"
    open(os.path.join(OUT, "paper.md"), "w", encoding="utf-8").write(txt)
    json.dump(paths, open(os.path.join(OUT, "paper_paths.json"), "w"))
    print(txt)


if __name__ == "__main__":
    main()
